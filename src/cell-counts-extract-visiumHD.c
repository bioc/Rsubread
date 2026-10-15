#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "subread.h"

#define INITIAL_CAPACITY (33554467) // prime number

// Structure for the hash map to mimic AWK's 'known' associative array
typedef struct Node {
    char *key; // Stores concatenated "r y"
    struct Node *next;
} Node;

typedef struct {
    Node **buckets;
    size_t capacity;
    size_t size;
} HashMap;

// Simple DJB2-like hash function for strings
unsigned long hash(const char *str, size_t capacity) {
    unsigned long h = 5381;
    int c;
    while ((c = *str++)) {
        h = ((h << 5) + h) + c;
    }
    return h % capacity;
}

HashMap* create_map() {
    HashMap *map = malloc(sizeof(HashMap));
    map->capacity = INITIAL_CAPACITY;
    map->size = 0;
    map->buckets = calloc(map->capacity, sizeof(Node*));
    return map;
}

// Inserts a key. Returns 1 if it's a newly inserted key, 0 if it already existed.
int insert_key(HashMap *map, const char *key) {
    unsigned long slot = hash(key, map->capacity);
    Node *curr = map->buckets[slot];
    
    // Check if key already exists
    while (curr != NULL) {
        if (strcmp(curr->key, key) == 0) {
            return 0; // Already known
        }
        curr = curr->next;
    }
    
    // Key is unique; insert it
    Node *new_node = malloc(sizeof(Node));
    new_node->key = strdup(key);
    new_node->next = map->buckets[slot];
    map->buckets[slot] = new_node;
    map->size++;
    
    return 1; // Successfully added new key
}

void free_map(HashMap *map) {
    for (size_t i = 0; i < map->capacity; i++) {
        Node *curr = map->buckets[i];
        while (curr != NULL) {
            Node *tmp = curr;
            curr = curr->next;
            free(tmp->key);
            free(tmp);
        }
    }
    free(map->buckets);
    free(map);
}


#define UMI_LENGTH_IN_1R 9	// skip the UMI and index the barcodes
#if defined(_WIN32) || defined(_WIN64)
    #define popen _popen
    #define pclose _pclose
#endif

#define BAM_RECORD_LINE_LENGTH (MAX_READ_NAME_LEN + MAX_READ_LENGTH*2+MAX_CHROMOSOME_NAME_LEN*2+8192)

static int shell_quote_argument(const char *argument, char **quoted_argument) {
    size_t quote_length = 2;
    const char *input;
    char *quoted, *output;

    for (input = argument; *input; input++) {
        size_t additional = (*input == '\'') ? 4 : 1;
        if (quote_length > SIZE_MAX - additional) return -1;
        quote_length += additional;
    }
    if (quote_length == SIZE_MAX) return -1;

    quoted = malloc(quote_length + 1);
    if (!quoted) return -1;

    output = quoted;
    *output++ = '\'';
    for (input = argument; *input; input++) {
        if (*input == '\'') {
            *output++ = '\'';
            *output++ = '\\';
            *output++ = '\'';
            *output++ = '\'';
        } else {
            *output++ = *input;
        }
    }
    *output++ = '\'';
    *output = '\0';
    *quoted_argument = quoted;
    return 0;
}

int cell_counts_extract_vHD_main(char * bamname, char * txtname){
    // Allocate a fixed buffer on the stack for fgets
    char line[BAM_RECORD_LINE_LENGTH]; 

    HashMap *known_map = create_map();
    char *quoted_bamname = NULL;
    char *command = NULL;
    size_t command_capacity;
    int command_length;

    if (shell_quote_argument(bamname, &quoted_bamname) != 0) {
        free_map(known_map);
        return -1;
    }
    if (strlen(quoted_bamname) > SIZE_MAX - sizeof("samtools view ")) {
        free(quoted_bamname);
        free_map(known_map);
        return -1;
    }
    command_capacity = sizeof("samtools view ") + strlen(quoted_bamname);
    command = malloc(command_capacity);
    if (!command) {
        free(quoted_bamname);
        free_map(known_map);
        return -1;
    }
    command_length = snprintf(command, command_capacity, "samtools view %s", quoted_bamname);
    if (command_length < 0 || (size_t)command_length >= command_capacity) {
        free(command);
        free(quoted_bamname);
        free_map(known_map);
        return -1;
    }
    free(quoted_bamname);
    
    FILE *samfp = popen(command, "r");
    free(command);
    FILE *txtfp = fopen(txtname, "w");

    if (!samfp || !txtfp) {
        // Basic error checking to prevent crashing if files fail to open
        if (samfp) pclose(samfp);
        if (txtfp) fclose(txtfp);
        free_map(known_map);
        return -1;
    }

    // Loop using fgets instead of getline
    while (fgets(line, sizeof(line), samfp) != NULL) {
        
        // Strip trailing newline character if present
        size_t len = strlen(line);
        if (len > 0 && line[len - 1] == '\n') {
            line[len - 1] = '\0';
        }

        // AWK: !/^@/
        if (line[0] == '@' || line[0] == '\0') {
            continue;
        }

        char *r = NULL;
        char *y = NULL;
        char *cb = NULL;

        char *line_ptr = line;
        char *token;
        char *saveptr;
        int field_idx = 1;

        // Tokenize line by tab character (\t)
        while ((token = strtok_r(line_ptr, "\t", &saveptr)) != NULL) {
            line_ptr = NULL; // subsequent calls to strtok_r need NULL

            // AWK: for(i=12; i<=NF; i++)
            if (field_idx >= 12) {
                // Check prefixes and extract substring from index 6 (0-indexed 5)
                if (strncmp(token, "1R:Z:", 5) == 0) {
                    r = token + 5 + UMI_LENGTH_IN_1R;
                } else if (strncmp(token, "1Y:Z:", 5) == 0) {
                    y = token + 5 + UMI_LENGTH_IN_1R;
                } else if (strncmp(token, "CB:Z:", 5) == 0 && !cb) {
                    cb = token + 5;
                } else if (strncmp(token, "sb:Z:", 5) == 0) {
                    cb = token + 5;
                }
            }
            field_idx++;
        }

        // AWK: if(r&&y&&cb)
        if (r && y && cb) {
            // AWK: known[r y] -> Create combined key
            size_t key_len = strlen(r) + strlen(y) + 1;
            char *combined_key = malloc(key_len);
            snprintf(combined_key, key_len, "%s%s", r, y);

            // AWK: if(!known[r y]) { print r, y, cb; known[r y]=1 }
            if (insert_key(known_map, combined_key)) {
                fprintf(txtfp, "%s\t%s\t%s\n", r, y, cb);
            }

            free(combined_key);
        }
    }
    
    // Crucial Windows Fix: streams opened with popen must be closed with pclose
    pclose(samfp); 
    fclose(txtfp);

    // Removed free(line) as 'line' is now a stack array, not dynamically allocated.
    free_map(known_map);
    return 0;
}
