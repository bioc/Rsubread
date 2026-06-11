#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
int cell_counts_extract_vHD_main(char * bamname, char * txtname){
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    
    HashMap *known_map = create_map();
    char command[100+MAX_FILE_NAME_LENGTH];
    sprintf(command, "samtools view %s", bamname);
    FILE *samfp = popen(command, "r");
    FILE *txtfp = fopen(txtname, "w");
    
    while ((read = getline(&line, &len, samfp)) != -1) {
        // Strip trailing newline character if present
        if (read > 0 && line[read - 1] == '\n') {
            line[read - 1] = '\0';
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
                } else if (strncmp(token, "CB:Z:", 5) == 0) {
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
    fclose(samfp);
    fclose(txtfp);

    free(line);
    free_map(known_map);
    return 0;
}
