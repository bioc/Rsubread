#include "subread.h"
#include "HelperFunctions.h"

int RSubread_parse_CIGAR_string(const char * CIGAR_Str, unsigned int * Staring_Points, unsigned short * Section_Length)
{
	unsigned int tmp_int=0;
	int cigar_cursor=0;
	unsigned short read_cursor=0;
	unsigned int chromosome_cursor=0;
	int ret=0;

	for(cigar_cursor=0; ; cigar_cursor++)
	{
		char ch = CIGAR_Str[cigar_cursor];

		if(ch >='0' && ch <= '9')
		{
			tmp_int=tmp_int*10+(ch - '0');
		}
		else
		{
			if(ch == 'M' || ch == 'D')
			{
				read_cursor += tmp_int;
				chromosome_cursor += tmp_int;
			}
			else if(ch == 'N' || ch == 0)
			{
				if(ret <6)
				{
					if(read_cursor>0)
					{
						Staring_Points[ret] = chromosome_cursor - read_cursor;
						Section_Length[ret] = read_cursor;
						ret ++;
					}
				}
				read_cursor = 0;

				if(ch == 'N') chromosome_cursor += tmp_int;
				else break;
			}
			//printf("C=%c, TV=%d, CC=%d, RC=%d\n", ch, tmp_int, chromosome_cursor, read_cursor);
			tmp_int = 0;
		}
		if(cigar_cursor>100) return -1;
	}

	return ret;
}

void display_sections(char * CIGAR_Str)
{
	unsigned int Staring_Points[6];
	unsigned short Section_Length[6];
	int retv = RSubread_parse_CIGAR_string(CIGAR_Str, Staring_Points, Section_Length);

	int x1;
	SUBREADprintf("Cigar=%s ; Sections=%d\n", CIGAR_Str, retv);
	for(x1=0; x1<retv; x1++)
	{
		SUBREADprintf("   Section #%d: offset=%u  length=%u\n",x1, Staring_Points[x1], Section_Length[x1]);
	}
	SUBREADprintf("\n");
	
}

#ifdef RSUBREAD_TEST_HELPER_FUNCTIONS
void main()
#else
void testi_helper_1_main()
#endif
{
	display_sections("");
	display_sections("*");
	display_sections("5S10M2D10M800N12M3I12M450N12M12D99M6S");
	display_sections("110M2I10M800N32M3I12M6S");
	display_sections("200S110M2I10M800N32M3I12M200N40M");
	display_sections("3M1663N61M1045N36M3D20M66N10M2D10M77N3M1663N61M1045N36M3D20M66N103M1663N61M1045N36M3D20M66N9M");
}
