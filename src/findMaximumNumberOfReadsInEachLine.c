# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include <math.h>
# include <stdbool.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

int main (int argc, char *argv[])
{
	/****************************************************************************************************************************************
	 * Variable declaration
	 ****************************************************************************************************************************************/
	char pass1_filename[FILENAME_LENGTH];

	char *buffer_for_pass1;
	char *rev_buffer_for_pass1;

	size_t len = 0;
	ssize_t line_len;

	FILE *fhr;

	int max_reads_in_a_line = 0;
	int reads_in_this_line = 0;
	int i;

	/****************************************************************************************************************************************/

	/****************************************************************************************************************************************
	 * Variable initialization
	 ****************************************************************************************************************************************/
	fhr = fopen (pass1_filename , "r");
	if ( fhr == NULL )
	{
		printf ("Error! File pass1_filename = %s not found" , pass1_filename);
		exit (1);
	}
	/****************************************************************************************************************************************/
	while ( ( line_len = getline ( &buffer_for_pass1 , &len , fhr) ) != -1 )
	{
		rev_buffer_for_pass1 = strrev (buffer_for_pass1);
		reads_in_this_line = 0;
		for ( i = 0 ; rev_buffer_for_pass1[i] != ' ' ; i++ )
			if ( rev_buffer_for_pass1[i] == ',' ) reads_in_this_line += 1;
		if ( reads_in_this_line > max_reads_in_a_line )
			max_reads_in_a_line = reads_in_this_line;
	}
	printf ("%d\n" , max_reads_in_a_line);
	return 0;
}
