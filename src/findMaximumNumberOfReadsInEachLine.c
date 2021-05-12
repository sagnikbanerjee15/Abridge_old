# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include <math.h>
# include <stdbool.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

void str_reverse (char *str)
{
	if ( !str ) return;

	int i, j;
	int length_of_string;
	char ch;

	length_of_string = strlen (str);
	for ( i = 0 ; i < length_of_string / 2 ; i++ )
	{
		ch = str[i];
		str[i] = str[length_of_string - i - 1];
		str[length_of_string - i - 1] = ch;
	}
}

int main (int argc, char *argv[])
{
	/****************************************************************************************************************************************
	 * Variable declaration
	 ****************************************************************************************************************************************/
	char pass1_filename[FILENAME_LENGTH];

	char *buffer_for_pass1;

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
	strcpy(pass1_filename , argv[1]);
	fhr = fopen (pass1_filename , "r");
	if ( fhr == NULL )
	{
		printf ("Error! File pass1_filename = %s not found" , pass1_filename);
		exit (1);
	}
	/****************************************************************************************************************************************/
	while ( ( line_len = getline ( &buffer_for_pass1 , &len , fhr) ) != -1 )
	{
		str_reverse (buffer_for_pass1);
		reads_in_this_line = 0;
		for ( i = 0 ; buffer_for_pass1[i] != ' ' ; i++ )
			if ( buffer_for_pass1[i] == ',' ) reads_in_this_line += 1;
		if ( reads_in_this_line > max_reads_in_a_line )
			max_reads_in_a_line = reads_in_this_line;
	}
	printf ("%d\n" , max_reads_in_a_line);
	return 0;
}
