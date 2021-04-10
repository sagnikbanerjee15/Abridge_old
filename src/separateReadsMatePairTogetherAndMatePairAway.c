# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include <math.h>
# include <stdbool.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

void writeToFile (char **split_line, FILE *fhw)
{
	fprintf (fhw , "%s" , split_line[0]);
	fprintf (fhw , "%s" , "\t");
	fprintf (fhw , "%s" , split_line[1]);
	//fflush (fhw);
}

void splitMappingInTwoPartsAndSetNHValue (char *line1, char **split_line1, int *NH_value1)
{
	int i, j0, j1, first_tab_found = 0, k, NH_string_index;
	char NH_string[100];
	char *temp;

	j0 = 0;
	j1 = 0;
	NH_string_index = 0;
	for ( i = 0 ; line1[i] != '\0' ; i++ )
	{
		if ( line1[i] == '\t' && first_tab_found == 0 )
		{
			first_tab_found = 1;
			continue;
		}
		if ( first_tab_found )
			split_line1[1][j1++ ] = line1[i];
		else split_line1[0][j0++ ] = line1[i];

		if ( line1[i - 1] == ':' && line1[i - 2] == 'i' && line1[i - 3] == ':' && line1[i - 4] == 'H' && line1[i - 5] == 'N' )
		{
			for ( k = i ; line1[k] != '\t' ; k++ )
				NH_string[NH_string_index++ ] = line1[k];
			NH_string[NH_string_index++ ] = '\0';
		}
	}
	( *NH_value1 ) = strtol (NH_string , &temp , 10);
	split_line1[0][j0] = '\0';
	split_line1[1][j1] = '\0';
}

void separateReads (char *input_samfilename, char *output_filename_mate_pairs_already_next_to_each_other, char *output_filename_mate_pairs_away_from_to_each_other)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char *line1 = NULL; // for reading each line1
	char **split_line1; // List of strings to store each element of a single alignment
	char *line2 = NULL; // for reading each line1
	char **split_line2; // List of strings to store each element of a single alignment

	int NH_value1, NH_value2;
	int i;

	size_t len1 = 0;
	ssize_t len2 = 0;
	ssize_t line1_len;
	ssize_t line2_len;

	FILE *fhr;
	FILE *fhw_already_next_to_each_other;
	FILE *fhw_away_from_to_each_other;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr = fopen (input_samfilename , "r");
	if ( fhr == NULL )
	{
		printf ("Error! File %s not found" , input_samfilename);
		exit (1);
	}
	fhw_already_next_to_each_other = fopen (output_filename_mate_pairs_already_next_to_each_other , "w");
	if ( fhw_already_next_to_each_other == NULL )
	{
		printf ("Error! File %s cannot be created" , fhw_already_next_to_each_other);
		exit (1);
	}
	fhw_away_from_to_each_other = fopen (output_filename_mate_pairs_away_from_to_each_other , "w");
	if ( fhw_away_from_to_each_other == NULL )
	{
		printf ("Error! File %s cannot be created" , fhw_away_from_to_each_other);
		exit (1);
	}

	split_line1 = ( char** ) malloc (sizeof(char*) * 10);
	for ( i = 0 ; i < 10 ; i++ )
		split_line1[i] = ( char* ) malloc (sizeof(char) * 1000);

	split_line2 = ( char** ) malloc (sizeof(char*) * 10);
	for ( i = 0 ; i < 10 ; i++ )
		split_line2[i] = ( char* ) malloc (sizeof(char) * 1000);
	/********************************************************************/

	while ( ( line1_len = getline ( &line1 , &len1 , fhr) ) != -1 )
	{
		if ( line1[0] == '@' )
		{
			fprintf (fhw_already_next_to_each_other , "%s" , line1);
			fprintf (fhw_away_from_to_each_other , "%s" , line1);
		}
		else break;
	}

	do
	{
		line2_len = getline ( &line2 , &len2 , fhr);
		if ( line2_len == -1 ) break;
		splitMappingInTwoPartsAndSetNHValue (line1 , split_line1 , &NH_value1);
		splitMappingInTwoPartsAndSetNHValue (line2 , split_line2 , &NH_value2);
		if ( strcmp (split_line1[0] , split_line2[0]) == 0 )
		{
			/*
			 * Same read names for consecutive reads
			 */
			writeToFile (split_line1 , fhw_already_next_to_each_other);
			writeToFile (split_line2 , fhw_already_next_to_each_other);
			line1_len = getline ( &line1 , &len1 , fhr);
			if ( line1_len == -1 ) break;
		}
		else
		{
			writeToFile (split_line1 , output_filename_mate_pairs_away_from_to_each_other);
			strcpy(line1 , line2);
		}
	} while ( 1 );
}

int main (int argc, char *argv[])
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	int read_length;
	int read_id[100];
	int i, j;

	char input_samfilename[FILENAME_LENGTH];
	char output_filename_mate_pairs_already_next_to_each_other[FILENAME_LENGTH];
	char output_filename_mate_pairs_away_from_to_each_other[FILENAME_LENGTH];

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	read_length = 0;
	strcpy(input_samfilename , argv[1]);
	strcpy(output_filename_mate_pairs_already_next_to_each_other , argv[2]);
	strcpy(output_filename_mate_pairs_away_from_to_each_other , argv[3]);

	/********************************************************************/
	separateReads (input_samfilename , output_filename_mate_pairs_already_next_to_each_other , output_filename_mate_pairs_away_from_to_each_other);

	return 0;
}
