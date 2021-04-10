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

void generateNextReadID (char *alphabets, int *read_id, int *read_length)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	int i;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/

	/********************************************************************/

	if ( *read_length == 0 )
	{
		read_id[0] = 0;
		( *read_length )++;
	}
	else
	{
		/*
		 * Check if the read is the last element of the maximum read_length
		 */
		for ( i = 0 ; i < *read_length ; i++ )
			if ( read_id[i] != strlen (alphabets) - 1 ) break;
		if ( i == *read_length )
		{
			( *read_length )++;
			for ( i = 0 ; i < *read_length ; i++ )
				read_id[i] = 0;
		}
		else
		{
			/*
			 * Increment the read_id
			 */
			for ( i = *read_length - 1 ; i >= 0 ; i-- )
			{
				if ( read_id[i] == ( strlen (alphabets) - 1 ) )
					read_id[i] = 0;
				else
				{
					read_id[i]++;
					break;
				}
			}
		}
	}
}

void convertReadIdToString (int *read_id, char *read_id_string, int read_length, char *alphabets)
{
	int i;
	for ( i = 0 ; i < read_length ; i++ )
		read_id_string[i] = alphabets[read_id[i]];
	read_id_string[i] = '\0';
}

void splitMappingInTwoPartsAndSetNHValue (char *line, char **split_line, int *NH_value)
{
	int i, j0, j1, first_tab_found = 0, k, NH_string_index;
	char NH_string[100];
	char *temp;

	j0 = 0;
	j1 = 0;
	NH_string_index = 0;
	for ( i = 0 ; line[i] != '\0' ; i++ )
	{
		if ( line[i] == '\t' && first_tab_found == 0 )
		{
			first_tab_found = 1;
			continue;
		}
		if ( first_tab_found )
			split_line[1][j1++ ] = line[i];
		else split_line[0][j0++ ] = line[i];

		if ( line[i - 1] == ':' && line[i - 2] == 'i' && line[i - 3] == ':' && line[i - 4] == 'H' && line[i - 5] == 'N' )
		{
			for ( k = i ; line[k] != '\t' ; k++ )
				NH_string[NH_string_index++ ] = line[k];
			NH_string[NH_string_index++ ] = '\0';
		}
	}
	( *NH_value ) = strtol (NH_string , &temp , 10);
	split_line[0][j0] = '\0';
	split_line[1][j1] = '\0';
}

void convertOldReadIdsToNewReadIds (char *input_samfilename_mate_pairs_next_to_each_other, char *input_samfilename_mate_pairs_away_from_each_other, char *output_samfilename)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char alphabets[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-='{}[]|?<>,.";
	char *temp; //Required for strtoi
	char *line1 = NULL; // for reading each line
	char *line2 = NULL; // for reading each line
	char **split_line1; // List of strings to store each element of a single alignment
	char **split_line2;
	char prev_old_read_id[1000];
	char read_id_string[100];

	int read_length;
	int i;
	int j;
	int read_id[100];
	int NH_value1;
	int NH_value2;
	int MAX_number_of_invalid_nodes_allowed = 10000;
	long long int read_number;

	size_t len1 = 0;
	ssize_t len2 = 0;
	ssize_t line1_len;
	ssize_t line2_len;

	FILE *fhr_mate_pairs_next_to_each_other;
	FILE *fhr_mate_pairs_away_from_each_other;
	FILE *fhw;

	struct Old_Read_ID_to_New_Read_ID_Circular_Linked_list *node_of_interest;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr_mate_pairs_next_to_each_other = fopen (input_samfilename_mate_pairs_next_to_each_other , "r");
	if ( fhr_mate_pairs_next_to_each_other == NULL )
	{
		printf ("Error! File %s not found" , input_samfilename_mate_pairs_next_to_each_other);
		exit (1);
	}
	fhr_mate_pairs_away_from_each_other = fopen (input_samfilename_mate_pairs_away_from_each_other , "r");
	if ( fhr_mate_pairs_away_from_each_other == NULL )
	{
		printf ("Error! File %s not found" , input_samfilename_mate_pairs_away_from_each_other);
		exit (1);
	}
	fhw = fopen (output_samfilename , "w");
	if ( fhw == NULL )
	{
		printf ("Error! File %s cannot be created" , output_samfilename);
		exit (1);
	}

	split_line1 = ( char** ) malloc (sizeof(char*) * 10);
	for ( i = 0 ; i < 10 ; i++ )
		split_line1[i] = ( char* ) malloc (sizeof(char) * 1000);

	split_line2 = ( char** ) malloc (sizeof(char*) * 10);
	for ( i = 0 ; i < 10 ; i++ )
		split_line2[i] = ( char* ) malloc (sizeof(char) * 1000);

	head = NULL;
	/********************************************************************/

	while ( ( line1_len = getline ( &line1 , &len1 , fhr_mate_pairs_next_to_each_other) ) != -1 )
	{
		if ( line1[0] == '@' )
			fprintf (fhw , "%s" , line1);
		else break;
	}
	while ( ( line1_len = getline ( &line1 , &len1 , fhr_mate_pairs_away_from_each_other) ) != -1 )
		if ( line1[0] != '@' ) break;

	read_number = 0;
	read_length = 0;
	do
	{
		line1_len = getline ( &line1 , &len1 , fhr_mate_pairs_next_to_each_other);
		if ( line1_len == -1 ) break;
		line2_len = getline ( &line2 , &len2 , fhr_mate_pairs_next_to_each_other);
		if ( line2_len == -1 ) break;

		splitMappingInTwoPartsAndSetNHValue (line1 , split_line1 , &NH_value1);
		splitMappingInTwoPartsAndSetNHValue (line2 , split_line2 , &NH_value2);

		generateNextReadID (alphabets , read_id , &read_length);
		convertReadIdToString (read_id , read_id_string , read_length , alphabets);
		strcpy(split_line1[0] , read_id_string);
		strcpy(split_line2[0] , read_id_string);

		writeToFile (split_line1 , fhw);
		writeToFile (split_line2 , fhw);
	} while ( 1 );

	do
	{
		line1_len = getline ( &line1 , &len1 , fhr_mate_pairs_away_from_each_other);
		if ( line1_len == -1 ) break;
		line2_len = getline ( &line2 , &len2 , fhr_mate_pairs_away_from_each_other);
		if ( line2_len == -1 ) break;

		splitMappingInTwoPartsAndSetNHValue (line1 , split_line1 , &NH_value1);
		splitMappingInTwoPartsAndSetNHValue (line2 , split_line2 , &NH_value2);

		generateNextReadID (alphabets , read_id , &read_length);
		convertReadIdToString (read_id , read_id_string , read_length , alphabets);
		strcpy(split_line1[0] , read_id_string);
		strcpy(split_line2[0] , read_id_string);

		writeToFile (split_line1 , fhw);
		writeToFile (split_line2 , fhw);
	} while ( 1 );

	fclose (fhr_mate_pairs_next_to_each_other);
	fclose (fhr_mate_pairs_away_from_each_other);
	fclose (fhw);
}

int main (int argc, char *argv[])
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	int read_length;
	int read_id[100];
	int i, j;

	char input_samfilename_mate_pairs_next_to_each_other[FILENAME_LENGTH];
	char input_samfilename_mate_pairs_away_from_each_other[FILENAME_LENGTH];
	char output_samfilename[FILENAME_LENGTH];

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	read_length = 0;
	strcpy(input_samfilename_mate_pairs_next_to_each_other , argv[1]);
	strcpy(input_samfilename_mate_pairs_away_from_each_other , argv[2]);
	strcpy(output_samfilename , argv[3]);

	/********************************************************************/
	convertOldReadIdsToNewReadIds (input_samfilename_mate_pairs_next_to_each_other , input_samfilename_mate_pairs_away_from_each_other , output_samfilename);
	return 0;
}
