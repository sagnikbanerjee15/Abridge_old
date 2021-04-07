# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include <math.h>
# include <stdbool.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

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

void writeToFile (char **split_line, FILE *fhw)
{
	fprintf (fhw , "%s" , split_line[0]);
	fprintf (fhw , "%s" , "\t");
	fprintf (fhw , "%s" , split_line[1]);
	fprintf (fhw , "%s" , "\n");
}

void insertNewEntryInMappingDictionary (char *new_read_name, char *old_read_name, int max_elements, struct Old_Read_ID_to_New_Read_ID **read_id_mapping, int NH_val)
{
	int i;
	for ( i = 0 ; i < max_elements ; i++ )
	{
		if ( read_id_mapping[i]->valid == 0 )
		{
			//printf ("\nAdding %s to mapping dictionary %s-->%s" , old_read_name , old_read_name , new_read_name);
			//fflush (stdout);
			strcpy(read_id_mapping[i]->new_read_id , new_read_name);
			strcpy(read_id_mapping[i]->old_read_id , old_read_name);
			read_id_mapping[i]->number_of_multi_maps = NH_val;
			read_id_mapping[i]->valid = 1;
			break;
		}
	}
}

int searchOldReadNameInMappingDictionary (char *old_read_name, int max_elements, struct Old_Read_ID_to_New_Read_ID **read_id_mapping)
{
	int i;
	for ( i = max_elements - 1 ; i >= 0 ; i-- )
		if ( read_id_mapping[i]->valid == 1 && strcmp (read_id_mapping[i]->old_read_id , old_read_name) == 0 )
			return i;

	return -1;
}

void splitMappingInTwoParts (char *line, char **split_line)
{
	int i, j0, j1, first_tab_found = 0, k;

	j0 = 0;
	j1 = 0;
	for ( i = 0 ; line[i] != '\0' ; i++ )
	{
		if ( first_tab_found )
			split_line[1][j1++ ] = line[i];
		else split_line[0][j0++ ] = line[i];
		if ( line[i] == '\t' && first_tab_found == 0 ) first_tab_found = 1;
	}
	split_line[0][j0] = '\0';
	split_line[1][j1] = '\0';
}

int findNumberOfValidMappings (struct Old_Read_ID_to_New_Read_ID **read_id_mapping, int num_elements_read_id_mapping_dictionary)
{
	int total = 0;
	int i;
	for ( i = 0 ; i < num_elements_read_id_mapping_dictionary ; i++ )
		if ( read_id_mapping[i]->valid == 1 ) total += 1;

	return total;
}

void convertOldReadIdsToNewReadIds (char *input_samfilename_um, char *input_samfilename_mm, char *output_samfilename)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char alphabets[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-='{}[]|?<>,.";
	char *temp; //Required for strtoi
	char *line = NULL; // for reading each line
	char **split_line; // List of strings to store each element of a single alignment
	char **split_tags;
	char prev_old_read[1000];
	char read_id_string[100];

	int read_length;
	int i;
	int j;
	int read_id[100];

	long long int read_number;

	size_t len = 0;
	ssize_t line_len;

	FILE *fhr_um;
	FILE *fhr_mm;
	FILE *fhw;

	struct Old_Read_ID_to_New_Read_ID **read_id_mapping;
	struct Sam_Alignment *curr_alignment;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr_um = fopen (input_samfilename_um , "r");
	if ( fhr_um == NULL )
	{
		printf ("Error! File %s not found" , input_samfilename_um);
		exit (1);
	}
	fhr_mm = fopen (input_samfilename_mm , "r");
	if ( fhr_mm == NULL )
	{
		printf ("Error! File %s not found" , input_samfilename_mm);
		exit (1);
	}
	fhw = fopen (output_samfilename , "w");
	if ( fhw == NULL )
	{
		printf ("Error! File %s cannot be created" , output_samfilename);
		exit (1);
	}

	split_line = ( char** ) malloc (sizeof(char*) * 2);
	for ( i = 0 ; i < 2 ; i++ )
		split_line[i] = ( char* ) malloc (sizeof(char) * 1000);

	/********************************************************************/

	while ( ( line_len = getline ( &line , &len , fhr_um) ) != -1 )
	{
		if ( line[0] == '@' )
			fprintf (fhw , "%s" , line);
		else break;
	}

	read_number = 0;

	do
	{
		read_number++;
		splitMappingInTwoParts (line , split_line);
		if ( read_number % 2 == 1 )
		{
			generateNextReadID (alphabets , read_id , &read_length);
			convertReadIdToString (read_id , read_id_string , read_length , alphabets);
			for ( i = 0 ; i < read_length ; i++ )
				printf ("\n%d " , read_id[i]);
			printf ("\nread_id_string = %s read_length = %d" , read_id_string , read_length);
		}
		strcpy(split_line[0] , read_id_string);
		writeToFile (split_line , fhw);

	} while ( ( line_len = getline ( &line , &len , fhr_um) ) != -1 );

	while ( ( line_len = getline ( &line , &len , fhr_mm) ) != -1 )
		if ( line[0] != '@' ) break;

	do
	{
		splitMappingInTwoParts (line , split_line);
		if ( strcmp (prev_old_read , split_line[0]) != 0 )
		{
			generateNextReadID (alphabets , read_id , &read_length);
			convertReadIdToString (read_id , read_id_string , read_length , alphabets);
			strcpy(prev_old_read , split_line[0]);
		}
		strcpy(split_line[0] , read_id_string);
		writeToFile (split_line , fhw);
	} while ( ( line_len = getline ( &line , &len , fhr_mm) ) != -1 );

	fclose (fhr_um);
	fclose (fhr_mm);
}

int main (int argc, char *argv[])
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	int read_length;
	int read_id[100];
	int i, j;

	char alphabets[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-='{}[]|?<>,.";
	char input_samfilename_um[FILENAME_LENGTH];
	char input_samfilename_mm[FILENAME_LENGTH];
	char output_samfilename[FILENAME_LENGTH];

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	read_length = 0;
	strcpy(input_samfilename_um , argv[1]);
	strcpy(input_samfilename_mm , argv[2]);
	strcpy(output_samfilename , argv[3]);

	/********************************************************************/
	convertOldReadIdsToNewReadIds (input_samfilename_um , input_samfilename_mm , output_samfilename);
	return 0;
}
