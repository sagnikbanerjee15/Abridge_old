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

void writeToFile (struct Sam_Alignment *curr_alignment, FILE *fhw)
{
	char str[1000];
	int i;

	fprintf (fhw , "%s" , curr_alignment->read_name);
	fprintf (fhw , "%s" , "\t");

	sprintf(str , "%ld" , curr_alignment->samflag);
	fprintf (fhw , "%s" , str);
	fprintf (fhw , "%s" , "\t");

	fprintf (fhw , "%s" , curr_alignment->reference_name);
	fprintf (fhw , "%s" , "\t");

	sprintf(str , "%ld" , curr_alignment->start_position);
	fprintf (fhw , "%s" , str);
	fprintf (fhw , "%s" , "\t");

	sprintf(str , "%ld" , curr_alignment->mapping_quality_score);
	fprintf (fhw , "%s" , str);
	fprintf (fhw , "%s" , "\t");

	fprintf (fhw , "%s" , curr_alignment->cigar);
	fprintf (fhw , "%s" , "\t");

	fprintf (fhw , "%s" , curr_alignment->reference_name_next_mate);
	fprintf (fhw , "%s" , "\t");

	sprintf(str , "%ld" , curr_alignment->start_position_next);
	fprintf (fhw , "%s" , str);
	fprintf (fhw , "%s" , "\t");

	sprintf(str , "%ld" , curr_alignment->template_length);
	fprintf (fhw , "%s" , str);
	fprintf (fhw , "%s" , "\t");

	fprintf (fhw , "%s" , curr_alignment->seq);
	fprintf (fhw , "%s" , "\t");

	for ( i = 0 ; curr_alignment->qual[i] != '\0' ; i++ )
		curr_alignment->qual[i] -= 90;
	fprintf (fhw , "%s" , curr_alignment->qual);
	fprintf (fhw , "%s" , "\t");

	for ( i = 0 ; i < curr_alignment->number_of_tag_items ; i++ )
	{
		fprintf (fhw , "%s" , curr_alignment->tags[i].name);
		fprintf (fhw , "%s" , ":");
		fprintf (fhw , "%s" , curr_alignment->tags[i].type);
		fprintf (fhw , "%s" , ":");
		fprintf (fhw , "%s" , curr_alignment->tags[i].val);
		fprintf (fhw , "%s" , "\t");
	}
}

void convertOldReadIdsToNewReadIds (char *input_samfilename, char *output_samfilename)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char alphabets[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-='{}[]|?<>,.";
	char *line = NULL; // for reading each line
	char **split_line; // List of strings to store each element of a single alignment
	char **split_tags;
	char read_id_string[100];

	int read_length;
	int i;
	int j;
	int num_elements_read_id_mapping_dictionary;
	int number_of_tags;
	int number_of_fields; // Number of fields in each sam alignment entry
	int NH_tag_index;
	int read_id[100];

	size_t len = 0;
	ssize_t line_len;

	FILE *fhr;
	FILE *fhw;

	struct Old_Read_ID_to_New_Read_ID **read_id_mapping;
	struct Sam_Alignment *curr_alignment;
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
	fhw = fopen (output_samfilename , "w");
	if ( fhw == NULL )
	{
		printf ("Error! File %s cannot be created" , output_samfilename);
		exit (1);
	}

	split_line = ( char** ) malloc (sizeof(char*) * ROWS);
	for ( i = 0 ; i < ROWS ; i++ )
		split_line[i] = ( char* ) malloc (sizeof(char) * COLS);

	split_tags = ( char** ) malloc (sizeof(char*) * ROWS);
	for ( i = 0 ; i < ROWS ; i++ )
		split_tags[i] = ( char* ) malloc (sizeof(char) * COLS);

	read_id_mapping = ( struct Old_Read_ID_to_New_Read_ID** ) malloc (sizeof(struct Old_Read_ID_to_New_Read_ID*) * ROWS);
	for ( i = 0 ; i < ROWS ; i++ )
		read_id_mapping[i] = allocateMemoryOld_Read_ID_to_New_Read_ID ();
	read_length = 0;
	num_elements_read_id_mapping_dictionary = 0;

	curr_alignment = allocateMemorySam_Alignment ();
	/********************************************************************/

	while ( ( line_len = getline ( &line , &len , fhr) ) != -1 )
	{
		if ( line[0] == '@' )
			fprintf (fhw , "%s" , line);
		else break;
	}

	do
	{
		number_of_fields = splitByDelimiter (line , '\t' , split_line);
		populateSamAlignmentInstance (curr_alignment , split_line , number_of_fields , split_tags);

		NH_tag_index = -1;
		for ( i = 0 ; i < curr_alignment->number_of_tag_items ; i++ )
			if ( strcmp (curr_alignment->tags[i].name , "NH") == 0 )
				NH_tag_index = i;

		if ( strcmp (curr_alignment->tags[NH_tag_index].val , "1") == 0 ) //Only a single occurance
		{
			generateNextReadID (alphabets , read_id , &read_length);
			convertReadIdToString (read_id , read_id_string , read_length , alphabets);
			strcpy(curr_alignment->read_name , read_id_string);
			writeToFile (curr_alignment , fhw);
		}
		else
		{

		}

	} while ( ( line_len = getline ( &line , &len , fhr) ) != -1 );
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
	char input_samfilename[FILENAME_LENGTH];
	char output_samfilename[FILENAME_LENGTH];

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	read_length = 0;
	strcpy(input_samfilename , argv[1]);
	strcpy(output_samfilename , argv[2]);

	/********************************************************************/
	convertOldReadIdsToNewReadIds (input_samfilename , output_samfilename);

	return 0;
	for ( i = 0 ; i < 15000000 ; i++ )
	{
		generateNextReadID (alphabets , read_id , &read_length);
		for ( j = 0 ; j < read_length ; j++ )
			printf ("%c" , alphabets[read_id[j]]);
		printf ("\n");
	}

	return 0;
}
