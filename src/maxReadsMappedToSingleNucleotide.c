# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

void findMaximumNumberOfReadsMappedToOneNucleotide (char *input_samfilename, char *output_filename, char *name_of_total_number_of_alignments_file)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	int i, j, k;
	int number_of_fields;

	long long int max_position, max_value;
	long long int prev_position, prev_value;
	long long int curr_position, curr_value;
	long long int total_number_of_alignments;

	FILE *fhr;
	FILE *fhw;
	FILE *fhw_tot_alignments;

	size_t len = 0;
	ssize_t line_len;

	char str[100];
	char *temp; //Useless
	char *line = NULL; // for reading each line
	char **split_line; // List of strings to store each element of a single alignment
	char **split_tags;

	//struct Sam_Alignment *curr_alignment;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	//curr_alignment = allocateMemorySam_Alignment ();
	fhr = fopen (input_samfilename , "r");
	if ( fhr == NULL )
	{
		printf ("Error! File %s not found" , input_samfilename);
		exit (1);
	}
	fhw = fopen (output_filename , "w");
	if ( fhw == NULL )
	{
		printf ("%s File cannot be created" , output_filename);
		exit (1);
	}
	fhw_tot_alignments = fopen (name_of_total_number_of_alignments_file , "w");
	if ( fhw_tot_alignments == NULL )
	{
		printf ("%s File cannot be created" , name_of_total_number_of_alignments_file);
		exit (1);
	}

	max_position = 0;
	max_value = 0;
	curr_position = 0;
	curr_value = 0;
	prev_position = 0;
	prev_value = 0;

	split_line = ( char** ) malloc (sizeof(char*) * ROWS);
	for ( i = 0 ; i < ROWS ; i++ )
		split_line[i] = ( char* ) malloc (sizeof(char) * COLS);

	split_tags = ( char** ) malloc (sizeof(char*) * ROWS);
	for ( i = 0 ; i < ROWS ; i++ )
		split_tags[i] = ( char* ) malloc (sizeof(char) * COLS);
	/********************************************************************/

	while ( ( line_len = getline ( &line , &len , fhr) ) != -1 )
		if ( line[0] != '@' ) break;

	total_number_of_alignments = 0;
	do
	{
		total_number_of_alignments += 1;
		number_of_fields = splitByDelimiter (line , '\t' , split_line);
		//populateSamAlignmentInstance ( curr_alignment , split_line , number_of_fields , split_tags );

		curr_position = strtol (split_line[3] , &temp , 10);
		if ( curr_position == 0 ) continue;
		if ( max_position == 0 )
		{
			max_position = curr_position;
			prev_position = curr_position;

			max_value = 1;
			curr_value = 1;
			prev_value = 1;
		}
		else
		{
			if ( prev_position == curr_position )
				prev_value += 1;
			else
			{
				if ( prev_value > max_value )
				{
					max_value = prev_value;
					max_position = prev_position;
				}
				prev_position = curr_position;
				prev_value = 1;
			}
		}

	} while ( ( line_len = getline ( &line , &len , fhr) ) != -1 );
	if ( prev_value > max_value )
	{
		max_value = prev_value;
		max_position = prev_position;
		prev_position = curr_position;
		prev_value = 1;
	}

	sprintf(str , "%lld" , max_value);
	strcat(str , "\n");
	fprintf (fhw , "%s" , str);

	sprintf(str , "%lld" , total_number_of_alignments);
	strcat(str , "\n");
	fprintf (fhw_tot_alignments , "%s" , str);

	fclose (fhw);
	fclose (fhw_tot_alignments);
	fclose (fhr);
}

int main (int argc, char *argv[])
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char input_samfilename[FILENAME_LENGTH];
	char output_filename[FILENAME_LENGTH];
	char name_of_total_number_of_alignments_file[FILENAME_LENGTH];

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/

	strcpy(input_samfilename , argv[1]);
	strcpy(output_filename , argv[2]);
	strcpy(name_of_total_number_of_alignments_file , argv[3]);

	/********************************************************************/
	findMaximumNumberOfReadsMappedToOneNucleotide (input_samfilename , output_filename , name_of_total_number_of_alignments_file);

}
