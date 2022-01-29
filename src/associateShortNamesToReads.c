# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

void generateNextShortReadId (char *short_read_id)
{
	int number_of_zeros = 0;
	if ( short_read_id[0] == '\0' ) //Empty string
	{
		short_read_id[0] = '0';
		short_read_id[1] = '\0';
	}
	else if ( strlen (short_read_id) == 1 && short_read_id[0] == 122 )
	{
		short_read_id[0] = 48;
		short_read_id[1] = 48;
		short_read_id[2] = '\0';
	}
	else
	{

		int i = strlen (short_read_id) - 1;
		short_read_id[i]++;
		i--;

		for ( ; i >= 0 && short_read_id[i + 1] == 123 ; i-- )
			short_read_id[i]++;

		for ( i = 0 ; i < strlen (short_read_id) ; i++ )
		{
			if ( short_read_id[i] == 58 )
				short_read_id[i] = 65;
			else if ( short_read_id[i] == 91 )
				short_read_id[i] = 97;
			else if ( short_read_id[i] == 123 )
			{
				short_read_id[i] = 48;
				number_of_zeros++;
			}
		}
		if ( number_of_zeros == strlen (short_read_id) )
		{
			short_read_id[strlen (short_read_id)] = '\0';
			short_read_id[strlen (short_read_id) - 1] = 48;
		}
	}
}

void assignShortenedReadsNames (char *inputfilename, char *outputfilename)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	FILE *fhr;
	FILE *fhw;

	int i;

	size_t len = 0;
	ssize_t line_len;

	char *temp; //Useless
	char *line = NULL; // for reading each line
	char *current_read_id;
	char *previous_read_id;

	char **split_on_tab;
	char *short_read_id;

	unsigned long long int line_number = 0;
	unsigned int number_of_characters_in_short_read_name = 1;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr = fopen (inputfilename , "r");
	if ( fhr == NULL )
	{
		printf ("\nCannot create file %s" , inputfilename);
		exit (0);
	}

	fhw = fopen (outputfilename , "w");
	if ( fhw == NULL )
	{
		printf ("\nCannot create file %s" , outputfilename);
		exit (0);
	}

	split_on_tab = ( char** ) malloc (sizeof(char*) * 5);
	for ( i = 0 ; i < 5 ; i++ )
		split_on_tab[i] = ( char* ) malloc (sizeof(char) * 1000);

	current_read_id = ( char* ) malloc (sizeof(char) * 1000);
	previous_read_id = ( char* ) malloc (sizeof(char) * 1000);
	short_read_id = ( char* ) malloc (sizeof(char) * 100);
	short_read_id[0] = '\0';
	current_read_id[0] = '\0';
	previous_read_id[0] = '\0';

	/********************************************************************/

	while ( ( line_len = getline ( &line , &len , fhr) ) != -1 )
	{
		splitByDelimiter (line , '\t' , split_on_tab);
		strcpy(current_read_id , split_on_tab[0]);
		if ( strcmp (previous_read_id , current_read_id) != 0 ) // Same read name
			generateNextShortReadId (short_read_id);

		fprintf (fhw , "%s" , split_on_tab[0]);
		fprintf (fhw , "%s" , "\t");
		fprintf (fhw , "%s" , split_on_tab[1]);
		fprintf (fhw , "%s" , "\t");
		fprintf (fhw , "%s" , split_on_tab[2]);
		fprintf (fhw , "%s" , "\t");
		fprintf (fhw , "%s" , short_read_id);
		fprintf (fhw , "%s" , "\n");

		strcpy(previous_read_id , current_read_id);
	}

	fclose (fhr);
	fclose (fhw);
}

int main (int argc, char *argv[])
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char read_name_to_line_numbers[FILENAME_LENGTH];
	char read_name_to_line_numbers_to_shortened_read_names[FILENAME_LENGTH];

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	strcpy(read_name_to_line_numbers , argv[1]);
	strcpy(read_name_to_line_numbers_to_shortened_read_names , argv[2]);

	/********************************************************************/

	assignShortenedReadsNames (read_name_to_line_numbers ,
			read_name_to_line_numbers_to_shortened_read_names);

	return 0;
}
