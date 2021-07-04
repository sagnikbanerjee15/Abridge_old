# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

void findMaximumNumberOfReadsInEachLine (char *pass1filename, char *output_filename)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	int i, j, k;
	int number_of_fields;
	int max_read_length;

	FILE *fhr;
	FILE *fhw;

	size_t len = 0;
	ssize_t line_len;

	unsigned long long int max_reads_in_each_line;
	unsigned long long int number_of_commas_in_each_line;
	unsigned long long int max_number_of_commas = 0;

	char *line = NULL; // for reading each line
	char str[100];
	char *temp; //Useless
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr = fopen (pass1filename , "r");
	if ( fhr == NULL )
	{
		printf ("Error! File %s not found" , pass1filename);
		exit (1);
	}
	fhw = fopen (output_filename , "w");
	if ( fhw == NULL )
	{
		printf ("%s File cannot be created" , output_filename);
		exit (1);
	}
	/********************************************************************/

	while ( ( line_len = getline ( &line , &len , fhr) ) != -1 )
		if ( line[0] != '@' ) break;

	max_reads_in_each_line = 0;
	do
	{
		number_of_commas_in_each_line = 1;
		for ( i = 0 ; line[i] != '\0' ; i++ )
			if ( line[i] == ',' ) number_of_commas_in_each_line++;
		if ( max_number_of_commas < number_of_commas_in_each_line )
			max_number_of_commas = number_of_commas_in_each_line;

	} while ( ( line_len = getline ( &line , &len , fhr) ) != -1 );

	sprintf(str , "%lld" , max_number_of_commas);
	strcat(str , "\n");
	fprintf (fhw , "%s" , str);

	fclose (fhr);
	fclose (fhw);
}

int main (int argc, char *argv[])
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char pass1filename[FILENAME_LENGTH];
	char output_filename[FILENAME_LENGTH];

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/

	strcpy(pass1filename , argv[1]);
	strcpy(output_filename , argv[2]);

	/********************************************************************/
	findMaximumNumberOfReadsInEachLine (pass1filename , output_filename);

}
