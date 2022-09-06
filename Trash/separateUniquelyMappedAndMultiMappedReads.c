# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

int isMappingUnique (char *line)
{
	int unique = 0;
	int i;

	for ( i = 0 ; line[i + 6] != '\0' ; i++ )
		if ( line[i] == 'N' && line[i + 1] == 'H' && line[i + 2] == ':' && line[i + 3] == 'i' && line[i + 4] == ':' && ( line[i + 5] == '1' || line[i + 5] == '0' ) && line[i + 6] == '\t' )
			unique = 1;

	return unique;
}

int main (int argc, char *argv[])
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char uniquely_mapped_reads_filename[FILENAME_LENGTH];
	char multi_mapped_reads_filename[FILENAME_LENGTH];
	char input_samfilename[FILENAME_LENGTH];
	char *line;

	size_t len = 0;
	ssize_t line_len;

	FILE *fhr;
	FILE *fhw_unique;
	FILE *fhw_multi;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	strcpy(input_samfilename , argv[1]);
	strcpy(uniquely_mapped_reads_filename , argv[2]);
	strcpy(multi_mapped_reads_filename , argv[3]);

	fhr = fopen (input_samfilename , "r");
	if ( fhr == NULL )
	{
		printf ("Error! File %s not found" , input_samfilename);
		exit (1);
	}
	fhw_unique = fopen (uniquely_mapped_reads_filename , "w");
	if ( fhw_unique == NULL )
	{
		printf ("%s File cannot be created" , uniquely_mapped_reads_filename);
		exit (1);
	}
	fhw_multi = fopen (multi_mapped_reads_filename , "w");
	if ( fhw_multi == NULL )
	{
		printf ("%s File cannot be created" , multi_mapped_reads_filename);
		exit (1);
	}
	/********************************************************************/

	while ( ( line_len = getline ( &line , &len , fhr) ) != -1 )
	{
		if ( line[0] == '@' )
		{
			fprintf (fhw_unique , "%s" , line);
			fprintf (fhw_multi , "%s" , line);
		}
		else break;
	}

	do
	{
		if ( isMappingUnique (line) )
			fprintf (fhw_unique , "%s" , line);
		else fprintf (fhw_multi , "%s" , line);
	} while ( ( line_len = getline ( &line , &len , fhr) ) != -1 );
	return 0;
}
