/*
 * Add NH tags to SAM files that do not have those
 * Input1: Inputsamfilename
 * Input2: Outputsamfilename
 * Input3: read_multi_mapping_filename
 * Input4: SE or PE
 */

# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

void findMultiMappedReadsFromSamFile (
		char *input_samfilename,
		char *output_samfilename)

{
	FILE *fhr;
	FILE *fhw;
	char temp[100];
	char faidx_filename[1000];
	char *buffer = NULL;
	char **split_on_tab;
	char line_to_be_written_to_file[10000];
	char **read_ids;
	char **multi_mapped_read_ids;
	unsigned long long int number_of_lines_in_samfile = 0;
	unsigned long long int read_ids_i, multi_mapped_read_ids_i;
	unsigned long long int i, j, k;

	size_t len = 0;
	ssize_t line_len;

	fhr = fopen (input_samfilename , "r");
	if ( fhr == NULL )
	{
		printf ("Error! File not found %s" , input_samfilename);
		exit (1);
	}

	split_on_tab = ( char** ) malloc (sizeof(char*) * 50);
	for ( i = 0 ; i < 50 ; i++ )
		split_on_tab[i] = ( char* ) malloc (sizeof(char) * COLS);

	while ( ( line_len = getline ( &buffer , &len , fhr) ) != -1 )
		number_of_lines_in_samfile += 1;

	read_ids = ( char** ) malloc (sizeof(char*) * number_of_lines_in_samfile);
	for ( k = 0 ; k < number_of_lines_in_samfile ; k++ )
		read_ids[k] = ( char* ) malloc (sizeof(char) * 50);

	multi_mapped_read_ids = ( char** ) malloc (sizeof(char*) * number_of_lines_in_samfile);
	for ( k = 0 ; k < number_of_lines_in_samfile ; k++ )
		multi_mapped_read_ids[k] = ( char* ) malloc (sizeof(char) * 50);

	rewind (fhr);
	while ( ( line_len = getline ( &buffer , &len , fhr) ) != -1 )
	{
		if ( buffer[0] == '@' )
			continue;
		else break;
	}

	read_ids_i = 0;
	do
	{
		splitByDelimiter (buffer , '\t' , split_on_tab);
		strcpy(read_ids[read_ids_i++ ] , split_on_tab[0]);
	} while ( ( line_len = getline ( &buffer , &len , fhr) ) != -1 );

	number_of_lines_in_samfile = read_ids_i;
	printf ("\nTotal number of reads %llu, Number of multi-mapped reads %llu" ,
			number_of_lines_in_samfile ,
			multi_mapped_read_ids_i);
	/*
	 * Find duplicated reads
	 */

	multi_mapped_read_ids_i = 0;
	for ( i = 0 ; i < number_of_lines_in_samfile ; i++ )
	{
		//if ( i % 10000 == 0 )
		printf ("\n%llu" , i);
		for ( j = i + 1 ; j < number_of_lines_in_samfile ; j++ )
		{
			if ( read_ids[i][0] == read_ids[j][0] )
			{
				if ( !strcmp (read_ids[i] , read_ids[j]) )
				{
					strcpy(multi_mapped_read_ids[multi_mapped_read_ids_i++ ] ,
							read_ids[i]);
					strcpy(read_ids[j] , "");
					break;
				}
			}
		}
		printf ("% llu" , j);
	}
	printf ("\nTotal number of reads %llu, Number of multi-mapped reads %llu" ,
			number_of_lines_in_samfile ,
			multi_mapped_read_ids_i);
}

int main (int argc, char *argv[])
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char input_samfilename[FILENAME_LENGTH];
	char output_samfilename[FILENAME_LENGTH];
	char read_multi_mapping_filename[FILENAME_LENGTH];
	char *temp;
	struct Read_Ids_to_NH *read_ids_to_nh;

	int paired;
	int number_of_lines_in_read_multi_mapping_filename;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	strcpy(input_samfilename , argv[1]);
	strcpy(output_samfilename , argv[2]);
	strcpy(read_multi_mapping_filename , argv[3]);
	paired = strtol (argv[4] , &temp , 10);
	/********************************************************************/

	/*
	 * Count the number of lines in read_multi_mapping_filename
	 */
//number_of_lines_in_read_multi_mapping_filename = countNumberOfLines (read_multi_mapping_filename);
	/*
	 * Create array of structures
	 */
	/*
	 read_ids_to_nh = ( struct read_ids_to_nh* ) malloc (sizeof(struct read_ids_to_nh) * number_of_lines_in_read_multi_mapping_filename);
	 readMultiMappingInformationInDS (read_ids_to_nh ,
	 number_of_lines_in_read_multi_mapping_filename ,
	 read_multi_mapping_filename);
	 */

	findMultiMappedReadsFromSamFile (input_samfilename , output_samfilename);

	return 0;
}
