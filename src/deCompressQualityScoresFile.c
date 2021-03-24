# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

void convertRLEtoQualValues ( char *input_qualityscore_filename, char *output_quality_score_filename )
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	FILE *fhr;
	FILE *fhw;

	size_t len = 0;
	ssize_t line_len;

	int i;
	int max_read_length;
	int rle_quality_scores_index = 0;
	int number_of_quality_scores_in_current_position = 0;

	char *line;

	struct RLE_Quality_Scores **rle_quality_scores;

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr = fopen ( input_qualityscore_filename , "r" );
	if ( fhr == NULL )
	{
		printf ( "Error! File %s not found" , input_qualityscore_filename );
		exit ( 1 );
	}
	fhw = fopen ( output_quality_score_filename , "w" );
	if ( fhw == NULL )
	{
		printf ( "%s File cannot be created" , output_quality_score_filename );
		exit ( 1 );
	}

	rle_quality_scores = ( struct RLE_Quality_Scores** ) malloc ( sizeof(struct RLE_Quality_Scores*) * MAX_SEQ_LEN );

	/********************************************************************/

	while ( ( line_len = getline ( &line , &len , fhr ) ) != -1 )
	{
		for ( i = 0 ; line[i] != '\0' ; i++ )
		{

		}
	}
}

int main ( int argc, char *argv[] )
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char *temp;
	char input_qualityscore_filename[FILENAME_LENGTH];
	char output_quality_score_filename[FILENAME_LENGTH];

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	strcpy( input_qualityscore_filename , argv[1] );
	strcpy( output_quality_score_filename , argv[2] );
	/********************************************************************/

	convertRLEtoQualValues ( input_qualityscore_filename , output_quality_score_filename );
}
