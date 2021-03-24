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
	int num;
	int max_read_length;
	int rle_quality_scores_index = 0;
	int number_of_quality_scores_in_current_position = 0;
	int number_of_quality_scores_in_current_position_index = 0;
	int checker_flag;
	int *quality_score_position_index;

	char *line;
	char *quality_score_of_read;

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

	max_read_length = 0;
	while ( ( line_len = getline ( &line , &len , fhr ) ) != -1 )
	{
		number_of_quality_scores_in_current_position = 0;
		for ( i = 0 ; line[i] != '\0' ; i++ )
		{
			if ( isdigit ( line[i] ) == 0 )
				number_of_quality_scores_in_current_position++;
		}
		rle_quality_scores[rle_quality_scores_index] = ( struct RLE_Quality_Scores* ) malloc ( sizeof(struct RLE_Quality_Scores) * number_of_quality_scores_in_current_position );
		number_of_quality_scores_in_current_position_index = 0;
		for ( i = 0 ; line[i] != '\0' ; i++ )
		{
			/*
			 * Check if the first element is a character or not
			 */
			if ( isdigit ( line[0] ) == 0 )
			{
				rle_quality_scores[rle_quality_scores_index][number_of_quality_scores_in_current_position_index].quality_score = line[0] - 30;
				rle_quality_scores[rle_quality_scores_index][number_of_quality_scores_in_current_position_index].frequency = 1;
				number_of_quality_scores_in_current_position_index++;
			}
			else
			{
				num = 0;
				while ( isidigit ( line[i] ) != 0 )
				{
					num = num * 10 + line[i] - 48;
					i++;
				}
				rle_quality_scores[rle_quality_scores_index][number_of_quality_scores_in_current_position_index].quality_score = line[i] - 30;
				rle_quality_scores[rle_quality_scores_index][number_of_quality_scores_in_current_position_index].frequency = num;
				number_of_quality_scores_in_current_position_index++;
			}
		}
		rle_quality_scores_index++;
	}
	max_read_length = rle_quality_scores_index;
	quality_score_position_index = ( int* ) malloc ( sizeof(int) * max_read_length );
	quality_score_of_read = ( char* ) malloc ( sizeof(char) * max_read_length );

	for ( i = 0 ; i < max_read_length ; i++ )
		quality_score_position_index[i] = 0;

	/*
	 * Start constructing the quality scores
	 */
	while ( 1 )
	{
		checker_flag = 0;
		for ( i = 0 ; i < max_read_length ; i++ )
		{
			quality_score_of_read[i] = rle_quality_scores[i][quality_score_position_index[i]].quality_score;
			rle_quality_scores[i][quality_score_position_index[i]].frequency--;
			if ( rle_quality_scores[i][quality_score_position_index[i]].frequency == 0 )
				quality_score_position_index[i]++;
			checker_flag += rle_quality_scores[i][quality_score_position_index[i]].frequency;
		}
		quality_score_of_read[i] = '\0';
		printf ( "\n%s" , quality_score_of_read );
		if ( checker_flag == 0 ) break;
	}

	fclose ( fhr );
	fclose ( fhw );
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
