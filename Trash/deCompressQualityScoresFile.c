# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

void convertRLEtoQualValues (
		char *input_qualityscore_filename,
		char *output_quality_score_filename)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	FILE *fhr;
	FILE *fhw;

	size_t len = 0;
	ssize_t line_len;

	short int line_number = 0;

	int i, j, k;
	int num;
	int max_read_length;
	int rle_quality_scores_index = 0;
	int number_of_quality_scores_in_current_position = 0;
	int number_of_quality_scores_in_current_position_index = 0;
	int checker_flag;
	int *quality_score_position_index;
	int *quality_score_position_max;

	char *line;
	char *quality_score_of_read;

	struct RLE_Quality_Scores **rle_quality_scores;

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr = fopen (input_qualityscore_filename , "r");
	if ( fhr == NULL )
	{
		printf ("Error! File %s not found" , input_qualityscore_filename);
		exit (1);
	}
	fhw = fopen (output_quality_score_filename , "w");
	if ( fhw == NULL )
	{
		printf ("%s File cannot be created" , output_quality_score_filename);
		exit (1);
	}

	rle_quality_scores = ( struct RLE_Quality_Scores** ) malloc (sizeof(struct RLE_Quality_Scores*) * MAX_SEQ_LEN);
	quality_score_position_max = ( int* ) malloc (sizeof(int) * MAX_SEQ_LEN);
	quality_score_position_index = ( int* ) malloc (sizeof(int) * MAX_SEQ_LEN);
	for ( i = 0 ; i < MAX_SEQ_LEN ; i++ )
	{
		quality_score_position_index[i] = 0;
		quality_score_position_max[i] = 0;
	}
	/********************************************************************/

	max_read_length = 0;
	while ( ( line_len = getline ( &line , &len , fhr) ) != -1 )
	{
		line_number++;
		//printf ("\nProcessing line %d" , line_number);
		number_of_quality_scores_in_current_position = 0;
		for ( i = 0 ; line[i] != '\0' ; i++ )
		{
			if ( isdigit (line[i]) == 0 )
				number_of_quality_scores_in_current_position++;
		}
		/*
		 printf ("\n%d %d" , rle_quality_scores_index , number_of_quality_scores_in_current_position);
		 fflush ( stdout);
		 */
		rle_quality_scores[rle_quality_scores_index] = ( struct RLE_Quality_Scores* ) malloc (sizeof(struct RLE_Quality_Scores) * ( number_of_quality_scores_in_current_position + 10 ));

		number_of_quality_scores_in_current_position_index = 0;
		for ( i = 0 ; line[i] != '\0' ; i++ )
		{
			/*
			 * Check if the first element is a character or not
			 */
			if ( line[i] >= 48 && line[i] <= 57 )
			{
				num = 0;
				while ( line[i] >= 48 && line[i] <= 57 )
				{
					num = num * 10 + line[i] - 48;
					i++;
				}
				rle_quality_scores[rle_quality_scores_index][number_of_quality_scores_in_current_position_index].quality_score = line[i] - 26;
				rle_quality_scores[rle_quality_scores_index][number_of_quality_scores_in_current_position_index].frequency = num;
				number_of_quality_scores_in_current_position_index++;
				//quality_score_position_max[i]++;

			}
			else
			{
				rle_quality_scores[rle_quality_scores_index][number_of_quality_scores_in_current_position_index].quality_score = line[i] - 26;
				rle_quality_scores[rle_quality_scores_index][number_of_quality_scores_in_current_position_index].frequency = 1;
				number_of_quality_scores_in_current_position_index++;
				//quality_score_position_max[i]++;
			}
		}

		quality_score_position_max[rle_quality_scores_index] = number_of_quality_scores_in_current_position_index;
		/*
		 printf ("\n%d quality_score_position_max %d i=%d" , rle_quality_scores_index , quality_score_position_max[rle_quality_scores_index] , i);
		 fflush ( stdout);
		 */
		num = 0;
		for ( j = 0 ; j < quality_score_position_max[rle_quality_scores_index] ;
				j++ )
			num += rle_quality_scores[rle_quality_scores_index][j].frequency;
		/*
		 *
		 *
		 * Calculate sum of each position
		 */

		//printf ("\nSum Position %d %lld" , rle_quality_scores_index , num);
		//fflush ( stdout);
		rle_quality_scores_index++;
	}
	//printf ("\nComing out of this loop");
	//fflush ( stdout);
	max_read_length = rle_quality_scores_index;
	quality_score_of_read = ( char* ) malloc (sizeof(char) * max_read_length);

	/*
	 * Print first line
	 */
	/*
	 for ( i = 0 ; i < max_read_length ; i++ )
	 {
	 for ( j = 0 ; j < quality_score_position_max[i] ; j++ )
	 printf ("\n%d %c" , rle_quality_scores[i][j].frequency , rle_quality_scores[i][j].quality_score);
	 }*/

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
			checker_flag += rle_quality_scores[i][quality_score_position_index[i]].frequency;
			//printf ("\n Frequency being added %d i=%d j=%d" , rle_quality_scores[i][quality_score_position_index[i]].frequency , i , quality_score_position_index[i]);
			if ( rle_quality_scores[i][quality_score_position_index[i]].frequency == 0 )
				quality_score_position_index[i]++;
			if ( quality_score_position_index[i] >= quality_score_position_max[i] )
			{
				printf ("\nTrouble %d %d" , quality_score_position_index[i] , quality_score_position_max[i]);
				fflush (stdout);
			}
		}
		//printf ("\nChecker flag %d" , checker_flag);
		quality_score_of_read[i] = '\0';
		//printf ( "\n%s" , quality_score_of_read );
		//fflush ( stdout );
		fprintf (fhw , "%s" , quality_score_of_read);
		fprintf (fhw , "%s" , "\n");
		if ( checker_flag == 0 ) break;
	}
	//printf ("\nReturning...");
	fflush ( stdout);
	fclose (fhr);
	fclose (fhw);
}

int main (int argc, char *argv[])
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
	strcpy(input_qualityscore_filename , argv[1]);
	strcpy(output_quality_score_filename , argv[2]);
	/********************************************************************/

	convertRLEtoQualValues (input_qualityscore_filename , output_quality_score_filename);
}
