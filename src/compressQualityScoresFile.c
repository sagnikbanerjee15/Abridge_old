# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

void performColumnWiseRLE ( char *input_qualityscore_filename, char *output_quality_score_filename, short int adjust_quality_scores )
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	FILE *fhr;
	FILE *fhw;

	struct Quality_Score_RLE **qsRLE;

	size_t len = 0;
	ssize_t line_len;

	char *line;
	char str[1000];
	char **lines_to_be_written_to_file;

	int i;
	int max_len_sequence = 0;
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
	qsRLE = ( struct Quality_Score_RLE** ) malloc ( sizeof(struct Quality_Score_RLE*) * MAX_SEQ_LEN );
	for ( i = 0 ; i < MAX_SEQ_LEN ; i++ )
		qsRLE[i] = ( struct Quality_Score_RLE* ) malloc ( sizeof(struct Quality_Score_RLE) );

	lines_to_be_written_to_file = ( char** ) malloc ( sizeof(char*) * MAX_SEQ_LEN );
	for ( i = 0 ; i < MAX_SEQ_LEN ; i++ )
	{
		lines_to_be_written_to_file[i] = ( char* ) malloc ( sizeof(char) * 100000 );
		lines_to_be_written_to_file[i][0] = '\0';
	}
	/********************************************************************/

	while ( ( line_len = getline ( &line , &len , fhr ) ) != -1 )
	{
		for ( i = 0 ; line[i] != '\0' ; i++ )
		{
			if ( i > max_len_sequence ) max_len_sequence = i;
			if ( qsRLE[i]->score_character == 'X' && qsRLE[i]->frequency == 0 )
			{
				qsRLE[i]->score_character = line[i];
				qsRLE[i]->frequency++;
			}
			else if ( line[i] == qsRLE[i]->score_character )
			{
				qsRLE[i]->frequency++;
			}
			else
			{
				sprintf( str , "%lld" , qsRLE[i]->frequency );
				strncat( str , &qsRLE[i]->score_character , 1 );
				strcat( lines_to_be_written_to_file[i] , str );
				str[0] = '\0';
				qsRLE[i]->frequency = 1;
				qsRLE[i]->score_character = line[i];
			}
		}
	}
	for ( i = 0 ; line[i] != '\0' ; i++ )
	{
		sprintf( str , "%lld" , qsRLE[i]->frequency );
		strncat( str , &qsRLE[i]->score_character , 1 );
		strcat( lines_to_be_written_to_file[i] , str );
		str[0] = '\0';
		qsRLE[i]->frequency = 0;
		qsRLE[i]->score_character = 'X';
	}
	for ( i = 0 ; i < max_len_sequence ; i++ )
	{
		fprintf ( fhw , "%s" , lines_to_be_written_to_file[i] );
		fprintf ( fhw , "%s" , "\n" );
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

	short int adjust_quality_scores;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	strcpy( input_qualityscore_filename , argv[1] );
	strcpy( output_quality_score_filename , argv[2] );
	adjust_quality_scores = strtol ( argv[3] , &temp , 10 );
	/********************************************************************/

	performColumnWiseRLE ( input_qualityscore_filename , output_quality_score_filename , adjust_quality_scores );
}
