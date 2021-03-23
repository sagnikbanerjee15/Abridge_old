# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

void expandMDString ( char *icigar, char *change_indicator )
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	int i;
	int num;
	int change_indicator_index = 0;

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/

	/********************************************************************/
	extractNHfromicigar ( icigar , strlen ( icigar ) );
	num = 0;
	for ( i = 0 ; icigar[i] != '\0' ; i++ )
	{
		if ( isCharacterInString ( "atgcn" , icigar[i] ) )
		{
			change_indicator[change_indicator_index++ ] = '0';
		}
		else if ( isdigit ( icigar[i] ) != 0 )
		{
			num = icigar[i] - 48;
			i++;
			while ( isdigit ( icigar[i] ) != 0 )
			{
				num = num * 10 + icigar[i] - 48;
				i++;
			}
			if ( icigar[i] != 'N' && icigar[i] != 'D' )
			{
				if ( icigar[i] == 'I' )
					while ( num-- )
						change_indicator[change_indicator_index++ ] = '0';
				else while ( num-- )
					change_indicator[change_indicator_index++ ] = '1';
			}
		}
		else if ( isCharacterInString ( insert_characters , icigar[i] ) || isCharacterInString ( mismatch_characters , icigar[i] ) )
			change_indicator[change_indicator_index++ ] = '0';
	}
	change_indicator[change_indicator_index++ ] = '\0';
}

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

	char str[1000];
	char *line;
	char **lines_to_be_written_to_file;
	char **split_on_tab;
	char *change_indicator;

	int i , j , k;
	int max_len_sequence = 0;
	int *lines_to_be_written_to_file_index;

	long long int number_of_records_read = 0;

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
		qsRLE[i] = allocateMemoryQuality_Score_RLE ();

	lines_to_be_written_to_file_index = ( int* ) malloc ( sizeof(int) * MAX_SEQ_LEN );
	lines_to_be_written_to_file = ( char** ) malloc ( sizeof(char*) * MAX_SEQ_LEN );
	for ( i = 0 ; i < MAX_SEQ_LEN ; i++ )
	{
		lines_to_be_written_to_file[i] = ( char* ) malloc ( sizeof(char) * 100000000 );
		lines_to_be_written_to_file[i][0] = '\0';
		lines_to_be_written_to_file_index[i] = 0;
	}

	split_on_tab = ( char** ) malloc ( sizeof(char*) * 3 );
	for ( i = 0 ; i < 3 ; i++ )
		split_on_tab[i] = ( char* ) malloc ( sizeof(char) * MAX_SEQ_LEN );
	change_indicator = ( char* ) malloc ( sizeof(char) * MAX_SEQ_LEN );

	/********************************************************************/

	line_len = getline ( &line , &len , fhr );
	splitByDelimiter ( line , '\t' , split_on_tab );
	for ( i = 0 ; split_on_tab[0][i] != '\0' ; i++ )
	{
		//if ( i > max_len_sequence ) max_len_sequence = i;
		qsRLE[i]->score_character = split_on_tab[0][i];
		qsRLE[i]->frequency++;
	}

	line_len = getline ( &line , &len , fhr );
	do
	{
		splitByDelimiter ( line , '\t' , split_on_tab );
		expandMDString ( split_on_tab[1] , change_indicator );
		printf ( "\n%s\t%s\n%s" , split_on_tab[0] , split_on_tab[1] , change_indicator );
		fflush ( stdout );
		for ( i = 0 ; split_on_tab[0][i] != '\0' ; i++ )
		{
			//if ( i > max_len_sequence ) max_len_sequence = i;
			if ( split_on_tab[0][i] == qsRLE[i]->score_character )
				qsRLE[i]->frequency++;
			else
			{
				if ( qsRLE[i]->frequency > 1 )
				{
					sprintf( str , "%lld" , qsRLE[i]->frequency );
					for ( j = 0 ; str[j] != '\0' ; j++ )
						lines_to_be_written_to_file[i][lines_to_be_written_to_file_index[i]++ ] = str[j];
				}
				lines_to_be_written_to_file[i][lines_to_be_written_to_file_index[i]++ ] = qsRLE[i]->score_character;
				lines_to_be_written_to_file[i][lines_to_be_written_to_file_index[i]] = '\0';
				str[0] = '\0';
				qsRLE[i]->frequency = 1;
				qsRLE[i]->score_character = split_on_tab[0][i];
			}
		}
	} while ( ( line_len = getline ( &line , &len , fhr ) ) != -1 );

	for ( i = 0 ; split_on_tab[0][i] != '\0' ; i++ )
	{
		if ( qsRLE[i]->frequency > 1 )
		{
			sprintf( str , "%lld" , qsRLE[i]->frequency );
			for ( j = 0 ; str[j] != '\0' ; j++ )
				lines_to_be_written_to_file[i][lines_to_be_written_to_file_index[i]++ ] = str[j];
		}
		lines_to_be_written_to_file[i][lines_to_be_written_to_file_index[i]++ ] = qsRLE[i]->score_character;
		lines_to_be_written_to_file[i][lines_to_be_written_to_file_index[i]] = '\0';
		str[0] = '\0';
		qsRLE[i]->frequency = 0;
		qsRLE[i]->score_character = 'X';
	}

	max_len_sequence = 151;
	for ( i = 0 ; i < max_len_sequence ; i++ )
	{
		fprintf ( fhw , "%s" , lines_to_be_written_to_file[i] );
		fprintf ( fhw , "%s" , "\n" );
		printf ( "\nPosition %d Index position %d Length of string %d" , i , lines_to_be_written_to_file_index[i] , strlen ( lines_to_be_written_to_file[i] ) );
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
