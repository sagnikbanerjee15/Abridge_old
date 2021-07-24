# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

void expandMDString (char *icigar, char *change_indicator, int val)
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
	//extractNHfromicigar (icigar , strlen (icigar));
	num = 0;
	for ( i = 0 ; icigar[i] != '\0' ; i++ )
	{
		if ( isCharacterInString ("atgcn" , icigar[i]) )
		{
			change_indicator[change_indicator_index++ ] = '0';
			/*
			 if ( val == 3540189 )
			 {
			 printf ("\n1.change_indicator_index = %d %s" , change_indicator_index , icigar);
			 fflush (stdout);
			 }
			 */
		}
		else if ( isdigit (icigar[i]) != 0 )
		{
			num = icigar[i] - 48;
			i++;
			while ( isdigit (icigar[i]) != 0 )
			{
				num = num * 10 + icigar[i] - 48;
				i++;
			}
			if ( icigar[i] != 'N' && icigar[i] != 'D' )
			{
				if ( icigar[i] == 'I' || icigar[i] == 'S' )
				{
					while ( num-- )
					{
						change_indicator[change_indicator_index++ ] = '0';
						/*
						 if ( val == 3540189 )
						 {
						 printf ("\n2.change_indicator_index = %d %s i=%d" , change_indicator_index , icigar , i);
						 fflush (stdout);
						 }
						 */
					}
				}
				else
				{
					while ( num-- )
					{
						change_indicator[change_indicator_index++ ] = '1';
						/*
						 if ( val == 3540189 )
						 {
						 printf ("\n3.change_indicator_index = %d %s i=%d" , change_indicator_index , icigar , i);
						 fflush (stdout);
						 }
						 */
					}
				}
			}
			/*
			 if ( val == 3540189 )
			 {
			 printf ("\n4.change_indicator_index = %d %s i=%d" , change_indicator_index , icigar , i);
			 fflush (stdout);
			 }
			 */
		}
		else if ( isCharacterInString (insert_characters , icigar[i]) || isCharacterInString (mismatch_characters , icigar[i]) )
		{
			change_indicator[change_indicator_index++ ] = '0';
			/*
			 if ( val == 3540189 )
			 {
			 printf ("\n5.change_indicator_index = %d %s i=%d" , change_indicator_index , icigar , i);
			 fflush (stdout);
			 }
			 */
		}
	}
	change_indicator[change_indicator_index++ ] = '\0';
	/*
	 if ( val == 3540189 )
	 {
	 printf ("\n6.change_indicator_index = %d %s i=%d" , change_indicator_index , icigar , i);
	 fflush (stdout);
	 printf ("\nReturning");
	 fflush (stdout);
	 }
	 */
}

void performColumnWiseRLE (
		char *input_qualityscore_filename,
		char *output_quality_score_filename,
		short int save_exact_quality_scores)
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
	char temp;

	int i, j, k;
	int max_len_sequence = 0;
	int *lines_to_be_written_to_file_index;
	int *count_max_reads_each_position;
	int line_number;

	long long int number_of_records_read = 0;

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

	qsRLE = ( struct Quality_Score_RLE** ) malloc (sizeof(struct Quality_Score_RLE*) * MAX_SEQ_LEN);
	for ( i = 0 ; i < MAX_SEQ_LEN ; i++ )
		qsRLE[i] = allocateMemoryQuality_Score_RLE ();

	lines_to_be_written_to_file_index = ( int* ) malloc (sizeof(int) * MAX_SEQ_LEN);
	lines_to_be_written_to_file = ( char** ) malloc (sizeof(char*) * MAX_SEQ_LEN);
	for ( i = 0 ; i < MAX_SEQ_LEN ; i++ )
	{
		lines_to_be_written_to_file[i] = ( char* ) malloc (sizeof(char) * 1000000);
		lines_to_be_written_to_file[i][0] = '\0';
		lines_to_be_written_to_file_index[i] = 0;
	}

	split_on_tab = ( char** ) malloc (sizeof(char*) * 5);
	for ( i = 0 ; i < 5 ; i++ )
		split_on_tab[i] = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	change_indicator = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	count_max_reads_each_position = ( int* ) malloc (sizeof(int) * MAX_SEQ_LEN);
	for ( i = 0 ; i < MAX_SEQ_LEN ; i++ )
		count_max_reads_each_position[i] = 0;
	/********************************************************************/

	line_len = getline ( &line , &len , fhr);
	splitByDelimiter (line , '\t' , split_on_tab);
	for ( i = 0 ; split_on_tab[0][i] != '\0' ; i++ )
	{
		//if ( i > max_len_sequence ) max_len_sequence = i;
		qsRLE[i]->score_character = split_on_tab[0][i];
		qsRLE[i]->frequency = 1;
	}

	//printf ( "\n%s\t%s\n%s" , split_on_tab[0] , split_on_tab[1] , change_indicator );
	//fflush ( stdout );
	line_number = 1;

	line_len = getline ( &line , &len , fhr);
	do
	{
		line_number++;
		/*
		 if ( line_number % 10000 == 0 )
		 printf ("\nProcessing line number %d" , line_number);
		 */
		splitByDelimiter (line , '\t' , split_on_tab);
		expandMDString (split_on_tab[1] , change_indicator , line_number);

		if ( strcmp (split_on_tab[2] , "2") == 0 ) // Reverse the change indicator
		{
			i = 0;
			j = strlen (change_indicator) - 1;
			while ( i < j )
			{
				temp = change_indicator[i];
				change_indicator[i++ ] = change_indicator[j];
				change_indicator[j-- ] = temp;
			}
		}
		continue;
		//printf ( "\n%s\t%s\n%s" , split_on_tab[0] , split_on_tab[1] , change_indicator );
		//fflush ( stdout );
		for ( i = 0 ; split_on_tab[0][i] != '\0' ; i++ )
		{
			if ( i > max_len_sequence ) max_len_sequence = i;
			//printf ( "\ni=%d split_on_tab[0][i]=%c qsRLE[i]->score_character=%c " , i , split_on_tab[0][i] , qsRLE[i]->score_character );
			if ( split_on_tab[0][i] == qsRLE[i]->score_character || ( save_exact_quality_scores == 0 && change_indicator[i] == '1' ) )
				qsRLE[i]->frequency++;
			else
			{
				if ( qsRLE[i]->frequency > 1 )
				{
					sprintf(str , "%lld" , qsRLE[i]->frequency);
					for ( j = 0 ; str[j] != '\0' ; j++ )
						lines_to_be_written_to_file[i][lines_to_be_written_to_file_index[i]++ ] = str[j];
				}
				count_max_reads_each_position[i] += qsRLE[i]->frequency;
				lines_to_be_written_to_file[i][lines_to_be_written_to_file_index[i]++ ] = qsRLE[i]->score_character + 30;// Done to avoid ambiguity between numeric quality scores and frequency value
				lines_to_be_written_to_file[i][lines_to_be_written_to_file_index[i]] = '\0';
				str[0] = '\0';
				qsRLE[i]->frequency = 1;
				qsRLE[i]->score_character = split_on_tab[0][i];
			}
		}
	} while ( ( line_len = getline ( &line , &len , fhr) ) != -1 );

	return;
	for ( i = 0 ; split_on_tab[0][i] != '\0' ; i++ )
	{
		if ( qsRLE[i]->frequency > 1 )
		{
			sprintf(str , "%lld" , qsRLE[i]->frequency);
			for ( j = 0 ; str[j] != '\0' ; j++ )
				lines_to_be_written_to_file[i][lines_to_be_written_to_file_index[i]++ ] = str[j];
		}
		count_max_reads_each_position[i] += qsRLE[i]->frequency;
		lines_to_be_written_to_file[i][lines_to_be_written_to_file_index[i]++ ] = qsRLE[i]->score_character + 30;// Done to avoid ambiguity between numeric quality scores and frequency value
		lines_to_be_written_to_file[i][lines_to_be_written_to_file_index[i]] = '\0';
		str[0] = '\0';
		qsRLE[i]->frequency = 0;
		qsRLE[i]->score_character = 'X';
	}
	max_len_sequence++;
	for ( i = 0 ; i < max_len_sequence ; i++ )
		printf ("\nMAX LEN %d %d" , i , count_max_reads_each_position[i]);

	//max_len_sequence = 151;
	for ( i = 0 ; i < max_len_sequence ; i++ )
	{
		fprintf (fhw , "%s" , lines_to_be_written_to_file[i]);
		fprintf (fhw , "%s" , "\n");
		printf ("\nPosition %d Index position %d Length of string %d" , i , lines_to_be_written_to_file_index[i] , strlen (lines_to_be_written_to_file[i]));
	}
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

	short int save_exact_quality_scores;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	strcpy(input_qualityscore_filename , argv[1]);
	strcpy(output_quality_score_filename , argv[2]);
	save_exact_quality_scores = strtol (argv[3] , &temp , 10);
	/********************************************************************/

	performColumnWiseRLE (input_qualityscore_filename , output_quality_score_filename , save_exact_quality_scores);
}
