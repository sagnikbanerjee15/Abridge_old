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
		switch ( icigar[i] )
		{
			case 'a':
			case 't':
			case 'g':
			case 'c':
			case 'n':
				change_indicator[change_indicator_index++ ] = '0';
				break;
		}
		if ( isdigit (icigar[i]) != 0 ) // icigar[i] is numeric
		{
			num = icigar[i] - 48;
			i++;
			while ( isdigit (icigar[i]) != 0 && icigar[i] != '\0' )
			{
				num = num * 10 + icigar[i] - 48;
				i++;
			}
			if ( icigar[i] == '\0' ) break;
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
		short int save_exact_quality_scores,
		int max_read_length)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	FILE *fhr;
	FILE *fhw;
	FILE **fhw_each_position;
	FILE **fhr_each_position;

	struct Quality_Score_RLE **qsRLE;

	size_t len = 0;
	ssize_t line_len;

	char str[1000];
	char *line;
	char *output_filename_for_each_position;
	//char **lines_to_be_written_to_file;
	char **split_on_tab;
	char *change_indicator;
	char *qual_score;
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
	output_filename_for_each_position = ( char* ) malloc (sizeof(char) * 1000);
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
	qual_score = ( char* ) malloc (sizeof(char) * ( max_read_length + 5 ));
	fhw_each_position = ( FILE** ) malloc (sizeof(FILE*) * max_read_length);
	fhr_each_position = ( FILE** ) malloc (sizeof(FILE*) * max_read_length);

	/*
	 * Create quality score files for each position of read
	 */
	for ( i = 0 ; i < max_read_length ; i++ )
	{
		sprintf(str , "%ld" , i);
		strcpy(output_filename_for_each_position , output_quality_score_filename);
		strcat(output_filename_for_each_position , str);
		fhw_each_position[i] = fopen (output_filename_for_each_position , "w");
		if ( fhw_each_position[i] == NULL )
		{
			printf ("%s File cannot be created" , strcat(output_quality_score_filename , str));
			exit (1);
		}
	}

	qsRLE = ( struct Quality_Score_RLE** ) malloc (sizeof(struct Quality_Score_RLE*) * ( max_read_length + 5 ));
	for ( i = 0 ; i < max_read_length ; i++ )
		qsRLE[i] = allocateMemoryQuality_Score_RLE ();

	/*
	 lines_to_be_written_to_file_index = ( int* ) malloc (sizeof(int) * max_read_length);
	 lines_to_be_written_to_file = ( char** ) malloc (sizeof(char*) * max_read_length);
	 for ( i = 0 ; i < max_read_length ; i++ )
	 {
	 lines_to_be_written_to_file[i] = ( char* ) malloc (sizeof(char) * 1000000);
	 lines_to_be_written_to_file[i][0] = '\0';
	 lines_to_be_written_to_file_index[i] = 0;
	 }*/

	split_on_tab = ( char** ) malloc (sizeof(char*) * 5);
	for ( i = 0 ; i < 5 ; i++ )
		split_on_tab[i] = ( char* ) malloc (sizeof(char) * ( max_read_length + 5 ));
	change_indicator = ( char* ) malloc (sizeof(char) * ( max_read_length + 5 ));
	count_max_reads_each_position = ( int* ) malloc (sizeof(int) * ( max_read_length + 5 ));
	for ( i = 0 ; i < max_read_length ; i++ )
		count_max_reads_each_position[i] = 0;
	/********************************************************************/

	line_number = 1;
	line_len = getline ( &line , &len , fhr);
	splitByDelimiter (line , '\t' , split_on_tab);
	strcpy(qual_score , split_on_tab[0]);
	for ( i = 0 ; qual_score[i] != '\0' ; i++ )
	{
		//if ( i > max_len_sequence ) max_len_sequence = i;
		qsRLE[i]->score_character = qual_score[i];
		qsRLE[i]->frequency = 1;
		count_max_reads_each_position[i]++;
	}

	while ( ( line_len = getline ( &line , &len , fhr) ) != -1 )
	{
		//if ( line_number == 100 ) break;
		line_number++;
		splitByDelimiter (line , '\t' , split_on_tab);
		strcpy(qual_score , split_on_tab[0]);
		//printf ("\nCIGAR = %s" , split_on_tab[1]);
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

		//printf ("\n%s %d %s" , change_indicator , line_number , split_on_tab[1]);

		/*
		 * Append blank spaces of short sized reads
		 */
		if ( strlen (qual_score) < max_read_length )
		{
			for ( i = strlen (qual_score) ; i < max_read_length ; i++ )
			{
				qual_score[i] = ' ';
				change_indicator[i] = '0';
			}
			qual_score[i] = '\0';
			change_indicator[i] = '\0';
		}

		for ( i = 0 ; qual_score[i] != '\0' ; i++ )
		{
			if ( qual_score[i] != ' ' && qsRLE[i]->score_character == ' ' )
				change_indicator[i] = '0';

			if ( qual_score[i] == ' ' && qsRLE[i]->score_character == ' ' )
				qsRLE[i]->frequency++;
			else if ( qual_score[i] != ' ' && ( qual_score[i] == qsRLE[i]->score_character || ( save_exact_quality_scores == 0 && change_indicator[i] == '1' ) ) )
				qsRLE[i]->frequency++;
			else
			{
				if ( qsRLE[i]->frequency > 1 )
				{
					sprintf(str , "%lld" , qsRLE[i]->frequency);
					fprintf (fhw_each_position[i] , "%s" , str);
				}
				count_max_reads_each_position[i] += qsRLE[i]->frequency;
				fputc (qsRLE[i]->score_character + 26 , fhw_each_position[i]);
				//printf ("\n %c %c" , qsRLE[i]->score_character , qsRLE[i]->score_character + 26);
				qsRLE[i]->frequency = 1;
				qsRLE[i]->score_character = qual_score[i];
			}
		}
	}

	for ( i = 0 ; i < max_read_length ; i++ )
	{
		if ( qsRLE[i]->frequency > 1 )
		{
			sprintf(str , "%lld" , qsRLE[i]->frequency);
			fprintf (fhw_each_position[i] , "%s" , str);
		}
		count_max_reads_each_position[i] += qsRLE[i]->frequency;
		fputc (qsRLE[i]->score_character + 26 , fhw_each_position[i]);
		fputc ('\n' , fhw_each_position[i]);

	}
	/*
	 for ( i = 0 ; i < max_read_length ; i++ )
	 printf ("\nMAX NUM READS IN POS %d %d" , i + 1 , count_max_reads_each_position[i]);
	 */
	for ( i = 0 ; i < max_read_length ; i++ )
		fclose (fhw_each_position[i]);
	fclose (fhr);
	fclose (fhw);

	/*
	 * Merge the RLE together and remove the intermediate files
	 */

	fhw = fopen (output_quality_score_filename , "w");
	for ( i = 0 ; i < max_read_length ; i++ )
	{
		sprintf(str , "%ld" , i);
		strcpy(output_filename_for_each_position , output_quality_score_filename);
		strcat(output_filename_for_each_position , str);
		fhr = fopen (output_filename_for_each_position , "r");
		if ( fhr == NULL )
		{
			printf ("%s File cannot be opened for reading" , strcat(output_filename_for_each_position , str));
			exit (1);
		}

		//Read the whole file
		line = NULL;
		len = 0;
		line_len = getline ( &line , &len , fhr);
		printf ("\nLine read: %s" , line);
		//strcat(line , "\n");
		fprintf (fhw , "%s" , line);
		fprintf (fhw , "%s" , "\n");
		fclose (fhr);
		remove (output_filename_for_each_position);
	}
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
	int max_read_length;

	short int save_exact_quality_scores;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	strcpy(input_qualityscore_filename , argv[1]);
	strcpy(output_quality_score_filename , argv[2]);
	max_read_length = strtol (argv[3] , &temp , 10);
	save_exact_quality_scores = strtol (argv[4] , &temp , 10);

	/********************************************************************/

	performColumnWiseRLE (input_qualityscore_filename , output_quality_score_filename , save_exact_quality_scores , max_read_length);
}
