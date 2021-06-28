# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include <math.h>
# include <stdbool.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

int main (int argc, char *argv[])
{
	/****************************************************************************************************************************************
	 * Variable declaration
	 ****************************************************************************************************************************************/
	char abridge_index_filename[FILENAME_LENGTH];
	char genome_filename[FILENAME_LENGTH];
	char genome_index_filename[FILENAME_LENGTH];
	char chromosome[1000];
	char pass1_filename[FILENAME_LENGTH];
	char output_sam_filename[FILENAME_LENGTH];
	char default_quality_value[10];
	char read_prefix[1000];
	char name_of_file_with_quality_scores[FILENAME_LENGTH];

	char *temp;
	char *buffer_for_pass1;
	char *buffer_for_qual;

	char **split_on_newline;
	char **split_on_tab;
	char **split_on_dash;
	char **split_on_comma;
	char **split_on_delimiter;
	char **split_on_newline_qual;

	struct Abridge_Index *abridge_index;
	struct Abridge_Index *genome_index;
	struct Whole_Genome_Sequence *single_genome_sequence;
	struct Sam_Alignment *sam_alignment;

	short int flag_ignore_mismatches;
	short int flag_ignore_soft_clippings;
	short int flag_ignore_unmapped_sequences;
	short int flag_ignore_quality_score;
	short int flag_ignore_sequence_information;
	short int flag_save_all_quality_scores;
	short int flag_save_exact_quality_scores;
	short int flag_save_scores;
	short int number_of_columns;
	short int first_record;
	short int read_names_stored;

	int i, j, k, l;
	int fseek_ret_val;
	int fread_ret_val;
	int number_of_entries_in_cluster;
	int sam_alignment_pool_index;
	int split_on_newline_qual_ROWS = ROWS;
	int split_on_newline_qual_COLS = COLS * 10;
	int number_of_newlines = 0;
	int number_of_distinct_cigars_in_a_line;
	int number_of_repititions_of_the_same_reads;
	int qual_pool_iterator;
	int total_quals;

	unsigned long long int read_number = 1;
	unsigned long long int from = -1;
	unsigned long long int to = -1;

	long long int start;
	long long int end;
	long long int abridge_match_start_index; // Start index of abridge index - need for multiple cluster overlaps
	long long int abridge_match_end_index; //
	long long int total_mapped_reads = 0;
	long long int curr_position;

	FILE *fhr_pass1;
	FILE *fhr_qual;
	FILE *fhw;
	/****************************************************************************************************************************************
	 * Variable initialization
	 ****************************************************************************************************************************************/

	strcpy(abridge_index_filename , argv[1]);
	strcpy(genome_filename , argv[2]);
	strcpy(genome_index_filename , argv[3]);
	strcpy(chromosome , argv[4]);
	strcpy(pass1_filename , argv[5]);
	start = strtol (argv[6] , &temp , 10);
	end = strtol (argv[7] , &temp , 10);
	strcpy(output_sam_filename , argv[8]);
	strcpy(default_quality_value , argv[9]);
	flag_ignore_sequence_information = strtol (argv[10] , &temp , 10);
	strcpy(read_prefix , argv[11]);
	strcpy(name_of_file_with_quality_scores , argv[12]);

	fhr_pass1 = fopen (pass1_filename , "rb");
	if ( fhr_pass1 == NULL )
	{
		printf ("Error! File %s not found" , pass1_filename);
		exit (1);
	}
	fhr_qual = fopen (name_of_file_with_quality_scores , "rb");
	if ( fhr_qual == NULL )
	{
		printf ("Error! File %s not found" , name_of_file_with_quality_scores);
		exit (1);
	}
	fhw = fopen (output_sam_filename , "w");
	if ( fhw == NULL )
	{
		printf ("Error! File %s cannot be opened for writing" , output_sam_filename);
		exit (1);
	}

	abridge_index = allocateMemoryAbridge_Index ();
	genome_index = allocateMemoryAbridge_Index ();

	split_on_newline_qual = ( char** ) malloc (sizeof(char*) * split_on_newline_qual_ROWS);
	for ( i = 0 ; i < split_on_newline_qual_ROWS ; i++ )
		split_on_newline_qual[i] = ( char* ) malloc (sizeof(char) * split_on_newline_qual_COLS);

	split_on_newline = ( char** ) malloc (sizeof(char*) * ROWS * 10);
	for ( i = 0 ; i < ROWS * 10 ; i++ )
		split_on_newline[i] = ( char* ) malloc (sizeof(char) * COLS * 10);

	split_on_tab = ( char** ) malloc (sizeof(char*) * 5);
	for ( i = 0 ; i < 5 ; i++ )
		split_on_tab[i] = ( char* ) malloc (sizeof(char) * COLS * 100);

	split_on_dash = ( char** ) malloc (sizeof(char*) * 5);
	for ( i = 0 ; i < 5 ; i++ )
		split_on_dash[i] = ( char* ) malloc (sizeof(char) * COLS);

	split_on_comma = ( char** ) malloc (sizeof(char*) * ROWS * 10);
	for ( i = 0 ; i < ROWS * 10 ; i++ )
		split_on_comma[i] = ( char* ) malloc (sizeof(char) * COLS * 10);

	split_on_delimiter = ( char** ) malloc (sizeof(char*) * ROWS);
	for ( i = 0 ; i < ROWS ; i++ )
		split_on_delimiter[i] = ( char* ) malloc (sizeof(char) * COLS);

	buffer_for_pass1 = ( char* ) malloc (sizeof(char) * MAX_BUFFER_SIZE_FOR_READING_PASS2_FILE);
	single_genome_sequence = ( struct Whole_Genome_Sequence* ) malloc (sizeof(struct Whole_Genome_Sequence));
	single_genome_sequence->number_of_reference_sequences = 0;
	single_genome_sequence->nucleotides = ( char** ) malloc (sizeof(char*) * 1);
	single_genome_sequence->reference_sequence_name = ( char** ) malloc (sizeof(char*) * 1);
	single_genome_sequence->reference_sequence_length = ( unsigned long long int* ) malloc (sizeof(unsigned long long int) * 1);
	sam_alignment = allocateMemorySam_Alignment ();

	/****************************************************************************************************************************************/
	readAbridgeIndex (abridge_index , abridge_index_filename , split_on_newline , &flag_ignore_mismatches , &flag_ignore_soft_clippings , &flag_ignore_unmapped_sequences , &flag_ignore_quality_score , &flag_save_all_quality_scores , &flag_save_exact_quality_scores , &flag_save_scores);
	readGenomeIndex (genome_index , genome_index_filename , split_on_newline);
	//readInGenomeSequenceSingleChromosome (single_genome_sequence , chromosome , genome_filename , genome_index);
	readInEachChromosome (genome_filename , single_genome_sequence , chromosome);
	first_record = findReadClusterFromAbridgeIndex (abridge_index , chromosome , start , end , &abridge_match_start_index , &abridge_match_end_index);
	writeSequenceHeaders (fhw , genome_filename , 0);

	buffer_for_qual = ( char* ) malloc (sizeof(char) * MAX_BUFFER_SIZE_FOR_READING_PASS2_FILE);
	from = start;
	to = end;
	/*
	 printf ("\n abridge_match_start_index %d abridge_match_end_index %d" , abridge_match_start_index , abridge_match_end_index);
	 fflush (stdout);
	 */
	for ( i = abridge_match_start_index ; i <= abridge_match_end_index ; i++ )
	{
		/*
		 printf ("\n%s %lld %lld %lld %lld %lld %lld" , abridge_index->chromosome[i] , abridge_index->start[i] , abridge_index->end[i] , abridge_index->start_byte[i] , abridge_index->end_byte[i] , abridge_index->start_byte_qual[i] , abridge_index->end_byte_qual[i]);
		 fflush (stdout);
		 */
		fseek_ret_val = fseek (fhr_pass1 , abridge_index->start_byte[i] , SEEK_SET);
		buffer_for_pass1[0] = '\0';
		fread_ret_val = fread (buffer_for_pass1 , 1 , abridge_index->end_byte[i] - abridge_index->start_byte[i] , fhr_pass1);
		//printf("\n fread_ret_val %d fseek_ret_val %d diff %lld", fread_ret_val, fseek_ret_val, MAX_BUFFER_SIZE_FOR_READING_PASS2_FILE - fread_ret_val);
		//fflush(stdout);
		buffer_for_pass1[fread_ret_val] = '\0';
		number_of_entries_in_cluster = splitByDelimiter (buffer_for_pass1 , '\n' , split_on_newline);
		number_of_entries_in_cluster--; //Last line is always empty

		//free (buffer_for_qual);
		//buffer_for_qual = ( char* ) malloc (sizeof(char) * ( abridge_index->end_byte_qual[i] - abridge_index->start_byte_qual[i] ));
		buffer_for_qual[0] = '\0';
		fseek_ret_val = fseek (fhr_qual , abridge_index->start_byte_qual[i] , SEEK_SET);
		fread_ret_val = fread (buffer_for_qual , 1 , abridge_index->end_byte_qual[i] - abridge_index->start_byte_qual[i] , fhr_qual);
		/*
		 for ( j = 0 ; j < split_on_newline_qual_ROWS ; j++ )
		 free (split_on_newline_qual[i]);
		 free (split_on_newline_qual);
		 */
		number_of_newlines = 0;
		for ( j = 0 ; buffer_for_qual[j] != '\0' ; j++ )
			if ( buffer_for_qual[j] == '\n' ) number_of_newlines += 1;
		split_on_newline_qual_ROWS = number_of_newlines + 10;
		split_on_newline_qual_COLS = MAX_SEQ_LEN;
		split_on_newline_qual = ( char** ) malloc (sizeof(char*) * split_on_newline_qual_ROWS);
		for ( j = 0 ; j < split_on_newline_qual_ROWS ; j++ )
			split_on_newline_qual[j] = ( char* ) malloc (sizeof(char) * split_on_newline_qual_COLS);

		total_quals = splitByDelimiter (buffer_for_qual , '\n' , split_on_newline_qual);
		/*
		 for ( j = 0 ; j < total_quals ; j++ )
		 printf ("\n%s" , split_on_newline_qual[j]);
		 */
		qual_pool_iterator = 0;
		/*
		 printf ("\nfread_ret_val %d" , fread_ret_val);
		 printf ("\nNew Record:%d\n" , i);
		 printf ("%s" , buffer_for_pass1);
		 printf ("\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
		 printf ("\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
		 fflush (stdout);
		 */
		curr_position = abridge_index->start[i];
		//if ( first_record == 1 ) curr_position = 0;
		/*
		 printf ("\ncurr_position First %d i=%d" , curr_position , i);
		 */
		for ( j = 0 ; j < number_of_entries_in_cluster ; j++ )
		{
			if ( strstr (split_on_newline[j] , "abridge_") )
				read_names_stored = 1;
			else read_names_stored = 0;

			number_of_columns = splitByDelimiter (split_on_newline[j] , '\t' , split_on_tab);

			if ( read_names_stored == 0 )
			{
				if ( number_of_columns == 1 )
					curr_position++;
				else if ( number_of_columns == 2 )
				{
					if ( j != 0 )
						curr_position += strtol (split_on_tab[0] , &temp , 10);
				}
			}
			else if ( read_names_stored == 1 )
			{
				if ( number_of_columns == 2 )
					curr_position++;
				else if ( number_of_columns == 2 )
				{
					if ( j != 0 )
						curr_position += strtol (split_on_tab[0] , &temp , 10);
				}
			}

			if ( ( number_of_columns == 1 && read_names_stored == 0 ) || ( number_of_columns == 2 && read_names_stored == 1 ) )
			{
				int max_number_of_commas = 0, number_of_commas = 0;
				for ( i = 0 ; split_on_tab[0][i] != '\0' ; i++ )
					if ( split_on_tab[0][i] == ',' ) number_of_commas++;
				number_of_distinct_cigars_in_a_line = splitByDelimiter (split_on_tab[0] , ',' , split_on_comma);
			}
			else if ( ( number_of_columns == 2 && read_names_stored == 0 ) || ( number_of_columns == 3 && read_names_stored == 1 ) )
			{
				int max_number_of_commas = 0, number_of_commas = 0;
				for ( i = 0 ; split_on_tab[1][i] != '\0' ; i++ )
					if ( split_on_tab[1][i] == ',' ) number_of_commas++;
				number_of_distinct_cigars_in_a_line = splitByDelimiter (split_on_tab[1] , ',' , split_on_comma);
			}

			/*
			 printf ("\ncurr_position %d\n" , curr_position);
			 */
			strcpy(sam_alignment->reference_name , chromosome);
			for ( k = 0 ; k < number_of_distinct_cigars_in_a_line ; k++ )
			{
				splitByDelimiter (split_on_comma[k] , '-' , split_on_dash);
				if ( flag_save_scores == 0 )
					number_of_repititions_of_the_same_reads = strtol (split_on_dash[1] , &temp , 10);
				else
				{
					if ( ! ( split_on_comma[j][1] == '-' && isalpha (split_on_dash[0][0]) != 0 ) )
					{
						sam_alignment->mapping_quality_score = strtol (split_on_dash[1] , &temp , 10);
						strcpy(sam_alignment->tags[3].val , split_on_dash[2]);
						number_of_repititions_of_the_same_reads = strtol (split_on_dash[3] , &temp , 10);
					}
				}

				sam_alignment->start_position = curr_position;
				if ( curr_position > end || curr_position < start ) continue;
				if ( split_on_comma[k][1] == '-' && isalpha (split_on_dash[0][0]) != 0 )
				{
					// Use the same cigar
					sprintf(temp , "%d" , read_number);
					read_number++;
					strcpy(sam_alignment->read_name , temp);
				}
				else
				{
					sprintf(temp , "%d" , read_number);
					read_number++;
					strcpy(sam_alignment->read_name , temp);
					strcpy(sam_alignment->icigar , split_on_dash[0]);
					//printSamAlignmentInstance (sam_alignment , 1);
					//fflush (stdout);
					convertIcigarToCigarandMDSingleEnded (single_genome_sequence , sam_alignment , chromosome , flag_ignore_mismatches , flag_ignore_soft_clippings , flag_ignore_unmapped_sequences , flag_ignore_quality_score , flag_ignore_sequence_information , default_quality_value);
					//printf ("\nconversion complete");
					//fflush (stdout);
					sprintf(temp , "%d" , read_number);
					read_number++;
					strcpy(sam_alignment->read_name , temp);
				}

				/*
				 * Write the alignments to stdout
				 */
				while ( number_of_repititions_of_the_same_reads-- )
				{
					printf ("%s" , read_prefix);
					printf ("%s_" , sam_alignment->read_name);
					printf ("%d" , number_of_repititions_of_the_same_reads);
					printf ("\t");

					printf ("%d" , sam_alignment->samflag);
					printf ("\t");

					printf ("%s" , sam_alignment->reference_name);
					printf ("\t");

					printf ("%d" , sam_alignment->start_position);
					printf ("\t");

					if ( flag_save_scores == 0 )
						printf ("255");
					else printf ("%d" , sam_alignment->mapping_quality_score);
					printf ("\t");

					printf ("%s" , sam_alignment->cigar);
					printf ("\t");

					printf ("*");
					printf ("\t");

					printf ("0");
					printf ("\t");

					printf ("0");
					printf ("\t");

					printf ("%s" , sam_alignment->seq);
					printf ("\t");

					if ( qual_pool_iterator >= total_quals )
					{
						printf ("\nQUAL SCORE ERROR");
						exit (2);
					}
					if ( flag_save_all_quality_scores == 1 )
					{
						/*
						 printf ("\nqual_pool_iterator %d total_quals %d\n " , qual_pool_iterator , total_quals);
						 fflush (stdout);
						 */
						printf ("%s" , split_on_newline_qual[qual_pool_iterator++ ]);
					}
					else printf ("%s" , sam_alignment->qual);
					printf ("\t");
					/*
					 printf ("\nWritten upto here");
					 fflush (stdout);
					 printf ("NH:i:%s" , sam_alignment->tags[0].val);
					 printf ("\t");
					 */
					if ( strcmp (sam_alignment->tags[1].val , ".") != 0 )
					{
						printf ("XS:A:%s" , sam_alignment->tags[1].val);
						printf ("\t");
					}
					printf ("MD:Z:%s" , sam_alignment->tags[2].val);
					if ( flag_save_scores == 1 )
					{
						printf ("\t");
						printf ("AS:i:%s" , sam_alignment->tags[3].val);
					}

					printf ("\n");
					fflush (stdout);
				}
				/*
				 printf ("\nFile write complete");
				 fflush (stdout);
				 */
			}
		}
	}
	return 0;
}
