# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include <math.h>
# include <stdbool.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

char findMatchCharacterIcigar (
		char *icigar,
		char samformatflag_replacer_characters[])
{
	int i;

	for ( i = 0 ; icigar[i] != '\0' ; i++ )
	{
		if ( strchr (samformatflag_replacer_characters , icigar[i]) != NULL )
			return icigar[i];
	}
	return ' ';
}

void generateNextReadID (char *alphabets, int *read_id, int *read_length)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	int i;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/

	/********************************************************************/

	if ( *read_length == 0 )
	{
		read_id[0] = 0;
		( *read_length )++;
	}
	else
	{
		/*
		 * Check if the read is the last element of the maximum read_length
		 */
		for ( i = 0 ; i < *read_length ; i++ )
			if ( read_id[i] != strlen (alphabets) - 1 ) break;
		if ( i == *read_length )
		{
			( *read_length )++;
			for ( i = 0 ; i < *read_length ; i++ )
				read_id[i] = 0;
		}
		else
		{
			/*
			 * Increment the read_id
			 */
			for ( i = *read_length - 1 ; i >= 0 ; i-- )
			{
				if ( read_id[i] == ( strlen (alphabets) - 1 ) )
					read_id[i] = 0;
				else
				{
					read_id[i]++;
					break;
				}
			}
		}
	}
}

void convertReadIdToString (
		int *read_id,
		char *read_id_string,
		int read_length,
		char *alphabets)
{
	int i;
	for ( i = 0 ; i < read_length ; i++ )
		read_id_string[i] = alphabets[read_id[i]];
	read_id_string[i] = '\0';
}

char returnDirection (
		char character,
		struct Paired_Ended_Flag_to_Single_Character *samflag_dictionary,
		int number_of_unique_samformatflags)
{
	int i;
	for ( i = 0 ; i < number_of_unique_samformatflags * 2 ; i++ )
	{
		if ( samflag_dictionary->character[i] == character )
			return samflag_dictionary->direction[i];
	}
	return ' ';
}

void writeToFile (
		short int flag_save_all_quality_scores,
		short int flag_save_exact_quality_scores,
		FILE *fhw_qual,
		FILE *fhw_pass1,
		struct Compressed_DS **compressed_ds_pool,
		int compressed_ds_pool_total,
		char *write_to_file_col1,
		char *write_to_file_col2,
		char *write_to_file_col3,
		char *encoded_string,
		long long int *count,
		char **qual_Scores,
		int quality_score_index,
		char samformatflag_replacer_characters[],
		int number_of_unique_samformatflags,
		struct Paired_Ended_Flag_to_Single_Character *samflag_dictionary,
		short int flag_ignore_soft_clippings,
		struct Cigar_Items *cigar_items_instance)
{
	int i, j, k, l;
	char str[1000];
	char qual[MAX_SEQ_LEN];
	char line_to_be_written_to_file[MAX_LINE_TO_BE_WRITTEN_TO_FILE];
	char list_of_read_names[MAX_LINE_TO_BE_WRITTEN_TO_FILE];
	int num_of_types;

	line_to_be_written_to_file[0] = '\0';
	list_of_read_names[0] = '\0';
	for ( i = 0 ; i < compressed_ds_pool_total ; i++ )
	{
		if ( i == 0 )
		{
			if ( compressed_ds_pool[i]->position != 1 )
			{
				sprintf(str , "%lld" , compressed_ds_pool[i]->position);
				strcat(line_to_be_written_to_file , str);
				strcat(line_to_be_written_to_file , "\t");
			}
			else str[0] = '\0'; // empty string
		}
		strcat(line_to_be_written_to_file , compressed_ds_pool[i]->icigar);
		strcat(line_to_be_written_to_file , "-");
		sprintf(str , "%ld" , compressed_ds_pool[i]->num_reads);
		strcat(line_to_be_written_to_file , str);
		if ( i != compressed_ds_pool_total - 1 )
		strcat(line_to_be_written_to_file , ",");

		for ( j = 0 ; j < compressed_ds_pool[i]->num_reads ; j++ )
		{
			strcat(list_of_read_names , compressed_ds_pool[i]->pointers_to_read_names[j]);
			if ( i != compressed_ds_pool_total - 1 || j != compressed_ds_pool[i]->num_reads - 1 )
			strcat(list_of_read_names , ",");
		}

		if ( flag_save_all_quality_scores == 1 )
		{
			for ( j = 0 ; j < compressed_ds_pool[i]->num_reads ; j++ )
			{
				if ( returnDirection (findMatchCharacterIcigar (compressed_ds_pool[i]->icigar , samformatflag_replacer_characters) , samflag_dictionary , number_of_unique_samformatflags) == '+' )
				{
					for ( k = 0 ;
							compressed_ds_pool[i]->pointers_to_qual_scores[j][k] != '\0' ;
							k++ )
						qual[k] = compressed_ds_pool[i]->pointers_to_qual_scores[j][k] - 90;
					qual[k] = '\0';
				}
				else if ( returnDirection (findMatchCharacterIcigar (compressed_ds_pool[i]->icigar , samformatflag_replacer_characters) , samflag_dictionary , number_of_unique_samformatflags) == '-' )
				{
					for ( k = strlen (compressed_ds_pool[i]->pointers_to_qual_scores[j]) - 1 ;
							k >= 0 ; k-- )
						qual[strlen (compressed_ds_pool[i]->pointers_to_qual_scores[j]) - 1 - k] = compressed_ds_pool[i]->pointers_to_qual_scores[j][k] - 90;
					qual[strlen (compressed_ds_pool[i]->pointers_to_qual_scores[j])] = '\0';
				}
				fprintf (fhw_qual , "%s" , qual);
				if ( flag_save_exact_quality_scores == 0 )
				{
					fprintf (fhw_qual , "%s" , "\t");
					if ( flag_ignore_soft_clippings == 1 )
					{
						splitCigar (compressed_ds_pool[i]->cigar , &num_of_types , cigar_items_instance);
						if ( cigar_items_instance[0].def == 'S' ) // Left soft clip exists
						{
							sprintf(str , "%ld" , cigar_items_instance[0].len);
							fprintf (fhw_qual , "%s" , str);
							fprintf (fhw_qual , "%s" , "S");
						}
					}
					if ( compressed_ds_pool[i]->icigar[1] != '\0' )
					{
						for ( k = 0 ;
								compressed_ds_pool[i]->icigar[k + 1] != '~' ;
								k++ )
							fputc (compressed_ds_pool[i]->icigar[k] , fhw_qual);
					}
					else
					{
						/*
						 * Find the icigar before the ith one which is not compressed
						 */
						l = i - 1;
						while ( l != 0 )
						{
							if ( compressed_ds_pool[l]->icigar[1] != '\0' )
								break;
							else l--;
						}

						for ( k = 0 ;
								compressed_ds_pool[l]->icigar[k + 1] != '~' ;
								k++ )
							fputc (compressed_ds_pool[l]->icigar[k] , fhw_qual);
					}
					if ( flag_ignore_soft_clippings == 1 )
					{
						//splitCigar (compressed_ds_pool[i]->cigar , &num_of_types , cigar_items_instance);
						if ( cigar_items_instance[num_of_types - 1].def == 'S' ) // Right soft clip exists
						{
							sprintf(str , "%ld" , cigar_items_instance[num_of_types - 1].len);
							fprintf (fhw_qual , "%s" , str);
							fprintf (fhw_qual , "%s" , "S");
						}
					}
					if ( returnDirection (findMatchCharacterIcigar (compressed_ds_pool[i]->icigar , samformatflag_replacer_characters) , samflag_dictionary , number_of_unique_samformatflags) == '-' )
						fprintf (fhw_qual , "%s" , "2");
					else fprintf (fhw_qual , "%s" , "1");
				}
				fprintf (fhw_qual , "%s" , "\n");
			}
		}
	}
	strcat(line_to_be_written_to_file , "\t");
	strcat(line_to_be_written_to_file , list_of_read_names);
	strcat(line_to_be_written_to_file , "\n");
	fprintf (fhw_pass1 , "%s" , line_to_be_written_to_file);
	*count = compressed_ds_pool_total;
}

void prepareIcigarForComparison (
		char *icigar1,
		char *icigar,
		char samformatflag_replacer_characters[])
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	int i, j, k;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	for ( i = 0 ; icigar[i] != '\0' ; i++ )
	{
		if ( strchr (samformatflag_replacer_characters , icigar[i]) == NULL )
			icigar1[i] = icigar[i];
		else icigar1[i] = 'M';
	}
	icigar1[i] = '\0';
}

void reModeliCIGARSPairedEnded (
		struct Compressed_DS **compressed_ds_pool,
		struct Compressed_DS **compressed_ds_pool_rearranged,
		short *already_processed,
		int compressed_ds_pool_index,
		char **modified_icigars,
		char samformatflag_replacer_characters[])
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	int i, j, k;
	int compressed_ds_pool_rearranged_index = 0;

	//char icigar1[MAX_SEQ_LEN];
	//char icigar2[MAX_SEQ_LEN];
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	for ( i = 0 ; i < compressed_ds_pool_index ; i++ )
		already_processed[i] = 0;

	/********************************************************************/
	for ( i = 0 ; i < compressed_ds_pool_index ; i++ )
		prepareIcigarForComparison (modified_icigars[i] , compressed_ds_pool[i]->icigar , samformatflag_replacer_characters);

	for ( i = 0 ; i < compressed_ds_pool_index ; i++ )
	{
		if ( already_processed[i] == 1 ) continue;
		already_processed[i] = 1;
		//prepareIcigarForComparison ( icigar1 , compressed_ds_pool[i]->icigar,samformatflag_replacer_characters );
		/*
		 * Copy the icigar entry into the rearranged pool
		 */
		//icigar1 = modified_icigars[i];
		strcpy(compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->icigar , compressed_ds_pool[i]->icigar);
		strcpy(compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->cigar , compressed_ds_pool[i]->cigar);
		compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->num_reads = compressed_ds_pool[i]->num_reads;
		compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->position = compressed_ds_pool[i]->position;
		for ( k = 0 ; k < compressed_ds_pool[i]->num_reads ; k++ )
		{
			compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->pointers_to_qual_scores[k] = compressed_ds_pool[i]->pointers_to_qual_scores[k];
			compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->pointers_to_read_names[k] = compressed_ds_pool[i]->pointers_to_read_names[k];
		}
		compressed_ds_pool_rearranged_index++;

		for ( j = i + 1 ; j < compressed_ds_pool_index ; j++ )
		{
			if ( already_processed[j] == 1 ) continue;
			//prepareIcigarForComparison ( icigar2 , compressed_ds_pool[j]->icigar );
			//icigar2 = modified_icigars[j];
			if ( strcmp (modified_icigars[i] , modified_icigars[j]) == 0 )
			{
				//printf ( "\n%s %s %d" , modified_icigars[i] , modified_icigars[j] , strcmp ( modified_icigars[i] , modified_icigars[j] ) );
				//printf ( "\n%s %s %d" , compressed_ds_pool[i]->icigar , compressed_ds_pool[j]->icigar , strcmp ( compressed_ds_pool[i]->icigar , compressed_ds_pool[j]->icigar ) );
				//strcpy( compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->icigar , compressed_ds_pool[i]->icigar );
				compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->icigar[0] = findMatchCharacterIcigar (compressed_ds_pool[j]->icigar , samformatflag_replacer_characters);
				compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->icigar[1] = '\0';
				compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->num_reads = compressed_ds_pool[j]->num_reads;
				compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->position = compressed_ds_pool[j]->position;
				for ( k = 0 ; k < compressed_ds_pool[j]->num_reads ; k++ )
				{
					compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->pointers_to_qual_scores[k] = compressed_ds_pool[j]->pointers_to_qual_scores[k];
					compressed_ds_pool_rearranged[compressed_ds_pool_rearranged_index]->pointers_to_read_names[k] = compressed_ds_pool[j]->pointers_to_read_names[k];
				}
				compressed_ds_pool_rearranged_index++;
				already_processed[j] = 1;
			}
		}
	}
	/*
	 printf ( "\n %d %d" , compressed_ds_pool_rearranged_index , compressed_ds_pool_index );
	 if ( compressed_ds_pool_index > 25000 )
	 {
	 for ( i = 0 ; i < compressed_ds_pool_index ; i++ )
	 {
	 printf ( "\n%s %d %s %d" , compressed_ds_pool[i]->icigar , compressed_ds_pool[i]->num_reads , compressed_ds_pool_rearranged[i]->icigar , compressed_ds_pool_rearranged[i]->num_reads );
	 }
	 printf ( "\n==============================================================================================================================================================================================" );
	 }
	 */
}

void compressPairedEndedAlignments (
		char *frequency_of_flags_filename,
		char *name_of_file_with_quality_scores,
		char *name_of_file_with_max_commas,
		char *input_samfilename,
		char *output_abridgefilename,
		char *unmapped_filename,
		char *genome_filename,
		short int flag_ignore_soft_clippings,
		short int flag_ignore_mismatches,
		short int flag_ignore_unmapped_sequences,
		short int flag_ignore_quality_score,
		short int run_diagnostics,
		long long int max_input_reads_in_a_single_nucl_loc,
		short int flag_save_all_quality_scores,
		short int flag_save_exact_quality_scores,
		long long int max_number_of_alignments,
		int max_read_length,
		char *dictionary_filename,
		short int flag_save_scores)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	FILE *fhr;
	FILE *fhw_pass1;
	FILE *fhw_unmapped;
	FILE *fhw_name_of_file_with_max_commas;
	FILE *fhw_qual;
	FILE *fhw_dictionary;
	FILE *fhr_freq_samflags;

	char **qual_scores;
	char **read_names;
	char **split_line; // List of strings to store each element of a single alignment
	char **split_tags; // List of strings to store tag information
	char **split_reference_info;
	char alphabets[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-='{}[]|?<>,.";
	char samformatflag_replacer_characters[] = "BEFHIJKLMOPQRSUVWXYZbdefhijklmopqrsuvwxyz";
	char *temp; //Useless
	char *line = NULL; // for reading each line
	char *entry_in_output_file; //entry in output file
	char *prev_reference_name;
	char *curr_reference_name;
	char *reference_id_quick_read;
	char *samflag_quick_read;
	char *write_to_file_col1;
	char *write_to_file_col2;
	char *write_to_file_col3;
	char *encoded_string;
	char **modified_icigars;
	char str[100];
	char read_id_string[100];
	char line_to_be_written_to_file[1000];

	size_t len = 0;
	ssize_t line_len;

	short *already_processed;
	short int NH_tag_index;

	int flag;
	int i, j, k; // Required in loops
	int number_of_tags;
	int number_of_fields; // Number of fields in each sam alignment entry
	int sam_tag_index;
	int tab_number;
	int num_items_in_alignment_pool = 0; // Items in pool
	int samflag_quick_read_index = 0;
	int compressed_ds_pool_index = 0;
	int quality_score_index = 0;
	int number_of_reference_sequences = 0;
	int reference_sequence_index = 0;
	int number_of_repetitions = 0;
	int read_length;
	int read_id[100];
	int number_of_unique_samformatflags;
	int current_length_of_line_to_be_written_to_file;

	long long int relative_position_to_previous_read_cluster;
	long long int previous_position = -1;
	long long int current_position;
	long long int number_of_records_written = 0;
	long long int number_of_records_read = 0;
	long long int num_pools_written = 0;
	long long int max_commas = 0;
	long long int curr_commas = 0;
	unsigned long long int mega_array_index;

	struct Sam_Alignment *prev_alignment;
	struct Sam_Alignment *curr_alignment;
	struct Sam_Alignment *sam_alignment_instance_diagnostics;
	struct Sam_Alignment *temp_alignment;
	struct Sam_Alignment **alignment_pool_same_position;
	struct Compressed_DS **compressed_ds_pool;
	struct Compressed_DS **compressed_ds_pool_rearranged;
	struct Reference_Sequence_Info **reference_info;
	struct Whole_Genome_Sequence *whole_genome;
	struct Paired_Ended_Flag_to_Single_Character *samflag_dictionary;
	struct Cigar_Items *cigar_items_instance;
	cigar_items_instance = ( struct Cigar_Items* ) malloc (sizeof(struct Cigar_Items) * 50);
	//printf ("\nmax_number_of_alignments %d" , max_number_of_alignments);

	//printf ("\nSize of mega array %d" , sizeof ( mega_array ));

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr = fopen (input_samfilename , "r");
	if ( fhr == NULL )
	{
		printf ("Error! File %s not found" , input_samfilename);
		exit (1);
	}
	fhr_freq_samflags = fopen (frequency_of_flags_filename , "r");
	if ( fhr_freq_samflags == NULL )
	{
		printf ("Error! File %s not found" , fhr_freq_samflags);
		exit (1);
	}
	fhw_pass1 = fopen (output_abridgefilename , "w");
	if ( fhw_pass1 == NULL )
	{
		printf ("%s File cannot be created" , output_abridgefilename);
		exit (1);
	}
	fhw_unmapped = fopen (unmapped_filename , "w");
	if ( fhw_unmapped == NULL )
	{
		printf ("%s File cannot be created" , unmapped_filename);
		exit (1);
	}
	fhw_name_of_file_with_max_commas = fopen (name_of_file_with_max_commas , "w");
	if ( fhw_name_of_file_with_max_commas == NULL )
	{
		printf ("%s File cannot be created" , name_of_file_with_max_commas);
		exit (1);
	}
	fhw_qual = fopen (name_of_file_with_quality_scores , "w");
	if ( fhw_qual == NULL )
	{
		printf ("%s File cannot be created" , name_of_file_with_quality_scores);
		exit (1);
	}
	fhw_dictionary = fopen (dictionary_filename , "w");
	if ( fhw_dictionary == NULL )
	{
		printf ("%s File cannot be created" , dictionary_filename);
		exit (1);
	}

	split_line = ( char** ) malloc (sizeof(char*) * ROWS);
	for ( i = 0 ; i < ROWS ; i++ )
		split_line[i] = ( char* ) malloc (sizeof(char) * COLS);

	split_tags = ( char** ) malloc (sizeof(char*) * ROWS);
	for ( i = 0 ; i < ROWS ; i++ )
		split_tags[i] = ( char* ) malloc (sizeof(char) * COLS);

	split_reference_info = ( char** ) malloc (sizeof(char*) * ROWS);
	for ( i = 0 ; i < ROWS ; i++ )
		split_reference_info[i] = ( char* ) malloc (sizeof(char) * COLS);

	already_processed = ( short* ) malloc (sizeof(short) * max_input_reads_in_a_single_nucl_loc);
	max_input_reads_in_a_single_nucl_loc += 5;
	compressed_ds_pool = ( struct Compressed_DS** ) malloc (sizeof(struct Compressed_DS*) * max_input_reads_in_a_single_nucl_loc);
	for ( i = 0 ; i < max_input_reads_in_a_single_nucl_loc ; i++ )
		compressed_ds_pool[i] = allocateMemoryCompressed_DS (max_input_reads_in_a_single_nucl_loc);

	compressed_ds_pool_rearranged = ( struct Compressed_DS** ) malloc (sizeof(struct Compressed_DS*) * max_input_reads_in_a_single_nucl_loc);
	for ( i = 0 ; i < max_input_reads_in_a_single_nucl_loc ; i++ )
		compressed_ds_pool_rearranged[i] = allocateMemoryCompressed_DS (max_input_reads_in_a_single_nucl_loc);

	write_to_file_col1 = ( char* ) malloc (sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE);
	write_to_file_col2 = ( char* ) malloc (sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE);
	write_to_file_col3 = ( char* ) malloc (sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE);
	encoded_string = ( char* ) malloc (sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE);
	write_to_file_col1[0] = '\0';
	write_to_file_col2[0] = '\0';
	write_to_file_col3[0] = '\0';

	reference_id_quick_read = ( char* ) malloc (sizeof(char) * 1000);
	samflag_quick_read = ( char* ) malloc (sizeof(char) * 1000);
	prev_reference_name = ( char* ) malloc (sizeof(char) * 1000);
	prev_reference_name[0] = '\0';
	curr_reference_name = ( char* ) malloc (sizeof(char) * 1000);
	curr_reference_name[0] = '\0';

	curr_alignment = allocateMemorySam_Alignment ();
	prev_alignment = allocateMemorySam_Alignment ();
	temp_alignment = allocateMemorySam_Alignment ();
	sam_alignment_instance_diagnostics = allocateMemorySam_Alignment ();
	reference_info = ( struct Reference_Sequence_Info** ) malloc (sizeof(struct Reference_Sequence_Info*) * MAX_REFERENCE_SEQUENCES);
	for ( i = 0 ; i < MAX_REFERENCE_SEQUENCES ; i++ )
		reference_info[i] = allocateMemoryReference_Sequence_Info ();

	temp = ( char* ) malloc (sizeof(char) * MAX_GENERAL_LEN);
	whole_genome = ( struct Whole_Genome_Sequence* ) malloc (sizeof(struct Whole_Genome_Sequence));
	qual_scores = ( char** ) malloc (sizeof(char*) * max_input_reads_in_a_single_nucl_loc);
	for ( i = 0 ; i < max_input_reads_in_a_single_nucl_loc ; i++ )
		qual_scores[i] = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	read_names = ( char** ) malloc (sizeof(char*) * max_input_reads_in_a_single_nucl_loc);
	for ( i = 0 ; i < max_input_reads_in_a_single_nucl_loc ; i++ )
		read_names[i] = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	modified_icigars = ( char** ) malloc (sizeof(char*) * max_input_reads_in_a_single_nucl_loc);
	for ( i = 0 ; i < max_input_reads_in_a_single_nucl_loc ; i++ )
		modified_icigars[i] = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	/********************************************************************/

	/*
	 * Write the first line in output file
	 */
	temp[0] = '\0';
	sprintf(str , "%lld" , flag_ignore_mismatches);
	strcat(temp , str);
	strcat(temp , "\t");
	sprintf(str , "%lld" , flag_ignore_soft_clippings);
	strcat(temp , str);
	strcat(temp , "\t");
	sprintf(str , "%lld" , flag_ignore_unmapped_sequences);
	strcat(temp , str);
	strcat(temp , "\t");
	sprintf(str , "%lld" , flag_ignore_quality_score);
	strcat(temp , str);
	strcat(temp , "\t");
	sprintf(str , "%lld" , flag_save_all_quality_scores);
	strcat(temp , str);
	strcat(temp , "\t");
	sprintf(str , "%lld" , flag_save_exact_quality_scores);
	strcat(temp , str);
	strcat(temp , "\t");
	sprintf(str , "%lld" , flag_save_scores);
	strcat(temp , str);
	strcat(temp , "\n");
	fprintf (fhw_pass1 , "%s" , temp);

	/*
	 * For diagnostics
	 */
	if ( run_diagnostics == 1 )
		readInTheEntireGenome (genome_filename , whole_genome);
	/*
	 *	Construct the mapping between samformatflag and the Single Character
	 */
	number_of_unique_samformatflags = 0;
	while ( ( line_len = getline ( &line , &len , fhr_freq_samflags) ) != -1 )
		number_of_unique_samformatflags++;
	samflag_dictionary = allocateMemoryPaired_Ended_Flag_to_Single_Character (number_of_unique_samformatflags);
	rewind (fhr_freq_samflags);
	i = 0;
	j = 0;
	while ( ( line_len = getline ( &line , &len , fhr_freq_samflags) ) != -1 )
	{
		if ( strtol (line , &temp , 10) == 77 || strtol (line , &temp , 10) == 141 )
		{
			number_of_unique_samformatflags--;
			continue;
		}
		samflag_dictionary->direction[i] = '+';
		samflag_dictionary->samflags[i] = strtol (line , &temp , 10);
		samflag_dictionary->character[i] = samformatflag_replacer_characters[j++ ];

		samflag_dictionary->direction[i + 1] = '-';
		samflag_dictionary->character[i + 1] = samformatflag_replacer_characters[j++ ];
		samflag_dictionary->samflags[i + 1] = strtol (line , &temp , 10);
		i += 2;
	}

	for ( i = 0 ; i < number_of_unique_samformatflags * 2 ; i++ )
	{
		//printf ("\n%c %d %c" , samflag_dictionary->direction[i] , samflag_dictionary->samflags[i] , samflag_dictionary->character[i]);
		line_to_be_written_to_file[0] = '\0';
		sprintf(temp , "%lld" , samflag_dictionary->samflags[i]);
		strcat(line_to_be_written_to_file , temp);
		strcat(line_to_be_written_to_file , "\t");
		current_length_of_line_to_be_written_to_file = strlen (line_to_be_written_to_file);
		line_to_be_written_to_file[current_length_of_line_to_be_written_to_file] = samflag_dictionary->direction[i];
		line_to_be_written_to_file[current_length_of_line_to_be_written_to_file + 1] = '\t';
		line_to_be_written_to_file[current_length_of_line_to_be_written_to_file + 2] = samflag_dictionary->character[i];
		line_to_be_written_to_file[current_length_of_line_to_be_written_to_file + 3] = '\n';
		line_to_be_written_to_file[current_length_of_line_to_be_written_to_file + 4] = '\0';
		fprintf (fhw_dictionary , "%s" , line_to_be_written_to_file);
	}
	/*
	 * Read in the reference sequence information
	 */
	while ( ( line_len = getline ( &line , &len , fhr) ) != -1 )
	{
		if ( line[0] == '@' )
		{
			if ( line[1] == 'S' && line[2] == 'Q' )
			{
				//printf("\n Reference: %s %d", line, strlen(line));
				//fflush(stdout);
				strcpy(reference_info[number_of_reference_sequences]->line , line);
				number_of_reference_sequences++;
			}
		}
		else break;
	}
	do
	{
		number_of_fields = splitByDelimiter (line , '\t' , split_line);
		populateSamAlignmentInstance (curr_alignment , split_line , number_of_fields , split_tags);
		strcpy(curr_reference_name , curr_alignment->reference_name);
		if ( curr_alignment->samflag == 77 || curr_alignment->samflag == 141 )
		{
			if ( flag_ignore_unmapped_sequences == 0 )
			{
				//Write the unmapped reads into file
				fprintf (fhw_unmapped , "%s" , curr_alignment->seq);
				fprintf (fhw_unmapped , "%s" , "\n");
				for ( i = 0 ; curr_alignment->qual[i] != '\0' ; i++ )
					curr_alignment->qual[i] -= 90;
				fprintf (fhw_unmapped , "%s" , curr_alignment->qual);
				fprintf (fhw_unmapped , "%s" , "\n");
			}
			continue;
		}
		current_position = curr_alignment->start_position;
		/*
		 * Change this function
		 */

		generateIntegratedCigarPairedEnded (curr_alignment , flag_save_scores , flag_ignore_soft_clippings , flag_ignore_mismatches , flag_ignore_unmapped_sequences , flag_ignore_quality_score , whole_genome , sam_alignment_instance_diagnostics , number_of_records_read , run_diagnostics , samflag_dictionary , number_of_unique_samformatflags , samformatflag_replacer_characters);
		NH_tag_index = -1;
		for ( i = 0 ; i < curr_alignment->number_of_tag_items ; i++ )
			if ( strcmp (curr_alignment->tags[i].name , "NH") == 0 )
				NH_tag_index = i;

		if ( strlen (prev_reference_name) == 0 ) // 1st chromosome - initialize stuffs
		{
			//printf("\n1. compressed_ds_pool_index %d", compressed_ds_pool_index);
			//fflush(stdout);
			previous_position = current_position;
			strcpy(prev_reference_name , curr_reference_name);
			strcpy(qual_scores[quality_score_index] , curr_alignment->qual);
			strcpy(read_names[quality_score_index] , curr_alignment->read_name);
			strcpy(compressed_ds_pool[compressed_ds_pool_index]->icigar , curr_alignment->icigar);
			strcpy(compressed_ds_pool[compressed_ds_pool_index]->cigar , curr_alignment->cigar);
			compressed_ds_pool[compressed_ds_pool_index]->num_reads = 1;
			compressed_ds_pool[compressed_ds_pool_index]->pointers_to_qual_scores[compressed_ds_pool[compressed_ds_pool_index]->num_reads - 1] = qual_scores[quality_score_index];
			compressed_ds_pool[compressed_ds_pool_index]->pointers_to_read_names[compressed_ds_pool[compressed_ds_pool_index]->num_reads - 1] = read_names[quality_score_index];
			compressed_ds_pool[compressed_ds_pool_index]->position = curr_alignment->start_position;
			//printf("\n1. Max_read_at_a_position %d chromosome %s position %d compressed_ds_pool_index %d", compressed_ds_pool[compressed_ds_pool_index]->num_reads, curr_alignment->reference_name, curr_alignment->start_position, compressed_ds_pool_index);
			quality_score_index++;
			compressed_ds_pool_index++;
			//printf("\n Writing Reference to file %s %d", reference_info[reference_sequence_index]->line, reference_sequence_index);
			//fflush(stdout);
			fprintf (fhw_pass1 , "%s" , reference_info[reference_sequence_index]->line);
			reference_sequence_index++;
		}
		else if ( strcmp (prev_reference_name , curr_reference_name) != 0 ) // New chromosome
		{
			//printf("\2. ncompressed_ds_pool_index %d", compressed_ds_pool_index);
			//fflush(stdout);
			reModeliCIGARSPairedEnded (compressed_ds_pool , compressed_ds_pool_rearranged , already_processed , compressed_ds_pool_index , modified_icigars , samformatflag_replacer_characters);
			writeToFile (flag_save_all_quality_scores , flag_save_exact_quality_scores , fhw_qual , fhw_pass1 , compressed_ds_pool_rearranged , compressed_ds_pool_index , write_to_file_col1 , write_to_file_col2 , write_to_file_col3 , encoded_string , &curr_commas , qual_scores , quality_score_index , samformatflag_replacer_characters , number_of_unique_samformatflags , samflag_dictionary , flag_ignore_soft_clippings , cigar_items_instance);
			if ( max_commas < curr_commas ) max_commas = curr_commas;
			//printf ( "\n%lld %lld" , curr_commas , max_commas );
			compressed_ds_pool_index = 0;
			previous_position = current_position;
			strcpy(prev_reference_name , curr_reference_name);
			strcpy(compressed_ds_pool[compressed_ds_pool_index]->icigar , curr_alignment->icigar);
			strcpy(compressed_ds_pool[compressed_ds_pool_index]->cigar , curr_alignment->cigar);
			compressed_ds_pool[compressed_ds_pool_index]->num_reads = 1;
			compressed_ds_pool[compressed_ds_pool_index]->position = curr_alignment->start_position;
			strcpy(qual_scores[quality_score_index] , curr_alignment->qual);
			strcpy(read_names[quality_score_index] , curr_alignment->read_name);
			quality_score_index++;
			//printf("\n2. Max_read_at_a_position %d chromosome %s position %d compressed_ds_pool_index %d", compressed_ds_pool[compressed_ds_pool_index]->num_reads, curr_alignment->reference_name, curr_alignment->start_position, compressed_ds_pool_index);
			compressed_ds_pool_index++;
			//printf("\n Writing Reference to file %s %d", reference_info[reference_sequence_index]->line, reference_sequence_index);
			//fflush(stdout);
			fprintf (fhw_pass1 , "%s" , reference_info[reference_sequence_index]->line);
			reference_sequence_index++;
		}
		else // Same chromosome
		{
			if ( previous_position == current_position )
			{
				//printf("\n3. compressed_ds_pool_index %d", compressed_ds_pool_index);
				//fflush(stdout);
				for ( i = 0 ; i < compressed_ds_pool_index ; i++ )
				{
					if ( strcmp (compressed_ds_pool[i]->icigar , curr_alignment->icigar) == 0 )
					{
						compressed_ds_pool[i]->num_reads++;
						strcpy(qual_scores[quality_score_index] , curr_alignment->qual);
						strcpy(read_names[quality_score_index] , curr_alignment->read_name);
						compressed_ds_pool[i]->pointers_to_qual_scores[compressed_ds_pool[i]->num_reads - 1] = qual_scores[quality_score_index];
						compressed_ds_pool[i]->pointers_to_read_names[compressed_ds_pool[i]->num_reads - 1] = read_names[quality_score_index];
						quality_score_index++;
						//printf("\n3. Max_read_at_a_position %d chromosome %s position %d compressed_ds_pool_index %d", compressed_ds_pool[i]->num_reads, curr_alignment->reference_name, curr_alignment->start_position, i);
						break;
					}
				}
				if ( i == compressed_ds_pool_index ) // New icigar encountered
				{
					strcpy(compressed_ds_pool[compressed_ds_pool_index]->icigar , curr_alignment->icigar);
					strcpy(compressed_ds_pool[compressed_ds_pool_index]->cigar , curr_alignment->cigar);
					compressed_ds_pool[compressed_ds_pool_index]->num_reads = 1;
					compressed_ds_pool[compressed_ds_pool_index]->position = compressed_ds_pool[0]->position;
					strcpy(qual_scores[quality_score_index] , curr_alignment->qual);
					strcpy(read_names[quality_score_index] , curr_alignment->read_name);
					compressed_ds_pool[compressed_ds_pool_index]->pointers_to_qual_scores[compressed_ds_pool[compressed_ds_pool_index]->num_reads - 1] = qual_scores[quality_score_index];
					compressed_ds_pool[compressed_ds_pool_index]->pointers_to_read_names[compressed_ds_pool[compressed_ds_pool_index]->num_reads - 1] = read_names[quality_score_index];
					quality_score_index++;
					//printf("\n4. Max_read_at_a_position %d chromosome %s position %d compressed_ds_pool_index %d", compressed_ds_pool[compressed_ds_pool_index]->num_reads, curr_alignment->reference_name, curr_alignment->start_position, compressed_ds_pool_index);
					compressed_ds_pool_index++;
				}
			}
			else
			{
				//printf("\n4. compressed_ds_pool_index %d", compressed_ds_pool_index);
				//fflush(stdout);
				reModeliCIGARSPairedEnded (compressed_ds_pool , compressed_ds_pool_rearranged , already_processed , compressed_ds_pool_index , modified_icigars , samformatflag_replacer_characters);
				writeToFile (flag_save_all_quality_scores , flag_save_exact_quality_scores , fhw_qual , fhw_pass1 , compressed_ds_pool_rearranged , compressed_ds_pool_index , write_to_file_col1 , write_to_file_col2 , write_to_file_col3 , encoded_string , &curr_commas , qual_scores , quality_score_index , samformatflag_replacer_characters , number_of_unique_samformatflags , samflag_dictionary , flag_ignore_soft_clippings , cigar_items_instance);
				if ( max_commas < curr_commas ) max_commas = curr_commas;
				//printf ( "\n%lld %lld" , curr_commas , max_commas );
				compressed_ds_pool_index = 0;
				quality_score_index = 0;
				strcpy(qual_scores[quality_score_index] , curr_alignment->qual);
				strcpy(read_names[quality_score_index] , curr_alignment->read_name);
				strcpy(compressed_ds_pool[compressed_ds_pool_index]->icigar , curr_alignment->icigar);
				strcpy(compressed_ds_pool[compressed_ds_pool_index]->cigar , curr_alignment->cigar);
				compressed_ds_pool[compressed_ds_pool_index]->num_reads = 1;
				compressed_ds_pool[compressed_ds_pool_index]->pointers_to_qual_scores[compressed_ds_pool[compressed_ds_pool_index]->num_reads - 1] = qual_scores[quality_score_index];
				compressed_ds_pool[compressed_ds_pool_index]->pointers_to_read_names[compressed_ds_pool[compressed_ds_pool_index]->num_reads - 1] = read_names[quality_score_index];
				compressed_ds_pool[compressed_ds_pool_index]->position = curr_alignment->start_position - previous_position;
				quality_score_index++;
				//printf("\n5. Max_read_at_a_position %d chromosome %s position %d compressed_ds_pool_index %d", compressed_ds_pool[compressed_ds_pool_index]->num_reads, curr_alignment->reference_name, curr_alignment->start_position, compressed_ds_pool_index);
				compressed_ds_pool_index++;
			}
			previous_position = current_position;
		}
		//printf("\nMax_read_at_a_position %d chromosome %s position %d compressed_ds_pool_index %d", compressed_ds_pool[compressed_ds_pool_index]->num_reads, curr_alignment->reference_name, curr_alignment->start_position, compressed_ds_pool_index);
		reInitializeSamAlignmentInstance (curr_alignment);
	} while ( ( line_len = getline ( &line , &len , fhr) ) != -1 );

	reModeliCIGARSPairedEnded (compressed_ds_pool , compressed_ds_pool_rearranged , already_processed , compressed_ds_pool_index , modified_icigars , samformatflag_replacer_characters);
	writeToFile (flag_save_all_quality_scores , flag_save_exact_quality_scores , fhw_qual , fhw_pass1 , compressed_ds_pool_rearranged , compressed_ds_pool_index , write_to_file_col1 , write_to_file_col2 , write_to_file_col3 , encoded_string , &curr_commas , qual_scores , quality_score_index , samformatflag_replacer_characters , number_of_unique_samformatflags , samflag_dictionary , flag_ignore_soft_clippings , cigar_items_instance);
	if ( max_commas < curr_commas ) max_commas = curr_commas;
	sprintf(temp , "%lld" , max_commas);
	strcat(temp , "\n");
	fprintf (fhw_name_of_file_with_max_commas , "%s" , temp);

	fclose (fhr);
	fclose (fhw_pass1);
	fclose (fhw_unmapped);
	fclose (fhw_name_of_file_with_max_commas);
	fclose (fhw_dictionary);
}

int main (int argc, char *argv[])
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char input_samfilename[FILENAME_LENGTH];
	char output_abridgefilename[FILENAME_LENGTH];
	char genome_filename[FILENAME_LENGTH];
	char unmapped_filename[FILENAME_LENGTH];
	char name_of_file_with_max_commas[FILENAME_LENGTH];
	char name_of_file_with_quality_scores[FILENAME_LENGTH];
	char frequency_of_flags_filename[FILENAME_LENGTH];
	char dictionary_filename[FILENAME_LENGTH];
	char *temp; //Required for strtoi

	short int flag_ignore_soft_clippings;
	short int flag_ignore_mismatches;
	short int flag_ignore_quality_score;
	short int flag_ignore_unmapped_sequences;
	short int run_diagnostics;
	short int save_all_quality_scores;
	short int save_exact_quality_scores;
	short int save_scores;

	long long int max_input_reads_in_a_single_nucl_loc;
	long long int max_number_of_alignments;

	int i;
	int max_read_length;

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	strcpy(genome_filename , argv[1]);
	flag_ignore_soft_clippings = strtol (argv[2] , &temp , 10);
	flag_ignore_mismatches = strtol (argv[3] , &temp , 10);
	flag_ignore_quality_score = strtol (argv[4] , &temp , 10);
	flag_ignore_unmapped_sequences = strtol (argv[5] , &temp , 10);
	strcpy(input_samfilename , argv[6]);
	strcpy(output_abridgefilename , argv[7]);
	strcpy(unmapped_filename , argv[8]);
	run_diagnostics = strtol (argv[9] , &temp , 10);
	max_input_reads_in_a_single_nucl_loc = strtol (argv[10] , &temp , 10);
	strcpy(name_of_file_with_max_commas , argv[11]);
	save_all_quality_scores = strtol (argv[12] , &temp , 10);
	save_exact_quality_scores = strtol (argv[13] , &temp , 10);
	strcpy(name_of_file_with_quality_scores , argv[14]);
	max_number_of_alignments = strtol (argv[15] , &temp , 10);
	max_read_length = strtol (argv[16] , &temp , 10);
	strcpy(frequency_of_flags_filename , argv[17]);
	strcpy(dictionary_filename , argv[18]);
	save_scores = strtol (argv[19] , &temp , 10);
	/********************************************************************/

	compressPairedEndedAlignments (frequency_of_flags_filename , name_of_file_with_quality_scores , name_of_file_with_max_commas , input_samfilename , output_abridgefilename , unmapped_filename , genome_filename , flag_ignore_soft_clippings , flag_ignore_mismatches , flag_ignore_unmapped_sequences , flag_ignore_quality_score , run_diagnostics , max_input_reads_in_a_single_nucl_loc , save_all_quality_scores , save_exact_quality_scores , max_number_of_alignments , max_read_length , dictionary_filename , save_scores);
	return 0;
}
