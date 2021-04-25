# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include <math.h>
# include <stdbool.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

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

void convertReadIdToString (int *read_id, char *read_id_string, int read_length, char *alphabets)
{
	int i;
	for ( i = 0 ; i < read_length ; i++ )
		read_id_string[i] = alphabets[read_id[i]];
	read_id_string[i] = '\0';
}

void compressPairedEndedAlignments (char *name_of_file_with_quality_scores, char *name_of_file_with_max_commas, char *input_samfilename, char *output_abridgefilename, char *unmapped_filename, char *genome_filename, short int flag_ignore_soft_clippings, short int flag_ignore_mismatches, short int flag_ignore_unmapped_sequences, short int flag_ignore_quality_score, short int run_diagnostics, long long int max_input_reads_in_a_single_nucl_loc, short int flag_save_all_quality_scores, short int flag_save_exact_quality_scores, long long int max_number_of_alignments, int max_read_length)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	FILE *fhr;
	FILE *fhw_pass1;
	FILE *fhw_unmapped;
	FILE *fhw_name_of_file_with_max_commas;
	FILE *fhw_qual;

	char **qual_scores;
	char **split_line; // List of strings to store each element of a single alignment
	char **split_tags; // List of strings to store tag information
	char **split_reference_info;
	char alphabets[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-='{}[]|?<>,.";
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
	strcat(temp , "\n");
	fprintf (fhw_pass1 , "%s" , temp);

	/*
	 * For diagnostics
	 */
	if ( run_diagnostics == 1 )
		readInTheEntireGenome (genome_filename , whole_genome);

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
		if ( curr_alignment->samflag == 4 )
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
		//generateIntegratedCigar (curr_alignment , flag_ignore_soft_clippings , flag_ignore_mismatches , flag_ignore_unmapped_sequences , flag_ignore_quality_score , whole_genome , sam_alignment_instance_diagnostics , number_of_records_read , run_diagnostics);
		NH_tag_index = -1;
		for ( i = 0 ; i < curr_alignment->number_of_tag_items ; i++ )
			if ( strcmp (curr_alignment->tags[i].name , "NH") == 0 )
				NH_tag_index = i;

	} while ( ( line_len = getline ( &line , &len , fhr) ) != -1 );

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
	char *temp; //Required for strtoi

	short int flag_ignore_soft_clippings;
	short int flag_ignore_mismatches;
	short int flag_ignore_quality_score;
	short int flag_ignore_unmapped_sequences;
	short int run_diagnostics;
	short int save_all_quality_scores;
	short int save_exact_quality_scores;

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

	/********************************************************************/

	compressPairedEndedAlignments (name_of_file_with_quality_scores , name_of_file_with_max_commas , input_samfilename , output_abridgefilename , unmapped_filename , genome_filename , flag_ignore_soft_clippings , flag_ignore_mismatches , flag_ignore_unmapped_sequences , flag_ignore_quality_score , run_diagnostics , max_input_reads_in_a_single_nucl_loc , save_all_quality_scores , save_exact_quality_scores , max_number_of_alignments , max_read_length);
	return 0;
}
