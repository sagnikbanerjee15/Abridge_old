# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include <math.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

long long int total_mapped_reads = 0;

void decompressFile (char *name_of_file_with_quality_scores, char *abridge_index_filename, char *genome_filename, char *output_sam_filename, char *pass2_filename, char *genome_prefix, char *unmapped_filename, char *default_quality_value, short int flag_ignore_sequence_information)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	FILE *fhr;
	FILE *fhw;
	FILE *fhr_qual;

	struct Abridge_Index *abridge_index;
	struct Sam_Alignment *sam_alignment;

	int i, j;
	int sam_alignment_pool_index;
	int fread_ret_val;
	int fseek_ret_val;
	int line_number;

	size_t len = 0;
	ssize_t line_len;

	short int flag_ignore_mismatches;
	short int flag_ignore_soft_clippings;
	short int flag_ignore_unmapped_sequences;
	short int flag_ignore_quality_score;

	long long int max_cluster_size;
	unsigned long long int line_num = 0;
	unsigned long long int read_number = 1;
	unsigned long long int from = -1;
	unsigned long long int to = -1;
	unsigned long long int number_of_newlines = 0;
	unsigned long long int number_of_commas_in_each_line = 0;
	unsigned long long int max_number_of_commas = 0;
	int number_of_entries_in_cluster;
	int number_of_elements_after_split_on_delimiter;

	char **split_on_newline;
	char **split_on_tab;
	char **split_on_dash;
	char **split_on_comma;
	char **split_on_delimiter;
	char *buffer;
	char **sequence_portions_from_reference;
	char *fasta_file_with_expressed_portions;
	char *cigar;
	char *md;
	char *output_prefix_without_path;
	char line_to_be_written_to_file[MAX_GENERAL_LEN];
	char temp[100];
	char read_prefix[10];

	struct Sam_Alignment **sam_alignment_pool;
	struct Sam_Alignment *sam_alignment_instance;
	struct Whole_Genome_Sequence *whole_genome;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr = fopen (pass2_filename , "rb");
	if ( fhr == NULL )
	{
		printf ("Error! File not found");
		exit (1);
	}
	fhr_qual = fopen (name_of_file_with_quality_scores , "r");
	if ( fhr_qual == NULL )
	{
		printf ("Error! File not found");
		exit (1);
	}
	fhw = fopen (output_sam_filename , "w");
	if ( fhw == NULL )
	{
		printf ("Error! File cannot be opened for writing");
		exit (1);
	}

	split_on_newline = ( char** ) malloc (sizeof(char*) * ROWS * 10);
	for ( i = 0 ; i < ROWS * 10 ; i++ )
		split_on_newline[i] = ( char* ) malloc (sizeof(char) * COLS * 10);

	split_on_tab = ( char** ) malloc (sizeof(char*) * ROWS);
	for ( i = 0 ; i < ROWS ; i++ )
		split_on_tab[i] = ( char* ) malloc (sizeof(char) * COLS);

	split_on_dash = ( char** ) malloc (sizeof(char*) * ROWS);
	for ( i = 0 ; i < ROWS ; i++ )
		split_on_dash[i] = ( char* ) malloc (sizeof(char) * COLS);

	split_on_comma = ( char** ) malloc (sizeof(char*) * ROWS * 10);
	for ( i = 0 ; i < ROWS * 10 ; i++ )
		split_on_comma[i] = ( char* ) malloc (sizeof(char) * COLS * 10);

	split_on_delimiter = ( char** ) malloc (sizeof(char*) * ROWS);
	for ( i = 0 ; i < ROWS ; i++ )
		split_on_delimiter[i] = ( char* ) malloc (sizeof(char) * COLS);

	output_prefix_without_path = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	sequence_portions_from_reference = ( char** ) malloc (sizeof(char*) * MAX_POOL_SIZE);
	fasta_file_with_expressed_portions = ( char* ) malloc (sizeof(char) * FILENAME_LENGTH);

	buffer = ( char* ) malloc (sizeof(char) * MAX_BUFFER_SIZE_FOR_READING_PASS2_FILE);
	abridge_index = allocateMemoryAbridge_Index ();
	sam_alignment = allocateMemorySam_Alignment ();
	whole_genome = ( struct Whole_Genome_Sequence* ) malloc (sizeof(struct Whole_Genome_Sequence));
	sam_alignment_instance = allocateMemorySam_Alignment ();
	read_prefix[0] = '\0'; // Empty string

	/********************************************************************/

	readAbridgeIndex (abridge_index , abridge_index_filename , split_on_newline , &flag_ignore_mismatches , &flag_ignore_soft_clippings , &flag_ignore_unmapped_sequences , &flag_ignore_quality_score);
	readInTheEntireGenome (genome_filename , whole_genome);
	writeSequenceHeaders (fhw , genome_filename);

	for ( i = 0 ; i < abridge_index->number_of_items ; i++ )
	{
		if ( i != 1178 ) continue;

		//printf("\nFile pointer at %ld", ftell(fhr));
		//printf("\nEntry in index file %s %d %d", abridge_index->chromosome[i], abridge_index->start_byte[i], abridge_index->end_byte[i]);
		fseek_ret_val = fseek (fhr , abridge_index->start_byte[i] , SEEK_SET);
		//printf("\nFile pointer moved to %ld", ftell(fhr));
		buffer[0] = '\0';
		fread_ret_val = fread (buffer , 1 , abridge_index->end_byte[i] - abridge_index->start_byte[i] , fhr);
		//printf("\n fread_ret_val %d fseek_ret_val %d diff %lld", fread_ret_val, fseek_ret_val, MAX_BUFFER_SIZE_FOR_READING_PASS2_FILE - fread_ret_val);
		//fflush(stdout);
		buffer[fread_ret_val] = '\0';
		/*
		 printf("\nfread_ret_val %d", fread_ret_val);
		 printf("\nNew Record:%d\n", i);
		 printf("%s", buffer);
		 printf("\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
		 printf("\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
		 fflush(stdout);
		 */
		printf ("\nBuffer content %s" , buffer);
		number_of_newlines = 0;
		max_number_of_commas = 0;
		number_of_commas_in_each_line = 0;
		for ( j = 0 ; buffer[j] != '\0' ; j++ )
		{
			if ( buffer[j] == '\n' )
			{
				number_of_newlines += 1;
				if ( max_number_of_commas < number_of_commas_in_each_line )
					max_number_of_commas = number_of_commas_in_each_line;
				number_of_commas_in_each_line = 0;
			}
			else if ( buffer[j] == ',' )
				number_of_commas_in_each_line += 1;
		}
		if ( number_of_newlines > ROWS * 10 )
		{
			number_of_newlines += 10;
			for ( j = 0 ; j < ROWS * 10 ; j++ )
				free (split_on_newline[j]);
			free (split_on_newline);
			split_on_newline = ( char** ) malloc (sizeof(char*) * number_of_newlines);
			for ( j = 0 ; j < number_of_newlines ; j++ )
				split_on_newline[j] = ( char* ) malloc (sizeof(char) * COLS * 10);
		}
		number_of_entries_in_cluster = splitByDelimiter (buffer , '\n' , split_on_newline);
		if ( i % 1000 == 0 )
		{
			//printf("\n%d record processed, Number of lines in cluster %d", i, number_of_entries_in_cluster);
			//fflush(stdout);
		}
		number_of_entries_in_cluster--; //Last line is always empty
		printf ("\nLine num %d Number of newlines %d Number of commas %d number_of_entries_in_cluster %d" , i , number_of_newlines , max_number_of_commas , number_of_entries_in_cluster);
		fflush (stdout);
		convertToAlignment (sam_alignment_instance , whole_genome , split_on_newline , sam_alignment , i , abridge_index , number_of_entries_in_cluster , split_on_tab , split_on_dash , split_on_comma , default_quality_value , flag_ignore_mismatches , flag_ignore_soft_clippings , flag_ignore_unmapped_sequences , flag_ignore_quality_score , flag_ignore_sequence_information , &read_number , &total_mapped_reads , read_prefix , from , to , fhw , fhr_qual);
		//if (i == 10) break;
	}
	return;
	/*
	 * Write all unmapped reads to samfile
	 */
	fhr = fopen (unmapped_filename , "r");
	if ( fhr == NULL )
	{
		printf ("Error! File not found");
		exit (1);
	}
	line_number = 1;
	free (buffer);
	buffer = NULL;
	while ( ( line_len = getline ( &buffer , &len , fhr) ) != -1 )
	{
		strcpy(sam_alignment->seq , buffer);
		line_len = getline ( &buffer , &len , fhr);
		strcpy(sam_alignment->qual , buffer);

		line_to_be_written_to_file[0] = '\0';
		strcat(line_to_be_written_to_file , "unmapped_");
		sprintf(temp , "%d" , line_number);
		strcat(line_to_be_written_to_file , temp);
		strcat(line_to_be_written_to_file , "\t");

		strcat(line_to_be_written_to_file , "4");
		strcat(line_to_be_written_to_file , "\t");

		strcat(line_to_be_written_to_file , "*");
		strcat(line_to_be_written_to_file , "\t");

		strcat(line_to_be_written_to_file , "0");
		strcat(line_to_be_written_to_file , "\t");

		strcat(line_to_be_written_to_file , "0");

		strcat(line_to_be_written_to_file , "\t");
		strcat(line_to_be_written_to_file , "*");

		strcat(line_to_be_written_to_file , "\t");
		strcat(line_to_be_written_to_file , "*");

		strcat(line_to_be_written_to_file , "\t");
		strcat(line_to_be_written_to_file , "0");

		strcat(line_to_be_written_to_file , "\t");
		strcat(line_to_be_written_to_file , "0");

		strcat(line_to_be_written_to_file , "\t");
		sam_alignment->seq[strlen (sam_alignment->seq) - 1] = '\0';
		strcat(line_to_be_written_to_file , sam_alignment->seq);

		strcat(line_to_be_written_to_file , "\t");
		sam_alignment->qual[strlen (sam_alignment->qual) - 1] = '\0';
		strcat(line_to_be_written_to_file , sam_alignment->qual);

		strcat(line_to_be_written_to_file , "\tNH:i:0\tHI:i:0\tnM:i:1\tuT:A:1");

		strcat(line_to_be_written_to_file , "\n");
		fprintf (fhw , "%s" , line_to_be_written_to_file);

		line_number++;
	}

	fclose (fhw);
	fclose (fhr);
}

int main (int argc, char *argv[])
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char abridge_index_filename[FILENAME_LENGTH];
	char genome_filename[FILENAME_LENGTH];
	char pass2_filename[FILENAME_LENGTH];
	char output_sam_filename[FILENAME_LENGTH];
	char genome_prefix[FILENAME_LENGTH];
	char name_of_file_with_quality_scores[FILENAME_LENGTH];
	char default_quality_value[10];
	char unmapped_filename[FILENAME_LENGTH];
	char *temp; //Required for strtoi

	short int flag_ignore_sequence_information;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	strcpy(abridge_index_filename , argv[1]);
	strcpy(genome_filename , argv[2]);
	strcpy(output_sam_filename , argv[3]);
	strcpy(pass2_filename , argv[4]);
	strcpy(genome_prefix , argv[5]);
	strcpy(default_quality_value , argv[6]);
	flag_ignore_sequence_information = strtol (argv[7] , &temp , 10);
	strcpy(unmapped_filename , argv[8]);
	strcpy(name_of_file_with_quality_scores , argv[9]);
	/********************************************************************/

	decompressFile (name_of_file_with_quality_scores , abridge_index_filename , genome_filename , output_sam_filename , pass2_filename , genome_prefix , unmapped_filename , default_quality_value , flag_ignore_sequence_information);
	//printf("\nTotal mapped reads %lld", total_mapped_reads);
	//printf("\n");
	return 0;
}
