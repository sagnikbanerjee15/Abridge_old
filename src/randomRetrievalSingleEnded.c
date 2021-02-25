# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include <math.h>
# include <stdbool.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

int main(int argc, char *argv[])
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char abridge_index_filename[FILENAME_LENGTH];
	char genome_filename[FILENAME_LENGTH];
	char genome_index_filename[FILENAME_LENGTH];
	char chromosome[1000];
	char pass2_filename[FILENAME_LENGTH];
	char output_sam_filename[FILENAME_LENGTH];
	char default_quality_value[10];
	char read_prefix[1000];

	char *temp;
	char *buffer;

	char **split_on_newline;
	char **split_on_tab;
	char **split_on_dash;
	char **split_on_comma;
	char **split_on_delimiter;

	struct Abridge_Index *abridge_index;
	struct Abridge_Index *genome_index;
	struct Whole_Genome_Sequence *single_genome_sequence;
	struct Sam_Alignment *sam_alignment_instance;

	short int flag_ignore_mismatches;
	short int flag_ignore_soft_clippings;
	short int flag_ignore_unmapped_sequences;
	short int flag_ignore_quality_score;
	short int flag_ignore_sequence_information;

	int i;
	int fseek_ret_val;
	int fread_ret_val;
	int number_of_entries_in_cluster;
	int sam_alignment_pool_index;

	unsigned long long int read_number = 1;
	unsigned long long int from = -1;
	unsigned long long int to = -1;

	long long int start;
	long long int end;
	long long int abridge_match_start_index; // Start index of abridge index - need for multiple cluster overlaps
	long long int abridge_match_end_index; //
	long long int total_mapped_reads = 0;

	FILE *fhr;
	FILE *fhw;
	/********************************************************************
	 * Variable initialization
	 ********************************************************************/

	strcpy(abridge_index_filename, argv[1]);
	strcpy(genome_filename, argv[2]);
	strcpy(genome_index_filename, argv[3]);
	strcpy(chromosome, argv[4]);
	strcpy(pass2_filename, argv[5]);
	start = strtol(argv[6], &temp, 10);
	end = strtol(argv[7], &temp, 10);
	strcpy(output_sam_filename, argv[8]);
	strcpy(default_quality_value, argv[9]);
	flag_ignore_sequence_information = strtol(argv[10], &temp, 10);
	strcpy(read_prefix, argv[11]);

	fhr = fopen(pass2_filename, "rb");
	if (fhr == NULL)
	{
		printf("Error! File not found");
		exit(1);
	}
	fhw = fopen(output_sam_filename, "w");
	if (fhw == NULL)
	{
		printf("Error! File cannot be opened for writing");
		exit(1);
	}
	abridge_index = allocateMemoryAbridge_Index();
	genome_index = allocateMemoryAbridge_Index();

	split_on_newline = (char**) malloc(sizeof(char*) * ROWS * 10);
	for (i = 0; i < ROWS * 10; i++)
		split_on_newline[i] = (char*) malloc(sizeof(char) * COLS * 10);

	split_on_tab = (char**) malloc(sizeof(char*) * ROWS);
	for (i = 0; i < ROWS; i++)
		split_on_tab[i] = (char*) malloc(sizeof(char) * COLS);

	split_on_dash = (char**) malloc(sizeof(char*) * ROWS);
	for (i = 0; i < ROWS; i++)
		split_on_dash[i] = (char*) malloc(sizeof(char) * COLS);

	split_on_comma = (char**) malloc(sizeof(char*) * ROWS * 10);
	for (i = 0; i < ROWS * 10; i++)
		split_on_comma[i] = (char*) malloc(sizeof(char) * COLS * 10);

	split_on_delimiter = (char**) malloc(sizeof(char*) * ROWS);
	for (i = 0; i < ROWS; i++)
		split_on_delimiter[i] = (char*) malloc(sizeof(char) * COLS);

	buffer = (char*) malloc(sizeof(char) * MAX_BUFFER_SIZE_FOR_READING_PASS2_FILE);
	single_genome_sequence = (struct Whole_Genome_Sequence*) malloc(sizeof(struct Whole_Genome_Sequence));
	sam_alignment_instance = allocateMemorySam_Alignment();

	/********************************************************************/

	readAbridgeIndex(abridge_index, abridge_index_filename, split_on_newline, &flag_ignore_mismatches, &flag_ignore_soft_clippings, &flag_ignore_unmapped_sequences, &flag_ignore_quality_score);
	readGenomeIndex(genome_index, genome_index_filename, split_on_newline);
	readInGenomeSequenceSingleChromosome(single_genome_sequence, chromosome, genome_filename, genome_index);
	findReadClusterFromAbridgeIndex(abridge_index, chromosome, start, end, &abridge_match_start_index, &abridge_match_end_index);
	writeSequenceHeaders(fhw, genome_filename);

	from = start;
	to = end;
	//printf("\n abridge_match_start_index %d abridge_match_end_index %d", abridge_match_start_index, abridge_match_end_index);
	//fflush(stdout);
	for (i = abridge_match_start_index; i <= abridge_match_end_index; i++)
	{
		//printf("\n%s %lld %lld %lld %lld", abridge_index->chromosome[i], abridge_index->start[i], abridge_index->end[i], abridge_index->start_byte[i], abridge_index->end_byte[i]);
		//fflush(stdout);

		fseek_ret_val = fseek(fhr, abridge_index->start_byte[i], SEEK_SET);
		buffer[0] = '\0';
		fread_ret_val = fread(buffer, 1, abridge_index->end_byte[i] - abridge_index->start_byte[i], fhr);
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
		number_of_entries_in_cluster = splitByDelimiter(buffer, '\n', split_on_newline);
		number_of_entries_in_cluster--; //Last line is always empty

		convertToAlignment(sam_alignment_instance, sam_alignment_pool_index, single_genome_sequence, split_on_newline, sam_alignment_instance, i, abridge_index, number_of_entries_in_cluster, split_on_tab, split_on_dash, split_on_comma, default_quality_value, flag_ignore_mismatches, flag_ignore_soft_clippings, flag_ignore_unmapped_sequences, flag_ignore_quality_score, flag_ignore_sequence_information, &read_number, &total_mapped_reads, read_prefix, from, to, fhw);
	}
	//printf("\n");

	return 0;
}
