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
	char *temp;

	long long int start;
	long long int end;
	long long int abridge_match_start_index; // Start index of abridge index - need for multiple cluster overlaps
	long long int abridge_match_end_index; //

	struct Abridge_Index *abridge_index;
	struct Abridge_Index *genome_index;
	struct Whole_Genome_Sequence *single_genome_sequence;

	char **split_on_newline;
	char **split_on_tab;
	char **split_on_dash;
	char **split_on_comma;
	char **split_on_delimiter;

	short int flag_ignore_mismatches;
	short int flag_ignore_soft_clippings;
	short int flag_ignore_unmapped_sequences;
	short int flag_ignore_quality_score;

	int i;
	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	strcpy(abridge_index_filename, argv[1]);
	strcpy(genome_filename, argv[2]);
	strcpy(genome_index_filename, argv[3]);
	strcpy(chromosome, argv[4]);
	start = strtol(argv[5], &temp, 10);
	end = strtol(argv[6], &temp, 10);

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

	single_genome_sequence = (struct Whole_Genome_Sequence*) malloc(sizeof(struct Whole_Genome_Sequence));

	/********************************************************************/

	readAbridgeIndex(abridge_index, abridge_index_filename, split_on_newline, &flag_ignore_mismatches, &flag_ignore_soft_clippings, &flag_ignore_unmapped_sequences, &flag_ignore_quality_score);
	readGenomeIndex(genome_index, genome_index_filename, split_on_newline);
	readInGenomeSequenceSingleChromosome(single_genome_sequence, chromosome, genome_filename, genome_index);
	findReadClusterFromAbridgeIndex(abridge_index, chromosome, start, end, &abridge_match_start_index, &abridge_match_end_index);

	for (i = abridge_match_start_index; i <= abridge_match_end_index; i++)
		printf("\n%s %lld %lld", abridge_index->chromosome[i], abridge_index->start[i], abridge_index->end[i]);
	printf("\n");
	return 0;
}
