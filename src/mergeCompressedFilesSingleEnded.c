# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

void seekFilePointersToAppropriatePosition(char *chromosome, FILE **fhr, int number_of_files_to_be_compressed, char **split_on_tab, struct Chromosome_Starting_Byte **starting_bytes, short int *chromosome_present)
{
	printf("\nAdjusting file pointer for chromosome %s", chromosome);
	fflush(stdout);
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	size_t len;
	ssize_t line_len;

	int i;
	int j;
	int index_of_chromosome;

	char *line = NULL;
	char temp[1000];
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/

	/********************************************************************/
	for (i = 0; i < number_of_files_to_be_compressed; i++)
	{
		index_of_chromosome = -1;
		for (j = 0; starting_bytes[i]->number_of_chromosomes; j++)
		{
			if (strcmp(chromosome, starting_bytes[i]->name[j]) == 0)
			{
				index_of_chromosome = j;
				chromosome_present[i] = 1;
				break;
			}
		}
		if (index_of_chromosome != -1) fseek(fhr[i], starting_bytes[i]->start_byte_in_pass2_file[index_of_chromosome], SEEK_SET);
	}

}
void mergeAbridgeCompressedFiles(char **pass2_filenames, FILE **fhr, int number_of_files_to_be_compressed, char **split_on_tab, struct Chromosome_Info *chromosome_info, struct Chromosome_Starting_Byte **starting_bytes)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	int i;
	int j;
	int test_eof_for_all_files;
	int test_for_end_of_chromsome_alignments;
	int number_of_columns;

	long long int lowest_position;

	unsigned long long int *line_numbers;

	short int *merge_lines_from_these_files;
	short int *reached_end_of_chromsome;
	short int *read_from_file;
	short int time_to_quit;
	short int *chromosome_present;

	size_t *len;
	ssize_t *line_len;

	char *temp; // For strtol
	char **line; // for reading each line

	struct Pass2_Compressed_DS **pass2_compressed_ds_instance;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	line_numbers = (unsigned long long int*) malloc(sizeof(unsigned long long int) * number_of_files_to_be_compressed);
	merge_lines_from_these_files = (short int*) malloc(sizeof(short int) * number_of_files_to_be_compressed);
	reached_end_of_chromsome = (short int*) malloc(sizeof(short int) * number_of_files_to_be_compressed);
	chromosome_present = (short int*) malloc(sizeof(short int) * number_of_files_to_be_compressed);
	len = (size_t*) malloc(sizeof(size_t) * number_of_files_to_be_compressed);
	line_len = (ssize_t*) malloc(sizeof(ssize_t) * number_of_files_to_be_compressed);
	line = (char**) malloc(sizeof(char*) * number_of_files_to_be_compressed);
	pass2_compressed_ds_instance = (struct Pass2_Compressed_DS**) malloc(sizeof(struct Pass2_Compressed_DS*) * number_of_files_to_be_compressed);

	for (i = 0; i < number_of_files_to_be_compressed; i++)
		pass2_compressed_ds_instance[i] = allocateMemoryPass2_Compressed_DS();
	/********************************************************************/
	for (j = 0; j < chromosome_info->number_of_chromosomes; j++)
	{
		//if (strcmp(chromosome_info->name[j], "KI270591.1") != 0) continue;
		printf("\nProcessing chromsosome %s", chromosome_info->name[j]);
		fflush(stdout);
		for (i = 0; i < number_of_files_to_be_compressed; i++)
		{
			fhr[i] = fopen(pass2_filenames[i], "rb");
			if (fhr[i] == NULL)
			{
				printf("Error! File %s not found", pass2_filenames[i]);
				exit(1);
			}
			line[i] = NULL;
			len[i] = 0;
			read_from_file[i] = 1;
			line_numbers[i] = 0;
			line_len[i] = 0;
			chromosome_present[i] = 0;
		}
		printf("\nInitialization complete");
		fflush(stdout);
		seekFilePointersToAppropriatePosition(chromosome_info->name[j], fhr, number_of_files_to_be_compressed, split_on_tab, starting_bytes, chromosome_present);
		printf("\nFile Pointers have been set");
		fflush(stdout);
		for (i = 0; i < number_of_files_to_be_compressed; i++)
		{
			pass2_compressed_ds_instance[i]->position = 0;
			//line_numbers[i] = 0;
			reached_end_of_chromsome[i] = 0;
		}
		while (1)
		{
			time_to_quit = 1;
			for (i = 0; i < number_of_files_to_be_compressed; i++)
			{
				if (read_from_file[i] == 1 && reached_end_of_chromsome[i] == 0 && chromosome_present[i] == 1)
				{
					line_len[i] = getline(&line[i], &len[i], fhr[i]);
					if (line_len[i] == -1)
					{
						reached_end_of_chromsome[i] = 1;
						continue;
					}
					if (line[i][0] == '@') // New chromosome encountered
					{
						//printf("\nLine read from file %s", line[i]);
						line_len[i] = -1;
						reached_end_of_chromsome[i] = 1;
						continue;
					}
					//line_numbers[i]++;
					number_of_columns = splitByDelimiter(line[i], '\t', split_on_tab);
					if (number_of_columns == 1)
					{
						pass2_compressed_ds_instance[i]->position++;
						strcpy(pass2_compressed_ds_instance[i]->col2, split_on_tab[0]);
						time_to_quit = 0;
					}
					else if (number_of_columns == 2)
					{
						pass2_compressed_ds_instance[i]->position += strtol(split_on_tab[0], &temp, 10);
						strcpy(pass2_compressed_ds_instance[i]->col2, split_on_tab[1]);
						time_to_quit = 0;
					}
				}
			}

			// Find lowest position
			lowest_position = -1;
			merge_lines_from_these_files[0] = 0;
			for (i = 0; i < number_of_files_to_be_compressed; i++)
			{
				if (lowest_position == -1) lowest_position = pass2_compressed_ds_instance[i]->position;
				if (lowest_position > pass2_compressed_ds_instance[i]->position) lowest_position = pass2_compressed_ds_instance[i]->position;
				merge_lines_from_these_files[i] = 0;
			}

			for (i = 0; i < number_of_files_to_be_compressed; i++)
			{
				if (pass2_compressed_ds_instance[i]->position == lowest_position)
				{
					read_from_file[i] = 1;
					merge_lines_from_these_files[i] = 1;
				}
				else read_from_file[i] = 0;
			}

			test_for_end_of_chromsome_alignments = 0;
			for (i = 0; i < number_of_files_to_be_compressed; i++)
				test_for_end_of_chromsome_alignments += reached_end_of_chromsome[i];
			if (test_for_end_of_chromsome_alignments == number_of_files_to_be_compressed) break;
			if (time_to_quit) break;
			printf("\n%s ", chromosome_info->name[j]);
			for (i = 0; i < number_of_files_to_be_compressed; i++)
				printf("%lld ", pass2_compressed_ds_instance[i]->position);
			for (i = 0; i < number_of_files_to_be_compressed; i++)
				printf("%d", merge_lines_from_these_files[i]);
			printf(" %lf", (double) lowest_position / (double) chromosome_info->length[j]);
			fflush(stdout);
		}
	}
}

int isChromosomePresent(struct Chromosome_Info *chromosome_info, char *chromosome)
{
	int i;

	for (i = 0; i < chromosome_info->number_of_chromosomes; i++)
		if (strcmp(chromosome_info->name[i], chromosome) == 0) return 1;
	return 0;
}

void collectChromosomeInformationFromAllFiles(char **pass2_filenames, int number_of_files_to_be_compressed, char **split_on_tab, struct Chromosome_Info *chromosome_info, struct Chromosome_Starting_Byte **starting_bytes)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	int i;
	int j;

	size_t len;
	ssize_t line_len;

	char temp[1000];
	char *temp_strtol;
	char *line = NULL; // for reading each line
	char *index_filename;

	FILE *fhr;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	index_filename = (char*) malloc(sizeof(char) * 10000);
	/********************************************************************/
	for (i = 0; i < number_of_files_to_be_compressed; i++)
	{
		fhr = fopen(pass2_filenames[i], "rb");
		if (fhr == NULL)
		{
			printf("Error! File %s not found", pass2_filenames[i]);
			exit(1);
		}
		while ((line_len = getline(&line, &len, fhr)) != -1)
		{
			if (line[0] == '@')
			{
				splitByDelimiter(line, '\t', split_on_tab);
				strcpy(temp, split_on_tab[1]);
				splitByDelimiter(temp, ':', split_on_tab);

				strcpy(starting_bytes[i]->name[starting_bytes[i]->number_of_chromosomes], split_on_tab[1]);
				starting_bytes[i]->start_byte_in_pass2_file[starting_bytes[i]->number_of_chromosomes] = ftell(fhr);
				starting_bytes[i]->number_of_chromosomes++;
				if (isChromosomePresent(chromosome_info, split_on_tab[1])) continue;

				strcpy(chromosome_info->name[chromosome_info->number_of_chromosomes], split_on_tab[1]);

				splitByDelimiter(line, '\t', split_on_tab);
				strcpy(temp, split_on_tab[2]);
				splitByDelimiter(temp, ':', split_on_tab);
				chromosome_info->length[chromosome_info->number_of_chromosomes] = strtol(split_on_tab[1], &temp_strtol, 10);
				chromosome_info->number_of_chromosomes++;

			}
		}
	}

	for (i = 0; i < number_of_files_to_be_compressed; i++)
		for (j = 0; j < starting_bytes[i]->number_of_chromosomes; j++)
			printf("\n%d %s %lld", i, starting_bytes[i]->name[j], starting_bytes[i]->start_byte_in_pass2_file[j]);
	//exit(1);
	fflush(stdout);

}

int main(int argc, char *argv[])
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char **pass2_filenames;
	char **split_on_tab;

	char abridge_index_filename[FILENAME_LENGTH];

	int number_of_files_to_be_compressed;
	int i;

	short int **flags;
	short int chromosome_order_check;

	struct Abridge_Index **abridge_index;
	struct Abrige_Index *super_index;
	struct Chromosome_Info *chromosome_info;
	struct Chromosome_Starting_Byte **starting_bytes;

	FILE **fhr;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	chromosome_info = allocateMemoryChromosome_Info();
	number_of_files_to_be_compressed = argc - 1;
	pass2_filenames = (char**) malloc(sizeof(char*) * number_of_files_to_be_compressed);
	flags = (short int**) malloc(sizeof(short int*) * number_of_files_to_be_compressed);
	starting_bytes = (struct Chromosome_Starting_Byte**) malloc(sizeof(struct Chromosome_Starting_Byte*) * number_of_files_to_be_compressed);
	abridge_index = (struct Abridge_Index**) malloc(sizeof(struct Abridge_Index*) * number_of_files_to_be_compressed);
	fhr = (FILE**) malloc(sizeof(FILE*) * number_of_files_to_be_compressed);
	for (i = 0; i < number_of_files_to_be_compressed; i++)
	{
		pass2_filenames[i] = (char*) malloc(sizeof(char) * FILENAME_LENGTH);
		strcpy(pass2_filenames[i], argv[i + 1]);
		flags[i] = (short int*) malloc(sizeof(short int) * 5);
		abridge_index[i] = allocateMemoryAbridge_Index();
		starting_bytes[i] = allocateMemoryChromosome_Starting_Byte();
	}

	split_on_tab = (char**) malloc(sizeof(char*) * ROWS * 10);
	for (i = 0; i < ROWS * 10; i++)
		split_on_tab[i] = (char*) malloc(sizeof(char) * COLS * 10);

	/********************************************************************/

	/*
	 chromosome_order_check = verifyChromosomeOrder(pass2_filenames, number_of_files_to_be_compressed, split_on_tab);
	 if (chromosome_order_check == -1)
	 {
	 printf("\nNumber of chromosomes is not same for all files");
	 exit(1);
	 }
	 else if (chromosome_order_check == -2)
	 {
	 printf("\nChromosomes mismatch");
	 exit(1);
	 }*/

	collectChromosomeInformationFromAllFiles(pass2_filenames, number_of_files_to_be_compressed, split_on_tab, chromosome_info, starting_bytes);
	mergeAbridgeCompressedFiles(pass2_filenames, fhr, number_of_files_to_be_compressed, split_on_tab, chromosome_info, starting_bytes);
	printf("\n");
	return 0;
}
