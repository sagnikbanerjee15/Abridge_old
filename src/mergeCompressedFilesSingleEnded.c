# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

void seekFilePointersToAppropriatePosition(char *chromosome, FILE **fhr, int number_of_files_to_be_compressed, char **split_on_tab, struct Chromosome_Starting_Byte **starting_bytes, short int *chromosome_present)
{
	//printf("\nAdjusting file pointer for chromosome %s - ", chromosome);
	//fflush(stdout);
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
		for (j = 0; j < starting_bytes[i]->number_of_chromosomes; j++)
		{
			//printf("\nSample=%d Chromosome_num=%d %s %s Number_of_chromosomes=%d", i, j, chromosome, starting_bytes[i]->name[j], starting_bytes[i]->number_of_chromosomes);
			//fflush(stdout);
			if (strcmp(chromosome, starting_bytes[i]->name[j]) == 0)
			{
				index_of_chromosome = j;
				chromosome_present[i] = 1;
				//printf("%d", i);
				//fflush(stdout);
				break;
			}
		}
		if (index_of_chromosome != -1)
		{
			fseek(fhr[i], starting_bytes[i]->start_byte_in_pass2_file[index_of_chromosome], SEEK_SET);
			//getline(&line, &len, fhr[i]);
			//printf("\nLine Read right after fseek to chromosome %s\n%s", chromosome, line);
			//fflush(stdout);
		}
	}
	//printf("\nExiting file pointer adjustment for chromosome %s ", chromosome);
}

void mergeEntriesFromMultipleFiles(struct Pass2_Compressed_DS **pass2_compressed_ds_instance, short int *merge_lines_from_these_files, int number_of_files_to_be_compressed, char **split_on_comma, char **split_on_dash, struct Merged_Compressed_DS *previous_alignment, struct Merged_Compressed_DS *current_alignment, char *line_to_be_written_to_file, FILE *fhw)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	int i;
	int j;
	int k;

	int number_of_icigars;
	int number_of_reads_for_icigar;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	char icigar[MAX_ICIGAR_LENGTH];
	char *temp_strtol;
	char temp_sprintf[100];

	/********************************************************************/
	for (i = 0; i < number_of_files_to_be_compressed; i++)			//ith file
	{
		if (merge_lines_from_these_files[i] == 1)
		{
			current_alignment->position = pass2_compressed_ds_instance[i]->position;
			number_of_icigars = splitByDelimiter(pass2_compressed_ds_instance[i]->col2, ',', split_on_comma);
			for (k = 0; k < number_of_icigars; k++)			//kth icigar in the ith file
			{
				splitByDelimiter(split_on_comma[k], '-', split_on_dash);
				strcpy(icigar, split_on_dash[0]);
				number_of_reads_for_icigar = strtol(split_on_dash[1], &temp_strtol, 10);
				for (j = 0; j < current_alignment->number_of_unique_cigars; j++)			//jth icigar in current alignment
				{
					if (strcmp(icigar, current_alignment->icigars[j]) == 0)
					{
						current_alignment->number_of_reads[j][i] = number_of_reads_for_icigar;
						break;
					}
				}
				if (j == current_alignment->number_of_unique_cigars)			//icigar not found
				{
					strcpy(current_alignment->icigars[current_alignment->number_of_unique_cigars], icigar);
					current_alignment->number_of_reads[current_alignment->number_of_unique_cigars][i] = number_of_reads_for_icigar;
					current_alignment->number_of_unique_cigars++;
				}
			}
		}
	}
	/*
	 * Create the merged col2 representation
	 */
	current_alignment->col2[0] = '\0';
	for (j = 0; j < current_alignment->number_of_unique_cigars; j++)
	{
		strcat(current_alignment->col2, current_alignment->icigars[j]);
		strcat(current_alignment->col2, "-");
		for (i = 0; i < number_of_files_to_be_compressed; i++)			//ith file
		{
			//printf("\ncurrent_alignment->icigars[j][i]=%lld", current_alignment->number_of_reads[j][i]);
			sprintf(temp_sprintf, "%d", current_alignment->number_of_reads[j][i]);
			strcat(current_alignment->col2, temp_sprintf);
			if (i != number_of_files_to_be_compressed - 1)
			strcat(current_alignment->col2, "|");
		}
		if (j != current_alignment->number_of_unique_cigars - 1)
		strcat(current_alignment->col2, ",");
	}
	//printf("\nPosition: %lld COMPRESSED_ICIGAR: %s", current_alignment->position, current_alignment->col2);

	/*
	 * Writing to file
	 */
	line_to_be_written_to_file[0] = '\0';
	if (previous_alignment->position == 0)
	{
		sprintf(temp_sprintf, "%d", current_alignment->position);
		strcat(line_to_be_written_to_file, temp_sprintf);
		strcat(line_to_be_written_to_file, "\t");
		strcat(line_to_be_written_to_file, current_alignment->col2);
		strcat(line_to_be_written_to_file, "\n");
		fprintf(fhw, "%s", line_to_be_written_to_file);
	}
	else
	{
		if (current_alignment->position - previous_alignment->position > 1)
		{
			sprintf(temp_sprintf, "%d", current_alignment->position - previous_alignment->position);
			strcat(line_to_be_written_to_file, temp_sprintf);
			strcat(line_to_be_written_to_file, "\t");
		}
		strcat(line_to_be_written_to_file, current_alignment->col2);
		strcat(line_to_be_written_to_file, "\n");
		fprintf(fhw, "%s", line_to_be_written_to_file);
	}

	/*
	 * Copy the position to previous alignment. Can skip rest of the fields
	 */
	previous_alignment->position = current_alignment->position;

	/*
	 * Reinitialize current alignment
	 */
	current_alignment->col1[0] = '\0';
	current_alignment->col2[0] = '\0';
	current_alignment->number_of_unique_cigars = 0;
	current_alignment->position = 0;
	for (i = 0; i < MAX_UNIQUE_CIGARS; i++)
	{
		for (j = 0; j < MAX_FILES_FOR_MERGING; j++)
			current_alignment->number_of_reads[i][j] = 0;
	}
	for (i = 0; i < MAX_UNIQUE_CIGARS; i++)
		current_alignment->icigars[i][0] = '\0';
}

void mergeAbridgeCompressedFiles(char **pass2_filenames, FILE **fhr, int number_of_files_to_be_compressed, char **split_on_tab, struct Chromosome_Info *chromosome_info, struct Chromosome_Starting_Byte **starting_bytes, char *merged_output_filename, char *chromosome_to_be_processed)
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
	char **split_on_comma;
	char **split_on_dash;
	char *line_to_be_written_to_file;
	char temp_sprintf[100];

	struct Pass2_Compressed_DS **pass2_compressed_ds_instance;
	struct Merged_Compressed_DS *previous_alignment;
	struct Merged_Compressed_DS *current_alignment;

	FILE *fhw;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	line_numbers = (unsigned long long int*) malloc(sizeof(unsigned long long int) * number_of_files_to_be_compressed);
	merge_lines_from_these_files = (short int*) malloc(sizeof(short int) * number_of_files_to_be_compressed);
	reached_end_of_chromsome = (short int*) malloc(sizeof(short int) * number_of_files_to_be_compressed);
	read_from_file = (short int*) malloc(sizeof(short int) * number_of_files_to_be_compressed);
	chromosome_present = (short int*) malloc(sizeof(short int) * number_of_files_to_be_compressed);
	len = (size_t*) malloc(sizeof(size_t) * number_of_files_to_be_compressed);
	line_len = (ssize_t*) malloc(sizeof(ssize_t) * number_of_files_to_be_compressed);
	line = (char**) malloc(sizeof(char*) * number_of_files_to_be_compressed);
	pass2_compressed_ds_instance = (struct Pass2_Compressed_DS**) malloc(sizeof(struct Pass2_Compressed_DS*) * number_of_files_to_be_compressed);
	split_on_comma = (char**) malloc(sizeof(char*) * ROWS * 10);
	split_on_dash = (char**) malloc(sizeof(char*) * ROWS * 10);
	line_to_be_written_to_file = (char*) malloc(sizeof(char) * MAX_GENERAL_LEN);

	for (i = 0; i < number_of_files_to_be_compressed; i++)
		pass2_compressed_ds_instance[i] = allocateMemoryPass2_Compressed_DS();

	split_on_comma = (char**) malloc(sizeof(char*) * ROWS * 10);
	for (i = 0; i < ROWS * 10; i++)
		split_on_comma[i] = (char*) malloc(sizeof(char) * COLS * 10);

	split_on_dash = (char**) malloc(sizeof(char*) * ROWS * 10);
	for (i = 0; i < ROWS * 10; i++)
		split_on_dash[i] = (char*) malloc(sizeof(char) * COLS * 10);

	previous_alignment = allocateMemoryMerged_Compressed_DS();
	current_alignment = allocateMemoryMerged_Compressed_DS();

	fhw = fopen(merged_output_filename, "w");
	if (fhw == NULL)
	{
		printf("\nError opening file %s for writing", merged_output_filename);
		exit(1);
	}

	/********************************************************************/
	for (j = 0; j < chromosome_info->number_of_chromosomes; j++)
	{
		if (strcmp(chromosome_to_be_processed, chromosome_info->name[j]) != 0) continue;
		//printf("\nProcessing chromsosome %s number_of_files_to_be_compressed %d", chromosome_info->name[j], number_of_files_to_be_compressed);
		//fflush(stdout);
		for (i = 0; i < number_of_files_to_be_compressed; i++)
		{
			//printf("\nInitialization ongoing for file no %d", i);
			//fflush(stdout);
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
		/*
		 * Write chromosome info to merged file
		 */
		line_to_be_written_to_file[0] = '\0';
		strcat(line_to_be_written_to_file, "@SQ");
		strcat(line_to_be_written_to_file, "\t");
		strcat(line_to_be_written_to_file, "SN:");
		strcat(line_to_be_written_to_file, chromosome_info->name[j]);
		strcat(line_to_be_written_to_file, "\t");
		strcat(line_to_be_written_to_file, "LN:");
		sprintf(temp_sprintf, "%d", chromosome_info->length[j]);
		strcat(line_to_be_written_to_file, temp_sprintf);
		strcat(line_to_be_written_to_file, "\n");
		fprintf(fhw, "%s", line_to_be_written_to_file);

		//printf("\nInitialization complete");
		//fflush(stdout);
		seekFilePointersToAppropriatePosition(chromosome_info->name[j], fhr, number_of_files_to_be_compressed, split_on_tab, starting_bytes, chromosome_present);
		/*printf(" ");
		 for (i = 0; i < number_of_files_to_be_compressed; i++)
		 {
		 printf("%d", chromosome_present[i]);
		 fflush(stdout);
		 }
		 */
		//continue;
		//printf("\nFile Pointers have been set");
		//fflush(stdout);
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
					//printf("%s", line[i]);
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
			test_for_end_of_chromsome_alignments = 0;
			for (i = 0; i < number_of_files_to_be_compressed; i++)
				test_for_end_of_chromsome_alignments += reached_end_of_chromsome[i];
			if (test_for_end_of_chromsome_alignments == number_of_files_to_be_compressed) break;
			if (time_to_quit) break;

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
			mergeEntriesFromMultipleFiles(pass2_compressed_ds_instance, merge_lines_from_these_files, number_of_files_to_be_compressed, split_on_comma, split_on_dash, previous_alignment, current_alignment, line_to_be_written_to_file, fhw);

			/*
			 printf("\n%s ", chromosome_info->name[j]);
			 for (i = 0; i < number_of_files_to_be_compressed; i++)
			 printf("%lld ", pass2_compressed_ds_instance[i]->position);
			 for (i = 0; i < number_of_files_to_be_compressed; i++)
			 printf("%d", merge_lines_from_these_files[i]);
			 printf(" %lf", (double) lowest_position / (double) chromosome_info->length[j]);
			 fflush(stdout);
			 */
		}
	}
	fclose(fhw);
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
		fclose(fhr);
	}
	/*for (i = 0; i < number_of_files_to_be_compressed; i++)
	 {
	 for (j = 0; j < starting_bytes[i]->number_of_chromosomes; j++)
	 {
	 printf("\n%d %s %lld %d", i, starting_bytes[i]->name[j], starting_bytes[i]->start_byte_in_pass2_file[j], starting_bytes[i]->number_of_chromosomes);
	 fflush(stdout);
	 }
	 }*/
	//exit(1);
}

int main(int argc, char *argv[])
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char **pass2_filenames;
	char **split_on_tab;

	char abridge_index_filename[FILENAME_LENGTH];
	char merged_output_filename[FILENAME_LENGTH];
	char chromosome_to_be_processed[1000];

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
	number_of_files_to_be_compressed = argc - 3;
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
	//printf("\n%d", argc);
	//fflush(stdout);
	strcpy(merged_output_filename, argv[i + 1]);
	//printf("\n%s", merged_output_filename);
	//fflush(stdout);
	strcpy(chromosome_to_be_processed, argv[i + 2]);
	//printf("\n%s", chromosome_to_be_processed);
	//fflush(stdout);
	//exit(1);

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
	mergeAbridgeCompressedFiles(pass2_filenames, fhr, number_of_files_to_be_compressed, split_on_tab, chromosome_info, starting_bytes, merged_output_filename, chromosome_to_be_processed);
	printf("\n");
	return 0;
}
