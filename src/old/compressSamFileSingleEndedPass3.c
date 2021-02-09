# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

/*
 * Usage:
 * compressSamFileSingleEnded \ argv[0]
 * <reference_sequence_id> \ argv[1]
 * <ignore-soft-clippings> \ argv[2]
 * <ignore-mismatches> \ argv[3]
 * <ignore-sequence-information> \ argv[4]
 * <ignore-unmapped-reads> \ argv[5]
 * <inputinputsamfilename> \ argv[6]
 * <outputfilename> argv[7]
 */

void writeToFile(struct Pass1_Compressed_DS **pass1_compressed_DS_pool, int pass1_compressed_DS_pool_index, FILE *fhw, char *write_to_file_col1, char *write_to_file_col3, char *line_to_be_written_to_file)
{
	int i;
	line_to_be_written_to_file[0] = '\0';
	write_to_file_col1[0] = '\0';
	write_to_file_col3[0] = '\0';

	if (pass1_compressed_DS_pool_index == 0) return;
	for (i = 0; i < pass1_compressed_DS_pool_index; i++)
	{
		strcat(write_to_file_col1, pass1_compressed_DS_pool[i]->col1);
		strcat(write_to_file_col1, ";");
		strcat(write_to_file_col3, pass1_compressed_DS_pool[i]->col3);
		strcat(write_to_file_col3, ";");
	}
	write_to_file_col1[strlen(write_to_file_col1) - 1] = '\0'; // Removing the last semi colon
	write_to_file_col3[strlen(write_to_file_col3) - 1] = '\0'; // Removing the last semi colon
	//printf("\n Col1: %s Col3: %s", write_to_file_col1, write_to_file_col3);

	strcat(line_to_be_written_to_file, write_to_file_col1);
	strcat(line_to_be_written_to_file, "\t");
	strcat(line_to_be_written_to_file, pass1_compressed_DS_pool[pass1_compressed_DS_pool_index - 1]->col2);
	strcat(line_to_be_written_to_file, "\t");
	strcat(line_to_be_written_to_file, write_to_file_col3);
	strcat(line_to_be_written_to_file, "\n");
	//printf("\n to file: %s", line_to_be_written_to_file);
	fprintf(fhw, "%s", line_to_be_written_to_file);

	line_to_be_written_to_file[0] = '\0';
	write_to_file_col1[0] = '\0';
	write_to_file_col3[0] = '\0';

}

void reModeliCIGARS(char **split_icigars, char **split_num_reads, int number_of_cigars, char **split_icigars_cp, char **split_num_reads_cp, char **split_icigars_final, char **split_num_reads_final, char *replacement_character)
{
	int i;
	int j;
	int final_index;
	/*
	 printf("\n Before: ");
	 for (i = 0; i < number_of_cigars; i++)
	 printf("%s,", split_icigars[i]);
	 */

	// Replace the special matching characters by M to compare different iCIGARs
	for (i = 0; i < number_of_cigars; i++)
		for (j = 0; j < split_icigars_cp[i][j] != '\0'; j++)
			switch (split_icigars_cp[i][j])
			{
				case 'B':
					replacement_character[i] = 'B';
					split_icigars_cp[i][j] = 'M';
					break;
				case 'C':
					replacement_character[i] = 'C';
					split_icigars_cp[i][j] = 'M';
					break;
				case 'E':
					replacement_character[i] = 'E';
					split_icigars_cp[i][j] = 'M';
					break;
				case 'F':
					replacement_character[i] = 'F';
					split_icigars_cp[i][j] = 'M';
					break;
				case 'H':
					replacement_character[i] = 'H';
					split_icigars_cp[i][j] = 'M';
					break;
				case 'J':
					replacement_character[i] = 'J';
					split_icigars_cp[i][j] = 'M';
					break;
				case 'K':
					replacement_character[i] = 'K';
					split_icigars_cp[i][j] = 'M';
					break;
				case 'L':
					replacement_character[i] = 'L';
					split_icigars_cp[i][j] = 'M';
					break;
				case 'O':
					replacement_character[i] = 'O';
					split_icigars_cp[i][j] = 'M';
					break;
				case 'P':
					replacement_character[i] = 'P';
					split_icigars_cp[i][j] = 'M';
					break;
				case 'Q':
					replacement_character[i] = 'Q';
					split_icigars_cp[i][j] = 'M';
					break;
				case 'R':
					replacement_character[i] = 'R';
					split_icigars_cp[i][j] = 'M';
					break;
			}

	final_index = 0;
	for (i = 0; i < number_of_cigars; i++)
	{
		while (strcmp(split_icigars_cp[i], "X") == 0)
			i++;
		strcpy(split_icigars_final[final_index], split_icigars[i]);
		strcpy(split_num_reads_final[final_index], split_num_reads[i]);
		final_index++;
		for (j = i + 1; j < number_of_cigars; j++)
		{
			if (strcmp(split_icigars_cp[i], split_icigars_cp[j]) == 0)
			{
				split_icigars_cp[j][0] = 'X';
				split_icigars_cp[j][1] = '\0';
				split_icigars_final[final_index][0] = replacement_character[j];
				split_icigars_final[final_index][1] = '\0';
				strcpy(split_num_reads_final[final_index], split_num_reads[j]);
				final_index++;
			}
		}
	}
	/*
	 if (number_of_cigars != final_index) printf("\n number_of_cigars: %d final_index: %d ", number_of_cigars, final_index);
	 printf("\n After: ");
	 for (i = 0; i < number_of_cigars; i++)
	 printf("%s,", split_icigars_final[i]);
	 */

}

void compressSimilarAlignments(char *input_filename, char *output_abridgefilename)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/

	FILE *fhr;
	FILE *fhw;

	int i;
	int number_of_cigars;

	size_t len = 0;
	ssize_t line_len;

	char *temp; //Useless
	char *line = NULL; // for reading each line
	char **split_line; // List of strings to store each element of a single alignment
	char *write_to_file_col1;
	char *write_to_file_col2;
	char *write_to_file_col3;
	char *line_to_be_written_to_file;
	char *replacement_character;
	char **split_icigars;
	char **split_num_reads;
	char **split_icigars_cp;
	char **split_num_reads_cp;
	char **split_icigars_final;
	char **split_num_reads_final;

	struct Pass1_Compressed_DS *pass1_compressed_DS;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	pass1_compressed_DS = (struct Pass1_Compressed_DS*) malloc(sizeof(struct Pass1_Compressed_DS));

	fhr = fopen(input_filename, "r");
	if (fhr == NULL)
	{
		printf("Error! File not found");
		exit(1);
	}
	fhw = fopen(output_abridgefilename, "w");
	if (fhw == NULL)
	{
		printf("File cannot be created");
		exit(1);
	}

	split_line = (char**) malloc(sizeof(char*) * ROWS);
	for (i = 0; i < ROWS; i++)
		split_line[i] = (char*) malloc(sizeof(char) * COLS);

	write_to_file_col1 = (char*) malloc(sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE);
	write_to_file_col3 = (char*) malloc(sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE);

	write_to_file_col1[0] = '\0';
	write_to_file_col3[0] = '\0';

	line_to_be_written_to_file = (char*) malloc(sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE);
	line_to_be_written_to_file[0] = '\0';

	replacement_character = (char*) malloc(sizeof(char) * MAX_POOL_SIZE);
	split_icigars = (char**) malloc(sizeof(char*) * MAX_POOL_SIZE);
	for (i = 0; i < MAX_POOL_SIZE; i++)
		split_icigars[i] = (char*) malloc(sizeof(char) * MAX_ICIGAR_LENGTH);

	split_icigars_cp = (char**) malloc(sizeof(char*) * MAX_POOL_SIZE);
	for (i = 0; i < MAX_POOL_SIZE; i++)
		split_icigars_cp[i] = (char*) malloc(sizeof(char) * MAX_ICIGAR_LENGTH);

	split_icigars_final = (char**) malloc(sizeof(char*) * MAX_POOL_SIZE);
	for (i = 0; i < MAX_POOL_SIZE; i++)
		split_icigars_final[i] = (char*) malloc(sizeof(char) * MAX_ICIGAR_LENGTH);

	split_num_reads = (char**) malloc(sizeof(char*) * MAX_POOL_SIZE);
	for (i = 0; i < MAX_POOL_SIZE; i++)
		split_num_reads[i] = (char*) malloc(sizeof(char) * 25);

	split_num_reads_cp = (char**) malloc(sizeof(char*) * MAX_POOL_SIZE);
	for (i = 0; i < MAX_POOL_SIZE; i++)
		split_num_reads_cp[i] = (char*) malloc(sizeof(char) * 25);

	split_num_reads_final = (char**) malloc(sizeof(char*) * MAX_POOL_SIZE);
	for (i = 0; i < MAX_POOL_SIZE; i++)
		split_num_reads_final[i] = (char*) malloc(sizeof(char) * 25);
	/********************************************************************/

	while ((line_len = getline(&line, &len, fhr)) != -1)
	{
		splitByDelimiter(line, '\t', split_line);
		/*printf("\nLine: %s Col1: %s Length: %d", line, split_line[0], strlen(split_line[0]));
		 fflush(stdout);*/
		/*printf("\n%d %d", pass1_compressed_DS_pool_index, MAX_POOL_SIZE);
		 fflush(stdout);*/
		if (isCommaInLine(split_line[1]) == 1)
		{
			splitByDelimiter(split_line[1], ',', split_icigars);
			number_of_cigars = splitByDelimiter(split_line[2], ',', split_num_reads);
			//printf("\n%s", line);
			//printf("\n icigar: %s number_of_cigars: %d", split_line[1], number_of_cigars);
			splitByDelimiter(split_line[1], ',', split_icigars_cp);
			splitByDelimiter(split_line[2], ',', split_num_reads_cp);
			reModeliCIGARS(split_icigars, split_num_reads, number_of_cigars, split_icigars_cp, split_num_reads_cp, split_icigars_final, split_num_reads_final, replacement_character);
			line_to_be_written_to_file[0] = '\0';
			if (strlen(split_line[0]) != 0)
			{
				strcpy(line_to_be_written_to_file, split_line[0]);
				strcat(line_to_be_written_to_file, "\t");
			}
			for (i = 0; i < number_of_cigars; i++)
			{
				strcat(line_to_be_written_to_file, split_icigars_final[i]);
				strcat(line_to_be_written_to_file, "-");
				strcat(line_to_be_written_to_file, split_num_reads_final[i]);
				if (i != number_of_cigars - 1) strcat(line_to_be_written_to_file, ",");
			}
			/*strcat(line_to_be_written_to_file, "\t");
			 for (i = 0; i < number_of_cigars; i++)
			 {
			 strcat(line_to_be_written_to_file, split_num_reads_final[i]);
			 if (i != number_of_cigars - 1) strcat(line_to_be_written_to_file, ",");
			 }*/
			strcat(line_to_be_written_to_file, "\n");
			fprintf(fhw, "%s", line_to_be_written_to_file);
		}
		else if (line[0] == '@') fprintf(fhw, "%s", line);
		else
		{
			line_to_be_written_to_file[0] = '\0';
			if (strlen(split_line[0]) != 0)
			{
				strcpy(line_to_be_written_to_file, split_line[0]);
				strcat(line_to_be_written_to_file, "\t");
			}
			strcat(line_to_be_written_to_file, split_line[1]);
			strcat(line_to_be_written_to_file, "-");
			strcat(line_to_be_written_to_file, split_line[2]);
			strcat(line_to_be_written_to_file, "\n");
			fprintf(fhw, "%s", line_to_be_written_to_file);
		}
	}
}

void swapCigarFrequenceItems(struct Cigar_Frequency **cigar_freq_pool, int left_pointer, int right_pointer)
{
	struct Cigar_Frequency *temp;
	temp = cigar_freq_pool[left_pointer];
	cigar_freq_pool[left_pointer] = cigar_freq_pool[right_pointer];
	cigar_freq_pool[right_pointer] = temp;
}

int partition(struct Cigar_Frequency **cigar_freq_pool, int left, int right, long int pivot_value)
{
	int left_pointer = left - 1;
	int right_pointer = right;
	while (1)
	{
		while (cigar_freq_pool[++left_pointer]->freq < pivot_value)
		{
			//do nothing
		}
		while (right_pointer > 0 && cigar_freq_pool[--right_pointer]->freq > pivot_value)
		{
			//do nothing
		}
		if (left_pointer >= right_pointer) break;
		else swapCigarFrequenceItems(cigar_freq_pool, left_pointer, right_pointer);
	}
	swapCigarFrequenceItems(cigar_freq_pool, left_pointer, right);
	return left_pointer;
}

void quickSort(struct Cigar_Frequency **cigar_freq_pool, int left, int right)
{
	if (right - left <= 0) return;
	else
	{
		long int pivot_value = cigar_freq_pool[right]->freq;
		int pivot_index = partition(cigar_freq_pool, left, right, pivot_value);
		quickSort(cigar_freq_pool, left, pivot_index - 1);
		quickSort(cigar_freq_pool, pivot_index + 1, right);
	}
}

void reverseCigarFrequencyPool(struct Cigar_Frequency **cigar_freq_pool, int N)
{
	int i;
	for (i = 0; i < N / 2; i++)
		swapCigarFrequenceItems(cigar_freq_pool, i, N - i - 1);
}

void quickSortCigarFrequencyPool(struct Cigar_Frequency **cigar_freq_pool, int N)
{
	quickSort(cigar_freq_pool, 0, N - 1);
	reverseCigarFrequencyPool(cigar_freq_pool, N);
}

int searchicigarInMappingDictionary(struct Pass3_Compression_Symbol_icigar_Mapping **symbol_icigar_mapping, int cigar_freq_pool_index, char *icigar)
{
	int index;
	int i;
	for (i = 0; i < cigar_freq_pool_index; i++)
		if (strcmp(icigar, symbol_icigar_mapping[i]->icigar) == 0) return i;
	return -1;
}

void replaceicigarWithSymbolAndWriteToFile(FILE *fhr, FILE *fhw, struct Pass3_Compression_Symbol_icigar_Mapping **symbol_icigar_mapping, int cigar_freq_pool_index)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	int i;
	int j;
	int number_of_cols;
	int number_of_cigars;

	char *temp; //Useless
	char *line = NULL; // for reading each line
	char **split_line; // List of strings to store each element of a single alignment
	char *line_to_be_written_to_file;
	char *icigar_field;
	char **split_icigars;

	size_t len = 0;
	ssize_t line_len;

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/

	icigar_field = (char*) malloc(sizeof(char) * MAX_ICIGAR_LENGTH);

	split_line = (char**) malloc(sizeof(char*) * ROWS);
	for (i = 0; i < ROWS; i++)
		split_line[i] = (char*) malloc(sizeof(char) * COLS);

	split_icigars = (char**) malloc(sizeof(char*) * MAX_POOL_SIZE);
	for (i = 0; i < MAX_POOL_SIZE; i++)
		split_icigars[i] = (char*) malloc(sizeof(char) * MAX_ICIGAR_LENGTH);

	line_to_be_written_to_file = (char*) malloc(sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE);
	line_to_be_written_to_file[0] = '\0';
	/********************************************************************/

	/*
	 * Write the entire mapping dictionary to file
	 */
	for (i = 0; i < cigar_freq_pool_index; i++)
	{
		strcat(line_to_be_written_to_file, symbol_icigar_mapping[i]->icigar);
		strcat(line_to_be_written_to_file, "~");
		strcat(line_to_be_written_to_file, symbol_icigar_mapping[i]->symbolic_icigar);
		strcat(line_to_be_written_to_file, " ");
	}
	strcat(line_to_be_written_to_file, "\n");
	fprintf(fhw, "%s", line_to_be_written_to_file);
	/********************************************************************/

	while ((line_len = getline(&line, &len, fhr)) != -1)
	{
		if (line[0] == '@') fprintf(fhw, "%s", line);
		else
		{
			number_of_cols = splitByDelimiter(line, '\t', split_line);
			if (number_of_cols == 1) strcpy(icigar_field, split_line[0]);
			else if (number_of_cols == 2) strcpy(icigar_field, split_line[1]);
			else
			{
				//Should never come here
			}
			line_to_be_written_to_file[0] = '\0';
			if (number_of_cols == 2)
			{
				strcat(line_to_be_written_to_file, split_line[0]);
				strcat(line_to_be_written_to_file, "\t");
			}

			number_of_cigars = splitByDelimiter(icigar_field, ',', split_icigars);
			for (i = 0; i < number_of_cigars; i++)
			{
				// icigar has a splice or deletions or soft clips
				if (strchr(split_icigars[i], 'N') != NULL || strchr(split_icigars[i], 'D') != NULL || strchr(split_icigars[i], 'A') != NULL || strchr(split_icigars[i], 'T') != NULL || strchr(split_icigars[i], 'G') != NULL || strchr(split_icigars[i], 'C') != NULL)
				{
					strcat(line_to_be_written_to_file, split_icigars[i]);
					if (i != number_of_cigars - 1) strcat(line_to_be_written_to_file, ",");
					continue;
				}

				for (j = 0; j < 5; j++)
					if (strchr(split_icigars[i], insert_characters[j]) != NULL) break;

				// icigar has insertions
				if (j < 5)
				{
					strcat(line_to_be_written_to_file, split_icigars[i]);
					if (i != number_of_cigars - 1) strcat(line_to_be_written_to_file, ",");
					continue;
				}
				for (j = 0; j < 5; j++)
					if (strchr(split_icigars[i], mismatch_characters[j]) != NULL) break;

				// icigar has mismatching characters
				if (j < 5)
				{
					strcat(line_to_be_written_to_file, split_icigars[i]);
					if (i != number_of_cigars - 1) strcat(line_to_be_written_to_file, ",");
					continue;
				}

				int index = searchicigarInMappingDictionary(symbol_icigar_mapping, cigar_freq_pool_index, split_icigars[i]);
				if (index == -1)
				{
					// No Match found. Write to file the same line as read from file
					strcat(line_to_be_written_to_file, split_icigars[i]);
					continue;
				}
				else
				{
					strcat(line_to_be_written_to_file, symbol_icigar_mapping[index]->symbolic_icigar);
				}
				if (i != number_of_cigars - 1) strcat(line_to_be_written_to_file, ",");
			}
			strcat(line_to_be_written_to_file, "\n");
			fprintf(fhw, "%s", line_to_be_written_to_file);
		}

	}
}

void reCodeiCIGARS(char *input_filename, char *output_abridgefilename)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/

	FILE *fhr;
	FILE *fhw;

	int i;
	int j;
	int number_of_cigars;
	int number_of_cols;
	int cigar_freq_pool_index = 0;

	size_t len = 0;
	ssize_t line_len;

	char *temp; //Useless
	char *line = NULL; // for reading each line
	char **split_line; // List of strings to store each element of a single alignment
	char *write_to_file_col1;
	char *write_to_file_col2;
	char *write_to_file_col3;
	char *line_to_be_written_to_file;
	char *replacement_character;
	char **split_icigars;
	char **split_num_reads;
	char **split_icigars_cp;
	char **split_num_reads_cp;
	char **split_icigars_final;
	char **split_num_reads_final;
	char *icigar_field;

	struct Pass2_Compressed_DS *pass2_compressed_DS;
	struct Pass3_Compression_Symbol_icigar_Mapping **symbol_icigar_mapping;
	struct Cigar_Frequency **cigar_freq_pool;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	pass2_compressed_DS = (struct Pass2_Compressed_DS*) malloc(sizeof(struct Pass2_Compressed_DS));
	icigar_field = (char*) malloc(sizeof(char) * MAX_ICIGAR_LENGTH);

	fhr = fopen(input_filename, "r");
	if (fhr == NULL)
	{
		printf("Error! File not found");
		exit(1);
	}
	fhw = fopen(output_abridgefilename, "w");
	if (fhw == NULL)
	{
		printf("File cannot be created");
		exit(1);
	}

	split_line = (char**) malloc(sizeof(char*) * ROWS);
	for (i = 0; i < ROWS; i++)
		split_line[i] = (char*) malloc(sizeof(char) * COLS);

	write_to_file_col1 = (char*) malloc(sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE);
	write_to_file_col3 = (char*) malloc(sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE);

	write_to_file_col1[0] = '\0';
	write_to_file_col3[0] = '\0';

	line_to_be_written_to_file = (char*) malloc(sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE);
	line_to_be_written_to_file[0] = '\0';

	replacement_character = (char*) malloc(sizeof(char) * MAX_POOL_SIZE);
	split_icigars = (char**) malloc(sizeof(char*) * MAX_POOL_SIZE);
	for (i = 0; i < MAX_POOL_SIZE; i++)
		split_icigars[i] = (char*) malloc(sizeof(char) * MAX_ICIGAR_LENGTH);

	split_icigars_cp = (char**) malloc(sizeof(char*) * MAX_POOL_SIZE);
	for (i = 0; i < MAX_POOL_SIZE; i++)
		split_icigars_cp[i] = (char*) malloc(sizeof(char) * MAX_ICIGAR_LENGTH);

	split_icigars_final = (char**) malloc(sizeof(char*) * MAX_POOL_SIZE);
	for (i = 0; i < MAX_POOL_SIZE; i++)
		split_icigars_final[i] = (char*) malloc(sizeof(char) * MAX_ICIGAR_LENGTH);

	split_num_reads = (char**) malloc(sizeof(char*) * MAX_POOL_SIZE);
	for (i = 0; i < MAX_POOL_SIZE; i++)
		split_num_reads[i] = (char*) malloc(sizeof(char) * 25);

	split_num_reads_cp = (char**) malloc(sizeof(char*) * MAX_POOL_SIZE);
	for (i = 0; i < MAX_POOL_SIZE; i++)
		split_num_reads_cp[i] = (char*) malloc(sizeof(char) * 25);

	split_num_reads_final = (char**) malloc(sizeof(char*) * MAX_POOL_SIZE);
	for (i = 0; i < MAX_POOL_SIZE; i++)
		split_num_reads_final[i] = (char*) malloc(sizeof(char) * 25);

	cigar_freq_pool = (struct Cigar_Frequency**) malloc(sizeof(struct Cigar_Frequency*) * MAX_POOL_SIZE);
	for (i = 0; i < MAX_POOL_SIZE; i++)
		cigar_freq_pool[i] = allocateMemoryCigar_Frequency();

	symbol_icigar_mapping = (struct Pass3_Compression_Symbol_icigar_Mapping**) malloc(sizeof(struct Pass3_Compression_Symbol_icigar_Mapping*) * (MAX_SYMBOLS_FOR_PASS3_COMPRESSION * MAX_SYMBOLS_FOR_PASS3_COMPRESSION));
	for (i = 0; i < (MAX_SYMBOLS_FOR_PASS3_COMPRESSION * MAX_SYMBOLS_FOR_PASS3_COMPRESSION); i++)
		symbol_icigar_mapping[i] = allocateMemoryPass3_Compression_Symbol_icigar_Mapping();

	initializePass3_Compression_Symbol_icigar_MappingPool(symbol_icigar_mapping);
	/********************************************************************/

	while ((line_len = getline(&line, &len, fhr)) != -1)
	{
		if (line[0] == '@') continue;
		else
		{
			number_of_cols = splitByDelimiter(line, '\t', split_line);
			if (number_of_cols == 1) strcpy(icigar_field, split_line[0]);
			else if (number_of_cols == 2) strcpy(icigar_field, split_line[1]);
			else
			{
				//Should never come here
			}

			number_of_cigars = splitByDelimiter(icigar_field, ',', split_icigars);
			for (i = 0; i < number_of_cigars; i++)
			{
				// icigar has a splice or deletions or soft clips
				if (strchr(split_icigars[i], 'N') != NULL || strchr(split_icigars[i], 'D') != NULL || strchr(split_icigars[i], 'A') != NULL || strchr(split_icigars[i], 'T') != NULL || strchr(split_icigars[i], 'G') != NULL || strchr(split_icigars[i], 'C') != NULL)
					continue;

				for (j = 0; j < 5; j++)
					if (strchr(split_icigars[i], insert_characters[j]) != NULL) break;

				// icigar has insertions
				if (j < 5) continue;
				for (j = 0; j < 5; j++)
					if (strchr(split_icigars[i], mismatch_characters[j]) != NULL) break;

				// icigar has mismatching characters
				if (j < 5) continue;

				// Check if cigar already present
				for (j = 0; j < cigar_freq_pool_index; j++)
					if (strcmp(cigar_freq_pool[j]->cigar, split_icigars[i]) == 0)
					{
						cigar_freq_pool[j]->freq++;
						break;
					}
				if (j == cigar_freq_pool_index)
				{
					strcpy(cigar_freq_pool[cigar_freq_pool_index]->cigar, split_icigars[i]);
					cigar_freq_pool[cigar_freq_pool_index]->freq = 1;
					cigar_freq_pool_index++;
				}
				//printf("\n cigar_freq_pool_index: %d %d %f j=%d %s freq=%d", cigar_freq_pool_index, MAX_POOL_SIZE, (float) cigar_freq_pool_index / (float) MAX_POOL_SIZE, j, split_icigars[i], cigar_freq_pool[j]->freq);
			}
		}
	}
	quickSortCigarFrequencyPool(cigar_freq_pool, cigar_freq_pool_index);
	/*
	 for (j = 0; j < cigar_freq_pool_index; j++)
	 printf("\n%d iCIGAR: %s ", cigar_freq_pool[j]->freq, cigar_freq_pool[j]->cigar);
	 */
	if (cigar_freq_pool_index > (MAX_SYMBOLS_FOR_PASS3_COMPRESSION * MAX_SYMBOLS_FOR_PASS3_COMPRESSION)) cigar_freq_pool_index = MAX_SYMBOLS_FOR_PASS3_COMPRESSION * MAX_SYMBOLS_FOR_PASS3_COMPRESSION;
	assignicigarsToSymbols(cigar_freq_pool, cigar_freq_pool_index, symbol_icigar_mapping);
	rewind(fhr);
	replaceicigarWithSymbolAndWriteToFile(fhr, fhw, symbol_icigar_mapping, cigar_freq_pool_index);
}

int main(int argc, char *argv[])
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char input_filename[FILENAME_LENGTH];
	char output_abridgefilename[FILENAME_LENGTH];
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	strcpy(input_filename, argv[1]);
	strcpy(output_abridgefilename, argv[2]);
	/********************************************************************/

	reCodeiCIGARS(input_filename, output_abridgefilename);

	return 0;
}
