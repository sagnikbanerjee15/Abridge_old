# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

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

void swapPointers(char **a, char **b)
{
	char *temp;
	temp = *a;
	*a = *b;
	*b = temp;

}

void reModeliCIGARS(char **split_icigars, char **split_num_reads, int number_of_cigars, char **split_icigars_cp, char **split_num_reads_cp, char **split_icigars_final, char **split_num_reads_final, char *replacement_character, int line_num)
{
	int i;
	int j;
	int final_index;

	// Replace the special matching characters by M to compare different iCIGARs
	for (i = 0; i < number_of_cigars; i++)
		for (j = 0; j < split_icigars_cp[i][j] != '\0'; j++)
			switch (split_icigars_cp[i][j])
			{
				case 'B':
					replacement_character[i] = 'B';
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
				case 'U':
					replacement_character[i] = 'U';
					split_icigars_cp[i][j] = 'M';
					break;
			}

	final_index = 0;
	for (i = 0; i < number_of_cigars; i++)
	{
		while (split_icigars_cp[i][0] == 'X' && i < number_of_cigars)
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
	if (line_num == -1)
	{
		for (i = 0; i < number_of_cigars; i++)
			printf("\n%d) Before: %s After: %s", i + 1, split_icigars[i], split_icigars_final[i]);
		printf("\n");
	}

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
	int line_num = 0;

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
		split_line[i] = (char*) malloc(sizeof(char) * MAX_ICIGAR_LENGTH);

	write_to_file_col1 = (char*) malloc(sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE);
	write_to_file_col3 = (char*) malloc(sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE);

	write_to_file_col1[0] = '\0';
	write_to_file_col3[0] = '\0';

	line_to_be_written_to_file = (char*) malloc(sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE);
	line_to_be_written_to_file[0] = '\0';

	replacement_character = (char*) malloc(sizeof(char) * MIN_POOL_SIZE * 10);
	split_icigars = (char**) malloc(sizeof(char*) * MAX_POOL_SIZE);
	for (i = 0; i < MIN_POOL_SIZE * 10; i++)
		split_icigars[i] = (char*) malloc(sizeof(char) * MAX_ICIGAR_LENGTH);

	split_icigars_cp = (char**) malloc(sizeof(char*) * MIN_POOL_SIZE * 10);
	for (i = 0; i < MIN_POOL_SIZE * 10; i++)
		split_icigars_cp[i] = (char*) malloc(sizeof(char) * MAX_ICIGAR_LENGTH);

	split_icigars_final = (char**) malloc(sizeof(char*) * MIN_POOL_SIZE * 10);
	for (i = 0; i < MIN_POOL_SIZE * 10; i++)
		split_icigars_final[i] = (char*) malloc(sizeof(char) * MAX_ICIGAR_LENGTH);

	split_num_reads = (char**) malloc(sizeof(char*) * MIN_POOL_SIZE * 10);
	for (i = 0; i < MIN_POOL_SIZE * 10; i++)
		split_num_reads[i] = (char*) malloc(sizeof(char) * 25);

	split_num_reads_cp = (char**) malloc(sizeof(char*) * MIN_POOL_SIZE * 10);
	for (i = 0; i < MIN_POOL_SIZE * 10; i++)
		split_num_reads_cp[i] = (char*) malloc(sizeof(char) * 25);

	split_num_reads_final = (char**) malloc(sizeof(char*) * MIN_POOL_SIZE * 10);
	for (i = 0; i < MIN_POOL_SIZE * 10; i++)
		split_num_reads_final[i] = (char*) malloc(sizeof(char) * 25);

	//line = (char*) malloc(sizeof(char) * MAX_ICIGAR_LENGTH);
	/********************************************************************/
	while ((line_len = getline(&line, &len, fhr)) != -1)
	{
		line_num++;
		if (line_num == 1)
		{
			fprintf(fhw, "%s", line);
			continue;
		}
		//printf("\nLine num: %d line_len %d len %d", line_num, line_len, len);
		//fflush(stdout);
		splitByDelimiter(line, '\t', split_line);
		//printf("\nLine: %s Col1: %s Length: %d", line, split_line[0], strlen(split_line[0]));
		//fflush(stdout);
		if (isCommaInLine(split_line[1]) == 1)
		{
			splitByDelimiter(split_line[1], ',', split_icigars);
			number_of_cigars = splitByDelimiter(split_line[2], ',', split_num_reads);
			//continue;
			//printf("\n%s", line);
			//printf("\n icigar: %s number_of_cigars: %d", split_line[1], number_of_cigars);
			splitByDelimiter(split_line[1], ',', split_icigars_cp);
			splitByDelimiter(split_line[2], ',', split_num_reads_cp);
			if (line_num == -1)
			{
				printf("\nLine read from file: %s", line);
				printf("\nMiddle column: %s", split_line[1]);
				printf("\nNumber of cigars: %d", number_of_cigars);
				for (i = 0; i < number_of_cigars; i++)
					printf("\n%d) Before function call: %s", i + 1, split_icigars[i]);
				printf("\n====================================================================================================================");
			}

			reModeliCIGARS(split_icigars, split_num_reads, number_of_cigars, split_icigars_cp, split_num_reads_cp, split_icigars_final, split_num_reads_final, replacement_character, line_num);
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
		free(line);
		line = NULL;
	}
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

	compressSimilarAlignments(input_filename, output_abridgefilename);
	return 0;
}
