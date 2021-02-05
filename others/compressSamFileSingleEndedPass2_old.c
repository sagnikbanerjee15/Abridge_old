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

void compressSimilarAlignments(char *input_filename, char *output_abridgefilename)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/

	FILE *fhr;
	FILE *fhw;

	int i;
	int pass1_compressed_DS_pool_index;

	size_t len = 0;
	ssize_t line_len;

	char *temp; //Useless
	char *line = NULL; // for reading each line
	char **split_line; // List of strings to store each element of a single alignment
	char *write_to_file_col1;
	char *write_to_file_col2;
	char *write_to_file_col3;
	char *line_to_be_written_to_file;

	struct Pass1_Compressed_DS **pass1_compressed_DS_pool;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	pass1_compressed_DS_pool = (struct Pass1_Compressed_DS**) malloc(sizeof(struct Pass1_Compressed_DS*) * MAX_POOL_SIZE);
	for (i = 0; i < MAX_POOL_SIZE; i++)
		pass1_compressed_DS_pool[i] = allocateMemoryPass1_Compressed_DS();

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
	/********************************************************************/

	pass1_compressed_DS_pool_index = 0;
	while ((line_len = getline(&line, &len, fhr)) != -1)
	{
		splitByDelimiter(line, '\t', split_line);
		/*printf("\n%d %d", pass1_compressed_DS_pool_index, MAX_POOL_SIZE);
		 fflush(stdout);*/
		if (isCommaInLine(split_line[0]) == 0)
		{
			//printf("\n No comma: %s %s %s", split_line[0], split_line[1], split_line[2]);
			if (pass1_compressed_DS_pool_index == 0)
			{
				strcpy(pass1_compressed_DS_pool[pass1_compressed_DS_pool_index]->col1, split_line[0]);
				strcpy(pass1_compressed_DS_pool[pass1_compressed_DS_pool_index]->col2, split_line[1]);
				strcpy(pass1_compressed_DS_pool[pass1_compressed_DS_pool_index]->col3, split_line[2]);
				pass1_compressed_DS_pool_index++;
			}
			else if (strcmp(pass1_compressed_DS_pool[pass1_compressed_DS_pool_index - 1]->col2, split_line[1]) == 0)
			{
				strcpy(pass1_compressed_DS_pool[pass1_compressed_DS_pool_index]->col1, split_line[0]);
				strcpy(pass1_compressed_DS_pool[pass1_compressed_DS_pool_index]->col2, split_line[1]);
				strcpy(pass1_compressed_DS_pool[pass1_compressed_DS_pool_index]->col3, split_line[2]);
				pass1_compressed_DS_pool_index++;
			}
			else
			{
				writeToFile(pass1_compressed_DS_pool, pass1_compressed_DS_pool_index, fhw, write_to_file_col1, write_to_file_col3, line_to_be_written_to_file);
				pass1_compressed_DS_pool_index = 0;
				strcpy(pass1_compressed_DS_pool[pass1_compressed_DS_pool_index]->col1, split_line[0]);
				strcpy(pass1_compressed_DS_pool[pass1_compressed_DS_pool_index]->col2, split_line[1]);
				strcpy(pass1_compressed_DS_pool[pass1_compressed_DS_pool_index]->col3, split_line[2]);
				pass1_compressed_DS_pool_index++;
			}
		}
		else
		{
			writeToFile(pass1_compressed_DS_pool, pass1_compressed_DS_pool_index, fhw, write_to_file_col1, write_to_file_col3, line_to_be_written_to_file);
			pass1_compressed_DS_pool_index = 0;
			// Write the line read from file without making any changes
			fprintf(fhw, "%s", line);
			//exit(1);
		}
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
