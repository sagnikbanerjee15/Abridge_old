# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

void compress(char *input_filename, char *output_abridgefilename)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	FILE *fhr;
	FILE *fhw;

	size_t len = 0;
	ssize_t line_len;

	char *temp; //Useless
	char *line = NULL; // for reading each line
	char **split_line; // List of strings to store each element of a single alignment
	char **split_cigars;

	int i;
	int j;
	int number_of_cigars;
	int cigar_frequency_index = 0;
	int number_of_records_processed = 0;

	struct Cigar_Frequency **cigar_frequency;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/

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

	split_line = (char**) malloc(sizeof(char*) * MAX_CIGAR_LENGTH);
	for (i = 0; i < MAX_CIGAR_LENGTH; i++)
		split_line[i] = (char*) malloc(sizeof(char) * MAX_CIGAR_LENGTH);

	split_cigars = (char**) malloc(sizeof(char*) * MAX_CIGAR_LENGTH);
	for (i = 0; i < MAX_CIGAR_LENGTH; i++)
		split_cigars[i] = (char*) malloc(sizeof(char) * MAX_CIGAR_LENGTH);

	cigar_frequency = (struct Cigar_Frequency**) malloc(sizeof(struct Cigar_Frequency*) * MAX_CIGAR_FREQ_SIZE);
	for (i = 0; i < MAX_CIGAR_FREQ_SIZE; i++)
		cigar_frequency[i] = allocateMemoryCigar_Frequency();

	/********************************************************************/
	while ((line_len = getline(&line, &len, fhr)) != -1)
	{
		number_of_records_processed += 1;
		splitByDelimiter(line, '\t', split_line);
		printf("\n%f %d %lld %d %d", (float) cigar_frequency_index / (float) MAX_CIGAR_FREQ_SIZE, cigar_frequency_index, MAX_CIGAR_FREQ_SIZE, strlen(split_line[1]), number_of_records_processed);
		number_of_cigars = 1;
		if (isCommaInLine(split_line[1])) number_of_cigars = splitByDelimiter(split_line[1], ',', split_cigars);
		else strcpy(split_cigars[0], split_line[1]);
		if (cigar_frequency_index == 0)
		{
			for (i = 0; i < number_of_cigars; i++)
			{
				strcpy(cigar_frequency[cigar_frequency_index]->cigar, split_cigars[i]);
				cigar_frequency[cigar_frequency_index]->freq = 1;
				cigar_frequency_index++;
			}
		}
		else
		{
			for (j = 0; j < number_of_cigars; j++)
			{
				for (i = 0; i < cigar_frequency_index; i++)
				{
					if (strcmp(cigar_frequency[i]->cigar, split_cigars[j]) == 0)
					{
						cigar_frequency[i]->freq++;
						printf("\n%d Match found!! CIGAR: %s", cigar_frequency[i]->freq, cigar_frequency[i]->cigar);
						break;
					}
				}
				if (i == cigar_frequency_index) // Cigar is not present
				{
					strcpy(cigar_frequency[cigar_frequency_index]->cigar, split_cigars[j]);
					cigar_frequency[cigar_frequency_index]->freq = 1;
					cigar_frequency_index++;
				}
			}
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

	compress(input_filename, output_abridgefilename);
	return 0;
}

