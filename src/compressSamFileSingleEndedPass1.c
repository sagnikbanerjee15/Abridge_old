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

void encodeByRunLength(struct Compressed_DS **compressed_ds_pool, int compressed_ds_pool_total, char *encoded_string)
{
	int i;
	int count;
	char str[1000];

	str[0] = '\0';
	for (i = 0; i < compressed_ds_pool_total; i++)
	{
		count = 1;
		while (i < compressed_ds_pool_total - 1 && compressed_ds_pool[i]->num_reads == compressed_ds_pool[i + 1]->num_reads)
		{
			count++;
			i++;
		}
		sprintf(str, "%d", count);
		strcat(encoded_string, str);
		strcat(encoded_string, "-");
		//printf("\n%s", encoded_string);
		sprintf(str, "%ld", compressed_ds_pool[i]->num_reads);
		strcat(encoded_string, str);
		strcat(encoded_string, ",");
		//printf("\n%s", encoded_string);
	}
	//printf("\n Encoded String: %s", encoded_string);
}

void writeToFile(FILE *fhw, struct Compressed_DS **compressed_ds_pool, int compressed_ds_pool_total, char *write_to_file_col1, char *write_to_file_col2, char *write_to_file_col3, char *encoded_string)
{
	int i;
	char str[1000];
	char line_to_be_written_to_file[MAX_LINE_TO_BE_WRITTEN_TO_FILE];

	//printf("\nInside writeToFile compressed_ds_pool_total: %d", compressed_ds_pool_total);
	for (i = 0; i < compressed_ds_pool_total; i++)
	{
		if (i == 0)
		{
			if (compressed_ds_pool[i]->position != 1) sprintf(str, "%lld", compressed_ds_pool[i]->position);
			else str[0] = '\0'; // empty string
			strcat(write_to_file_col1, str);
		}

		//strcat(write_to_file_col1, ",");

		strcat(write_to_file_col2, compressed_ds_pool[i]->icigar);
		strcat(write_to_file_col2, ",");

		sprintf(str, "%ld", compressed_ds_pool[i]->num_reads);
		strcat(write_to_file_col3, str);
		strcat(write_to_file_col3, ",");

	}

	//write_to_file_col1[strlen(write_to_file_col1) - 1] = '\0'; // Removing the last comma
	write_to_file_col2[strlen(write_to_file_col2) - 1] = '\0'; // Removing the last comma
	write_to_file_col3[strlen(write_to_file_col3) - 1] = '\0'; // Removing the last comma

	strcpy(line_to_be_written_to_file, write_to_file_col1);
	strcat(line_to_be_written_to_file, "\t");
	strcat(line_to_be_written_to_file, write_to_file_col2);
	strcat(line_to_be_written_to_file, "\t");
	if (compressed_ds_pool_total == 1)
	{
		strcat(line_to_be_written_to_file, write_to_file_col3);
	}
	else
	{
		strcat(line_to_be_written_to_file, write_to_file_col3);
		/*
		 encoded_string[0] = '\0';
		 encodeByRunLength(compressed_ds_pool, compressed_ds_pool_total, encoded_string);
		 strcat(line_to_be_written_to_file, encoded_string);
		 */
	}
	strcat(line_to_be_written_to_file, "\n");
	fprintf(fhw, "%s", line_to_be_written_to_file);

	// Reinitialize for next iteration
	write_to_file_col1[0] = '\0';
	write_to_file_col2[0] = '\0';
	write_to_file_col3[0] = '\0';

}

void readAlignmentsAndCompress(char *input_samfilename, char *output_abridgefilename, char *reference_id, short int flag_ignore_soft_clippings, short int flag_ignore_mismatches, short int flag_ignore_unmapped_sequences)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	FILE *fhr;
	FILE *fhw;

	char **split_line; // List of strings to store each element of a single alignment
	char **split_tags; // List of strings to store tag information
	char **split_reference_info;
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

	size_t len = 0;
	ssize_t line_len;

	int flag;
	int i, j, k; // Required in loops
	int number_of_tags;
	int number_of_fields; // Number of fields in each sam alignment entry
	int sam_tag_index;
	int tab_number;
	int num_items_in_alignment_pool = 0; // Items in pool
	int samflag_quick_read_index = 0;
	int compressed_ds_pool_index = 0;
	int number_of_reference_sequences = 0;
	int reference_sequence_index = 0;

	long long int relative_position_to_previous_read_cluster;
	long long int previous_position = -1;
	long long int current_position;
	long long int number_of_records_written = 0;
	long long int number_of_records_read = 0;
	long long int num_pools_written = 0;

	struct Sam_Alignment *prev_alignment;
	struct Sam_Alignment *curr_alignment;
	struct Sam_Alignment *temp_alignment;
	struct Sam_Alignment **alignment_pool_same_position;
	struct Compressed_DS **compressed_ds_pool;
	struct Reference_Sequence_Info **reference_info;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr = fopen(input_samfilename, "r");
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

	split_tags = (char**) malloc(sizeof(char*) * ROWS);
	for (i = 0; i < ROWS; i++)
		split_tags[i] = (char*) malloc(sizeof(char) * COLS);

	split_reference_info = (char**) malloc(sizeof(char*) * ROWS);
	for (i = 0; i < ROWS; i++)
		split_reference_info[i] = (char*) malloc(sizeof(char) * COLS);

	compressed_ds_pool = (struct Compressed_DS**) malloc(sizeof(struct Compressed_DS*) * MAX_POOL_SIZE);
	for (i = 0; i < MAX_POOL_SIZE; i++)
		compressed_ds_pool[i] = allocateMemoryCompressed_DS();

	write_to_file_col1 = (char*) malloc(sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE);
	write_to_file_col2 = (char*) malloc(sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE);
	write_to_file_col3 = (char*) malloc(sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE);
	encoded_string = (char*) malloc(sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE);
	write_to_file_col1[0] = '\0';
	write_to_file_col2[0] = '\0';
	write_to_file_col3[0] = '\0';
	reference_id_quick_read = (char*) malloc(sizeof(char) * 1000);
	samflag_quick_read = (char*) malloc(sizeof(char) * 1000);
	prev_reference_name = (char*) malloc(sizeof(char) * 1000);
	prev_reference_name[0] = '\0';
	curr_reference_name = (char*) malloc(sizeof(char) * 1000);
	curr_reference_name[0] = '\0';

	curr_alignment = allocateMemorySam_Alignment();
	prev_alignment = allocateMemorySam_Alignment();
	temp_alignment = allocateMemorySam_Alignment();
	reference_info = (struct Reference_Sequence_Info**) malloc(sizeof(struct Reference_Sequence_Info*) * MAX_REFERENCE_SEQUENCES);
	for (i = 0; i < MAX_REFERENCE_SEQUENCES; i++)
		reference_info[i] = allocateMemoryReference_Sequence_Info();

	printf("\n Starting to read file");
	fflush(stdout);
	/********************************************************************/

	/*
	 * Read in the reference sequence information
	 */
	while ((line_len = getline(&line, &len, fhr)) != -1)
	{
		if (line[0] == '@')
		{
			if (line[1] == 'S' && line[2] == 'Q')
			{
				//printf("\n Reference: %s %d", line, strlen(line));
				//fflush(stdout);
				strcpy(reference_info[number_of_reference_sequences]->line, line);
				number_of_reference_sequences++;
			}
		}
		else break;
	}

	do
	{
		number_of_records_read += 1;
		if (number_of_records_read % 1000000 == 0)
		{
			printf("\nNumber of records read: %lld Million %s", number_of_records_read / 1000000, output_abridgefilename);
			fflush(stdout);
		}
		number_of_fields = splitByDelimiter(line, '\t', split_line);
		populateSamAlignmentInstance(curr_alignment, split_line, number_of_fields, split_tags);
		strcpy(curr_reference_name, curr_alignment->reference_name);
		if (curr_alignment->samflag == 4) continue;
		current_position = curr_alignment->start_position;
		//printSamAlignmentInstance(curr_alignment);
		generateIntegratedCigar(curr_alignment, flag_ignore_soft_clippings, flag_ignore_mismatches, flag_ignore_unmapped_sequences);
		//printf("\n Position:%lld iCIGAR: %s", curr_alignment->start_position, curr_alignment->icigar);
		if (strlen(prev_reference_name) == 0) // 1st chromosome - initialize stuffs
		{
			previous_position = current_position;
			strcpy(prev_reference_name, curr_reference_name);
			strcpy(compressed_ds_pool[compressed_ds_pool_index]->icigar, curr_alignment->icigar);
			compressed_ds_pool[compressed_ds_pool_index]->num_reads = 1;
			compressed_ds_pool[compressed_ds_pool_index]->position = curr_alignment->start_position;
			compressed_ds_pool_index++;
			//printf("\n Writing Reference to file %s %d", reference_info[reference_sequence_index]->line, reference_sequence_index);
			//fflush(stdout);
			fprintf(fhw, "%s", reference_info[reference_sequence_index]->line);
			reference_sequence_index++;
		}

		else if (strcmp(prev_reference_name, curr_reference_name) != 0) // New chromosome
		{
			writeToFile(fhw, compressed_ds_pool, compressed_ds_pool_index, write_to_file_col1, write_to_file_col2, write_to_file_col3, encoded_string);
			compressed_ds_pool_index = 0;
			previous_position = current_position;
			strcpy(prev_reference_name, curr_reference_name);
			strcpy(compressed_ds_pool[compressed_ds_pool_index]->icigar, curr_alignment->icigar);
			compressed_ds_pool[compressed_ds_pool_index]->num_reads = 1;
			compressed_ds_pool[compressed_ds_pool_index]->position = curr_alignment->start_position;
			compressed_ds_pool_index++;
			//printf("\n Writing Reference to file %s %d", reference_info[reference_sequence_index]->line, reference_sequence_index);
			//fflush(stdout);
			fprintf(fhw, "%s", reference_info[reference_sequence_index]->line);
			reference_sequence_index++;
		}
		else // Same chromosome
		{
			if (previous_position == current_position)
			{
				for (i = 0; i < compressed_ds_pool_index; i++)
				{
					if (strcmp(compressed_ds_pool[i]->icigar, curr_alignment->icigar) == 0)
					{
						compressed_ds_pool[i]->num_reads++;
						break;
					}
				}
				if (i == compressed_ds_pool_index) // New icigar encountered
				{
					strcpy(compressed_ds_pool[compressed_ds_pool_index]->icigar, curr_alignment->icigar);
					compressed_ds_pool[compressed_ds_pool_index]->num_reads = 1;
					compressed_ds_pool[compressed_ds_pool_index]->position = compressed_ds_pool[0]->position;
					compressed_ds_pool_index++;
				}
			}
			else
			{
				writeToFile(fhw, compressed_ds_pool, compressed_ds_pool_index, write_to_file_col1, write_to_file_col2, write_to_file_col3, encoded_string);
				compressed_ds_pool_index = 0;
				strcpy(compressed_ds_pool[compressed_ds_pool_index]->icigar, curr_alignment->icigar);
				compressed_ds_pool[compressed_ds_pool_index]->num_reads = 1;
				compressed_ds_pool[compressed_ds_pool_index]->position = curr_alignment->start_position - previous_position;
				compressed_ds_pool_index++;
			}
			previous_position = current_position;
		}
	} while ((line_len = getline(&line, &len, fhr)) != -1);
	//Write final data to file
	writeToFile(fhw, compressed_ds_pool, compressed_ds_pool_index, write_to_file_col1, write_to_file_col2, write_to_file_col3, encoded_string);
	/*printf("\n Execution complete");
	 fflush(stdout);*/
	/*for (i = 0; i < MAX_POOL_SIZE; i++)
	 freeSamAlignmentInstance(alignment_pool_same_position[i]);
	 free(alignment_pool_same_position);*/
	printf("\n");
}

int main(int argc, char *argv[])
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char input_samfilename[FILENAME_LENGTH];
	char output_abridgefilename[FILENAME_LENGTH];
	char reference_id[FILENAME_LENGTH];
	char *temp; //Required for strtoi

	short int flag_ignore_soft_clippings;
	short int flag_ignore_mismatches;
	short int flag_ignore_sequence_information;
	short int flag_ignore_unmapped_sequences;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	strcpy(reference_id, argv[1]);
	flag_ignore_soft_clippings = strtol(argv[2], &temp, 10);
	flag_ignore_mismatches = strtol(argv[3], &temp, 10);
	flag_ignore_sequence_information = strtol(argv[4], &temp, 10);
	flag_ignore_unmapped_sequences = strtol(argv[5], &temp, 10);
	strcpy(input_samfilename, argv[6]);
	strcpy(output_abridgefilename, argv[7]);
	/********************************************************************/

	/*
	 * If user requests no sequence information then everything else is also ignored
	 */
	if (flag_ignore_sequence_information == 1)
	{
		flag_ignore_mismatches = 0;
		flag_ignore_soft_clippings = 0;
		flag_ignore_unmapped_sequences = 0;
	}
	readAlignmentsAndCompress(input_samfilename, output_abridgefilename, reference_id, flag_ignore_soft_clippings, flag_ignore_mismatches, flag_ignore_unmapped_sequences);
	return 0;
}
