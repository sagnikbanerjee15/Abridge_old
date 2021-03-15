# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

void addTagToFile(char *sam_inputfilename, char *sam_outputfilename, char *tag_name, char *tag_value)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	FILE *fhr, *fhw;

	size_t len = 0;
	ssize_t line_len;

	int i, j, k;

	char *line = NULL; // for reading each line
	char *tag_and_value;
	char *line_to_be_written_to_file;

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/

	tag_and_value = (char*) malloc(sizeof(char) * (strlen(tag_name) + strlen(tag_value) + 5));
	strcpy(tag_and_value, tag_name);
	strcat(tag_and_value, ":");
	strcat(tag_and_value, tag_value);

	line_to_be_written_to_file = (char*) malloc(sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE);

	fhr = fopen(sam_inputfilename, "r");
	if (fhr == NULL)
	{
		printf("\nError opening file %s for reading", sam_inputfilename);
		exit(1);
	}

	fhw = fopen(sam_outputfilename, "w");
	if (fhw == NULL)
	{
		printf("\nError opening file %s for writing", sam_outputfilename);
		exit(1);
	}
	/********************************************************************/
	while ((line_len = getline(&line, &len, fhr)) != -1)
	{
		if (line[0] == '@') fprintf(fhw, "%s", line);
		else break;
	}

	do
	{
		for (i = 0; line[i] != '\n'; i++)
			line_to_be_written_to_file[i] = line[i];
		line_to_be_written_to_file[i] = '\0';
		strcat(line_to_be_written_to_file, "\t");
		strcat(line_to_be_written_to_file, tag_and_value);
		strcat(line_to_be_written_to_file, "\n");
		fprintf(fhw, "%s", line_to_be_written_to_file);

	} while ((line_len = getline(&line, &len, fhr)) != -1);

	fclose(fhr);
	fclose(fhw);
}

int main(int argc, char *argv[])
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char *sam_inputfilename;
	char *sam_outputfilename;
	char *tag_name;
	char *tag_value;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	sam_inputfilename = (char*) malloc(sizeof(char) * 1000);
	sam_outputfilename = (char*) malloc(sizeof(char) * 1000);
	tag_name = (char*) malloc(sizeof(char) * 100);
	tag_value = (char*) malloc(sizeof(char) * 10);

	strcpy(sam_inputfilename, argv[1]);
	strcpy(sam_outputfilename, argv[2]);
	strcpy(tag_name, argv[3]);
	strcpy(tag_value, argv[4]);
	/********************************************************************/

	addTagToFile(sam_inputfilename, sam_outputfilename, tag_name, tag_value);
}
