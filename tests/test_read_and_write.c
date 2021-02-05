/*
 * test_read_and_write.c
 *
 *  Created on: Jan 16, 2021
 *      Author: sagnik
 */
# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>

int main(int argc, char *argv[])
{
	char inputfilename[1000];
	char outputfilename[1000];
	char *line; // for reading each line

	size_t len = 0;
	ssize_t line_len;
	FILE *fhr, *fhw;

	strcpy(inputfilename, argv[1]);
	strcpy(outputfilename, argv[2]);

	fhr = fopen(inputfilename, "r");
	if (fhr == NULL)
	{
		printf("Error! File not found");
		exit(1);
	}
	fhw = fopen(outputfilename, "w");
	if (fhw == NULL)
	{
		printf("File cannot be created");
		exit(1);
	}
	exit(1);

	while ((line_len = getline(&line, &len, fhr)) != -1)
	{
		continue;
		//fprintf(fhw, "%s", line);
	}

	fclose(fhr);
	fclose(fhw);
	return 0;
}
