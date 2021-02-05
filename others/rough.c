# include <stdio.h>
# include <string.h>
# include <time.h>

# define MAX_ITERATIONS 10000

int main()
{
	// Some absolute rubbish!!

	FILE *fhr;
	long long int i;
	long long int start_line = 5598448;
	long long int end_line = 5605621;

	long long int start_byte = 143212037;
	long long int end_byte = 150509071;
	char *buffer;
	size_t len = 0;
	ssize_t line_len;
	char *line = NULL; // for reading each line
	int line_num;

	buffer = (char*) malloc(sizeof(char) * 1000000);

	clock_t start, end;
	double cpu_time_used;
	char indexfilename[] = "/project/maizegdb/sagnik/ABRIDGE/analysis/developing_abridge/SRR13009993_STAR_Aligned.sortedByCoord.out.abridge.pass2.index";

	//File SEEK operations
	start = clock();
	fhr = fopen(indexfilename, "r");
	for (i = 0; i < MAX_ITERATIONS; i++)
	{
		fseek(fhr, start_byte, SEEK_SET);
		fread(buffer, end_byte - start_byte, 1, fhr);
		rewind(fhr);
	}
	fclose(fhr);
	end = clock();
	printf("\n File seek time: %lf", ((double) (end - start)) / CLOCKS_PER_SEC);

	start = clock();
	fhr = fopen(indexfilename, "r");
	for (i = 0; i < MAX_ITERATIONS; i++)
	{
		line_num = 0;
		buffer[0] = '\0';
		while ((line_len = getline(&line, &len, fhr)) != -1)
		{
			line_num++;
			if (line_num >= start_line && line_num <= end_line)
			{
				strcpy(buffer, line);
			}
		}
		rewind(fhr);
	}
	fclose(fhr);
	end = clock();
	printf("\n File read time: %lf", ((double) (end - start)) / CLOCKS_PER_SEC);

	return 0;
}
