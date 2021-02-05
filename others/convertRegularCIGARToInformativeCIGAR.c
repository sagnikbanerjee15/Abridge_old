# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include <time.h>

// ATGCN - ASCII Codes 33-37
char insert_characters[] =
{ '!', '"', '#', '$', '%' };
// ATGCN - ASCII Codes 38-42
char mismatch_characters[] =
{ '&', '\'', '(', ')', '*' };

struct Cigar_Items
{
	/*
	 * Template to store cigar items
	 */
	char def;
	int len;
};

void splitCigar(char *cigar, int *num_of_types, struct Cigar_Items *cigar_items_instance)
{
	/********************************************************************
	 * Variable declarations
	 ********************************************************************/
	//struct Cigar_Items **cigar_items_instance;
	int i, j = 0;
	int current_length = 0;
	/********************************************************************/

	//printf("\nInside splitCigar. CIGAR: %s", cigar);
	for (i = 0; cigar[i] != '\0'; i++)
	{
		if (isdigit(cigar[i]) != 0) //cigar[i] is a digit
		{
			current_length = current_length * 10 + cigar[i] - 48;
		}
		else if (isalpha(cigar[i]) != 0) //cigar[i] is a character
		{
			cigar_items_instance[j].def = cigar[i];
			cigar_items_instance[j].len = current_length;
			current_length = 0;
			j += 1;
		}
		else
		{
			//Should not come here
		}
	}
	*num_of_types = j;
}

void insertCharacterInString(char *str, int *str_length, char ins, int loc, int number_of_insertions_to_be_made)
{
	/*
	 * Inserts a character. Assumes that a large
	 */
	int i = (*str_length) + number_of_insertions_to_be_made;
	for (; i > loc; i--)
		str[i] = str[i - number_of_insertions_to_be_made];
	(*str_length) += number_of_insertions_to_be_made;
	while (number_of_insertions_to_be_made--)
		str[i++] = ins;
}

void replaceNucleotidesInSeqWithInsertSymbols(char *seq, int *seq_length, struct Cigar_Items *cigar_items_instance, int number_of_cigar_items)
{
	int i;
	int j;
	int sequence_index_to_be_replaced = 0;
	for (i = 0; i <= number_of_cigar_items; i++)
	{
		switch (cigar_items_instance[i].def)
		{
			case 'M':
				sequence_index_to_be_replaced += cigar_items_instance[i].len;
				break;
			case 'I':
				for (j = 0; j < cigar_items_instance[i].len; j++)
				{
					switch (seq[sequence_index_to_be_replaced])
					{
						case 'A':
							seq[sequence_index_to_be_replaced] = insert_characters[0];
							break;
						case 'T':
							seq[sequence_index_to_be_replaced] = insert_characters[1];
							break;
						case 'G':
							seq[sequence_index_to_be_replaced] = insert_characters[2];
							break;
						case 'C':
							seq[sequence_index_to_be_replaced] = insert_characters[3];
							break;
						case 'N':
							seq[sequence_index_to_be_replaced] = insert_characters[4];
							break;
					}
					sequence_index_to_be_replaced++;
				}
				sequence_index_to_be_replaced += cigar_items_instance[i].len;
				break;
		}
	}
}

void convertRegularCIGARToStringRepresentation(struct Cigar_Items *cigar_items_instance, int number_of_cigar_items, char *cigar_extended)
{
	int i;
	int j;
	int cigar_extended_index = 0;
	for (i = 0; i < number_of_cigar_items; i++)
	{
		switch (cigar_items_instance[i].def)
		{
			case 'M':
				for (j = 0; j < cigar_items_instance[i].len; j++)
					cigar_extended[cigar_extended_index++] = 'M';
				break;
			case 'I':
				for (j = 0; j < cigar_items_instance[i].len; j++)
					cigar_extended[cigar_extended_index++] = 'I';
				break;
		}
	}
}

void designInformativeCIGAR(char *md, char *seq, int *seq_length, struct Cigar_Items *cigar_items_instance, int number_of_cigar_items, char *cigar, char *cigar_extended, char *md_extended, char *icigar, char **splices)
{
	int i;
	int num;
	int j;
	int md_extended_index;
	int splice_number_index;
	int md_extended_length;
	/*char *cigar_extended;
	 char *md_extended;
	 char **splices;
	 char *icigar;*/
	int splices_index;
	char temp_convert_int_to_string[100];
	int print_outputs = 1;

	int intron_splice_index;
	/*static clock_t times[2];*/

	md_extended_index = 0;
	i = 0;
	num = 0;
	/*times[0] = clock();*/
	/*cigar_extended = (char*) malloc(sizeof(char) * (*seq_length * 2));
	 md_extended = (char*) malloc(sizeof(char) * (*seq_length * 2));
	 icigar = (char*) malloc(sizeof(char) * (*seq_length * 2));
	 splices = (char**) malloc(sizeof(char*) * 100);
	 for (i = 0; i < 100; i++)
	 splices[i] = (char*) malloc(sizeof(char) * 100);*/
	/*times[0] = clock() - times[0];*/

	splitCigar(cigar, &number_of_cigar_items, cigar_items_instance);
	if (print_outputs) printf("\n%s", seq);
	replaceNucleotidesInSeqWithInsertSymbols(seq, seq_length, cigar_items_instance, number_of_cigar_items);
	if (print_outputs) printf("\n%s", seq);
	convertRegularCIGARToStringRepresentation(cigar_items_instance, number_of_cigar_items, cigar_extended);

	/*
	 * Expand the md string
	 */
	i = 0;
	while (md[i] != '\0')
	{
		if (isdigit(md[i]) != 0) // md[i] is a digit
		num = num * 10 + md[i] - 48;
		else if (md[i] == '^')
		{
			for (j = 0; j < num; j++)
				md_extended[md_extended_index++] = '=';
			num = 0;
			i += 1;
			while (isdigit(md[i]) == 0) // Iterate till you pick up all the deletions and ignore them for now
				i += 1;
			i -= 1;
		}
		else
		{
			for (j = 0; j < num; j++)
				md_extended[md_extended_index++] = '=';
			md_extended[md_extended_index++] = 'X';
			num = 0;
		}
		i += 1;
	}
	for (j = 0; j < num; j++)
		md_extended[md_extended_index++] = '=';
	md_extended[md_extended_index++] = '\0';
	md_extended_length = md_extended_index - 1;

	if (print_outputs) printf("\n%s %d", md_extended, md_extended_length);

	/*
	 * Add insert symbols to md_extended
	 */
	for (i = 0; i < *seq_length; i++)
		if (seq[i] >= 33 && seq[i] <= 36) insertCharacterInString(md_extended, &md_extended_length, seq[i], i, 1);
	/*printf("\n%s %d", md_extended, md_extended_length);*/

	/*
	 * Add mismatch characters to md_extended
	 */

	for (i = 0; i < md_extended_length; i++)
	{
		if (md_extended[i] == 'X')
		{
			switch (seq[i])
			{
				case 'A':
					md_extended[i] = mismatch_characters[0];
					break;
				case 'T':
					md_extended[i] = mismatch_characters[1];
					break;
				case 'G':
					md_extended[i] = mismatch_characters[2];
					break;
				case 'C':
					md_extended[i] = mismatch_characters[3];
					break;
			}
		}
	}
	if (print_outputs) printf("\n%s %d", md_extended, md_extended_length);

	/*
	 * Encode deletions and splice junction info in md_extended
	 */
	md_extended_index = 0;
	splice_number_index = 0;
	for (i = 0; i < number_of_cigar_items; i++)
	{
		switch (cigar_items_instance[i].def)
		{
			case 'M':
				md_extended_index += cigar_items_instance[i].len;
				break;
			case 'N':
				md_extended[md_extended_index] = splice_number_index++ + 48;
				break;
			case 'I':
				md_extended_index += cigar_items_instance[i].len;
				break;
			case 'D':
				insertCharacterInString(md_extended, &md_extended_length, '-', md_extended_index, cigar_items_instance[i].len);
				md_extended_index += cigar_items_instance[i].len;
				break;
		}
	}
	if (print_outputs) printf("\n%s %d", md_extended, md_extended_length);

	/*
	 * Generate the informative CIGAR
	 */

	/*
	 * Combine the splice items into an array
	 */
	splices_index = 0;
	for (i = 0; i < number_of_cigar_items; i++)
	{
		if (cigar_items_instance[i].def == 'N')
		{
			sprintf(temp_convert_int_to_string, "%d", cigar_items_instance[i].len);
			strcpy(splices[splices_index], temp_convert_int_to_string);
			strcat(splices[splices_index], "N");
			splices_index++;
		}
	}
	num = 0;
	icigar[0] = '\0';
	intron_splice_index = 0;
	for (i = 0; i < md_extended_length; i++)
	{
		if (md_extended[i] == '=') num++;
		else
		{
			sprintf(temp_convert_int_to_string, "%d", num);
			if (num > 0)
			{
				strcat(icigar, temp_convert_int_to_string);
				strcat(icigar, "M");
			}
			num = 0;
			if ((md_extended[i] >= 33 && md_extended[i] <= 37) || (md_extended[i] >= 38 && md_extended[i] <= 42)) // Either a mismatch or an insertion
			{
				temp_convert_int_to_string[1] = '\0';
				temp_convert_int_to_string[0] = md_extended[i];
				strcat(icigar, temp_convert_int_to_string);
			}
			else if (md_extended[i] >= '0' && md_extended[i] <= '9')
			{
				num = 0;
				intron_splice_index = 0;
				while (md_extended[i] >= '0' && md_extended[i] <= '9')
				{
					intron_splice_index = intron_splice_index * 10 + md_extended[i] - 48;
					i++;
					num++;
				}
				strcat(icigar, splices[intron_splice_index]);
				i--;
			}
			else if (md_extended[i] == '-') // Deletion
			{
				while (md_extended[i] == '-')
				{
					num++;
					i++;
				}
				i--;
				sprintf(temp_convert_int_to_string, "%d", num);
				strcat(icigar, temp_convert_int_to_string);
				strcat(icigar, "D");
				num = 0;
			}
		}
	}
	sprintf(temp_convert_int_to_string, "%d", num);
	strcat(icigar, temp_convert_int_to_string);
	strcat(icigar, "M");

	/*times[1] = clock();*/
	if (print_outputs) printf("\n%s", icigar);
	/*free(cigar_extended);
	 free(md_extended);
	 for (i = 0; i < 100; i++)
	 free(splices[i]);
	 free(splices);
	 free(icigar);*/
	/*times[1] = clock() - times[1];*/
	/*return times;*/
}

/*
 *
 * 26M1I16M620N41M1D16M
 * 28A49A4^A16
 * AGAGGTGGCCACCATGACGAGGAGATAGAGGCGTTCAAAGACACAATATGAAGAGAGCAGATGAGAAATTTGTTGAGGAGATACAAAAAAACGTAACGCC
 *
 * MMMMMMMMMMMMMMMMMMMMMMMMMMIMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM-MMMMMMMMMMMMMMMM
 * AGAGGTGGCCACCATGACGAGGAGATAGAGGCGTTCAAAGACACAATATGAAGAGAGCAGATGAGAAATTTGTTGAGGAGATAC-AAAAAAACGTAACGCC
 * ==========================I==X=================================================X====-================
 *
 * MMMMMMMMMMMMMMMMMMMMMMMMMMIMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
 * AGAGGTGGCCACCATGACGAGGAGATAGAGGCGTTCAAAGACACAATATGAAGAGAGCAGATGAGAAATTTGTTGAGGAGATACAAAAAAACGTAACGCC
 * AGAGGTGGCCACCATGACGAGGAGAT!GAGGCGTTCAAAGACACAATATGAAGAGAGCAGATGAGAAATTTGTTGAGGAGATACAAAAAAACGTAACGCC
 * ==========================I==X=================================================X====================
 * ==========================!==*=================================================*==================== Replace with special characters
 * ==========================!================1===================================*====-================
 *
 * MMMMMMMMMMMMMMMMMMMMMMMMMMIMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
 * AGAGGTGGCCACCATGACGAGGAGAT!GAGGCGTTCAAAGACACAATATGAAGAGAGCAGATGAGAAATTTGTTGAGGAGATACAAAAAAACGTAACGCC
 * ============================X=================================================X====================
 * ==========================!==X=================================================X====================
 * ==========================!==*=================================================*====================
 *
 *
 * TGAAAAAAAAAAAAAAAGTTTGGGTTTTTAAAAGAAGTCTTCTTTTTTTTGGGTTTATGGGCTGAAAAAAATTGTTGTGTTTTTTGTATTTTATTTTTT
 * TGAAAAAAAAAAAAAAAGTTTGGGTTTTTAAAAGAAGTCTTCTTTTTTTTGGGTTTATGGGCTGAAA!AAATTGTTGTGTTTTTTGTATTTTATTTTTT
 * ======================================X==X========================================================
 *
 *

 */

//=========================================0=============
int main()
{
	char *regular_cigar = "100M";
	char seq[2500] = "CCCGGAAACCCTTGAAATCGGACGGTTTAGTGAAAATGGAGGATCAAGTTGGGTTTGGGTTCCGTCCGAACGACGAGGAGCTCGTTGGTCACTATCTCCG";
	char *md = "10A89";
	struct Cigar_Items cigar_items_instance[100];
	int number_of_cigar_items = 0;
	int seq_length;
	int i;

	char *cigar_extended;
	char *md_extended;
	char *icigar;
	char **splices;

	seq_length = strlen(seq);

	cigar_extended = (char*) malloc(sizeof(char) * (seq_length * 2));
	md_extended = (char*) malloc(sizeof(char) * (seq_length * 2));
	icigar = (char*) malloc(sizeof(char) * (seq_length * 2));
	splices = (char**) malloc(sizeof(char*) * 100);
	for (i = 0; i < 100; i++)
		splices[i] = (char*) malloc(sizeof(char) * 100);

	/*printf("\n Seq length: %i", seq_length);
	 fflush(stdout);*/
	/*expanded_md = (char*) malloc(sizeof(char) * (seq_length * 2));*/
	/*printf("\nRegular CIGAR: %s", regular_cigar);*/

	/*
	 printf("\nRegular CIGAR: %s", regular_cigar);
	 fflush(stdout);
	 return 0;
	 */

	designInformativeCIGAR(md, seq, &seq_length, cigar_items_instance, number_of_cigar_items, regular_cigar, cigar_extended, md_extended, icigar, splices);
	printf("\n");
	return 0;

}
