# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include <math.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

void readInTheEntireGenome(char *genome_filename, struct Whole_Genome_Sequence *whole_genome)
{
	FILE *fhr;
	char *buffer = NULL;
	size_t len = 0;
	ssize_t line_len;
	int i;

	//buffer = (char*) malloc(sizeof(char) * pow(2, 32));
	fhr = fopen(genome_filename, "rb");
	if (fhr == NULL)
	{
		printf("Error! File not found");
		exit(1);
	}

	whole_genome->nucleotides = (char**) malloc(sizeof(char*) * MAX_REFERENCE_SEQUENCES);
	whole_genome->reference_sequence_name = (char**) malloc(sizeof(char*) * MAX_REFERENCE_SEQUENCES);
	whole_genome->reference_sequence_length = (unsigned long long int*) malloc(sizeof(unsigned long long int) * MAX_REFERENCE_SEQUENCES);
	whole_genome->number_of_reference_sequences = 0;

	while ((line_len = getline(&buffer, &len, fhr)) != -1)
	{
		if (buffer[0] == '>')
		{
			whole_genome->reference_sequence_name[whole_genome->number_of_reference_sequences] = (char*) malloc(sizeof(char) * (len + 1));
			for (i = 0; buffer[i] != '\0'; i++)
				if (buffer[i] == 32)
				{
					buffer[i] = '\0';
					break;
				}
			for (i = 1; buffer[i] != '\0'; i++)
				buffer[i - 1] = buffer[i];
			buffer[i - 1] = buffer[i];
			strcpy(whole_genome->reference_sequence_name[whole_genome->number_of_reference_sequences], buffer);
		}
		else
		{
			whole_genome->nucleotides[whole_genome->number_of_reference_sequences] = (char*) malloc(sizeof(char) * (len + 1));
			strcpy(whole_genome->nucleotides[whole_genome->number_of_reference_sequences], buffer);
			whole_genome->reference_sequence_length[whole_genome->number_of_reference_sequences] = strlen(buffer);
			whole_genome->number_of_reference_sequences++;
		}
		//free(buffer);
		//buffer = NULL;
	}
	/*
	 for (i = 0; i < whole_genome->number_of_reference_sequences; i++)
	 printf("\n %d %s %lld %lld", i + 1, whole_genome->reference_sequence_name[i], strlen(whole_genome->nucleotides[i]), whole_genome->reference_sequence_length[i]);
	 */
	fclose(fhr);
}

void readAbridgeIndex(struct Abridge_Index *abridge_index, char *abridge_index_filename, char **split_line, short int *flag_ignore_mismatches, short int *flag_ignore_soft_clippings, short int *flag_ignore_unmapped_sequences, short int *flag_ignore_quality_score)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	FILE *fhr;

	unsigned long long int line_num = 0;

	size_t len = 0;
	ssize_t line_len;

	char *temp; //Useless
	char *line = NULL; // for reading each line
	char str[100];

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr = fopen(abridge_index_filename, "r");
	if (fhr == NULL)
	{
		printf("Error! File not found");
		exit(1);
	}
	/********************************************************************/
	while ((line_len = getline(&line, &len, fhr)) != -1)
	{
		line_num++;
		splitByDelimiter(line, '\t', split_line);
		if (line_num == 1)
		{
			*flag_ignore_mismatches = strtol(split_line[0], &temp, 10);
			*flag_ignore_soft_clippings = strtol(split_line[1], &temp, 10);
			*flag_ignore_unmapped_sequences = strtol(split_line[2], &temp, 10);
			*flag_ignore_quality_score = strtol(split_line[3], &temp, 10);
			continue;
		}
		//printf("\n%d %s", abridge_index->number_of_items, split_line[0]);
		//fflush(stdout);
		strcpy(abridge_index->chromosome[abridge_index->number_of_items], split_line[0]);
		abridge_index->start[abridge_index->number_of_items] = strtol(split_line[1], &temp, 10);
		abridge_index->end[abridge_index->number_of_items] = strtol(split_line[2], &temp, 10);
		abridge_index->start_byte[abridge_index->number_of_items] = strtol(split_line[3], &temp, 10);
		abridge_index->end_byte[abridge_index->number_of_items] = strtol(split_line[4], &temp, 10);
		abridge_index->number_of_items++;
		//printf("\n %lld", MAX_POOL_SIZE - abridge_index->number_of_items);
	}
	fclose(fhr);
}

int maxClusterSize(struct Abridge_Index *abridge_index)
{
	long long int i;
	long long int max_cluster_size = -1;
	for (i = 0; i < abridge_index->number_of_items; i++)
		if (max_cluster_size < (abridge_index->end[i] - abridge_index->start[i])) max_cluster_size = (abridge_index->end[i] - abridge_index->start[i]);
	return max_cluster_size;

}

void loadSequencesFromFile(char **sequences, char *filename)
{
	FILE *fhr;
	long long int num_sequence = 0;
	size_t len = 0;
	ssize_t line_len;
	char *line = NULL; // for reading each line

	fhr = fopen(filename, "rb");
	if (fhr == NULL)
	{
		printf("Error! File not found");
		exit(1);
	}
	while ((line_len = getline(&line, &len, fhr)) != -1)
		if (line[0] != '>') strcpy(sequences[num_sequence++], line);

}

long long int reverseInteger(long long int val)
{
	long long int temp = 0;
	while (val)
	{
		temp = temp * 10 + (val % 10);
		val /= 10;
	}
	return temp;
}

long long int extractNHfromicigar(char *icigar, int icigar_length)
{
	char nh_val_string[15];
	char nh_val_string_rev[15];
	long long int nh_val = 0;
	int i;
	int j;
	int k;
	char *t;

	j = 0;
	for (i = icigar_length - 1; i >= 0; i--)
	{
		if (isdigit(icigar[i]) == 0) break;
		//nh_val = nh_val * 10 + icigar[i] - 48;
		nh_val_string[j++] = icigar[i];
	}
	icigar[i + 1] = '\0';
	nh_val_string[j] = '\0';

	k = j - 1;
	for (i = 0; i < j; i++)
		nh_val_string_rev[i] = nh_val_string[k--];
	nh_val_string_rev[j] = '\0';
	nh_val = strtol(nh_val_string_rev, &t, 10);
	return nh_val;
}

void replaceCharacterInString(char *str, char ch_to_be_replaced, char replace_with, int str_length)
{
	int i;
	for (i = 0; i < str_length; i++)
		if (str[i] == ch_to_be_replaced) str[i] = replace_with;
}

int findSamFormatFlag(char *icigar, int icigar_length, char *XS)
{
	int i;
	int samformatflag;

	for (i = 0; i < icigar_length; i++)
	{
		switch (icigar[i])
		{
			case 'B':
				//replaceCharacterInString(icigar, 'B', 'M', icigar_length);
				samformatflag = 0;
				XS[0] = '.';
				break;
			case 'E':
				//replaceCharacterInString(icigar, 'E', 'M', icigar_length);
				samformatflag = 16;
				XS[0] = '.';
				break;
			case 'F':
				//replaceCharacterInString(icigar, 'F', 'M', icigar_length);
				samformatflag = 256;
				XS[0] = '.';
				break;
			case 'H':
				//replaceCharacterInString(icigar, 'H', 'M', icigar_length);
				samformatflag = 272;
				XS[0] = '.';
				break;
			case 'J':
				//replaceCharacterInString(icigar, 'J', 'M', icigar_length);
				samformatflag = 0;
				XS[0] = '+';
				break;
			case 'K':
				//replaceCharacterInString(icigar, 'K', 'M', icigar_length);
				samformatflag = 16;
				XS[0] = '+';
				break;
			case 'L':
				//replaceCharacterInString(icigar, 'L', 'M', icigar_length);
				samformatflag = 256;
				XS[0] = '+';
				break;
			case 'O':
				//replaceCharacterInString(icigar, 'O', 'M', icigar_length);
				samformatflag = 272;
				XS[0] = '+';
				break;
			case 'P':
				//replaceCharacterInString(icigar, 'P', 'M', icigar_length);
				samformatflag = 0;
				XS[0] = '-';
				break;
			case 'Q':
				//replaceCharacterInString(icigar, 'Q', 'M', icigar_length);
				samformatflag = 16;
				XS[0] = '-';
				break;
			case 'R':
				//replaceCharacterInString(icigar, 'R', 'M', icigar_length);
				samformatflag = 256;
				XS[0] = '-';
				break;
			case 'U':
				//replaceCharacterInString(icigar, 'U', 'M', icigar_length);
				samformatflag = 272;
				XS[0] = '-';
				break;
		}
	}
	return samformatflag;
}

void generateReadSequenceAndMDString(struct Sam_Alignment *sam_alignment_instance, struct Whole_Genome_Sequence *whole_genome)
{
	int i;
	int j;
	int k;
	char read_from_genome[MAX_SEQ_LEN];
	char MD[MAX_SEQ_LEN];

	for (i = 0; i < sam_alignment_instance->number_of_cigar_items; i++)
		if (sam_alignment_instance->cigar_items[i].def == 'N') break;

	if (i != sam_alignment_instance->number_of_cigar_items) //Spliced read
	{

	}
	else // Not spliced read
	{
		if (sam_alignment_instance->samflag == 0 || sam_alignment_instance->samflag == 256) // Forward read - simply pick up from genome
		{
			for (i = 0; i < whole_genome->number_of_reference_sequences; i++)
			{
				if (strcmp(sam_alignment_instance->reference_name, whole_genome->reference_sequence_name[i]) == 0)
				{
					int start = sam_alignment_instance->start_position - 1;
					int end = start + strlen(sam_alignment_instance->seq);
					k = 0;
					for (j = start; j <= end; j++)
						read_from_genome[k++] = whole_genome->nucleotides[i][j];
					read_from_genome[k] = '\0';
					break;
				}
			}
		}
	}
}

void convertIcigarToCigarandMD(struct Whole_Genome_Sequence *whole_genome, struct Sam_Alignment *sam_alignment_instance, char *chromosome, short int flag_ignore_mismatches, short int flag_ignore_soft_clippings, short int flag_ignore_unmapped_sequences, short int flag_ignore_quality_score, short int flag_ignore_sequence_information, char *default_quality_value)
{
	/*
	 * Coverts the iCIGAR into CIGAR and MD
	 * Returns the samformatflag
	 */

	int samformatflag = -1;
	int NH_value;
	int icigar_length;
	int i;
	int j;
	int processing_left_soft_clip = 1;
	int left_soft_clip_index = 0;
	int right_soft_clip_index = 0;
	int cigar_index = 0;
	int num = 0;
	int cigar_temp_index = 0;
	int cigar_items_instance_index;
	int flag = 0;
	int soft_clips_removed_qual_index;
	int soft_clips_removed_seq_index;

	char XS[2];

	struct Cigar_Items cigar_items_instance[MAX_SEQ_LEN];

	char temp[100];
	char cigar_temp[MAX_SEQ_LEN];

//printf("\niCIGAR at the beginning %s ", sam_alignment_instance->icigar);

	icigar_length = strlen(sam_alignment_instance->icigar);
	NH_value = extractNHfromicigar(sam_alignment_instance->icigar, icigar_length);

	XS[1] = '\0';
	samformatflag = findSamFormatFlag(sam_alignment_instance->icigar, icigar_length, XS);
	sam_alignment_instance->samflag = samformatflag;
	/*
	 * Construct the Cigar string, MD string, Soft Clips, etc.
	 */
	splitCigar(sam_alignment_instance->icigar, &cigar_items_instance_index, cigar_items_instance);
//printf("\nICIGAR %s Number of items: %d", sam_alignment_instance->icigar, cigar_items_instance_index);

	/*
	 for (i = 0; i < cigar_items_instance_index; i++)
	 printf("\n%d %c ", cigar_items_instance[i].len, cigar_items_instance[i].def);
	 fflush(stdout);
	 */
	strcpy(sam_alignment_instance->reference_name, chromosome);
	sam_alignment_instance->cigar[0] = '\0';
	sam_alignment_instance->temp[0] = '\0';
	sam_alignment_instance->soft_clippings.left[0] = '\0';
	sam_alignment_instance->soft_clippings.right[0] = '\0';
	sam_alignment_instance->soft_clippings.left_qual[0] = '\0';
	sam_alignment_instance->soft_clippings.right_qual[0] = '\0';
	sam_alignment_instance->soft_clips_removed_qual[0] = '\0';
	sam_alignment_instance->soft_clips_removed_seq[0] = '\0';
	soft_clips_removed_qual_index = 0;
	soft_clips_removed_seq_index = 0;

	strcpy(sam_alignment_instance->tags[0].name, "NH");
	strcpy(sam_alignment_instance->tags[0].type, "i");
	sprintf(temp, "%d", NH_value);
	strcpy(sam_alignment_instance->tags[0].val, temp);

	strcpy(sam_alignment_instance->tags[1].name, "XS");
	strcpy(sam_alignment_instance->tags[1].type, ".");
	strcpy(sam_alignment_instance->tags[1].val, XS);

	strcpy(sam_alignment_instance->tags[2].name, "MD");
	strcpy(sam_alignment_instance->tags[2].type, "Z");
	strcpy(sam_alignment_instance->tags[2].val, "dummy");

	/*
	 for (i = 0; i < cigar_items_instance_index; i++)
	 printf("\n%d %c", cigar_items_instance[i].len, cigar_items_instance[i].def);
	 */
	for (i = 0; i < cigar_items_instance_index; i++)
	{
		if (processing_left_soft_clip == 1 && isCharacterInString("atgcn", cigar_items_instance[i].def))
		{
			sam_alignment_instance->soft_clippings.left[left_soft_clip_index] = cigar_items_instance[i].def - 32;
			sam_alignment_instance->soft_clippings.left[left_soft_clip_index + 1] = '\0';
			if (flag_ignore_quality_score == 0)
			{
				i++;
				sam_alignment_instance->soft_clippings.left_qual[left_soft_clip_index] = cigar_items_instance[i].def - 90;
				sam_alignment_instance->soft_clippings.left_qual[left_soft_clip_index + 1] = '\0';
			}
			left_soft_clip_index++;
		}
		else if (processing_left_soft_clip == 0 && isCharacterInString("atgcn", cigar_items_instance[i].def))
		{
			sam_alignment_instance->soft_clippings.right[right_soft_clip_index] = cigar_items_instance[i].def - 32;
			sam_alignment_instance->soft_clippings.right[right_soft_clip_index + 1] = '\0';
			if (flag_ignore_quality_score == 0)
			{
				i++;
				sam_alignment_instance->soft_clippings.right_qual[right_soft_clip_index] = cigar_items_instance[i].def - 90;
				sam_alignment_instance->soft_clippings.right_qual[right_soft_clip_index + 1] = '\0';
			}
			right_soft_clip_index++;
		}
		else // Capital letters & special characters - cannot be soft clips
		{
			flag = 0;
			processing_left_soft_clip = 0;
			if (left_soft_clip_index > 0) //There were some left soft clips
			{
				sprintf(temp, "%d", left_soft_clip_index);
				strcat(sam_alignment_instance->cigar, temp);
				strcat(sam_alignment_instance->cigar, "S");
				num = 0;
				left_soft_clip_index = 0;
			}
			if (cigar_items_instance[i].def == 'N')
			{
				sprintf(temp, "%d", cigar_items_instance[i].len);
				strcat(sam_alignment_instance->cigar, temp);
				strcat(sam_alignment_instance->cigar, "N");
				num = 0;
			}
			else if (cigar_items_instance[i].def == 'D')
			{
				sprintf(temp, "%d", cigar_items_instance[i].len);
				strcat(sam_alignment_instance->cigar, temp);
				strcat(sam_alignment_instance->cigar, "D");
				num = 0;
			}
			num = 0;
			// Merging consecutive matches, mismatches, insertions
			while (i < cigar_items_instance_index && (isCharacterInString("BEFHJKLOPQRU", cigar_items_instance[i].def) || isCharacterInString(mismatch_characters, cigar_items_instance[i].def) || isCharacterInString(insert_characters, cigar_items_instance[i].def)))
			{
				if (isCharacterInString("BEFHJKLOPQRU", cigar_items_instance[i].def))
				{
					for (j = 0; j < cigar_items_instance[i].len; j++)
					{
						sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++] = '-';
						sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index] = '\0';
					}

					if (flag_ignore_quality_score == 0)
					{
						for (j = 0; j < cigar_items_instance[i].len; j++)
						{
							sam_alignment_instance->soft_clips_removed_qual[soft_clips_removed_qual_index] = default_quality_value[0];
							sam_alignment_instance->soft_clips_removed_qual[soft_clips_removed_qual_index + 1] = '\0';
							soft_clips_removed_qual_index++;
						}
					}
					num += cigar_items_instance[i].len;
				}
				else if (isCharacterInString(mismatch_characters, cigar_items_instance[i].def) || isCharacterInString(insert_characters, cigar_items_instance[i].def))
				{
					if (isCharacterInString(mismatch_characters, cigar_items_instance[i].def))
					{
						switch (cigar_items_instance[i].def)
						{
							case '&':
								sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++] = 'A';
								break;
							case '\'':
								sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++] = 'T';
								break;
							case '(':
								sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++] = 'G';
								break;
							case ')':
								sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++] = 'C';
								break;
							case '*':
								sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++] = 'N';
								break;
						}
					}

					if (isCharacterInString(insert_characters, cigar_items_instance[i].def))
					{
						switch (cigar_items_instance[i].def)
						{
							case '!':
								sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++] = 'A';
								break;
							case '"':
								sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++] = 'T';
								break;
							case '#':
								sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++] = 'G';
								break;
							case '$':
								sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++] = 'C';
								break;
							case '%':
								sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++] = 'N';
								break;
						}
					}

					if (flag_ignore_quality_score == 0)
					{
						i++;
						sam_alignment_instance->soft_clips_removed_qual[soft_clips_removed_qual_index] = cigar_items_instance[i].def - 90;
						sam_alignment_instance->soft_clips_removed_qual[soft_clips_removed_qual_index + 1] = '\0';
						soft_clips_removed_qual_index++;
					}
					num++;
				}
				i++;
				flag = 1;
			}
			if (flag)
			{
				sprintf(temp, "%d", num);
				strcat(sam_alignment_instance->cigar, temp);
				strcat(sam_alignment_instance->cigar, "M");
				num = 0;
			}
			if (flag && i < cigar_items_instance_index) i--;
		}
	}
	if (right_soft_clip_index > 0) //There were some right soft clips
	{
		sprintf(temp, "%d", right_soft_clip_index);
		strcat(sam_alignment_instance->cigar, temp);
		strcat(sam_alignment_instance->cigar, "S");
		num = 0;
	}

	strcpy(sam_alignment_instance->seq, sam_alignment_instance->soft_clippings.left);
	strcat(sam_alignment_instance->seq, sam_alignment_instance->soft_clips_removed_seq);
	strcat(sam_alignment_instance->seq, sam_alignment_instance->soft_clippings.right);

	strcpy(sam_alignment_instance->qual, sam_alignment_instance->soft_clippings.left_qual);
	strcat(sam_alignment_instance->qual, sam_alignment_instance->soft_clips_removed_qual);
	strcat(sam_alignment_instance->qual, sam_alignment_instance->soft_clippings.right_qual);

	splitCigar(sam_alignment_instance->cigar, &sam_alignment_instance->number_of_cigar_items, sam_alignment_instance->cigar_items);
	if (flag_ignore_sequence_information == 0) generateReadSequenceAndMDString(sam_alignment_instance, whole_genome);
	else
	{
		splitCigar(sam_alignment_instance->cigar, &cigar_items_instance_index, cigar_items_instance);
		int length_of_read = 0;
		for (i = 0; i < cigar_items_instance_index; i++)
			if (cigar_items_instance[cigar_items_instance_index].def == 'M') length_of_read += cigar_items_instance[cigar_items_instance_index].len;
		for (i = 0; i < length_of_read; i++)
		{
			sam_alignment_instance->seq[i] = 'A';
			sam_alignment_instance->qual[i] = default_quality_value[0];
		}
		sam_alignment_instance->seq[i] = '\0';
		sam_alignment_instance->qual[i] = '\0';
	}

	/*
	 //printf("\nsoft_clips_removed_qual_index:%d ", soft_clips_removed_qual_index);
	 //printf("\nflag_ignore_mismatches %d flag_ignore_soft_clippings %d flag_ignore_unmapped_sequences %d flag_ignore_quality_score %d", flag_ignore_mismatches, flag_ignore_soft_clippings, flag_ignore_unmapped_sequences, flag_ignore_quality_score);
	 *
	 printf("\nICIGAR %s ", sam_alignment_instance->icigar);
	 printf("\nCIGAR %s ", sam_alignment_instance->cigar);
	 printf("\nSEQ  %s ", sam_alignment_instance->seq);
	 printf("\nQUAL %s ", sam_alignment_instance->qual);
	 printf("\nNH:%d \nsamformat_flag: %d", NH_value, samformatflag);
	 printf("\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
	 printf("\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
	 */

}

void writeAlignmentToFile(struct Sam_Alignment *sam_alignment, short int flag_ignore_sequence_information, int number_of_repititions_of_the_same_reads, FILE *fhw)
{
	int i;
	char line_to_be_written_to_file[MAX_GENERAL_LEN];
	char temp[100];
	for (i = 0; i < number_of_repititions_of_the_same_reads; i++)
	{
		line_to_be_written_to_file[0] = '\0';
		strcat(line_to_be_written_to_file, sam_alignment->read_name);
		sprintf(temp, "%d", i + 1);
		strcat(line_to_be_written_to_file, "_");
		strcat(line_to_be_written_to_file, temp);
		strcat(line_to_be_written_to_file, "\t");

		sprintf(temp, "%d", sam_alignment->samflag);
		strcat(line_to_be_written_to_file, temp);

		strcat(line_to_be_written_to_file, "\t");
		strcat(line_to_be_written_to_file, sam_alignment->reference_name);

		strcat(line_to_be_written_to_file, "\t");
		sprintf(temp, "%d", sam_alignment->start_position);
		strcat(line_to_be_written_to_file, temp);

		strcat(line_to_be_written_to_file, "\t");
		strcat(line_to_be_written_to_file, "255");

		strcat(line_to_be_written_to_file, "\t");
		strcat(line_to_be_written_to_file, sam_alignment->cigar);

		strcat(line_to_be_written_to_file, "\t");
		strcat(line_to_be_written_to_file, "*");

		strcat(line_to_be_written_to_file, "\t");
		strcat(line_to_be_written_to_file, "*");

		strcat(line_to_be_written_to_file, "\t");
		strcat(line_to_be_written_to_file, "*");

		strcat(line_to_be_written_to_file, "\t");
		strcat(line_to_be_written_to_file, sam_alignment->seq);

		strcat(line_to_be_written_to_file, "\t");
		strcat(line_to_be_written_to_file, sam_alignment->qual);

		strcat(line_to_be_written_to_file, "\t");
		//Tags
		strcat(line_to_be_written_to_file, "NH:i:");
		strcat(line_to_be_written_to_file, sam_alignment->tags[0].val);
		strcat(line_to_be_written_to_file, "\t");

		if (strcmp(sam_alignment->tags[1].val, ".") != 0)
		{
			strcat(line_to_be_written_to_file, "XS:");
			strcat(line_to_be_written_to_file, sam_alignment->tags[1].val);
			strcat(line_to_be_written_to_file, "\t");
		}

		if (flag_ignore_sequence_information != 0)
		{
			strcat(line_to_be_written_to_file, "MD:i");
			strcat(line_to_be_written_to_file, sam_alignment->tags[2].val);
			strcat(line_to_be_written_to_file, "\t");
		}

		strcat(line_to_be_written_to_file, "\n");
		fprintf(fhw, "%s", line_to_be_written_to_file);
	}
}

void convertToAlignment(struct Whole_Genome_Sequence *whole_genome, struct Sam_Alignment *sam_alignment_instance, int sam_alignment_pool_index, struct Whole_Genome_Sequence *whole_genome, char **split_on_newline, struct Sam_Alignment *sam_alignment, int cluster_index, struct Abridge_Index *abridge_index, int number_of_entries_in_cluster, char **split_on_tab, char **split_on_dash, char **split_on_comma, char *default_quality_value, short int flag_ignore_mismatches, short int flag_ignore_soft_clippings, short int flag_ignore_unmapped_sequences, short int flag_ignore_quality_score, short int flag_ignore_sequence_information, unsigned long long int *read_number, FILE *fhw)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	long long int start_position;

	int number_of_columns;
	int number_of_distinct_cigars_in_a_line;
	int number_of_repititions_of_the_same_reads;
	int samformatflag;
	int i;
	int j;

	char *temp; //Useless
	char *distinct_icigars_in_a_line;
	char *icigar;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	distinct_icigars_in_a_line = (char*) malloc(sizeof(char) * MAX_GENERAL_LEN);

	/********************************************************************/
	for (i = 0; i < number_of_entries_in_cluster; i++)
	{
		//printf("\nProcessing this line %s", split_on_newline[i]);
		number_of_columns = splitByDelimiter(split_on_newline[i], '\t', split_on_tab);
		//printf("\nNumber of columns %d", number_of_columns);

		if (number_of_columns == 1)
		{
			strcpy(distinct_icigars_in_a_line, split_on_tab[0]);
			start_position++;
		}
		else if (number_of_columns == 2)
		{
			if (i == 0) start_position = abridge_index->start[cluster_index];
			else start_position += strtol(split_on_tab[0], &temp, 10);
			strcpy(distinct_icigars_in_a_line, split_on_tab[1]);
		}
		else
		{
			// Should never enter here
			printf("\nTrouble");
		}
		number_of_distinct_cigars_in_a_line = splitByDelimiter(distinct_icigars_in_a_line, ',', split_on_comma);
		sam_alignment_pool_index = 0;
		for (j = 0; j < number_of_distinct_cigars_in_a_line; j++)
		{
			splitByDelimiter(split_on_comma[j], '-', split_on_dash);
			number_of_repititions_of_the_same_reads = strtol(split_on_dash[1], &temp, 10);
			sam_alignment_instance->start_position = start_position;
			if (split_on_dash[0][1] == '-' && isalpha(split_on_dash[0][0]) != 0)
			{
				// Use the same cigar
				sprintf(temp, "%d", *read_number);
				(*read_number)++;
				strcpy(sam_alignment_instance->read_name, temp);
			}
			else
			{
				//printf("\nSplit_on_dash %s %d sam_alignment_pool_index %d ", split_on_dash[0], strlen(split_on_dash[0]), sam_alignment_pool_index);
				strcpy(sam_alignment_instance->icigar, split_on_dash[0]);
				//printf("\nsam_alignment_pool[i]->icigar %s %d Original_icigar %s %d", sam_alignment_pool[sam_alignment_pool_index]->icigar, strlen(sam_alignment_pool[sam_alignment_pool_index]->icigar), split_on_dash[0], strlen(split_on_dash[0]));
				convertIcigarToCigarandMD(whole_genome, sam_alignment_instance, abridge_index->chromosome[cluster_index], flag_ignore_mismatches, flag_ignore_soft_clippings, flag_ignore_unmapped_sequences, flag_ignore_quality_score, flag_ignore_sequence_information, default_quality_value);
				sprintf(temp, "%d", *read_number);
				(*read_number)++;
				strcpy(sam_alignment_instance->read_name, temp);
			}
			writeAlignmentToFile(sam_alignment_instance, flag_ignore_sequence_information, number_of_repititions_of_the_same_reads, fhw);
			sam_alignment_pool_index++;
		}
	}
}

void decompressFile(char *abridge_index_filename, char *genome_filename, char *output_sam_filename, char *pass2_filename, char *genome_prefix, char *default_quality_value, short int flag_ignore_sequence_information)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	FILE *fhr;
	FILE *fhw;

	struct Abridge_Index *abridge_index;
	struct Sam_Alignment *sam_alignment;

	int i;
	int sam_alignment_pool_index;
	int fread_ret_val;
	int fseek_ret_val;

	short int flag_ignore_mismatches;
	short int flag_ignore_soft_clippings;
	short int flag_ignore_unmapped_sequences;
	short int flag_ignore_quality_score;

	long long int max_cluster_size;
	unsigned long long int line_num = 0;
	unsigned long long int read_number = 1;
	int number_of_entries_in_cluster;
	int number_of_elements_after_split_on_delimiter;

	char **split_on_newline;
	char **split_on_tab;
	char **split_on_dash;
	char **split_on_comma;
	char **split_on_delimiter;
	char *buffer;
	char **sequence_portions_from_reference;
	char *fasta_file_with_expressed_portions;
	char *cigar;
	char *md;
	char *output_prefix_without_path;

	struct Sam_Alignment **sam_alignment_pool;
	struct Sam_Alignment *sam_alignment_instance;
	struct Whole_Genome_Sequence *whole_genome;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr = fopen(pass2_filename, "rb");
	if (fhr == NULL)
	{
		printf("Error! File not found");
		exit(1);
	}
	fhw = fopen(output_sam_filename, "w");
	if (fhw == NULL)
	{
		printf("Error! File cannot be opened for writing");
		exit(1);
	}

	split_on_newline = (char**) malloc(sizeof(char*) * ROWS * 10);
	for (i = 0; i < ROWS * 10; i++)
		split_on_newline[i] = (char*) malloc(sizeof(char) * COLS * 10);

	split_on_tab = (char**) malloc(sizeof(char*) * ROWS);
	for (i = 0; i < ROWS; i++)
		split_on_tab[i] = (char*) malloc(sizeof(char) * COLS);

	split_on_dash = (char**) malloc(sizeof(char*) * ROWS);
	for (i = 0; i < ROWS; i++)
		split_on_dash[i] = (char*) malloc(sizeof(char) * COLS);

	split_on_comma = (char**) malloc(sizeof(char*) * ROWS * 10);
	for (i = 0; i < ROWS * 10; i++)
		split_on_comma[i] = (char*) malloc(sizeof(char) * COLS * 10);

	split_on_delimiter = (char**) malloc(sizeof(char*) * ROWS);
	for (i = 0; i < ROWS; i++)
		split_on_delimiter[i] = (char*) malloc(sizeof(char) * COLS);

	output_prefix_without_path = (char*) malloc(sizeof(char) * MAX_SEQ_LEN);
	sequence_portions_from_reference = (char**) malloc(sizeof(char*) * MAX_POOL_SIZE);
	fasta_file_with_expressed_portions = (char*) malloc(sizeof(char) * FILENAME_LENGTH);

	buffer = (char*) malloc(sizeof(char) * MAX_BUFFER_SIZE_FOR_READING_PASS2_FILE);
	abridge_index = allocateMemoryAbridge_Index();
	sam_alignment = allocateMemorySam_Alignment();
	whole_genome = (struct Whole_Genome_Sequence*) malloc(sizeof(struct Whole_Genome_Sequence));
	sam_alignment_instance = allocateMemorySam_Alignment();

	/********************************************************************/

	readAbridgeIndex(abridge_index, abridge_index_filename, split_on_newline, &flag_ignore_mismatches, &flag_ignore_soft_clippings, &flag_ignore_unmapped_sequences, &flag_ignore_quality_score);
	readInTheEntireGenome(genome_filename, whole_genome);

	printf("\nProcessing each line from index file");
	fflush(stdout);
	for (i = 0; i < abridge_index->number_of_items; i++)
	{
		if (i % 1000 == 0)
		{
			printf("\n%d record processed", i);
			fflush(stdout);
		}
		//printf("\nFile pointer at %ld", ftell(fhr));
		//printf("\nEntry in index file %s %d %d", abridge_index->chromosome[i], abridge_index->start_byte[i], abridge_index->end_byte[i]);
		fseek_ret_val = fseek(fhr, abridge_index->start_byte[i], SEEK_SET);
		//printf("\nFile pointer moved to %ld", ftell(fhr));
		//continue;
		buffer[0] = '\0';
		fread_ret_val = fread(buffer, 1, abridge_index->end_byte[i] - abridge_index->start_byte[i], fhr);
		//printf("\n fread_ret_val %d fseek_ret_val %d diff %lld", fread_ret_val, fseek_ret_val, MAX_BUFFER_SIZE_FOR_READING_PASS2_FILE - fread_ret_val);
		//fflush(stdout);
		buffer[fread_ret_val] = '\0';
		/*
		 printf("\nfread_ret_val %d", fread_ret_val);
		 printf("\nNew Record:%d\n", i);
		 printf("%s", buffer);
		 printf("\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
		 printf("\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
		 fflush(stdout);
		 */
		number_of_entries_in_cluster = splitByDelimiter(buffer, '\n', split_on_newline);
		number_of_entries_in_cluster--; //Last line is always empty
		convertToAlignment(whole_genome, sam_alignment_instance, sam_alignment_pool_index, whole_genome, split_on_newline, sam_alignment, i, abridge_index, number_of_entries_in_cluster, split_on_tab, split_on_dash, split_on_comma, default_quality_value, flag_ignore_mismatches, flag_ignore_soft_clippings, flag_ignore_unmapped_sequences, flag_ignore_quality_score, flag_ignore_sequence_information, &read_number, fhw);
		//if (i == 10) break;
	}
}

int main(int argc, char *argv[])
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char abridge_index_filename[FILENAME_LENGTH];
	char genome_filename[FILENAME_LENGTH];
	char pass2_filename[FILENAME_LENGTH];
	char output_sam_filename[FILENAME_LENGTH];
	char genome_prefix[FILENAME_LENGTH];
	char default_quality_value[10];
	char *temp; //Required for strtoi

	short int flag_ignore_sequence_information;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	strcpy(abridge_index_filename, argv[1]);
	strcpy(genome_filename, argv[2]);
	strcpy(output_sam_filename, argv[3]);
	strcpy(pass2_filename, argv[4]);
	strcpy(genome_prefix, argv[5]);
	strcpy(default_quality_value, argv[6]);
	flag_ignore_sequence_information = strtol(argv[7], &temp, 10);
	/********************************************************************/

	decompressFile(abridge_index_filename, genome_filename, output_sam_filename, pass2_filename, genome_prefix, default_quality_value, flag_ignore_sequence_information);
	return 0;
}

