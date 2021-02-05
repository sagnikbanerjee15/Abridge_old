#ifndef _ABRIDGE_SRC_FUNCTIONS_DEFINITIONS_H_
#define _ABRIDGE_SRC_FUNCTIONS_DEFINITIONS_H_

struct Sam_Alignment* allocateMemorySam_Alignment()
{
	/********************************************************************
	 * Variable declarations
	 ********************************************************************/
	struct Sam_Alignment *s;
	int i;
	/********************************************************************/
	s = (struct Sam_Alignment*) malloc(sizeof(struct Sam_Alignment));

	s->read_name = (char*) malloc(sizeof(char) * MAX_SEQ_LEN);
	s->reference_name = (char*) malloc(sizeof(char) * MAX_SEQ_LEN);
	s->cigar = (char*) malloc(sizeof(char) * MAX_SEQ_LEN);
	s->reference_name_next_mate = (char*) malloc(sizeof(char) * MAX_SEQ_LEN);
	s->seq = (char*) malloc(sizeof(char) * MAX_SEQ_LEN);
	s->qual = (char*) malloc(sizeof(char) * MAX_SEQ_LEN);
	s->soft_clippings.left = (char*) malloc(sizeof(char) * MAX_SEQ_LEN);
	s->soft_clippings.right = (char*) malloc(sizeof(char) * MAX_SEQ_LEN);
	s->cigar_items = (struct Cigar_Items*) malloc(sizeof(struct Cigar_Items) * MAX_CIGAR_ITEMS);
	s->number_of_cigar_items = 0;
	s->number_of_tag_items = 0;

	for (i = 0; i < 50; i++)
	{
		s->tags[i].name = (char*) malloc(sizeof(char) * MAX_TAG_LENGTH);
		s->tags[i].val = (char*) malloc(sizeof(char) * MAX_TAG_LENGTH);
		s->tags[i].type = (char*) malloc(sizeof(char) * MAX_TAG_LENGTH);
		s->tags[i].name[0] = '\0';
		s->tags[i].val[0] = '\0';
		s->tags[i].type[0] = '\0';
	}

	s->cigar_extended = (char*) malloc(sizeof(char) * (MAX_SEQ_LEN * 2));
	s->md_extended = (char*) malloc(sizeof(char) * (MAX_SEQ_LEN * 2));
	s->icigar = (char*) malloc(sizeof(char) * (MAX_SEQ_LEN * 2));
	s->splices = (char**) malloc(sizeof(char*) * 100);
	for (i = 0; i < 100; i++)
		s->splices[i] = (char*) malloc(sizeof(char) * 50);
	s->soft_clips_removed_seq = (char*) malloc(sizeof(char) * MAX_SEQ_LEN);
	s->temp = (char*) malloc(sizeof(char) * MAX_SEQ_LEN);

	s->temp[0] = '\0';
	s->cigar_extended[0] = '\0';
	s->icigar[0] = '\0';
	s->md_extended[0] = '\0';
	s->soft_clips_removed_seq[0] = '\0';
	s->soft_clippings.left[0] = '\0';
	s->soft_clippings.right[0] = '\0';

	return s;
}

struct Pass3_Compression_Symbol_icigar_Mapping* allocateMemoryPass3_Compression_Symbol_icigar_Mapping()
{
	struct Pass3_Compression_Symbol_icigar_Mapping *s;
	s = (struct Pass3_Compression_Symbol_icigar_Mapping*) malloc(sizeof(struct Pass3_Compression_Symbol_icigar_Mapping));
	s->icigar = (char*) malloc(sizeof(char) * MAX_ICIGAR_LENGTH);
	s->symbolic_icigar = (char*) malloc(sizeof(char) * MAX_SYMBOLIC_ICIGAR_LENGTH);
	return s;
}

struct Compressed_DS* allocateMemoryCompressed_DS()
{
	struct Compressed_DS *s;
	s = (struct Compressed_DS*) malloc(sizeof(struct Compressed_DS));
	s->icigar = (char*) malloc(sizeof(char) * (MAX_SEQ_LEN * 2));
	s->num_reads = 0;
	s->position = 0;
	return s;
}

struct Pass1_Compressed_DS* allocateMemoryPass1_Compressed_DS()
{
	struct Pass1_Compressed_DS *s;
	s = (struct Pass1_Compressed_DS*) malloc(sizeof(struct Pass1_Compressed_DS));
	s->col1 = (char*) malloc(sizeof(char) * (MAX_SEQ_LEN * 2));
	s->col2 = (char*) malloc(sizeof(char) * MAX_ICIGAR_LENGTH_PASS1_COL2);
	s->col3 = (char*) malloc(sizeof(char) * (MAX_SEQ_LEN * 2));
	return s;
}

struct Pass2_Compressed_DS* allocateMemoryPass2_Compressed_DS()
{
	struct Pass2_Compressed_DS *s;
	s = (struct Pass2_Compressed_DS*) malloc(sizeof(struct Pass2_Compressed_DS));
	s->col1 = (char*) malloc(sizeof(char) * (MAX_SEQ_LEN * 2));
	s->col2 = (char*) malloc(sizeof(char) * MAX_ICIGAR_LENGTH_PASS1_COL2);
	return s;
}

struct Cigar_Frequency* allocateMemoryCigar_Frequency()
{
	struct Cigar_Frequency *s;
	s = (struct Cigar_Frequency*) malloc(sizeof(struct Cigar_Frequency));
	s->cigar = (char*) malloc(sizeof(char) * MAX_CIGAR_LENGTH);
	s->freq = 0;
	return s;
}

struct Reference_Sequence_Info* allocateMemoryReference_Sequence_Info()
{
	struct Reference_Sequence_Info *s;
	s = (struct Reference_Sequence_Info*) malloc(sizeof(struct Reference_Sequence_Info));
	s->line = (char*) malloc(sizeof(char) * MAX_SEQ_LEN);
	return s;
}

struct Abridge_Index* allocateMemoryAbridge_Index()
{
	int i;
	struct Abridge_Index *s;
	s = (struct Abridge_Index*) malloc(sizeof(struct Abridge_Index));
	s->chromosome = (char**) malloc(sizeof(char*) * MAX_POOL_SIZE);
	for (i = 0; i < MAX_POOL_SIZE; i++)
		s->chromosome[i] = (char*) malloc(sizeof(char) * MAX_REFERENCE_SEQ_LENGTH);
	s->start = (long long int*) malloc(sizeof(long long int) * MAX_POOL_SIZE);
	s->end = (long long int*) malloc(sizeof(long long int) * MAX_POOL_SIZE);
	s->start_byte = (long long int*) malloc(sizeof(long long int) * MAX_POOL_SIZE);
	s->end_byte = (long long int*) malloc(sizeof(long long int) * MAX_POOL_SIZE);
	s->number_of_items = 0;
	return s;
}

int isCommaInLine(char *line)
{
	int i;

	for (i = 0; line[i] != '\0'; i++)
		if (line[i] == ',') return 1;

	return 0;
}

int splitByDelimiter(char *line, char delimiter, char **new_string)
{
	/*********************************************************************
	 * Splits a string one a requested character and
	 * returns the number of splits
	 *********************************************************************/

	/********************************************************************
	 * Variable declarations
	 ********************************************************************/
	int i = 0, j = 0, ctr = 0;
	/********************************************************************/

	for (i = 0; line[i] != '\0'; i++)
	{
		// if space or NULL found, assign NULL into new_string[ctr]
		if (line[i] == delimiter)
		{
			new_string[ctr][j] = '\0';
			ctr++; //for next word
			j = 0; //for next word, init index to 0
		}
		else if (line[i] == '\n') continue;
		else
		{
			new_string[ctr][j] = line[i];
			j++;
		}
	}
	new_string[ctr][j] = '\0';
	return ctr + 1;
}

int locateSamTags(char *tag)
{
	/********************************************************************
	 * Variable declarations
	 ********************************************************************/
	int i = 0;
	/********************************************************************/

	for (i = 0; i < sizeof(sam_tags) / sizeof(sam_tags[0]); i++)
	{
		if (strcmp(sam_tags[i], tag) == 0) return i;
	}
	return -1;

}

void freeSamAlignmentInstance(struct Sam_Alignment *s)
{
	int i;
	for (i = 0; i < 100; i++)
	{
		free(s->tags[i].name);
		free(s->tags[i].type);
		free(s->tags[i].val);
	}
	free(s->read_name);
	free(s->reference_name);
	free(s->reference_name_next_mate);
	free(s->cigar);
	free(s->seq);
	free(s->qual);
	free(s->soft_clippings.left);
	free(s->soft_clippings.right);
	free(s->cigar_items);
	free(s->cigar_extended);
	free(s->md_extended);
	free(s->icigar);
	for (i = 0; i < ROWS; i++)
		free(s->splices[i]);
	free(s->temp);
}

void populateSamAlignmentInstance(struct Sam_Alignment *dest, char **src, int number_of_fields, char **split_tags)
{
	/********************************************************************
	 * Variable declarations
	 ********************************************************************/
	char *temp;
	int i;
	/********************************************************************/

// Assign the first 11 mandatory sam fields
	/*
	 dest->read_name = src[0];
	 dest->samflag = strtol(src[1], &temp, 10);
	 dest->reference_name = src[2];
	 dest->start_position = strtol(src[3], &temp, 10);
	 dest->mapping_quality_score = strtol(src[4], &temp, 10);
	 dest->cigar = src[5];
	 dest->reference_name_next_mate = src[6];
	 dest->start_position_next = strtol(src[7], &temp, 10);
	 dest->template_length = strtol(src[8], &temp, 10);
	 dest->seq = src[9];
	 dest->qual = src[10];
	 */

	strcpy(dest->read_name, src[0]);
	dest->samflag = strtol(src[1], &temp, 10);
	strcpy(dest->reference_name, src[2]);
	dest->start_position = strtol(src[3], &temp, 10);
	dest->mapping_quality_score = strtol(src[4], &temp, 10);
	strcpy(dest->cigar, src[5]);
	strcpy(dest->reference_name_next_mate, src[6]);
	dest->start_position_next = strtol(src[7], &temp, 10);
	dest->template_length = strtol(src[8], &temp, 10);
	strcpy(dest->seq, src[9]);
	strcpy(dest->qual, src[10]);

// Assign SAM tags
	for (i = 11; i < number_of_fields; i++)
	{
		//printf("%s\n",src[i]);
		//fflush ( stdout );
		splitByDelimiter(src[i], ':', split_tags);
		//sam_tag_index = locateSamTags(split_tags[0]);
		//dest->tags[i - 11].name = sam_tags[sam_tag_index];
		if (i == number_of_fields - 1) split_tags[2][strlen(split_tags[2])] = '\0';
		strcpy(dest->tags[i - 11].name, split_tags[0]);
		strcpy(dest->tags[i - 11].type, split_tags[1]);
		strcpy(dest->tags[i - 11].val, split_tags[2]);
		//printf("\n Tags %s Parts of the tag %s %s %s ", src[i], split_tags[0], split_tags[1], split_tags[2]);
	}
	dest->read_seq_len = strlen(dest->seq);
	dest->number_of_tag_items = number_of_fields - 11;
}

void printSamAlignmentInstance(struct Sam_Alignment *s)
{
	int i;
	printf("\n");
	printf("^^^^^^^^^^^^^^^^^SAM ALIGNMENT^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
	printf("\n");
	printf("Read name: %s", s->read_name);
	printf("\n");
	printf("SamFlag: %d", s->samflag);
	printf("\n");
	printf("Reference name: %s", s->reference_name);
	printf("\n");
	printf("Start position: %lld", s->start_position);
	printf("\n");
	printf("Mapping quality score: %d", s->mapping_quality_score);
	printf("\n");
	printf("Cigar: %s", s->cigar);
	printf("\n");
	printf("Reference name next mate: %s", s->reference_name_next_mate);
	printf("\n");
	printf("Start position next: %lld", s->start_position_next);
	printf("\n");
	printf("Template length: %d", s->template_length);
	printf("\n");
	printf("Sequence: %s %d ", s->seq, s->read_seq_len);
	printf("\n");
	printf("Quality scores: %s", s->qual);
	printf("\n");
	for (i = 0; i < s->number_of_tag_items; i++)
	{
		printf("Tag item %s Tag value %s", s->tags[i].name, s->tags[i].val);
		printf("\n");
	}
	printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
	fflush (stdout);
}

int isSequenceSoftClipped(char *cigar)
{
	int i;
	for (i = 0; cigar[i] != '\0'; i++)
		if (cigar[i] == 'S') return 1;
	return 0;
}

void extractSubString(char *str, char *substr, int start_index, int end_index)
{
	/********************************************************************
	 * Variable declarations
	 ********************************************************************/
	int length = strlen(str);
	int i;
	int j = 0;

	/********************************************************************/

	if (length < (end_index - start_index + 1))
	{
		return;
	}
	for (i = start_index; i <= end_index; i++)
		substr[j++] = str[i];
	substr[j] = '\0';

}

void splitCigar(char *cigar, int *num_of_types, struct Cigar_Items *cigar_items_instance)
{
	/********************************************************************
	 * Variable declarations
	 ********************************************************************/

	int i, j = 0;
	int current_length = 0;
	/********************************************************************/

	for (i = 0; cigar[i] != '\0'; i++)
	{
		if (isdigit(cigar[i]) != 0) //cigar[i] is a digit
		{
			current_length = current_length * 10 + cigar[i] - 48;
		}
		else //cigar[i] is a character
		{
			cigar_items_instance[j].def = cigar[i];
			cigar_items_instance[j].len = current_length;
			current_length = 0;
			j += 1;
		}
	}
	*num_of_types = j;
}

int isCharacterInString(char *str, char character)
{
	int i;
	for (i = 0; str[i] != '\0'; i++)
		if (character == str[i]) return 1;
	return 0;
}

int isNumber(char *str)
{
	int i;
	for (i = 0; str[i] != '\0'; i++)
		if (str[i] < 48 || str[i] > 57) return 0;
	return 1;
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

void designIntegratedCIGAR(char *md, char *seq, int *seq_length, struct Cigar_Items *cigar_items_instance, int number_of_cigar_items, char *cigar, char *cigar_extended, char *md_extended, char *icigar, char **splices)
{
	int i;
	int num;
	int j;
	int md_extended_index;
	int splice_number_index;
	int md_extended_length;
	int splices_index;
	char temp_convert_int_to_string[100];
	int intron_splice_index;
	int print_outputs = 0;

	md_extended_index = 0;
	i = 0;
	num = 0;
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
		if (seq[i] >= 33 && seq[i] <= 37) insertCharacterInString(md_extended, &md_extended_length, seq[i], i, 1);
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
				case 'N':
					md_extended[i] = mismatch_characters[4];
					break;
			}
		}
	}

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
}

int isAlignmentPerfect(char *cigar, struct Sam_Tags *tags, int MD_tag_index, int NM_tag_index, int nM_tag_index)
{
	int i;
	if (NM_tag_index != -1 && nM_tag_index != -1)
	{
		if (tags[NM_tag_index].val[0] == '0' && tags[nM_tag_index].val[0] == '0') return 1;
		else return 0;
	}

	/*
	 * Checking for insertions
	 */
	for (i = 0; cigar[i] != '\0'; i++)
		if (cigar[i] == 'I') return 0;

	/*
	 * Checking for deletions
	 */
	for (i = 0; tags[MD_tag_index].val[i] != '\0'; i++)
		if (tags[MD_tag_index].val[i] == '^' || (tags[MD_tag_index].val[i] >= 65 && tags[MD_tag_index].val[i] <= 90) || (tags[MD_tag_index].val[i] >= 97 && tags[MD_tag_index].val[i] <= 122)) return 0;

	return 1;
}

void generateIntegratedCigar(struct Sam_Alignment *curr_alignment, short int flag_ignore_soft_clippings, short int flag_ignore_mismatches, short int flag_ignore_unmapped_sequences)
{
	/*
	 * Creates the integrated cigar
	 */
	/********************************************************************
	 * Variable declarations
	 ********************************************************************/
	int left_soft_clip_point = 0;
	int right_soft_clip_point = 0;
	int i;
	int number_of_mismatches_in_each_pair = -1;
	int MD_tag_index;
	int nM_tag_index;
	int XS_tag_index;
	int NH_tag_index;
	int NM_tag_index;
	int print_outputs = 0;
	int perfect_alignment_indicator = 0;
	int spliced_alignment_indicator = 0;
	char str[1000];
	char temp_str[5];
	char M_replacement_character;
	/********************************************************************/
	splitCigar(curr_alignment->cigar, &curr_alignment->number_of_cigar_items, curr_alignment->cigar_items);

	for (i = 0; i < curr_alignment->number_of_cigar_items; i++)
	{
		if (curr_alignment->cigar_items[i].def == 'N')
		{
			spliced_alignment_indicator = 1;
			break;
		}
	}
	/*
	 * Process each alignment to extract soft clipped portion of reads
	 */
	if (isSequenceSoftClipped(curr_alignment->cigar) == 1)
	{
		if (curr_alignment->cigar_items[0].def == 'S')
		{
			left_soft_clip_point = curr_alignment->cigar_items[0].len;
			extractSubString(curr_alignment->seq, curr_alignment->soft_clippings.left, 0, left_soft_clip_point - 1);
			if (print_outputs > 1) printf("\nLSC: %d %s", curr_alignment->cigar_items[0].len, curr_alignment->soft_clippings.left);
			curr_alignment->soft_clips_removed_seq_len = curr_alignment->read_seq_len - curr_alignment->cigar_items[0].len;
			for (i = 0; i <= left_soft_clip_point - 1; i++)
				//curr_alignment->soft_clippings.left[i] =(curr_alignment->soft_clippings.left[i] >= 65 && curr_alignment->soft_clippings.left[i] <= 90) ? curr_alignment->soft_clippings.left[i] + 32 : curr_alignment->soft_clippings.left[i];
				curr_alignment->soft_clippings.left[i] =
						(curr_alignment->soft_clippings.left[i] >= 65 && curr_alignment->soft_clippings.left[i] <= 90) ? curr_alignment->soft_clippings.left[i] + 32 : curr_alignment->soft_clippings.left[i];
		}
		if (curr_alignment->cigar_items[curr_alignment->number_of_cigar_items - 1].def == 'S')
		{
			right_soft_clip_point = curr_alignment->cigar_items[curr_alignment->number_of_cigar_items - 1].len;
			extractSubString(curr_alignment->seq, curr_alignment->soft_clippings.right, curr_alignment->read_seq_len - right_soft_clip_point, curr_alignment->read_seq_len);
			if (print_outputs > 1) printf("\nRSC: %d %s", curr_alignment->cigar_items[curr_alignment->number_of_cigar_items - 1].len, curr_alignment->soft_clippings.right);
			curr_alignment->soft_clips_removed_seq_len = curr_alignment->read_seq_len - curr_alignment->cigar_items[curr_alignment->number_of_cigar_items - 1].len;

			for (i = 0; i <= right_soft_clip_point; i++)
				//curr_alignment->soft_clippings.right[i] =(curr_alignment->soft_clippings.right[i] >= 65 && curr_alignment->soft_clippings.right[i] <= 90) ? curr_alignment->soft_clippings.right[i] + 32 : curr_alignment->soft_clippings.right[i];
				curr_alignment->soft_clippings.right[i] =
						(curr_alignment->soft_clippings.right[i] >= 65 && curr_alignment->soft_clippings.right[i] <= 90) ? curr_alignment->soft_clippings.right[i] + 32 : curr_alignment->soft_clippings.right[i];
		}

		/*
		 * Remove the soft-clips and construct new sequence
		 */
		if (left_soft_clip_point != 0 || right_soft_clip_point != 0)
		{
			for (i = left_soft_clip_point; i <= curr_alignment->read_seq_len - right_soft_clip_point; i++)
				curr_alignment->soft_clips_removed_seq[i] = curr_alignment->seq[i];
			if (print_outputs > 1) printf("\nSEQ: %s newSEQ: %s", curr_alignment->seq, curr_alignment->soft_clips_removed_seq);
		}
		else
		{
			curr_alignment->soft_clips_removed_seq = curr_alignment->seq;
			if (print_outputs > 1) printf("\nSEQ: %s newSEQ: %s", curr_alignment->seq, curr_alignment->soft_clips_removed_seq);
		}
	}
	else
	{
		curr_alignment->soft_clips_removed_seq = curr_alignment->seq;
		if (print_outputs > 1) printf("\nSEQ: %s newSEQ: %s", curr_alignment->seq, curr_alignment->soft_clips_removed_seq);
	}

	/*
	 * Find the different tags
	 */
	XS_tag_index = -1;
	NH_tag_index = -1;
	NM_tag_index = -1;
	nM_tag_index = -1;
	MD_tag_index = -1;
	for (i = 0; i < curr_alignment->number_of_tag_items; i++)
	{
		if (strcmp(curr_alignment->tags[i].name, "MD") == 0) MD_tag_index = i;
		if (strcmp(curr_alignment->tags[i].name, "XS") == 0) XS_tag_index = i;
		if (strcmp(curr_alignment->tags[i].name, "NH") == 0) NH_tag_index = i;
		if (strcmp(curr_alignment->tags[i].name, "NM") == 0) NM_tag_index = i;
		if (strcmp(curr_alignment->tags[i].name, "nM") == 0) nM_tag_index = i;
	}

	perfect_alignment_indicator = isAlignmentPerfect(curr_alignment->cigar, curr_alignment->tags, MD_tag_index, NM_tag_index, nM_tag_index);

	if (perfect_alignment_indicator == 0) designIntegratedCIGAR(curr_alignment->tags[MD_tag_index].val, curr_alignment->soft_clips_removed_seq, &curr_alignment->soft_clips_removed_seq_len, curr_alignment->cigar_items, curr_alignment->number_of_cigar_items, curr_alignment->cigar, curr_alignment->cigar_extended, curr_alignment->md_extended, curr_alignment->icigar, curr_alignment->splices);
	else
	{
		if (curr_alignment->cigar_items[0].def == 'S' && curr_alignment->cigar_items[curr_alignment->number_of_cigar_items - 1].def == 'S')
		{
			strcpy(curr_alignment->icigar, "");
			for (i = 1; i < curr_alignment->number_of_cigar_items - 1; i++)
			{
				sprintf(str, "%d", curr_alignment->cigar_items[i].len);
				strcat(curr_alignment->icigar, str);
				temp_str[0] = curr_alignment->cigar_items[i].def;
				temp_str[1] = '\0';
				strcat(curr_alignment->icigar, temp_str);
			}
		}
		else if (curr_alignment->cigar_items[0].def != 'S' && curr_alignment->cigar_items[curr_alignment->number_of_cigar_items - 1].def == 'S')
		{
			strcpy(curr_alignment->icigar, "");
			for (i = 0; i < curr_alignment->number_of_cigar_items - 1; i++)
			{
				sprintf(str, "%d", curr_alignment->cigar_items[i].len);
				strcat(curr_alignment->icigar, str);
				temp_str[0] = curr_alignment->cigar_items[i].def;
				temp_str[1] = '\0';
				strcat(curr_alignment->icigar, temp_str);
			}
		}
		else if (curr_alignment->cigar_items[0].def == 'S' && curr_alignment->cigar_items[curr_alignment->number_of_cigar_items - 1].def != 'S')
		{
			strcpy(curr_alignment->icigar, "");
			for (i = 1; i < curr_alignment->number_of_cigar_items; i++)
			{
				sprintf(str, "%d", curr_alignment->cigar_items[i].len);
				strcat(curr_alignment->icigar, str);
				temp_str[0] = curr_alignment->cigar_items[i].def;
				temp_str[1] = '\0';
				strcat(curr_alignment->icigar, temp_str);
			}
		}
		else strcpy(curr_alignment->icigar, curr_alignment->cigar);
	}
	if (print_outputs > 1) printf("\nPrevious iCIGAR: %s", curr_alignment->icigar);

	if (left_soft_clip_point != 0 && right_soft_clip_point == 0)
	{
		strcpy(curr_alignment->temp, curr_alignment->soft_clippings.left);
		strcat(curr_alignment->temp, curr_alignment->icigar);
	}
	else if (left_soft_clip_point == 0 && right_soft_clip_point != 0)
	{
		strcpy(curr_alignment->temp, curr_alignment->icigar);
		strcat(curr_alignment->temp, curr_alignment->soft_clippings.right);
	}
	else if (left_soft_clip_point != 0 && right_soft_clip_point != 0)
	{
		strcpy(curr_alignment->temp, curr_alignment->soft_clippings.left);
		strcat(curr_alignment->temp, curr_alignment->icigar);
		strcat(curr_alignment->temp, curr_alignment->soft_clippings.right);
	}
	else
	{
		strcpy(curr_alignment->temp, curr_alignment->icigar);
	}

	strcpy(curr_alignment->icigar, curr_alignment->temp);

	/*
	 * Add NH tag
	 */
	/*if (XS_tag_index != -1) strcat(curr_alignment->icigar, curr_alignment->tags[XS_tag_index].val);*/
	if (NH_tag_index != -1) strcat(curr_alignment->icigar, curr_alignment->tags[NH_tag_index].val);

	/*
	 * Change the iCIGAR representation to reflect the samformatflag and XS tag
	 */
	if (spliced_alignment_indicator == 0)
	{
		switch (curr_alignment->samflag)
		{
			case 0:
				M_replacement_character = 'B';
				break;
			case 16:
				M_replacement_character = 'E';
				break;
			case 256:
				M_replacement_character = 'F';
				break;
			case 272:
				M_replacement_character = 'H';
				break;
				break;
		}
		for (i = 0; i < strlen(curr_alignment->icigar); i++)
			if (curr_alignment->icigar[i] == 'M') curr_alignment->icigar[i] = M_replacement_character;
	}
	else
	{
		if (XS_tag_index == -1)
		{
			switch (curr_alignment->samflag)
			{
				case 0:
					M_replacement_character = 'B';
					break;
				case 16:
					M_replacement_character = 'E';
					break;
				case 256:
					M_replacement_character = 'F';
					break;
				case 272:
					M_replacement_character = 'H';
					break;
			}
		}
		else
		{
			if (strcmp(curr_alignment->tags[XS_tag_index].val, "+") == 0)
			{
				switch (curr_alignment->samflag)
				{
					case 0:
						M_replacement_character = 'J';
						break;
					case 16:
						M_replacement_character = 'K';
						break;
					case 256:
						M_replacement_character = 'L';
						break;
					case 272:
						M_replacement_character = 'O';
						break;
				}
			}
			else if (strcmp(curr_alignment->tags[XS_tag_index].val, "-") == 0)
			{
				switch (curr_alignment->samflag)
				{
					case 0:
						M_replacement_character = 'P';
						break;
					case 16:
						M_replacement_character = 'Q';
						break;
					case 256:
						M_replacement_character = 'R';
						break;
					case 272:
						M_replacement_character = 'U';
						break;
				}
			}

		}
		for (i = 0; i < strlen(curr_alignment->icigar); i++)
			if (curr_alignment->icigar[i] == 'M') curr_alignment->icigar[i] = M_replacement_character;
	}
	/*
	 switch (curr_alignment->samflag)
	 {
	 case 0:
	 strcat(curr_alignment->icigar, samformat_characters_SE[0]);
	 break;
	 case 16:
	 strcat(curr_alignment->icigar, samformat_characters_SE[1]);
	 break;
	 case 256:
	 strcat(curr_alignment->icigar, samformat_characters_SE[2]);
	 break;
	 case 272:
	 strcat(curr_alignment->icigar, samformat_characters_SE[3]);
	 break;
	 }*/

	/*
	 if (curr_alignment->samflag == 256)
	 {
	 printSamAlignmentInstance(curr_alignment);
	 printf("\nIntegrated CIGAR: %s CIGAR: %s MD: %s perfect_alignment_indicator: %d SEQ: %s READ_NAME: %s", curr_alignment->icigar, curr_alignment->cigar, curr_alignment->tags[MD_tag_index].val, perfect_alignment_indicator, curr_alignment->seq, curr_alignment->read_name);
	 }*/

	if (print_outputs > 0)
		printf("\nIntegrated CIGAR: %s CIGAR: %s MD: %s perfect_alignment_indicator: %d SEQ: %s READ_NAME: %s", curr_alignment->icigar, curr_alignment->cigar, curr_alignment->tags[MD_tag_index].val, perfect_alignment_indicator, curr_alignment->seq, curr_alignment->read_name);
	if (print_outputs > 1) printf("\n\n\n");
}

void initializePass3_Compression_Symbol_icigar_MappingPool(struct Pass3_Compression_Symbol_icigar_Mapping **symbol_icigar_mapping)
{
	int i;
	int j;
	int symbol_icigar_mapping_index = 0;
	char *individual_symbols = "+-0123456789abcdefghijklmnopqrstuvwxyzUVWXYZ";
	for (i = 0; i < MAX_SYMBOLS_FOR_PASS3_COMPRESSION; i++)
	{
		for (j = 0; j < MAX_SYMBOLS_FOR_PASS3_COMPRESSION; j++)
		{
			symbol_icigar_mapping[symbol_icigar_mapping_index]->symbolic_icigar[0] = individual_symbols[i];
			symbol_icigar_mapping[symbol_icigar_mapping_index]->symbolic_icigar[1] = individual_symbols[j];
			symbol_icigar_mapping[symbol_icigar_mapping_index]->symbolic_icigar[2] = '\0';
			symbol_icigar_mapping_index++;
		}
	}
	/*
	 printf("\n%d", symbol_icigar_mapping_index);
	 for (i = 0; i < symbol_icigar_mapping_index; i++)
	 printf("\n%s", symbol_icigar_mapping[i]->symbolic_icigar);
	 */
}

void assignicigarsToSymbols(struct Cigar_Frequency **cigar_freq_pool, int cigar_freq_pool_index, struct Pass3_Compression_Symbol_icigar_Mapping **symbol_icigar_mapping)
{
	int i;
	for (i = 0; i < cigar_freq_pool_index; i++)
		strcpy(symbol_icigar_mapping[i]->icigar, cigar_freq_pool[i]->cigar);
	/*
	 for (i = 0; i < cigar_freq_pool_index; i++)
	 printf("\n icigar: %s symbol: %s", symbol_icigar_mapping[i]->icigar, symbol_icigar_mapping[i]->symbolic_icigar);
	 */
}

#endif /* ABRIDGE_SRC_FUNCTIONS_DEFINITIONS_H_ */

