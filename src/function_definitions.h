#ifndef _ABRIDGE_SRC_FUNCTIONS_DEFINITIONS_H_
#define _ABRIDGE_SRC_FUNCTIONS_DEFINITIONS_H_

#include <stdbool.h>

struct Sam_Alignment* allocateMemorySam_Alignment ()
{
	/********************************************************************
	 * Variable declarations
	 ********************************************************************/
	struct Sam_Alignment *s;
	int i;
	/********************************************************************/
	s = ( struct Sam_Alignment* ) malloc (sizeof(struct Sam_Alignment));

	s->read_name = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->reference_name = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->cigar = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->reference_name_next_mate = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->seq = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->qual = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->soft_clippings.left = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->soft_clippings.right = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->soft_clippings.left_qual = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->soft_clippings.right_qual = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->cigar_items = ( struct Cigar_Items* ) malloc (sizeof(struct Cigar_Items) * MAX_CIGAR_ITEMS);
	s->number_of_cigar_items = 0;
	s->number_of_tag_items = 0;

	for ( i = 0 ; i < 50 ; i++ )
	{
		s->tags[i].name = ( char* ) malloc (sizeof(char) * MAX_TAG_LENGTH);
		s->tags[i].val = ( char* ) malloc (sizeof(char) * MAX_TAG_LENGTH);
		s->tags[i].type = ( char* ) malloc (sizeof(char) * MAX_TAG_LENGTH);
		s->tags[i].name[0] = '\0';
		s->tags[i].val[0] = '\0';
		s->tags[i].type[0] = '\0';
	}

	s->cigar_extended = ( char* ) malloc (sizeof(char) * ( MAX_SEQ_LEN * 2 ));
	s->md_extended = ( char* ) malloc (sizeof(char) * ( MAX_SEQ_LEN * 2 ));
	s->icigar = ( char* ) malloc (sizeof(char) * ( MAX_SEQ_LEN * 2 ));
	s->splices = ( char** ) malloc (sizeof(char*) * 100);
	for ( i = 0 ; i < 100 ; i++ )
		s->splices[i] = ( char* ) malloc (sizeof(char) * 50);
	s->soft_clips_removed_seq = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->soft_clips_removed_qual = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->selected_qual = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->temp = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);

	s->temp[0] = '\0';
	s->cigar_extended[0] = '\0';
	s->icigar[0] = '\0';
	s->md_extended[0] = '\0';
	s->soft_clips_removed_seq[0] = '\0';
	s->soft_clippings.left[0] = '\0';
	s->soft_clippings.right[0] = '\0';
	s->soft_clippings.left_qual[0] = '\0';
	s->soft_clippings.right_qual[0] = '\0';
	s->soft_clips_removed_qual[0] = '\0';
	s->soft_clips_removed_seq[0] = '\0';

	return s;
}

struct Paired_Ended_Flag_to_Single_Character* allocateMemoryPaired_Ended_Flag_to_Single_Character (int size)
{
	struct Paired_Ended_Flag_to_Single_Character *s;
	s = ( struct Paired_Ended_Flag_to_Single_Character* ) malloc (sizeof(struct Paired_Ended_Flag_to_Single_Character));
	s->character = ( char* ) malloc (sizeof(char) * size * 2);
	s->samflags = ( int* ) malloc (sizeof(int) * size * 2);
	s->direction = ( char* ) malloc (sizeof(int) * size * 2);
	return s;
}

struct Quality_Score_RLE* allocateMemoryQuality_Score_RLE ()
{
	struct Quality_Score_RLE *s;
	s = ( struct Quality_Score_RLE* ) malloc (sizeof(struct Quality_Score_RLE));
	s->score_character = 'X';
	s->frequency = 0;
	return s;
}

void reInitializeSamAlignmentInstance (struct Sam_Alignment *s)
{
	s->number_of_cigar_items = 0;
	s->number_of_tag_items = 0;
	s->temp[0] = '\0';
	s->cigar_extended[0] = '\0';
	s->icigar[0] = '\0';
	s->md_extended[0] = '\0';
	s->soft_clips_removed_seq[0] = '\0';
	s->soft_clippings.left[0] = '\0';
	s->soft_clippings.right[0] = '\0';

}

struct Pass3_Compression_Symbol_icigar_Mapping* allocateMemoryPass3_Compression_Symbol_icigar_Mapping ()
{
	struct Pass3_Compression_Symbol_icigar_Mapping *s;
	s = ( struct Pass3_Compression_Symbol_icigar_Mapping* ) malloc (sizeof(struct Pass3_Compression_Symbol_icigar_Mapping));
	s->icigar = ( char* ) malloc (sizeof(char) * MAX_ICIGAR_LENGTH);
	s->symbolic_icigar = ( char* ) malloc (sizeof(char) * MAX_SYMBOLIC_ICIGAR_LENGTH);
	return s;
}

struct Compressed_DS* allocateMemoryCompressed_DS (long long int max_input_reads_in_a_single_nucl_loc)
{
	struct Compressed_DS *s;
	int i;
	s = ( struct Compressed_DS* ) malloc (sizeof(struct Compressed_DS));
	s->icigar = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->pointers_to_qual_scores = ( char** ) malloc (sizeof(char*) * max_input_reads_in_a_single_nucl_loc);
	s->pointers_to_read_names = ( char** ) malloc (sizeof(char*) * max_input_reads_in_a_single_nucl_loc);
	s->num_reads = 0;
	s->position = 0;
	return s;
}

struct Pass1_Compressed_DS* allocateMemoryPass1_Compressed_DS ()
{
	struct Pass1_Compressed_DS *s;
	s = ( struct Pass1_Compressed_DS* ) malloc (sizeof(struct Pass1_Compressed_DS));
	s->col1 = ( char* ) malloc (sizeof(char) * ( MAX_SEQ_LEN * 2 ));
	s->col2 = ( char* ) malloc (sizeof(char) * MAX_ICIGAR_LENGTH_PASS1_COL2);
	s->col3 = ( char* ) malloc (sizeof(char) * ( MAX_SEQ_LEN * 2 ));
	return s;
}

struct Pass2_Compressed_DS* allocateMemoryPass2_Compressed_DS ()
{
	struct Pass2_Compressed_DS *s;
	s = ( struct Pass2_Compressed_DS* ) malloc (sizeof(struct Pass2_Compressed_DS));
	s->col1 = ( char* ) malloc (sizeof(char) * ( MAX_SEQ_LEN * 2 ));
	s->col2 = ( char* ) malloc (sizeof(char) * MAX_ICIGAR_LENGTH_PASS1_COL2);
	s->position = 0;
	return s;
}

struct Merged_Compressed_DS* allocateMemoryMerged_Compressed_DS ()
{
	struct Merged_Compressed_DS *s;
	int i;
	int j;

	s = ( struct Merged_Compressed_DS* ) malloc (sizeof(struct Merged_Compressed_DS));
	s->col1 = ( char* ) malloc (sizeof(char) * ( MAX_SEQ_LEN * 2 ));
	s->col2 = ( char* ) malloc (sizeof(char) * MAX_ICIGAR_LENGTH_PASS1_COL2);
	s->position = 0;
	s->icigars = ( char** ) malloc (sizeof(char*) * MAX_UNIQUE_CIGARS);
	for ( i = 0 ; i < MAX_UNIQUE_CIGARS ; i++ )
		s->icigars[i] = ( char* ) malloc (sizeof(char) * 1000);
	s->number_of_reads = ( int** ) malloc (sizeof(int*) * MAX_UNIQUE_CIGARS);
	for ( i = 0 ; i < MAX_UNIQUE_CIGARS ; i++ )
	{
		s->number_of_reads[i] = ( int* ) malloc (sizeof(int) * MAX_FILES_FOR_MERGING);
		for ( j = 0 ; j < MAX_FILES_FOR_MERGING ; j++ )
			s->number_of_reads[i][j] = 0;
	}
	s->number_of_unique_cigars = 0;
	return s;
}

struct Chromosome_Starting_Byte* allocateMemoryChromosome_Starting_Byte ()
{
	struct Chromosome_Starting_Byte *s;
	int i;

	s = ( struct Chromosome_Starting_Byte* ) malloc (sizeof(struct Chromosome_Starting_Byte));
	s->name = ( char** ) malloc (sizeof(char*) * MAX_REFERENCE_SEQUENCES);
	s->start_byte_in_pass2_file = ( long int* ) malloc (sizeof(long int) * MAX_REFERENCE_SEQUENCES);
	for ( i = 0 ; i < MAX_REFERENCE_SEQUENCES ; i++ )
		s->name[i] = ( char* ) malloc (sizeof(char) * 1000);
	s->number_of_chromosomes = 0;
	return s;
}

struct Chromosome_Info* allocateMemoryChromosome_Info ()
{
	struct Chromosome_Info *s;
	int i;

	s = ( struct Chromosome_Info* ) malloc (sizeof(struct Chromosome_Info));
	s->length = ( long long int* ) malloc (sizeof(long long int) * MAX_REFERENCE_SEQUENCES);
	s->number_of_chromosomes = 0;
	s->name = ( char** ) malloc (sizeof(char*) * MAX_REFERENCE_SEQUENCES);
	for ( i = 0 ; i < MAX_REFERENCE_SEQUENCES ; i++ )
		s->name[i] = ( char* ) malloc (sizeof(char) * 1000);
	return s;
}

struct Cigar_Frequency* allocateMemoryCigar_Frequency ()
{
	struct Cigar_Frequency *s;
	s = ( struct Cigar_Frequency* ) malloc (sizeof(struct Cigar_Frequency));
	s->cigar = ( char* ) malloc (sizeof(char) * MAX_CIGAR_LENGTH);
	s->freq = 0;
	return s;
}

struct Reference_Sequence_Info* allocateMemoryReference_Sequence_Info ()
{
	struct Reference_Sequence_Info *s;
	s = ( struct Reference_Sequence_Info* ) malloc (sizeof(struct Reference_Sequence_Info));
	s->line = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	return s;
}

struct Abridge_Index* allocateMemoryAbridge_Index ()
{
	int i;
	struct Abridge_Index *s;
	s = ( struct Abridge_Index* ) malloc (sizeof(struct Abridge_Index));
	s->chromosome = ( char** ) malloc (sizeof(char*) * MAX_POOL_SIZE);
	for ( i = 0 ; i < MAX_POOL_SIZE ; i++ )
		s->chromosome[i] = ( char* ) malloc (sizeof(char) * MAX_REFERENCE_SEQ_LENGTH);
	s->start = ( long long int* ) malloc (sizeof(long long int) * MAX_POOL_SIZE);
	s->end = ( long long int* ) malloc (sizeof(long long int) * MAX_POOL_SIZE);
	s->start_byte = ( long long int* ) malloc (sizeof(long long int) * MAX_POOL_SIZE);
	s->end_byte = ( long long int* ) malloc (sizeof(long long int) * MAX_POOL_SIZE);
	s->start_byte_qual = ( long long int* ) malloc (sizeof(long long int) * MAX_POOL_SIZE);
	s->end_byte_qual = ( long long int* ) malloc (sizeof(long long int) * MAX_POOL_SIZE);
	s->number_of_items = 0;
	return s;
}

struct Old_Read_ID_to_New_Read_ID_Circular_Linked_list* allocateMemoryOld_Read_ID_to_New_Read_ID_Circular_Linked_list ()
{
	struct Old_Read_ID_to_New_Read_ID_Circular_Linked_list *s;
	s = ( struct Old_Read_ID_to_New_Read_ID_Circular_Linked_list* ) malloc (sizeof(struct Old_Read_ID_to_New_Read_ID_Circular_Linked_list));
	s->old_read_id = ( char* ) malloc (sizeof(char) * 100);
	s->new_read_id = ( char* ) malloc (sizeof(char) * 100);
	s->number_of_multi_maps = 0;
	s->valid = 0;
	s->prev = NULL;
	s->next = NULL;
	return s;
}

struct All_Relevant_Info_PE_per_Alignment* allocateMemoryAll_Relevant_Info_PE_per_Alignment (int max_read_length)
{
	struct All_Relevant_Info_PE_per_Alignment *s;
	s = ( struct All_Relevant_Info_PE_per_Alignment* ) malloc (sizeof(struct All_Relevant_Info_PE_per_Alignment));
	s->NH_value = 0;
	s->icigar = ( char* ) malloc (sizeof(char) * max_read_length);
	s->new_read_id = ( char* ) malloc (sizeof(char) * 10);
	s->old_read_id = ( char* ) malloc (sizeof(char) * 50);
	s->position = 0;
	s->new_read_id_assigned = 0;
	return s;
}

long long int countNumberOfCharatersInString (char *line, char needle)
{
	/********************************************************************
	 * Variable declarations
	 ********************************************************************/
	long long int count = 0;

	int i;

	for ( i = 0 ; line[i] != '\0' ; i++ )
		if ( line[i] == needle ) count += 1;

	return count;

}

int isCommaInLine (char *line)
{
	int i;

	for ( i = 0 ; line[i] != '\0' ; i++ )
		if ( line[i] == ',' ) return 1;

	return 0;
}

int splitByDelimiter (char *line, char delimiter, char **new_string)
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

	for ( i = 0 ; line[i] != '\0' ; i++ )
	{
		// if space or NULL found, assign NULL into new_string[ctr]
		if ( line[i] == delimiter )
		{
			new_string[ctr][j] = '\0';
			ctr++; //for next word
			j = 0; //for next word, init index to 0
		}
		else if ( line[i] == '\n' )
			continue;
		else
		{
			new_string[ctr][j] = line[i];
			j++;
		}
	}
	new_string[ctr][j] = '\0';
	return ctr + 1;
}

int locateSamTags (char *tag)
{
	/********************************************************************
	 * Variable declarations
	 ********************************************************************/
	int i = 0;
	/********************************************************************/

	for ( i = 0 ; i < sizeof ( sam_tags ) / sizeof ( sam_tags[0] ) ; i++ )
	{
		if ( strcmp (sam_tags[i] , tag) == 0 ) return i;
	}
	return -1;

}

void freeSamAlignmentInstance (struct Sam_Alignment *s)
{
	int i;
	for ( i = 0 ; i < 100 ; i++ )
	{
		free (s->tags[i].name);
		free (s->tags[i].type);
		free (s->tags[i].val);
	}
	free (s->read_name);
	free (s->reference_name);
	free (s->reference_name_next_mate);
	free (s->cigar);
	free (s->seq);
	free (s->qual);
	free (s->soft_clippings.left);
	free (s->soft_clippings.right);
	free (s->cigar_items);
	free (s->cigar_extended);
	free (s->md_extended);
	free (s->icigar);
	for ( i = 0 ; i < ROWS ; i++ )
		free (s->splices[i]);
	free (s->temp);
}

void populateSamAlignmentInstance (struct Sam_Alignment *dest, char **src, int number_of_fields, char **split_tags)
{
	/********************************************************************
	 * Variable declarations
	 ********************************************************************/
	char *temp;
	int i;
	/********************************************************************/

	// Assign the first 11 mandatory sam fields
	strcpy (dest->read_name , src[0]);
	dest->samflag = strtol (src[1] , &temp , 10);
	strcpy (dest->reference_name , src[2]);
	dest->start_position = strtol (src[3] , &temp , 10);
	dest->mapping_quality_score = strtol (src[4] , &temp , 10);
	strcpy (dest->cigar , src[5]);
	strcpy (dest->reference_name_next_mate , src[6]);
	dest->start_position_next = strtol (src[7] , &temp , 10);
	dest->template_length = strtol (src[8] , &temp , 10);
	strcpy (dest->seq , src[9]);
	strcpy (dest->qual , src[10]);
	for ( i = 0 ; dest->qual[i] != '\0' ; i++ )
		dest->qual[i] += 90;

	// Assign SAM tags
	for ( i = 11 ; i < number_of_fields ; i++ )
	{
		//printf("%s\n",src[i]);
		//fflush ( stdout );
		splitByDelimiter (src[i] , ':' , split_tags);
		//sam_tag_index = locateSamTags(split_tags[0]);
		//dest->tags[i - 11].name = sam_tags[sam_tag_index];
		if ( i == number_of_fields - 1 )
			split_tags[2][strlen (split_tags[2])] = '\0';
		strcpy (dest->tags[i - 11].name , split_tags[0]);
		strcpy (dest->tags[i - 11].type , split_tags[1]);
		strcpy (dest->tags[i - 11].val , split_tags[2]);
		//printf("\n Tags %s Parts of the tag %s %s %s ", src[i], split_tags[0], split_tags[1], split_tags[2]);
	}
	dest->read_seq_len = strlen (dest->seq);
	dest->number_of_tag_items = number_of_fields - 11;
}

void printSamAlignmentInstance (struct Sam_Alignment *s, short int print_everything)
{
	int i;
	printf ("\n");
	printf ("^^^^^^^^^^^^^^^^^SAM ALIGNMENT^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
	printf ("\n");
	printf ("Read name: %s" , s->read_name);
	printf ("\n");
	printf ("SamFlag: %d" , s->samflag);
	printf ("\n");
	printf ("Reference name: %s" , s->reference_name);
	printf ("\n");
	printf ("Start position: %lld" , s->start_position);
	printf ("\n");
	printf ("Mapping quality score: %d" , s->mapping_quality_score);
	printf ("\n");
	printf ("CIGAR: %s" , s->cigar);
	printf ("\n");
	printf ("iCIGAR: %s" , s->icigar);
	printf ("\n");
	printf ("Reference name next mate: %s" , s->reference_name_next_mate);
	printf ("\n");
	printf ("Start position next: %lld" , s->start_position_next);
	printf ("\n");
	printf ("Template length: %d" , s->template_length);
	printf ("\n");
	printf ("Seq:  %s %d " , s->seq , strlen (s->seq));
	printf ("\n");
	printf ("Qual: %s");
	for ( i = 0 ; s->qual[i] != '\0' ; i++ )
		printf ("%c" , s->qual[i] - 90);
	printf ("\n");
	for ( i = 0 ; i < s->number_of_tag_items ; i++ )
	{
		printf ("Tag item %s Tag value %s" , s->tags[i].name , s->tags[i].val);
		printf ("\n");
	}

	if ( print_everything == 1 )
	{
		printf ("Left Soft clipped portion");
		printf ("\n");
		printf ("%s" , s->soft_clippings.left);
		printf ("\n");
		printf ("Right Soft clipped portion");
		printf ("\n");
		printf ("%s" , s->soft_clippings.right);
		printf ("\n");
		printf ("MD_extended: ");
		printf ("\n");
		printf ("CIGAR_extended: ");
		printf ("\n");
		printf ("soft clip removed sequence: ");
		printf ("\n");
		printf ("\n%s\n%s\n%s" , s->md_extended , s->cigar_extended , s->soft_clips_removed_seq);
		printf ("\n");
	}
	printf ("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
}

int isSequenceSoftClipped (char *cigar)
{
	int i;
	for ( i = 0 ; cigar[i] != '\0' ; i++ )
		if ( cigar[i] == 'S' ) return 1;
	return 0;
}

void convertIcigarToCigarandMDPairedEnded (struct Whole_Genome_Sequence *whole_genome, struct Sam_Alignment *sam_alignment_instance, char *chromosome, short int flag_ignore_mismatches, short int flag_ignore_soft_clippings, short int flag_ignore_unmapped_sequences, short int flag_ignore_quality_score, short int flag_ignore_sequence_information, char *default_quality_value, struct Paired_Ended_Flag_to_Single_Character *samflag_dictionary, int number_of_unique_samformatflags, char samformatflag_replacer_characters[])
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
	int icigar_items_instance_index;
	int cigar_items_instance_index;
	int flag = 0;
	int soft_clips_removed_qual_index;
	int soft_clips_removed_seq_index;
	int MD_index = 0;
	int MD_extended_index = 0;

	char XS[2];
	char MD[1000];
	char character_to_be_replaced;

	struct Cigar_Items cigar_items_instance[MAX_SEQ_LEN];

	char temp[100];
	char cigar_temp[MAX_SEQ_LEN];

	icigar_length = strlen (sam_alignment_instance->icigar);
	NH_value = extractNHfromicigar (sam_alignment_instance->icigar , icigar_length);

	XS[1] = '\0';
	samformatflag = findSamFormatFlagPairedEnded (sam_alignment_instance->icigar , icigar_length , XS , samflag_dictionary , number_of_unique_samformatflags , samformatflag_replacer_characters , &character_to_be_replaced);
	sam_alignment_instance->samflag = samformatflag;

	/*
	 * Construct the Cigar string, MD string, Soft Clips, etc.
	 */
	splitCigar (sam_alignment_instance->icigar , &icigar_items_instance_index , cigar_items_instance);

	strcpy (sam_alignment_instance->reference_name , chromosome);
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

	strcpy (sam_alignment_instance->tags[0].name , "NH");
	strcpy (sam_alignment_instance->tags[0].type , "i");
	sprintf (temp , "%d" , NH_value);
	strcpy (sam_alignment_instance->tags[0].val , temp);

	strcpy (sam_alignment_instance->tags[1].name , "XS");
	strcpy (sam_alignment_instance->tags[1].type , ".");
	strcpy (sam_alignment_instance->tags[1].val , XS);

	strcpy (sam_alignment_instance->tags[2].name , "MD");
	strcpy (sam_alignment_instance->tags[2].type , "Z");
	strcpy (sam_alignment_instance->tags[2].val , "dummy");

	sam_alignment_instance->number_of_tag_items = 3;
	sam_alignment_instance->md_extended[0] = '\0';
	MD_extended_index = 0;

	for ( i = 0 ; i < icigar_items_instance_index ; i++ )
	{
		if ( processing_left_soft_clip == 1 && isCharacterInString ("atgcn" , cigar_items_instance[i].def) )
		{
			sam_alignment_instance->soft_clippings.left[left_soft_clip_index] = cigar_items_instance[i].def - 32;
			sam_alignment_instance->soft_clippings.left[left_soft_clip_index + 1] = '\0';
			if ( flag_ignore_quality_score == 0 )
			{
				i++;
				sam_alignment_instance->soft_clippings.left_qual[left_soft_clip_index] = cigar_items_instance[i].def - 90;
				sam_alignment_instance->soft_clippings.left_qual[left_soft_clip_index + 1] = '\0';
			}
			left_soft_clip_index++;
		}
		else if ( processing_left_soft_clip == 0 && isCharacterInString ("atgcn" , cigar_items_instance[i].def) )
		{
			sam_alignment_instance->soft_clippings.right[right_soft_clip_index] = cigar_items_instance[i].def - 32;
			sam_alignment_instance->soft_clippings.right[right_soft_clip_index + 1] = '\0';
			if ( flag_ignore_quality_score == 0 )
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
			if ( left_soft_clip_index > 0 ) //There were some left soft clips
			{
				sprintf (temp , "%d" , left_soft_clip_index);
				strcat (sam_alignment_instance->cigar , temp);
				strcat (sam_alignment_instance->cigar , "S");
				num = 0;
				left_soft_clip_index = 0;
			}
			if ( cigar_items_instance[i].def == 'N' )
			{
				sprintf (temp , "%d" , cigar_items_instance[i].len);
				strcat (sam_alignment_instance->cigar , temp);
				strcat (sam_alignment_instance->cigar , "N");
				num = 0;
			}
			else if ( cigar_items_instance[i].def == 'D' )
			{
				sprintf (temp , "%d" , cigar_items_instance[i].len);
				strcat (sam_alignment_instance->cigar , temp);
				strcat (sam_alignment_instance->cigar , "D");
				num = 0;

				for ( j = 0 ; j < cigar_items_instance[i].len ; j++ )
					sam_alignment_instance->md_extended[MD_extended_index++ ] = '-';
			}
			else if ( isCharacterInString (insert_characters , cigar_items_instance[i].def) )
			{
				num = 0;
				while ( i < icigar_items_instance_index && isCharacterInString (insert_characters , cigar_items_instance[i].def) )
				{
					switch ( cigar_items_instance[i].def )
					{
						case '!':
							sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++ ] = 'A';
							sam_alignment_instance->md_extended[MD_extended_index++ ] = '!';
							break;
						case '"':
							sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++ ] = 'T';
							sam_alignment_instance->md_extended[MD_extended_index++ ] = '"';
							break;
						case '#':
							sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++ ] = 'G';
							sam_alignment_instance->md_extended[MD_extended_index++ ] = '#';
							break;
						case '$':
							sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++ ] = 'C';
							sam_alignment_instance->md_extended[MD_extended_index++ ] = '$';
							break;
						case '%':
							sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++ ] = 'N';
							sam_alignment_instance->md_extended[MD_extended_index++ ] = '%';
							break;
					}
					if ( flag_ignore_quality_score == 0 )
					{
						i++;
						sam_alignment_instance->soft_clips_removed_qual[soft_clips_removed_qual_index] = cigar_items_instance[i].def - 90;
						sam_alignment_instance->soft_clips_removed_qual[soft_clips_removed_qual_index + 1] = '\0';
						soft_clips_removed_qual_index++;
					}
					num++;
					i++;
				}
				i--;
				sprintf (temp , "%d" , num);
				strcat (sam_alignment_instance->cigar , temp);
				strcat (sam_alignment_instance->cigar , "I");
				num = 0;

			}
			// Merging consecutive matches and mismatches
			if ( isCharacterInString (samformatflag_replacer_characters , cigar_items_instance[i].def) || isCharacterInString (mismatch_characters , cigar_items_instance[i].def) )
			{
				num = 0;
				while ( i < icigar_items_instance_index && ( isCharacterInString (samformatflag_replacer_characters , cigar_items_instance[i].def) || isCharacterInString (mismatch_characters , cigar_items_instance[i].def) ) )
				{
					if ( isCharacterInString (samformatflag_replacer_characters , cigar_items_instance[i].def) )
					{
						for ( j = 0 ; j < cigar_items_instance[i].len ; j++ )
						{
							sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++ ] = '-';
							sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index] = '\0';
							sam_alignment_instance->md_extended[MD_extended_index++ ] = '=';
						}

						if ( flag_ignore_quality_score == 0 )
						{
							for ( j = 0 ; j < cigar_items_instance[i].len ;
									j++ )
							{
								sam_alignment_instance->soft_clips_removed_qual[soft_clips_removed_qual_index] = default_quality_value[0];
								sam_alignment_instance->soft_clips_removed_qual[soft_clips_removed_qual_index + 1] = '\0';
								soft_clips_removed_qual_index++;
							}
						}
						num += cigar_items_instance[i].len;
					}
					else if ( isCharacterInString (mismatch_characters , cigar_items_instance[i].def) )
					{
						if ( isCharacterInString (mismatch_characters , cigar_items_instance[i].def) )
						{
							switch ( cigar_items_instance[i].def )
							{
								case '&':
									sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++ ] = 'A';
									sam_alignment_instance->md_extended[MD_extended_index++ ] = '&';
									break;
								case '\'':
									sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++ ] = 'T';
									sam_alignment_instance->md_extended[MD_extended_index++ ] = '\'';
									break;
								case '(':
									sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++ ] = 'G';
									sam_alignment_instance->md_extended[MD_extended_index++ ] = '(';
									break;
								case ')':
									sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++ ] = 'C';
									sam_alignment_instance->md_extended[MD_extended_index++ ] = ')';
									break;
								case '*':
									sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++ ] = 'N';
									sam_alignment_instance->md_extended[MD_extended_index++ ] = '*';
									break;
							}
						}
						if ( flag_ignore_quality_score == 0 )
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
				i--;
				sprintf (temp , "%d" , num);
				strcat (sam_alignment_instance->cigar , temp);
				strcat (sam_alignment_instance->cigar , "M");
				num = 0;
			}
		}
	}
	if ( right_soft_clip_index > 0 ) //There were some right soft clips
	{
		sprintf (temp , "%d" , right_soft_clip_index);
		strcat (sam_alignment_instance->cigar , temp);
		strcat (sam_alignment_instance->cigar , "S");
		num = 0;
	}
	sam_alignment_instance->md_extended[MD_extended_index++ ] = '\0';

	strcpy (sam_alignment_instance->seq , sam_alignment_instance->soft_clippings.left);
	strcat (sam_alignment_instance->seq , sam_alignment_instance->soft_clips_removed_seq);
	strcat (sam_alignment_instance->seq , sam_alignment_instance->soft_clippings.right);

	strcpy (sam_alignment_instance->qual , sam_alignment_instance->soft_clippings.left_qual);
	strcat (sam_alignment_instance->qual , sam_alignment_instance->soft_clips_removed_qual);
	strcat (sam_alignment_instance->qual , sam_alignment_instance->soft_clippings.right_qual);

	splitCigar (sam_alignment_instance->cigar , &sam_alignment_instance->number_of_cigar_items , sam_alignment_instance->cigar_items);
	if ( flag_ignore_sequence_information == 0 )
		generateReadSequenceAndMDString (sam_alignment_instance , whole_genome);
	else
	{
		splitCigar (sam_alignment_instance->cigar , &cigar_items_instance_index , cigar_items_instance);
		int length_of_read = 0;
		for ( i = 0 ; i < cigar_items_instance_index ; i++ )
			if ( cigar_items_instance[cigar_items_instance_index].def == 'M' )
				length_of_read += cigar_items_instance[cigar_items_instance_index].len;
		for ( i = 0 ; i < length_of_read ; i++ )
		{
			sam_alignment_instance->seq[i] = 'A';
			sam_alignment_instance->qual[i] = default_quality_value[0];
		}
		sam_alignment_instance->seq[i] = '\0';
		sam_alignment_instance->qual[i] = '\0';
	}

}

void convertIcigarToCigarandMDSingleEnded (struct Whole_Genome_Sequence *whole_genome, struct Sam_Alignment *sam_alignment_instance, char *chromosome, short int flag_ignore_mismatches, short int flag_ignore_soft_clippings, short int flag_ignore_unmapped_sequences, short int flag_ignore_quality_score, short int flag_ignore_sequence_information, char *default_quality_value)
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
	int icigar_items_instance_index;
	int cigar_items_instance_index;
	int flag = 0;
	int soft_clips_removed_qual_index;
	int soft_clips_removed_seq_index;
	int MD_index = 0;
	int MD_extended_index = 0;

	char XS[2];
	char MD[1000];

	struct Cigar_Items cigar_items_instance[MAX_SEQ_LEN];

	char temp[100];
	char cigar_temp[MAX_SEQ_LEN];

	icigar_length = strlen (sam_alignment_instance->icigar);
	NH_value = extractNHfromicigar (sam_alignment_instance->icigar , icigar_length);

	XS[1] = '\0';
	samformatflag = findSamFormatFlagSingleEnded (sam_alignment_instance->icigar , icigar_length , XS);
	sam_alignment_instance->samflag = samformatflag;
	/*
	 * Construct the Cigar string, MD string, Soft Clips, etc.
	 */
	splitCigar (sam_alignment_instance->icigar , &icigar_items_instance_index , cigar_items_instance);

	strcpy (sam_alignment_instance->reference_name , chromosome);
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

	strcpy (sam_alignment_instance->tags[0].name , "NH");
	strcpy (sam_alignment_instance->tags[0].type , "i");
	sprintf (temp , "%d" , NH_value);
	strcpy (sam_alignment_instance->tags[0].val , temp);

	strcpy (sam_alignment_instance->tags[1].name , "XS");
	strcpy (sam_alignment_instance->tags[1].type , ".");
	strcpy (sam_alignment_instance->tags[1].val , XS);

	strcpy (sam_alignment_instance->tags[2].name , "MD");
	strcpy (sam_alignment_instance->tags[2].type , "Z");
	strcpy (sam_alignment_instance->tags[2].val , "dummy");

	sam_alignment_instance->number_of_tag_items = 3;
	sam_alignment_instance->md_extended[0] = '\0';
	MD_extended_index = 0;

	for ( i = 0 ; i < icigar_items_instance_index ; i++ )
	{
		if ( processing_left_soft_clip == 1 && isCharacterInString ("atgcn" , cigar_items_instance[i].def) )
		{
			sam_alignment_instance->soft_clippings.left[left_soft_clip_index] = cigar_items_instance[i].def - 32;
			sam_alignment_instance->soft_clippings.left[left_soft_clip_index + 1] = '\0';
			if ( flag_ignore_quality_score == 0 )
			{
				i++;
				sam_alignment_instance->soft_clippings.left_qual[left_soft_clip_index] = cigar_items_instance[i].def - 90;
				sam_alignment_instance->soft_clippings.left_qual[left_soft_clip_index + 1] = '\0';
			}
			left_soft_clip_index++;
		}
		else if ( processing_left_soft_clip == 0 && isCharacterInString ("atgcn" , cigar_items_instance[i].def) )
		{
			sam_alignment_instance->soft_clippings.right[right_soft_clip_index] = cigar_items_instance[i].def - 32;
			sam_alignment_instance->soft_clippings.right[right_soft_clip_index + 1] = '\0';
			if ( flag_ignore_quality_score == 0 )
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
			if ( left_soft_clip_index > 0 ) //There were some left soft clips
			{
				sprintf (temp , "%d" , left_soft_clip_index);
				strcat (sam_alignment_instance->cigar , temp);
				strcat (sam_alignment_instance->cigar , "S");
				num = 0;
				left_soft_clip_index = 0;
			}
			if ( cigar_items_instance[i].def == 'N' )
			{
				sprintf (temp , "%d" , cigar_items_instance[i].len);
				strcat (sam_alignment_instance->cigar , temp);
				strcat (sam_alignment_instance->cigar , "N");
				num = 0;
			}
			else if ( cigar_items_instance[i].def == 'D' )
			{
				sprintf (temp , "%d" , cigar_items_instance[i].len);
				strcat (sam_alignment_instance->cigar , temp);
				strcat (sam_alignment_instance->cigar , "D");
				num = 0;

				for ( j = 0 ; j < cigar_items_instance[i].len ; j++ )
					sam_alignment_instance->md_extended[MD_extended_index++ ] = '-';
			}
			else if ( isCharacterInString (insert_characters , cigar_items_instance[i].def) )
			{
				num = 0;
				while ( i < icigar_items_instance_index && isCharacterInString (insert_characters , cigar_items_instance[i].def) )
				{
					switch ( cigar_items_instance[i].def )
					{
						case '!':
							sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++ ] = 'A';
							sam_alignment_instance->md_extended[MD_extended_index++ ] = '!';
							break;
						case '"':
							sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++ ] = 'T';
							sam_alignment_instance->md_extended[MD_extended_index++ ] = '"';
							break;
						case '#':
							sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++ ] = 'G';
							sam_alignment_instance->md_extended[MD_extended_index++ ] = '#';
							break;
						case '$':
							sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++ ] = 'C';
							sam_alignment_instance->md_extended[MD_extended_index++ ] = '$';
							break;
						case '%':
							sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++ ] = 'N';
							sam_alignment_instance->md_extended[MD_extended_index++ ] = '%';
							break;
					}
					if ( flag_ignore_quality_score == 0 )
					{
						i++;
						sam_alignment_instance->soft_clips_removed_qual[soft_clips_removed_qual_index] = cigar_items_instance[i].def - 90;
						sam_alignment_instance->soft_clips_removed_qual[soft_clips_removed_qual_index + 1] = '\0';
						soft_clips_removed_qual_index++;
					}
					num++;
					i++;
				}
				i--;
				sprintf (temp , "%d" , num);
				strcat (sam_alignment_instance->cigar , temp);
				strcat (sam_alignment_instance->cigar , "I");
				num = 0;

			}
			// Merging consecutive matches and mismatches
			if ( isCharacterInString ("BEFHJKLOPQRU" , cigar_items_instance[i].def) || isCharacterInString (mismatch_characters , cigar_items_instance[i].def) )
			{
				num = 0;
				while ( i < icigar_items_instance_index && ( isCharacterInString ("BEFHJKLOPQRU" , cigar_items_instance[i].def) || isCharacterInString (mismatch_characters , cigar_items_instance[i].def) ) )
				{
					if ( isCharacterInString ("BEFHJKLOPQRU" , cigar_items_instance[i].def) )
					{
						for ( j = 0 ; j < cigar_items_instance[i].len ; j++ )
						{
							sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++ ] = '-';
							sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index] = '\0';
							sam_alignment_instance->md_extended[MD_extended_index++ ] = '=';
						}

						if ( flag_ignore_quality_score == 0 )
						{
							for ( j = 0 ; j < cigar_items_instance[i].len ;
									j++ )
							{
								sam_alignment_instance->soft_clips_removed_qual[soft_clips_removed_qual_index] = default_quality_value[0];
								sam_alignment_instance->soft_clips_removed_qual[soft_clips_removed_qual_index + 1] = '\0';
								soft_clips_removed_qual_index++;
							}
						}
						num += cigar_items_instance[i].len;
					}
					else if ( isCharacterInString (mismatch_characters , cigar_items_instance[i].def) )
					{
						if ( isCharacterInString (mismatch_characters , cigar_items_instance[i].def) )
						{
							switch ( cigar_items_instance[i].def )
							{
								case '&':
									sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++ ] = 'A';
									sam_alignment_instance->md_extended[MD_extended_index++ ] = '&';
									break;
								case '\'':
									sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++ ] = 'T';
									sam_alignment_instance->md_extended[MD_extended_index++ ] = '\'';
									break;
								case '(':
									sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++ ] = 'G';
									sam_alignment_instance->md_extended[MD_extended_index++ ] = '(';
									break;
								case ')':
									sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++ ] = 'C';
									sam_alignment_instance->md_extended[MD_extended_index++ ] = ')';
									break;
								case '*':
									sam_alignment_instance->soft_clips_removed_seq[soft_clips_removed_seq_index++ ] = 'N';
									sam_alignment_instance->md_extended[MD_extended_index++ ] = '*';
									break;
							}
						}
						if ( flag_ignore_quality_score == 0 )
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
				i--;
				sprintf (temp , "%d" , num);
				strcat (sam_alignment_instance->cigar , temp);
				strcat (sam_alignment_instance->cigar , "M");
				num = 0;
			}
		}
	}
	if ( right_soft_clip_index > 0 ) //There were some right soft clips
	{
		sprintf (temp , "%d" , right_soft_clip_index);
		strcat (sam_alignment_instance->cigar , temp);
		strcat (sam_alignment_instance->cigar , "S");
		num = 0;
	}
	sam_alignment_instance->md_extended[MD_extended_index++ ] = '\0';

	strcpy (sam_alignment_instance->seq , sam_alignment_instance->soft_clippings.left);
	strcat (sam_alignment_instance->seq , sam_alignment_instance->soft_clips_removed_seq);
	strcat (sam_alignment_instance->seq , sam_alignment_instance->soft_clippings.right);

	strcpy (sam_alignment_instance->qual , sam_alignment_instance->soft_clippings.left_qual);
	strcat (sam_alignment_instance->qual , sam_alignment_instance->soft_clips_removed_qual);
	strcat (sam_alignment_instance->qual , sam_alignment_instance->soft_clippings.right_qual);

	splitCigar (sam_alignment_instance->cigar , &sam_alignment_instance->number_of_cigar_items , sam_alignment_instance->cigar_items);
	if ( flag_ignore_sequence_information == 0 )
		generateReadSequenceAndMDString (sam_alignment_instance , whole_genome);
	else
	{
		splitCigar (sam_alignment_instance->cigar , &cigar_items_instance_index , cigar_items_instance);
		int length_of_read = 0;
		for ( i = 0 ; i < cigar_items_instance_index ; i++ )
			if ( cigar_items_instance[cigar_items_instance_index].def == 'M' )
				length_of_read += cigar_items_instance[cigar_items_instance_index].len;
		for ( i = 0 ; i < length_of_read ; i++ )
		{
			sam_alignment_instance->seq[i] = 'A';
			sam_alignment_instance->qual[i] = default_quality_value[0];
		}
		sam_alignment_instance->seq[i] = '\0';
		sam_alignment_instance->qual[i] = '\0';
	}
}

void extractSubString (char *str, char *substr, int start_index, int end_index)
{
	/********************************************************************
	 * Variable declarations
	 ********************************************************************/
	int length = strlen (str);
	int i;
	int j = 0;

	/********************************************************************/

	if ( length < ( end_index - start_index + 1 ) )
	{
		return;
	}
	for ( i = start_index ; i <= end_index ; i++ )
		substr[j++ ] = str[i];
	substr[j] = '\0';

}

void splitCigar (char *cigar, int *num_of_types, struct Cigar_Items *cigar_items_instance)
{
	/********************************************************************
	 * Variable declarations
	 ********************************************************************/

	int i, j = 0;
	int current_length = 0;
	/********************************************************************/

	for ( i = 0 ; cigar[i] != '\0' ; i++ )
	{
		if ( isdigit (cigar[i]) != 0 ) //cigar[i] is a digit
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

int isCharacterInString (char *str, char character)
{
	int i;
	for ( i = 0 ; str[i] != '\0' ; i++ )
		if ( character == str[i] ) return 1;
	return 0;
}

int isNumber (char *str)
{
	int i;
	for ( i = 0 ; str[i] != '\0' ; i++ )
		if ( str[i] < 48 || str[i] > 57 ) return 0;
	return 1;
}

void insertCharacterInString (char *str, int *str_length, char ins, int loc, int number_of_insertions_to_be_made)
{
	/*
	 * Inserts a character. Assumes that a large
	 */
	int i = ( *str_length ) + number_of_insertions_to_be_made;
	for ( ; i > loc ; i-- )
		str[i] = str[i - number_of_insertions_to_be_made];
	( *str_length ) += number_of_insertions_to_be_made;
	while ( number_of_insertions_to_be_made-- )
		str[i++ ] = ins;
}

void replaceNucleotidesInSeqWithInsertSymbols (char *seq, int *seq_length, struct Cigar_Items *cigar_items_instance, int number_of_cigar_items)
{
	int i;
	int j;
	int sequence_index_to_be_replaced = 0;
	for ( i = 0 ; i <= number_of_cigar_items ; i++ )
	{
		switch ( cigar_items_instance[i].def )
		{
			case 'M':
				sequence_index_to_be_replaced += cigar_items_instance[i].len;
				break;
			case 'I':
				for ( j = 0 ; j < cigar_items_instance[i].len ; j++ )
				{
					switch ( seq[sequence_index_to_be_replaced] )
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
				//sequence_index_to_be_replaced += cigar_items_instance[i].len;
				break;
		}
	}
}

void convertRegularCIGARToStringRepresentation (struct Cigar_Items *cigar_items_instance, int number_of_cigar_items, char *cigar_extended)
{
	int i;
	int j;
	int cigar_extended_index = 0;
	cigar_extended[0] = '\0';
	for ( i = 0 ; i < number_of_cigar_items ; i++ )
	{
		switch ( cigar_items_instance[i].def )
		{
			case 'M':
				for ( j = 0 ; j < cigar_items_instance[i].len ; j++ )
					cigar_extended[cigar_extended_index++ ] = 'M';
				break;
			case 'I':
				for ( j = 0 ; j < cigar_items_instance[i].len ; j++ )
					cigar_extended[cigar_extended_index++ ] = 'I';
				break;
		}
	}
	cigar_extended[cigar_extended_index++ ] = '\0';
}

void designIntegratedCIGARSingleEnded (char *md, char *seq, int *seq_length, char *soft_clips_removed_qual, struct Cigar_Items *cigar_items_instance, int number_of_cigar_items, char *cigar, char *cigar_extended, char *md_extended, char *icigar, char **splices, short int flag_ignore_quality_score, short int flag_ignore_mismatches)
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
	int soft_clips_removed_qual_len = 0;

	md_extended_index = 0;
	i = 0;
	num = 0;
	replaceNucleotidesInSeqWithInsertSymbols (seq , seq_length , cigar_items_instance , number_of_cigar_items);
	convertRegularCIGARToStringRepresentation (cigar_items_instance , number_of_cigar_items , cigar_extended);
	soft_clips_removed_qual_len = strlen (soft_clips_removed_qual);

	/*
	 * Expand the md string
	 */
	i = 0;
	while ( md[i] != '\0' )
	{
		if ( isdigit (md[i]) != 0 ) // md[i] is a digit
			num = num * 10 + md[i] - 48;
		else if ( md[i] == '^' )
		{
			for ( j = 0 ; j < num ; j++ )
				md_extended[md_extended_index++ ] = '=';
			num = 0;
			i += 1;
			while ( isdigit (md[i]) == 0 ) // Iterate till you pick up all the deletions and ignore them for now
				i += 1;
			i -= 1;
		}
		else
		{
			for ( j = 0 ; j < num ; j++ )
				md_extended[md_extended_index++ ] = '=';
			md_extended[md_extended_index++ ] = 'X';
			num = 0;
		}
		i += 1;
	}
	for ( j = 0 ; j < num ; j++ )
		md_extended[md_extended_index++ ] = '=';
	md_extended[md_extended_index++ ] = '\0';
	md_extended_length = md_extended_index - 1;

	/*
	 * Add insert symbols to md_extended
	 */
	j = 0;
	for ( i = 0 ; i < *seq_length ; i++ )
	{
		if ( seq[i] >= 33 && seq[i] <= 37 )
		{
			insertCharacterInString (md_extended , &md_extended_length , seq[i] , i , 1);
			j = -1;
		}
	}

	if ( strlen (seq) != strlen (md_extended) )
	{
		printf ("\nInside here");
		printf ("\nLengths do not match %d %d %d %s" , strlen (soft_clips_removed_qual) , strlen (md_extended) , strlen (seq) , cigar);
		printf ("\n%s" , seq);
		printf ("\n");
		for ( i = 0 ; soft_clips_removed_qual[i] != '\0' ; i++ )
			printf ("%c" , soft_clips_removed_qual[i] - 90);
		printf ("\n%s" , md_extended);
	}
	/*printf("\n%s %d", md_extended, md_extended_length);*/

	/*
	 * Add mismatch characters to md_extended
	 */

	for ( i = 0 ; i < md_extended_length ; i++ )
	{
		if ( md_extended[i] == 'X' )
		{
			switch ( seq[i] )
			{
				case 'A':
					md_extended[i] =
							flag_ignore_mismatches == 0 ? mismatch_characters[0] : '=';
					break;
				case 'T':
					md_extended[i] =
							flag_ignore_mismatches == 0 ? mismatch_characters[1] : '=';
					break;
				case 'G':
					md_extended[i] =
							flag_ignore_mismatches == 0 ? mismatch_characters[2] : '=';
					break;
				case 'C':
					md_extended[i] =
							flag_ignore_mismatches == 0 ? mismatch_characters[3] : '=';
					break;
				case 'N':
					md_extended[i] =
							flag_ignore_mismatches == 0 ? mismatch_characters[4] : '=';
					break;
			}
		}
	}

	/*
	 * Encode deletions and splice junction info in md_extended
	 */
	md_extended_index = 0;
	splice_number_index = 0;
	for ( i = 0 ; i < number_of_cigar_items ; i++ )
	{
		switch ( cigar_items_instance[i].def )
		{
			case 'M':
				md_extended_index += cigar_items_instance[i].len;
				break;
			case 'N':
				//md_extended[md_extended_index] = splice_number_index++ + 48;
				insertCharacterInString (md_extended , &md_extended_length , splice_number_index++ + 48 , md_extended_index , 1);
				insertCharacterInString (soft_clips_removed_qual , &soft_clips_removed_qual_len , 'X' , md_extended_index , 1);
				md_extended_index += 1;
				break;
			case 'I':
				md_extended_index += cigar_items_instance[i].len;
				break;
			case 'D':
				insertCharacterInString (md_extended , &md_extended_length , '-' , md_extended_index , cigar_items_instance[i].len);
				insertCharacterInString (soft_clips_removed_qual , &soft_clips_removed_qual_len , '-' , md_extended_index , cigar_items_instance[i].len);
				md_extended_index += cigar_items_instance[i].len;
				break;
		}
	}

	/*
	 * Generate the informative CIGAR
	 */

	/*
	 * Combine the splice items into an array
	 */
	splices_index = 0;
	for ( i = 0 ; i < number_of_cigar_items ; i++ )
	{
		if ( cigar_items_instance[i].def == 'N' )
		{
			sprintf (temp_convert_int_to_string , "%d" , cigar_items_instance[i].len);
			strcpy (splices[splices_index] , temp_convert_int_to_string);
			strcat (splices[splices_index] , "N");
			splices_index++;
		}
	}

	/*
	 * Prepare the icigar string
	 */
	num = 0;
	icigar[0] = '\0';
	intron_splice_index = 0;
	for ( i = 0 ; md_extended[i] != '\0' ; i++ )
	{
		if ( md_extended[i] == '=' )
			num++;
		else
		{
			sprintf (temp_convert_int_to_string , "%d" , num);
			if ( num > 0 )
			{
				strcat (icigar , temp_convert_int_to_string);
				strcat (icigar , "M");
			}
			num = 0;
			if ( ( md_extended[i] >= 33 && md_extended[i] <= 37 ) || ( md_extended[i] >= 38 && md_extended[i] <= 42 ) ) // Either a mismatch or an insertion
			{

				temp_convert_int_to_string[0] = md_extended[i];
				if ( flag_ignore_quality_score == 0 )
				{
					temp_convert_int_to_string[1] = soft_clips_removed_qual[i];
					temp_convert_int_to_string[2] = '\0';
				}
				else temp_convert_int_to_string[1] = '\0';
				strcat (icigar , temp_convert_int_to_string);
			}
			else if ( md_extended[i] >= '0' && md_extended[i] <= '9' )
			{
				num = 0;
				intron_splice_index = 0;
				while ( md_extended[i] >= '0' && md_extended[i] <= '9' )
				{
					intron_splice_index = intron_splice_index * 10 + md_extended[i] - 48;
					i++;
					num++;
				}
				num = 0;
				strcat (icigar , splices[intron_splice_index]);
				i--;
			}
			else if ( md_extended[i] == '-' ) // Deletion
			{
				while ( md_extended[i] == '-' )
				{
					num++;
					i++;
				}
				i--;
				sprintf (temp_convert_int_to_string , "%d" , num);
				strcat (icigar , temp_convert_int_to_string);
				strcat (icigar , "D");
				num = 0;
			}
		}
	}
	sprintf (temp_convert_int_to_string , "%d" , num);
	strcat (icigar , temp_convert_int_to_string);
	strcat (icigar , "M");
}

int isAlignmentPerfect (char *cigar, struct Sam_Tags *tags, int MD_tag_index, int NM_tag_index, int nM_tag_index)
{
	int i;
	if ( NM_tag_index != -1 && nM_tag_index != -1 )
	{
		if ( tags[NM_tag_index].val[0] == '0' && tags[nM_tag_index].val[0] == '0' )
			return 1;
		else return 0;
	}

	/*
	 * Checking for insertions
	 */
	for ( i = 0 ; cigar[i] != '\0' ; i++ )
		if ( cigar[i] == 'I' ) return 0;

	/*
	 * Checking for deletions
	 */
	for ( i = 0 ; tags[MD_tag_index].val[i] != '\0' ; i++ )
		if ( tags[MD_tag_index].val[i] == '^' || ( tags[MD_tag_index].val[i] >= 65 && tags[MD_tag_index].val[i] <= 90 ) || ( tags[MD_tag_index].val[i] >= 97 && tags[MD_tag_index].val[i] <= 122 ) )
			return 0;

	return 1;
}

void generateIntegratedCigarSingleEnded (struct Sam_Alignment *curr_alignment, short int flag_ignore_soft_clippings, short int flag_ignore_mismatches, short int flag_ignore_unmapped_sequences, short int flag_ignore_quality_score, struct Whole_Genome_Sequence *whole_genome, struct Sam_Alignment *sam_alignment_instance_diagnostics, long long int number_of_records_read, short int run_diagnostics)
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
	int flag;
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
	char dummy;
	/********************************************************************/

	splitCigar (curr_alignment->cigar , &curr_alignment->number_of_cigar_items , curr_alignment->cigar_items);

	for ( i = 0 ; i < curr_alignment->number_of_cigar_items ; i++ )
	{
		if ( curr_alignment->cigar_items[i].def == 'N' )
		{
			spliced_alignment_indicator = 1;
			break;
		}
	}
	/*
	 * Process each alignment to extract soft clipped portion of reads
	 */
	if ( isSequenceSoftClipped (curr_alignment->cigar) == 1 )
	{
		if ( curr_alignment->cigar_items[0].def == 'S' )
		{
			left_soft_clip_point = curr_alignment->cigar_items[0].len;
			extractSubString (curr_alignment->seq , curr_alignment->soft_clippings.left , 0 , left_soft_clip_point - 1);
			extractSubString (curr_alignment->qual , curr_alignment->soft_clippings.left_qual , 0 , left_soft_clip_point - 1);
			curr_alignment->soft_clips_removed_seq_len = curr_alignment->read_seq_len - curr_alignment->cigar_items[0].len;
			for ( i = 0 ; i <= left_soft_clip_point - 1 ; i++ )
				curr_alignment->soft_clippings.left[i] =
						( curr_alignment->soft_clippings.left[i] >= 65 && curr_alignment->soft_clippings.left[i] <= 90 ) ? curr_alignment->soft_clippings.left[i] + 32 : curr_alignment->soft_clippings.left[i];
		}
		if ( curr_alignment->cigar_items[curr_alignment->number_of_cigar_items - 1].def == 'S' )
		{
			right_soft_clip_point = curr_alignment->cigar_items[curr_alignment->number_of_cigar_items - 1].len;
			extractSubString (curr_alignment->seq , curr_alignment->soft_clippings.right , curr_alignment->read_seq_len - right_soft_clip_point , curr_alignment->read_seq_len - 1);
			extractSubString (curr_alignment->qual , curr_alignment->soft_clippings.right_qual , curr_alignment->read_seq_len - right_soft_clip_point , curr_alignment->read_seq_len - 1);
			curr_alignment->soft_clips_removed_seq_len = curr_alignment->read_seq_len - curr_alignment->cigar_items[curr_alignment->number_of_cigar_items - 1].len;
			for ( i = 0 ; i <= right_soft_clip_point ; i++ )
				curr_alignment->soft_clippings.right[i] =
						( curr_alignment->soft_clippings.right[i] >= 65 && curr_alignment->soft_clippings.right[i] <= 90 ) ? curr_alignment->soft_clippings.right[i] + 32 : curr_alignment->soft_clippings.right[i];
		}

		/*
		 * Remove the soft-clips and construct new sequence
		 */

		int j = 0;
		for ( i = left_soft_clip_point ;
				i < curr_alignment->read_seq_len - right_soft_clip_point ; i++ )
		{
			curr_alignment->soft_clips_removed_seq[j] = curr_alignment->seq[i];
			curr_alignment->soft_clips_removed_qual[j] = curr_alignment->qual[i];
			j++;
		}
		curr_alignment->soft_clips_removed_seq[j] = '\0';
		curr_alignment->soft_clips_removed_qual[j] = '\0';
	}
	else
	{
		strcpy (curr_alignment->soft_clips_removed_seq , curr_alignment->seq);
		strcpy (curr_alignment->soft_clips_removed_qual , curr_alignment->qual);
		curr_alignment->soft_clips_removed_seq_len = strlen (curr_alignment->soft_clips_removed_seq);
	}

	/*
	 * Find the different tags
	 */
	XS_tag_index = -1;
	NH_tag_index = -1;
	NM_tag_index = -1;
	nM_tag_index = -1;
	MD_tag_index = -1;
	for ( i = 0 ; i < curr_alignment->number_of_tag_items ; i++ )
	{
		if ( strcmp (curr_alignment->tags[i].name , "MD") == 0 )
			MD_tag_index = i;
		if ( strcmp (curr_alignment->tags[i].name , "XS") == 0 )
			XS_tag_index = i;
		if ( strcmp (curr_alignment->tags[i].name , "NH") == 0 )
			NH_tag_index = i;
		if ( strcmp (curr_alignment->tags[i].name , "NM") == 0 )
			NM_tag_index = i;
		if ( strcmp (curr_alignment->tags[i].name , "nM") == 0 )
			nM_tag_index = i;
	}

	perfect_alignment_indicator = isAlignmentPerfect (curr_alignment->cigar , curr_alignment->tags , MD_tag_index , NM_tag_index , nM_tag_index);

	if ( perfect_alignment_indicator == 0 )
		designIntegratedCIGARSingleEnded (curr_alignment->tags[MD_tag_index].val , curr_alignment->soft_clips_removed_seq , &curr_alignment->soft_clips_removed_seq_len , curr_alignment->soft_clips_removed_qual , curr_alignment->cigar_items , curr_alignment->number_of_cigar_items , curr_alignment->cigar , curr_alignment->cigar_extended , curr_alignment->md_extended , curr_alignment->icigar , curr_alignment->splices , flag_ignore_quality_score , flag_ignore_mismatches);
	else
	{
		if ( curr_alignment->cigar_items[0].def == 'S' && curr_alignment->cigar_items[curr_alignment->number_of_cigar_items - 1].def == 'S' )
		{
			strcpy (curr_alignment->icigar , "");
			for ( i = 1 ; i < curr_alignment->number_of_cigar_items - 1 ; i++ )
			{
				sprintf (str , "%d" , curr_alignment->cigar_items[i].len);
				strcat (curr_alignment->icigar , str);
				temp_str[0] = curr_alignment->cigar_items[i].def;
				temp_str[1] = '\0';
				strcat (curr_alignment->icigar , temp_str);
			}
		}
		else if ( curr_alignment->cigar_items[0].def != 'S' && curr_alignment->cigar_items[curr_alignment->number_of_cigar_items - 1].def == 'S' )
		{
			strcpy (curr_alignment->icigar , "");
			for ( i = 0 ; i < curr_alignment->number_of_cigar_items - 1 ; i++ )
			{
				sprintf (str , "%d" , curr_alignment->cigar_items[i].len);
				strcat (curr_alignment->icigar , str);
				temp_str[0] = curr_alignment->cigar_items[i].def;
				temp_str[1] = '\0';
				strcat (curr_alignment->icigar , temp_str);
			}
		}
		else if ( curr_alignment->cigar_items[0].def == 'S' && curr_alignment->cigar_items[curr_alignment->number_of_cigar_items - 1].def != 'S' )
		{
			strcpy (curr_alignment->icigar , "");
			for ( i = 1 ; i < curr_alignment->number_of_cigar_items ; i++ )
			{
				sprintf (str , "%d" , curr_alignment->cigar_items[i].len);
				strcat (curr_alignment->icigar , str);
				temp_str[0] = curr_alignment->cigar_items[i].def;
				temp_str[1] = '\0';
				strcat (curr_alignment->icigar , temp_str);
			}
		}
		else strcpy (curr_alignment->icigar , curr_alignment->cigar);
		//printf("\nPartial icigar %s left soft clip %d right soft clip %d", curr_alignment->icigar, left_soft_clip_point, right_soft_clip_point);
	}

	/*
	 * Prepend and append the soft clips to icigar
	 */
	curr_alignment->temp[0] = '\0';
	if ( flag_ignore_soft_clippings == 0 ) // DO NOT ignore soft clippings
	{
		if ( left_soft_clip_point != 0 && right_soft_clip_point == 0 )
		{
			if ( flag_ignore_quality_score == 0 )
			{
				int j = 0;
				for ( i = 0 ;
						curr_alignment->soft_clippings.left_qual[i] != '\0' ;
						i++ )
				{
					curr_alignment->temp[j++ ] = curr_alignment->soft_clippings.left[i];
					curr_alignment->temp[j++ ] = curr_alignment->soft_clippings.left_qual[i];
				}
				curr_alignment->temp[j] = '\0';
			}
			else strcpy (curr_alignment->temp , curr_alignment->soft_clippings.left);
			strcat (curr_alignment->temp , curr_alignment->icigar);
		}
		else if ( left_soft_clip_point == 0 && right_soft_clip_point != 0 )
		{
			strcpy (curr_alignment->temp , curr_alignment->icigar);
			if ( flag_ignore_quality_score == 1 )
				strcat (curr_alignment->temp , curr_alignment->soft_clippings.right);
			else
			{
				int j = strlen (curr_alignment->temp);
				for ( i = 0 ;
						curr_alignment->soft_clippings.right_qual[i] != '\0' ;
						i++ )
				{
					curr_alignment->temp[j++ ] = curr_alignment->soft_clippings.right[i];
					curr_alignment->temp[j++ ] = curr_alignment->soft_clippings.right_qual[i];
				}
				curr_alignment->temp[j] = '\0';
			}
		}
		else if ( left_soft_clip_point != 0 && right_soft_clip_point != 0 )
		{
			if ( flag_ignore_quality_score == 0 )
			{
				int j = 0;
				for ( i = 0 ;
						curr_alignment->soft_clippings.left_qual[i] != '\0' ;
						i++ )
				{
					curr_alignment->temp[j++ ] = curr_alignment->soft_clippings.left[i];
					curr_alignment->temp[j++ ] = curr_alignment->soft_clippings.left_qual[i];
				}
				curr_alignment->temp[j] = '\0';
			}
			else strcpy (curr_alignment->temp , curr_alignment->soft_clippings.left);
			strcat (curr_alignment->temp , curr_alignment->icigar);
			if ( flag_ignore_quality_score == 1 )
				strcat (curr_alignment->temp , curr_alignment->soft_clippings.right);
			else
			{
				int j = strlen (curr_alignment->temp);
				for ( i = 0 ;
						curr_alignment->soft_clippings.right_qual[i] != '\0' ;
						i++ )
				{
					curr_alignment->temp[j++ ] = curr_alignment->soft_clippings.right[i];
					curr_alignment->temp[j++ ] = curr_alignment->soft_clippings.right_qual[i];
				}
				curr_alignment->temp[j] = '\0';
			}
		}
		else if ( left_soft_clip_point == 0 && right_soft_clip_point == 0 )
		{
			strcpy (curr_alignment->temp , curr_alignment->icigar);
		}
	}
	else // DO ignore soft clippings
	{
		strcpy (curr_alignment->temp , curr_alignment->icigar);
	}

	strcpy (curr_alignment->icigar , curr_alignment->temp);

	/*
	 * Add NH tag
	 */
	/*if (XS_tag_index != -1) strcat(curr_alignment->icigar, curr_alignment->tags[XS_tag_index].val);*/
	if ( NH_tag_index != -1 )
		strcat (curr_alignment->icigar , curr_alignment->tags[NH_tag_index].val);

	/*
	 * Change the iCIGAR representation to reflect the samformatflag and XS tag
	 */
	if ( spliced_alignment_indicator == 0 )
	{
		switch ( curr_alignment->samflag )
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
		for ( i = 0 ; i < strlen (curr_alignment->icigar) ; i++ )
			if ( curr_alignment->icigar[i] == 'M' )
				curr_alignment->icigar[i] = M_replacement_character;
	}
	else
	{
		if ( XS_tag_index == -1 )
		{
			switch ( curr_alignment->samflag )
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
			if ( strcmp (curr_alignment->tags[XS_tag_index].val , "+") == 0 )
			{
				switch ( curr_alignment->samflag )
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
			else if ( strcmp (curr_alignment->tags[XS_tag_index].val , "-") == 0 )
			{
				switch ( curr_alignment->samflag )
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
		for ( i = 0 ; i < strlen (curr_alignment->icigar) ; i++ )
			if ( curr_alignment->icigar[i] == 'M' )
				curr_alignment->icigar[i] = M_replacement_character;
	}

	/*
	 * For diagnostics
	 */

	if ( run_diagnostics == 1 )
	{
		dummy = 'Z';
		strcpy (sam_alignment_instance_diagnostics->icigar , curr_alignment->icigar);
		sam_alignment_instance_diagnostics->start_position = curr_alignment->start_position;
		convertIcigarToCigarandMDSingleEnded (whole_genome , sam_alignment_instance_diagnostics , curr_alignment->reference_name , flag_ignore_mismatches , flag_ignore_soft_clippings , flag_ignore_unmapped_sequences , flag_ignore_quality_score , 0 , &dummy);
		/************************************************************************
		 * Remove this section later
		 *************************************************************************/
		/*if (isCharacterInString(sam_alignment_instance_diagnostics->cigar, 'I') && isCharacterInString(sam_alignment_instance_diagnostics->cigar, 'D') && isCharacterInString(sam_alignment_instance_diagnostics->cigar, 'N') && isCharacterInString(sam_alignment_instance_diagnostics->cigar, 'S'))
		 {
		 printSamAlignmentInstance(curr_alignment, 1);
		 printSamAlignmentInstance(sam_alignment_instance_diagnostics, 1);
		 exit(1);
		 }
		 */
		/************************************************************************/

		flag = 0;
		if ( strcmp (sam_alignment_instance_diagnostics->cigar , curr_alignment->cigar) != 0 )
		{
			printf ("\nCIGAR Mismatch");
			flag = 1;
		}
		if ( strcmp (sam_alignment_instance_diagnostics->seq , curr_alignment->seq) != 0 )
		{
			printf ("\nSEQ Mismatch");
			flag = 1;
		}
		if ( strcmp (sam_alignment_instance_diagnostics->tags[2].val , curr_alignment->tags[MD_tag_index].val) != 0 )
		{
			printf ("\nMD Mismatch Actual: %s Deciphered: %s" , curr_alignment->tags[MD_tag_index].val , sam_alignment_instance_diagnostics->tags[2].val);
			flag = 1;
		}

		if ( flag )
		{
			printf ("\nRecord Number : %lld" , number_of_records_read);
			printSamAlignmentInstance (curr_alignment , 1);
			printSamAlignmentInstance (sam_alignment_instance_diagnostics , 1);
			exit (1);
		}
	}
}

char findReplamentCharacterForPairedEndedReads (int samflag, struct Paired_Ended_Flag_to_Single_Character *samflag_dictionary, int number_of_unique_samformatflags, int XS_tag_index, struct Sam_Tags *tags)
{
	int i, j, k;
	char replacement_character;
	for ( i = 0 ; i < number_of_unique_samformatflags * 2 ; i++ )
	{
		if ( XS_tag_index == -1 ) // Either splice does not exist or direction could not be determined from the alignment
		{
			if ( samflag_dictionary->direction[i] == '+' && samflag == samflag_dictionary->samflags[i] )
			{
				replacement_character = samflag_dictionary->character[i];
				break;
			}
		}
		else
		{
			if ( strcmp (tags[XS_tag_index].val , "+") == 0 && samflag_dictionary->direction[i] == '+' && samflag == samflag_dictionary->samflags[i] )
			{
				replacement_character = samflag_dictionary->character[i];
				break;
			}
			else if ( strcmp (tags[XS_tag_index].val , "-") == 0 && samflag_dictionary->direction[i] == '-' && samflag == samflag_dictionary->samflags[i] )
			{
				replacement_character = samflag_dictionary->character[i];
				break;
			}
		}
	}

	return replacement_character;

}

void generateIntegratedCigarPairedEnded (struct Sam_Alignment *curr_alignment, short int flag_ignore_soft_clippings, short int flag_ignore_mismatches, short int flag_ignore_unmapped_sequences, short int flag_ignore_quality_score, struct Whole_Genome_Sequence *whole_genome, struct Sam_Alignment *sam_alignment_instance_diagnostics, long long int number_of_records_read, short int run_diagnostics, struct Paired_Ended_Flag_to_Single_Character *samflag_dictionary, int number_of_unique_samformatflags, char samformatflag_replacer_characters[])
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
	int flag;
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
	char dummy;
	/********************************************************************/

	splitCigar (curr_alignment->cigar , &curr_alignment->number_of_cigar_items , curr_alignment->cigar_items);

	for ( i = 0 ; i < curr_alignment->number_of_cigar_items ; i++ )
	{
		if ( curr_alignment->cigar_items[i].def == 'N' )
		{
			spliced_alignment_indicator = 1;
			break;
		}
	}
	/*
	 * Process each alignment to extract soft clipped portion of reads
	 */
	if ( isSequenceSoftClipped (curr_alignment->cigar) == 1 )
	{
		if ( curr_alignment->cigar_items[0].def == 'S' )
		{
			left_soft_clip_point = curr_alignment->cigar_items[0].len;
			extractSubString (curr_alignment->seq , curr_alignment->soft_clippings.left , 0 , left_soft_clip_point - 1);
			extractSubString (curr_alignment->qual , curr_alignment->soft_clippings.left_qual , 0 , left_soft_clip_point - 1);
			curr_alignment->soft_clips_removed_seq_len = curr_alignment->read_seq_len - curr_alignment->cigar_items[0].len;
			for ( i = 0 ; i <= left_soft_clip_point - 1 ; i++ )
				curr_alignment->soft_clippings.left[i] =
						( curr_alignment->soft_clippings.left[i] >= 65 && curr_alignment->soft_clippings.left[i] <= 90 ) ? curr_alignment->soft_clippings.left[i] + 32 : curr_alignment->soft_clippings.left[i];
		}
		if ( curr_alignment->cigar_items[curr_alignment->number_of_cigar_items - 1].def == 'S' )
		{
			right_soft_clip_point = curr_alignment->cigar_items[curr_alignment->number_of_cigar_items - 1].len;
			extractSubString (curr_alignment->seq , curr_alignment->soft_clippings.right , curr_alignment->read_seq_len - right_soft_clip_point , curr_alignment->read_seq_len - 1);
			extractSubString (curr_alignment->qual , curr_alignment->soft_clippings.right_qual , curr_alignment->read_seq_len - right_soft_clip_point , curr_alignment->read_seq_len - 1);
			curr_alignment->soft_clips_removed_seq_len = curr_alignment->read_seq_len - curr_alignment->cigar_items[curr_alignment->number_of_cigar_items - 1].len;
			for ( i = 0 ; i <= right_soft_clip_point ; i++ )
				curr_alignment->soft_clippings.right[i] =
						( curr_alignment->soft_clippings.right[i] >= 65 && curr_alignment->soft_clippings.right[i] <= 90 ) ? curr_alignment->soft_clippings.right[i] + 32 : curr_alignment->soft_clippings.right[i];
		}

		/*
		 * Remove the soft-clips and construct new sequence
		 */

		int j = 0;
		for ( i = left_soft_clip_point ;
				i < curr_alignment->read_seq_len - right_soft_clip_point ; i++ )
		{
			curr_alignment->soft_clips_removed_seq[j] = curr_alignment->seq[i];
			curr_alignment->soft_clips_removed_qual[j] = curr_alignment->qual[i];
			j++;
		}
		curr_alignment->soft_clips_removed_seq[j] = '\0';
		curr_alignment->soft_clips_removed_qual[j] = '\0';
	}
	else
	{
		strcpy (curr_alignment->soft_clips_removed_seq , curr_alignment->seq);
		strcpy (curr_alignment->soft_clips_removed_qual , curr_alignment->qual);
		curr_alignment->soft_clips_removed_seq_len = strlen (curr_alignment->soft_clips_removed_seq);
	}

	/*
	 * Find the different tags
	 */
	XS_tag_index = -1;
	NH_tag_index = -1;
	NM_tag_index = -1;
	nM_tag_index = -1;
	MD_tag_index = -1;
	for ( i = 0 ; i < curr_alignment->number_of_tag_items ; i++ )
	{
		if ( strcmp (curr_alignment->tags[i].name , "MD") == 0 )
			MD_tag_index = i;
		if ( strcmp (curr_alignment->tags[i].name , "XS") == 0 )
			XS_tag_index = i;
		if ( strcmp (curr_alignment->tags[i].name , "NH") == 0 )
			NH_tag_index = i;
		if ( strcmp (curr_alignment->tags[i].name , "NM") == 0 )
			NM_tag_index = i;
		if ( strcmp (curr_alignment->tags[i].name , "nM") == 0 )
			nM_tag_index = i;
	}

	perfect_alignment_indicator = isAlignmentPerfect (curr_alignment->cigar , curr_alignment->tags , MD_tag_index , NM_tag_index , nM_tag_index);

	if ( perfect_alignment_indicator == 0 )
		designIntegratedCIGARSingleEnded (curr_alignment->tags[MD_tag_index].val , curr_alignment->soft_clips_removed_seq , &curr_alignment->soft_clips_removed_seq_len , curr_alignment->soft_clips_removed_qual , curr_alignment->cigar_items , curr_alignment->number_of_cigar_items , curr_alignment->cigar , curr_alignment->cigar_extended , curr_alignment->md_extended , curr_alignment->icigar , curr_alignment->splices , flag_ignore_quality_score , flag_ignore_mismatches);
	else
	{
		if ( curr_alignment->cigar_items[0].def == 'S' && curr_alignment->cigar_items[curr_alignment->number_of_cigar_items - 1].def == 'S' )
		{
			strcpy (curr_alignment->icigar , "");
			for ( i = 1 ; i < curr_alignment->number_of_cigar_items - 1 ; i++ )
			{
				sprintf (str , "%d" , curr_alignment->cigar_items[i].len);
				strcat (curr_alignment->icigar , str);
				temp_str[0] = curr_alignment->cigar_items[i].def;
				temp_str[1] = '\0';
				strcat (curr_alignment->icigar , temp_str);
			}
		}
		else if ( curr_alignment->cigar_items[0].def != 'S' && curr_alignment->cigar_items[curr_alignment->number_of_cigar_items - 1].def == 'S' )
		{
			strcpy (curr_alignment->icigar , "");
			for ( i = 0 ; i < curr_alignment->number_of_cigar_items - 1 ; i++ )
			{
				sprintf (str , "%d" , curr_alignment->cigar_items[i].len);
				strcat (curr_alignment->icigar , str);
				temp_str[0] = curr_alignment->cigar_items[i].def;
				temp_str[1] = '\0';
				strcat (curr_alignment->icigar , temp_str);
			}
		}
		else if ( curr_alignment->cigar_items[0].def == 'S' && curr_alignment->cigar_items[curr_alignment->number_of_cigar_items - 1].def != 'S' )
		{
			strcpy (curr_alignment->icigar , "");
			for ( i = 1 ; i < curr_alignment->number_of_cigar_items ; i++ )
			{
				sprintf (str , "%d" , curr_alignment->cigar_items[i].len);
				strcat (curr_alignment->icigar , str);
				temp_str[0] = curr_alignment->cigar_items[i].def;
				temp_str[1] = '\0';
				strcat (curr_alignment->icigar , temp_str);
			}
		}
		else strcpy (curr_alignment->icigar , curr_alignment->cigar);
		//printf("\nPartial icigar %s left soft clip %d right soft clip %d", curr_alignment->icigar, left_soft_clip_point, right_soft_clip_point);
	}

	/*
	 * Prepend and append the soft clips to icigar
	 */
	curr_alignment->temp[0] = '\0';
	if ( flag_ignore_soft_clippings == 0 ) // DO NOT ignore soft clippings
	{
		if ( left_soft_clip_point != 0 && right_soft_clip_point == 0 )
		{
			if ( flag_ignore_quality_score == 0 )
			{
				int j = 0;
				for ( i = 0 ;
						curr_alignment->soft_clippings.left_qual[i] != '\0' ;
						i++ )
				{
					curr_alignment->temp[j++ ] = curr_alignment->soft_clippings.left[i];
					curr_alignment->temp[j++ ] = curr_alignment->soft_clippings.left_qual[i];
				}
				curr_alignment->temp[j] = '\0';
			}
			else strcpy (curr_alignment->temp , curr_alignment->soft_clippings.left);
			strcat (curr_alignment->temp , curr_alignment->icigar);
		}
		else if ( left_soft_clip_point == 0 && right_soft_clip_point != 0 )
		{
			strcpy (curr_alignment->temp , curr_alignment->icigar);
			if ( flag_ignore_quality_score == 1 )
				strcat (curr_alignment->temp , curr_alignment->soft_clippings.right);
			else
			{
				int j = strlen (curr_alignment->temp);
				for ( i = 0 ;
						curr_alignment->soft_clippings.right_qual[i] != '\0' ;
						i++ )
				{
					curr_alignment->temp[j++ ] = curr_alignment->soft_clippings.right[i];
					curr_alignment->temp[j++ ] = curr_alignment->soft_clippings.right_qual[i];
				}
				curr_alignment->temp[j] = '\0';
			}
		}
		else if ( left_soft_clip_point != 0 && right_soft_clip_point != 0 )
		{
			if ( flag_ignore_quality_score == 0 )
			{
				int j = 0;
				for ( i = 0 ;
						curr_alignment->soft_clippings.left_qual[i] != '\0' ;
						i++ )
				{
					curr_alignment->temp[j++ ] = curr_alignment->soft_clippings.left[i];
					curr_alignment->temp[j++ ] = curr_alignment->soft_clippings.left_qual[i];
				}
				curr_alignment->temp[j] = '\0';
			}
			else strcpy (curr_alignment->temp , curr_alignment->soft_clippings.left);
			strcat (curr_alignment->temp , curr_alignment->icigar);
			if ( flag_ignore_quality_score == 1 )
				strcat (curr_alignment->temp , curr_alignment->soft_clippings.right);
			else
			{
				int j = strlen (curr_alignment->temp);
				for ( i = 0 ;
						curr_alignment->soft_clippings.right_qual[i] != '\0' ;
						i++ )
				{
					curr_alignment->temp[j++ ] = curr_alignment->soft_clippings.right[i];
					curr_alignment->temp[j++ ] = curr_alignment->soft_clippings.right_qual[i];
				}
				curr_alignment->temp[j] = '\0';
			}
		}
		else if ( left_soft_clip_point == 0 && right_soft_clip_point == 0 )
		{
			strcpy (curr_alignment->temp , curr_alignment->icigar);
		}
	}
	else // DO ignore soft clippings
	{
		strcpy (curr_alignment->temp , curr_alignment->icigar);
	}

	strcpy (curr_alignment->icigar , curr_alignment->temp);

	/*
	 * Add NH tag
	 */
	/*if (XS_tag_index != -1) strcat(curr_alignment->icigar, curr_alignment->tags[XS_tag_index].val);*/
	if ( NH_tag_index != -1 )
		strcat (curr_alignment->icigar , curr_alignment->tags[NH_tag_index].val);

	/*
	 * Change the iCIGAR representation to reflect the samformatflag and XS tag
	 */
	/*if ( spliced_alignment_indicator == 0 )
	 {
	 M_replacement_character = findReplamentCharacterForPairedEndedReads (curr_alignment->samflag , samflag_dictionary , number_of_unique_samformatflags , XS_tag_index , curr_alignment->tags);
	 for ( i = 0 ; i < strlen (curr_alignment->icigar) ; i++ )
	 if ( curr_alignment->icigar[i] == 'M' )
	 curr_alignment->icigar[i] = M_replacement_character;
	 }
	 else
	 {
	 M_replacement_character = findReplamentCharacterForPairedEndedReads (curr_alignment->samflag , samflag_dictionary , number_of_unique_samformatflags , XS_tag_index , curr_alignment->tags);
	 for ( i = 0 ; i < strlen (curr_alignment->icigar) ; i++ )
	 if ( curr_alignment->icigar[i] == 'M' )
	 curr_alignment->icigar[i] = M_replacement_character;
	 }
	 */
	M_replacement_character = findReplamentCharacterForPairedEndedReads (curr_alignment->samflag , samflag_dictionary , number_of_unique_samformatflags , XS_tag_index , curr_alignment->tags);
	for ( i = 0 ; i < strlen (curr_alignment->icigar) ; i++ )
		if ( curr_alignment->icigar[i] == 'M' )
			curr_alignment->icigar[i] = M_replacement_character;
	/*
	 printf ("\nCIGAR %s iCIGAR: %s" , curr_alignment->cigar , curr_alignment->icigar);
	 fflush (stdout);
	 */

	/*
	 * For diagnostics
	 */

	if ( run_diagnostics == 1 )
	{
		dummy = 'Z';
		strcpy (sam_alignment_instance_diagnostics->icigar , curr_alignment->icigar);
		sam_alignment_instance_diagnostics->start_position = curr_alignment->start_position;
		convertIcigarToCigarandMDPairedEnded (whole_genome , sam_alignment_instance_diagnostics , curr_alignment->reference_name , flag_ignore_mismatches , flag_ignore_soft_clippings , flag_ignore_unmapped_sequences , flag_ignore_quality_score , 0 , &dummy , samflag_dictionary , number_of_unique_samformatflags , samformatflag_replacer_characters);
		/************************************************************************
		 * Remove this section later
		 *************************************************************************/
		/*if (isCharacterInString(sam_alignment_instance_diagnostics->cigar, 'I') && isCharacterInString(sam_alignment_instance_diagnostics->cigar, 'D') && isCharacterInString(sam_alignment_instance_diagnostics->cigar, 'N') && isCharacterInString(sam_alignment_instance_diagnostics->cigar, 'S'))
		 {
		 printSamAlignmentInstance(curr_alignment, 1);
		 printSamAlignmentInstance(sam_alignment_instance_diagnostics, 1);
		 exit(1);
		 }
		 */
		/************************************************************************/

		flag = 0;
		if ( strcmp (sam_alignment_instance_diagnostics->cigar , curr_alignment->cigar) != 0 )
		{
			printf ("\nCIGAR Mismatch");
			flag = 1;
		}
		if ( strcmp (sam_alignment_instance_diagnostics->seq , curr_alignment->seq) != 0 )
		{
			printf ("\nSEQ Mismatch");
			flag = 1;
		}
		if ( strcmp (sam_alignment_instance_diagnostics->tags[2].val , curr_alignment->tags[MD_tag_index].val) != 0 )
		{
			printf ("\nMD Mismatch Actual: %s Deciphered: %s" , curr_alignment->tags[MD_tag_index].val , sam_alignment_instance_diagnostics->tags[2].val);
			flag = 1;
		}

		if ( flag )
		{
			printf ("\nRecord Number : %lld" , number_of_records_read);
			printSamAlignmentInstance (curr_alignment , 1);
			printSamAlignmentInstance (sam_alignment_instance_diagnostics , 1);
			exit (1);
		}
	}
}

void initializePass3_Compression_Symbol_icigar_MappingPool (struct Pass3_Compression_Symbol_icigar_Mapping **symbol_icigar_mapping)
{
	int i;
	int j;
	int symbol_icigar_mapping_index = 0;
	char *individual_symbols = "+-0123456789abcdefghijklmnopqrstuvwxyzUVWXYZ";
	for ( i = 0 ; i < MAX_SYMBOLS_FOR_PASS3_COMPRESSION ; i++ )
	{
		for ( j = 0 ; j < MAX_SYMBOLS_FOR_PASS3_COMPRESSION ; j++ )
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

void assignicigarsToSymbols (struct Cigar_Frequency **cigar_freq_pool, int cigar_freq_pool_index, struct Pass3_Compression_Symbol_icigar_Mapping **symbol_icigar_mapping)
{
	int i;
	for ( i = 0 ; i < cigar_freq_pool_index ; i++ )
		strcpy (symbol_icigar_mapping[i]->icigar , cigar_freq_pool[i]->cigar);
	/*
	 for (i = 0; i < cigar_freq_pool_index; i++)
	 printf("\n icigar: %s symbol: %s", symbol_icigar_mapping[i]->icigar, symbol_icigar_mapping[i]->symbolic_icigar);
	 */
}

void readInEachChromosome (char *genome_filename, struct Whole_Genome_Sequence *whole_genome, char *chromosome)
{
	FILE *fhr;
	char *buffer = NULL;
	size_t len = 0;
	ssize_t line_len;
	int i;
	int j;

//buffer = (char*) malloc(sizeof(char) * pow(2, 32));
	fhr = fopen (genome_filename , "r");
	if ( fhr == NULL )
	{
		printf ("Error! File not found");
		exit (1);
	}

//printf ("\nwhole_genome->number_of_reference_sequences %d" , whole_genome->number_of_reference_sequences);
	if ( whole_genome->number_of_reference_sequences == 1 )
	{
		int total = strlen (whole_genome->reference_sequence_name[whole_genome->number_of_reference_sequences - 1]) + strlen (whole_genome->nucleotides[whole_genome->number_of_reference_sequences - 1]);
		printf ("\nFreeing %f KBytes" , ( float ) total / ( float ) 1024);
		free (whole_genome->reference_sequence_name[whole_genome->number_of_reference_sequences - 1]);
		free (whole_genome->nucleotides[whole_genome->number_of_reference_sequences - 1]);
	}
	whole_genome->number_of_reference_sequences = 0;
//printf ("\n Loading chromosome %s" , chromosome);
//fflush (stdout);
	while ( ( line_len = getline ( &buffer , &len , fhr) ) != -1 )
	{
		//printf(“\n%lld”, strlen(buffer));
		if ( strlen (buffer) <= 1 ) continue;
		if ( buffer[0] == '>' )
		{
			for ( i = 1 ; buffer[i] != 32 ; i++ )
				buffer[i - 1] = buffer[i];
			buffer[i - 1] = '\0';
			if ( strcmp (buffer , chromosome) == 0 )
			{
				/*
				 printf ("\n Found chromosome %s" , chromosome);
				 fflush (stdout);
				 */
				whole_genome->reference_sequence_name[whole_genome->number_of_reference_sequences] = ( char* ) malloc (sizeof(char) * ( line_len + 1 ));
				strcpy (whole_genome->reference_sequence_name[whole_genome->number_of_reference_sequences] , buffer);
				line_len = getline ( &buffer , &len , fhr);

				whole_genome->nucleotides[whole_genome->number_of_reference_sequences] = ( char* ) malloc (sizeof(char) * ( line_len + 1 ));
				for ( i = 0 ; buffer[i] != '\0' ; i++ )
					if ( buffer[i] >= 90 && buffer[i] <= 122 ) buffer[i] -= 32;
				strcpy (whole_genome->nucleotides[whole_genome->number_of_reference_sequences] , buffer);
				whole_genome->reference_sequence_length[whole_genome->number_of_reference_sequences] = strlen (buffer);
				whole_genome->number_of_reference_sequences = 1;
				break;
			}
		}

	}
	fclose (fhr);
	/*
	 for ( i = 0 ; i < whole_genome->number_of_reference_sequences ; i++ )
	 printf ("\nChromosome in whole genome %s" , whole_genome->reference_sequence_name[i]);
	 */
}

void readInTheEntireGenome (char *genome_filename, struct Whole_Genome_Sequence *whole_genome)
{
	FILE *fhr;
	char *buffer = NULL;
	size_t len = 0;
	ssize_t line_len;
	int i;
	int j;

//buffer = (char*) malloc(sizeof(char) * pow(2, 32));
	fhr = fopen (genome_filename , "rb");
	if ( fhr == NULL )
	{
		printf ("Error! File not found");
		exit (1);
	}

	whole_genome->nucleotides = ( char** ) malloc (sizeof(char*) * MAX_REFERENCE_SEQUENCES);
	whole_genome->reference_sequence_name = ( char** ) malloc (sizeof(char*) * MAX_REFERENCE_SEQUENCES);
	whole_genome->reference_sequence_length = ( unsigned long long int* ) malloc (sizeof(unsigned long long int) * MAX_REFERENCE_SEQUENCES);
	whole_genome->number_of_reference_sequences = 0;

	while ( ( line_len = getline ( &buffer , &len , fhr) ) != -1 )
	{
		//printf("\n%lld", strlen(buffer));
		if ( strlen (buffer) <= 1 ) continue;
		if ( buffer[0] == '>' )
		{
			whole_genome->reference_sequence_name[whole_genome->number_of_reference_sequences] = ( char* ) malloc (sizeof(char) * ( line_len + 1 ));
			j = 0;
			for ( i = 1 ; buffer[i] != 32 ; i++ )
				whole_genome->reference_sequence_name[whole_genome->number_of_reference_sequences][j++ ] = buffer[i];
			whole_genome->reference_sequence_name[whole_genome->number_of_reference_sequences][j++ ] = '\0';
		}
		else
		{
			whole_genome->nucleotides[whole_genome->number_of_reference_sequences] = ( char* ) malloc (sizeof(char) * ( line_len + 1 ));
			for ( i = 0 ; buffer[i] != '\0' ; i++ )
				if ( buffer[i] >= 90 && buffer[i] <= 122 ) buffer[i] -= 32;
			strcpy (whole_genome->nucleotides[whole_genome->number_of_reference_sequences] , buffer);
			whole_genome->reference_sequence_length[whole_genome->number_of_reference_sequences] = strlen (buffer);
			whole_genome->number_of_reference_sequences++;
		}

	}
//free(buffer);
//buffer = NULL;
	/*
	 for (i = 0; i < whole_genome->number_of_reference_sequences; i++)
	 printf("\n %d %s %lld %lld", i + 1, whole_genome->reference_sequence_name[i], strlen(whole_genome->nucleotides[i]), whole_genome->reference_sequence_length[i]);
	 */
	fclose (fhr);
}

void reverseComplement (char *seq, char *rev_seq)
{
	int i;
	int j;

	j = 0;
	for ( i = strlen (seq) - 1 ; i >= 0 ; i-- )
	{
		switch ( seq[i] )
		{
			case 'A':
				rev_seq[j++ ] = 'T';
				break;
			case 'T':
				rev_seq[j++ ] = 'A';
				break;
			case 'G':
				rev_seq[j++ ] = 'C';
				break;
			case 'C':
				rev_seq[j++ ] = 'G';
				break;
			case 'N':
				rev_seq[j++ ] = 'N';
				break;
		}
	}
	rev_seq[j++ ] = '\0';
}

void readInGenomeSequenceSingleChromosome (struct Whole_Genome_Sequence *single_genome_sequence, char *chromosome, char *genome_filename, struct Abridge_Index *genome_index)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	FILE *fhr;

	char *buffer;

	size_t len = 0;
	ssize_t line_len;

	int i;
	int j;
	int fseek_ret_val;
	int fread_ret_val;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr = fopen (genome_filename , "rb");
	if ( fhr == NULL )
	{
		printf ("Error! File not found");
		exit (1);
	}
	/********************************************************************/

	buffer = ( char* ) malloc (sizeof(char) * MAX_REFERENCE_SEQ_LEN);
	single_genome_sequence->nucleotides = ( char** ) malloc (sizeof(char*) * MAX_REFERENCE_SEQUENCES);
	single_genome_sequence->reference_sequence_name = ( char** ) malloc (sizeof(char*) * MAX_REFERENCE_SEQUENCES);
	single_genome_sequence->reference_sequence_length = ( unsigned long long int* ) malloc (sizeof(unsigned long long int) * MAX_REFERENCE_SEQUENCES);
	single_genome_sequence->number_of_reference_sequences = 1;

	for ( i = 0 ; i < genome_index->number_of_items ; i++ )
		if ( strcmp (chromosome , genome_index->chromosome[i]) == 0 ) break;

	fseek_ret_val = fseek (fhr , genome_index->start_byte[i] , SEEK_SET);
	fread_ret_val = fread (buffer , 1 , genome_index->end_byte[i] - genome_index->start_byte[i] , fhr);
	if ( buffer == NULL ) printf ("\nBuffer is null chromosome number %d" , i);
//printf("\nReading from %lld Num bytes read %lld fseek_ret_val %lld fread_ret_val %lld ferror %lld feof %lld Genome filename %s\n", genome_index->start_byte[i], genome_index->end_byte[i] - genome_index->start_byte[i], fseek_ret_val, fread_ret_val, ferror(fhr), feof(fhr), genome_filename);
//fflush (stdout);
	buffer[genome_index->end_byte[i] - genome_index->start_byte[i] + 1] = '\0';
	single_genome_sequence->reference_sequence_name[0] = ( char* ) malloc (sizeof(char) * ( 100 ));
	strcpy (single_genome_sequence->reference_sequence_name[0] , chromosome);
	single_genome_sequence->nucleotides[0] = ( char* ) malloc (sizeof(char) * strlen (buffer));
	strcpy (single_genome_sequence->nucleotides[0] , buffer);
	single_genome_sequence->reference_sequence_length[0] = genome_index->end_byte[i] - genome_index->start_byte[i];
	fclose (fhr);
}

int findReadClusterFromAbridgeIndex (struct Abridge_Index *abridge_index, char *chromosome, long long int start, long long int end, long long int *abridge_match_start_index, long long int *abridge_match_end_index)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	int i;
	short int first_record = 0;
	long long int start1, start2;
	long long int end1, end2;
	bool condition1;
	bool condition2;
	bool condition3;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	*abridge_match_start_index = -1;
	*abridge_match_end_index = -1;
	start1 = start;
	end1 = end;
	/********************************************************************/
	for ( i = 0 ; i < abridge_index->number_of_items ; i++ )
	{
		if ( strcmp (chromosome , abridge_index->chromosome[i]) == 0 )
		{
			start2 = abridge_index->start[i];
			end2 = abridge_index->end[i];
			condition1 = start1 >= start2 && start1 <= end2 && end1 >= start1 && end1 <= end2;
			condition2 = start1 <= start2 && start1 <= end2 && end1 >= start1 && end1 <= end2;
			condition3 = start1 >= start2 && start1 <= end2 && end1 >= start2 && end1 >= end2;

			condition1 = end1 < start2;
			condition2 = end2 < start1;
			//printf("\n%d %d %d %d %d %d %d", start1, end1, start2, end2, condition1, condition2, !(condition1 || condition2));
			if ( ! ( condition1 || condition2 ) )
			{
				if ( *abridge_match_end_index == -1 )
				{
					*abridge_match_start_index = i;
					if ( i == 0 )
						first_record = 1;
					else
					{
						if ( strcmp (abridge_index->chromosome[i] , abridge_index->chromosome[i - 1]) != 0 )
							first_record = 1;
					}
				}
				*abridge_match_end_index = i;
			}
		}
	}
	return first_record;
}

void readGenomeIndex (struct Abridge_Index *genome_index, char *genome_index_filename, char **split_line)
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
	fhr = fopen (genome_index_filename , "r");
	if ( fhr == NULL )
	{
		printf ("Error! File not found");
		exit (1);
	}
	/********************************************************************/

	genome_index->number_of_items = 0;
	while ( ( line_len = getline ( &line , &len , fhr) ) != -1 )
	{
		line_num++;
		splitByDelimiter (line , '\t' , split_line);
		strcpy (genome_index->chromosome[genome_index->number_of_items] , split_line[0]);
		genome_index->start_byte[genome_index->number_of_items] = strtol (split_line[2] , &temp , 10);
		genome_index->end_byte[genome_index->number_of_items] = genome_index->start_byte[genome_index->number_of_items] + strtol (split_line[3] , &temp , 10);
		genome_index->number_of_items++;
	}
}

void readAbridgeIndex (struct Abridge_Index *abridge_index, char *abridge_index_filename, char **split_line, short int *flag_ignore_mismatches, short int *flag_ignore_soft_clippings, short int *flag_ignore_unmapped_sequences, short int *flag_ignore_quality_score, short int *flag_save_all_quality_scores, short int *flag_save_exact_quality_scores)
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

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr = fopen (abridge_index_filename , "r");
	if ( fhr == NULL )
	{
		printf ("Error! File not found");
		exit (1);
	}
	/********************************************************************/
	while ( ( line_len = getline ( &line , &len , fhr) ) != -1 )
	{
		line_num++;
		splitByDelimiter (line , '\t' , split_line);
		if ( line_num == 1 )
		{
			//printf ("\n First line %s" , line);
			*flag_ignore_mismatches = strtol (split_line[0] , &temp , 10);
			*flag_ignore_soft_clippings = strtol (split_line[1] , &temp , 10);
			*flag_ignore_unmapped_sequences = strtol (split_line[2] , &temp , 10);
			*flag_ignore_quality_score = strtol (split_line[3] , &temp , 10);
			*flag_save_all_quality_scores = strtol (split_line[4] , &temp , 10);
			*flag_save_exact_quality_scores = strtol (split_line[5] , &temp , 10);
			continue;
		}
		//printf("\n%d %s", abridge_index->number_of_items, split_line[0]);
		//fflush(stdout);
		strcpy (abridge_index->chromosome[abridge_index->number_of_items] , split_line[0]);
		abridge_index->start[abridge_index->number_of_items] = strtol (split_line[1] , &temp , 10);
		abridge_index->end[abridge_index->number_of_items] = strtol (split_line[2] , &temp , 10);
		abridge_index->start_byte[abridge_index->number_of_items] = strtol (split_line[3] , &temp , 10);
		abridge_index->end_byte[abridge_index->number_of_items] = strtol (split_line[4] , &temp , 10);
		abridge_index->start_byte_qual[abridge_index->number_of_items] = strtol (split_line[5] , &temp , 10);
		abridge_index->end_byte_qual[abridge_index->number_of_items] = strtol (split_line[6] , &temp , 10);
		abridge_index->number_of_items++;
		//printf("\n %lld", MAX_POOL_SIZE - abridge_index->number_of_items);
	}
	fclose (fhr);
}

int maxClusterSize (struct Abridge_Index *abridge_index)
{
	long long int i;
	long long int max_cluster_size = -1;
	for ( i = 0 ; i < abridge_index->number_of_items ; i++ )
		if ( max_cluster_size < ( abridge_index->end[i] - abridge_index->start[i] ) )
			max_cluster_size = ( abridge_index->end[i] - abridge_index->start[i] );
	return max_cluster_size;

}

void loadSequencesFromFile (char **sequences, char *filename)
{
	FILE *fhr;
	long long int num_sequence = 0;
	size_t len = 0;
	ssize_t line_len;
	char *line = NULL; // for reading each line

	fhr = fopen (filename , "rb");
	if ( fhr == NULL )
	{
		printf ("Error! File not found");
		exit (1);
	}
	while ( ( line_len = getline ( &line , &len , fhr) ) != -1 )
		if ( line[0] != '>' ) strcpy (sequences[num_sequence++ ] , line);

}

long long int reverseInteger (long long int val)
{
	long long int temp = 0;
	while ( val )
	{
		temp = temp * 10 + ( val % 10 );
		val /= 10;
	}
	return temp;
}

long long int extractNHfromicigar (char *icigar, int icigar_length)
{
	char nh_val_string[15];
	char nh_val_string_rev[15];
	long long int nh_val = 0;
	int i;
	int j;
	int k;
	char *t;

	j = 0;
	for ( i = icigar_length - 1 ; i >= 0 ; i-- )
	{
		if ( isdigit (icigar[i]) == 0 ) break;
		//nh_val = nh_val * 10 + icigar[i] - 48;
		nh_val_string[j++ ] = icigar[i];
	}
	icigar[i + 1] = '\0';
	nh_val_string[j] = '\0';

	k = j - 1;
	for ( i = 0 ; i < j ; i++ )
		nh_val_string_rev[i] = nh_val_string[k-- ];
	nh_val_string_rev[j] = '\0';
	nh_val = strtol (nh_val_string_rev , &t , 10);
	return nh_val;
}

void replaceCharacterInString (char *str, char ch_to_be_replaced, char replace_with, int str_length)
{
	int i;
	for ( i = 0 ; i < str_length ; i++ )
		if ( str[i] == ch_to_be_replaced ) str[i] = replace_with;
}

int findSamFormatFlagPairedEnded (char *icigar, int icigar_length, char *XS, struct Paired_Ended_Flag_to_Single_Character *samflag_dictionary, int number_of_unique_samformatflags, char samformatflag_replacer_characters[], char *character_to_be_replaced)
{
	int i, j;
	int samformatflag = -1;
	XS[0] = '-1';
	for ( i = 0 ; i < icigar_length ; i++ )
	{
		if ( ( icigar[i] != 'a' || icigar[i] != 't' || icigar[i] != 'g' || icigar[i] != 'c' ) && ! ( icigar[i] >= 48 && icigar[i] <= 57 ) && strchr (samformatflag_replacer_characters , icigar[i]) != NULL )
		{
			for ( j = 0 ; j < number_of_unique_samformatflags * 2 ; j++ )
			{
				if ( icigar[i] == samflag_dictionary->character[j] )
				{
					XS[0] = samflag_dictionary->direction[j];
					samformatflag = samflag_dictionary->samflags[j];
					( *character_to_be_replaced ) = icigar[i];
					printf ("\nicigar_character %c direction %c samformatflag %d " , icigar[i] , samflag_dictionary->direction[j] , samflag_dictionary->samflags[j]);
					break;
				}
			}
		}
		if ( samformatflag != -1 ) break;
	}
	return samformatflag;
}

int findSamFormatFlagSingleEnded (char *icigar, int icigar_length, char *XS)
{
	int i;
	int samformatflag;

	for ( i = 0 ; i < icigar_length ; i++ )
	{
		switch ( icigar[i] )
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

void generateReadSequenceAndMDString (struct Sam_Alignment *sam_alignment_instance, struct Whole_Genome_Sequence *whole_genome)
{
	int i;
	int j;
	int k;
	int chromosome_index;
	int read_from_genome_index;
	int distance_from_start_pos;
	int MD_index = 0;
	int splices_index = 0;
	int read_from_genome_including_deletions_index = 0;
	int num;
	int temp_index = 0;

	char read_from_genome[MAX_SEQ_LEN];
	char read_from_genome_including_deletions[MAX_SEQ_LEN];
	char MD[MAX_SEQ_LEN];
	char temp[100];

	/*
	 * Find the chromosome
	 */
	for ( chromosome_index = 0 ;
			chromosome_index < whole_genome->number_of_reference_sequences ;
			chromosome_index++ )
		if ( strcmp (sam_alignment_instance->reference_name , whole_genome->reference_sequence_name[chromosome_index]) == 0 )
			break;

	if ( chromosome_index == whole_genome->number_of_reference_sequences )
		printf ("\n Chromosome not found %s" , sam_alignment_instance->reference_name);

	/*
	 * Construct Sequence from CIGAR
	 */
	read_from_genome_index = 0;
	distance_from_start_pos = 0;
	MD[0] = '\0';
	/*
	 printf ("\nsam_alignment_instance->start_position %d" , sam_alignment_instance->start_position);
	 for ( i = 0 ; i < sam_alignment_instance->number_of_cigar_items ; i++ )
	 {
	 printf ("\n%d %c" , sam_alignment_instance->cigar_items[i].len , sam_alignment_instance->cigar_items[i].def);
	 }*/
	for ( i = 0 ; i < sam_alignment_instance->number_of_cigar_items ; i++ )
	{
		if ( sam_alignment_instance->cigar_items[i].def == 'S' )
		{
			//Soft clips - not present in genome. Hence ignore
			for ( j = 0 ; j < sam_alignment_instance->cigar_items[i].len ; j++ )
			{
				read_from_genome[read_from_genome_index++ ] = 'S';
				//read_from_genome_including_deletions[read_from_genome_including_deletions_index++] = 'S';
				/*
				 printf ("\nSoft clips j=%d" , j);
				 fflush (stdout);
				 */
			}
		}
		else if ( sam_alignment_instance->cigar_items[i].def == 'M' )
		{
			for ( j = distance_from_start_pos ;
					j < distance_from_start_pos + sam_alignment_instance->cigar_items[i].len ;
					j++ )
			{
				read_from_genome[read_from_genome_index++ ] = whole_genome->nucleotides[chromosome_index][sam_alignment_instance->start_position + j - 1];
				read_from_genome_including_deletions[read_from_genome_including_deletions_index++ ] = whole_genome->nucleotides[chromosome_index][sam_alignment_instance->start_position + j - 1];
				/*
				 printf ("\nPosition %d Length of genome %d read_from_genome_including_deletions_index %d j=%d" , sam_alignment_instance->start_position + j - 1 , strlen (whole_genome->nucleotides[chromosome_index]) , read_from_genome_including_deletions_index , j);
				 printf ("\nM whole_genome->nucleotides[chromosome_index][sam_alignment_instance->start_position + j - 1] %c" , whole_genome->nucleotides[chromosome_index][sam_alignment_instance->start_position + j - 1]);
				 fflush (stdout);
				 */
			}
			distance_from_start_pos += sam_alignment_instance->cigar_items[i].len;
		}
		else if ( sam_alignment_instance->cigar_items[i].def == 'N' )
		{
			// Genome sequence is not consumed but pointer will move ahead
			distance_from_start_pos += sam_alignment_instance->cigar_items[i].len;
		}
		else if ( sam_alignment_instance->cigar_items[i].def == 'I' )
		{
			for ( j = 0 ; j < sam_alignment_instance->cigar_items[i].len ; j++ )
			{
				read_from_genome[read_from_genome_index++ ] = 'I';
				read_from_genome_including_deletions[read_from_genome_including_deletions_index++ ] = 'I';
			}
		}
		else if ( sam_alignment_instance->cigar_items[i].def == 'D' )
		{

			for ( j = distance_from_start_pos ;
					j < distance_from_start_pos + sam_alignment_instance->cigar_items[i].len ;
					j++ )
			{
				read_from_genome_including_deletions[read_from_genome_including_deletions_index++ ] = whole_genome->nucleotides[chromosome_index][sam_alignment_instance->start_position + j - 1];
				/*
				 printf ("\nD whole_genome->nucleotides[chromosome_index][sam_alignment_instance->start_position + j - 1] %c" , whole_genome->nucleotides[chromosome_index][sam_alignment_instance->start_position + j - 1]);
				 */
			}
			// Genome sequence is not consumed but pointer will move ahead
			distance_from_start_pos += sam_alignment_instance->cigar_items[i].len;
		}
	}
	/*
	 printf ("\n");
	 for ( j = 0 ; j < read_from_genome_including_deletions_index ; j++ )
	 printf ("%c" , read_from_genome_including_deletions[j]);
	 printf ("\n");
	 fflush (stdout);
	 */
	/*
	 printf ("\nread_from_genome_including_deletions_index=%d" , read_from_genome_including_deletions_index);
	 fflush (stdout);*/
	read_from_genome[read_from_genome_index++ ] = '\0';
	read_from_genome_including_deletions[read_from_genome_including_deletions_index++ ] = '\0';
	for ( i = 0 ; sam_alignment_instance->seq[i] != '\0' ; i++ )
		if ( sam_alignment_instance->seq[i] == '-' )
			sam_alignment_instance->seq[i] = read_from_genome[i];

	strcpy (sam_alignment_instance->cigar_extended , read_from_genome_including_deletions);
	if ( strlen (sam_alignment_instance->md_extended) != strlen (read_from_genome_including_deletions) )
	{
		printSamAlignmentInstance (sam_alignment_instance , 1);
		printf ("\nCIGAR %s iCIGAR %s" , sam_alignment_instance->cigar , sam_alignment_instance->icigar);
		printf ("\nLengths do not match");
		printf ("\n%d %d" , strlen (sam_alignment_instance->md_extended) , strlen (read_from_genome_including_deletions));
		printf ("\n%s\n%s" , sam_alignment_instance->md_extended , read_from_genome_including_deletions);
		exit (1);
	}

	/*
	 * Construct MD String
	 */
	MD_index = 0;
	num = 0;
	sam_alignment_instance->tags[2].val[0] = '\0';
	for ( i = 0 ; sam_alignment_instance->md_extended[i] != '\0' ; i++ )
	{
		if ( sam_alignment_instance->md_extended[i] == '=' )
			num += 1;
		else if ( sam_alignment_instance->md_extended[i] == '-' )
		{
			if ( num > 0 )
			{
				sprintf (temp , "%d" , num);
				strcat (sam_alignment_instance->tags[2].val , temp);
			}
			num = 0;

			strcat (sam_alignment_instance->tags[2].val , "^");
			temp_index = 0;
			while ( sam_alignment_instance->md_extended[i] == '-' )
			{
				temp[temp_index++ ] = read_from_genome_including_deletions[i];
				i++;
			}
			temp[temp_index++ ] = '\0';
			strcat (sam_alignment_instance->tags[2].val , temp);
			i--;
		}
		else if ( sam_alignment_instance->md_extended[i] >= 38 && sam_alignment_instance->md_extended[i] <= 42 ) // Mismatches
		{
			sprintf (temp , "%d" , num);
			strcat (sam_alignment_instance->tags[2].val , temp);
			num = 0;
			temp[0] = read_from_genome_including_deletions[i];
			temp[1] = '\0';
			strcat (sam_alignment_instance->tags[2].val , temp);
		}
	}

	for ( i = 0 ; sam_alignment_instance->tags[2].val[i] != '\0' ; i++ )
		if ( sam_alignment_instance->tags[2].val[i] >= 'A' && sam_alignment_instance->tags[2].val[i] <= 'Z' && sam_alignment_instance->tags[2].val[i] != 'A' && sam_alignment_instance->tags[2].val[i] != 'T' && sam_alignment_instance->tags[2].val[i] != 'G' && sam_alignment_instance->tags[2].val[i] != 'C' )
			sam_alignment_instance->tags[2].val[i] = 'N';

	if ( num > 0 )
	{
		sprintf (temp , "%d" , num);
		strcat (sam_alignment_instance->tags[2].val , temp);
	}
	/*
	 printf("\n%d %d", read_from_genome_index, read_from_genome_including_deletions_index);
	 printf("\n%s\n%s\n%s", sam_alignment_instance->soft_clips_removed_seq, read_from_genome_including_deletions, read_from_genome);
	 printf("\n");
	 */
}

void writeAlignmentToFileSingleEnded (struct Sam_Alignment *sam_alignment, short int flag_ignore_sequence_information, int number_of_repititions_of_the_same_reads, char *read_prefix, FILE *fhw, FILE *fhr_qual, short int flag_save_all_quality_scores)
{
	int i;

	char line_to_be_written_to_file[MAX_GENERAL_LEN];
	char temp[100];
	char *buffer;

	size_t len = 0;
	ssize_t line_len;

	for ( i = 0 ; i < number_of_repititions_of_the_same_reads ; i++ )
	{
		line_to_be_written_to_file[0] = '\0';
		strcat (line_to_be_written_to_file , read_prefix);
		strcat (line_to_be_written_to_file , sam_alignment->read_name);
		sprintf (temp , "%d" , i + 1);
		strcat (line_to_be_written_to_file , "_");
		strcat (line_to_be_written_to_file , temp);
		strcat (line_to_be_written_to_file , "\t");

		sprintf (temp , "%d" , sam_alignment->samflag);
		strcat (line_to_be_written_to_file , temp);

		strcat (line_to_be_written_to_file , "\t");
		strcat (line_to_be_written_to_file , sam_alignment->reference_name);

		strcat (line_to_be_written_to_file , "\t");
		sprintf (temp , "%d" , sam_alignment->start_position);
		strcat (line_to_be_written_to_file , temp);

		strcat (line_to_be_written_to_file , "\t");
		strcat (line_to_be_written_to_file , "255");

		strcat (line_to_be_written_to_file , "\t");
		strcat (line_to_be_written_to_file , sam_alignment->cigar);

		strcat (line_to_be_written_to_file , "\t");
		strcat (line_to_be_written_to_file , "*");

		strcat (line_to_be_written_to_file , "\t");
		strcat (line_to_be_written_to_file , "0");

		strcat (line_to_be_written_to_file , "\t");
		strcat (line_to_be_written_to_file , "0");

		strcat (line_to_be_written_to_file , "\t");
		strcat (line_to_be_written_to_file , sam_alignment->seq);

		strcat (line_to_be_written_to_file , "\t");
		if ( flag_save_all_quality_scores == 1 && ( line_len = getline ( &buffer , &len , fhr_qual) ) != -1 )
		{
			buffer[strlen (buffer) - 1] = '\0';
			strcat (line_to_be_written_to_file , buffer);
		}
		else strcat (line_to_be_written_to_file , sam_alignment->qual);

		strcat (line_to_be_written_to_file , "\t");
		//Tags
		strcat (line_to_be_written_to_file , "NH:i:");
		strcat (line_to_be_written_to_file , sam_alignment->tags[0].val);
		strcat (line_to_be_written_to_file , "\t");

		if ( strcmp (sam_alignment->tags[1].val , ".") != 0 )
		{
			strcat (line_to_be_written_to_file , "XS:A:");
			strcat (line_to_be_written_to_file , sam_alignment->tags[1].val);
			strcat (line_to_be_written_to_file , "\t");
		}

		if ( flag_ignore_sequence_information == 0 )
		{
			strcat (line_to_be_written_to_file , "MD:Z:");
			strcat (line_to_be_written_to_file , sam_alignment->tags[2].val);
			strcat (line_to_be_written_to_file , "\t");
		}

		strcat (line_to_be_written_to_file , "\n");
		fprintf (fhw , "%s" , line_to_be_written_to_file);
	}
}

void writeAlignmentToFilePairedEnded (struct Sam_Alignment *sam_alignment, short int flag_ignore_sequence_information, int number_of_repititions_of_the_same_reads, FILE *fhw, FILE *fhr_qual, short int flag_save_all_quality_scores, char **read_names, int *read_names_index)
{
	int i;

	char line_to_be_written_to_file[MAX_GENERAL_LEN];
	char temp[100];
	char *buffer;

	size_t len = 0;
	ssize_t line_len;

	for ( i = 0 ; i < number_of_repititions_of_the_same_reads ; i++ )
	{
		line_to_be_written_to_file[0] = '\0';
		strcat (line_to_be_written_to_file , read_names[ *read_names_index]);
		( *read_names_index )++;
		//sprintf (temp , "%d" , i + 1);
		//strcat (line_to_be_written_to_file , "_");
		//strcat (line_to_be_written_to_file , temp);
		strcat (line_to_be_written_to_file , "\t");

		sprintf (temp , "%d" , sam_alignment->samflag);
		strcat (line_to_be_written_to_file , temp);

		strcat (line_to_be_written_to_file , "\t");
		strcat (line_to_be_written_to_file , sam_alignment->reference_name);

		strcat (line_to_be_written_to_file , "\t");
		sprintf (temp , "%d" , sam_alignment->start_position);
		strcat (line_to_be_written_to_file , temp);

		strcat (line_to_be_written_to_file , "\t");
		strcat (line_to_be_written_to_file , "255");

		strcat (line_to_be_written_to_file , "\t");
		strcat (line_to_be_written_to_file , sam_alignment->cigar);

		strcat (line_to_be_written_to_file , "\t");
		strcat (line_to_be_written_to_file , "*");

		strcat (line_to_be_written_to_file , "\t");
		strcat (line_to_be_written_to_file , "0");

		strcat (line_to_be_written_to_file , "\t");
		strcat (line_to_be_written_to_file , "0");

		strcat (line_to_be_written_to_file , "\t");
		strcat (line_to_be_written_to_file , sam_alignment->seq);

		strcat (line_to_be_written_to_file , "\t");
		if ( flag_save_all_quality_scores == 1 && ( line_len = getline ( &buffer , &len , fhr_qual) ) != -1 )
		{
			buffer[strlen (buffer) - 1] = '\0';
			strcat (line_to_be_written_to_file , buffer);
		}
		else strcat (line_to_be_written_to_file , sam_alignment->qual);

		strcat (line_to_be_written_to_file , "\t");
		//Tags
		strcat (line_to_be_written_to_file , "NH:i:");
		strcat (line_to_be_written_to_file , sam_alignment->tags[0].val);
		strcat (line_to_be_written_to_file , "\t");

		if ( strcmp (sam_alignment->tags[1].val , ".") != 0 )
		{
			strcat (line_to_be_written_to_file , "XS:A:");
			strcat (line_to_be_written_to_file , sam_alignment->tags[1].val);
			strcat (line_to_be_written_to_file , "\t");
		}

		if ( flag_ignore_sequence_information == 0 )
		{
			strcat (line_to_be_written_to_file , "MD:Z:");
			strcat (line_to_be_written_to_file , sam_alignment->tags[2].val);
			strcat (line_to_be_written_to_file , "\t");
		}

		strcat (line_to_be_written_to_file , "\n");
		fprintf (fhw , "%s" , line_to_be_written_to_file);
	}
}

void convertToAlignmentPairedEnded (struct Sam_Alignment *sam_alignment_instance, struct Whole_Genome_Sequence *whole_genome, char **split_on_tab, char **split_on_dash, char **split_on_comma, char **read_names, char *default_quality_value, short int flag_ignore_mismatches, short int flag_ignore_soft_clippings, short int flag_ignore_unmapped_sequences, short int flag_ignore_quality_score, short int flag_ignore_sequence_information, unsigned long long int *read_number, unsigned long long int *total_mapped_reads, FILE *fhw, FILE *fhr_qual, short int flag_save_all_quality_scores, short int number_of_columns, unsigned long long int curr_position, char *chromosome, struct Paired_Ended_Flag_to_Single_Character *samflag_dictionary, int number_of_unique_samformatflags, char samformatflag_replacer_characters[])
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/

	int number_of_distinct_cigars_in_a_line;
	int number_of_repititions_of_the_same_reads;
	int samformatflag;

	int i, j;
	int read_names_index = 0;

	char *temp; //Useless
	char *distinct_icigars_in_a_line;
	char *icigar;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/

	/********************************************************************/

	if ( number_of_columns == 2 )
	{
		int max_number_of_commas = 0, number_of_commas = 0;
		for ( i = 0 ; split_on_tab[0][i] != '\0' ; i++ )
			if ( split_on_tab[0][i] == ',' ) number_of_commas++;
		number_of_distinct_cigars_in_a_line = splitByDelimiter (split_on_tab[0] , ',' , split_on_comma);
		splitByDelimiter (split_on_tab[1] , ',' , read_names);
	}
	else if ( number_of_columns == 3 )
	{
		int max_number_of_commas = 0, number_of_commas = 0;
		for ( i = 0 ; split_on_tab[1][i] != '\0' ; i++ )
			if ( split_on_tab[1][i] == ',' ) number_of_commas++;
		number_of_distinct_cigars_in_a_line = splitByDelimiter (split_on_tab[1] , ',' , split_on_comma);
		splitByDelimiter (split_on_tab[2] , ',' , read_names);
	}
	for ( j = 0 ; j < number_of_distinct_cigars_in_a_line ; j++ )
	{
		splitByDelimiter (split_on_comma[j] , '-' , split_on_dash);
		number_of_repititions_of_the_same_reads = strtol (split_on_dash[1] , &temp , 10);
		sam_alignment_instance->start_position = curr_position;

		if ( split_on_comma[j][1] == '-' && isalpha (split_on_dash[0][0]) != 0 )
		{
			// Use the same cigar
			sprintf (temp , "%d" , *read_number);
			( *read_number )++;
			strcpy (sam_alignment_instance->read_name , temp);
			/*if ( sam_alignment_instance->start_position == 27381 && strcmp (sam_alignment_instance->reference_name , "1") == 0 )
			 {
			 printf ("\nSame iCIGAR");
			 printf ("\nsplit_on_comma[j] = %s" , split_on_comma[j]);
			 printf ("\nWeird Location cigar %s number_of_repititions_of_the_same_reads %d" , sam_alignment_instance->icigar , number_of_repititions_of_the_same_reads);
			 printf ("\nMD String: %s" , sam_alignment_instance->tags[2].val);
			 printf ("\n==============================================================================================================================");
			 fflush (stdout);
			 }*/
		}
		else
		{
			strcpy (sam_alignment_instance->icigar , split_on_dash[0]);
			//printf ("\nj=%d number_of_distinct_cigars_in_a_line=%d Inside ICIGAR %s" , j , number_of_distinct_cigars_in_a_line , sam_alignment_instance->icigar);
			//fflush (stdout);
			convertIcigarToCigarandMDPairedEnded (whole_genome , sam_alignment_instance , chromosome , flag_ignore_mismatches , flag_ignore_soft_clippings , flag_ignore_unmapped_sequences , flag_ignore_quality_score , flag_ignore_sequence_information , default_quality_value , samflag_dictionary , number_of_unique_samformatflags , samformatflag_replacer_characters);

			/*if ( sam_alignment_instance->start_position == 27381 && strcmp (sam_alignment_instance->reference_name , "1") == 0 )
			 {
			 printf ("\nsplit_on_comma[j] = %s" , split_on_comma[j]);
			 printf ("\nWeird Location cigar %s number_of_repititions_of_the_same_reads %d" , sam_alignment_instance->icigar , number_of_repititions_of_the_same_reads);
			 printf ("\nMD String: %s" , sam_alignment_instance->tags[2].val);
			 printf ("\n split_on_dash[0] %s split_on_dash[1] %s" , split_on_dash[0] , split_on_dash[1]);
			 printf ("\nsplit_on_dash[0][1] == '-' %d isalpha (split_on_dash[0][0])  %d" , split_on_dash[0][1] == '-' , isalpha (split_on_dash[0][0]));
			 printf ("\n==============================================================================================================================");
			 fflush (stdout);
			 }*/
			sprintf (temp , "%d" , *read_number);
			( *read_number )++;
			strcpy (sam_alignment_instance->read_name , temp);
		}
		writeAlignmentToFilePairedEnded (sam_alignment_instance , flag_ignore_sequence_information , number_of_repititions_of_the_same_reads , fhw , fhr_qual , flag_save_all_quality_scores , read_names , &read_names_index);
		( *total_mapped_reads ) += number_of_repititions_of_the_same_reads;
	}
}

void convertToAlignmentSingleEnded (struct Sam_Alignment *sam_alignment_instance, struct Whole_Genome_Sequence *whole_genome, char **split_on_tab, char **split_on_dash, char **split_on_comma, char *default_quality_value, short int flag_ignore_mismatches, short int flag_ignore_soft_clippings, short int flag_ignore_unmapped_sequences, short int flag_ignore_quality_score, short int flag_ignore_sequence_information, unsigned long long int *read_number, unsigned long long int *total_mapped_reads, char *read_prefix, FILE *fhw, FILE *fhr_qual, short int flag_save_all_quality_scores, short int number_of_columns, unsigned long long int curr_position, char *chromosome)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/

	int number_of_distinct_cigars_in_a_line;
	int number_of_repititions_of_the_same_reads;
	int samformatflag;

	int i, j;

	char *temp; //Useless
	char *distinct_icigars_in_a_line;
	char *icigar;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/

	/********************************************************************/

	if ( number_of_columns == 1 )
	{
		int max_number_of_commas = 0, number_of_commas = 0;
		for ( i = 0 ; split_on_tab[0][i] != '\0' ; i++ )
			if ( split_on_tab[0][i] == ',' ) number_of_commas++;
		number_of_distinct_cigars_in_a_line = splitByDelimiter (split_on_tab[0] , ',' , split_on_comma);
	}
	else if ( number_of_columns == 2 )
	{
		int max_number_of_commas = 0, number_of_commas = 0;
		for ( i = 0 ; split_on_tab[1][i] != '\0' ; i++ )
			if ( split_on_tab[1][i] == ',' ) number_of_commas++;
		number_of_distinct_cigars_in_a_line = splitByDelimiter (split_on_tab[1] , ',' , split_on_comma);
	}
	for ( j = 0 ; j < number_of_distinct_cigars_in_a_line ; j++ )
	{
		splitByDelimiter (split_on_comma[j] , '-' , split_on_dash);
		number_of_repititions_of_the_same_reads = strtol (split_on_dash[1] , &temp , 10);
		sam_alignment_instance->start_position = curr_position;

		if ( split_on_comma[j][1] == '-' && isalpha (split_on_dash[0][0]) != 0 )
		{
			// Use the same cigar
			sprintf (temp , "%d" , *read_number);
			( *read_number )++;
			strcpy (sam_alignment_instance->read_name , temp);
			/*if ( sam_alignment_instance->start_position == 27381 && strcmp (sam_alignment_instance->reference_name , "1") == 0 )
			 {
			 printf ("\nSame iCIGAR");
			 printf ("\nsplit_on_comma[j] = %s" , split_on_comma[j]);
			 printf ("\nWeird Location cigar %s number_of_repititions_of_the_same_reads %d" , sam_alignment_instance->icigar , number_of_repititions_of_the_same_reads);
			 printf ("\nMD String: %s" , sam_alignment_instance->tags[2].val);
			 printf ("\n==============================================================================================================================");
			 fflush (stdout);
			 }*/
		}
		else
		{
			strcpy (sam_alignment_instance->icigar , split_on_dash[0]);
			//printf ("\nj=%d number_of_distinct_cigars_in_a_line=%d Inside ICIGAR %s" , j , number_of_distinct_cigars_in_a_line , sam_alignment_instance->icigar);
			//fflush (stdout);

			convertIcigarToCigarandMDSingleEnded (whole_genome , sam_alignment_instance , chromosome , flag_ignore_mismatches , flag_ignore_soft_clippings , flag_ignore_unmapped_sequences , flag_ignore_quality_score , flag_ignore_sequence_information , default_quality_value);

			/*if ( sam_alignment_instance->start_position == 27381 && strcmp (sam_alignment_instance->reference_name , "1") == 0 )
			 {
			 printf ("\nsplit_on_comma[j] = %s" , split_on_comma[j]);
			 printf ("\nWeird Location cigar %s number_of_repititions_of_the_same_reads %d" , sam_alignment_instance->icigar , number_of_repititions_of_the_same_reads);
			 printf ("\nMD String: %s" , sam_alignment_instance->tags[2].val);
			 printf ("\n split_on_dash[0] %s split_on_dash[1] %s" , split_on_dash[0] , split_on_dash[1]);
			 printf ("\nsplit_on_dash[0][1] == '-' %d isalpha (split_on_dash[0][0])  %d" , split_on_dash[0][1] == '-' , isalpha (split_on_dash[0][0]));
			 printf ("\n==============================================================================================================================");
			 fflush (stdout);
			 }*/
			sprintf (temp , "%d" , *read_number);
			( *read_number )++;
			strcpy (sam_alignment_instance->read_name , temp);
		}
		writeAlignmentToFileSingleEnded (sam_alignment_instance , flag_ignore_sequence_information , number_of_repititions_of_the_same_reads , read_prefix , fhw , fhr_qual , flag_save_all_quality_scores);
		( *total_mapped_reads ) += number_of_repititions_of_the_same_reads;
	}
}

void writeSequenceHeaders (FILE *fhw, char *genome_filename)
{
	FILE *fhr;

	char temp[100];
	char faidx_filename[1000];
	char *buffer = NULL;
	char **split_on_tab;
	char line_to_be_written_to_file[10000];

	size_t len = 0;
	ssize_t line_len;

	int i;

	strcpy (faidx_filename , genome_filename);
	strcat (faidx_filename , ".fai");

	fhr = fopen (faidx_filename , "r");
	if ( fhr == NULL )
	{
		printf ("Error! File not found");
		exit (1);
	}

	split_on_tab = ( char** ) malloc (sizeof(char*) * ROWS);
	for ( i = 0 ; i < ROWS ; i++ )
		split_on_tab[i] = ( char* ) malloc (sizeof(char) * COLS);

	strcpy (line_to_be_written_to_file , "@HD");
	strcat (line_to_be_written_to_file , "\t");
	strcat (line_to_be_written_to_file , "VN:1.4");
	strcat (line_to_be_written_to_file , "\t");
	strcat (line_to_be_written_to_file , "SO:coordinate");
	strcat (line_to_be_written_to_file , "\n");
//fprintf (fhw , "%s" , line_to_be_written_to_file);
	printf ("%s" , line_to_be_written_to_file);
	while ( ( line_len = getline ( &buffer , &len , fhr) ) != -1 )
	{
		splitByDelimiter (buffer , '\t' , split_on_tab);
		strcpy (line_to_be_written_to_file , "@SQ");
		strcat (line_to_be_written_to_file , "\t");
		strcat (line_to_be_written_to_file , "SN:");
		strcat (line_to_be_written_to_file , split_on_tab[0]);
		strcat (line_to_be_written_to_file , "\t");
		strcat (line_to_be_written_to_file , "LN:");
		strcat (line_to_be_written_to_file , split_on_tab[1]);
		strcat (line_to_be_written_to_file , "\t");
		strcat (line_to_be_written_to_file , "\n");
		//fprintf (fhw , "%s" , line_to_be_written_to_file);
		printf ("%s" , line_to_be_written_to_file);
	}

	fclose (fhr);
}

#endif /* ABRIDGE_SRC_FUNCTIONS_DEFINITIONS_H_ */

