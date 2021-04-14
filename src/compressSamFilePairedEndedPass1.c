# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include <math.h>
# include <stdbool.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

int main (int argc, char *argv[])
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char input_samfilename[FILENAME_LENGTH];
	char output_abridgefilename[FILENAME_LENGTH];
	char genome_filename[FILENAME_LENGTH];
	char unmapped_filename[FILENAME_LENGTH];
	char name_of_file_with_max_commas[FILENAME_LENGTH];
	char name_of_file_with_quality_scores[FILENAME_LENGTH];
	char *temp; //Required for strtoi

	short int flag_ignore_soft_clippings;
	short int flag_ignore_mismatches;
	short int flag_ignore_quality_score;
	short int flag_ignore_unmapped_sequences;
	short int run_diagnostics;
	short int save_all_quality_scores;
	short int save_exact_quality_scores;

	long long int max_input_reads_in_a_single_nucl_loc;
	long long int max_number_of_alignments;

	int i;
	int max_read_length;

	struct All_Relevant_Info_PE_per_Alignment **mega_array;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	strcpy(genome_filename , argv[1]);
	flag_ignore_soft_clippings = strtol (argv[2] , &temp , 10);
	flag_ignore_mismatches = strtol (argv[3] , &temp , 10);
	flag_ignore_quality_score = strtol (argv[4] , &temp , 10);
	flag_ignore_unmapped_sequences = strtol (argv[5] , &temp , 10);
	strcpy(input_samfilename , argv[6]);
	strcpy(output_abridgefilename , argv[7]);
	strcpy(unmapped_filename , argv[8]);
	run_diagnostics = strtol (argv[9] , &temp , 10);
	max_input_reads_in_a_single_nucl_loc = strtol (argv[10] , &temp , 10);
	strcpy(name_of_file_with_max_commas , argv[11]);
	save_all_quality_scores = strtol (argv[12] , &temp , 10);
	save_exact_quality_scores = strtol (argv[13] , &temp , 10);
	strcpy(name_of_file_with_quality_scores , argv[14]);
	max_number_of_alignments = strtol (argv[15] , &temp , 10);
	max_read_length = strtol (argv[16] , &temp , 10);

	printf ("\nmax_number_of_alignments %d" , max_number_of_alignments);
	mega_array = ( struct All_Relevant_Info_PE_per_Alignment** ) malloc (sizeof(struct All_Relevant_Info_PE_per_Alignment*) * max_number_of_alignments);
	for ( i = 0 ; i < max_number_of_alignments ; i++ )
		mega_array[i] = allocateMemoryAll_Relevant_Info_PE_per_Alignment (max_read_length);
	printf ("\nSize of mega array %d" , sizeof ( mega_array ));
	/********************************************************************/

	return 0;
}
