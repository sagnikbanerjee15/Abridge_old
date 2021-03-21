# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

long long int findEndPointOfiCIGAR ( char *icigar )
{
	long long int end_point = 0;
	long long int num = 0;
	int start = -1;
	int end;
	int i , j , k;

	// Removing Left Soft clips
	for ( i = 0 ; icigar[i] != '\0' ; i++ )
		if ( isdigit ( icigar[i] ) != 0 && start == -1 )
		{
			start = i;
			break;
		}
	// Ignore NH values
	for ( j = strlen ( icigar ) - 1 ; j >= 0 ; j-- )
		if ( isdigit ( icigar[j] ) == 0 ) break;

	// Removing Right Soft clips
	for ( k = j ; k >= 0 ; k-- )
		if ( icigar[k] != 'A' && icigar[k] != 'T' && icigar[k] != 'G' && icigar[k] != 'C' && icigar[k] != 'N' )
			break;
	end = k;

	for ( i = start ; i <= end ; i++ )
	{
		if ( icigar[i] == '!' || icigar[i] == '"' || icigar[i] == '#' || icigar[i] == '$' || icigar[i] == '%' ) // Insertion characters
		{
			end_point += num;
			num = 0;
		}
		else if ( icigar[i] == '&' || icigar[i] == '\'' || icigar[i] == '(' || icigar[i] == ')' || icigar[i] == '*' ) // Mismatch characters
		{
			end_point += 1;
			num = 0;
		}
		else if ( isdigit ( icigar[i] ) != 0 )
		{
			num = num * 10 + icigar[i] - 48;
		}
		else if ( isalpha ( icigar[i] ) != 0 )
		{
			end_point += num;
			num = 0;
		}
	}
	/*if (start != 0 && end != j)
	 {
	 printf("\nstart=%d end=%d j=%d icigar: %s, end_point= %lld icigar_portion: ", start, end, j, icigar, end_point);
	 for (i = start; i <= end; i++)
	 printf("%c", icigar[i]);
	 }*/

	return end_point;
}

void findFarthestMapping ( long long int *start, long long int *end, char **split_icigar_field, int number_of_fields, char **split_icigar_and_num_reads, int print )
{
	int i;
	*end = -1;
	/*if (print)
	 {
	 printf("\n New Line");
	 printf("\n Number of fields: %d", number_of_fields);
	 for (i = 0; i < number_of_fields; i++)
	 printf("\n %d. %s", i, split_icigar_field[i]);

	 }*/
	for ( i = 0 ; i < number_of_fields ; i++ )
	{
		splitByDelimiter ( split_icigar_field[i] , '-' , split_icigar_and_num_reads ); // will always have 2 fields
		if ( strlen ( split_icigar_and_num_reads[0] ) != 1 )
		{
			if ( print )
				printf ( "\n Start:%lld Length of read:%d End:%lld icigar:%s" , *start , findEndPointOfiCIGAR ( split_icigar_and_num_reads[0] ) , *start + findEndPointOfiCIGAR ( split_icigar_and_num_reads[0] ) - 1 , split_icigar_field[i] );
			if ( *start + findEndPointOfiCIGAR ( split_icigar_and_num_reads[0] ) - 1 > *end )
				*end = *start + findEndPointOfiCIGAR ( split_icigar_and_num_reads[0] ) - 1;

		}
	}
}

void findContinousClusters ( char *input_filename, char *output_filename, long long int max_input_reads_in_a_single_nucl_loc )
{
	/*
	 *
	 */

	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	FILE *fhr;
	FILE *fhw;

	int i;
	int number_of_fields;

	long long int start_position_of_cluster = -1;
	long long int end_position_of_cluster = -1;
	long long int start_position_of_read;
	long long int end_position_of_read;
	long long int num_lines_read = 0;
	long long int byte_number_start_cluster = -1;
	long long int byte_number_end_cluster = -1;
	long int file_position;

	size_t len = 0;
	ssize_t line_len;
	ssize_t line_len_previous;

	char *temp; //Useless
	char *line = NULL; // for reading each line
	char *current_reference_sequence;
	char *line_to_be_written_to_file;
	char **split_on_tab;
	char **split_on_colon;
	char **split_icigar_field;
	char **split_icigar_and_num_reads;

	struct Pass2_Compressed_DS *pass2_compressed_ds_instance;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr = fopen ( input_filename , "r" );
	if ( fhr == NULL )
	{
		printf ( "Error! File not found" );
		exit ( 1 );
	}
	fhw = fopen ( output_filename , "w" );
	if ( fhw == NULL )
	{
		printf ( "File cannot be created" );
		exit ( 1 );
	}
	//printf ( "\n max_input_reads_in_a_single_nucl_loc = %d\n" , max_input_reads_in_a_single_nucl_loc );
	//fflush ( stdout );

	max_input_reads_in_a_single_nucl_loc += 5;
	split_on_tab = ( char** ) malloc ( sizeof(char*) * 5 );
	for ( i = 0 ; i < 5 ; i++ )
		split_on_tab[i] = ( char* ) malloc ( sizeof(char) * max_input_reads_in_a_single_nucl_loc * 10000 );
	split_on_colon = ( char** ) malloc ( sizeof(char*) * 5 );
	for ( i = 0 ; i < 5 ; i++ )
		split_on_colon[i] = ( char* ) malloc ( sizeof(char) * 1000 );

	split_icigar_field = ( char** ) malloc ( sizeof(char*) * max_input_reads_in_a_single_nucl_loc );
	for ( i = 0 ; i < max_input_reads_in_a_single_nucl_loc ; i++ )
		split_icigar_field[i] = ( char* ) malloc ( sizeof(char) * 10000 );

	split_icigar_and_num_reads = ( char** ) malloc ( sizeof(char*) * 5 );
	for ( i = 0 ; i < 5 ; i++ )
		split_icigar_and_num_reads[i] = ( char* ) malloc ( sizeof(char) * 100000 );

	line_to_be_written_to_file = ( char* ) malloc ( sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE );
	pass2_compressed_ds_instance = allocateMemoryPass2_Compressed_DS ();
	current_reference_sequence = ( char* ) malloc ( sizeof(char) * MAX_REFERENCE_SEQ_LENGTH );
	temp = ( char* ) malloc ( sizeof(char) * MAX_GENERAL_LEN );
	line_to_be_written_to_file[0] = '\0';
	/********************************************************************/

	line_len_previous = -1;
	while ( ( line_len = getline ( &line , &len , fhr ) ) != -1 )
	{
		num_lines_read++;
		if ( num_lines_read == 1 )
		{
			fprintf ( fhw , "%s" , line );
			continue;
		}
		if ( line[0] == '@' )
		{
			if ( start_position_of_cluster != -1 && end_position_of_cluster != -1 )
			{
				byte_number_end_cluster = ftell ( fhr ) - line_len;
				//printf("\n%s\t%lld\t%lld\t%ld", current_reference_sequence, start_position_of_cluster, end_position_of_cluster, ftell(fhr));
				strcpy( line_to_be_written_to_file , current_reference_sequence );
				strcat( line_to_be_written_to_file , "\t" );

				sprintf( temp , "%d" , start_position_of_cluster );
				strcat( line_to_be_written_to_file , temp );
				strcat( line_to_be_written_to_file , "\t" );

				sprintf( temp , "%d" , end_position_of_cluster );
				strcat( line_to_be_written_to_file , temp );
				strcat( line_to_be_written_to_file , "\t" );

				sprintf( temp , "%ld" , byte_number_start_cluster );
				strcat( line_to_be_written_to_file , temp );
				strcat( line_to_be_written_to_file , "\t" );

				sprintf( temp , "%ld" , byte_number_end_cluster );
				strcat( line_to_be_written_to_file , temp );
				strcat( line_to_be_written_to_file , "\t" );

				strcat( line_to_be_written_to_file , "\n" );
				fprintf ( fhw , "%s" , line_to_be_written_to_file );

				file_position = -1;

				byte_number_start_cluster = ftell ( fhr ) - line_len;
			}
			// New chromosome encountered;
			splitByDelimiter ( line , '\t' , split_on_tab );
			splitByDelimiter ( split_on_tab[1] , ':' , split_on_colon );
			strcpy( current_reference_sequence , split_on_colon[1] );
			file_position = ftell ( fhr );

			line_len = getline ( &line , &len , fhr );
			num_lines_read++;
			number_of_fields = splitByDelimiter ( line , '\t' , split_on_tab );

			if ( number_of_fields == 1 )
			{

				start_position_of_read = 1;
				number_of_fields = splitByDelimiter ( split_on_tab[0] , ',' , split_icigar_field );
				//continue;
				/*
				 if (strcmp(current_reference_sequence, "MT") == 0) findFarthestMapping(&start_position_of_read, &end_position_of_read, split_icigar_field, number_of_fields, split_icigar_and_num_reads, 1);
				 else findFarthestMapping(&start_position_of_read, &end_position_of_read, split_icigar_field, number_of_fields, split_icigar_and_num_reads, 0);
				 */
				findFarthestMapping ( &start_position_of_read , &end_position_of_read , split_icigar_field , number_of_fields , split_icigar_and_num_reads , 0 );
				start_position_of_cluster = start_position_of_read;
				end_position_of_cluster = end_position_of_read;
				byte_number_start_cluster = file_position;
				//printf("\n Trouble %s", line);
			}
			else
			{
				start_position_of_read = strtol ( split_on_tab[0] , &temp , 10 );
				number_of_fields = splitByDelimiter ( split_on_tab[1] , ',' , split_icigar_field );
				/*
				 if (strcmp(current_reference_sequence, "MT") == 0) findFarthestMapping(&start_position_of_read, &end_position_of_read, split_icigar_field, number_of_fields, split_icigar_and_num_reads, 1);
				 else findFarthestMapping(&start_position_of_read, &end_position_of_read, split_icigar_field, number_of_fields, split_icigar_and_num_reads, 0);
				 */
				findFarthestMapping ( &start_position_of_read , &end_position_of_read , split_icigar_field , number_of_fields , split_icigar_and_num_reads , 0 );
				start_position_of_cluster = start_position_of_read;
				end_position_of_cluster = end_position_of_read;
				byte_number_start_cluster = file_position;
			}
		}
		else
		{
			//continue;
			number_of_fields = splitByDelimiter ( line , '\t' , split_on_tab );
			if ( number_of_fields == 1 )
			{
				start_position_of_read++;
				number_of_fields = splitByDelimiter ( split_on_tab[0] , ',' , split_icigar_field );
			}
			else
			{
				start_position_of_read += strtol ( split_on_tab[0] , &temp , 10 );
				number_of_fields = splitByDelimiter ( split_on_tab[1] , ',' , split_icigar_field );
			}

			/*
			 if (strcmp(current_reference_sequence, "MT") == 0) findFarthestMapping(&start_position_of_read, &end_position_of_read, split_icigar_field, number_of_fields, split_icigar_and_num_reads, 1);
			 else findFarthestMapping(&start_position_of_read, &end_position_of_read, split_icigar_field, number_of_fields, split_icigar_and_num_reads, 0);
			 */
			findFarthestMapping ( &start_position_of_read , &end_position_of_read , split_icigar_field , number_of_fields , split_icigar_and_num_reads , 0 );
			// New cluster is found
			if ( start_position_of_read > end_position_of_cluster )
			{
				byte_number_end_cluster = ftell ( fhr ) - line_len;
				//printf("\n%s\t%lld\t%lld\t%ld", current_reference_sequence, start_position_of_cluster, end_position_of_cluster, ftell(fhr));
				strcpy( line_to_be_written_to_file , current_reference_sequence );
				strcat( line_to_be_written_to_file , "\t" );

				sprintf( temp , "%d" , start_position_of_cluster );
				strcat( line_to_be_written_to_file , temp );
				strcat( line_to_be_written_to_file , "\t" );

				sprintf( temp , "%d" , end_position_of_cluster );
				strcat( line_to_be_written_to_file , temp );
				strcat( line_to_be_written_to_file , "\t" );

				sprintf( temp , "%ld" , byte_number_start_cluster );
				strcat( line_to_be_written_to_file , temp );
				strcat( line_to_be_written_to_file , "\t" );

				sprintf( temp , "%ld" , byte_number_end_cluster );
				strcat( line_to_be_written_to_file , temp );
				strcat( line_to_be_written_to_file , "\t" );

				strcat( line_to_be_written_to_file , "\n" );
				fprintf ( fhw , "%s" , line_to_be_written_to_file );

				start_position_of_cluster = start_position_of_read;
				end_position_of_cluster = end_position_of_read;
				file_position = ftell ( fhr ) - line_len;

				byte_number_start_cluster = ftell ( fhr ) - line_len;

			}
			else if ( end_position_of_read > end_position_of_cluster )
			{
				end_position_of_cluster = end_position_of_read;
				//if (strcmp(current_reference_sequence, "MT") == 0) printf("\nEnd_position_of_cluster: %d", end_position_of_cluster);
			}
		}
		line_len_previous = line_len;
		//if (strcmp(current_reference_sequence, "MT") == 0) printf("\nLine read from file:%s", line);

	}
	byte_number_end_cluster = ftell ( fhr ) - line_len - line_len_previous;
	//printf("\n%s\t%lld\t%lld\t%ld", current_reference_sequence, start_position_of_cluster, end_position_of_cluster, ftell(fhr));
	strcpy( line_to_be_written_to_file , current_reference_sequence );
	strcat( line_to_be_written_to_file , "\t" );

	sprintf( temp , "%d" , start_position_of_cluster );
	strcat( line_to_be_written_to_file , temp );
	strcat( line_to_be_written_to_file , "\t" );

	sprintf( temp , "%d" , end_position_of_cluster );
	strcat( line_to_be_written_to_file , temp );
	strcat( line_to_be_written_to_file , "\t" );

	sprintf( temp , "%ld" , byte_number_start_cluster );
	strcat( line_to_be_written_to_file , temp );
	strcat( line_to_be_written_to_file , "\t" );

	sprintf( temp , "%ld" , byte_number_end_cluster );
	strcat( line_to_be_written_to_file , temp );
	strcat( line_to_be_written_to_file , "\t" );

	strcat( line_to_be_written_to_file , "\n" );
	fprintf ( fhw , "%s" , line_to_be_written_to_file );

	start_position_of_cluster = start_position_of_read;
	end_position_of_cluster = end_position_of_read;
	file_position = ftell ( fhr ) - line_len;
}

int main ( int argc, char *argv[] )
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char input_filename[FILENAME_LENGTH];
	char output_filename[FILENAME_LENGTH];
	char *temp;

	long long int max_input_reads_in_a_single_nucl_loc;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	strcpy( input_filename , argv[1] );
	strcpy( output_filename , argv[2] );
	max_input_reads_in_a_single_nucl_loc = strtol ( argv[3] , &temp , 10 );
	/********************************************************************/

	//printf ( "\n max_input_reads_in_a_single_nucl_loc = %d\n" , max_input_reads_in_a_single_nucl_loc );
	//fflush ( stdout );
	findContinousClusters ( input_filename , output_filename , max_input_reads_in_a_single_nucl_loc );
	return 0;
}
