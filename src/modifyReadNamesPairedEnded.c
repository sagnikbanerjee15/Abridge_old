# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include <math.h>
# include <stdbool.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

void writeToFile (char **split_line, FILE *fhw)
{
	fprintf (fhw , "%s" , split_line[0]);
	fprintf (fhw , "%s" , "\t");
	fprintf (fhw , "%s" , split_line[1]);
}

void generateNextReadID (char *alphabets, int *read_id, int *read_length)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	int i;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/

	/********************************************************************/

	if ( *read_length == 0 )
	{
		read_id[0] = 0;
		( *read_length )++;
	}
	else
	{
		/*
		 * Check if the read is the last element of the maximum read_length
		 */
		for ( i = 0 ; i < *read_length ; i++ )
			if ( read_id[i] != strlen (alphabets) - 1 ) break;
		if ( i == *read_length )
		{
			( *read_length )++;
			for ( i = 0 ; i < *read_length ; i++ )
				read_id[i] = 0;
		}
		else
		{
			/*
			 * Increment the read_id
			 */
			for ( i = *read_length - 1 ; i >= 0 ; i-- )
			{
				if ( read_id[i] == ( strlen (alphabets) - 1 ) )
					read_id[i] = 0;
				else
				{
					read_id[i]++;
					break;
				}
			}
		}
	}
}

void convertReadIdToString (int *read_id, char *read_id_string, int read_length, char *alphabets)
{
	int i;
	for ( i = 0 ; i < read_length ; i++ )
		read_id_string[i] = alphabets[read_id[i]];
	read_id_string[i] = '\0';
}

void splitMappingInTwoPartsAndSetNHValue (char *line, char **split_line, int *NH_value)
{
	int i, j0, j1, first_tab_found = 0, k, NH_string_index;
	char NH_string[100];
	char *temp;

	j0 = 0;
	j1 = 0;
	NH_string_index = 0;
	for ( i = 0 ; line[i] != '\0' ; i++ )
	{
		if ( first_tab_found )
			split_line[1][j1++ ] = line[i];
		else split_line[0][j0++ ] = line[i];
		if ( line[i] == '\t' && first_tab_found == 0 ) first_tab_found = 1;
		if ( line[i - 1] == ':' && line[i - 2] == 'i' && line[i - 3] == ':' && line[i - 4] == 'H' && line[i - 5] == 'N' )
		{
			for ( k = i ; line[k] != '\t' ; k++ )
				NH_string[NH_string_index++ ] = line[k];
			NH_string[NH_string_index++ ] = '\0';
		}
	}
	( *NH_value ) = strtol (NH_string , &temp , 10);
	split_line[0][j0] = '\0';
	split_line[1][j1] = '\0';
}

struct Old_Read_ID_to_New_Read_ID_Circular_Linked_list* updateNodeInCircularLinkedList (struct Old_Read_ID_to_New_Read_ID_Circular_Linked_list *head, char *old_read_id, int *number_of_invalid_nodes)
{
	/*
	 * Scan the entire Circular linked list
	 * to find the old_read_id
	 * If old_read_id is found then reduce the value of number_of_multi_maps field by 1
	 * If number_of_multi_maps becomes zero, delete the node
	 */
	struct Old_Read_ID_to_New_Read_ID_Circular_Linked_list *temp = NULL;
	if ( head == NULL ) return temp;
	temp = head;
	while ( temp->next != head )
	{
		if ( strcmp (temp->old_read_id , old_read_id) == 0 )
		{
			temp->number_of_multi_maps--;
			if ( temp->number_of_multi_maps == 0 ) //Remove this node and add it before head - Do not Delete it
			{
				temp->prev = head->prev;
				head->prev = temp;
				temp->next = head;
				temp->prev->next = temp;
				temp->valid = 0;
				( *number_of_invalid_nodes )++;
			}
			return temp;
		}
		temp = temp->next;
	}
	return temp;
}

void insertNodeInCircularLinkedList (struct Old_Read_ID_to_New_Read_ID_Circular_Linked_list *head, char *old_read_id, char *new_read_id, int NH_value)
{
	if ( head == NULL )
	{
		head = allocateMemoryOld_Read_ID_to_New_Read_ID_Circular_Linked_list ();
		head->prev = head;
		strcpy(head->new_read_id , new_read_id);
		strcpy(head->old_read_id , old_read_id);
		head->valid = 1;
		head->number_of_multi_maps = NH_value * 2 - 1;
		head->next = head;
	}
	else
	{
		struct Old_Read_ID_to_New_Read_ID_Circular_Linked_list *temp;
		struct Old_Read_ID_to_New_Read_ID_Circular_Linked_list *node;
		temp = head;
		if ( temp->prev == head ) //Only one node in linked list
			node = allocateMemoryOld_Read_ID_to_New_Read_ID_Circular_Linked_list ();
		else
		{
			while ( temp->prev->prev != head && temp->prev->valid == 0 )
				temp = temp->prev;
		}
		if ( temp == head )
			node = allocateMemoryOld_Read_ID_to_New_Read_ID_Circular_Linked_list ();

		temp = node;
		temp->next = allocateMemoryOld_Read_ID_to_New_Read_ID_Circular_Linked_list ();
		temp->next->next = head;
		temp->next->prev = temp;
		head->prev = temp->next;
		temp = temp->next;
		strcpy(temp->new_read_id , new_read_id);
		strcpy(temp->old_read_id , old_read_id);
		temp->number_of_multi_maps = NH_value * 2 - 1;
		temp->valid = 1;
	}
}

void deleteEntireLinkedList (struct Old_Read_ID_to_New_Read_ID_Circular_Linked_list *head)
{
	if ( head == NULL ) return;
	struct Old_Read_ID_to_New_Read_ID_Circular_Linked_list *temp;
	struct Old_Read_ID_to_New_Read_ID_Circular_Linked_list *node_to_be_removed;
	temp = head;
	while ( temp->next != NULL )
	{
		node_to_be_removed = temp->next;
		temp->next->next->prev = temp;
		temp->next = temp->next->next;
		node_to_be_removed->prev = NULL;
		node_to_be_removed->next = NULL;
		free (node_to_be_removed->new_read_id);
		free (node_to_be_removed->old_read_id);
		free (node_to_be_removed);
	}

	node_to_be_removed = head;
	node_to_be_removed->prev = NULL;
	node_to_be_removed->next = NULL;
	free (node_to_be_removed->new_read_id);
	free (node_to_be_removed->old_read_id);
	free (node_to_be_removed);
}

void deleteInvalidNodesFromCircularLinkedList (struct Old_Read_ID_to_New_Read_ID_Circular_Linked_list *head, int MAX_number_of_invalid_nodes_allowed)
{
	struct Old_Read_ID_to_New_Read_ID_Circular_Linked_list *temp;
	struct Old_Read_ID_to_New_Read_ID_Circular_Linked_list *node_to_be_removed;

	while ( MAX_number_of_invalid_nodes_allowed-- )
	{
		temp = head->prev;
		head->prev = temp->prev;
		head->prev->next = head;
		temp->prev = NULL;
		temp->next = NULL;
		free (temp->old_read_id);
		free (temp->new_read_id);
		free (temp);
	}
}

void convertOldReadIdsToNewReadIds (char *input_samfilename, char *output_samfilename)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char alphabets[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-='{}[]|?<>,.";
	char *temp; //Required for strtoi
	char *line = NULL; // for reading each line
	char **split_line; // List of strings to store each element of a single alignment
	char prev_old_read_id[1000];
	char read_id_string[100];

	int read_length;
	int i;
	int j;
	int read_id[100];
	int NH_value;
	int number_of_invalid_nodes = 0;
	int MAX_number_of_invalid_nodes_allowed = 10000;

	long long int read_number;

	size_t len = 0;
	ssize_t line_len;

	FILE *fhr;
	FILE *fhw;

	struct Old_Read_ID_to_New_Read_ID_Circular_Linked_list *read_mapping_head = NULL;
	struct Old_Read_ID_to_New_Read_ID_Circular_Linked_list *locate_node_of_interest;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr = fopen (input_samfilename , "r");
	if ( fhr == NULL )
	{
		printf ("Error! File %s not found" , input_samfilename);
		exit (1);
	}
	fhw = fopen (output_samfilename , "w");
	if ( fhw == NULL )
	{
		printf ("Error! File %s cannot be created" , output_samfilename);
		exit (1);
	}

	split_line = ( char** ) malloc (sizeof(char*) * 100);
	for ( i = 0 ; i < 100 ; i++ )
		split_line[i] = ( char* ) malloc (sizeof(char) * 1000);

	/********************************************************************/

	while ( ( line_len = getline ( &line , &len , fhr) ) != -1 )
	{
		if ( line[0] == '@' )
			fprintf (fhw , "%s" , line);
		else break;
	}

	read_number = 0;
	prev_old_read_id[0] = '\0';
	do
	{
		read_number++;
		splitMappingInTwoPartsAndSetNHValue (line , split_line , &NH_value);
		if ( read_number % 2 == 0 && strcmp (prev_old_read_id , split_line[0]) == 0 ) //Read id is exactly same as previous line no need to update or consult dictionary
		{
			strcpy(split_line[0] , read_id_string);
			updateNodeInCircularLinkedList (read_mapping_head , split_line[0] , &number_of_invalid_nodes);
		}
		else
		{
			strcpy(prev_old_read_id , split_line[0]);
			locate_node_of_interest = updateNodeInCircularLinkedList (read_mapping_head , split_line[0] , &number_of_invalid_nodes);
			if ( locate_node_of_interest != NULL )
			{
				strcpy(split_line[0] , locate_node_of_interest->new_read_id);
			}
			else
			{
				generateNextReadID (alphabets , read_id , &read_length);
				convertReadIdToString (read_id , read_id_string , read_length , alphabets);
				insertNodeInCircularLinkedList (read_mapping_head , split_line[0] , read_id_string , NH_value);
			}
		}
		writeToFile (split_line , fhw);
		if ( number_of_invalid_nodes > MAX_number_of_invalid_nodes_allowed )
			deleteInvalidNodesFromCircularLinkedList (read_mapping_head , MAX_number_of_invalid_nodes_allowed);
	} while ( ( line_len = getline ( &line , &len , fhr) ) != -1 );

	deleteEntireLinkedList (read_mapping_head);
	fclose (fhr);
	fclose (fhw);
}

int main (int argc, char *argv[])
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	int read_length;
	int read_id[100];
	int i, j;

	char alphabets[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-='{}[]|?<>,.";
	char input_samfilename[FILENAME_LENGTH];
	char output_samfilename[FILENAME_LENGTH];

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	read_length = 0;
	strcpy(input_samfilename , argv[1]);
	strcpy(output_samfilename , argv[2]);

	/********************************************************************/
	convertOldReadIdsToNewReadIds (input_samfilename , output_samfilename);
	return 0;
}
