#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>

float HUMAN_THRESH = 1.0/3.0;  /*lower threshold for a contig to be considered non-human*/
int RUNS = 100; /*maximum number of rounds of extending a contig in one direction*/
int SUBLEN = 24; /*kmer size*/
int MAX_SIZE = 2560000000;
//int NUMEL = 0;
//float TOTSC = 0;

/*Takes as input string s, which is trusted to be of length <= len
 *allocates a new string of size 2* len +1 and copies s to the beginning of
 *the new string*/
char *double_string_right(char *s, int len)
{
    char *s2 = malloc(2 * len + 1);
    strncpy(s2, s, len+1); //Assumes accuracy of len
	//s2[len] = '\0';
    return s2;
}

/*structure to return a string and some ancillary information
  c is the string; len is the length without the null terminator; numel is number of scores summed; totsc is the total score*/
struct ret {
    char *c;
    int len;
    int flen;
    float totbac;
    float totvir;
};

/*Implements != with wildcard N*/
bool neq(char a, char b) {
	return a != b && a != 'N' && b != 'N';
}
/*Naive strstr with single-character wildcard N*/
//https://codereview.stackexchange.com/questions/35396/strstr-implementation
char *strstrn(char *haystack, char *needle) {
    const char *a = haystack, *b = needle;
    for (;;) {
        if      (!*b)          return haystack;
        else if (!*a)          return NULL;
        else if (neq(*a++, *b++)) { //gratuitous pointer math
			a = ++haystack; b = needle;
		}
    }
}

/*find_sub_right tries to extend a contig on the right
 * reads is a two-dimensional array of all reads
 * sb0 is the end piece to extend
 * status is the current read status array (see below)
 * num_reads is the length of reads
 * kmer_contains_n is a flag for whether sb0 contains an N
 * num_reads is the number of reads, equivalently the size of the major dimension of **reads
 * Returns the read with the earliest occurence of sb0 in one exists, specifically
 *   the index of the read, the index of sb0 in the read, and the location of an N in the read
 * Returns by setting the first three values of int array cmin
 * Marks all reads containing sb0 as used
 */
void find_sub_right(char *reads[], char sb0[], int *status, int num_reads, bool kmer_contains_n, int *cmin) {
	int pos;
	char *ptr;
	//static int cmin[3];
	cmin[0] = 100;
	cmin[1] = -1;
	cmin[2] = -1;
	//printf(" starting values: [%d,%d,%d] ", cmin[0], cmin[1], cmin[2]);
	//printf("kmer: %s\n", sb0);
	for (int i=0; i < num_reads; i++) {
		if (status[i]) {
			ptr = (kmer_contains_n || status[i] > 0) ? strstrn(reads[i], sb0) : strstr(reads[i], sb0);
			if (ptr) {
			    pos = ptr-reads[i];
			    if (pos < cmin[0] && pos >= 0) { //>=0 check seems unneeded but harmless
			        cmin[0] = pos;
			        cmin[1] = i;
			        cmin[2] = status[i]-1;
			    }
			    status[i] = 0;
			}
		}
	}
}

/*find_sub_left  tries to extend a contig on the left
 *Same idea as find_sub_self, except finds the read with the first occurence of sb0
 */
void find_sub_left(char *reads[], char sb0[], int *status, int num_reads, bool kmer_contains_n, int *cmax) {
	int pos;
	char *ptr;
	//static int cmax[3] = {-1,-1,-1};
	cmax[0] = -1;
	cmax[1] = -1;
	cmax[2] = -1;
	for (int i=0; i < num_reads; i++) {
		if (status[i]) {
			ptr = (kmer_contains_n || status[i] > 0) ? strstrn(reads[i], sb0) : strstr(reads[i], sb0);
			if (ptr) {
			  	pos = ptr-reads[i];
			  	if (pos > cmax[0] && pos>=0) {
					cmax[0] = pos;
					cmax[1] = i;
					cmax[2] = status[i]-1;
				}
				status[i] = 0;
			}
		}
	}
}

/*
 *assemble_right builds a contig starting from the seed read and extending to the right; read_list is the two-dimensional array of reads
 *bact_scores and viral_scores are the model predictions of each read; status is the current read status array (see below) 
 *num_reads is the number of reads or the size of upper dimension of read_list, SEGMENT_LENGTH is the read length
 */
//char* assemble_read(float *bact_scores, float *viral_scores, char *ch_arr[], int *status, float SEGMENT_LENGTH, int num_reads, int i){ [old signature]
struct ret assemble_right(int seed, char *read_list[], float *bact_scores, float* viral_scores, int *status, int SEGMENT_LENGTH, int num_reads){

	float TOTBAC = bact_scores[seed] * SEGMENT_LENGTH; //This stores the sum of the score of each base
	float TOTVIR = viral_scores[seed] * SEGMENT_LENGTH;
	char *read = read_list[seed];
	int clen = SEGMENT_LENGTH;
	int ccap = 2*SEGMENT_LENGTH; /*current upper bound on the length of the contig array*/
	char *contig = double_string_right(read, SEGMENT_LENGTH);

	char sb0[SUBLEN+1];
	strncpy(sb0,&read[SEGMENT_LENGTH-SUBLEN], SUBLEN);
	sb0[SUBLEN] = '\0';   //need to add a termination character because strncpy does not
	int cnt = 0;
	int *poss = malloc(3 * sizeof(int));
	int last_n = status[seed]-1; //Key note: this is an offset from the SEGMENT_LENGTH-th last position, -1 = no n in last S_L bases
	status[seed] = 0;
	int new_read_offset, amt_added;
	float  total_score = TOTBAC+TOTVIR;
	while (cnt < RUNS) { //keep trying to extend the contig as long as there is a possible extension and average non-human score per element stays above 0.5
		find_sub_right(read_list, sb0, status, num_reads, (last_n >= SEGMENT_LENGTH-SUBLEN), poss); //returns pos in read, index of read
		//Excludes 2 cases: poss[0] == S_L-SUB: best match has kmer at the far right. read not extended at all
		//or poss[0] == 100 (default value): no read found
		if (poss[0] >= SEGMENT_LENGTH-SUBLEN) { 
			break;
		}

		//printf("right extension w/ %d (new status: %d) ", poss[1], status[poss[1]]);
		new_read_offset = poss[0] + SUBLEN;
		amt_added = SEGMENT_LENGTH - new_read_offset;
		total_score += (bact_scores[poss[1]] + viral_scores[poss[1]]) * amt_added;
		if (total_score / (clen + amt_added) <= HUMAN_THRESH) {
			break;
		}

		//If the length of the contig is getting close to ccap, then double the available length and free the memory for the shorter array
		if (clen + amt_added +1 > ccap) {
			ccap = ccap * 2;
			char *new_contig = double_string_right(contig, ccap);
			free(contig);
			contig = new_contig;
		}
		//clen += SEGMENT_LENGTH;  //Why was this segment_length?

		clen += amt_added;
		last_n -= amt_added; //now an offset in new contig
		//Append the rest of the contig
		strcat(contig, &read_list[poss[1]][new_read_offset]); //only copies  part after kmer

		if (last_n != poss[2]) { //1) if N's are in same spot or 2) neither contig nor new read have N's, nothing to do
			if (last_n >= 0 && last_n < new_read_offset) {
				contig[clen-SEGMENT_LENGTH+last_n] = read_list[poss[1]][last_n]; //Overwrite N in contig with nt from read
			}
			last_n = poss[2] >= new_read_offset ? poss[2] : -1; //Save new N location if it was in the new part of the read
		}

		strncpy(sb0, &contig[clen-SUBLEN], SUBLEN);
		sb0[SUBLEN] = '\0'; //need to add a termination character because strncpy does not
		TOTBAC += bact_scores[poss[1]] * amt_added;
		TOTVIR += viral_scores[poss[1]] * amt_added;
		cnt++;	//update number of loop iterations
	}
	free(poss);
	struct ret rres = { contig, clen, clen, TOTBAC, TOTVIR };
	return rres;
}

/*assemble_left builds a contig starting from read and extending to the right;
 *arguments are the same as assemble_right, with three numbers describing the intermedate results from assemble_right
 */
struct ret assemble_left(int seed, char *read_list[], float *bact_scores, float* viral_scores, int *status, int SEGMENT_LENGTH, int num_reads, int FULLLEN, float TOTBAC, float TOTVIR){

	int clen = 0;
	//int ccap = SEGMENT_LENGTH;
	char* contig = calloc(1, sizeof(char));

	char sb0[SUBLEN+1];
	strncpy(sb0, read_list[seed], SUBLEN);
	sb0[SUBLEN] = '\0'; //need to add a termination character because strncpy does not
	int cnt = 0;
	int *poss = malloc(3*sizeof(int));
	int first_n = -1; //Offset from start; can't have an N in an empty contig
	int amt_added; //== new_read_offset
	float total_score = TOTBAC+TOTVIR;
	while (cnt<RUNS)  //keep trying to extend the contig as long as there is a possible extension and average score per element stays above 0.5
	{
		find_sub_left(read_list, sb0, status, num_reads, (first_n >= 0 && first_n < SUBLEN), poss);
		//Excludes 2 cases: poss[0] == 0: best match has kmer at the far left. read not extended at all
		//or poss[0] == -1 (default value): no read found		
		if (poss[0] <= 0){ 
			break;
		}
		//printf("left extension w/ %d (new status: %d) ", poss[1], status[poss[1]]);
		//new_read_offset = poss[0];
		amt_added = poss[0];
		total_score += (bact_scores[poss[1]] + viral_scores[poss[1]]) * amt_added;
		if (total_score / (FULLLEN + amt_added) <= HUMAN_THRESH) {
			break;
		}
		//Allocate new space for the contig because there is no easy way to prepend the added read on the left
		char *new_contig = malloc((amt_added + clen + 1) * sizeof(char));
		//Copy the new piece on the left
		strncpy(new_contig, read_list[poss[1]], amt_added);

		//new_contig[amt_added] = '\0'; // add termination character because strncpy does not

		clen += amt_added;
		FULLLEN += amt_added;
		first_n += amt_added;
		if (first_n >= SEGMENT_LENGTH) {
			first_n = -1; //Nothing to be done about that N
		}
		//Copy the old contig on the right
		strcpy(new_contig+amt_added, contig);

		free(contig);
		contig = new_contig;

		if (first_n != poss[2]) { //1) if N's are in same spot or 2) neither contig nor new read have N's, nothing to do
			if (amt_added <= first_n) {
				contig[first_n] = read_list[poss[1]][first_n]; //Overwrite N in contig with nt from read
			}
			first_n = (poss[2] < amt_added) ? poss[2] : -1; //Save new N location if it was in the new part of the read (could be -1)
		}
		
		strncpy(sb0, contig, SUBLEN);
		if (clen < SUBLEN) { //If we've added less than SUBLEN bases so far, we need some of the original read
			strncpy(sb0+clen, read_list[seed], SUBLEN-clen);
		}
		sb0[SUBLEN] = '\0'; // add termination character because strncpy does not
		//TOTSC+=poss[1]; //update  total score
		TOTBAC += bact_scores[poss[1]] * amt_added;
		TOTVIR += viral_scores[poss[1]] * amt_added;
		cnt++;	//update number of loop iterations

	}
	free(poss);
	struct ret lres = {contig, clen, FULLLEN, TOTBAC, TOTVIR };
	return lres;
}


/*
 *Done by python: identify seeds and determine sorted order,
 */

//                     bact scores,        viral scores,        reads,          sorted list of indicies of reads to used as seeds
int assemble_read_loop(float *bact_scores, float *viral_scores, char *ch_arr[], int *seed_indices, int SEGMENT_LENGTH, int num_reads, int num_seeds, char *fname){
	/*Status array tracks whether a read is used (0) or not, and also the nosition of an N in an unused read:
	 * -1 = unused, no N; value x between 1 and SEGMENT_LENGTH = N at position x (1-indexed)
	 */
	int* status = malloc(num_reads * sizeof(int));
	int i, j, seed, count=0, count2=0;
	char *contig, *ptr;

	for (i = 0; i < num_reads; i++) { //fill in status
		ptr = strchr(ch_arr[i], 'N');
		status[i] = ptr ? ptr - ch_arr[i] + 1 : -1;
	}

	FILE *fp;
	fp = fopen(fname, "w+");
	fprintf(fp, "#Recieved %d reads, of which %d are seeds.\n", num_reads, num_seeds);
	for(j = 0; j < num_seeds ; j++) {
		seed = seed_indices[j];
		if (status[seed]) { //includes N-reads. If there's an N in the seed covered by both left & right reads, the right read will determine the nt Noam: Not sure I understand
			//status[seed] = 0;
			struct ret rres = assemble_right(seed, ch_arr, bact_scores, viral_scores, status, SEGMENT_LENGTH, num_reads);
			//The current TOTSC & NUMEL values are passed as fields in rres to match previous behavior
			struct ret lres = assemble_left(seed, ch_arr, bact_scores, viral_scores, status, SEGMENT_LENGTH, num_reads, rres.flen, rres.totbac, rres.totvir);
			//fprintf(stdout,"Got left contig %s with expected length %d and actual length %d\n",lres.c, lres.len, strlen(lres.c));
			//fprintf(stdout,"Got right contig %s with expected length %d and actual length %d\n",rres.c, rres.len, strlen(rres.c));
			contig = malloc((lres.len + rres.len + 1) * sizeof(char));
			strcpy(contig, lres.c);
			strcpy(contig + lres.len, rres.c);
			free(lres.c);
			free(rres.c);
			float bac_score = (lres.totbac)/(lres.flen); 
			float vir_score = (lres.totvir)/(lres.flen);
			char *contig_n = strchr(contig, 'N'); //Not the most efficient, but want to check so we can use strstr if there is no N
			//fprintf(stdout,"About to write %s with expected length %d and actual length %d\n",contig, lres.len + rres.len, strlen(contig));

			//remove reads that occur in the contig, even if they weren't detected in extension
		   	for(i = 0; i < num_reads; i++) {
				if (status[i]){
					ptr = (contig_n || status[i] > 0) ? strstrn(contig, ch_arr[i]) : strstr(contig, ch_arr[i]);
					if (ptr) {
						status[i] = 0;
					}
				}
			}
			fprintf(fp,">contig_%d_bact:%f_virus:%f\n", seed, bac_score, vir_score);
			fprintf(fp, contig);
			fprintf(fp,"\n");
			count ++;
		} else { //if the seed was skipped b/c already marked as used
			count2 ++;
		}


	}
	fprintf(fp, "#Done. Wrote %d contigs and skipped %d used seeds\n", count, count2);
	fclose(fp);
	return 0;
}

