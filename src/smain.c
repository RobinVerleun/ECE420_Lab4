
#define LAB4_EXTEND

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../util/Lab4_IO.h"
#include "../util/timer.h"

#define EPSILON 0.00001
#define DAMPING_FACTOR 0.85


int main( int argc, char** argv) {

	//Operational data
	struct node *nodehead;
	int nodecount, *in_links, *out_links;
	double damp_const, *r, *r_pre;
	double start, end;
	int iteration = 0;
	int i, j;


	// ================== Initialize and Read In Matrix ======================
	get_node_stat(&nodecount, &in_links, &out_links);
	
	// Load the inputs using Lab4_IO.c method.
	get_node_stat(&nodecount, &in_links, &out_links);

	
	// Initialize the nodes
	node_init(&nodehead, in_links, out_links, 0, nodecount);

	
	// Allocate space for r and r_pre 
	r = malloc(nodecount * sizeof(double));
	r_pre = malloc(nodecount * sizeof(double));

	/*
	 * Assign initial values to r - formula 1 from lab handout.
	 * Calculate dampining constant.
	 * Calculate the next iteration probabilities until error gap is met.
	 */
	for(i = 0; i < nodecount; ++i) 
		r[i] = 1.0 / nodecount;

	damp_const = (1.0 - DAMPING_FACTOR) / nodecount;
	
	GET_TIME(start);

	do {
		++iteration;
		vec_cp(r, r_pre, nodecount);
		for( i = 0; i < nodecount; ++i ) {
			r[i] = 0;
			for ( j = 0; j < nodehead[i].num_in_links; ++j ) 
				r[i] += r_pre[nodehead[i].inlinks[j]] / out_links[nodehead[i].inlinks[j]];
			r[i] *= DAMPING_FACTOR;
			r[i] += damp_const;
		}
	}while(rel_error(r, r_pre, nodecount) >= EPSILON);

	GET_TIME(end);

	Lab4_saveoutput(r, nodecount, end-start);
	printf("Time elapsed: %f\n", end-start);

	/*
	 * Free data.
	 */
	node_destroy(nodehead, nodecount);
	free(r);
	free(r_pre); 

	return 0;
}
