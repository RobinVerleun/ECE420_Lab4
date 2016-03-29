
#define LAB4_EXTEND

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "../util/timer.h"
#include "../util/Lab4_IO.h"

#define EPSILON 0.00001
#define DAMPING_FACTOR 0.85
#define PRIMARY 0

int main( int argc, char** argv) {

	//Initialization data
	struct node *nodehead;
	int nodecount, *in_links, *out_links;
	int local_start, local_end, local_chunk;
	double damp_const, *r, *r_pre, *local_work;
	
	//Timing data
	double startTime, endTime;
	
	//Iterators
	int i, j;

	//MPI variables
	int threads, rank;
	MPI_Comm comm;

	MPI_Init(&argc, &argv);

	//Get MPI meta_data
	comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &threads);
	MPI_Comm_rank(comm, &rank);


	// ================== Initialize and Load Matrix ======================
	// Load the inputs using Lab4_IO.c method.
	get_node_stat(&nodecount, &in_links, &out_links);
	node_init(&nodehead, in_links, out_links, 0, nodecount);
	
	
	// Allocate space for r and r_pre 
	r = malloc(nodecount * sizeof(double));
	
	// Assign initial values to r - formula 1 from lab handout. 
	for(i = 0; i < nodecount; ++i){
		r[i] = 1.0 / nodecount;
	}
	//if(rank == PRIMARY){
		r_pre = malloc(nodecount * sizeof(double));
	//}

	// Calculate the dampening constant.
	damp_const = (1.0 - DAMPING_FACTOR) / nodecount;
	
	// Local values for each process
	local_chunk = nodecount/threads;
	local_start = local_chunk * rank;
	local_end = local_start + local_chunk;
	
	// Space for local copy of algorithm
	local_work = malloc(local_chunk*(sizeof(double)));
	
	// ========================== Page Rank Work ==========================
	
	GET_TIME(startTime);
	//Error to compare to epsilon
	double error = 100;

	do {
		//if(rank == PRIMARY){
			vec_cp(r, r_pre, nodecount);
		//}
		for( i = local_start; i < local_end; ++i ) { 
			local_work[i-local_start] = 0;
			for ( j = 0; j < nodehead[i].num_in_links; ++j ) {
				local_work[i-local_start] += r_pre[nodehead[i].inlinks[j]] / out_links[nodehead[i].inlinks[j]];
			}
			local_work[i-local_start] *= DAMPING_FACTOR;
			local_work[i-local_start] += damp_const;
		}
		
		MPI_Allgather(local_work, local_chunk, MPI_DOUBLE, r, local_chunk, MPI_DOUBLE, comm);
		
		if(rank == PRIMARY){
			error = rel_error(r, r_pre, nodecount);
		}
		
		MPI_Bcast(&error, 1, MPI_DOUBLE, PRIMARY, comm);
	}while(error >= EPSILON);

	GET_TIME(endTime);
	
	if(rank == PRIMARY){
		Lab4_saveoutput(r, nodecount, endTime-startTime);
		printf("Time elapsed: %f\n", endTime-startTime);
	}
	
	//Free Data
	node_destroy(nodehead, nodecount);
	free(local_work);
	free(r);
	free(r_pre); 
	MPI_Finalize();

	return 0;
}


