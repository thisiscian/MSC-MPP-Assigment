#ifndef __PARALLEL_RECONSTRUCT_H__
#define __PARALLEL_RECONSTRUCT_H__
#include "pgmio.h"
#include "stdio.h"
#include "reconstruct.h"
#include <mpi.h>

void p_reconstruct(char *inpath, char *infile, char *outpath, int limit_type, float *limits);
void swap_halos(int n[], int dims[], int n_range[], void *segment, int up, int down, int left, int right, int cart_comm);
void scatter_2D(void *buf, void *edge, int rank, int nproc, int cart_comm, int ndims, int n[], int n_range[], int dims[], MPI_Datatype segblock[]);
void gather_2D(void *buf_in, void *out_in, int rank, int nproc, int cart_comm, int ndims, int n[], int n_range[], int dims[]);
#endif
