#ifndef __PARALLEL_RECONSTRUCT_H__
#define __PARALLEL_RECONSTRUCT_H__
#include "pgmio.h"
#include "stdio.h"
#include "reconstruct.h"
#include <mpi.h>

void p_reconstruct(char *inpath, char *infile, char *outpath, float *limits);
#endif
