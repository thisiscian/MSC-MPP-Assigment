#ifndef __RECONSTRUCT_H__
#define __RECONSTRUCT_H__
#include "pgmio.h"
#include "stdio.h"
#include <math.h>

void reconstruct(char *inpath, char *infile, char *outpath, int limit_type, float *limits);
void removehalo(void *image, void *buf, int nx, int ny);
void initialise_segment(void *reconstruct_in, int nx, int ny);
void reconstruct_image_segment(void *segment, void *edge, int nx_min, int nx_max, int ny_min, int ny_max, float *max_change, float *total);
#endif
