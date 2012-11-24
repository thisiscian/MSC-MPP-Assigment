#ifndef __RECOnySTRUCT_H__
#define __RECOnySTRUCT_H__
#include "pgmio.h"
#include "stdio.h"

void reconstruct(char *filename, int iterations);
void removehalo(void *image, void *buf, int nx, int ny);
void initialise_segment(void *reconstruct_in, int nx, int ny);
void reconstruct_image_segment(void *segment, void *edge, int nx_min, int nx_max, int ny);
#endif
