#include "edgeview.h"

/*	Reads a pgm file and immediately writes it */
void edgeview(char *infilename, int M, int N)
{
	double image[M][N];
	char outfilename[256];
	sprintf(outfilename, "../output_files/edgeview_%dx%d.pgm", M, N);
	pgmread(infilename, image, M, N);
	pgmwrite(outfilename, image, M, N);
	return;
}
