#include "reconstruct.h"

void reconstruct(char *infilename, int iterations)
{
	int i;
	int nx, ny;
	char outfilename[256];
	float max_change = 0;
	pgmsize(infilename, &nx, &ny);
	
	float buf[nx][ny];
	float segment[nx+2][ny+2];
	pgmread(infilename, buf, nx, ny);
	initialise_segment(segment, nx+2, ny+2);
	for(i=0; i<iterations;i++)
	{
		reconstruct_image_segment(segment, buf, 0, nx, 0, ny, &max_change);
	}
	removehalo(segment, buf, nx, ny);
	sprintf(outfilename, "../output_files/reconstruct_%dx%d_%d-iterations.pgm", nx, ny, iterations);
	pgmwrite(outfilename, buf, nx, ny);
}

void initialise_segment(void *segment_in, int nx, int ny)
{
	int i,j;
	float *segment = (float *) segment_in;
	for(j=0;j<ny;j++)
	{
		for(i=0;i<nx;i++)
		{
			segment[i*ny+j] = 255.0;
		}
	}	
}

void reconstruct_image_segment(void *segment, void *edge_in, int nx_min, int nx_max, int ny_min, int ny_max, float *max_change)
{
	int i,j;
	int nx_range = nx_max-nx_min;
	int ny_range = ny_max-ny_min;
	float *old = (float *) segment;
	float new[nx_range+2][ny_range+2];
	float *edge = (float *) edge_in;
	*max_change = 0;
	for(j=1;j<ny_range+1;j++)
	{
		for(i=1;i<nx_range+1;i++)
		{
			new[i][j] = 0.25*(old[(i-1)*(ny_range+2)+j]+old[(i+1)*(ny_range+2)+j]+old[i*(ny_range+2)+(j-1)]+old[i*(ny_range+2)+(j+1)]-edge[(i-1)*ny_range+j-1]);
		}
	}	
	for(j=1;j<ny_range+1;j++)
	{
		for(i=1;i<nx_range+1;i++)
		{
			if(fabs(old[i*(ny_range+2)+j]-new[i][j]) > *max_change) *max_change = fabs(old[i*(ny_range+2)+j]-new[i][j]);
			old[i*(ny_range+2)+j] = new[i][j];
		}
	}	
}

void removehalo(void *image_in, void *buf_out, int nx, int ny)
{
	int i,j;
	float *buf = (float *) buf_out;
	float *image = (float *) image_in;
	for(j=0;j<ny;j++)
	{
		for(i=0;i<nx;i++)
		{
			buf[i*ny+j] = image[(i+1)*(ny+2)+(j+1)];
		}
	}
}
