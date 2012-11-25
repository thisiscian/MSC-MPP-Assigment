#include "parallel_reconstruct.h"

void p_reconstruct(char *infilename, int type, float *limits)
{
	setbuf(stdout, NULL);
	int rank=0, nproc=1;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	float change[nproc], max_change;
	int nx, ny, nx_min, nx_max, nx_range;
	int i=0, j=0, keep_going=1;
	int cartesian_communicator;
	int ndims = 1;
	int dims[ndims];
	dims[0] = 0;
	int period[ndims];
	int left;
	int right;
	period[0] = 0;
	MPI_Dims_create(nproc, ndims, dims);
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, period, 0, &cartesian_communicator);
	MPI_Cart_shift(cartesian_communicator, 0, 1,  &left, &right);
	pgmsize(infilename, &nx, &ny);
	nx_min = (int)((double)rank/(double)nproc*nx);
	nx_max = (int)((double)(rank+1)/(double)nproc*nx);
	if(nx_max > nx) nx_max = nx;
	nx_range = nx_max-nx_min;

	float buf[nx][ny];
	float edge[nx_range][ny];
	float out[nx_range][ny];
	float segment[nx_range+2][ny+2];
	if(rank == 0)
	{	
		pgmread(infilename, buf, nx, ny);
	}
	MPI_Scatter(buf, ny*nx_range, MPI_FLOAT, edge, ny*nx_range, MPI_FLOAT, 0, MPI_COMM_WORLD);
	initialise_segment(segment, nx_range+2, ny+2);
	while(keep_going == 1)
	{
		MPI_Request request[2];
		MPI_Status status[2];
		reconstruct_image_segment(segment, edge, nx_min, nx_max, ny, &max_change); 
		MPI_Issend(&segment[1][1],ny,MPI_FLOAT,left,0,cartesian_communicator, &request[0]);
		MPI_Issend(&segment[nx_range][1],ny,MPI_FLOAT,right,1,cartesian_communicator, &request[1]);
		MPI_Recv(&segment[nx_range+1][1],ny,MPI_FLOAT,right,0,cartesian_communicator, &status[0]);
		MPI_Recv(&segment[0][1],ny,MPI_FLOAT,left,1,cartesian_communicator, &status[1]);
		i++;
		switch (type)
		{
			case 0:
				if(i >= *limits) keep_going = 0;
				break;
			case 1:
				MPI_Gather(&max_change, 1, MPI_FLOAT, change, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
				if(rank == 0)
				{
					for(j=1; j<nproc; j++)
					{
						if(max_change < change[j]) max_change = change[j];
					}
					if(max_change <= *limits) keep_going = 0;
					for(j=1;j<nproc;j++)
					{
						MPI_Issend(&keep_going, 1, MPI_INT, j, 0, MPI_COMM_WORLD, &request[0]);
					}
					printf("iteration=%d (max_change = %f)\n", i, max_change);
				}
				else
				{
					MPI_Recv(&keep_going, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status[0]);
				}
				break;
		}
	}
	removehalo(segment, out, nx_range, ny);
	MPI_Gather(out, ny*nx_range, MPI_FLOAT, buf, ny*nx_range, MPI_FLOAT, 0, MPI_COMM_WORLD);
	if(rank == 0)
	{
		char outfilename[256];
		switch (type)
		{
			case 0:
				sprintf(outfilename, "../output_files/p%d_reconstruct_%dx%d_%.0fi.pgm", nproc, nx, ny, *limits);
				break;
			case 1:
				sprintf(outfilename, "../output_files/p%d_reconstruct_%dx%d_%fl.pgm", nproc, nx, ny, *limits);
				break;
		}
		pgmwrite(outfilename, buf, nx, ny);
		printf("Program did %d iterations.\n", i);
	}
	MPI_Finalize();
}
