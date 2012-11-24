#include "parallel_reconstruct.h"

void p_reconstruct(char *infilename, int iterations)
{
	setbuf(stdout, NULL);
	int i,j;
	int nx, ny;
	int nx_min, nx_max, nx_range;
	int rank=0, nproc=1;

	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
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
	MPI_Datatype contiguous;
	MPI_Type_contiguous(ny, MPI_FLOAT, &contiguous);
	MPI_Type_commit(&contiguous);
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
	for(i=0; i<iterations;i++)
	{
		MPI_Request request[2];
		MPI_Status status[2];
		//if(rank == 0) printf("on iteration %d\n", i);
		reconstruct_image_segment(segment, edge, nx_min, nx_max, ny); 
		MPI_Issend(&segment[1][1],ny,MPI_FLOAT,left,0,cartesian_communicator, &request[0]);
		MPI_Issend(&segment[nx_range][1],ny,MPI_FLOAT,right,1,cartesian_communicator, &request[1]);

		int recv = 0;
		while(recv == 0)
		{
			MPI_Irecv(&segment[nx_range+1][1],ny,MPI_FLOAT,right,0,cartesian_communicator, &status[0]);
			MPI_Irecv(&segment[0][1],ny,MPI_FLOAT,left,1,cartesian_communicator, &status[1]);
			MPI_Testall(2, request, &recv, status);
		}
	}
	removehalo(segment, out, nx_range, ny);
	MPI_Gather(out, ny*nx_range, MPI_FLOAT, buf, ny*nx_range, MPI_FLOAT, 0, MPI_COMM_WORLD);
	if(rank == 0)
	{
		char outfilename[256];
		sprintf(outfilename, "../output_files/p%d_reconstruct_%dx%d_%di.pgm", nproc, nx, ny, iterations);
		pgmwrite(outfilename, buf, nx, ny);
	}
	MPI_Finalize();
}
