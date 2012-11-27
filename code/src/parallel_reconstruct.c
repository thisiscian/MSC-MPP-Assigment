#include "parallel_reconstruct.h"

void p_reconstruct(char *inpath, char *infile, char *outpath, float *limits)
{
	char outfilename[256], infilename[256];

	/* initialises MPI
 	 * sets 'rank' to the id of the current thread
 	 * sets 'nproc' to the total number of threads */
	int rank=0, nproc=1;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);


	/* integers that dicate sizes of the image
   * and sizes of the subsegments of the image */
	int n[2], n_min[2], n_max[2], n_range[2];
	int m_min[2], m_max[2], m_range[2];

	/* integers required for the various loops
   * if the limit is the number of iterations, 
   * and the limit is less than one, don't do an[1] iterations*/
	int i=0, j=0, keep_going=1;

	/* floats used when the limit is to a minimum change in the pixels */
	float change[nproc], max_change; // floats that 
	
	/* things required for setting the cartesian communicator */
	int cartesian_communicator;
	int ndims = 2;
	int dims[ndims], period[ndims], coords[ndims];
	int left, right, up, down;


	for(i=0;i<ndims; i++)
	{
		dims[i] = 0; 	 // automatically choose processers per dim
		period[i] = 0; // do not want periodicity
	}
	MPI_Dims_create(nproc, ndims, dims);
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, period, 0, &cartesian_communicator);
	MPI_Cart_shift(cartesian_communicator, 0, 1,  &left, &right);
	MPI_Cart_shift(cartesian_communicator, 1, 1,  &up, &down);
	MPI_Cart_coords(cartesian_communicator, rank, ndims, coords);
	
	/* this section sets the sizes of the file, according to the image */
	sprintf(infilename, "%s/%s", inpath, infile);
	pgmsize(infilename, &n[0], &n[1]);
	for(i=0; i<2; i++)
	{
		n_min[i] = (int)((double)coords[i]/(double)dims[i]*n[i]);
		n_max[i] = (int)((double)(coords[i]+1)/(double)dims[i]*n[i]);
		if(n_max[i] > n[i]) n_max[i] = n[i]; // 'n_max[0]' shouldn't be greater than 'n[0]'
		n_range[i] = n_max[i]-n_min[i];
	}
	/* the arrays that actually contain the data to do with the reconstruction */
	float buf[n[0]][n[1]];
	float edge[n_range[0]][n_range[1]];
	float out[n_range[0]][n_range[1]];
	float segment[n_range[0]+2][n_range[1]+2]; //+2 because it contains 1px halos

	/* make the master thread read in the file
 	 * and then gives the other threads the appropriate sections */
	if(rank == 0)
	{	
		pgmread(infilename, buf, n[0], n[1]);
	}

	/* things required for the datatype that splits the image up into segments */
	MPI_Datatype segblock_main, segblock_right, segblock_down, segblock_rightdown;
	int n_range_right = n[0]-(int)((double)(dims[0]-1)/(double)dims[0]*n[0]);
	int n_range_down = n[1]-(int)((double)(dims[1]-1)/(double)dims[1]*n[1]);
	
	MPI_Type_vector(n_range[0], n_range[1], n[1], MPI_FLOAT, &segblock_main);
	MPI_Type_commit(&segblock_main);

	MPI_Type_vector(n_range[0], n_range_down, n[1], MPI_FLOAT, &segblock_down);
	MPI_Type_commit(&segblock_down);

	MPI_Type_vector(n_range_right, n_range[1], n[1], MPI_FLOAT, &segblock_right);
	MPI_Type_commit(&segblock_right);

	MPI_Type_vector(n_range_right, n_range_down, n[1], MPI_FLOAT, &segblock_rightdown);
	MPI_Type_commit(&segblock_rightdown);
	
	/* required for MPI_Issend and MPI_Recv */
	MPI_Request request[4];
	MPI_Status status[4];

	/* sending out the data to the various segments	
   * possibly writing to unsafe memory locations...(with uneven segment sizes)*/
	if(rank ==0)
	{
		for(j=0; j<nproc; j++)
		{
			MPI_Cart_coords(cartesian_communicator, j, ndims, coords);
			for(i=0; i<ndims; i++)
			{
				m_min[i] = (int)((double)coords[i]/(double)dims[i]*n[i]);
				m_max[i] = (int)((double)(coords[i]+1)/(double)dims[i]*n[i]);
				if(m_max[i] > n[i]) m_max[i] = n[i]; // 'n_max[0]' shouldn't be greater than 'n[0]'
				m_range[i] = m_max[i]-m_min[i];
			}
			if(m_range[0] != n_range[0] && m_range[1] != n_range[1])
				MPI_Issend(&buf[m_min[0]][m_min[1]], 1, segblock_rightdown, j, 0, cartesian_communicator, &request[0]);
			else if(m_range[0] != n_range[0])
				MPI_Issend(&buf[m_min[0]][m_min[1]], 1, segblock_right, j, 0, cartesian_communicator, &request[0]);
			else if(m_range[1] != n_range[1])
				MPI_Issend(&buf[m_min[0]][m_min[1]], 1, segblock_down, j, 0, cartesian_communicator, &request[0]);
			else 
				MPI_Issend(&buf[m_min[0]][m_min[1]], 1, segblock_main, j, 0, cartesian_communicator, &request[0]);
		}
	}
	MPI_Recv(edge, n_range[0]*n_range[1], MPI_FLOAT, 0, 0, cartesian_communicator, &status[0]);
	
	/* sets all the values of the thread's 'segment' to 255, including the halo */
	initialise_segment(segment, n_range[0]+2, n_range[1]+2);

	/* needed to send noncontiguous data */
	MPI_Datatype noncontig, noncontig_right;
	MPI_Type_vector(n_range[0], 1, n_range[1]+2, MPI_FLOAT, &noncontig);
	MPI_Type_commit(&noncontig);
	MPI_Type_vector(n_range_right, 1, n_range[1]+2, MPI_FLOAT, &noncontig_right);
	MPI_Type_commit(&noncontig_right);
	
	/* reconstruct the image
   * send the halo data
   * then check if the limit has been achieved */
	i=0;
	while(keep_going == 1)
	{
		reconstruct_image_segment(segment, edge, n_min[0], n_max[0], n_min[1], n_max[1], &max_change);
		if(n[0]/dims[0] != n_range[0] && n[1]/dims[1] != n_range[1])
		{
			MPI_Issend(&segment[1][1]							,1						,noncontig_right,up		,0,cartesian_communicator, &request[0]);
			MPI_Issend(&segment[1][n_range_down]	,1						,noncontig_right,down	,1,cartesian_communicator, &request[1]);
			MPI_Issend(&segment[1][1]							,n_range_down	,MPI_FLOAT			,left	,2,cartesian_communicator, &request[2]);
			MPI_Issend(&segment[n_range_right][1]	,n_range_down	,MPI_FLOAT			,right,3,cartesian_communicator, &request[3]);
		
			MPI_Recv(&segment[1][n_range_down+1]	,1						,noncontig_right,down	,0,cartesian_communicator, &status[0]);
			MPI_Recv(&segment[1][0]								,1						,noncontig_right,up		,1,cartesian_communicator, &status[1]);
			MPI_Recv(&segment[n_range_right+1][1]	,n_range_down	,MPI_FLOAT			,right,2,cartesian_communicator, &status[2]);
			MPI_Recv(&segment[0][1]								,n_range_down	,MPI_FLOAT			,left	,3,cartesian_communicator, &status[3]);
		}
		else if(n[1]/dims[1] != n_range[1])
		{
			MPI_Issend(&segment[1][1]							,1						,noncontig			,up		,0,cartesian_communicator, &request[0]);
			MPI_Issend(&segment[1][n_range_down]	,1						,noncontig			,down	,1,cartesian_communicator, &request[1]);
			MPI_Issend(&segment[1][1]							,n_range_down	,MPI_FLOAT			,left	,2,cartesian_communicator, &request[2]);
			MPI_Issend(&segment[n_range[0]][1]		,n_range_down	,MPI_FLOAT			,right,3,cartesian_communicator, &request[3]);
			
			MPI_Recv(&segment[1][n_range_down+1]	,1						,noncontig			,down	,0,cartesian_communicator, &status[0]);
			MPI_Recv(&segment[1][0]								,1						,noncontig			,up		,1,cartesian_communicator, &status[1]);
			MPI_Recv(&segment[n_range[0]+1][1]		,n_range_down	,MPI_FLOAT			,right,2,cartesian_communicator, &status[2]);
			MPI_Recv(&segment[0][1]								,n_range_down	,MPI_FLOAT			,left	,3,cartesian_communicator, &status[3]);
		}
		else if(n[0]/dims[0] != n_range[0])
		{
			MPI_Issend(&segment[1][1]							,1						,noncontig_right,up		,0,cartesian_communicator, &request[0]);
			MPI_Issend(&segment[1][n_range[1]]		,1						,noncontig_right,down	,1,cartesian_communicator, &request[1]);
			MPI_Issend(&segment[1][1]							,n_range[1]		,MPI_FLOAT			,left	,2,cartesian_communicator, &request[2]);
			MPI_Issend(&segment[n_range_right][1]	,n_range[1]		,MPI_FLOAT			,right,3,cartesian_communicator, &request[3]);
			
			MPI_Recv(&segment[1][n_range[1]+1]		,1						,noncontig_right,down	,0,cartesian_communicator, &status[0]);
			MPI_Recv(&segment[1][0]								,1						,noncontig_right,up		,1,cartesian_communicator, &status[1]);
			MPI_Recv(&segment[n_range_right+1][1]	,n_range[1]		,MPI_FLOAT			,right,2,cartesian_communicator, &status[2]);
			MPI_Recv(&segment[0][1]								,n_range[1]		,MPI_FLOAT			,left	,3,cartesian_communicator, &status[3]);
		}
		else
		{
			MPI_Issend(&segment[1][1]							,1						,noncontig			,up		,0,cartesian_communicator, &request[0]);
			MPI_Issend(&segment[1][n_range[1]]		,1						,noncontig			,down	,1,cartesian_communicator, &request[1]);
			MPI_Issend(&segment[1][1]							,n_range[1]		,MPI_FLOAT			,left	,2,cartesian_communicator, &request[2]);
			MPI_Issend(&segment[n_range[0]][1]		,n_range[1]		,MPI_FLOAT			,right,3,cartesian_communicator, &request[3]);
			
			MPI_Recv(&segment[1][n_range[1]+1]		,1						,noncontig			,down	,0,cartesian_communicator, &status[0]);
			MPI_Recv(&segment[1][0]								,1						,noncontig			,up		,1,cartesian_communicator, &status[1]);
			MPI_Recv(&segment[n_range[0]+1][1]		,n_range[1]		,MPI_FLOAT			,right,2,cartesian_communicator, &status[2]);
			MPI_Recv(&segment[0][1]								,n_range[1]		,MPI_FLOAT			,left	,3,cartesian_communicator, &status[3]);
		}

		i++;
		if(i >= *limits && *limits >= 1) keep_going = 0;
		else if(*limits < 1)
		{
			MPI_Gather(&max_change, 1, MPI_FLOAT, change, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
			if(rank == 0)
			{
				for(j=1; j<nproc; j++)
				{
					if(max_change < change[j]) max_change = change[j];
				}
				if(max_change <= *limits ) keep_going = 0;
				for(j=1;j<nproc;j++)
				{
					MPI_Issend(&keep_going, 1, MPI_INT, j, 0, MPI_COMM_WORLD, &request[0]);
				}
			}
			else
			{
				MPI_Recv(&keep_going, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status[0]);
			}
			break;
		}
	}
	/* sets 'out' to the nonhalo region of 'segment'*/
	removehalo(segment, out, n_range[0], n_range[1]);

	/* have the master thread gather all 'out' arrays
   * 	and set it to the 'buf' array */
	MPI_Issend(out, n_range[0]*n_range[1], MPI_FLOAT, 0, rank, cartesian_communicator, &request[0]);
	if(rank ==0)
	{
		printf("Program did %d iterations.\n", i);
		for(j=0; j<nproc; j++)
		{
			MPI_Cart_coords(cartesian_communicator, j, ndims, coords);
			for(i=0; i<ndims; i++)
			{
				m_min[i] = (int)((double)coords[i]/(double)dims[i]*n[i]);
				m_max[i] = (int)((double)(coords[i]+1)/(double)dims[i]*n[i]);
				if(m_max[i] > n[i]) m_max[i] = n[i]; // 'n_max[0]' shouldn't be greater than 'n[0]'
				m_range[i] = m_max[i]-m_min[i];
			}
			if(m_range[0] != n_range[0] && m_range[1] != n_range[1])
			{
				MPI_Recv(&buf[m_min[0]][m_min[1]], 1, segblock_rightdown, j, j, cartesian_communicator, &status[0]);
			}
			else if(m_range[0] != n_range[0])
			{
				MPI_Recv(&buf[m_min[0]][m_min[1]], 1, segblock_right, j, j, cartesian_communicator, &status[0]);
			}
			else if(m_range[1] != n_range[1])
			{
				MPI_Recv(&buf[m_min[0]][m_min[1]], 1, segblock_down, j, j, cartesian_communicator, &status[0]);
			}
			else
			{
				MPI_Recv(&buf[m_min[0]][m_min[1]], 1, segblock_main, j, j, cartesian_communicator, &status[0]);
			}
		}
		/* have the master thread write the data to a .pgm file
	   * the name of the pgm is give by:
	   * NAME_reconstruct_pN_{i,l}L.pgm
	   * where
	   *  NAME is the name of the input file
	   * 	N is the number of threads
	   * 	{i,l} is the type of limit
	   * 	and L is the limiting number */
		if(*limits >= 1)
			sprintf(outfilename, "%s/%s_reconstruct_p%d_i%0.f.pgm", outpath, infile, nproc, *limits);
		if(*limits < 1)
			sprintf(outfilename, "%s/%s_reconstruct_p%d_l%f.pgm", outpath, infile, nproc, *limits);
		pgmwrite(outfilename, buf, n[0], n[1]);
	}
}
