#include "parallel_reconstruct.h"

void p_reconstruct(char *inpath, char *infile, char *outpath, int limit_type, float *limits)
{
	setbuf(stdout, NULL);
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

	/* integers required for the various loops */
	int i=0, j=0, keep_going=1;

	/* floats used when the limit is to a minimum change in the pixels */
	float max_change, max, total;
	
	/* things required for setting the cartesian communicator */
	int cart_comm; // this should be MPI_Comm, but it seems to complain about using that
	int ndims = 2;
	int dims[ndims], period[ndims], coords[ndims];
	int left, right, up, down; //the neighbouring process's rank

	for(i=0;i<ndims; i++)
	{
		dims[i] = 0; 	 // automatically choose processers per dim
		period[i] = 0; // do not want periodicity
	}
	MPI_Dims_create(nproc, ndims, dims);
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, period, 0, &cart_comm);
	MPI_Cart_shift(cart_comm, 0, 1,  &left, &right);
	MPI_Cart_shift(cart_comm, 1, 1,  &up, &down);
	MPI_Cart_coords(cart_comm, rank, ndims, coords);
	
	/* this section sets the sizes of the file, according to the image */
	sprintf(infilename, "%s/%s", inpath, infile);
	pgmsize(infilename, &n[0], &n[1]);
	for(i=0; i<2; i++)
	{
		n_min[i] = n[i]*coords[i]/dims[i];
		n_max[i] = n[i]*(coords[i]+1)/dims[i];
		if(n_max[i] > n[i]) n_max[i] = n[i]; // 'n_max[0]' shouldn't be greater than 'n[0]'
		n_range[i] = n_max[i]-n_min[i];
	}

	/* the arrays that actually contain the data to do with the reconstruction */
	float buf[n[0]][n[1]];
	float edge[n_range[0]][n_range[1]];
	float out[n_range[0]][n_range[1]];
	float segment[n_range[0]+2][n_range[1]+2]; //+2 because it contains 1px halos


	/* things required for the datatype that splits the image up into segments.
	 * segblock is in 4, because there are two values that any given segment
	 * might have for width or height, due to the floor of integer division
	 * [0] = "inner" segments, which have the same size as the 0th segment
	 * [1] = "down" segments, which have a different height
	 * [2] = "right" segments, which have a different width
	 * [3] = bottom right segment, which has both a different height and width
	 *
	 * Of course, there are "inner" segments with different heights,
	 * it's just easier to visualise what is happening this way
	 * */
	MPI_Datatype segblock[4];
	int n_range_right = n[0]-n[0]*(dims[0]-1)/dims[0];
	int n_range_down = n[1]-n[1]*(dims[1]-1)/dims[1];
	
	MPI_Type_vector(n_range[0], n_range[1], n[1], MPI_FLOAT, &segblock[0]);
	MPI_Type_commit(&segblock[0]);

	MPI_Type_vector(n_range[0], n_range_down, n[1], MPI_FLOAT, &segblock[1]);
	MPI_Type_commit(&segblock[1]);

	MPI_Type_vector(n_range_right, n_range[1], n[1], MPI_FLOAT, &segblock[2]);
	MPI_Type_commit(&segblock[2]);

	MPI_Type_vector(n_range_right, n_range_down, n[1], MPI_FLOAT, &segblock[3]);
	MPI_Type_commit(&segblock[3]);
	
	if(rank == 0)	pgmread(infilename, buf, n[0], n[1]);
	scatter_2D(buf, edge, rank, nproc, cart_comm, ndims, n, n_range, dims, segblock);
	initialise_segment(segment, n_range[0]+2, n_range[1]+2);

	i=0;
	while(keep_going == 1)
	{
		i++;
		reconstruct_image_segment(segment, edge, n_min[0], n_max[0], n_min[1], n_max[1], &max_change, &total);
		swap_halos(n, dims, n_range, segment, up, down, left, right, cart_comm);

		/* choose when to print average pixel value */
		if(i % 100 == 0)
		{
			float aver;
			MPI_Reduce(&total, &aver, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
			if(rank == 0) printf("Average Pixel Value at iteration %d = %f\n", i, aver/(n[0]*n[1]));
			
		}
		
		if(i >= *limits && limit_type == 0) keep_going = 0;
		else if(limit_type == 1)
		{
			MPI_Allreduce(&max_change, &max, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
			if(max <= *limits ) keep_going = 0;
		}
	}

	removehalo(segment, out, n_range[0], n_range[1]);
	gather_2D(buf, out, rank, nproc, cart_comm, ndims, n, n_range, dims);
	
	/* have the master thread write the data to a .pgm file
	  * the name of the pgm is give by:
	  * INPUTFILENAME_reconstruct_pNPROC_{i,l}LIMIT.pgm */
	MPI_Reduce(&max_change, &max, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
	if(rank == 0)
	{
		printf("\nProgram did %d iterations.\n\nGreatest change on this iteration was: %f\n", i, max);
		if(limit_type == 0)
			sprintf(outfilename, "%s/%s_reconstruct_p%d_i%0.f.pgm", outpath, infile, nproc, *limits);
		else if(limit_type == 1)
			sprintf(outfilename, "%s/%s_reconstruct_p%d_l%f.pgm", outpath, infile, nproc, *limits);
		pgmwrite(outfilename, buf, n[0], n[1]);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}
	
void swap_halos(int n[], int dims[], int n_range[], void *seg, int up, int down, int left, int right, int cart_comm)
{
		/* the 'segment' array containing the 
		 * image w/ halo for this iteration */
		float *segment = (float *) seg;
		int seg_width = (n_range[1]+2);

		/* these are here to simulate the values of the ranges of the 
		 * outer segments */
		int n_range_right = n[0]-((dims[0]-1)*n[0])/dims[0];
		int n_range_down = n[1]-((dims[1]-1)*n[1])/dims[1];

		
		/* needed to send noncontiguous data */
		MPI_Datatype noncontig, noncontig_right, vect_type;
		MPI_Type_vector(n_range[0], 1, n_range[1]+2, MPI_FLOAT, &noncontig);
		MPI_Type_commit(&noncontig);
		MPI_Type_vector(n_range_right, 1, n_range[1]+2, MPI_FLOAT, &noncontig_right);
		MPI_Type_commit(&noncontig_right);

		MPI_Request request[4];
		MPI_Status	status[4];

		vect_type = noncontig;
		int point_location_up = seg_width+1; // brings us to the first non-halo point
		int point_location_down = seg_width+n_range[1]; //brings us to the start of the bottom row
		int point_location_left = seg_width+1; // brings us to the first non-halo point
		int point_location_right = n_range[0]*seg_width+1; //brings us the start of the right row
		int range = n_range[1];
		
		/* this part is neccessary because the outer edges of the segments
		 * can have different sizes to the inner parts */
		if(n[1]/dims[1] != n_range[1])
		{
			range = n_range_down;
			point_location_down = seg_width+n_range_down;	
		}
		if(n[0]/dims[0] != n_range[0])
		{
			vect_type = noncontig_right;
			point_location_right = n_range_right*seg_width+1;
		}

		/* send and recieve the 4 lines that make up the halo */
		MPI_Issend(&segment[point_location_up]		,1		,vect_type,up		,0,cart_comm, &request[0]);
		MPI_Issend(&segment[point_location_down]	,1		,vect_type,down	,1,cart_comm, &request[1]);
		MPI_Issend(&segment[point_location_left]	,range,MPI_FLOAT,left	,2,cart_comm, &request[2]);
		MPI_Issend(&segment[point_location_right]	,range,MPI_FLOAT,right,3,cart_comm, &request[3]);
			
		MPI_Recv(&segment[point_location_down+1]					,1		,vect_type,down	,0,cart_comm, &status[0]);
		MPI_Recv(&segment[seg_width]											,1		,vect_type,up		,1,cart_comm, &status[1]);
		MPI_Recv(&segment[point_location_right+seg_width]	,range,MPI_FLOAT,right,2,cart_comm, &status[2]);
		MPI_Recv(&segment[1]															,range,MPI_FLOAT,left	,3,cart_comm, &status[3]);
}

void scatter_2D(void *buf_in, void *edge_in, int rank, int nproc, int cart_comm, int ndims, int n[], int n_range[], int dims[], MPI_Datatype segblock[])
{
	float *buf = (float *) buf_in;
	float *edge = (float *) edge_in;
	int i,j;
	int m_min[ndims], m_max[ndims], m_range[ndims];
	int coords[ndims];
	MPI_Request request;
	MPI_Status status;
	if(rank == 0)
	{
		for(j=0; j<nproc; j++)
		{
			MPI_Cart_coords(cart_comm, j, ndims, coords);
			for(i=0; i<ndims; i++)
			{
				m_min[i] = n[i]*coords[i]/dims[i];
				m_max[i] = n[i]*(coords[i]+1)/dims[i];
				if(m_max[i] > n[i]) m_max[i] = n[i]; // 'n_max[0]' shouldn't be greater than 'n[0]'
				m_range[i] = m_max[i]-m_min[i];
			}
			if(m_range[0] != n_range[0] && m_range[1] != n_range[1])
				MPI_Issend(&buf[m_min[0]*n[1]+m_min[1]], 1, segblock[3], j, 0, cart_comm, &request);
			else if(m_range[0] != n_range[0])
				MPI_Issend(&buf[m_min[0]*n[1]+m_min[1]], 1, segblock[2], j, 0, cart_comm, &request);
			else if(m_range[1] != n_range[1])
				MPI_Issend(&buf[m_min[0]*n[1]+m_min[1]], 1, segblock[1], j, 0, cart_comm, &request);
			else 
				MPI_Issend(&buf[m_min[0]*n[1]+m_min[1]], 1, segblock[0], j, 0, cart_comm, &request);
		}
	}
	MPI_Recv(edge, n_range[0]*n_range[1], MPI_FLOAT, 0, 0, cart_comm, &status);
}

void gather_2D(void *buf_in, void *out_in, int rank, int nproc, int cart_comm, int ndims, int n[], int n_range[], int dims[])
{
	float *buf = (float *) buf_in;
	float *out = (float *) out_in;
	int i,j;
	int m_min[ndims], m_max[ndims], m_range[ndims];
	int coords[ndims];
	MPI_Request request;
	MPI_Status status;

	MPI_Issend(out, n_range[0]*n_range[1], MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &request);
	if(rank ==0)
	{
		for(j=0; j<nproc; j++)
		{
			MPI_Cart_coords(cart_comm, j, ndims, coords);
			/* imitate the values of the other processors */
			for(i=0; i<ndims; i++)
			{
				m_min[i] = n[i]*coords[i]/dims[i];
				m_max[i] = n[i]*(coords[i]+1)/dims[i];
				if(m_max[i] > n[i]) m_max[i] = n[i]; // 'n_max[0]' shouldn't be greater than 'n[0]'
				m_range[i] = m_max[i]-m_min[i];
			}
			float hold[m_range[0]][m_range[1]];
			MPI_Recv(hold, m_range[0]*m_range[1], MPI_FLOAT, j, 0, MPI_COMM_WORLD, &status);
			int x,y;
			for(x=0;x<m_range[0];x++)
			{
				for(y=0;y<m_range[1];y++)
				{
					buf[(m_min[0]+x)*n[1]+m_min[1]+y] = hold[x][y];
				}
			}
		}
	}
}
