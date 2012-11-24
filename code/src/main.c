#include <stdlib.h>
#include <stdio.h>
#include <math.h>
int main()
{
	int iterations = 1000;
	//char infilename[256] = "../input_files/edge192x128.pgm";
	char infilename[256] = "../input_files/edge256x192.pgm";
	//char infilename[256] = "../input_files/edge512x384.pgm";
	//char infilename[256] = "../input_files/edge768x768.pgm";
	//char infilename[256] = "../input_files/edge1024x1408.pgm";
	
	
	p_reconstruct(infilename, iterations);
	//reconstruct(infilename, iterations);
	return 0;
}
