#include <iostream>
#include <iomanip>
#include <cmath>
#include <mpi.h>

int main(int argc, char *argv[]){ int i,j,k,size,index;
	int index1,index2;
	int mynode, totalnodes;
	double alpha,gamma;
	const int numrows = 5;
	
	MPI_Status status;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	
	size = (int) pow(2,log2(totalnodes+1)+1)-1;
	double ** A = new double*[numrows]; for(i=0;i<numrows;i++){
  A[i] = new double[size+1];
  for(j=0;j<size+1;j++)
	  A[i][j] = 0.0;
	}
	if(mynode==0){
		A[0][0] = -2.0; A[0][1] = 1.0;
		A[1][0] = 1.0; A[1][1] = -2.0; A[1][2] = A[2][1] = 1.0; A[2][2] = -2.0; A[2][3] =
	}
	else if(mynode==(totalnodes-1)){
		index = 2*mynode;
		A[0][index-1] = 1.0; A[0][index] = -2.0; A[0][index+1] = 1.0;
		index = 2*mynode+1;
		A[1][index-1] = 1.0; A[1][index] = -2.0; A[1][index+1] = 1.0;
		A[2][size-2] = 1.0; A[2][size-1] = -2.0;
	} else{
  for(i=0;i<3;i++){
	  index = i + 2*mynode;
	  1.0; 1.0;
  } }
	A[i][index-1]
	A[i][index]
	A[i][index+1]
	= 1.0;
	= -2.0;
	= 1.0;
	for(i=0;i<3;i++) A[i][size] = 2*mynode+i;
	int numactivep = totalnodes;
	int * activep = new int[totalnodes]; for(j=0;j<numactivep;j++)
  activep[j] = j;
	for(j=0;j<size+1;j++){
  A[3][j] = A[0][j];
  A[4][j] = A[2][j];
	}