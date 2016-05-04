//
//  main.cpp
//  Rosenbrock-MPI
//
//  Created by 雪竜 on 16/4/27.
//
//

#include <iostream>
#include <mpi.h>
#include "MatrixUtility.h"
#include "LDLUtility.h"

int dimension = 10;
int processor_id, processor_count; //Worker ID, aka rank

#define ___PROBLEM_SIZE___ dimension
#define ___MASTER_PROCESSOR___ (processor_count - 1)
#define ___IS_MASTER_PROCESSOR___ processor_id == ___MASTER_PROCESSOR___
#define ___IS_PARALLEL___ processor_count > 1
#define ___TAG___ 574
#define ___TOLERANCE___ 0.0001

double F_xn(const double* v, const int n) {
	
	int res = 0;
	
	if (n < ___PROBLEM_SIZE___ - 1)
		res = 400 * ( v[n] * v[n] - v[n+1]) * v[n] + 2 * (v[n] - 1);
	
	if(n > 0)
		res += -200 * (v[n - 1]*v[n - 1] - v[n]);
	
	return res;
}

double F_xn_xn(const double* v, const int n) {
	int res = 0;
	
	if (n < ___PROBLEM_SIZE___ - 1)
		res = 1200 * v[n] * v[n] - 400 * v[n+1] + 2;
	
	if(n > 0)
		res += 200;
	
	return res;
}

double F_xn_xn_minus_1(const double* v, const int n) {
	return -400*v[n-1];
}

double F_xn_xn_plus_1(const double* v, const int n) {
	return -400*v[n];
}

bool is_my_row(int row_number, int processor_count, int processor_id){
	return (row_number % processor_count) == processor_id;
}

void array_copy(double* dst, double* src, int size){
	for (int i  = 0; i < size; i++) {
		dst[i] = src[i];
	}
}

int main(int argc, char** argv) {
	
	if(argc > 1)
		___PROBLEM_SIZE___ = atoi(argv[1]);
	
	// MPI code starts here...
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &processor_count);
	MPI_Comm_rank(MPI_COMM_WORLD, &processor_id);
	MPI_Status mpi_stat;
	
	double** A = MatrixUtility::InitializeEmptyMatrix(___PROBLEM_SIZE___, ___PROBLEM_SIZE___);
	double** A_res = MatrixUtility::InitializeEmptyMatrix(___PROBLEM_SIZE___, ___PROBLEM_SIZE___);
	double* B = MatrixUtility::InitializeEmptyVector(___PROBLEM_SIZE___);
	double* A_row_k = MatrixUtility::InitializeEmptyVector(___PROBLEM_SIZE___);
	double* X = MatrixUtility::InitializeEmptyVector(___PROBLEM_SIZE___);
	
	for (int i = 0; i < ___PROBLEM_SIZE___; i++) {
		X[i] = i+1;
	}
	
	/*
	
	A[0][0] = 2;
	A[0][1] = -1;
	A[0][2] = 0;
	A[1][0] = -1;
	A[1][1] = 2;
	A[1][2] = -1;
	A[2][0] = 0;
	A[2][1] = -1;
	A[2][2] = 1;
	
	B[0] = 1;
	B[1] = 0;
	B[2] = 0;
	 */
	
	double residual = 0;
	
	do {
		
		int row_per_processor = ___PROBLEM_SIZE___ / processor_count;
		
		int row_range_min = processor_id * row_per_processor;
		int row_range_max = row_range_min + row_per_processor - 1;
		
		if (___IS_MASTER_PROCESSOR___) {
			row_range_max = ___PROBLEM_SIZE___ - 1;
			
			for (int i = 0; i < ___PROBLEM_SIZE___; i++) {
				B[i] = -F_xn(X, i);
				
				A[i][i] = F_xn_xn(X, i);
				
				
				if(i < ___PROBLEM_SIZE___ - 1)
				{
					A[i][i + 1] = F_xn_xn_plus_1(X, i);
				}
				
				if(i > 0)
				{
					A[i][i - 1] = F_xn_xn_minus_1(X, i);
				}
			}
			
			if (___PROBLEM_SIZE___ <= 30) {
				printf("\nB:\n");
				for (int i = 0; i < ___PROBLEM_SIZE___; i++) {
					std::cout << std::setw(12) << B[i];
				}
				
				printf("\nA:\n");
				MatrixUtility::PrintMatrix(A, ___PROBLEM_SIZE___, ___PROBLEM_SIZE___);
			}
		}
		
		/*------------GET LDLT FACTOR------------*/
		
		for (int k = 0; k < ___PROBLEM_SIZE___; k++) {
			MPI_Bcast(A[k], ___PROBLEM_SIZE___, MPI_DOUBLE, ___MASTER_PROCESSOR___, MPI_COMM_WORLD);
		}
		
		for (int k = 0; k < ___PROBLEM_SIZE___; k++) {
			if (k >= row_range_min && k <= row_range_max) {
				for (int j = k + 1; j < ___PROBLEM_SIZE___; j++) {
					A[k][j] /= A[k][k];
				}
				
				if (___IS_PARALLEL___) {
					for (int i = processor_id; i < processor_count; i++) {
						if (i != processor_id) {
							MPI_Send(A[k], ___PROBLEM_SIZE___, MPI_DOUBLE, i, ___TAG___, MPI_COMM_WORLD);
						}
					}
				}
			}
			
			if (___IS_PARALLEL___){
				if (k <= row_range_max ) {
					if (k < row_range_min) {
						MPI_Recv(A_row_k, ___PROBLEM_SIZE___, MPI_DOUBLE, MPI_ANY_SOURCE, ___TAG___, MPI_COMM_WORLD, &mpi_stat);
					}
					else
					{
						array_copy(A_row_k, A[k], ___PROBLEM_SIZE___);
					}
					
					if (___IS_MASTER_PROCESSOR___) {
						array_copy(A_res[k], A_row_k, ___PROBLEM_SIZE___);
					}
				}
			}
			
			MPI_Barrier(MPI_COMM_WORLD);
			
			for (int i = ( ( (k+1) > row_range_min) ? (k+1) : row_range_min); i<=row_range_max; i++){
				for(int j = k + 1; j < ___PROBLEM_SIZE___; j++){
					if(___IS_PARALLEL___)
						A[i][j] = A[i][j] - A[i][k] * A_row_k[j];
					else
						A[i][j] = A[i][j] - A[i][k] * A[k][j];
				}
			}
		}
		
		if (___IS_PARALLEL___ && ___PROBLEM_SIZE___ <= 30)
		{
			if (___IS_MASTER_PROCESSOR___)
			{
				printf("\nA_res:\n");
				MatrixUtility::PrintMatrix(A_res, ___PROBLEM_SIZE___, ___PROBLEM_SIZE___);
			}
		}
		else
		{
			if (___PROBLEM_SIZE___ <= 30) {
				printf("\nA_res:\n");
				MatrixUtility::PrintMatrix(A, ___PROBLEM_SIZE___, ___PROBLEM_SIZE___);
			}
			
		}
		
		/*------------GET LDLT FACTOR FINISHED------------*/
		
		/*------------BACKWARD SUBSTITUTION------------*/
		
		if (___IS_MASTER_PROCESSOR___)
		{
			double* D;
			double** L_trans;
			if (___IS_PARALLEL___)
			{
				D = MatrixUtility::ExtractDiagnoal(A_res, ___PROBLEM_SIZE___);
				L_trans = MatrixUtility::GetMatrixTranspose(A_res, ___PROBLEM_SIZE___, ___PROBLEM_SIZE___);
			}
			else
			{
				D = MatrixUtility::ExtractDiagnoal(A, ___PROBLEM_SIZE___);
				L_trans = MatrixUtility::GetMatrixTranspose(A, ___PROBLEM_SIZE___, ___PROBLEM_SIZE___);
			}
			
			double* Y = new double[___PROBLEM_SIZE___];
			double* Z = new double[___PROBLEM_SIZE___];
			double* X_delta = MatrixUtility::InitializeEmptyVector(___PROBLEM_SIZE___);
			
			LDLUtility::SolveXWithLD(L_trans, ___PROBLEM_SIZE___, D, B, Z, Y, X_delta);
			//MatrixUtility::PrintVector(X_delta, ___PROBLEM_SIZE___);
			
			for (int i = 0; i < ___PROBLEM_SIZE___; i++) {
				X[i] += X_delta[i];
			}
			
			if (___PROBLEM_SIZE___ <= 30) {
				printf("\nX_delta:\n");
				for (int i = 0; i < ___PROBLEM_SIZE___; i++) {
					std::cout << std::setw(12) << X_delta[i];
				}
				printf("\nX:\n");
				for (int i = 0; i < ___PROBLEM_SIZE___; i++) {
					std::cout << std::setw(12) << X[i];
				}
				
				std::cout << std::endl;
				std::cout << std::endl;
			}
			residual = MatrixUtility::ComputeResidual(X_delta, ___PROBLEM_SIZE___);
			
			printf("Residual: %f\n", residual);
			
			delete [] Y;
			delete [] Z;
			delete [] X_delta;
		}
		
		MPI_Bcast(&residual, 1, MPI_DOUBLE, ___MASTER_PROCESSOR___, MPI_COMM_WORLD);
		
	} while (residual > ___TOLERANCE___);
	
	/*
	
	if (___IS_MASTER_PROCESSOR___) {
		
		printf("\nB:\n");
		for (int i = 0; i < ___PROBLEM_SIZE___; i++) {
			std::cout << std::setw(12) << B[i];
		}
		
		printf("\nA:\n");
		MatrixUtility::PrintMatrix(A, ___PROBLEM_SIZE___, ___PROBLEM_SIZE___);
		
		printf("\nA_res:\n");
		MatrixUtility::PrintMatrix(A, ___PROBLEM_SIZE___, ___PROBLEM_SIZE___);
		
		printf("\nX:\n");
		for (int i = 0; i < ___PROBLEM_SIZE___; i++) {
			std::cout << std::setw(12) << X[i];
		}
		
		std::cout << std::endl;
	}
	 */
	
	//printf("done\n");
	
	/*------------BACKWARD SUBSTITUTION FINISHED------------*/
	
	MPI_Finalize();
	
	return 0;
}
