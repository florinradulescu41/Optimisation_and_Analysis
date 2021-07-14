/*
 * Tema 2 ASC
 * 2021 Spring
 */
#include "utils.h"
#include <stdlib.h>
#include <string.h>
#include "cblas.h"

/* BLAS implementation */
// op(X) is one of op(X) = X   or   op(X) = X**T (X transposed)

// C = A × B × Bt + At × A = ABBt + AtA
double* my_solver(int N, double *A, double *B) {
	printf("BLAS SOLVER\n");

	double *AtA = calloc(N * N, sizeof(double));
	double *C = calloc(N * N, sizeof(double));

	// C = B * Bt
	// A will be multiplied next step for ABBt and AtA will be added later
	cblas_dgemm(				// Performs C = alpha * op(A) * op(B) + beta * C
			CblasRowMajor,	// Layout of the matrix
			CblasNoTrans,		// Is the 1st matrix used in its trasnposed form?
			CblasTrans,			// Is the 2nd matrix used in its trasnposed form?
			N,							// Row_Size of the 1st matrix
			N,							// Column_Size of the 1st matrix
			N,							// Column_Size of the 2nd matrix
			1.0,						// Multiplication quoeficient (alpha)
			B,							// Pointer to the 1st matrix
			N,							// First dimension of the 1st matrix
			B,							// Pointer to the 2nd matrix
			N,							// First dimension of the 2nd matrix
			0.0,						// 3rd matrix quoeficient (beta)
			C,							// Pointer to the 3rd matrix
			N								// First dimension of the 3rd matrix
	);

	// C = A * C (where C has the previous calculated value of BBt)
	// AtA will be added after it is calculated next step
	cblas_dtrmm(				// B = alpha * op(A) * B or B = alpha * B * op(A)
			CblasRowMajor,	// Layout of the matrix
			CblasLeft,			// Side where 1st matrix is multiplied to 2nd matrix
			CblasUpper,			// Triangularity type of the 1st matrix
			CblasNoTrans,		// Is the 1st matrix used in its transposed form?
			CblasNonUnit,		// Is the 1st matrix supposed to have unitary diagonal?
			N,							// Row_Size of the 2nd matrix
			N,							// Column_Size of the 2nd matrix
			1.0,						// Multiplication quoeficient (alpha)
			A,							// Pointer to the 1st matrix
			N,							// First dimension of the 1st matrix
			C,							// Pointer to the 2nd matrix
			N								// First dimension of the 2nd matrix
	);

	// AtA = A (At will be multiplied next step)
	memcpy(AtA, A, N * N * sizeof(*AtA));

	// AtA = At * AtA (where AtA has the previous copied value of initial A)
	// Check the parameters descriptions of the previous cblas_dtrmm (above)
	cblas_dtrmm(
			CblasRowMajor,
			CblasLeft,
			CblasUpper,
			CblasTrans,
			CblasNonUnit,
			N,
			N,
			1.0,
			A,
			N,
			AtA,
			N
	);

	int i, j;
	// C = C + AtA (where C has the previous calculated value of ABBt)
	for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {
					C[i * N + j] += AtA[i * N + j];
			}
	}

	free(AtA);

	return C;
}
