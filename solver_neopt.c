/*
 * Tema 2 ASC
 * 2021 Spring
 */
#include "utils.h"

// Function to get the minimum value of two variables
double minim(double i, double j) {
		if (i <= j) return i;
		else return j;
}

/* Unoptimized implementation */

// C = A × B × Bt + At × A
double* my_solver(int N, double *A, double* B) {
		printf("NEOPT SOLVER\n");

		double *BBt = calloc(N * N, sizeof(double));
		double *ABBt = calloc(N * N, sizeof(double));
		double *AtA = calloc(N * N, sizeof(double));
		double *C = calloc(N * N, sizeof(double));

		int i, j, k;

		// BBt = B × Bt
		// Bt will be considered just B iterated on columns instead of on lines
		for (i = 0; i < N; i++) {
				for (j = 0; j < N; j++) {
						for (k = 0; k < N; k++) {
								BBt[i * N + j] += B[i * N + k] * B[j * N + k];
						}
				}
		}

		// ABBt = A × BBt
		// Considering A a superior triangular matrix, the first i - 1 elements
		// on each row of the result will be 0, so iterate k = {1...N} interval
		for (i = 0; i < N; i++) {
				for (j = 0; j < N; j++) {
						for (k = i; k < N; k++) {
								ABBt[i * N + j] += A[i * N + k] * BBt[k * N + j];
						}
				}
		}

		// AtA = At × A
		// Considering both A and At as triangulars, A superior and At inferior
		// At will be considered just A iterated on columns instead of on lines
		for (i = 0; i < N; i++) {
				for (j = 0; j < N; j++) {
						for (k = 0; k < minim(i, j) + 1; k++) {
								AtA[i * N + j] += A[k * N + i] * A[k * N + j];
						}
				}
		}

		// C = ABBt + AtA
		for (i = 0; i < N; i++) {
				for (j = 0; j < N; j++) {
						C[i * N + j] = ABBt[i * N + j] + AtA[i * N + j];
				}
		}

		free(BBt);
		free(ABBt);
		free(AtA);

		return C;
}
