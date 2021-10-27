# Matrix_Optimisation
Optimised Matrix Calculus and Analysis

Repository structure:

Makefile, solver_blas.c [1], solver_neopt [2], solver_opt [3], grafic [4]
for comparative evaluation and 3 separate individual graphs [5] [6] [7].
The mentioned .c files contain different implementations of the operation in the theme,
C = A × B × Bt + At × A = ABBt + AtA, where the result matrix C is considered
and the two matrices given as a parameter of the my_solver function A and B, the three
being square matrices of size N * N with double type elements, and A
being a superior triangular matrix. Note Xt transposed to the matrix X.

[1] solver_blas.c

Implementation description:
Functions in the cblas.h library are used to perform various operations
on the matrices (multiplications and additions) in order to achieve the result.
The result matrix C will store the partial results of the first sequence of
multiplications, and the second sequence will be added at the end. So,
using the cblas_dgemm function (parameters detailed in the comments), C arrives
to have the value B * Bt (transposed to B is not calculated explicitly, but is
uses a function parameter to announce the calculation mode).
Then, on C (which now has the value B * Bt) apply a multiplication with
matrix A using the cblas_dtrmm function (parameters explained in the comments),
operation a trace whose value of C will be A * B * Bt, ie the first term
of the assembly from the initial function. Separately, using the same cblas_dtrmm function
At * A is calculated, multiplying the transposed matrix At to the left of the matrix
AtA, which contains copied the initial value of A, resulting in At * A.
The two terms, C and AtA are summed at the end in the final matrix C.

[2] solver_neopt.c

Implementation description:
The non-optimized variant contains the classical calculation of the matrix multiplication in C,
using iterations through for loops. The first loop is used to scroll
the lines of the matrix B, then the columns in both matrices and the lines in the second
matrix. However, given that the second matrix in the first nested loop
it is a transpose, its passage in the last forum is done by reversing the positions
line and column indices (both matrices thus appear to be traversed "identically").
The difference can be seen at the next nesting of for type clauses, where
multiplication is done classically. These first two nests calculate, in turn,
B * Bt in the matrix BBt and A * BBt in the matrix ABBt. The AtA matrix therefore contains
the result of the multiplication between the matrix A and its transpose, in which case, the same
as in the previous A * BBt multiplication, the described case is specially treated
in enunciation in the form of the property of the superior triangular matrix which
has the matrix A. Therefore, in both iterations, the last for loop is
modified, iterating ignoring null elements (in the first case, it starts on
each row in the column with index immediately above that row, and in al
the second case stopping at the minimum of the line index in the matrix A and of the
of the column in the matrix At, which in turn is triangular, although inferior).
The last nesting is used to add the two terms, ABBt and AtA in
as a result of the returned matrix, C.

[3] solver_opt.c

Implementation description:
The variant follows the same approach as solver_neopt.c, but with an efficiency
brought to the calculation process that makes the program run faster.
Starting from the observation that denotes that in the last for loop of
nests used for matrix multiplication, the result can be treated
found at a certain index as a constant. As a result, none will be added
value to it in the iterations, but only after establishing the sum of all
the coefficients involved. This temporary amount is kept in a variable
register type to speed up memory access and writing
within it. How the change affects the execution time of
this approach can be verified in the graph in point [4], in
which shows a speedup compared to the initial, non-optimized version. A mention
the important thing is that the optimizations used were not enough to reduce
the execution time of this implementation below the threshold described as satisfactory
in the theme statement (12 seconds on the ibm-nehalem.q queue for the test with
matrix size of N * N = 1200 * 1200 elements).

[4] graphic

The graph is used as a basis for time-based benchmarking
running between the three implementation variants mentioned above. So,
for each of the variants a number was made on the same architecture
of 6 tests, each test having a different value for the overall size
of the calculation matrices N * N, gradually varying the value of N from 400 to 1400
with a step of 200. The following results were obtained for the base of the graph:

BLAS, NEOPT, OPT = run using solver_blas.c, solver_neopt.c, solver_opt.c.
Time is the time expressed in seconds elapsed to obtain the correct result.

N = 400:
BLAS SOLVER Time = 0.050631
NEOPT SOLVER Time = 1.112644
OPT SOLVER Time = 0.746751

N = 600:
BLAS SOLVER Time = 0.107851
NEOPT SOLVER Time = 3.570956
OPT SOLVER Time = 2.337827

N = 800:
BLAS SOLVER Time = 0.202301
NEOPT SOLVER Time = 8.536873
OPT SOLVER Time = 5.608161

N = 1000:
BLAS SOLVER Time = 0.383041
NEOPT SOLVER Time = 16.422619
OPT SOLVER Time = 10.776759

N = 1200:
BLAS SOLVER Time = 0.652319
NEOPT SOLVER Time = 28.680351
OPT SOLVER Time = 18.785624

N = 1400:
BLAS SOLVER Time = 1.029703
NEOPT SOLVER Time = 46.510609
OPT SOLVER Time = 30.587818

Conclusions: It can be seen that the runs using the functions in cblas.h hold
considerably less than those using the classic, optimized or variants
no, and this difference is all the more visible on the chart. Growth in
time of running the BLAS type variants for a progression of dimensions of
The matrix is ​​almost a linear function, being only almost double with each
step compared to the previous version (the value being very small, the doubling factor
it is almost imperceptible). On the other hand, on average, the growth of others
methods in the graph (NEOPT AND OPT_M) can be approximated to a double linearity,
but on a considerably high spectrum of values, which makes them seem more
much the equivalent of exponential functions. In reality, you can see how
each increase is achieved with similar uniformities in the other three
separate graphics included for each implementation variant [5] [6] [7].
