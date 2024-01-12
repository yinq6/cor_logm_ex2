#include <stdlib.h>
#include <stdio.h>
#include"mkl.h"
#include<math.h>
#include "mkl_dss.h"
#include "mkl_types.h"
#include "mkl_spblas.h"
//This function solve the linear system problem by using Intel-MKL library. The idea is Ax=b ->  x=inv(A)b
void full_lu(double* C, double** A, double* Near, int col, int row)
{
    //A rownumber=row; A column number=col;Near=B;
    //C=INV(A)*B;

    int rowA, columnA;
    double alpha, beta;
    rowA = row; columnA = col; 
    double* MatrixA;


    alpha = 1.0; beta = 0.00;
  
    MatrixA = (double*)mkl_malloc(rowA * columnA * sizeof(double), 64);
    /*
    if (MatrixA == NULL ) {
        printf("\n ERROR: Can't allocate memory for matrices. Aborting... \n\n");
        mkl_free(MatrixA);
        //return 1;
    }*/
    
    for (int i = 0; i < rowA; i++)
        for (int j = 0; j < columnA; j++)
        {
            //MatrixA[i * columnA + j] = *((double*)A + i * columnA + j);
            
            if (MatrixA != NULL) {
                MatrixA[i * columnA + j] = A[i][j];
            }
        }
    /*
    printf("matrixA:\n");
    for (int i = 0; i < rowA; i++) {
        for (int j = 0; j < columnA; j++)
        {
            printf("%lf ", MatrixA[i * columnA + j]);
        }
        printf("\n");
    }
    */
    
    int InfoHelp, * VectorHelp;
    VectorHelp = (int*)mkl_malloc(columnA * sizeof(int), 64);
    for (int i = 0; i < columnA; i++)
        VectorHelp[i] = 0;
    //invA
    InfoHelp = LAPACKE_dgetrf(CblasRowMajor, columnA, columnA, MatrixA, columnA, VectorHelp);
    InfoHelp = LAPACKE_dgetri(CblasRowMajor, columnA, MatrixA, columnA, VectorHelp);
    /*
    * printf("invmatrixA:\n");
    for (int i = 0; i < rowA; i++) {
        for (int j = 0; j < columnA; j++)
        {
            printf("%lf ", MatrixA[i * columnA + j]);
        }
        printf("\n");
    }
    */
    
        
    double* BB, * INV_A_B;
    BB = (double*)mkl_malloc(columnA * 1 * sizeof(double), 64);
    INV_A_B = (double*)mkl_malloc(columnA * 1 * sizeof(double), 64);
    for (int i = 0; i < columnA; i++)
    {
        BB[i] = Near[i];
    }
    for (int i = 0; i < columnA; i++)
    {
        INV_A_B[i] = 0.00;
    }

    //INV_A_B=invA*B;
    cblas_dgemv(CblasRowMajor, CblasNoTrans, rowA, columnA, alpha, MatrixA, columnA, BB, 1, beta, INV_A_B, 1);


    C[0] = 0.00;
    for (int i = 0; i < rowA; i++)
    {
        C[i] = INV_A_B[i];
    }
    mkl_free(MatrixA); mkl_free(BB);  mkl_free(INV_A_B); mkl_free(VectorHelp);

}





int sp_ds(double* solution, MKL_INT* ia, MKL_INT* ja, double* values, double* rhs, MKL_INT N)
{
    MKL_INT nrhs = 1;
    // Descriptor of main sparse matrix properties
    struct matrix_descr descrA;
    // Structure with sparse matrix stored in CSR format
    sparse_matrix_t       csrA;
    sparse_operation_t    transA = SPARSE_OPERATION_NON_TRANSPOSE;
    MKL_INT i, j, ivar;
    MKL_INT nnonzero = 5 * N;
    /* Allocate storage for the solver handle and the right-hand side. */
    _MKL_DSS_HANDLE_t handle;
    //_INTEGER_t error;
    MKL_INT error;
    _CHARACTER_t statIn[] = "determinant", * uplo;
    _DOUBLE_PRECISION_t statOut[5], eps = 1e-8;
    MKL_INT opt = MKL_DSS_DEFAULTS+ MKL_DSS_ZERO_BASED_INDEXING, opt1;
    MKL_INT sym = MKL_DSS_NON_SYMMETRIC;
    MKL_INT type = MKL_DSS_INDEFINITE;
    /* --------------------- */
    /* Initialize the solver */
    /* --------------------- */
    double* residual = new double[N];
    double dvar;
    ivar = N;

    error = dss_create(handle, opt);
    if (error != MKL_DSS_SUCCESS)
        goto printError;
    /* ------------------------------------------- */
    /* Define the non-zero structure of the matrix */
    /* ------------------------------------------- */
    error = dss_define_structure(handle,sym, ia, ivar,  ivar, ja, nnonzero);
    if (error != MKL_DSS_SUCCESS)
        goto printError;
    /* ------------------ */
    /* Reorder the matrix */
    /* ------------------ */
    opt = MKL_DSS_DEFAULTS;
    error = dss_reorder(handle, opt, 0);
    if (error != MKL_DSS_SUCCESS)
        goto printError;
    /* ------------------ */
    /* Factor the matrix  */
    /* ------------------ */
    error = dss_factor_real(handle, type, values);
    if (error != MKL_DSS_SUCCESS)
        goto printError;
    /* ------------------------ */
    /* Get the solution vector for Ax=b and ATx=b and check correctness */
    /* ------------------------ */
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
    descrA.mode = SPARSE_FILL_MODE_UPPER;
    descrA.diag = SPARSE_DIAG_NON_UNIT;
    //mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ONE, nRows, nCols, rowIndex, rowIndex + 1, columns, values);
    mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ZERO, N, N, ia, ia + 1, ja, values);
    
   opt = MKL_DSS_DEFAULTS;
  
        // Nullify solution on entry (for sure)
        for (j = 0; j < N; j++)
            solution[j] = rhs[j];

        // Apply trans or non-trans option, solve system
        
        error = dss_solve_real(handle, opt, rhs, nrhs, solution);
        if (error != MKL_DSS_SUCCESS)
            goto printError;
        

     mkl_sparse_d_mv(transA, 1.0, csrA, descrA, solution, 0.0, residual);
     dvar = -1.0E0;
     i = 1;
     //cblas_daxpy(N, dvar, b, 1, residual, 1);  
     daxpy(&N, &dvar, rhs, &i, residual, &i);// residual = residual-b;

     dvar = dnrm2(&N, residual, &i);

            if (dvar>eps)
            {
                
                printf("Incorrect solution\n");
                error = 1000 + i;
                goto printError;
            }
        /*
          printf("Print solution array: ");
        for (j = 0; j < N; j++)
            printf(" %g", solution[j]);

        printf("\n");
    */
        
      
    /* -------------------------- */
    /* Deallocate solver storage  */
    /* -------------------------- */
    mkl_sparse_destroy(csrA);
    error = dss_delete(handle, opt);
    if (error != MKL_DSS_SUCCESS)
        goto printError;
    /* ---------------------- */
    /* Print solution vector  */
    /* ---------------------- */
    delete[] residual;
    //printf("\nExample successfully PASSED!\n");
    return 0;
printError:
    delete[] residual;
    printf("Solver returned error code %d\n", error);
    return 1;
}