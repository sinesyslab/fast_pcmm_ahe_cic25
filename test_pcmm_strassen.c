// This test is to profile plaintext-ciphertext matrix multiplication using Schoolbook and Strassen approaches.

// -------------------------------------------------------------------- Test Begins -----------------------------------------------------------------------

#include <stdio.h> // for printf
#include <stdlib.h> // for malloc and free
#include <string.h>
#include <time.h> // To measure time and do profiling. Required for also random number generation.
#include <math.h> // To use "round" function for rounding the vector lengths after each iteration to nearest integer.
#include "randapi.h" // For CREATE_CSPRNG function call

#include "test_pcmm_strassen.h"

// Include NIST256 curve specific header files
#include "ecp_NIST256.h"
#include "fp_NIST256.h"

#define TIME_PROFILE

// --------------------------------------------------------- Helper Functions Begin --------------------------------------------------------------------------

// This function is used in time profiling
long timespec_diff (struct timespec start, struct timespec end)
{
    struct timespec temp;
    if ((end.tv_nsec-start.tv_nsec) < 0) {
        temp.tv_sec = end.tv_sec-start.tv_sec-1;
        temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    } else {
        temp.tv_sec = end.tv_sec-start.tv_sec;
        temp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }
 
    return ((1e9*temp.tv_sec) + temp.tv_nsec);
}

// Function to allocate memory for a plaintext matrix
int** allocateMatrix(int size) { 
    // Allocate memory for the array of row pointers
    int** matrix = (int**)malloc(size * sizeof(int*)); 
    if (matrix == NULL) {
        perror("Failed to allocate memory for matrix rows");
        exit(EXIT_FAILURE);
    }
    
    // Allocate memory for each row
    for (int i = 0; i < size; i++) {
        matrix[i] = (int*)malloc(size * sizeof(int)); 
        if (matrix[i] == NULL) {
            perror("Failed to allocate memory for matrix row");
            exit(EXIT_FAILURE);
        }
    }
    return matrix;
}

// Function to allocate memory for a ciphertext matrix
ECP_NIST256** allocate_ECP_Matrix(int size) { 
    // Allocate memory for the array of row pointers
    ECP_NIST256** matrix = (ECP_NIST256**)malloc(size * sizeof(ECP_NIST256*)); 
    if (matrix == NULL) {
        perror("Failed to allocate memory for matrix rows");
        exit(EXIT_FAILURE);
    }
    
    // Allocate memory for each row
    for (int i = 0; i < size; i++) {
        matrix[i] = (ECP_NIST256*)malloc(size * sizeof(ECP_NIST256)); 
        if (matrix[i] == NULL) {
            perror("Failed to allocate memory for matrix row");
            exit(EXIT_FAILURE);
        }
    }
    return matrix;
}

// Function to free memory of a matrix
void freeMatrix(int** matrix, int size) {
    for (int i = 0; i < size; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

// Function to free memory of a ECP matrix
void free_ECP_Matrix(ECP_NIST256** matrix, int size) {
    for (int i = 0; i < size; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

// Function to add two matrices
void addMatrices(int** A, int** B, int** C, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}

// Function to add two ECP matrices
void add_ECP_Matrices(ECP_NIST256** A, ECP_NIST256** B, ECP_NIST256** C, int size) {
    ECP_NIST256 add_temp1, add_temp2;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            add_temp1 = A[i][j];
            add_temp2 = B[i][j];
            // Adds ECP instance Q to ECP instance P - extern void ECP_NIST256_add(ECP_NIST256 *P, ECP_NIST256 *Q);
            ECP_NIST256_add(&add_temp1/* ECP_NIST256 *P */, &add_temp2/* ECP_NIST256 *Q */);
            C[i][j] = add_temp1;
        }
    }
}

// Function to subtract two matrices
void subtractMatrices(int** A, int** B, int** C, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
}

// Function to subtract two ECP matrices
void subtract_ECP_Matrices(ECP_NIST256** A, ECP_NIST256** B, ECP_NIST256** C, int size) {
    ECP_NIST256 sub_temp1, sub_temp2;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            sub_temp1 = A[i][j];
            sub_temp2 = B[i][j];
            // Subtracts ECP instance Q from ECP instance P - on exit =P-Q - extern void ECP_ZZZ_sub(ECP_ZZZ *P, ECP_ZZZ *Q);
            ECP_NIST256_sub(&sub_temp1/* ECP_ZZZ *P */, &sub_temp2/* ECP_ZZZ *Q */);            
            C[i][j] = sub_temp1;
        }
    }
}

// Strassen's Algorithm
void strassen(int** A, ECP_NIST256** B, ECP_NIST256** C, int size, int actual_size) {
    if (size == 1) {
        ECP_NIST256 strassen_temp1;
        strassen_temp1 = B[0][0];
        if (A[0][0] > 0) {
            // Multiplies an ECP instance P by a small integer, side-channel resistant - extern void ECP_ZZZ_pinmul(ECP_ZZZ *P, int i, int b); 
            ECP_NIST256_pinmul(&strassen_temp1/* ECP_NIST256 *P */, A[0][0]/* int i */, (int) (BIT_WIDTH_EACH_PLAINTEXT + log2(actual_size))/* int b */); // "(BIT_WIDTH_EACH_PLAINTEXT + 1)" bits considering (2 ** BIT_WIDTH_EACH_PLAINTEXT) although (2 ** BIT_WIDTH_EACH_PLAINTEXT) can not be the random number generated.
        } else {
            // Multiplies an ECP instance P by a small integer, side-channel resistant - extern void ECP_ZZZ_pinmul(ECP_ZZZ *P, int i, int b); 
            ECP_NIST256_pinmul(&strassen_temp1/* ECP_NIST256 *P */, (A[0][0] * (-1))/* int i */, (int) (BIT_WIDTH_EACH_PLAINTEXT + log2(actual_size))/* int b */); // "(BIT_WIDTH_EACH_PLAINTEXT + 1)" bits considering (2 ** BIT_WIDTH_EACH_PLAINTEXT) although (2 ** BIT_WIDTH_EACH_PLAINTEXT) can not be the random number generated.
            // Negation of an ECP point - extern void ECP_ZZZ_neg(ECP_ZZZ *P);
            ECP_NIST256_neg(&strassen_temp1/* ECP_ZZZ *P */);
        }
        C[0][0] = strassen_temp1;
        return;
    }

    int newSize = size / 2;
    int** A11 = allocateMatrix(newSize);
    int** A12 = allocateMatrix(newSize);
    int** A21 = allocateMatrix(newSize);
    int** A22 = allocateMatrix(newSize);

    ECP_NIST256** B11 = allocate_ECP_Matrix(newSize);
    ECP_NIST256** B12 = allocate_ECP_Matrix(newSize);
    ECP_NIST256** B21 = allocate_ECP_Matrix(newSize);
    ECP_NIST256** B22 = allocate_ECP_Matrix(newSize);

    ECP_NIST256** C11 = allocate_ECP_Matrix(newSize);
    ECP_NIST256** C12 = allocate_ECP_Matrix(newSize);
    ECP_NIST256** C21 = allocate_ECP_Matrix(newSize);
    ECP_NIST256** C22 = allocate_ECP_Matrix(newSize);

    ECP_NIST256** M1 = allocate_ECP_Matrix(newSize);
    ECP_NIST256** M2 = allocate_ECP_Matrix(newSize);
    ECP_NIST256** M3 = allocate_ECP_Matrix(newSize);
    ECP_NIST256** M4 = allocate_ECP_Matrix(newSize);
    ECP_NIST256** M5 = allocate_ECP_Matrix(newSize);
    ECP_NIST256** M6 = allocate_ECP_Matrix(newSize);
    ECP_NIST256** M7 = allocate_ECP_Matrix(newSize);

    int** temp1 = allocateMatrix(newSize);
    ECP_NIST256** temp2 = allocate_ECP_Matrix(newSize);
    ECP_NIST256** temp3 = allocate_ECP_Matrix(newSize);

    // Dividing matrices into 4 submatrices
    for (int i = 0; i < newSize; i++) {
        for (int j = 0; j < newSize; j++) {
            A11[i][j] = A[i][j];
            A12[i][j] = A[i][j + newSize];
            A21[i][j] = A[i + newSize][j];
            A22[i][j] = A[i + newSize][j + newSize];

            B11[i][j] = B[i][j];
            B12[i][j] = B[i][j + newSize];
            B21[i][j] = B[i + newSize][j];
            B22[i][j] = B[i + newSize][j + newSize];
        }
    }
         
    // M1 = (A11 + A22) * (B11 + B22)
    addMatrices(A11, A22, temp1, newSize);
    add_ECP_Matrices(B11, B22, temp2, newSize); 
    strassen(temp1, temp2, M1, newSize, actual_size);

    // M2 = (A21 + A22) * B11
    addMatrices(A21, A22, temp1, newSize);
    strassen(temp1, B11, M2, newSize, actual_size);

    // M3 = A11 * (B12 - B22)
    subtract_ECP_Matrices(B12, B22, temp2, newSize);
    strassen(A11, temp2, M3, newSize, actual_size);

    // M4 = A22 * (B21 - B11)
    subtract_ECP_Matrices(B21, B11, temp2, newSize);
    strassen(A22, temp2, M4, newSize, actual_size);

    // M5 = (A11 + A12) * B22
    addMatrices(A11, A12, temp1, newSize);
    strassen(temp1, B22, M5, newSize, actual_size);

    // M6 = (A21 - A11) * (B11 + B12)
    subtractMatrices(A21, A11, temp1, newSize);
    add_ECP_Matrices(B11, B12, temp2, newSize);
    strassen(temp1, temp2, M6, newSize, actual_size);

    // M7 = (A12 - A22) * (B21 + B22)
    subtractMatrices(A12, A22, temp1, newSize);
    add_ECP_Matrices(B21, B22, temp2, newSize);
    strassen(temp1, temp2, M7, newSize, actual_size);

    // Calculating C11, C12, C21, C22
    add_ECP_Matrices(M1, M4, temp2, newSize);
    subtract_ECP_Matrices(temp2, M5, temp3, newSize);
    add_ECP_Matrices(temp3, M7, C11, newSize);

    add_ECP_Matrices(M3, M5, C12, newSize);

    add_ECP_Matrices(M2, M4, C21, newSize);

    add_ECP_Matrices(M1, M3, temp2, newSize);
    subtract_ECP_Matrices(temp2, M2, temp3, newSize);
    add_ECP_Matrices(temp3, M6, C22, newSize);
    
    // Combining C11, C12, C21, C22 into C
    for (int i = 0; i < newSize; i++) {
        for (int j = 0; j < newSize; j++) {
            C[i][j] = C11[i][j];
            C[i][j + newSize] = C12[i][j];
            C[i + newSize][j] = C21[i][j];
            C[i + newSize][j + newSize] = C22[i][j];
        }
    }

    // Freeing allocated memory
    freeMatrix(A11, newSize);
    freeMatrix(A12, newSize);
    freeMatrix(A21, newSize);
    freeMatrix(A22, newSize);
    free_ECP_Matrix(B11, newSize);
    free_ECP_Matrix(B12, newSize);
    free_ECP_Matrix(B21, newSize);
    free_ECP_Matrix(B22, newSize);
    free_ECP_Matrix(C11, newSize);
    free_ECP_Matrix(C12, newSize);
    free_ECP_Matrix(C21, newSize);
    free_ECP_Matrix(C22, newSize);
    free_ECP_Matrix(M1, newSize);
    free_ECP_Matrix(M2, newSize);
    free_ECP_Matrix(M3, newSize);
    free_ECP_Matrix(M4, newSize);
    free_ECP_Matrix(M5, newSize);
    free_ECP_Matrix(M6, newSize);
    free_ECP_Matrix(M7, newSize);
    freeMatrix(temp1, newSize);
    free_ECP_Matrix(temp2, newSize);
    free_ECP_Matrix(temp3, newSize);
}

// --------------------------------------------------------- Helper Functions End --------------------------------------------------------------------------

// Main function to test the algorithm
int main() {

    BIG_256_56 CURVE_Order_NIST256_minus_1;
 
    // Copy BIG to another BIG (Constant Time) - extern void BIG_256_56_copy(BIG_256_56 x, BIG_256_56 y); - y is the num to be copied 
    BIG_256_56_copy(CURVE_Order_NIST256_minus_1/* BIG_256_56 x */, CURVE_Order_NIST256/* BIG_256_56 y */);

    // Decrement BIG by a small integer - output not normalised (Constant Time) - extern void BIG_256_56_dec(BIG_256_56 x, int i); 
    BIG_256_56_dec(CURVE_Order_NIST256_minus_1/* BIG_256_56 x */, 1/* int i */);
    
    BIG_256_56 max_256_bit_number, modulus_for_random_num_gen;

    // Copy BIG to another BIG (Constant Time) - extern void BIG_256_56_copy(BIG_256_56 x, BIG_256_56 y); - y is the num to be copied 
    BIG_256_56_copy(max_256_bit_number/* BIG_256_56 x */, CURVE_Order_NIST256_minus_1/* BIG_256_56 y */);
    BIG_256_56_copy(modulus_for_random_num_gen/* BIG_256_56 x */, CURVE_Order_NIST256/* BIG_256_56 y */);

    // Below code w.r.t. CSPRNG is inspired from testecc.c
    int i;
    unsigned long ran;
    char raw[100];
    
    octet RAW = {0, sizeof(raw), raw};
    csprng RNG;                // Crypto Strong RNG
    
    time((time_t *)&ran);
    RAW.len = 100;              // fake random seed source
    RAW.val[0] = ran;
    RAW.val[1] = ran >> 8;
    RAW.val[2] = ran >> 16;
    RAW.val[3] = ran >> 24;
    
    for (i = 4; i < 100; i++) RAW.val[i] = i;
    
    CREATE_CSPRNG(&RNG, &RAW);  // initialise strong RNG
    
    ECP_NIST256 point_generator, random_point;

#if CURVETYPE_NIST256==MONTGOMERY
    // Set ECP to point(x,[y]) given x - extern int ECP_NIST256_set(ECP_NIST256 *P, BIG_256_56 x);
#else  
    // Set ECP to point(x,y) given x and y - ECP_NIST256_set(ECP_NIST256 *P, BIG_256_56 x, BIG_256_56 y)
    ECP_NIST256_set(&point_generator, CURVE_Gx_NIST256, CURVE_Gy_NIST256);
#endif
   
// In this test, size refers to matrix size and BIT_WIDTH_EACH_PLAINTEXT refers to the bitsize of each element of the matrix
// Change "size" variable in the below line accordingly to profile across specific matrix sizes.
for (int size = 8; size <= 512; size *= 2) {
    // Note that to populate various matrix size measurements corresponding to each bit size in separate excel files, file should be opened in "w" mode when profiling the first matrix size corresponding to a particular bit size and then, in rest of the profiling runs for that particular bit size, file should be opened in "a" mode
    FILE *csv_output_file;
    if (size == 8) {
        if (BIT_WIDTH_EACH_PLAINTEXT == 4) {
            csv_output_file = fopen("mm_strassen_size04.csv", "w");
        } else if (BIT_WIDTH_EACH_PLAINTEXT == 8) {
            csv_output_file = fopen("mm_strassen_size08.csv", "w");
        } else if (BIT_WIDTH_EACH_PLAINTEXT == 12) {
            csv_output_file = fopen("mm_strassen_size12.csv", "w");
        } else if (BIT_WIDTH_EACH_PLAINTEXT == 16) {
            csv_output_file = fopen("mm_strassen_size16.csv", "w");
        }
    } else {
        if (BIT_WIDTH_EACH_PLAINTEXT == 4) {
            csv_output_file = fopen("mm_strassen_size04.csv", "a");
        } else if (BIT_WIDTH_EACH_PLAINTEXT == 8) {
            csv_output_file = fopen("mm_strassen_size08.csv", "a");
        } else if (BIT_WIDTH_EACH_PLAINTEXT == 12) {
            csv_output_file = fopen("mm_strassen_size12.csv", "a");
        } else if (BIT_WIDTH_EACH_PLAINTEXT == 16) {
            csv_output_file = fopen("mm_strassen_size16.csv", "a");
        }
    }

    if (csv_output_file == NULL) {
        printf("Error opening file!\n");
        return 1;
    }
    
    // Since all the matrix size measurements corresponding to a particular bit size will be in one (i.e., same) excel file, header need to be written for only the first matrix size profiling run corresponding to the particular bit size
    if (size == 8) {
        // Write header
        fprintf(csv_output_file, "matrix size, bit size, schoolbook (us), strassen (us), Schoolbook / Strassen\n");
    } 

    unsigned int test_num;
    int num_successful_tests = 0;
#ifdef TIME_PROFILE
    struct timespec time_start, time_end;
    long time_elapsed;
    long time_TOTAL_schoolbook = 0;
    long time_TOTAL_strassen = 0;
#endif
   
    for (test_num = 0; test_num < NUM_TESTS; test_num++)
    {        
        BIG_256_56 random_num;
        
        srand(time(0) + test_num); // Seed the PRNG for rand()

//        printf("\nTest: %d\n", test_num);
        
        int** A = allocateMatrix(size);
        ECP_NIST256** B = allocate_ECP_Matrix(size);
        ECP_NIST256** B1 = allocate_ECP_Matrix(size);

// -----------------------------------------------------------------------------------------------------------------------------------------------------------
// -------------------------------------- Populating input plaintext matrix (A) and two ECP point matrices (B, B1) begin -------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------------------------

        int mask;
        if (BIT_WIDTH_EACH_PLAINTEXT == 4) {
            mask = 0xF;
        } else if (BIT_WIDTH_EACH_PLAINTEXT == 8) {
            mask = 0xFF;
        } else if (BIT_WIDTH_EACH_PLAINTEXT == 12) {
            mask = 0xFFF;
        } else if (BIT_WIDTH_EACH_PLAINTEXT == 16) {
            mask = 0xFFFF;
        }
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                // Populating input vector of 'BIT_WIDTH_EACH_PLAINTEXT'-bit unsigned numbers
                A[i][j] = rand() & mask;
              
                // Create an unbiased random BIG from a random number generator, reduced with respect to a modulus (Constant Time as used) - extern void BIG_256_56_randomnum(BIG_256_56 x, BIG_256_56 n, csprng *r);
                BIG_256_56_randomnum(random_num/* BIG_256_56 x */, modulus_for_random_num_gen/* BIG_256_56 n */, &RNG/* csprng *r */);
                random_point = point_generator;
                // Multiplies an ECP instance P by a BIG, side-channel resistant - extern void ECP_NIST256_clmul(ECP_NIST256 *P, BIG_256_56 e, BIG_256_56 maxe);
                ECP_NIST256_clmul(&random_point /* ECP_NIST256 *P */, random_num /* BIG_256_56 e */, max_256_bit_number/* BIG_256_56 maxe */);
                B[i][j] = random_point;
                
                BIG_256_56_randomnum(random_num/* BIG_256_56 x */, modulus_for_random_num_gen/* BIG_256_56 n */, &RNG/* csprng *r */);
                random_point = point_generator;
                // Multiplies an ECP instance P by a BIG, side-channel resistant - extern void ECP_NIST256_clmul(ECP_NIST256 *P, BIG_256_56 e, BIG_256_56 maxe);
                ECP_NIST256_clmul(&random_point /* ECP_NIST256 *P */, random_num /* BIG_256_56 e */, max_256_bit_number/* BIG_256_56 maxe */);
                B1[i][j] = random_point;
            }
        }

// -------------------------------------- Populating input plaintext matrix (A) and two ECP point matrices (B, B1) end -------------------------------------

// -----------------------------------------------------------------------------------------------------------------------------------------------------------
// ------------------------------ Profiling timing for plaintext-ciphertext matrix multiplication using Schoolbook Approach begin ----------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------------------------

#ifdef TIME_PROFILE
//        printf("Schoolbook Approach (ns): ");
        clock_gettime (CLOCK_PROCESS_CPUTIME_ID, &time_start);
#endif
        ECP_NIST256** result_schoolbook1 = allocate_ECP_Matrix(size);
        ECP_NIST256** result_schoolbook2 = allocate_ECP_Matrix(size);
       
        ECP_NIST256 schoolbook_temp1, schoolbook_temp2;
        schoolbook_temp1 = point_generator;

        // Calculating output matrix for first ECP point matrix in schoolbook approach.
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                // Set ECP to point-at-infinity - extern void ECP_ZZZ_inf(ECP_ZZZ *P);
                ECP_NIST256_inf(&schoolbook_temp2/* ECP_ZZZ *P */);
                                
                for (int k = 0; k < size; k++) {
                    schoolbook_temp1 = B[k][j];
                    // Multiplies an ECP instance P by a small integer, side-channel resistant - extern void ECP_NIST256_pinmul(ECP_NIST256 *P, int i, int b);
                    ECP_NIST256_pinmul(&schoolbook_temp1/* ECP_NIST256 *P */, A[i][k]/* int i */, (BIT_WIDTH_EACH_PLAINTEXT)/* int b */); // "(BIT_WIDTH_EACH_PLAINTEXT + 1)" bits considering (2 ** BIT_WIDTH_EACH_PLAINTEXT) although (2 ** BIT_WIDTH_EACH_PLAINTEXT) can not be the random number generated.
                    // Adds ECP instance Q to ECP instance P - extern void ECP_NIST256_add(ECP_NIST256 *P, ECP_NIST256 *Q);
                    ECP_NIST256_add(&schoolbook_temp2/* ECP_NIST256 *P */, &schoolbook_temp1/* ECP_NIST256 *Q */);                       
                }
                result_schoolbook1[i][j] = schoolbook_temp2;
            }
        }
        // Calculating output matrix for second ECP point matrix in schoolbook approach.
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                // Set ECP to point-at-infinity - extern void ECP_ZZZ_inf(ECP_ZZZ *P);
                ECP_NIST256_inf(&schoolbook_temp2/* ECP_ZZZ *P */);
                                
                for (int k = 0; k < size; k++) {
                    schoolbook_temp1 = B1[k][j];
                    // Multiplies an ECP instance P by a small integer, side-channel resistant - extern void ECP_NIST256_pinmul(ECP_NIST256 *P, int i, int b);
                    ECP_NIST256_pinmul(&schoolbook_temp1/* ECP_NIST256 *P */, A[i][k]/* int i */, (BIT_WIDTH_EACH_PLAINTEXT)/* int b */); // "(BIT_WIDTH_EACH_PLAINTEXT + 1)" bits considering (2 ** BIT_WIDTH_EACH_PLAINTEXT) although (2 ** BIT_WIDTH_EACH_PLAINTEXT) can not be the random number generated.
                    // Adds ECP instance Q to ECP instance P - extern void ECP_NIST256_add(ECP_NIST256 *P, ECP_NIST256 *Q);
                    ECP_NIST256_add(&schoolbook_temp2/* ECP_NIST256 *P */, &schoolbook_temp1/* ECP_NIST256 *Q */);                       
                }
                result_schoolbook2[i][j] = schoolbook_temp2;
            }
        } 
#ifdef TIME_PROFILE
        clock_gettime (CLOCK_PROCESS_CPUTIME_ID, &time_end);
        time_elapsed = timespec_diff (time_start, time_end); // time_elapsed = (time_end.tv_nsec - time_start.tv_nsec);
   
//        printf("%ld\n", time_elapsed);
        time_TOTAL_schoolbook += time_elapsed;
#endif

// ------------------------------ Profiling timing for plaintext-ciphertext matrix multiplication using Schoolbook Approach end ----------------------------


// -----------------------------------------------------------------------------------------------------------------------------------------------------------
// -------------------------------- Profiling timing for plaintext-ciphertext matrix multiplication using Strassen Approach begin ----------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------------------------

#ifdef TIME_PROFILE
//        printf("Strassen Approach (ns): ");
        clock_gettime (CLOCK_PROCESS_CPUTIME_ID, &time_start);
#endif
        ECP_NIST256** result_strassen1 = allocate_ECP_Matrix(size);
        ECP_NIST256** result_strassen2 = allocate_ECP_Matrix(size);
        // Perform plaintext-ciphertext matrix multiplication using Strassen algorithm for first ECP point matrix
        strassen(A, B, result_strassen1, size, size);

        // Perform plaintext-ciphertext matrix multiplication using Strassen algorithm for Second ECP point matrix
        strassen(A, B1, result_strassen2, size, size);
#ifdef TIME_PROFILE
        clock_gettime (CLOCK_PROCESS_CPUTIME_ID, &time_end);
        time_elapsed = timespec_diff (time_start, time_end); // time_elapsed = (time_end.tv_nsec - time_start.tv_nsec);
   
//        printf("%ld\n", time_elapsed);
        time_TOTAL_strassen += time_elapsed;
#endif    
// -------------------------------- Profiling timing for plaintext-ciphertext matrix multiplication using Strassen Approach end ----------------------------

// ------------------------------------------------------ Post Processing begin ----------------------------------------------------------------------------

        int count_correct_output = 0;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                // Tests for equality of two ECPs - return 1 if P=Q, else returns 0 - extern int ECP_NIST256_equals(ECP_NIST256 *P, ECP_NIST256 *Q);
                if(ECP_NIST256_equals(&(result_schoolbook1[i][j])/* ECP_NIST256 *P */, &(result_strassen1[i][j])/* ECP_NIST256 *Q */) == 1) {  
                    if (ECP_NIST256_equals(&(result_schoolbook2[i][j])/* ECP_NIST256 *P */, &(result_strassen2[i][j])/* ECP_NIST256 *Q */) == 1) {
                        count_correct_output += 1;  
                    }               
                }
            }
        }
        
        if (count_correct_output == (size * size)) {
            printf("\nTest %d is successful for matrix size of %d\n", test_num, size);
            num_successful_tests += 1;
        } else {
            printf("\nTest %d is not successful for matrix size of %d\n", test_num, size);
        }
    
        // Free allocated memory
        freeMatrix(A, size);
        free_ECP_Matrix(B, size);
        free_ECP_Matrix(B1, size);
        free_ECP_Matrix(result_schoolbook1, size);
        free_ECP_Matrix(result_schoolbook2, size);
        free_ECP_Matrix(result_strassen1, size);
        free_ECP_Matrix(result_strassen2, size);
    }
    // So, comparing "num_successful_tests" with "(NUM_TESTS)" in below part of the code
    if (num_successful_tests == (NUM_TESTS)) {
        printf("\nAll %d tests successful for matrix size of %d!\n", (NUM_TESTS), size);
    }

#ifdef TIME_PROFILE
    time_TOTAL_schoolbook /= NUM_TESTS;
    time_TOTAL_strassen /= NUM_TESTS;
#endif

    printf("Average Time (us) for schoolbook approach: %.3f\n", (time_TOTAL_schoolbook*1e-3));
    printf("Average Time (us) for strassen approach: %.3f\n", time_TOTAL_strassen*1e-3);
    printf("Timing in schoolbook approach / Timing in strassen approach: %.3f\n", (time_TOTAL_schoolbook*1e-3)/(time_TOTAL_strassen*1e-3));
    // Populate each row corresponding to each vector length
    fprintf(csv_output_file, "%d, %d, %.3f, %.3f, %.3f\n", size, BIT_WIDTH_EACH_PLAINTEXT, (time_TOTAL_schoolbook*1e-3), (time_TOTAL_strassen*1e-3), ((time_TOTAL_schoolbook*1e-3)/(time_TOTAL_strassen*1e-3)));
    // Close the file
    fclose(csv_output_file);

// ------------------------------------------------------ Post Processing end ----------------------------------------------------------------------------

  }
    
  printf("Data written to csv file\n");
  return 0;
}
