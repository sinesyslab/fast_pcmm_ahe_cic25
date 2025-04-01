// This test is to profile plaintext-ciphertext matrix multiplication using Schoolbook and Proposed spproaches (the Proposed approach is inspired by Cussen's algorithm and extends Cussen's algorithm from plaintext to encrypted setting)

// -------------------------------------------------------------------- Test Begins -----------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h> // To measure time and do profiling. Required for also random number generation.
#include <math.h> // To use "round" function for rounding the vector lengths after each iteration to nearest integer.
#include "randapi.h" // For CREATE_CSPRNG function call

#include "test_pcmm_proposed.h"

// Include NIST256 curve specific header files
#include "ecp_NIST256.h"
#include "fp_NIST256.h"

// To profile for "entire cussen" scenario, comment both the lines below. To profile for "only forward path" / "only backward path" scenarios, uncomment only the corresponding define accordingly
//#define ONLY_FORWARD_PATH // Define to time profile "only forward path of cussen's algorithm"
//#define ONLY_BACKWARD_PATH // Define to time profile "only backward path of cussen's algorithm"

#define TIME_PROFILE

// --------------------------------------------------------- Helper Functions Begin --------------------------------------------------------------------------

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

// Function to free memory of an ECP matrix
void free_ECP_Matrix(ECP_NIST256** matrix, int size) {
    for (int i = 0; i < size; i++) {
        free(matrix[i]);
    }
    free(matrix);
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

void copy_matrix(ECP_NIST256** A, ECP_NIST256** B, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            A[i][j] = B[i][j];
        }
    }
}

// Comparison function
// The comparison function’s return value informs qsort how to order two elements.
// A negative return value means the first element precedes the second.
// Zero means the elements are equal.
// A positive return value means the first element follows the second.
int compare(const void* a, const void* b)
{
    return (*(int*)a - *(int*)b); // For sorting in ascending order
}

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

// --------------------------------------------------------- Helper Functions End --------------------------------------------------------------------------

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

// In this test, "num_elements" refers to matrix size and BIT_WIDTH_EACH_PLAINTEXT refers to bitsize of elements of the plaintext matrix
// Change "num_elements" variable in the below line accordingly to profile across specific matrix sizes
 for (int num_elements = 8; num_elements <= 512; (num_elements *= 2)) {
    // Note that to populate various matrix size measurements corresponding to each bit size in separate excel files, file should be opened in "w" mode when profiling the first matrix size corresponding to a particular bit size and then, in rest of the profiling runs for that particular bit size, file should be opened in "a" mode
    FILE *csv_output_file;
    if (num_elements == 8) {
        if (BIT_WIDTH_EACH_PLAINTEXT == 4) {
            csv_output_file = fopen("mm_cussen_size04.csv", "w");
        } else if (BIT_WIDTH_EACH_PLAINTEXT == 8) {
            csv_output_file = fopen("mm_cussen_size08.csv", "w");
        } else if (BIT_WIDTH_EACH_PLAINTEXT == 12) {
            csv_output_file = fopen("mm_cussen_size12.csv", "w");
        } else if (BIT_WIDTH_EACH_PLAINTEXT == 16) {
            csv_output_file = fopen("mm_cussen_size16.csv", "w");
        }
    } else {
        if (BIT_WIDTH_EACH_PLAINTEXT == 4) {
            csv_output_file = fopen("mm_cussen_size04.csv", "a");
        } else if (BIT_WIDTH_EACH_PLAINTEXT == 8) {
            csv_output_file = fopen("mm_cussen_size08.csv", "a");
        } else if (BIT_WIDTH_EACH_PLAINTEXT == 12) {
            csv_output_file = fopen("mm_cussen_size12.csv", "a");
        } else if (BIT_WIDTH_EACH_PLAINTEXT == 16) {
            csv_output_file = fopen("mm_cussen_size16.csv", "a");
        }
    }

    if (csv_output_file == NULL) {
        printf("Error opening file!\n");
        return 1;
    }

#if !defined(ONLY_FORWARD_PATH) && !defined(ONLY_BACKWARD_PATH)
  // Since all the matrix size measurements corresponding to a particular bit size will be in one (i.e., same) excel file, header need to be written for only the first matrix size profiling run corresponding to the particular bit size
    if (num_elements == 8) {
        // Write header
        fprintf(csv_output_file, "matrix size, bit size, schoolbook (us), CUssen (us), Schoolbook / Cussen\n");
    }
#else
    if (num_elements == 8) {
        // Write header
        fprintf(csv_output_file, "Matrix size, bit size, Cussen (us)\n");
    }
#endif

    ECP_NIST256 point_generator, random_point;

#if CURVETYPE_NIST256==MONTGOMERY
    // Set ECP to point(x,[y]) given x - extern int ECP_NIST256_set(ECP_NIST256 *P, BIG_256_56 x);
#else  
    // Set ECP to point(x,y) given x and y - ECP_NIST256_set(ECP_NIST256 *P, BIG_256_56 x, BIG_256_56 y)
    ECP_NIST256_set(&point_generator, CURVE_Gx_NIST256, CURVE_Gy_NIST256);
#endif
   
    random_point = point_generator;

    unsigned int test_num;
    int num_successful_tests = 0;
#ifdef TIME_PROFILE
    struct timespec time_start, time_end;
    long time_TOTAL_schoolbook = 0;
    long time_elapsed;
    long time_TOTAL_cussen = 0;
#endif
    for (test_num = 0; test_num < NUM_TESTS; test_num++)
    {        
        BIG_256_56 random_num;
        
        srand(time(0) + test_num); // Seed the PRNG for rand()

//        printf("\nTest: %d\n", test_num);
        
        int** A = allocateMatrix(num_elements);

        ECP_NIST256** B = allocate_ECP_Matrix(num_elements);
        ECP_NIST256** B1 = allocate_ECP_Matrix(num_elements);

// -----------------------------------------------------------------------------------------------------------------------------------------------------------
// -------------------------------- Populating an input plaintext matrix (A) and two input ECP point matrices (B, B1) begin ----------------------------------
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
        for (int i = 0; i < num_elements; i++) {
            for (int j = 0; j < num_elements; j++) {
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

// -------------------------------- Populating an input plaintext matrix (A) and two input ECP point matrices (B, B1) end ----------------------------------

// -----------------------------------------------------------------------------------------------------------------------------------------------------------
// ----------------------------- Profiling timing for plaintext-ciphertext matrix multiplication using Schoolbook Approach begin -----------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------------------------

#if !defined(ONLY_FORWARD_PATH) && !defined(ONLY_BACKWARD_PATH)
// Profiling Schoolbook Approach timing
#ifdef TIME_PROFILE
//        printf("Schoolbook Approach (ns): ");
        clock_gettime (CLOCK_PROCESS_CPUTIME_ID, &time_start);
#endif
        ECP_NIST256** result1_schoolbook = allocate_ECP_Matrix(num_elements);
        ECP_NIST256** result2_schoolbook = allocate_ECP_Matrix(num_elements);
       
        ECP_NIST256 schoolbook_temp1, schoolbook_temp2;
        schoolbook_temp1 = point_generator;

        // Calculating output matrix for first ECP point matrix in schoolbook approach.
        for (int i = 0; i < num_elements; i++) {
            for (int j = 0; j < num_elements; j++) {
                // Set ECP to point-at-infinity - extern void ECP_ZZZ_inf(ECP_ZZZ *P);
                ECP_NIST256_inf(&schoolbook_temp2/* ECP_ZZZ *P */);
                                
                for (int k = 0; k < num_elements; k++) {
                    schoolbook_temp1 = B[k][j];
                    // Multiplies an ECP instance P by a small integer, side-channel resistant - extern void ECP_NIST256_pinmul(ECP_NIST256 *P, int i, int b);
                    ECP_NIST256_pinmul(&schoolbook_temp1/* ECP_NIST256 *P */, A[i][k]/* int i */, (BIT_WIDTH_EACH_PLAINTEXT)/* int b */); // "(BIT_WIDTH_EACH_PLAINTEXT + 1)" bits considering (2 ** BIT_WIDTH_EACH_PLAINTEXT) although (2 ** BIT_WIDTH_EACH_PLAINTEXT) can not be the random number generated.
                    // Adds ECP instance Q to ECP instance P - extern void ECP_NIST256_add(ECP_NIST256 *P, ECP_NIST256 *Q);
                    ECP_NIST256_add(&schoolbook_temp2/* ECP_NIST256 *P */, &schoolbook_temp1/* ECP_NIST256 *Q */);                       
                }
                result1_schoolbook[i][j] = schoolbook_temp2;
            }
        }

        // Calculating output matrix for Second ECP point matrix in schoolbook approach.
        for (int i = 0; i < num_elements; i++) {
            for (int j = 0; j < num_elements; j++) {
                // Set ECP to point-at-infinity - extern void ECP_ZZZ_inf(ECP_ZZZ *P);
                ECP_NIST256_inf(&schoolbook_temp2/* ECP_ZZZ *P */);
                                
                for (int k = 0; k < num_elements; k++) {
                    schoolbook_temp1 = B1[k][j];
                    // Multiplies an ECP instance P by a small integer, side-channel resistant - extern void ECP_NIST256_pinmul(ECP_NIST256 *P, int i, int b);
                    ECP_NIST256_pinmul(&schoolbook_temp1/* ECP_NIST256 *P */, A[i][k]/* int i */, (BIT_WIDTH_EACH_PLAINTEXT)/* int b */); // "(BIT_WIDTH_EACH_PLAINTEXT + 1)" bits considering (2 ** BIT_WIDTH_EACH_PLAINTEXT) although (2 ** BIT_WIDTH_EACH_PLAINTEXT) can not be the random number generated.
                    // Adds ECP instance Q to ECP instance P - extern void ECP_NIST256_add(ECP_NIST256 *P, ECP_NIST256 *Q);
                    ECP_NIST256_add(&schoolbook_temp2/* ECP_NIST256 *P */, &schoolbook_temp1/* ECP_NIST256 *Q */);                       
                }
                result2_schoolbook[i][j] = schoolbook_temp2;
            }
        }
       
#ifdef TIME_PROFILE
        clock_gettime (CLOCK_PROCESS_CPUTIME_ID, &time_end);
        time_elapsed = timespec_diff (time_start, time_end); // time_elapsed = (time_end.tv_nsec - time_start.tv_nsec);
   
//        printf("%ld\n", time_elapsed);
        time_TOTAL_schoolbook += time_elapsed;
#endif
#endif

// ----------------------------- Profiling timing for plaintext-ciphertext matrix multiplication using Schoolbook Approach end -----------------------------

// -----------------------------------------------------------------------------------------------------------------------------------------------------------
// -------------------------------- Profiling timing for plaintext-ciphertext matrix multiplication using Cussen Approach begin ------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------------------------

#ifndef ONLY_BACKWARD_PATH
// Profiling Cussen Approach timing
// Comment below 3 lines for time profiling of "only backward path of cussen algorithm". Uncomment them for "entire cussen algorithm" / "only forward path" time profiling.
#ifdef TIME_PROFILE
        clock_gettime (CLOCK_PROCESS_CPUTIME_ID, &time_start);
#endif
#endif

    int length_of_vector_after_iter[num_elements];
    int length_of_vector_before_elim_dupl_iter1[num_elements], length_of_vector_before_elim_dupl_iter2[num_elements], length_of_vector_before_elim_dupl_iter3[num_elements], length_of_vector_before_elim_dupl_iter4[num_elements];
    
    int** pointer_array_iter1 = allocateMatrix(num_elements);
    int** pointer_array_iter2 = allocateMatrix(num_elements);
    int** pointer_array_iter3 = allocateMatrix(num_elements);
    int** pointer_array_iter4 = allocateMatrix(num_elements);
            
    int** A_post_compression = allocateMatrix(num_elements);

    for (int i = 0; i < num_elements; i++) {
        length_of_vector_after_iter[i] = num_elements; // Initializing to 'num_elements' for the first iteration (for all the columns of the input plaintext matrix)
    }

    int length_of_vector_before_elim_dupl;

// -----------------------------------------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------- Compression Phase begin -----------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------------------------

// Cussen's algorithm (No alignment scenario and four iterations) 
// Forward path - Taking differences, sorting, eliminating duplicates and populating tracking pointers array  
    for (int col_index = 0; col_index < num_elements; col_index ++) {
        int input_vector[num_elements];

        for (int i = 0; i < num_elements; i++) {
            input_vector[i] = A[i][col_index];
        }

        int input_vector_after_diff[num_elements];
        int input_vector_after_diff_copy[num_elements];
        int vector_after_iter[num_elements];
       
        for (int iteration = 0; iteration < NUM_ITERATIONS; iteration ++) {
        // Start of Iteration (Forward path)            
            if (iteration != 0) {
                input_vector_after_diff[0] = input_vector[0];
                for (int index = 1; index < length_of_vector_after_iter[col_index]; index ++) {
                    input_vector_after_diff[index] = input_vector[index] - input_vector[index - 1];
                }
            } else {
                for (int index = 0; index < length_of_vector_after_iter[col_index]; index ++) {
                    input_vector_after_diff[index] = input_vector[index];
                }
            }
   
            for(int index = 0; index < length_of_vector_after_iter[col_index]; index ++) {
                input_vector_after_diff_copy[index] = input_vector_after_diff[index];
            }

            // Quick sort function in C - void qsort(void* base, size_t num, size_t size, int (*compare)(const void*, const void*));
            qsort(input_vector_after_diff, length_of_vector_after_iter[col_index], sizeof(unsigned int), compare);
   
            if (iteration == 0) {
                length_of_vector_before_elim_dupl_iter1[col_index] = length_of_vector_after_iter[col_index];
                length_of_vector_before_elim_dupl = length_of_vector_before_elim_dupl_iter1[col_index];
            } else if (iteration == 1) {
                length_of_vector_before_elim_dupl_iter2[col_index] = length_of_vector_after_iter[col_index];
                length_of_vector_before_elim_dupl = length_of_vector_before_elim_dupl_iter2[col_index];
            } else if (iteration == 2) {
                length_of_vector_before_elim_dupl_iter3[col_index] = length_of_vector_after_iter[col_index];
                length_of_vector_before_elim_dupl = length_of_vector_before_elim_dupl_iter3[col_index];
            } else if (iteration == 3) {
                length_of_vector_before_elim_dupl_iter4[col_index] = length_of_vector_after_iter[col_index];
                length_of_vector_before_elim_dupl = length_of_vector_before_elim_dupl_iter4[col_index];
            }
   
            // Eliminating Duplicates
            int index_vector_after_iter = 0;
            for (int index = 0; index < length_of_vector_after_iter[col_index]; index++) {
                if (index == 0) {
                    vector_after_iter[0] = input_vector_after_diff[0];
                    index_vector_after_iter = 1;
                } else if (vector_after_iter[index_vector_after_iter - 1] != input_vector_after_diff[index]) {
                    vector_after_iter[index_vector_after_iter ++] = input_vector_after_diff[index];
                }
            }
           
            length_of_vector_after_iter[col_index] = index_vector_after_iter;
            
            // Vector of pointers for tracking
            for (int i = 0; i < length_of_vector_before_elim_dupl; i++) {
                for (int j = 0; j < length_of_vector_after_iter[col_index]; j++) {
                    if (input_vector_after_diff_copy[i] == vector_after_iter[j]) {
                        if (iteration == 0) {
                            pointer_array_iter1[i][col_index] = j;
                        } else if (iteration == 1) {
                            pointer_array_iter2[i][col_index] = j;
                        } else if (iteration == 2) {
                            pointer_array_iter3[i][col_index] = j;
                        } else if (iteration == 3) {
                            pointer_array_iter4[i][col_index] = j;
                        }
                        break;
                    }
                }
            }

            // Update the input vector for the next iteration
            for (int index = 0; index < length_of_vector_after_iter[col_index]; index ++) {
                input_vector[index] = vector_after_iter[index];
            }
        // End of Iteration (Forward path)
        }
        for (int index = 0; index < length_of_vector_after_iter[col_index]; index++) {
            A_post_compression[index][col_index] = input_vector[index];
        }
     }

// ---------------------------------------------------------------- Compression Phase end -----------------------------------------------------------------

#ifdef ONLY_FORWARD_PATH
// Uncomment below 3 lines for "only forward path" time profiling of cussen algorithm. Comment them for "entire cussen algorithm" / "only backward path" time profiling.
#ifdef TIME_PROFILE
        clock_gettime (CLOCK_PROCESS_CPUTIME_ID, &time_end);
#endif
#endif

// -----------------------------------------------------------------------------------------------------------------------------------------------------------

// -----------------------------------------------------------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------ Reconstruction phase begin -----------------------------------------------------  
// -----------------------------------------------------------------------------------------------------------------------------------------------------------

#ifdef ONLY_BACKWARD_PATH
// Comment below 3 lines for "entire cussen algorithm" / "only forward path" time profiling. Uncomment them for "only backward path" time profiling
#ifdef TIME_PROFILE
        clock_gettime (CLOCK_PROCESS_CPUTIME_ID, &time_start);
#endif
#endif
    ECP_NIST256** result1_cussen = allocate_ECP_Matrix(num_elements);
    ECP_NIST256** result2_cussen = allocate_ECP_Matrix(num_elements);
    ECP_NIST256** result_outer_product = allocate_ECP_Matrix(num_elements);
    
    ECP_NIST256 output_vector_cussen_approach_post_iter1[num_elements];
    ECP_NIST256 output_vector_cussen_approach_post_iter2[num_elements];
    ECP_NIST256 output_vector_cussen_approach_post_iter3[num_elements];
    ECP_NIST256 output_vector_cussen_approach_post_iter4[num_elements];

    ECP_NIST256 output_vector_cussen_approach[num_elements];

// -----------------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------- Reconstruction phase for first ECP point matrix in the ciphertext begin -------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------------------------

    for (int col_index = 0; col_index < num_elements; col_index ++) {
        for (int col_ct_matrix = 0; col_ct_matrix < num_elements; col_ct_matrix ++) {
            for (int index = 0; index < length_of_vector_after_iter[col_index]; index ++) {
                output_vector_cussen_approach[index] = B[col_index][col_ct_matrix];
                // Multiplies an ECP instance P by a small integer, side-channel resistant - extern void ECP_NIST256_pinmul(ECP_NIST256 *P, int i, int b);
                ECP_NIST256_pinmul(&output_vector_cussen_approach[index]/* ECP_NIST256 *P */, A_post_compression[index][col_index]/* int i */, (BIT_WIDTH_EACH_PLAINTEXT)/* int b */); // "(BIT_WIDTH_EACH_PLAINTEXT + 1)" bits considering (2 ** BIT_WIDTH_EACH_PLAINTEXT) although (2 ** BIT_WIDTH_EACH_PLAINTEXT) can not be the random number generated.
            }
            
// Backward path - Point Multiplication, Copying duplicates using pointers into appropriate places, Point Additions.        
            for (int iteration = (NUM_ITERATIONS - 1); iteration >= 0; iteration --) {
            // Start of Iteration (Backward Path)
                if (iteration == 0) {
                    length_of_vector_before_elim_dupl = length_of_vector_before_elim_dupl_iter1[col_index];
                } else if (iteration == 1) {
                    length_of_vector_before_elim_dupl = length_of_vector_before_elim_dupl_iter2[col_index];
                } else if (iteration == 2) {
                    length_of_vector_before_elim_dupl = length_of_vector_before_elim_dupl_iter3[col_index];
                } else if (iteration == 3) {
                    length_of_vector_before_elim_dupl = length_of_vector_before_elim_dupl_iter4[col_index];
                }
               
                for (int index = 0; index < length_of_vector_before_elim_dupl; index ++) {
                    if (iteration == 0) {
                        output_vector_cussen_approach_post_iter1[index] = output_vector_cussen_approach[(pointer_array_iter1[index][col_index])];
                    } else if (iteration == 1) {
                        output_vector_cussen_approach_post_iter2[index] = output_vector_cussen_approach[(pointer_array_iter2[index][col_index])];
                    } else if (iteration == 2) {
                        output_vector_cussen_approach_post_iter3[index] = output_vector_cussen_approach[(pointer_array_iter3[index][col_index])];
                    } else if (iteration == 3) {
                        output_vector_cussen_approach_post_iter4[index] = output_vector_cussen_approach[(pointer_array_iter4[index][col_index])];
                    }
                }
               
                if (iteration != 0) {
                   // Point Additions
                   for (int index = 1; index < length_of_vector_before_elim_dupl; index ++) {
                      // Adds ECP instance Q to ECP instance P - extern void ECP_NIST256_add(ECP_NIST256 *P, ECP_NIST256 *Q);
                      if (iteration == 3) {
                          ECP_NIST256_add(&output_vector_cussen_approach_post_iter4[index]/* ECP_NIST256 *P */, &output_vector_cussen_approach_post_iter4[index - 1]/* ECP_NIST256 *Q */);
                      } else if (iteration == 2) {
                          ECP_NIST256_add(&output_vector_cussen_approach_post_iter3[index]/* ECP_NIST256 *P */, &output_vector_cussen_approach_post_iter3[index - 1]/* ECP_NIST256 *Q */);
                      } else if (iteration == 1) {
                          ECP_NIST256_add(&output_vector_cussen_approach_post_iter2[index]/* ECP_NIST256 *P */, &output_vector_cussen_approach_post_iter2[index - 1]/* ECP_NIST256 *Q */);
                      }
                   }
                }
             
                for (int index = 0; index < length_of_vector_before_elim_dupl; index ++) {
                   if (iteration == 3) {
                       output_vector_cussen_approach[index] = output_vector_cussen_approach_post_iter4[index];
                   } else if (iteration == 2) {
                       output_vector_cussen_approach[index] = output_vector_cussen_approach_post_iter3[index];
                   } else if (iteration == 1) {
                       output_vector_cussen_approach[index] = output_vector_cussen_approach_post_iter2[index];
                   } else if (iteration == 0) {
                       output_vector_cussen_approach[index] = output_vector_cussen_approach_post_iter1[index];
                   }
                }
            // End of iteration (Backward Path)                                
            } 
            // Write into corresponding column of the outer product
            for (int index = 0; index < num_elements; index ++) {
                result_outer_product[index][col_ct_matrix] = output_vector_cussen_approach[index];
            }
          }
         if (col_index != 0) {
             add_ECP_Matrices(result1_cussen/* ECP_NIST256** A */, result_outer_product/* ECP_NIST256** B */, result1_cussen/* ECP_NIST256** C */, num_elements/* int size */);
         } else {
             copy_matrix(result1_cussen/* ECP_NIST256** A */, result_outer_product/* ECP_NIST256** B */, num_elements/* int size */); // A = B - copy all matrix elements
         }  
      } 

// --------------------------------------------- Reconstruction phase for first ECP point matrix in the ciphertext end -------------------------------------

// -----------------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------- Reconstruction phase for second ECP point matrix in the ciphertext begin ------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------------------------

    for (int col_index = 0; col_index < num_elements; col_index ++) {
        for (int col_ct_matrix = 0; col_ct_matrix < num_elements; col_ct_matrix ++) {
            for (int index = 0; index < length_of_vector_after_iter[col_index]; index ++) {
                output_vector_cussen_approach[index] = B1[col_index][col_ct_matrix];
                // Multiplies an ECP instance P by a small integer, side-channel resistant - extern void ECP_NIST256_pinmul(ECP_NIST256 *P, int i, int b);
                ECP_NIST256_pinmul(&output_vector_cussen_approach[index]/* ECP_NIST256 *P */, A_post_compression[index][col_index]/* int i */, (BIT_WIDTH_EACH_PLAINTEXT)/* int b */); // "(BIT_WIDTH_EACH_PLAINTEXT + 1)" bits considering (2 ** BIT_WIDTH_EACH_PLAINTEXT) although (2 ** BIT_WIDTH_EACH_PLAINTEXT) can not be the random number generated.
            }
            
// Backward path - Point Multiplication, Copying duplicates using pointers into appropriate places, Point Additions.        
            for (int iteration = (NUM_ITERATIONS - 1); iteration >= 0; iteration --) {
            // Start of Iteration (Backward Path)
                if (iteration == 0) {
                    length_of_vector_before_elim_dupl = length_of_vector_before_elim_dupl_iter1[col_index];
                } else if (iteration == 1) {
                    length_of_vector_before_elim_dupl = length_of_vector_before_elim_dupl_iter2[col_index];
                } else if (iteration == 2) {
                    length_of_vector_before_elim_dupl = length_of_vector_before_elim_dupl_iter3[col_index];
                } else if (iteration == 3) {
                    length_of_vector_before_elim_dupl = length_of_vector_before_elim_dupl_iter4[col_index];
                }
               
                for (int index = 0; index < length_of_vector_before_elim_dupl; index ++) {
                    if (iteration == 0) {
                        output_vector_cussen_approach_post_iter1[index] = output_vector_cussen_approach[(pointer_array_iter1[index][col_index])];
                    } else if (iteration == 1) {
                        output_vector_cussen_approach_post_iter2[index] = output_vector_cussen_approach[(pointer_array_iter2[index][col_index])];
                    } else if (iteration == 2) {
                        output_vector_cussen_approach_post_iter3[index] = output_vector_cussen_approach[(pointer_array_iter3[index][col_index])];
                    } else if (iteration == 3) {
                        output_vector_cussen_approach_post_iter4[index] = output_vector_cussen_approach[(pointer_array_iter4[index][col_index])];
                    }
                }
               
                if (iteration != 0) {
                   // Point Additions
                   for (int index = 1; index < length_of_vector_before_elim_dupl; index ++) {
                      // Adds ECP instance Q to ECP instance P - extern void ECP_NIST256_add(ECP_NIST256 *P, ECP_NIST256 *Q);
                      if (iteration == 3) {
                          ECP_NIST256_add(&output_vector_cussen_approach_post_iter4[index]/* ECP_NIST256 *P */, &output_vector_cussen_approach_post_iter4[index - 1]/* ECP_NIST256 *Q */);
                      } else if (iteration == 2) {
                          ECP_NIST256_add(&output_vector_cussen_approach_post_iter3[index]/* ECP_NIST256 *P */, &output_vector_cussen_approach_post_iter3[index - 1]/* ECP_NIST256 *Q */);
                      } else if (iteration == 1) {
                          ECP_NIST256_add(&output_vector_cussen_approach_post_iter2[index]/* ECP_NIST256 *P */, &output_vector_cussen_approach_post_iter2[index - 1]/* ECP_NIST256 *Q */);
                      }
                   }
                }
             
                for (int index = 0; index < length_of_vector_before_elim_dupl; index ++) {
                   if (iteration == 3) {
                       output_vector_cussen_approach[index] = output_vector_cussen_approach_post_iter4[index];
                   } else if (iteration == 2) {
                       output_vector_cussen_approach[index] = output_vector_cussen_approach_post_iter3[index];
                   } else if (iteration == 1) {
                       output_vector_cussen_approach[index] = output_vector_cussen_approach_post_iter2[index];
                   } else if (iteration == 0) {
                       output_vector_cussen_approach[index] = output_vector_cussen_approach_post_iter1[index];
                   }
                }
            // End of iteration (Backward Path)                                
            } 
            // Write into corresponding column of the outer product
            for (int index = 0; index < num_elements; index ++) {
                result_outer_product[index][col_ct_matrix] = output_vector_cussen_approach[index];
            }
          }
         if (col_index != 0) {
             add_ECP_Matrices(result2_cussen/* ECP_NIST256** A */, result_outer_product/* ECP_NIST256** B */, result2_cussen/* ECP_NIST256** C */, num_elements/* int size */);
         } else {
             copy_matrix(result2_cussen/* ECP_NIST256** A */, result_outer_product/* ECP_NIST256** B */, num_elements/* int size */); // A = B - copy all matrix elements
         }  
      }

// --------------------------------------------- Reconstruction phase for second ECP point matrix in the ciphertext end ------------------------------------

#ifndef ONLY_FORWARD_PATH
// Comment below 3 lines for "only forward path" time profiling of cussen algorithm. Uncomment them for "entire cussen algorithm" / "only backward path" time profiling.
#ifdef TIME_PROFILE
        clock_gettime (CLOCK_PROCESS_CPUTIME_ID, &time_end);
#endif
#endif

// ------------------------------------------------------------ Reconstruction phase end --------------------------------------------------------------------- 

// -------------------------------- Profiling timing for plaintext-ciphertext matrix multiplication using Cussen Approach end --------------------------------

// ------------------------------------------------------------ Post Processing begin ------------------------------------------------------------------------

#ifdef TIME_PROFILE
        time_elapsed = timespec_diff (time_start, time_end); // time_elapsed = (time_end.tv_nsec - time_start.tv_nsec);
//        printf("\nCussen Approach (ns): ");
//        printf("%ld\n", time_elapsed);
        time_TOTAL_cussen += time_elapsed;
#endif  
              
        int count_correct_output = 0;
        for (int i = 0; i < num_elements; i++) {
            for (int j = 0; j < num_elements; j++) {
                // Tests for equality of two ECPs - return 1 if P=Q, else returns 0 - extern int ECP_NIST256_equals(ECP_NIST256 *P, ECP_NIST256 *Q);
                if(ECP_NIST256_equals(&(result1_schoolbook[i][j])/* ECP_NIST256 *P */, &(result1_cussen[i][j])/* ECP_NIST256 *Q */) == 1) {  
                    if (ECP_NIST256_equals(&(result2_schoolbook[i][j])/* ECP_NIST256 *P */, &(result2_cussen[i][j])/* ECP_NIST256 *Q */) == 1) {
                        count_correct_output += 1; 
                    }                
                }
            }
        }

        if (count_correct_output == (num_elements * num_elements)) {
            printf("\nTest %d is successful for matrix size of %d\n", test_num, num_elements);
            num_successful_tests += 1;
        } else {
            printf("\nTest %d is not successful for matrix size of %d\n", test_num, num_elements);
        }
    
        // Free allocated memory
        freeMatrix(A, num_elements);
        freeMatrix(pointer_array_iter1, num_elements);
        freeMatrix(pointer_array_iter2, num_elements);
        freeMatrix(pointer_array_iter3, num_elements);
        freeMatrix(pointer_array_iter4, num_elements);
        free_ECP_Matrix(B, num_elements);
        free_ECP_Matrix(B1, num_elements);
        free_ECP_Matrix(result1_schoolbook, num_elements);
        free_ECP_Matrix(result2_schoolbook, num_elements);
        free_ECP_Matrix(result_outer_product, num_elements);
        free_ECP_Matrix(result1_cussen, num_elements);         
        free_ECP_Matrix(result2_cussen, num_elements);         
    }
   
    if (num_successful_tests == NUM_TESTS) {
        printf("\nAll tests successful for num_elements = %d!\n", num_elements);
    }
#ifdef TIME_PROFILE
    time_TOTAL_schoolbook /= NUM_TESTS;
#endif

#ifdef TIME_PROFILE
    time_TOTAL_cussen /= NUM_TESTS;
#if !defined(ONLY_FORWARD_PATH) && !defined(ONLY_BACKWARD_PATH)
    printf("Average Time (us) for schoolbook approach: %.3f\n", (time_TOTAL_schoolbook*1e-3));
    printf("Average Time (us) for cussen approach: %.3f\n", time_TOTAL_cussen*1e-3);
    printf("Timing in schoolbook approach / Timing in cussen approach: %.3f\n", (time_TOTAL_schoolbook*1e-3)/(time_TOTAL_cussen*1e-3));
#else
  #ifdef ONLY_FORWARD_PATH
    printf("Average Time (us) for only forward path of cussen approach: %.3f\n", time_TOTAL_cussen*1e-3);
  #endif
  #ifdef ONLY_BACKWARD_PATH
    printf("Average Time (us) for only backward path of cussen approach: %.3f\n", time_TOTAL_cussen*1e-3);
  #endif
#endif
#endif

#if !defined(ONLY_FORWARD_PATH) && !defined(ONLY_BACKWARD_PATH)
    // Populate each row corresponding to each vector length
    fprintf(csv_output_file, "%d, %d, %.3f, %.3f, %.3f\n", num_elements, BIT_WIDTH_EACH_PLAINTEXT, (time_TOTAL_schoolbook*1e-3), (time_TOTAL_cussen*1e-3), ((time_TOTAL_schoolbook*1e-3)/(time_TOTAL_cussen*1e-3)));
#else
    // Populate each row corresponding to each vector length
    fprintf(csv_output_file, "%d, %d, %.3f\n", num_elements, BIT_WIDTH_EACH_PLAINTEXT, (time_TOTAL_cussen*1e-3));
#endif
    // Close the file
    fclose(csv_output_file);

// ------------------------------------------------------------ Post Processing end ------------------------------------------------------------------------

  }

  printf("Data written to csv file\n");
  return 0;
}
