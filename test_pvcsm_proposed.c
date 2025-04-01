// This test is to profile plaintext vector and ciphertext scalar multiplication using Schoolbook and Proposed approaches (the Proposed approach is inspired by Cussen's algorithm and extends Cussen's algorithm from plaintext to encrypted setting)

// -------------------------------------------------------------------- Test Begins -----------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h> // To measure time and do profiling. Required for also random number generation.
#include <math.h> // To use "round" function for rounding the vector lengths after each iteration to nearest integer.
#include "randapi.h" // For CREATE_CSPRNG function call

#include "test_pvcsm_proposed.h"

// Include NIST256 curve specific header files
#include "ecp_NIST256.h"
#include "fp_NIST256.h"

// To profile for "both compression and reconstruction phase", comment both the lines below. To profile for "only compression phase" / "only reconstruction phase" scenarios, uncomment only the corresponding define accordingly.
//#define ONLY_COMPRESSION_PHASE // Define to time profile "only compression phase of cussen's algorithm"
//#define ONLY_RECONSTRUCTION_PHASE // Define to time profile "only reconstruction phase of cussen's algorithm"

#define TIME_PROFILE

// --------------------------------------------------------- Helper Functions Begin --------------------------------------------------------------------------

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

 for (int num_elements = 8; num_elements <= 512; (num_elements *= 2)) {
    FILE *csv_output_file;
    if (num_elements == 8) {
        if (BIT_WIDTH_EACH_PLAINTEXT == 4) {
            csv_output_file = fopen("sv_ct_size04.csv", "w");
        } else if (BIT_WIDTH_EACH_PLAINTEXT == 8) {
            csv_output_file = fopen("sv_ct_size08.csv", "w");
        } else if (BIT_WIDTH_EACH_PLAINTEXT == 12) {
            csv_output_file = fopen("sv_ct_size12.csv", "w");
        } else if (BIT_WIDTH_EACH_PLAINTEXT == 16) {
            csv_output_file = fopen("sv_ct_size16.csv", "w");
        }
    } else {
        if (BIT_WIDTH_EACH_PLAINTEXT == 4) {
            csv_output_file = fopen("sv_ct_size04.csv", "a");
        } else if (BIT_WIDTH_EACH_PLAINTEXT == 8) {
            csv_output_file = fopen("sv_ct_size08.csv", "a");
        } else if (BIT_WIDTH_EACH_PLAINTEXT == 12) {
            csv_output_file = fopen("sv_ct_size12.csv", "a");
        } else if (BIT_WIDTH_EACH_PLAINTEXT == 16) {
            csv_output_file = fopen("sv_ct_size16.csv", "a");
        }
    }

    if (csv_output_file == NULL) {
        printf("Error opening file!\n");
        return 1;
    }

#if !defined(ONLY_COMPRESSION_PHASE) && !defined(ONLY_RECONSTRUCTION_PHASE)
  // Since all the matrix size measurements corresponding to a particular bit size will be in one (i.e., same) excel file, header need to be written for only the first matrix size profiling run corresponding to the particular bit size
    if (num_elements == 8) {
        // Write header
        fprintf(csv_output_file, "Vector Length, L_1, L_2, L_3, L_4, Schoolbook (us), Cussen (us), Schoolbook / Cussen\n");
    }
#else
 // Write header
    if (num_elements == 8) {
        // Write header
        fprintf(csv_output_file, "Vector Length, Cussen (us)\n");
    }
#endif
 
    ECP_NIST256 point_generator, random_point1, random_point2;

#if CURVETYPE_NIST256==MONTGOMERY
    // Set ECP to point(x,[y]) given x - extern int ECP_NIST256_set(ECP_NIST256 *P, BIG_256_56 x);
#else  
    // Set ECP to point(x,y) given x and y - ECP_NIST256_set(ECP_NIST256 *P, BIG_256_56 x, BIG_256_56 y)
    ECP_NIST256_set(&point_generator, CURVE_Gx_NIST256, CURVE_Gy_NIST256);
#endif
   
    random_point1 = point_generator;
    random_point2 = point_generator;

    unsigned int test_num;
    int num_successful_tests = 0;
#ifdef TIME_PROFILE
    struct timespec time_start, time_end;
    long time_elapsed;
    long time_TOTAL_schoolbook = 0;
    long time_TOTAL_cussen = 0;
    float length_of_vector_after_iter1 = 0;
    float length_of_vector_after_iter2 = 0;
    float length_of_vector_after_iter3 = 0;
    float length_of_vector_after_iter4 = 0;
#endif
     
    for (test_num = 0; test_num < NUM_TESTS; test_num++)
    {
        // Create an unbiased random BIG from a random number generator, reduced with respect to a modulus (Constant Time as used) - extern void BIG_256_56_randomnum(BIG_256_56 x, BIG_256_56 n, csprng *r);
        BIG_256_56 random_num;

// -----------------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------- Two random ECP Points for ciphertext begin ------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------------------------

        BIG_256_56_randomnum(random_num/* BIG_256_56 x */, modulus_for_random_num_gen/* BIG_256_56 n */, &RNG/* csprng *r */);
        // Multiplies an ECP instance P by a BIG, side-channel resistant - extern void ECP_NIST256_clmul(ECP_NIST256 *P, BIG_256_56 e, BIG_256_56 maxe);
        ECP_NIST256_clmul(&random_point1 /* ECP_NIST256 *P */, random_num /* BIG_256_56 e */, max_256_bit_number/* BIG_256_56 maxe */);
        
        BIG_256_56_randomnum(random_num/* BIG_256_56 x */, modulus_for_random_num_gen/* BIG_256_56 n */, &RNG/* csprng *r */);
        // Multiplies an ECP instance P by a BIG, side-channel resistant - extern void ECP_NIST256_clmul(ECP_NIST256 *P, BIG_256_56 e, BIG_256_56 maxe);
        ECP_NIST256_clmul(&random_point2 /* ECP_NIST256 *P */, random_num /* BIG_256_56 e */, max_256_bit_number/* BIG_256_56 maxe */);

// --------------------------------------------------- Two random ECP Points for ciphertext end ------------------------------------------------------------

        srand(time(0) + test_num); // Seed the PRNG
 
//      printf("\nTest: %d\n", test_num);
       
        int input_vector[num_elements];
        ECP_NIST256 output_vector1_schoolbook_approach[num_elements], output_vector2_schoolbook_approach[num_elements];

// -----------------------------------------------------------------------------------------------------------------------------------------------------------
// -------------------------------------------- Random plaintext vector elements generation begin ------------------------------------------------------------
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
        // Populating input vector of 'BIT_WIDTH_EACH_PLAINTEXT'-bit unsigned numbers and calculating output vector in schoolbook approach.
        for (int i = 0; i < num_elements; i++) {
            input_vector[i] = rand() & mask;        
        }

// -------------------------------------------- Random plaintext vector elements generation end ------------------------------------------------------------

// -----------------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------- Output in Schoolbook method begin ---------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------------------------

#if !defined(ONLY_COMPRESSION_PHASE) && !defined(ONLY_RECONSTRUCTION_PHASE)
// Profiling Schoolbook Approach timing
#ifdef TIME_PROFILE
//        printf("Schoolbook Approach (ns): ");
        clock_gettime (CLOCK_PROCESS_CPUTIME_ID, &time_start);
#endif
        for (int i = 0; i < num_elements; i++) {
            // Calculating output vector in schoolbook approach.
            // Considering the cipher text in vector-ciphertext multiplication as the earlier generated "random_point1".
            output_vector1_schoolbook_approach[i] = random_point1;
       
            // Multiplies an ECP instance P by a small integer, side-channel resistant - extern void ECP_NIST256_pinmul(ECP_NIST256 *P, int i, int b);
            ECP_NIST256_pinmul(&output_vector1_schoolbook_approach[i]/* ECP_NIST256 *P */, input_vector[i]/* int i */, (BIT_WIDTH_EACH_PLAINTEXT)/* int b */);
            
            // Considering the cipher text in vector-ciphertext multiplication as the earlier generated "random_point2".
            output_vector2_schoolbook_approach[i] = random_point2;
       
            // Multiplies an ECP instance P by a small integer, side-channel resistant - extern void ECP_NIST256_pinmul(ECP_NIST256 *P, int i, int b);
            ECP_NIST256_pinmul(&output_vector2_schoolbook_approach[i]/* ECP_NIST256 *P */, input_vector[i]/* int i */, (BIT_WIDTH_EACH_PLAINTEXT)/* int b */);
        }
#ifdef TIME_PROFILE
        clock_gettime (CLOCK_PROCESS_CPUTIME_ID, &time_end);
        time_elapsed = timespec_diff (time_start, time_end); // time_elapsed = (time_end.tv_nsec - time_start.tv_nsec);
   
//      printf("%ld\n", time_elapsed);
        time_TOTAL_schoolbook += time_elapsed;
#endif
#endif

// --------------------------------------------------- Output in Schoolbook method end ---------------------------------------------------------------------

// -----------------------------------------------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------------- Output in Cussen approach begin ---------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------------------------

#ifndef ONLY_RECONSTRUCTION_PHASE
#ifdef TIME_PROFILE
        clock_gettime (CLOCK_PROCESS_CPUTIME_ID, &time_start);
#endif
#endif
        int input_vector_after_diff[num_elements];
        int input_vector_after_diff_copy[num_elements];
        int vector_after_iter[num_elements];
        int pointer_array_iter1[num_elements];
        int pointer_array_iter2[num_elements];
        int pointer_array_iter3[num_elements];
        int pointer_array_iter4[num_elements];
       
        int length_of_vector_before_elim_dupl_iter1, length_of_vector_before_elim_dupl_iter2, length_of_vector_before_elim_dupl_iter3, length_of_vector_before_elim_dupl_iter4;
        int length_of_vector_before_elim_dupl;

        ECP_NIST256 output_vector_cussen_approach_post_iter1[num_elements];
        ECP_NIST256 output_vector_cussen_approach_post_iter2[num_elements];
        ECP_NIST256 output_vector_cussen_approach_post_iter3[num_elements];
        ECP_NIST256 output_vector_cussen_approach_post_iter4[num_elements];

        ECP_NIST256 output_vector1_cussen_approach[num_elements], output_vector2_cussen_approach[num_elements];
 
        int length_of_vector_after_iter = num_elements; // Initializing to 'num_elements' for the first iteration

       
// -----------------------------------------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------- Compression Phase begin -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------------------------

// Cussen's algorithm (No alignment scenario and four iterations)
// Compression phase - Taking differences, sorting, eliminating duplicates and populating tracking pointers array  
        for (int iteration = 0; iteration < NUM_ITERATIONS; iteration ++) {
        // Start of Iteration (Compression phase)            
            if (iteration != 0) {
                input_vector_after_diff[0] = input_vector[0];
                for (int index = 1; index < length_of_vector_after_iter; index ++) {
                    input_vector_after_diff[index] = input_vector[index] - input_vector[index - 1];
                }
            } else {
                for (int index = 0; index < length_of_vector_after_iter; index ++) {
                    input_vector_after_diff[index] = input_vector[index];
                }
            }
   
            for(int index = 0; index < length_of_vector_after_iter; index ++) {
                input_vector_after_diff_copy[index] = input_vector_after_diff[index];
            }

            // Quick sort function in C - void qsort(void* base, size_t num, size_t size, int (*compare)(const void*, const void*));
            qsort(input_vector_after_diff, length_of_vector_after_iter, sizeof(unsigned int), compare);
   
            if (iteration == 0) {
                length_of_vector_before_elim_dupl_iter1 = length_of_vector_after_iter;
                length_of_vector_before_elim_dupl = length_of_vector_before_elim_dupl_iter1;
            } else if (iteration == 1) {
                length_of_vector_before_elim_dupl_iter2 = length_of_vector_after_iter;
                length_of_vector_before_elim_dupl = length_of_vector_before_elim_dupl_iter2;
            } else if (iteration == 2) {
                length_of_vector_before_elim_dupl_iter3 = length_of_vector_after_iter;
                length_of_vector_before_elim_dupl = length_of_vector_before_elim_dupl_iter3;
            } else if (iteration == 3) {
                length_of_vector_before_elim_dupl_iter4 = length_of_vector_after_iter;
                length_of_vector_before_elim_dupl = length_of_vector_before_elim_dupl_iter4;
            }
   
            // Eliminating Duplicates
            int index_vector_after_iter = 0;
            for (int index = 0; index < length_of_vector_after_iter; index++) {
                if (index == 0) {
                    vector_after_iter[0] = input_vector_after_diff[0];
                    index_vector_after_iter = 1;
                } else if (vector_after_iter[index_vector_after_iter - 1] != input_vector_after_diff[index]) {
                    vector_after_iter[index_vector_after_iter ++] = input_vector_after_diff[index];
                }
            }
           
            length_of_vector_after_iter = index_vector_after_iter;
            
            // Vector of pointers for tracking
            for (int i = 0; i < length_of_vector_before_elim_dupl; i++) {
                for (int j = 0; j < length_of_vector_after_iter; j++) {
            if (input_vector_after_diff_copy[i] == vector_after_iter[j]) {
                        if (iteration == 0) {
                            pointer_array_iter1[i] = j;
                        } else if (iteration == 1) {
                            pointer_array_iter2[i] = j;
                        } else if (iteration == 2) {
                            pointer_array_iter3[i] = j;
                        } else if (iteration == 3) {
                            pointer_array_iter4[i] = j;
                        }
                        break;
                    }
                }
            }

            if (test_num != 0) {
                if (iteration == 0) {
                    length_of_vector_after_iter1 += length_of_vector_after_iter;
                } else if (iteration == 1) {
                    length_of_vector_after_iter2 += length_of_vector_after_iter;
                } else if (iteration == 2) {
                    length_of_vector_after_iter3 += length_of_vector_after_iter;
                } else if (iteration == 3) {
                    length_of_vector_after_iter4 += length_of_vector_after_iter;
                }
            }

            // Update the input vector for the next iteration
            for (int index = 0; index < length_of_vector_after_iter; index ++) {
                input_vector[index] = vector_after_iter[index];
            }
        // End of Iteration (Compression phase)
        }

#ifdef ONLY_COMPRESSION_PHASE
#ifdef TIME_PROFILE
        clock_gettime (CLOCK_PROCESS_CPUTIME_ID, &time_end);
#endif
#endif

// ---------------------------------------------------- Compression Phase end -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------------------------------------------------------------------------------------
// ------------------------------------- Reconstruction phase for first ECP point in the ciphertext begin ----------------------------------------------------  
// -----------------------------------------------------------------------------------------------------------------------------------------------------------
// Reconstruction phase - Point Multiplication, Copying duplicates using pointers into appropriate places, Point Additions.        

#ifdef ONLY_RECONSTRUCTION_PHASE
#ifdef TIME_PROFILE
        clock_gettime (CLOCK_PROCESS_CPUTIME_ID, &time_start);
#endif
#endif

        for (int index = 0; index < length_of_vector_after_iter; index ++) {
            output_vector1_cussen_approach[index] = random_point1;
            // Multiplies an ECP instance P by a small integer, side-channel resistant - extern void ECP_NIST256_pinmul(ECP_NIST256 *P, int i, int b);
            ECP_NIST256_pinmul(&output_vector1_cussen_approach[index]/* ECP_NIST256 *P */, input_vector[index]/* int i */, (BIT_WIDTH_EACH_PLAINTEXT)/* int b */);
        }

        for (int iteration = (NUM_ITERATIONS - 1); iteration >= 0; iteration --) {
        // Start of Iteration (Backward Path)
            if (iteration == 0) {
                length_of_vector_before_elim_dupl = length_of_vector_before_elim_dupl_iter1;
            } else if (iteration == 1) {
                length_of_vector_before_elim_dupl = length_of_vector_before_elim_dupl_iter2;
            } else if (iteration == 2) {
                length_of_vector_before_elim_dupl = length_of_vector_before_elim_dupl_iter3;
            } else if (iteration == 3) {
                length_of_vector_before_elim_dupl = length_of_vector_before_elim_dupl_iter4;
            }
           
            for (int index = 0; index < length_of_vector_before_elim_dupl; index ++) {
                if (iteration == 0) {
                    output_vector_cussen_approach_post_iter1[index] = output_vector1_cussen_approach[pointer_array_iter1[index]];
                } else if (iteration == 1) {
                    output_vector_cussen_approach_post_iter2[index] = output_vector1_cussen_approach[pointer_array_iter2[index]];
                } else if (iteration == 2) {
                    output_vector_cussen_approach_post_iter3[index] = output_vector1_cussen_approach[pointer_array_iter3[index]];
                } else if (iteration == 3) {
                    output_vector_cussen_approach_post_iter4[index] = output_vector1_cussen_approach[pointer_array_iter4[index]];
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
                   output_vector1_cussen_approach[index] = output_vector_cussen_approach_post_iter4[index];
               } else if (iteration == 2) {
                   output_vector1_cussen_approach[index] = output_vector_cussen_approach_post_iter3[index];
               } else if (iteration == 1) {
                   output_vector1_cussen_approach[index] = output_vector_cussen_approach_post_iter2[index];
               } else if (iteration == 0) {
                   output_vector1_cussen_approach[index] = output_vector_cussen_approach_post_iter1[index];
               }
            }
        // End of iteration (Backward Path)
        } 

// ------------------------------------- Reconstruction phase for first ECP point in the ciphertext end ----------------------------------------------------  

// -----------------------------------------------------------------------------------------------------------------------------------------------------------
// ------------------------------------------- Reconstruction phase for second ECP point in the ciphertext begin ---------------------------------------------  
// -----------------------------------------------------------------------------------------------------------------------------------------------------------
// Reconstruction phase - Point Multiplication, Copying duplicates using pointers into appropriate places, Point Additions.        

        for (int index = 0; index < length_of_vector_after_iter; index ++) {
            output_vector2_cussen_approach[index] = random_point2;
            // Multiplies an ECP instance P by a small integer, side-channel resistant - extern void ECP_NIST256_pinmul(ECP_NIST256 *P, int i, int b);
            ECP_NIST256_pinmul(&output_vector2_cussen_approach[index]/* ECP_NIST256 *P */, input_vector[index]/* int i */, (BIT_WIDTH_EACH_PLAINTEXT)/* int b */);
        }

        for (int iteration = (NUM_ITERATIONS - 1); iteration >= 0; iteration --) {
        // Start of Iteration (Backward Path)
            if (iteration == 0) {
                length_of_vector_before_elim_dupl = length_of_vector_before_elim_dupl_iter1;
            } else if (iteration == 1) {
                length_of_vector_before_elim_dupl = length_of_vector_before_elim_dupl_iter2;
            } else if (iteration == 2) {
                length_of_vector_before_elim_dupl = length_of_vector_before_elim_dupl_iter3;
            } else if (iteration == 3) {
                length_of_vector_before_elim_dupl = length_of_vector_before_elim_dupl_iter4;
            }
           
            for (int index = 0; index < length_of_vector_before_elim_dupl; index ++) {
                if (iteration == 0) {
                    output_vector_cussen_approach_post_iter1[index] = output_vector2_cussen_approach[pointer_array_iter1[index]];
                } else if (iteration == 1) {
                    output_vector_cussen_approach_post_iter2[index] = output_vector2_cussen_approach[pointer_array_iter2[index]];
                } else if (iteration == 2) {
                    output_vector_cussen_approach_post_iter3[index] = output_vector2_cussen_approach[pointer_array_iter3[index]];
                } else if (iteration == 3) {
                    output_vector_cussen_approach_post_iter4[index] = output_vector2_cussen_approach[pointer_array_iter4[index]];
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
                   output_vector2_cussen_approach[index] = output_vector_cussen_approach_post_iter4[index];
               } else if (iteration == 2) {
                   output_vector2_cussen_approach[index] = output_vector_cussen_approach_post_iter3[index];
               } else if (iteration == 1) {
                   output_vector2_cussen_approach[index] = output_vector_cussen_approach_post_iter2[index];
               } else if (iteration == 0) {
                   output_vector2_cussen_approach[index] = output_vector_cussen_approach_post_iter1[index];
               }
            }
        // End of iteration (Backward Path)
        }
#ifndef ONLY_COMPRESSION_PHASE
#ifdef TIME_PROFILE
        clock_gettime (CLOCK_PROCESS_CPUTIME_ID, &time_end);
#endif
#endif

// ------------------------------------------- Reconstruction phase for second ECP point in the ciphertext end ---------------------------------------------  

// ----------------------------------------------------- Output in Cussen approach end ---------------------------------------------------------------------

// ----------------------------------------------------- Post Processing Start -----------------------------------------------------------------------------

#ifdef TIME_PROFILE
        time_elapsed = timespec_diff (time_start, time_end); // time_elapsed = (time_end.tv_nsec - time_start.tv_nsec);
       
//        printf("\nCussen Approach (ns): ");
//        printf("%ld\n", time_elapsed);
        time_TOTAL_cussen += time_elapsed;
#endif  
        int count_correct_output = 0;
        for (int index = 0; index < num_elements; index++) {
            // Tests for equality of two ECPs - return 1 if P=Q, else returns 0 - extern int ECP_NIST256_equals(ECP_NIST256 *P, ECP_NIST256 *Q);
            if(ECP_NIST256_equals(&output_vector1_cussen_approach[index]/* ECP_NIST256 *P */, &output_vector1_schoolbook_approach[index]/* ECP_NIST256 *Q */) == 1) {
                if(ECP_NIST256_equals(&output_vector1_cussen_approach[index]/* ECP_NIST256 *P */, &output_vector1_schoolbook_approach[index]/* ECP_NIST256 *Q */) == 1) {
                    count_correct_output += 1;
                }
            }
        }
        if (count_correct_output == num_elements) {
            num_successful_tests += 1;
        }
    }
   
    if (num_successful_tests == NUM_TESTS) {
        printf("\nAll tests successful for num_elements = %d!\n", num_elements);
    }
#ifdef TIME_PROFILE
    time_TOTAL_schoolbook /= NUM_TESTS;
#endif

#ifdef TIME_PROFILE
    time_TOTAL_cussen /= NUM_TESTS;
#if !defined(ONLY_COMPRESSION_PHASE) && !defined(ONLY_RECONSTRUCTION_PHASE)
    length_of_vector_after_iter1 /= NUM_TESTS;
    length_of_vector_after_iter2 /= NUM_TESTS;
    length_of_vector_after_iter3 /= NUM_TESTS;
    length_of_vector_after_iter4 /= NUM_TESTS;
    // round() returns the nearest integer (returns in double ="%f" format) to the given value. It rounds halfway cases away from zero.
    printf("Average length of vector after iteration 1 is %d\n", (int) round(length_of_vector_after_iter1));
    printf("Average length of vector after iteration 2 is %d\n", (int) round(length_of_vector_after_iter2));
    printf("Average length of vector after iteration 3 is %d\n", (int) round(length_of_vector_after_iter3));
    printf("Average length of vector after iteration 4 is %d\n", (int) round(length_of_vector_after_iter4));
    printf("Average Time (us) for schoolbook approach: %.3f\n", (time_TOTAL_schoolbook*1e-3));
    printf("Average Time (us) for cussen approach: %.3f\n", time_TOTAL_cussen*1e-3);
    printf("Timing in schoolbook approach / Timing in cussen approach: %.3f\n", (time_TOTAL_schoolbook*1e-3)/(time_TOTAL_cussen*1e-3));
#else
  #ifdef ONLY_COMPRESSION_PHASE
    printf("Average Time (us) for only compression phase of cussen approach: %.3f\n", time_TOTAL_cussen*1e-3);
  #endif
  #ifdef ONLY_RECONSTRUCTION_PHASE
    printf("Average Time (us) for only reconstruction phase of cussen approach: %.3f\n", time_TOTAL_cussen*1e-3);
  #endif
#endif
#endif

#if !defined(ONLY_COMPRESSION_PHASE) && !defined(ONLY_RECONSTRUCTION_PHASE)
    // Populate each row corresponding to each vector length
    fprintf(csv_output_file, "%d, %d, %d, %d, %d, %.3f, %.3f, %.3f\n", num_elements, (int) round(length_of_vector_after_iter1), (int) round(length_of_vector_after_iter2), (int) round(length_of_vector_after_iter3), (int) round(length_of_vector_after_iter4), (time_TOTAL_schoolbook*1e-3), (time_TOTAL_cussen*1e-3), ((time_TOTAL_schoolbook*1e-3)/(time_TOTAL_cussen*1e-3)));
#else
    // Populate each row corresponding to each vector length
    fprintf(csv_output_file, "%d, %.3f\n", num_elements, (time_TOTAL_cussen*1e-3));
#endif
    // Close the file
    fclose(csv_output_file);

// ----------------------------------------------------- Post Processing End -----------------------------------------------------------------------------
  }

    printf("Data written to csv file\n");
    return 0;
}

