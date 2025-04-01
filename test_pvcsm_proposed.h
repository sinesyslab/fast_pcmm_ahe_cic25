// Corresponding output csv file will get populated on running this test!
	// bitsize of 04 -> sv_ct_size04.csv
	// bitsize of 08 -> sv_ct_size08.csv
	// bitsize of 12 -> sv_ct_size12.csv
	// bitsize of 16 -> sv_ct_size16.csv

#define BIT_WIDTH_EACH_PLAINTEXT 12 // Bit width of each plaintext vector element
#define NUM_TESTS 100     // Number of random trials to be averaged for profiling.
#define NUM_ITERATIONS 4  // Number of iterations in Cussen's Algorithm.
