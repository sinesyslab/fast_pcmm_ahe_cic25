# Fast Plaintext-Ciphertext Matrix Multiplication from Additively Homomorphic Encryption

This repository provides reference code for the software implementation presented in the following paper: K. S. T. Ramapragada and U. Banerjee, “Fast Plaintext-Ciphertext Matrix Multiplication from Additively Homomorphic Encryption,” IACR Communications in Cryptology, vol. 2, no. 1, pp. 1-33, April 2025.

### Prerequisites

The software implementation requires use of the [MIRACL core](https://github.com/miracl/core/tree/master/c) C library.

Following steps are necessary before running the PC-MM software tests:
1. Clone the MIRACL core repository: https://github.com/miracl/core.git
2. Setup the MIRACL core with the NIST P256 curve, install the library and do the basic tests as described in the [MIRACL documentation](https://github.com/miracl/core/blob/master/c/readme.md)
3. Copy the tests from this repository to the ```MIRACL/core/c``` directory

For the paper, the following tests were performed on Raspberry Pi 5 edge computing platform:
1. ```test_pvcsm_proposed.c```: Profile and compare plaintext vector and ciphertext scalar multiplication using Schoolbook and Proposed approaches (the Proposed approach is inspired by Cussen's algorithm and extends Cussen's algorithm from plaintext to encrypted setting)
2. ```test_pcmm_strassen.c```: Profile and compare plaintext-ciphertext matrix multiplication using Schoolbook and Strassen approaches
3. ```test_pcmm_proposed.c```: Profile and compare plaintext-ciphertext matrix multiplication using Schoolbook and Proposed approaches (the Proposed approach is inspired by Cussen's algorithm and extends Cussen's algorithm from plaintext to encrypted setting)

Note: Set the appropriate ```#define``` macros in the corresponding header file to run the test for various scenarios.

Example: The following command can be used to compile the tests for the proposed fast PC-MM: ```gcc -O3 test_pcmm_proposed.c core.a -o test_pcmm_proposed -lm```

