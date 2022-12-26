#include <iostream>
#include <vector>
#include <stack>
#include <ctime>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <gmpxx.h>

#include "fermat.h"
#include "rho.h"
#include "qs.h"

// Input size threshold below which we resort to trial division.
const static uint32_t TRIAL_THRESHOLD = 1000000000;

// Maximum input size we can handle (bits).
const static uint32_t MAX_DIGITS = 300;

// GMP random number generator.
extern gmp_randclass rng;

int main(int argc, char *argv[]) {
    // Determine the algorithm to use.
    std::string argN;
    int threadNum = 1;
    bool badArg = false;
    if (argc == 2) {
        argN = std::string(argv[1]);
    }else if(argc ==3){
        argN = std::string(argv[1]);
        threadNum = atoi(argv[2]);
        if(threadNum <1 || threadNum >12){
            badArg = true;
        }
    }else{
        badArg = true;
    }

    if (badArg) {
        std::cerr << "Usage: quadratic {N} {threadNum}" << std::endl;
        std::cerr << "N stands for the Big Number of RSA" << std::endl;
        std::cerr << "threadNum must greater than 1 and less than 12" << std::endl;
        std::cerr << "threadNum will be set to 1 if argc == 2" << std::endl;

        return 1;
    }

    // Seed the GMP random number generator.
    rng.seed(time(0));

    // Seed the standard library random number generator.
    srand(time(0));

    // Find some primes for trial division.
    std::vector<uint32_t> primes;
    uint32_t max = ceil(sqrt(TRIAL_THRESHOLD)) + 1;
    std::vector<bool> sieve(max, false);
    for (uint32_t p = 2; p < max; ++p) {
        if (sieve[p])
            continue;
        primes.push_back(p);
        for (uint32_t i = p; i < max; i += p)
            sieve[i] = true;
    }

    // Factor each of the 100 integers.
    mpz_class N(argN);
    if (mpz_sizeinbase(N.get_mpz_t(), 2) > MAX_DIGITS) {
        std::cerr << "fail" << std::endl << std::endl; // Too many digits.
        return 1;
    }

    if (mpz_probab_prime_p(N.get_mpz_t(), 10)) {
        // N is prime.
        std::cout << N << std::endl << std::endl;
        return 1;
    }

    std::stack<mpz_class> factors;
    factors.push(N);

    while (!factors.empty()) {
        mpz_class factor = factors.top();
        factors.pop();

        // N is prime.
        if (mpz_probab_prime_p(factor.get_mpz_t(), 10)) {
            std::cout << factor << " ";
            continue;
        }
        mpz_class result = quadraticSieve(threadNum, factor);
        factors.push(result);
        factors.push(factor / result);
    }
    std::cout << std::endl;
    return 0;
}

