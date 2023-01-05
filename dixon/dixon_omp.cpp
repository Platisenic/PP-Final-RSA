#include <stdio.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <stdint.h>
#include <gmp.h>
#include <gmpxx.h>
#include <omp.h>
#include "CycleTimer.h"
#include "matrix.h"
#include <algorithm>
using namespace std;
// Minimal smoothness bound.
const static uint32_t MINIMAL_BOUND = 300;

// Sieving interval length.
const static uint32_t interval =65536;
std::vector<uint32_t> generateFactorBase(const mpz_class& n, uint32_t B) {
    std::vector<uint32_t> factorBase;
    /*
     * This is essentially the sieve of Eratosthenes, with the addition
     * that it also checks for (n/p) = 1.
     */
    std::vector<bool> sieve(B + 1, false);
    for (uint32_t p = 2; p <= B; ++p) {
        if (sieve[p])
            continue;

        // Add p to factor base if n is a quadratic residue modulo p.
        if (mpz_legendre(n.get_mpz_t(), mpz_class(p).get_mpz_t()) == 1)
            factorBase.push_back(p);

        // Add multiples of p to sieve.
        for (uint32_t i = p; i <= B; i += p)
            sieve[i] = true;
    }

    return factorBase;
}

int main(int argc,char** argv){
    
    string num ;
    int thread_count;
    bool badArg = false;
    if(argc == 2 ){
        num = argv[1];
        thread_count = omp_get_max_threads();
    }else if(argc ==3){
        num =  argv[1];
        thread_count = atoi(argv[2]);
        if(thread_count <1 || thread_count>12) badArg = true;
    }else{
        badArg = true;
    }
    if (badArg) {
        std::cerr << "Usage: ./dixon_omp {N} {threadNum}" << std::endl;
        std::cerr << "N stands for the Big Number of RSA" << std::endl;
        std::cerr << "threadNum must greater than 1 and less than 12" << std::endl;
        std::cerr << "threadNum will be set to 1 if argc == 2" << std::endl;

        return 1;
    }
    
    omp_set_num_threads(thread_count);
    mpz_class n (num);
    const float logN = mpz_sizeinbase(n.get_mpz_t(), 2) * std::log(2); // Approx.
    const float loglogN = std::log(logN);
    const mpz_class sqrtN = sqrt(n);

    // Smoothness bound B.
    const uint32_t B = MINIMAL_BOUND + std::ceil(std::exp(0.55*std::sqrt(logN * loglogN)));
    const std::vector<uint32_t> factorBase = generateFactorBase(n, B);

    vector<uint32_t> smooth;
    vector<vector<uint32_t>> smoothFactors;
    
    omp_lock_t lck;
    omp_init_lock(&lck);
    double start = CycleTimer::currentSeconds();
#pragma omp parallel
{

    uint32_t id = omp_get_thread_num();
    uint32_t z = id*interval;
    while(smooth.size()<factorBase.size()+20){
        uint32_t a,i ;
        
        for(i=0,a = z;i<interval;i++,a++){
            // test if b-smooth
            mpz_class z_square = ((a+sqrtN)*(a+sqrtN))%n; 
            vector<uint32_t> factors;

            double p_now = 0;
            for(uint32_t j=0;j<factorBase.size();++j){ 
                const int p = factorBase[j];
                
                while(mpz_divisible_ui_p(z_square.get_mpz_t(), p)){
                    mpz_divexact_ui(z_square.get_mpz_t(), z_square.get_mpz_t(), p);
                    factors.push_back(j);
                } 
            }
            //check whether z_sqaure is b_smooth
            if(z_square ==1){
                omp_set_lock(&lck);
                smooth.push_back(a); //need change
                smoothFactors.push_back(factors);
                omp_unset_lock(&lck);
            }
            
        }
        if(smooth.size()>factorBase.size()+20){ 
                break;
            }

        z+= thread_count*interval;
    }
}   

    Matrix M(factorBase.size(),smoothFactors.size()+1);

    for (uint32_t i = 0; i < smoothFactors.size(); ++i) {
        for (uint32_t j = 0; j < smoothFactors[i].size(); ++j) {
            M(smoothFactors[i][j], i).flip();
        }
    }
    M.reduce();


    mpz_class a;
    mpz_class b;

    do {
        std::vector<uint32_t> x = M.solve();

        a = 1;
        b = 1;

        /*
         * Calculate b = product(smooth[i] + sqrt(n)).
         *
         * Also calculate the the power of each prime in a's decomposition on the
         * factor base.
         */
        std::vector<uint32_t> decomp(factorBase.size(), 0);
        for (uint32_t i = 0; i < smoothFactors.size(); ++i) {
            if (x[i] == 1) {
                for(uint32_t p = 0; p < smoothFactors[i].size(); ++p)
                    ++decomp[smoothFactors[i][p]];
                b *= (smooth[i] + sqrtN);
            }
        }

        /*
         * Calculate a = sqrt(product(factorBase[p])).
         */
        for(uint32_t p = 0; p < factorBase.size(); ++p) {
            for(uint32_t i = 0; i < (decomp[p] / 2); ++i)
                a *= factorBase[p];
        }

        // a = +/- b (mod n) means we have a trivial factor :(
    } while (a % n == b % n || a % n == (- b) % n + n);



    /************************************
     *                                  *
     * STAGE 4: Success!                *
     *                                  *
     ***********************************/
    mpz_class factor;
    mpz_gcd(factor.get_mpz_t(), mpz_class(b - a).get_mpz_t(), n.get_mpz_t());
    cout<<factor<<" "<<n/factor<<" "<<factor*(n/factor)<<endl;
    double conduct_time = CycleTimer::currentSeconds()-start;
    printf("time : %f \n",conduct_time);

}

