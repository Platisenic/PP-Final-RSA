#include <iostream>
#include <vector>
#include <cmath>
#include <gmp.h>
#include <gmpxx.h>

#include "CycleTimer.h"

int main() {
    mpz_class n("914758754609077333");
    mpz_class tmp = sqrt(n);
    long sqrtN = tmp.get_si();

    std::vector<bool> primeBool(sqrtN, true);
    primeBool[0] = primeBool[1] = false;

    double start = CycleTimer::currentSeconds();
    for (long i=2; i<sqrtN; i++) {
        if (primeBool[i]) {
            if (n%i == 0) {
                std::cout << i << " " << n/i << " " <<  n << std::endl;
                break;
            }
            for (long j=i+i; j<sqrtN; j+=i) {
                primeBool[j] = false;
            }
        }
    }
    double end = CycleTimer::currentSeconds();

    std::cout << "time: " << end-start << std::endl;

    return 0;
}
