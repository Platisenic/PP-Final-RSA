#include <iostream>
#include <gmpxx.h>

// https://github.com/jumpmanmv/gmp

int main() {
  mpz_class x("7612058254738945");
  mpz_class y("9263591128439081");
  mpz_class xy = x * y;

  std::cout << "\n    " << x << "\n*\n    " << y;
  std::cout << "\n--------------------\n" << xy << "\n\n";
}

