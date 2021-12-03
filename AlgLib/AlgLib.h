#ifndef ALGLIB_H
#define ALGLIB_H

#include "BigInt.h"
#include <map>

namespace alglib {

    al::BigInt gcd(const al::BigInt& x, const al::BigInt& y) noexcept;

    al::BigInt pollardRhoDivision(const al::BigInt& n, size_t nIt = 100000);
    auto pollardRhoFactorisation(const al::BigInt& n, size_t nIt = 100000) -> std::map<al::BigInt, al::BigInt>;

    al::BigInt babyStepGiantStepLog(const al::BigInt& a, const al::BigInt& b, const al::BigInt& p) noexcept;

    al::BigInt euler(const al::BigInt& n);
    al::BigInt mobius(const al::BigInt& n);

    /// undefined behavior if p is not prime number
    al::BigInt legendre(const al::BigInt& a, const al::BigInt& p);
    al::BigInt jacobi(const al::BigInt& a, const al::BigInt& p);

    auto cipolla(al::BigInt a, al::BigInt p, int countRetry = 10) -> std::pair<al::BigInt, al::BigInt>;

    bool testPrimeMillerRabin(const al::BigInt& a, size_t nIt = 10);

}

#endif //ALGLIB_H
