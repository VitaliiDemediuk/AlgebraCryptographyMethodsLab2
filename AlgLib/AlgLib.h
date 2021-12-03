#ifndef ALGLIB_H
#define ALGLIB_H

#include "BigInt.h"

namespace alglib {

    al::BigInt gcd(const al::BigInt& x, const al::BigInt& y) noexcept;

    al::BigInt pollardRhoDivision(const al::BigInt& n, size_t nIt = 100000);
    std::vector<al::BigInt> pollardRhoFactorisation(const al::BigInt& n, size_t nIt = 100000);

    al::BigInt babyStepGiantStepLog(const al::BigInt& a, const al::BigInt& b, const al::BigInt& p) noexcept;

}

#endif //ALGLIB_H
