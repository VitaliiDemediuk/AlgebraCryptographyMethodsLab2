/// My header
#include "AlgLib.h"

/// STL
#include <map>
#include <random>

using namespace alglib;

al::BigInt alglib::gcd(const al::BigInt& x, const al::BigInt& y) noexcept {
    return x % y == 0 ? y : alglib::gcd(y, x%y);
}

al::BigInt alglib::pollardRhoDivision(const al::BigInt &n, size_t nIt) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> distrib(1, (n % al::BigInt(1000000)).toSizeT());

    al::BigInt b0 = al::BigInt(distrib(gen)) % n, b1 = b0, g;
    b1 = al::multiplyMod(b1, b1, n);
    if (++b1 == n)
        b1 = 0;
    g = gcd (abs (b1 - b0), n);
    for (size_t i=0; i < nIt && (g == 1 || g == n); ++i) {
        b0 = al::multiplyMod(b0, b0, n);
        if (++b0 == n) {
            b0 = 0;
        }
        b1 = al::multiplyMod(b1, b1, n);
        ++b1;
        b1 = al::multiplyMod(b1, b1, n);
        if (++b1 == n) {
            b1 = 0;
        }
        g = gcd (abs (b1 - b0), n);
    }
    return g;
}

auto alglib::pollardRhoFactorisation(const al::BigInt& n, size_t nIt) -> std::map<al::BigInt, al::BigInt> {
    auto nCopy = n;
    std::map<al::BigInt, al::BigInt> res;
    al::BigInt p = pollardRhoDivision(nCopy, nIt);
    while (p != 1) {
        ++res[p];
        nCopy /= p;
        p = pollardRhoDivision(nCopy, nIt);
    }
    if (nCopy != 1) {
        ++res[nCopy];
    }
    return res;
}

al::BigInt alglib::babyStepGiantStepLog(const al::BigInt& a,
                const al::BigInt& b, const al::BigInt& p) noexcept
{
    const al::BigInt n = al::sqrt(p) + 1;
    const al::BigInt an = al::powMod(a, n, p);

    std::map<al::BigInt, al::BigInt> np;
    for (al::BigInt i = 1, tmp = an; i <= n; ++i) {
        if (!np.count(tmp)) {
            np[tmp] = i;
        }
        tmp = al::multiplyMod(tmp, an, p);
    }

    for (al::BigInt i = 0, tmp = b; i <= n; ++i) {
        if (np.count(tmp)) {
            al::BigInt ans = np[tmp] * n - i;
            if (ans < p) {
                return ans;
            }
        }
        tmp = al::multiplyMod(tmp, a, p);
    }

    return -1;
}

al::BigInt alglib::euler(const al::BigInt& n) {
    al::BigInt res{1};
    for (const auto& [p, count] : pollardRhoFactorisation(n)) {
        res *= al::pow(p-1, count);
    }
    return res;
}

al::BigInt alglib::mobius(const al::BigInt& n) {
    if (n == 1) { return 1; }
    auto fact = alglib::pollardRhoFactorisation(n);
    for (const auto& [p, count] : fact) {
        if (count > 1) {
            return 0;
        }
    }
    return fact.size() % 2 == 0 ? 1 : -1;
}

al::BigInt alglib::legendre(const al::BigInt& a, const al::BigInt& p) {
    return al::powMod(a, (p-1)/2, p);
}

al::BigInt alglib::jacobi(const al::BigInt& a, const al::BigInt& p) {
    al::BigInt res{1};
    for (const auto& [pi, count] : pollardRhoFactorisation(p)) {
        res *= al::pow(alglib::legendre(a, pi), count);
    }
    return res;
}