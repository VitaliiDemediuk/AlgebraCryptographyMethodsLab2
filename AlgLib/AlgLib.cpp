/// My header
#include "AlgLib.h"

/// STL
#include <map>
#include <random>

/// Boost
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/random.hpp>

using namespace alglib;

namespace {

    al::BigInt randBigInt(const al::BigInt& f, const al::BigInt& l) noexcept {
        using namespace boost::multiprecision;

        std::random_device rd;
        std::mt19937 gen(rd());
        boost::random::uniform_int_distribution<cpp_int> uid(cpp_int(f.toStdString()), cpp_int(l.toStdString()));
        return uid(gen).convert_to<std::string>();
    }

    al::BigInt randBigInt(const al::BigInt& n)  {
        using namespace boost::multiprecision;

        std::random_device rd;
        std::mt19937 gen(rd());
        boost::random::uniform_int_distribution<cpp_int> uid(cpp_int(2), cpp_int(n.toStdString()));
        return uid(gen).convert_to<std::string>();
    }

} // anon namespace

al::BigInt alglib::gcd(const al::BigInt& x, const al::BigInt& y) noexcept {
    return x % y == 0 ? y : alglib::gcd(y, x%y);
}

al::BigInt alglib::pollardRhoDivision(const al::BigInt &n, size_t nIt) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> distrib(1, (n % al::BigInt(1000000)).toSizeT());

    al::BigInt b0 = al::BigInt(distrib(gen)) % n;
    al::BigInt b1 = b0;
    al::BigInt g;
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

    while (nCopy % 2 == 0) {
        nCopy = nCopy / 2;
        if (res.empty()) {
            res[2] = 1;
        } else {
            ++res[2];
        }
    }

    al::BigInt p = pollardRhoDivision(nCopy, nIt);
    while (p != 1) {
        nCopy /= p;
        al::BigInt count = 1;
        while (nCopy % p == 0) {
            ++count;
            nCopy /= p;
        }
        res[p] = count;
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

bool alglib::testPrimeMillerRabin(const al::BigInt& a, size_t nIt) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> distrib(1, (a % al::BigInt(1000000)).toSizeT());

    if (a == 2 || a == 3) {
        return true;
    }
    if (a == 1 || a % 2 == 0 || a % 3 == 0) {
        return false;
    }

    al::BigInt q = a - 1;
    al::BigInt t = q;
    int n = 0;
    while (q % 2 == 0) {
        q = q / 2;
        n++;
    }
    for (long long i = 0; i < nIt; i++) {
        al::BigInt k = distrib(gen);
        al::BigInt b = al::powMod(k, q, a);
        if (b == 1 || b == t) {
            continue;
        }
        al::BigInt j = 0;
        while (j < n) {
            auto two = al::BigInt(2);
            b = al::powMod(b, two, a);
            if (b == t) {
                break;
            }
            j++;
        }
        if (j < n) {
            continue;
        }
        else {
            return false;
        }
    }
    return true;
}

namespace {

    class SqrtPolynomial {
    public:
        SqrtPolynomial(al::BigInt x, al::BigInt y, al::BigInt w)
                : x_{std::move(x)}, y_{std::move(y)}, w_{std::move(w)} {}

        SqrtPolynomial pow(const al::BigInt& n, const al::BigInt& mod) {
            if (n == 0) {
                return {1, 0, w_};
            }
            if (n == 1) {
                return *this;
            }
            SqrtPolynomial tmp = this->pow(n / 2, mod);
            SqrtPolynomial tmp1((tmp.x_ * tmp.x_ % mod + tmp.y_ * tmp.y_ * w_ % mod) % mod, tmp.x_ * 2 * tmp.y_ % mod, w_);
            if (n % 2 == 0) {
                return tmp1;
            } else {
                SqrtPolynomial tmp2((tmp1.x_ * x_ % mod + tmp1.y_ * y_ * w_ % mod) % mod,
                                     (tmp1.x_ * y_ % mod + tmp1.y_ * x_ % mod) % mod, w_);
                return tmp2;
            }
        }

        const al::BigInt& getX() const noexcept {
            return y_;
        }
    private:
        const al::BigInt x_;
        const al::BigInt y_;
        const al::BigInt w_;
    };

} // anon namespace

auto alglib::cipolla(al::BigInt a, al::BigInt p, int countRetry) -> std::pair<al::BigInt, al::BigInt> {
    if (countRetry == 0) {
        return { -1, -1 };
    }

    al::BigInt b = randBigInt(p);
    al::BigInt w = (b * b - a) % p;

    if (w == 0) {
        al::BigInt c = p - b;
        if (c < b) {
            std::swap(c, b);
        }
        return { b, c };
    }
    else {
        if (legendre(w, p) == 1) {
            return cipolla(a, p, countRetry - 1);
        }
        else {
            SqrtPolynomial tmp(b, 1, w);
            SqrtPolynomial ans1_ = tmp.pow((p + 1) / 2, p);
            al::BigInt ans1 = ans1_.getX() % p;
            al::BigInt ans2 = p - ans1;
            if (ans2 < ans1) {
                std::swap(ans1, ans2);
            }
            return { ans1, ans2 };
        }
    }
}