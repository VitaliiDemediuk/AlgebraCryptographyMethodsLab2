/// My header
#include "AlgLib.h"

/// STL
#include <map>

using namespace alglib;

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