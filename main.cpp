#include <iostream>

#include "AlgLib.h"

template <typename K, typename V>
std::ostream& operator<<(std::ostream& out, const std::map<K, V>& v) {
    for (size_t i = 0; const auto& [p, count] : v) {
        for(V j = 0; j < count; ++j) {
            out << p << " ";
        }
        if (i < v.size() - 1) {
            out << " ";
        }
        ++i;
    }
    return out;
}

int main() {
    std::cout << alglib::babyStepGiantStepLog(5, 2, 9) << std::endl; /// 2
    std::cout << alglib::pollardRhoDivision(10) << std::endl; /// 2 or 5
    std::cout << alglib::pollardRhoFactorisation(12343422) << std::endl; /// 7 26 13 111 47
    std::cout << alglib::euler(33) << std::endl; /// 20
    std::cout << alglib::mobius(16) << " " << alglib::mobius(30) << " " << alglib::mobius(210) << std::endl; /// 0 -1 1

    return 0;
}
