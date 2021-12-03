/// My header
#include "ElGamal.h"

/// STL
#include <random>

/// Boost
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/random.hpp>

using namespace alglib;

namespace {
    al::BigInt randBigInt(const al::BigInt& n)  {
        using namespace boost::multiprecision;

        std::random_device rd;
        std::mt19937 gen(rd());
        boost::random::uniform_int_distribution<cpp_int> uid(cpp_int(2), cpp_int(n.toStdString()));
        return uid(gen).convert_to<std::string>();
    }

}

/// Point --------------------------------------------------------------------------------------------------------------

Point::Point(al::BigInt x, al::BigInt y) noexcept
    : x_{std::move(x)}, y_{std::move(y)}, isZero_{false} {}

Point::Point(bool isZero) noexcept
    : x_{0}, y_{0}, isZero_{isZero} {}

Point Point::operator-() const noexcept {
    if (isZero_) {
        return *this;
    } else {
        return {x_, -y_};
    }
}

const al::BigInt& Point::x() const noexcept {
    return x_;
}

const al::BigInt& Point::y() const noexcept {
    return y_;
}

bool Point::operator== (const Point& rhs) const noexcept {
    return x_ == rhs.x_ && y_ == rhs.y_;
}

bool Point::operator!= (const Point& rhs) const noexcept {
    return !(*this == rhs);
}

/// ElGamal ------------------------------------------------------------------------------------------------------------

ElGamal::ElGamal() noexcept
    : p_{(al::pow(al::BigInt(2), 128) - 3) / 76439},
      a_{hexToBigInt("DB7C2ABF62E35E668076BEAD2088")},
      b_{hexToBigInt("659EF8BA043916EEDE8911702B22")},
      P_{hexToBigInt("09487239995A5EE76B55F9C2F098"), hexToBigInt("A89CE5AF8724C0A23E0E0FF77500")},
      n_{hexToBigInt("DB7C2ABF62E35E7628DFAC6561C5")} {}

ElGamal::ElGamal(al::BigInt p, al::BigInt a, al::BigInt b, al::BigInt n, Point P) noexcept
        : p_{std::move(p)},
          a_{std::move(a)},
          b_{std::move(b)},
          P_{std::move(P)},
          n_{std::move(n)} {}

std::pair<Point, Point> ElGamal::encrypt(Point M, Point Y) {
    const al::BigInt r = randBigInt(n_ - 1) + 1;
    const Point d = mulPoint(r, Y);
    const Point g = mulPoint(r, P_);
    const Point h = addPoints(M, d);
    return { g, h };
}

Point ElGamal::decrypt(std::pair<Point, Point> c, al::BigInt k) {
    const Point s = mulPoint(k, c.first);
    Point s1 = -s;
    Point M = addPoints(s1, c.second);
    return M;
}

Point ElGamal::addPoints(Point p1, Point p2) {
    if (p1.isZero_) {
        return p2;
    }
    if (p2.isZero_) {
        return p1;
    }
    if (p1 != p2 && p1.x_ == p2.x_) {
        return Point(true);
    }
    if (p1 == p2 && p1.y_ == 0) {
        return Point(true);
    }

    al::BigInt lambda;
    if (p1 != p2 && p1.x_ != p2.x_) {
        al::BigInt x2x1 = al::minusMod(p2.x_, p1.x_, p_);
        lambda = al::divisionMod(al::minusMod(p2.y_, p1.y_, p_),x2x1, p_);
    }
    else {
        lambda = al::divisionMod(((p1.x_ * 3) % p_ * p1.x_) % p_ + a_, p1.y_ * 2, p_);
    }

    al::BigInt xNew = (lambda * lambda - (p1.x_ + p2.x_)) % p_;
    al::BigInt yNew = (lambda * (p1.x_ - xNew) - p1.y_) % p_;
    return Point(xNew, yNew);
}


Point ElGamal::mulPoint(al::BigInt k, Point p1) {
    if (p1.isZero_ || k == 0) {
        return Point(true);
    }
    if (k == 1) {
        return p1;
    }
    Point tmp = mulPoint(k / 2, p1);

    Point ans = addPoints(tmp, tmp);

    if (k % 2 == 1) {
        return addPoints(ans, p1);
    }
    return ans;
}


Point ElGamal::getRandomPointOnCurve(int countRetry = 5) {
    if (countRetry == 0) {
        return {-1, -1};
    }
    al::BigInt x = randBigInt(p_);
    al::BigInt v = (al::powMod(x, al::BigInt(3), p_) + a_ * x % p_ + b_) % p_;
    if (jacobi(v, p_) != 1) {
        return getRandomPointOnCurve(countRetry - 1);
    }
    auto y2 = cipolla(v, p_);
    if (y2.first != -1) {
        return {x, y2.first};
    }
    else {
        return getRandomPointOnCurve(countRetry - 1);
    }
}

Point ElGamal::getP() {
    return P_;
}

void ElGamal::testPoint(const Point& p1, std::ostream& out) {
    out << "testing of point" << std::endl;
    out << "x                = " << p1.x_ << std::endl;
    out << "y                = " << p1.y_ << std::endl;
    out << "y^2              = " << p1.y_ * p1.y_ % p_ << std::endl;
    // cout << "x^3 + a*x + b    = " << (p1.x.pow(3, p) + a * p1.x % p + b) % p << endl;
}

al::BigInt ElGamal::getRandomSecretKey() const noexcept {
    return randBigInt(p_);
}

/// functions ----------------------------------------------------------------------------------------------------------

std::ostream& operator << (std::ostream& os, const Point& p1) {
    os << "(" << p1.x() << ", " << p1.y() << ")";
    return os;
}

/// ostream operator << ------------------------------------------------------------------------------------------------

al::BigInt alglib::hexToBigInt(const std::string& s) {
    al::BigInt tmp = 1, ans = 0;
    for (int64_t i = static_cast<int64_t>(s.size()) - 1; i >= 0; i--) {
        const char ch = s[i];
        long long v = 0;
        if (ch >= '0' && ch <= '9') {
            v = ch - '0';
        } else {
            v = ch - 'A' + 10;
        }
        ans = ans + (tmp * v);
        tmp = tmp * 16;
    }
    return ans;
}