#ifndef ACM2_ELGAMAL_H
#define ACM2_ELGAMAL_H

#include "AlgLib.h"

namespace alglib {

    class Point {
        friend class ElGamal;
    public:
        explicit Point() noexcept = default;
        Point(al::BigInt x, al::BigInt y) noexcept;
        explicit Point(bool isZero) noexcept;

        [[nodiscard]] const al::BigInt& x() const noexcept;
        [[nodiscard]] const al::BigInt& y() const noexcept;

        Point operator-() const noexcept;
        bool operator == (const Point& rhs) const noexcept;
        bool operator != (const Point& rhs) const noexcept;

    private:
        const al::BigInt x_{0};
        const al::BigInt y_{0};
        const bool isZero_{true};
    };

    class ElGamal {
    public:
        explicit ElGamal() noexcept;
        explicit ElGamal(al::BigInt p, al::BigInt a, al::BigInt b, al::BigInt n, Point P) noexcept;

        std::pair<Point, Point> encrypt(Point M, Point Y);
        Point decrypt(std::pair<Point, Point> c, al::BigInt k);
        Point addPoints(Point p1, Point p2);
        Point mulPoint(al::BigInt k, Point p1);
        Point getP();
        void testPoint(const Point& p1, std::ostream& out);
        Point getRandomPointOnCurve(int countRetry);
        [[nodiscard]] al::BigInt getRandomSecretKey() const noexcept;
    private:
        const al::BigInt p_;
        const al::BigInt a_, b_;
        const al::BigInt n_;
        Point P_;
    };

    al::BigInt hexToBigInt(const std::string& s);

}

std::ostream& operator << (std::ostream& os, alglib::Point bigInt);

#endif //ACM2_ELGAMAL_H
