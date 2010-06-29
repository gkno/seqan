/* Base calling for the 454 reads.
 */

#ifndef SIMULATE_454_BASE_CALLING_H_
#define SIMULATE_454_BASE_CALLING_H_

#include <cmath>

using namespace seqan;

class ThresholdMatrix
{
public:
    // The scaling parameter k.
    double _k;
    // Whether or not to use the sqrt for the std deviation computation.
    bool _useSqrt;
    // The edge length of the matrix.
    mutable unsigned _size;
    // The data of the matrix.
    mutable String<double> _data;

    ThresholdMatrix()
            : _k(0), _useSqrt(false), _size(0)
    {}
    
    ThresholdMatrix(double k, bool useSqrt)
            : _k(k), _useSqrt(useSqrt), _size(0)
    {}
};

inline double
normalDensityF(double x, double mu, double sigma)
{
//     std::cout << "normalDistF(" << x << ", " << mu << ", " << sigma << ")" << std::endl;
    const double PI = 3.14159265;
    double sigma2 = sigma * sigma;
    return exp(- (x - mu) * (x - mu) / (2 * sigma2)) / sqrt(2 * PI * sigma2);
}

inline double
lognormalDensityF(double x, double mu, double sigma)
{
//     std::cout << "lognormalDistF(" << x << ", " << mu << ", " << sigma << ")" << std::endl;
    if (x <= 0)
        return 0;
    const double PI = 3.14159265;
    double sigma2 = sigma * sigma;
    double log_mu2 = (log(x) - mu) * (log(x) - mu);
    return exp(-log_mu2 / (2 * sigma2)) / (x * sigma * sqrt(2 * PI));
}

inline double
dispatchDensityFunction(ThresholdMatrix const & matrix, unsigned r, double x)
{
    if (r == 0)
//         return lognormalDensityF(x, 0.23, 0.15);
        return normalDensityF(x, 0, 0.15);  // TODO(holtgrew): Using this for now until Huson answers.
    else
        return normalDensityF(x, r, (matrix._useSqrt ? sqrt(r) : r));
}

inline double
computeThreshold(ThresholdMatrix const & matrix, unsigned r1, unsigned r2)
{
//     std::cout << "computeThreshold(matrix, r1=" << r1 << ", r2=" << r2 << ")" << std::endl;
//     std::cout << "matrix._useSqrt == " << matrix._useSqrt << std::endl;
    if (r1 > r2)
        return computeThreshold(matrix, r2, r1);
    // The epsilon we use for convergence detection.
    const double EPSILON = 0.00001;

    // In i, we will count the number of iterations so we can limit the maximal
    // number of iterations.
    unsigned i = 0;
    (void) i;  // In case assertions are disabled.

    // f1 is the density function for r1 and f2 the density function for r2.

    // Pick left such that f1(left) > f2(left).
    double left = r1;
    if (left == 0) left = 0.23;
    while (dispatchDensityFunction(matrix, r1, left) <= dispatchDensityFunction(matrix, r2, left)) {
//         std::cout << "r1 = " << r1 << ", left = " << left << ", r2 = " << r2 << std::endl;
//         std::cout << "dispatchDensityFunction(matrix, r1, left)  == " << dispatchDensityFunction(matrix, r1, left) << ", dispatchDensityFunction(matrix, r2, left) == " <<  dispatchDensityFunction(matrix, r2, left) << std::endl;
        left /= 2.0;
    }
    // And pick right such that f1(right) < f2(right).
    double right = r2;
    if (right == 0) right = 0.5;
    while (dispatchDensityFunction(matrix, r1, right) >= dispatchDensityFunction(matrix, r2, right)) {
//         std::cout << "r1 = " << r1 << ", left = " << right << ", r2 = " << r2 << std::endl;
//         std::cout << "dispatchDensityFunction(matrix, r1, right)  == " << dispatchDensityFunction(matrix, r1, right) << ", dispatchDensityFunction(matrix, r2, right) == " <<  dispatchDensityFunction(matrix, r2, right) << std::endl;
        right *= 2.;
    }

    // Now, search for the intersection point.
    while (true) {
        SEQAN_ASSERT_LT_MSG(i++, 1000u, "Too many iterations (%u)! r1 = %u, r2 = %u.", i, r1, r2);
//         std::cout << "i == " << i << std::endl;

        double center = (left + right) / 2;
//         std::cout << "k == " << matrix._k << std::endl;
//         std::cout << "left = " << left << ", right = " << right << ", center = " << center << std::endl;
        double fCenter1 = dispatchDensityFunction(matrix, r1, center);
        double fCenter2 = dispatchDensityFunction(matrix, r2, center);
//         std::cout << "fCenter1 == " << fCenter1 << ", fCenter2 == " << fCenter2 << std::endl;
        double delta = fabs(fCenter1 - fCenter2);
        if (delta < EPSILON) {
//             std::cout << delta << " done" << std::endl;
            return center;
        }

        if (fCenter1 < fCenter2)
            right = center;
        else
            left = center;
    }
}

inline void
extendThresholds(ThresholdMatrix const & matrix, unsigned dim)
{
    // Allocate new data array for matrix.  Then compute values or copy
    // over existing ones.
    String<double> newData;
    resize(newData, dim * dim, Exact());
    for (unsigned i = 0; i < dim; ++i) {
        for (unsigned j = 0; j < dim; ++j) {
            if (i == j)
                continue;
            if (i < matrix._size && j < matrix._size)
                newData[i * dim + j] = matrix._data[i * matrix._size + j];
            else
                newData[i * dim + j] = computeThreshold(matrix, i, j);
        }
    }
    // Update matrix.
    assign(matrix._data, newData);
    matrix._size = dim;
}

inline double
getThreshold(ThresholdMatrix const & matrix, unsigned r1, unsigned r2)
{
//     std::cout << "getThreshold(matrix, r1=" << r1 << ", r2=" << r2 << ")" << std::endl;
    if (matrix._size <= r1 || matrix._size <= r2)
        extendThresholds(matrix, _max(r1, r2) + 1);
    return matrix._data[r1 * matrix._size + r2];
}

inline void
setK(ThresholdMatrix & matrix, double k)
{
    matrix._k = k;
}

inline void
setUseSqrt(ThresholdMatrix & matrix, bool useSqrt)
{
    matrix._useSqrt = useSqrt;
}

#endif  // SIMULATE_454_BASE_CALLING_H_
