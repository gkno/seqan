/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de 
 ===========================================================================
  Copyright (C) 2010
  
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.
  
  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
 
 ===========================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
 ===========================================================================
  Code for normally distributed random number generation.
 ===========================================================================
*/

#ifndef SEQAN_RANDOM_RANDOM_NORMAL_H_
#define SEQAN_RANDOM_RANDOM_NORMAL_H_

namespace seqan {

// ===========================================================================
// Forwards, Tags.
// ===========================================================================

// Specialization Tag for normal distribution.
struct Normal {};

// ===========================================================================
// Classes
// ===========================================================================

/**
.Spec.Normal Pdf
..signature:Pdf<Normal>
..general:Class.Pdf
..summary:Normal probability density function.
..cat:Random
..include:seqan/random.h
*/
template <>
class Pdf<Normal>
{
public:
    double _mu;
    double _sigma;

/**
.Memfunc.Normal Pdf#Pdf
..class:Spec.Normal Pdf
..summary:Constructor for normal Pdf.
..signature:Pdf<Normal>(mu, sigma)
..param.mu:Mean of the normal distribution.
...type:nolink:double
..param.sigma:Standard deviation of the normal distribution.
...type:nolink:double
*/
    Pdf(double mu, double sigma)
            : _mu(mu), _sigma(sigma)
    {
        SEQAN_CHECKPOINT;
    }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

template <>
struct Value<Pdf<Normal> >
{
    typedef double Type;
};

template <>
struct Value<const Pdf<Normal> > : Value<Pdf<Normal> > {};

// ===========================================================================
// Functions
// ===========================================================================

/*
..summary:Pick a normally distributed random number.
*/
template <typename TRNG>
inline
typename Value<Pdf<Normal> >::Type
pickRandomNumber(TRNG & rng, Pdf<Normal> const & pdf)
{
    SEQAN_CHECKPOINT;

    // Normal Distribution Heuristics, ported from Python.
    //
    // Kinderman and Monahan method. Reference: Kinderman, A.J. and
    // Monahan, J.F., "Computer generation of random variables using
    // the ratio of uniform deviates", ACM Trans Math Software, 3,
    // (1977), pp257-260.

    double z;
    Pdf<Uniform<double> > pdfUniform(0, 1);
    while (true) {
        double u1 = pickRandomNumber(rng, pdfUniform);
        double u2 = 1 - pickRandomNumber(rng, pdfUniform);
        z = SEQAN_NV_MAGICCONST * (u1 - 0.5) / u2;
        double zz = z * z / 4.0;
        if (zz < -::std::log10(u2))
            break;
    }
    return pdf._mu + z * pdf._sigma;
}

}  // namespace seqan

#endif  // SEQAN_RANDOM_RANDOM_NORMAL_H_
