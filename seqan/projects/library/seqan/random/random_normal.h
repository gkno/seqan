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

struct Normal {};

// ===========================================================================
// Classes
// ===========================================================================

template <>
class PDF<Normal>
{
public:
    double _mu;
    double _sigma;

    PDF(double mu, double sigma)
            : _mu(mu), _sigma(sigma)
    {
        SEQAN_CHECKPOINT;
    }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

template <>
struct Value<PDF<Normal> >
{
    typedef double Type;
};

template <>
struct Value<const PDF<Normal> > : Value<PDF<Normal> > {};

// ===========================================================================
// Functions
// ===========================================================================

/*
..summary:Pick a normally distributed random number.
*/
template <typename TRNG>
inline
typename Value<PDF<Normal> >::Type
pickRandomNumber(TRNG & rng, PDF<Normal> & pdf)
{
    SEQAN_CHECKPOINT;

    // Normal Distribution Heuristics, ported from Python.
    //
    // Kinderman and Monahan method. Reference: Kinderman, A.J. and
    // Monahan, J.F., "Computer generation of random variables using
    // the ratio of uniform deviates", ACM Trans Math Software, 3,
    // (1977), pp257-260.

    double z;
    while (true) {
        // TODO(holtgrew): Maybe better use a uniform distribution here?
//         double r1 = static_cast<double>(pickRandomNumber(rng) - InfimumValue<TRNG>::VALUE) / static_cast<double>(SupremumValue<TRNG>::VALUE - InfimumValue<TRNG>::VALUE);
//         double r2 = static_cast<double>(pickRandomNumber(rng) - InfimumValue<TRNG>::VALUE) / static_cast<double>(SupremumValue<TRNG>::VALUE - InfimumValue<TRNG>::VALUE);
        double r1 = static_cast<double>(pickRandomNumber(rng)) / static_cast<double>(SupremumValue<TRNG>::VALUE);
        double r2 = static_cast<double>(pickRandomNumber(rng)) / static_cast<double>(SupremumValue<TRNG>::VALUE);

        double u1 = r1;
        double u2 = 1 - r2;
        z = SEQAN_NV_MAGICCONST * (u1 - 0.5) / u2;
        double zz = z * z / 4.0;
        if (zz < -::std::log10(u2))
            break;
    }
    return pdf._mu + z * pdf._sigma;
}

}  // namespace seqan

#endif  // SEQAN_RANDOM_RANDOM_NORMAL_H_
