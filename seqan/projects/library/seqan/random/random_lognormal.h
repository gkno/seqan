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
  Code for log-normally distributed random number generation.
  ===========================================================================
*/

#ifndef SEQAN_RANDOM_RANDOM_LOGNORMAL_H_
#define SEQAN_RANDOM_RANDOM_LOGNORMAL_H_

namespace seqan {

// ===========================================================================
// Forwards, Tags.
// ===========================================================================

struct Normal;

struct LogNormal {};
struct MuSigma {};
struct MeanStdDev {};

// ===========================================================================
// Classes
// ===========================================================================

template <>
class PDF<LogNormal>
{
public:
    PDF<Normal> _normalDist;

    PDF(double mu, double sigma, MuSigma const &)
            : _normalDist(mu, sigma)
    {
        SEQAN_CHECKPOINT;
    }

    PDF(double mean, double stddev, MeanStdDev const &)
            : _normalDist(::std::log(mean) - 0.5 * ::std::log(1.0 + stddev * stddev / mean / mean),
                          ::std::sqrt(::std::log(1.0 + stddev * stddev / mean / mean)))
    {
        SEQAN_CHECKPOINT;
    }

    PDF(double mu, double sigma)
            : _normalDist(mu, sigma)
    {
        SEQAN_CHECKPOINT;
    }
};

// ===========================================================================
// Metafunctions
// ===========================================================================

template <>
struct Value<PDF<LogNormal> >
{
    typedef double Type;
};

template <>
struct Value<const PDF<LogNormal> > : Value<PDF<LogNormal> > {};

// ===========================================================================
// Functions
// ===========================================================================

/*
..summary:Pick a log-normally distributed random number.
*/
template <typename TRandomNumberGenerator>
inline
typename Value<PDF<LogNormal> >::Type
pickRandomNumber(TRandomNumberGenerator & rng, PDF<LogNormal> & pdf)
{
    SEQAN_CHECKPOINT;
    return exp(pickRandomNumber(rng, pdf._normalDist));
}

}  // namespace seqan

#endif  // SEQAN_RANDOM_RANDOM_LOGNORMAL_H_
