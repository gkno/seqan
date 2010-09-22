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

// Forward-declarations.
struct Normal;

// Specialization Tag for log-normal distribution.
struct LogNormal {};

/**
.Tag.MuSigma
..signature:MuSigma
..summary:Tag to specify that the given parameters are mu and sigma of the underlying normal distribution for lognormal distributions.
..cat:Random
..include:seqan/random.h
..see:Spec.Log-Normal PDF
..see:Tag.MeanStdDev
*/
struct MuSigma {};
/**
.Tag.MeanStdDev
..signature:MeanStdDev
..summary:Tag to specify that the given parameters are mean an standard deviation of the lognormal distribution.
..cat:Random
..include:seqan/random.h
..see:Spec.Log-Normal PDF
..see:Tag.MuSigma
*/
struct MeanStdDev {};

// ===========================================================================
// Classes
// ===========================================================================

/**
.Spec.Log-Normal PDF
..general:Class.PDF
..summary:Log-normal probability density function.
..remark:Note that you can construct this either with mu/sigma of the underlying normal distribution (default) or with the mean and standard deviation of the log-normal distribution.
..cat:Random
..include:seqan/random.h
*/
template <>
class PDF<LogNormal>
{
public:
    PDF<Normal> _normalDist;

/**
.Memfunc.Log-Normal PDF#PDF
..class:Spec.Log-Normal PDF
..summary:Constructor for log-normal PDF.
Log-normal PDFs can either be initialized by the mean and standard deviation of the underlying normal distribution or directly of the log-normal distribution.
..signature:PDF<LogNormal>(mu, sigma[, MuSigma()])
..signature:PDF<LogNormal>(mean, stdDev, MeanStdDev())
..param.mu:Mean of the underlying normal distribution.
...type:nolink:double
..param.sigma:Standard deviation of the underlying normal distribution.
...type:nolink:double
..param.mean:Mean of the log-normal distribution.
...type:nolink:double
..param.stdDev:Standard deviation of the log-normal distribution.
...type:nolink:double
..see:Tag.MuSigma
..see:Tag.MeanStdDev
*/
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

template <typename TRandomNumberGenerator>
inline
typename Value<PDF<LogNormal> >::Type
pickRandomNumber(TRandomNumberGenerator & rng, PDF<LogNormal> const & pdf)
{
    SEQAN_CHECKPOINT;
    return exp(pickRandomNumber(rng, pdf._normalDist));
}

}  // namespace seqan

#endif  // SEQAN_RANDOM_RANDOM_LOGNORMAL_H_
