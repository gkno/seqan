// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Jonathan Goeke <goeke@molgen.mpg.de>
// ==========================================================================
// This file contains helper functions to count words in sequences and to 
// calculate probabilities and variances of word occurences.
// ==========================================================================
#ifndef SANDBOX_ALIGNMENT_FREE_INCLUDE_SEQAN_ALIGNMENT_FREE_KMER_FUNCTIONS_H_
#define SANDBOX_ALIGNMENT_FREE_INCLUDE_SEQAN_ALIGNMENT_FREE_KMER_FUNCTIONS_H_

namespace seqan {

template <typename TAlphabet>
struct _UnmaskedAlphabet
{
    typedef TAlphabet Type;
};

template <>
struct _UnmaskedAlphabet<Dna5>
{
    typedef Dna Type;
};

template <typename TAlphabet>
struct _UnmaskedAlphabet<const TAlphabet>
{
    typedef const typename _UnmaskedAlphabet<TAlphabet>::Type Type;
};

/*
 * Function to count kmers, Ns are not considered
 */
template <typename TString>
void countKmers(String<unsigned int> & kmerCounts, TString const & sequence, unsigned int const k)
{
    typedef typename Value<TString>::Type TAlphabet;
    typedef typename _UnmaskedAlphabet<TAlphabet>::Type TUnmaskedAlphabet;

    typedef typename Iterator<TString>::Type                    TIterator;
    typedef typename Position<TIterator>::Type                      TPosition;
    typedef Shape<TUnmaskedAlphabet, SimpleShape>                           TShape;

    // Declare variables
    TShape myShape;     // Shape, length can be changed (kmer_length)

    // int kmer_length=5;
    resize(myShape, k);

    // Calculate number of kmers/ length of count vector, respectively background vector
    int kmerNumber = _intPow((unsigned)ValueSize<TUnmaskedAlphabet>::VALUE, weight(myShape));
// arrayFill(kmerCounts, kmerNumber, 0);
    resize(kmerCounts, kmerNumber, 0);
    TIterator itSequence = begin(sequence);
    int counterN = 0;

    // Check for any N that destroys the first kmers
    unsigned int j = (k - 1);
    for (TPosition i = position(itSequence); i <= j; ++i)
    {
        if (_repeatMaskValue(sequence[i]))
        {
            counterN = i + 1;
        }
    }

    for (; itSequence <= (end(sequence) - k); ++itSequence)
    {
        // Check if there is a "N" at the end of the new kmer
        if (_repeatMaskValue(value(itSequence + (k - 1))))
            counterN = k; // do not consider any kmer covering this "N"

        // If there is no "N" overlapping with the current kmer, count it
        if (counterN <= 0)
        {
            unsigned hashValue = hash(myShape, itSequence);
            ++kmerCounts[hashValue];
        }
        counterN--;
    }
}

/*
 * Function to count kmers and background nucleotide frequencies, Ns are not considered
 * (for zero order background model)
 */
template <typename TValueBG, typename TString>
void countKmers(String<unsigned int> & kmerCounts, String<TValueBG> & backgroundFrequencies, TString const & sequence, unsigned int k)
{
    // besser TString, dann wird alles abgedeckt, z.b. String<Dna,array>, oder segment<> wie infix,...
    typedef typename Value<TString>::Type TAlphabet;
    typedef typename _UnmaskedAlphabet<TAlphabet>::Type     TUnmaskedAlphabet;
    // typedef typename Value<TStringBG>::Type               TValueBG;

    typedef typename Iterator<TString>::Type                   TIterator;
    // typedef typename Iterator<TStringBG >::Type		TIteratorTStringBG;
    typedef typename Iterator<String<TValueBG> >::Type     TIteratorTStringBG;

    typedef typename Position<TIterator>::Type                      TPosition;

    typedef Shape<TUnmaskedAlphabet, SimpleShape>           TShape;

    unsigned int alphabetSize = ValueSize<TUnmaskedAlphabet>::VALUE;

    // Declare variables
    TShape myShape;     // Shape, length can be changed (kmer_length)
    TShape myShapeBG;       // Shape for background, set to markovlen+1, here zero order only

    resize(myShape, k);
    resize(myShapeBG, 1);    // Markov model of zero order (count background frequencies)

    // Calculate number of kmers/ length of count vector, respectively background vector
    unsigned int kmerNumber = _intPow(alphabetSize, k);
    unsigned int kmerNumberBG = alphabetSize; // zero order model for DNA sequences

    // arrayFill(kmerCounts, kmerNumber, 0);
    // arrayFill(backgroundFrequencies, kmerNumberBG, (TValueBG) 0);
    resize(kmerCounts, kmerNumber, 0);
    resize(backgroundFrequencies, kmerNumberBG, (TValueBG) 0);

    TIterator itSequence = begin(sequence);

    int counterN = 0; // counter that counts how many kmers are effected by a N
    int counterNbg = 0;   // count for bg (different shape size)

    // Check for any N that destroys the first kmers
    unsigned int j = (k - 1);
    for (TPosition i = position(itSequence); i <= j; ++i)
    {
        if (_repeatMaskValue(sequence[i]))
            counterN = i + 1;
    }

    int sumBG = 0;    // count number of nucleotides for frequency calculation (Ns are not considered anymore!)
    for (; itSequence <= (end(sequence) - k); ++itSequence)
    {
        // Check if there is a "N" at the end of the new kmer
        if (_repeatMaskValue(value(itSequence + (k - 1))))
        {
            counterN = k; // do not consider any kmer covering this "N"
        }
        // If there is no "N" overlapping with the current kmer, count it
        if (counterN <= 0)
        {
            unsigned hashValue = hash(myShape, itSequence);
            ++kmerCounts[hashValue];
        }

        // Check if there is a "N" at the end of the new background word, here single letters only
        if (_repeatMaskValue(value(itSequence)))
        {
            counterNbg = 1;
        }
        if (counterNbg <= 0)
        {
            unsigned hashValueBG = hash(myShapeBG, itSequence);
            backgroundFrequencies[hashValueBG] += 1.0;
            ++sumBG;
        }
        counterN--;
        counterNbg--;
    }
    // The background counts are updated until the last base is covered
    for (; itSequence < end(sequence); ++itSequence)
    {
        if (_repeatMaskValue(value(itSequence)))
        {
            counterNbg = 1;
        }
        if (counterNbg <= 0)
        {
            unsigned hashValueBG = hash(myShapeBG, itSequence);
            ++backgroundFrequencies[hashValueBG];
            ++sumBG;
        }
        counterNbg--;
    }

    // Normalise the background counts to frequencies
    TIteratorTStringBG itBackground = begin(backgroundFrequencies);
    for (; itBackground < end(backgroundFrequencies); ++itBackground)
    {
        value(itBackground) /= ((TValueBG) sumBG);
    }
}

/*
 * Function to count kmers and build a background markov model with masked sequences
 *
 */
template <typename TString, typename TAlphabetBG, typename TValue>
void countKmers(String<unsigned int> & kmerCounts, MarkovModel<TAlphabetBG, TValue> & backgroundModel, TString & sequence, unsigned int k)
{
    typedef typename Value<TString>::Type TAlphabet;
    typedef typename _UnmaskedAlphabet<TAlphabet>::Type     TUnmaskedAlphabet;

    typedef typename Iterator<TString, Rooted>::Type        TIterator;
    typedef typename Iterator<String<int>, Rooted>::Type        TIteratorInt;
    typedef typename Position<TIterator>::Type                  TPosition;
    typedef Shape<TAlphabetBG, SimpleShape> TShape;

    // Declare variables
    TShape myShape;     // Shape, length can be changed (kmer_length)

    // int length=length(sequence);
    resize(myShape, k);

    // we only consider kmers without N
    int kmerNumber = _intPow((unsigned)ValueSize<TAlphabetBG>::VALUE, weight(myShape));
    // arrayFill(kmerCounts, kmerNumber, 0);
    resize(kmerCounts, kmerNumber, 0);

    // Create sequence set for the markov model, if Ns occur, the sequence is split and Ns are removed
    StringSet<String<TAlphabetBG> > seqSetMM;

    TIterator itSeq = begin(sequence);

    // Check for any N that destroys the first kmers
    unsigned int j = (k - 1);
    for (TPosition i = position(itSeq); i <= j; ++i)
    {
        // if(_repeatMaskValue(sequence[i]))

        if (_repeatMaskValue(sequence[i]))
        {
            // std::cout<<"N detected in begin of sequence";
            if ((i - position(itSeq)) > 0)
            {
                appendValue(seqSetMM, infix(sequence, position(itSeq), i));
            }
            goFurther(itSeq, i + 1 - position(itSeq));
            j = i + k - 1;
        }
    }

    int counterN = 0;

    TPosition startSplitSequence = position(itSeq); // to split sequences, position of possible start of a sequence after NNs

    for (; itSeq <= (end(sequence) - k); ++itSeq)
    {
        if (_repeatMaskValue(value(itSeq + (k - 1))))
        {
            counterN = k;

            if (((position(itSeq) + k - 1) > startSplitSequence))
                appendValue(seqSetMM, infix(sequence, startSplitSequence, (position(itSeq) + k - 1))); // infix(sequence2,position(itSeq2),i));

            startSplitSequence = (position(itSeq) + k); // position after n, possible start
        }
        if (counterN <= 0)
        {
            unsigned hashValue = hash(myShape, itSeq);
            ++kmerCounts[hashValue];
        }

        counterN--;
    }
    // string set for markov model
    if ((position(itSeq) + k - 1) > startSplitSequence)
    {
        appendValue(seqSetMM, infix(sequence, startSplitSequence, (position(itSeq) + k - 1))); // infix(sequence2,position(itSeq2),i));
    }

    // Build background model
    buildMarkovModel(backgroundModel, seqSetMM);


}

template <typename TValue, typename TString, typename TStringBG>
void calculateProbability(TValue & probability, TString & sequence, TStringBG & backgroundFrequencies)
{
    typedef typename Value<TString>::Type TAlphabet;
    typedef typename Iterator<TString, Rooted>::Type    TIteratorTString;
    TIteratorTString itSequence = begin(sequence);
    probability = (TValue) 1;
    for (; itSequence < end(sequence); ++itSequence)
        probability *= backgroundFrequencies[ordValue(*itSequence)];
}

template <typename TValue, typename TString, typename TStringBG>
void calculateVariance(TValue & variance, TString & word, TStringBG & backgroundFrequencies, int n)
{
    typedef typename Value<TString>::Type TAlphabet;
    typedef typename Value<TStringBG>::Type TValueBG;

    // Probability for word1
    typedef typename Iterator<String<int>, Rooted>::Type        TIteratorInt;
    typedef typename Iterator<String<TAlphabet>, Rooted>::Type      TIterator;


    int l = length(word);
    TValueBG p_w;
    calculateProbability(p_w, word, backgroundFrequencies);

    String<int> periodicity;
    calculatePeriodicity(periodicity, word, word);
    variance = (TValue) (n - l + 1) * p_w;
    for (TIteratorInt i = begin(periodicity); i < end(periodicity); ++i)
    {
        TValueBG p_clump;
        TValueBG p_tmp;
        calculateProbability(p_tmp, word, backgroundFrequencies);
        String<TAlphabet> wordPrefix = prefix(word, value(i));
        calculateProbability(p_clump, wordPrefix, backgroundFrequencies);
        p_clump *= p_tmp;
        variance += (TValue) 2 * (n - l + 1 - value(i)) * p_clump;

        // std::cout<<"\n"<<p_clump;
    }
    // std::cout<<"\np_w:"<<p_w;
    variance += (TValue) p_w * p_w * (n - 2 * n * l + 3 * l * l - 4 * l + 1);
    // std::cout<<word<<"\np_w:"<<p_w<<", var: "<<variance<<"\n";
}

template <typename TValue, typename TString, typename TStringBG>
void calculateCovariance(TValue & covariance, TString & word1, TString & word2, TStringBG & backgroundFrequencies, int n)
{
    if (word1 == word2)
    {
        calculateVariance(covariance, word1, backgroundFrequencies, n);
        return;
    }

    typedef typename Value<TString>::Type TAlphabet;
    typedef typename Value<TStringBG>::Type TValueBG;

    // Probability for word1
    typedef typename Iterator<String<int>, Rooted>::Type        TIteratorInt;
    typedef typename Iterator<String<TAlphabet>, Rooted>::Type      TIterator;


    covariance = 0;

    int l1 = length(word1);
    TValueBG p_w1;
    calculateProbability(p_w1, word1, backgroundFrequencies);

    String<int> periodicity1;
    calculatePeriodicity(periodicity1, word1, word2);
    for (TIteratorInt i = begin(periodicity1); i < end(periodicity1); ++i)
    {
        TValueBG p_clump;
        TValueBG p_tmp;
        calculateProbability(p_tmp, word2, backgroundFrequencies);
        String<TAlphabet> wordPrefix = prefix(word1, value(i));
        calculateProbability(p_clump, wordPrefix, backgroundFrequencies);
        p_clump *= p_tmp;
        covariance += (TValue) (n - l1 + 1 - value(i)) * p_clump;
        // std::cout<<"\n"<<p_clump;
    }

    int l2 = length(word2);
    TValueBG p_w2;
    calculateProbability(p_w2, word2, backgroundFrequencies);

    String<int> periodicity2;
    calculatePeriodicity(periodicity2, word2, word1);
    for (TIteratorInt i = begin(periodicity2); i < end(periodicity2); ++i)
    {

        TValueBG p_clump;
        TValueBG p_tmp;
        calculateProbability(p_tmp, word1, backgroundFrequencies);
        String<TAlphabet> wordPrefix = prefix(word2, value(i));
        calculateProbability(p_clump, wordPrefix, backgroundFrequencies);
        p_clump *= p_tmp;
        covariance += (TValue) (n - l2 + 1 - value(i)) * p_clump;
        // std::cout<<"\n"<<p_clump;
    }

    // std::cout<<"\np_w:"<<p_w;
    covariance += (TValue) p_w1 * p_w2 * (n - 2 * n * l1 + 3 * l1 * l1 - 4 * l1 + 1);
    // std::cout<<word1<<"\np_w:"<<p_w1<<", covar: "<<covariance<<"\n";
}

template <typename TValue, typename TSpec, typename TAlphabet>
void calculateVariance(TValue & variance, String<TAlphabet, TSpec> & word, MarkovModel<TAlphabet, TValue> & bgModel, int n)
{
    // typedef typename Value<TString>::Type TAlphabet;
    // typedef typename Value<TStringBG>::Type TValueBG;

    // Probability for word1
    typedef typename Iterator<String<int>, Rooted>::Type        TIteratorInt;
    typedef typename Iterator<String<TAlphabet>, Rooted>::Type      TIterator;

    int l = length(word);
    TValue p_w;
    p_w = emittedProbability(bgModel, word);
// std::cout<<"wordp:"<<p_w;
    String<int> periodicity;
    calculatePeriodicity(periodicity, word, word);
    variance = (TValue) (n - l + 1) * p_w;
    for (TIteratorInt i = begin(periodicity); i < end(periodicity); ++i)
    {
// std::cout<<"\ni:"<<value(i)<<"\n";
        TValue p_clump;
        // TValue p_tmp;

// String<TAlphabet> clump=word;

        // calculateProbability(p_clum,word,backgroundFrequencies);
//			String<TAlphabet> wordPrefix = prefix(word,value(i));
// append(wordPrefix,clump);
        String<TAlphabet> clump = prefix(word, value(i));
        append(clump, word);
// std::cout<<"\nclump:"<<clump<<"\n";
        // calculateProbability(p_clump,wordPrefix,backgroundFrequencies);
        // p_clump*=p_tmp;
        p_clump = emittedProbability(bgModel, clump);
        variance += (TValue) 2 * (n - l + 1 - value(i)) * p_clump;
        // std::cout<<"\n"<<clump<<p_clump;
    }
    // std::cout<<"\np_w:"<<p_w;
    variance += (TValue) p_w * p_w * (n - 2 * n * l + 3 * l * l - 4 * l + 1);
    // std::cout<<word<<"\np_w:"<<p_w<<", var: "<<variance<<"\n";
}
/*
// like the r code except the third term includes the differences for word occurences when they are so close that the stationary distribution cannot be assumed
template <typename TValue, typename TSpec, typename TAlphabet>
void calculateVarianceRobinBook(TValue & variance, String<TAlphabet, TSpec> & word, MarkovModel<TAlphabet, TValue> & bgModel, int l)
{
    typedef typename Iterator<String<int>, Rooted>::Type        TIteratorInt;
    typedef typename Iterator<String<TAlphabet>, Rooted>::Type      TIterator;


    unsigned k = length(word);
    TValue p_w;
    p_w = emittedProbability(bgModel, word);
    String<int> eps;
    calculateOverlapIndicator(eps, word, word);
    variance = (TValue) (l - k + 1) * p_w * (1 - p_w);
    for (unsigned d = 1; d <= k - 1; ++d)
    {
        String<TAlphabet> wordSuffix = suffix(word, k - d);
        TValue p_suffix = 1.0;
        Shape<TAlphabet, SimpleShape> orderShape;
        resize(orderShape, bgModel.order);
        for (unsigned j = (k - d + 1); j <= k; ++j)
        {
            p_suffix = p_suffix * value(bgModel.transition, hash(orderShape, begin(word) + j - 2), hash(orderShape, begin(word) + j - 1)); //1,1);//(infix(word,j-2,j-1)),(infix(word,j-1,j)));
        }

        variance += (TValue) 2.0 * p_w * (l - k - d + 1) * (eps[k - d - 1] * p_suffix - p_w);

    }
    // std::cout<<word<<"\np_w:"<<p_w<<", var: "<<variance<<"\n";
}*/

template <typename TValue, typename TSpec, typename TAlphabet>
void calculateCovariance(TValue & covariance, String<TAlphabet, TSpec> & word1, String<TAlphabet, TSpec> & word2, MarkovModel<TAlphabet, TValue> & bgModel, int n)
{
    if (word1 == word2)
    {
        calculateVariance(covariance, word1, bgModel, n);
        return;
    }

    // Probability for word1
    typedef typename Iterator<String<int>, Rooted>::Type        TIteratorInt;
    typedef typename Iterator<String<TAlphabet>, Rooted>::Type      TIterator;

    covariance = 0;

    int l1 = length(word1);
    TValue p_w1;
    p_w1 = emittedProbability(bgModel, word1);

    String<int> periodicity1;
    calculatePeriodicity(periodicity1, word1, word2);  // word2 is right
    for (TIteratorInt i = begin(periodicity1); i < end(periodicity1); ++i)
    {
        TValue p_clump;
        String<TAlphabet> clump = prefix(word1, value(i));
        append(clump, word2);

        p_clump = emittedProbability(bgModel, clump);

        covariance += (TValue) (n - l1 + 1 - value(i)) * p_clump;
    }
    TValue p_w2;
    p_w2 = emittedProbability(bgModel, word2);


    String<int> periodicity2;
    calculatePeriodicity(periodicity2, word2, word1);
    for (TIteratorInt i = begin(periodicity2); i < end(periodicity2); ++i)
    {
        TValue p_clump;
        String<TAlphabet> clump = prefix(word2, value(i));
        append(clump, word1);

        p_clump = emittedProbability(bgModel, clump);

        covariance += (TValue) (n - l1 + 1 - value(i)) * p_clump;
        // std::cout<<"\n"<<p_clump;
    }

    // std::cout<<"\np_w:"<<p_w;
    covariance += (TValue) p_w1 * p_w2 * (n - 2 * n * l1 + 3 * l1 * l1 - 4 * l1 + 1);
    // std::cout<<word1<<"\np_w:"<<p_w1<<", covar: "<<covariance<<"\n";
}
/*
// like the R code without term 3, everything as in robin book
template <typename TValue, typename TSpec, typename TAlphabet>
void calculateCovarianceRobinBook(TValue & covariance, String<TAlphabet, TSpec> & word1, String<TAlphabet, TSpec> & word2, MarkovModel<TAlphabet, TValue> & bgModel, int l)
{
    if (word1 == word2)
    {
        calculateVarianceRobinBook(covariance, word1, bgModel, l);
        return;
    }

    // Probability for word1
    typedef typename Iterator<String<int>, Rooted>::Type        TIteratorInt;
    typedef typename Iterator<String<TAlphabet>, Rooted>::Type      TIterator;

    covariance = 0;

    int k = length(word1);
    TValue p_w1;
    p_w1 = emittedProbability(bgModel, word1);
    TValue p_w2;
    p_w2 = emittedProbability(bgModel, word2);

    String<int> eps1_2;
    calculateOverlapIndicator(eps1_2, word1, word2);
    String<int> eps2_1;
    calculateOverlapIndicator(eps2_1, word2, word1);

    covariance = 0;

    // First term
    for (int d = 1; d <= k - 1; ++d)
    {
        String<TAlphabet> wordSuffix = suffix(word2, k - d);
        TValue p_suffix = 1.0;
        Shape<TAlphabet, SimpleShape> orderShape;
        resize(orderShape, bgModel.order);
        for (int j = (k - d + 1); j <= k; ++j)
        {
            p_suffix = p_suffix * value(bgModel.transition, hash(orderShape, begin(word2) + j - 2), hash(orderShape, begin(word2) + j - 1)); //1,1);//(infix(word,j-2,j-1)),(infix(word,j-1,j)));
        }

        covariance += (TValue) p_w1 * (l - k - d + 1) * (eps1_2[k - d - 1] * p_suffix - p_w2);

    }
    // std::cout<<"\nterm1:"<<covariance;
    // Second term
    for (int d = 1; d <= k - 1; ++d)
    {
        String<TAlphabet> wordSuffix = suffix(word1, k - d);
        TValue p_suffix = 1.0;
        Shape<TAlphabet, SimpleShape> orderShape;
        resize(orderShape, bgModel.order);
        for (int j = (k - d + 1); j <= k; ++j)
        {
            p_suffix = p_suffix * value(bgModel.transition, hash(orderShape, begin(word1) + j - 2), hash(orderShape, begin(word1) + j - 1)); //1,1);//(infix(word,j-2,j-1)),(infix(word,j-1,j)));
        }

        covariance += (TValue) p_w2 * (l - k - d + 1) * (eps2_1[k - d - 1] * p_suffix - p_w1);

    }
    // std::cout<<"\nterm1+term2:"<<covariance;
    covariance -= (l - k + 1) * p_w1 * p_w2;
}*/

// calculate word peridicity (indicator for overlaps)
// P(w1,w2)= all p where w2_j=w1_{j+p} for all j=1...k-p
template <typename TString>
void calculatePeriodicity(String<int> & periodicity, TString & word1, TString & word2)
{
    typedef typename Value<TString>::Type TAlphabet;
    // simple case for variance only: word1==word2
    typedef typename Iterator<TString, Rooted>::Type        TIterator;
    typedef typename Size<TString>::Type TSize;
    TIterator itWord1 = begin(word1);
    TSize length1 = length(word1);
    TSize length2 = length(word2);

    for (TSize i = 1; i < length1; ++i)
    {
        String<TAlphabet> my_suffix = suffix(word1, i);       // overlap of suffix of word1 with prefix of word2
        TSize my_min = std::min(length2, (length1 - i));
        String<TAlphabet> my_prefix = prefix(word2, my_min);
        if (my_suffix == my_prefix)
        {
            appendValue(periodicity, i);
            // std::cout<<"\nperiodicity: "<<i<<periodicity;
        }
    }
}

// calculate word overlap indiciator epsilon (robin et al) (indicator for overlaps)
// eps(w1,w2)= 1 where w2_j=w1_{j+p} for all j=1...k-p
template <typename TString>
void calculateOverlapIndicator(String<int> & epsilon, TString & word1, TString & word2)
{

    typedef typename Value<TString>::Type TAlphabet;
    // simple case for variance only: word1==word2
    typedef typename Iterator<TString, Rooted>::Type        TIterator;
    typedef typename Size<TString>::Type TSize;
    TIterator itWord1 = begin(word1);
    TSize length1 = length(word1);
    TSize length2 = length(word2);
    // fill(epsilon,length1,0);
    resize(epsilon, length1, 0);
    for (TSize i = 0; i < length1; ++i)
    {
        String<TAlphabet> my_suffix = suffix(word1, length1 - i - 1);     // overlap of suffix of word1 with prefix of word2
        TSize my_min = min(length2, i + 1);
        String<TAlphabet> my_prefix = prefix(word2, my_min);
        if (my_suffix == my_prefix)
            epsilon[i] = 1;
    }
}

template <typename TString>
void
stringToStringSet(StringSet<TString> & stringSet, TString const & sequence)
{
    resize(stringSet, 1);
    stringSet[0] = sequence;
}

void
stringToStringSet(StringSet<String<Dna> > & dnaStringSet, String<Dna5> const & sequence)
{
    typedef Iterator<String<Dna5> const, Rooted>::Type       TIterator;
    typedef Iterator<String<int>, Rooted>::Type     TIteratorInt;
    typedef Position<TIterator>::Type                   TPosition;
    typedef Shape<Dna, SimpleShape> TShape;

    int length1;
    length1 = length(sequence);

    TIterator itSeq = begin(sequence);
    // kmer_length=1;

    // Check for any N that destroys the first kmers
    unsigned int j = 0;
    for (TPosition i = position(itSeq); i <= j; ++i)
    {
        if (sequence[i] == 'N')
        {
            // std::cout<<"N detected in begin of sequence";
            if ((i - position(itSeq)) > 0)
                appendValue(dnaStringSet, infix(sequence, position(itSeq), i));
            goFurther(itSeq, i + 1 - position(itSeq));
            j = i;
        }
    }

    int counterN = 0;

    TPosition startSplitSequence = position(itSeq);     // to split sequences, position of possible start of a sequence after NNs

    for (; itSeq <= (end(sequence) - 1); ++itSeq)
    {
        if (value(itSeq) == 'N')
        {
            counterN = 1;

            if (((position(itSeq)) > startSplitSequence))
            {
                appendValue(dnaStringSet, infix(sequence, startSplitSequence, position(itSeq))); // infix(sequence2,position(itSeq2),i));
            }
            startSplitSequence = (position(itSeq) + 1); // position after n, possible start

        }
        counterN--;
    }
    // string set for markov model
    if (position(itSeq) > startSplitSequence)
        appendValue(dnaStringSet, infix(sequence, startSplitSequence, position(itSeq))); // infix(sequence2,position(itSeq2),i));
}

void
cutNs(String<Dna5> & sequenceCut, String<Dna5> const & sequence)
{
    typedef Iterator<String<Dna5> const, Rooted>::Type       TIterator;
    typedef Iterator<String<int>, Rooted>::Type     TIteratorInt;
    typedef Position<TIterator>::Type                   TPosition;
    typedef Shape<Dna5, SimpleShape>    TShape;

    sequenceCut = "";
    int length1;
    length1 = length(sequence);

    TIterator itSeq = begin(sequence);
    // kmer_length=1;

    // Check for any N that destroys the first kmers
    unsigned int j = 0;
    for (TPosition i = position(itSeq); i <= j; ++i)
    {
        if (sequence[i] == 'N')
        {
            // std::cout<<"N detected in begin of sequence";
            if ((i - position(itSeq)) > 0)
                sequenceCut += infix(sequence, position(itSeq), i); // appendValue(dnaStringSet,infix(sequence,position(itSeq),i));
            goFurther(itSeq, i + 1 - position(itSeq));
            j = i;
        }
    }

    int counterN = 0;

    TPosition startSplitSequence = position(itSeq);     // to split sequences, position of possible start of a sequence after NNs

    for (; itSeq <= (end(sequence) - 1); ++itSeq)
    {
        if (value(itSeq) == 'N')
        {
            counterN = 1;
            if (((position(itSeq)) > startSplitSequence))
                sequenceCut += infix(sequence, startSplitSequence, position(itSeq)); // appendValue(dnaStringSet,infix(sequence,startSplitSequence,position(itSeq)));// infix(sequence2,position(itSeq2),i));
            startSplitSequence = (position(itSeq) + 1); // position after n, possible start
        }
        counterN--;
    }
    // string set for markov model
    if (position(itSeq) > startSplitSequence)
    {
        sequenceCut += infix(sequence, startSplitSequence, position(itSeq)); // appendValue(dnaStringSet,infix(sequence,startSplitSequence,position(itSeq)));// infix(sequence2,position(itSeq2),i));
    }
}

}  // namespace seqan

#endif  // SANDBOX_ALIGNMENT_FREE_INCLUDE_SEQAN_ALIGNMENT_FREE_KMER_FUNCTIONS_H_
