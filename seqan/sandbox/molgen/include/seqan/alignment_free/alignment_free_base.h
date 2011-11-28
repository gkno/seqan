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
// Definition of all AFScore and the specialisations, D2, D2star, D2z and N2
// ==========================================================================

#ifndef SANDBOX_ALIGNMENT_FREE_INCLUDE_SEQAN_ALIGNMENT_FREE_ALIGNMENT_FREE_BASE_H_
#define SANDBOX_ALIGNMENT_FREE_INCLUDE_SEQAN_ALIGNMENT_FREE_ALIGNMENT_FREE_BASE_H_

namespace seqan {

/**
.Class.AFScore:
..cat:Miscellaneous
..summary:Used to specify parameters and methods for alignment-free sequence comparison
..signature:AFScore<TScoreType>
..param.TScoreType:Method to use
*/

template <typename TSpec>
struct AFScore;


/**
.Spec.D2
..cat:Miscellaneous
..summary:D2 computes the inner product of the kmer count vectors
..signature:AFScore<D2>
..general:Class.AFScore
..remarks:D2 can be used for alignment-free sequence comparison
.Memvar.D2#kmerSize:
..class:Spec.D2
..summary:Size of the kmers
*/

/**.Memfunc.D2#AFScore<D2>:
..class:Spec.D2
..summary:Constructor
..signature:AFScore<D2>(kmerSize)
..param.kmerSize: Size of kmer
...type:nolink:unsigned
..remarks:
...text:No remarks
*/

struct D2_;     // Inner product of k-mer counts, d2 score
typedef Tag<D2_> const D2;

template <>
struct AFScore<D2>
{
    unsigned kmerSize;
    bool verbose;
    AFScore<D2>(unsigned k, bool verbose_ = false) : kmerSize(k), verbose(verbose_)
    {
    }

};


/**
.Spec.D2Star
..cat:Miscellaneous
..summary:D2Star computes the inner product of the standardised kmer count vectors
..signature:AFScore<D2Star>
..general:Class.AFScore
..remarks:D2Star can be used for alignment-free sequence comparison, this version calculates the background model on the concatenation of both sequences
..remarks:Reinert, G.; Chew, D.; Sun, F. & Waterman, M. S. Alignment-Free Sequence Comparison (I): Statistics and Power. J Comput Biol, 2009
.Memvar.D2Star#kmerSize:
..class:Spec.D2Star
..summary:Size of the kmers
.Memvar.D2Star#bgModelOrder:
..class:Spec.D2Star
..summary:Order of the background model
.Memvar.D2Star#outputFile:
..class:Spec.D2Star
..summary: When specified, all kmerWeights will be written to this file, for every sequence, and for every sequence comparison
*/


/**.Memfunc.D2Star#AFScore<D2Star>:
..class:Spec.D2Star
..summary:Constructor
..signature:AFScore<D2Star>(kmerSize, bgModelOrder)
..signature:AFScore<D2Star>(kmerSize, bgModelOrder, outputFile)
..param.kmerSize: Size of kmer
...type:nolink:unsigned
..param.bgModelOrder: Order of the background Markov model
...type:nolink:unsigned
...default:0
count variance
...type:nolink:bool
...default:false
..param.outputFile:When specified, all kmerWeights will be written to this file, for every sequence, and for every sequence comparison. The file can be very large.
...type:nolink:String<char>
...default:""
..remarks:
...text:No remarks
*/

struct D2Star_;        // Reinert and Waterman, D2 with centralised and standardised counts
typedef Tag<D2Star_> const D2Star;

template <>
struct AFScore<D2Star>
{
    unsigned kmerSize;
    unsigned bgModelOrder;
    bool verbose;


    AFScore<D2Star>(unsigned k, unsigned m, bool verbose_ = false) : kmerSize(k), bgModelOrder(m), verbose(verbose_)
    {
    }
};

/**
.Spec.N2
..cat:Miscellaneous
..summary:N2 computes the inner product of the standardised neighbourhood kmer count vectors
..signature:AFScore<N2>
..general:Class.AFScore
..remarks:N2 can be used for alignment-free sequence comparison
..remarks:Goeke et al. to appear
.Memvar.N2#kmerSize:
..class:Spec.N2
..summary:Size of the kmers
.Memvar.N2#bgModelOrder:
..class:Spec.N2
..summary:Order of the background model
.Memvar.N2#revCom:
..class:Spec.N2
..summary:Scoring of reverse complements words [max/min/mean/bothStrands]
.Memvar.N2#mismatches:
..class:Spec.N2
..summary:Approximate word matches [0(exact)/1(one mismatch)]
.Memvar.N2#mismatchWeight:
..class:Spec.N2
..summary:Weight for approximate word matches
.Memvar.D2star_old#outputFile:
..class:Spec.N2
*/


/**.Memfunc.N2#AFScore<N2>:
..class:Spec.N2
..summary:Constructor
..signature:AFScore<N2>(kmerSize, bgModelOrder,outputFile)
..param.kmerSize: Size of kmer
...type:nolink:unsigned
..param.bgModelOrder: Order of the background Markov model
...type:nolink:unsigned
...default:0
..param.outputFile:When specified, all normalised and standardised kmer neighbourhood counts will be written to this file for every sequence
...default:""
..remarks:
...text:No remarks
*/

struct N2_;     // Reinert and Waterman, D2 with centralised and standardised counts
typedef Tag<N2_> const N2;

template <>
struct AFScore<N2>
{
    unsigned kmerSize;
    unsigned bgModelOrder;
    String<char> revCom;    // Count reverse complement words?
    // revCom="";"mean","max","bothStrands"
    unsigned mismatches;    // Currently 0 or 1
    double mismatchWeight;  // Weight of words in the mismatch neighbourhood (0-1)
    bool verbose;
    bool norm;          // Normalize score? needed to provide a proper similarity measure
    String<char> outputFile;    // output of all pairwise kmerSores for pairwise visualisation
    
	//Constructor for the simple case with only exact word counts (N2*)
    AFScore<N2>(unsigned k, unsigned m, String<char> kmerWeightsFile = "", bool verbose_ = false)
    {
        kmerSize = k;
        bgModelOrder = m;
        outputFile = kmerWeightsFile;
	verbose = verbose_;
        revCom = "";
        mismatches = 0;
        mismatchWeight = 1.0;
        norm = true;

    };
	//Constructor for the case with exact word counts and reverse complement (N2rc)
    AFScore<N2>(unsigned k, unsigned m, String<char> revCom_, String<char> kmerWeightsFile = "", bool verbose_ = false)
    {
        kmerSize = k;
        bgModelOrder = m;
        revCom = revCom_;
	outputFile = kmerWeightsFile;
	verbose = verbose_;
        mismatches = 0;
        mismatchWeight = 1.0;
        norm = true;
    };
    	//Constructor for the case with mismatch-neighbourhood word counts and reverse complement (N2mmrc)
    AFScore<N2>(unsigned k, unsigned m, String<char> revCom_, unsigned mm, double mmw, String<char> kmerWeightsFile = "", bool verbose_ = false)
    {
        kmerSize = k;
        bgModelOrder = m;
        revCom = revCom_;
        mismatches = mm;
        mismatchWeight = mmw;
	outputFile = kmerWeightsFile;
        verbose = verbose_;
	norm = true;
    };
};


/**
.Spec.D2z
..cat:Miscellaneous
..summary:D2z computes a z-score of the inner product of kmer count vectors
..signature:AFScore<D2z>
..general:Class.AFScore
..remarks:D2z can be used for alignment-free sequence comparison. The algorithm differs from the original implementation by the way masked sequences are handled
..remarks:Kantorovitz, M. R.; Robinson, G. E. & Sinha, S. A statistical method for alignment-free comparison of regulatory sequences. Bioinformatics, 2007
.Memvar.D2z#kmerSize:
..class:Spec.D2z
..summary:Size of the kmers
.Memvar.D2z#bgModelOrder:
..class:Spec.D2z
..summary:Order of the background model
*/

/**.Memfunc.D2z#AFScore<D2z>:
..class:Spec.D2z
..summary:Constructor
..signature:AFScore<D2z>(kmerSize, bgModelOrder)
..param.kmerSize: Size of kmer
..param.bgModelOrder: Order of the background Markov model
..remarks:
...text:nolink:No remarks
*/

struct D2z_;        // Inner product of k-mer counts, d2 score with z-score
typedef Tag<D2z_> const D2z;

template <>
struct AFScore<D2z>
{
    unsigned kmerSize;
    unsigned bgModelOrder;
    bool verbose;
    AFScore<D2z>(unsigned k, unsigned m, bool verbose_ = false) : kmerSize(k), bgModelOrder(m), verbose(verbose_)
    {
    }
};

}  // namespace seqan

#endif  // SANDBOX_ALIGNMENT_FREE_INCLUDE_SEQAN_ALIGNMENT_FREE_ALIGNMENT_FREE_BASE_H_