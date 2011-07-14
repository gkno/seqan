// ==========================================================================
//                               alignmentFree
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

#ifndef SANDBOX_ALIGNMENTFREE_INCLUDE_SEQAN_ALIGNMENTFREE_ALIGNMENTFREE_BASE_H_
#define SANDBOX_ALIGNMENTFREE_INCLUDE_SEQAN_ALIGNMENTFREE_ALIGNMENTFREE_BASE_H_

namespace seqan
{


/**
.Class.AF_Score:
..cat:Miscellaneous
..summary:Used to specify parameters and methods for alignment-free sequence comparison
..signature:AF_Score<TScoreType>
..param.TScoreType:Method to use
*/

template <typename TSpec>
struct AF_Score;




/**
.Spec.D2
..cat:Miscellaneous
..summary:D2 computes the inner product of the kmer count vectors
..signature:AF_Score<D2>
..general:Class.AF_Score
..remarks:D2 can be used for alignment-free sequence comparison
.Memvar.D2#kmerSize:
..class:Spec.D2
..summary:Size of the kmers
*/

/**.Memfunc.D2#AF_Score<D2>:
..class:Spec.D2
..summary:Constructor
..signature:AF_Score<D2>(kmerSize)
..param.kmerSize: Size of kmer
...type:nolink:unsigned
..remarks:
...text:No remarks
*/

struct D2_;     //Inner product of k-mer counts, d2 score
typedef Tag<D2_> const D2;

template <>
struct AF_Score<D2>
{
    unsigned kmerSize;

    AF_Score<D2>(unsigned k)
    {
        kmerSize = k;
    };

};

/**
.Spec.D2star
..cat:Miscellaneous
..summary:D2star computes the inner product of the standardised kmer count vectors
..signature:AF_Score<D2star>
..general:Class.AF_Score
..remarks:D2star can be used for alignment-free sequence comparison, this version calculates the background model on the individual sequences
..remarks:Reinert, G.; Chew, D.; Sun, F. & Waterman, M. S. Alignment-Free Sequence Comparison (I): Statistics and Power. J Comput Biol, 2009
.Memvar.D2star#kmerSize:
..class:Spec.D2star
..summary:Size of the kmers
.Memvar.D2star#bgModelOrder:
..class:Spec.D2star
..summary:Order of the background model
.Memvar.D2star#outputFile:
..class:Spec.D2star
..summary: When specified, all kmerWeights will be written to this file, for every sequence, and for every sequence comparison
*/


/**.Memfunc.D2star#AF_Score<D2star>:
..class:Spec.D2star
..summary:Constructor
..signature:AF_Score<D2star>(kmerSize, bgModelOrder)
..signature:AF_Score<D2star>(kmerSize, bgModelOrder, revCom)
..signature:AF_Score<D2star>(kmerSize, bgModelOrder, outputFile)
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
struct D2star_;     //Reinert and Waterman, D2 with centralised and standardised counts
typedef Tag<D2star_> const D2star;

template <>
struct AF_Score<D2star>
{
    unsigned kmerSize;
    unsigned bgModelOrder;
    String<char> revCom;    //Count reverse complement words?
    //revCom="";"mean","max", "min"
    String<char> outputFile;    //output of all pairwise kmerSores for pairwise visualisation
    AF_Score<D2star>(unsigned k, unsigned m, String<char> kmerWeightsFile = "")
    {
        kmerSize = k;
        bgModelOrder = m;
        outputFile = kmerWeightsFile;
        revCom = "";
    };
    AF_Score<D2star>(unsigned k, unsigned m, String<char> revComOption, String<char> kmerWeightsFile = "")
    {
        kmerSize = k;
        bgModelOrder = m;
        outputFile = kmerWeightsFile;
        revCom = revComOption;
    };
};

/**
.Spec.D2star_original
..cat:Miscellaneous
..summary:D2star_original computes the inner product of the standardised kmer count vectors
..signature:AF_Score<D2star_original>
..general:Class.AF_Score
..remarks:D2star_original can be used for alignment-free sequence comparison, this version calculates the background model on the concatenation of both sequences
..remarks:Reinert, G.; Chew, D.; Sun, F. & Waterman, M. S. Alignment-Free Sequence Comparison (I): Statistics and Power. J Comput Biol, 2009
.Memvar.D2star_original#kmerSize:
..class:Spec.D2star_original
..summary:Size of the kmers
.Memvar.D2star_original#bgModelOrder:
..class:Spec.D2star_original
..summary:Order of the background model
.Memvar.D2star_original#outputFile:
..class:Spec.D2star_original
..summary: When specified, all kmerWeights will be written to this file, for every sequence, and for every sequence comparison
*/


/**.Memfunc.D2star_original#AF_Score<D2star_original>:
..class:Spec.D2star_original
..summary:Constructor
..signature:AF_Score<D2star_original>(kmerSize, bgModelOrder)
..signature:AF_Score<D2star_original>(kmerSize, bgModelOrder, outputFile)
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

struct D2star_original_;        //Reinert and Waterman, D2 with centralised and standardised counts
typedef Tag<D2star_original_> const D2star_original;

template <>
struct AF_Score<D2star_original>
{
    unsigned kmerSize;
    unsigned bgModelOrder;
    bool norm; //shoud be set to false, however for testpurpose
    //String<char> revCom;	//Count reverse complement words?
    //revCom="";"mean","max", "min"
    String<char> outputFile;    //output of all pairwise kmerSores for pairwise visualisation
    AF_Score<D2star_original>(unsigned k, unsigned m, String<char> kmerWeightsFile = "")
    {
        kmerSize = k;
        bgModelOrder = m;
        outputFile = kmerWeightsFile;
        norm = false;
    };
};

/**
.Spec.N2
..cat:Miscellaneous
..summary:N2 computes the inner product of the standardised neighbourhood kmer count vectors
..signature:AF_Score<N2>
..general:Class.AF_Score
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
.Memvar.D2star#outputFile:
..class:Spec.N2
*/


/**.Memfunc.N2#AF_Score<N2>:
..class:Spec.N2
..summary:Constructor
..signature:AF_Score<N2>(kmerSize, bgModelOrder,outputFile)
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
struct N2_;     //Reinert and Waterman, D2 with centralised and standardised counts
typedef Tag<N2_> const N2;

template <>
struct AF_Score<N2>
{
    unsigned kmerSize;
    unsigned bgModelOrder;
    String<char> revCom;    //Count reverse complement words?
    //revCom="";"mean","max","bothStrands"
    unsigned mismatches;    //Currently 0 or 1
    double mismatchWeight;  //Weight of words in the mismatch neighbourhood (0-1)
    bool norm;          //Normaize score? needed to provide a proper similarity measure
    bool verbose;
    bool onlyPositiveKmers;     //should all kmers be considered or only positively scored kmers?
    String<char> outputFile;    //output of all pairwise kmerSores for pairwise visualisation
    AF_Score<N2>(unsigned k, unsigned m, String<char> kmerWeightsFile = "")
    {
        kmerSize = k;
        bgModelOrder = m;
        //exactVariance=exact;
        outputFile = kmerWeightsFile;
        revCom = "";
        mismatches = 0;
        onlyPositiveKmers = false;
        mismatchWeight = 1.0;
        norm = false;
        verbose = true;
    };
    AF_Score<N2>(unsigned k, unsigned m, String<char> revComOption, unsigned mm, double mmw, bool v, bool n, String<char> kmerWeightsFile = "")
    {
        kmerSize = k;
        bgModelOrder = m;
        //exactVariance=exact;
        outputFile = kmerWeightsFile;
        revCom = revComOption;
        mismatches = mm;
        onlyPositiveKmers = false;
        mismatchWeight = mmw;
        norm = n;
        verbose = v;
    };

};


/**
.Spec.D2z
..cat:Miscellaneous
..summary:D2z computes a z-score of the inner product of kmer count vectors
..signature:AF_Score<D2z>
..general:Class.AF_Score
..remarks:D2z can be used for alignment-free sequence comparison. The algorithm differs from the original implementation by the way masked sequences are handled
..remarks:Kantorovitz, M. R.; Robinson, G. E. & Sinha, S. A statistical method for alignment-free comparison of regulatory sequences. Bioinformatics, 2007
.Memvar.D2z#kmerSize:
..class:Spec.D2z
..summary:Size of the kmers
.Memvar.D2z#bgModelOrder:
..class:Spec.D2z
..summary:Order of the background model
*/

/**.Memfunc.D2z#AF_Score<D2z>:
..class:Spec.D2z
..summary:Constructor
..signature:AF_Score<D2z>(kmerSize, bgModelOrder)
..param.kmerSize: Size of kmer
..param.bgModelOrder: Order of the background Markov model
..remarks:
...text:nolink:No remarks
*/

struct D2z_;        //Inner product of k-mer counts, d2 score with z-score
typedef Tag<D2z_> const D2z;

template <>
struct AF_Score<D2z>
{
    unsigned kmerSize;
    unsigned bgModelOrder;
    AF_Score<D2z>(unsigned k, unsigned m)
    {
        kmerSize = k;
        bgModelOrder = m;
    };



};



/**
.Spec.MplusD
..cat:Miscellaneous
..summary:MplusD combines a Markov model with kmerCounts for alignment-free sequence comparison (Distance measure)
..signature:AF_Score<MplusD>
..general:Class.AF_Score
..remarks:MplusD can be used for alignment-free sequence comparison, this is a distance measure. The Algorithm is different from the program that was published since not all details are provided, however there is reason to believe that this implementation is more trustworthy.
..remarks:Dai, Q.; Yang, Y. & Wang, T. Markov model plus k-word distributions: a synergy that produces novel statistical measures for sequence comparison. Bioinformatics, 2008, 24, 2296-2302
.Memvar.MplusD#kmerSize:
..class:Spec.MplusD
..summary:Size of the kmers
.Memvar.MplusD#bgModelOrder:
..class:Spec.MplusD
..summary:Order of the background model
*/

/**.Memfunc.MplusD#AF_Score<MplusD>:
..class:Spec.MplusD
..summary:Constructor
..signature:AF_Score<MplusD>(kmerSize, bgModelOrder)
..param.kmerSize: Size of kmer
...type:nolink:unsigned
..param.bgModelOrder: Order of the background Markov model
...type:nolink:unsigned
..remarks:
...text:No remarks
*/
struct MplusD_;
typedef Tag<MplusD_> const MplusD;

template <>
struct AF_Score<MplusD>
{
    unsigned kmerSize;
    unsigned bgModelOrder;
    AF_Score<MplusD>(unsigned k, unsigned m)
    {
        kmerSize = k;
        bgModelOrder = m;
    };


};



}  // namespace seqan

#endif  // SANDBOX_ALIGNMENTFREE_INCLUDE_SEQAN_ALIGNMENTFREE_ALIGNMENTFREE_BASE_H_