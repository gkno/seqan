#ifndef SEQAN_HEADER_FIND_MOTIF_BASE_H
#define SEQAN_HEADER_FIND_MOTIF_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

///////////////////////////////////////////////////////////////////////////////////////////////////

/**
.Class.MotifFinder:
..summary:Holds the algorithm parameter values and the motif instance(s) found by the apprppriate
          motif discovery algorithm.
..cat:Motif Finding
..signature:MotifFinder<TValue, TSpec>
..param.TValue:The type of sequences to be analyzed.
...type:Spec.Dna
...type:Spec.AminoAcid
..param.TSpec:The motif finding algorithm to search with.
...type:Spec.Projection
...type:Spec.EPatternBranching
...type:Spec.PMS1
...type:Spec.PMSP
*/

template <typename TValue, typename TSpec>
class MotifFinder;

//////////////////////////////////////////////////////////////////////////////
//Metafunctions
//////////////////////////////////////////////////////////////////////////////

/*
.Metafunction.Value:
..summary:Returns the sequence type of a @Class.MotifFinder@ type.
..cat:Motif Finding
..signature:Value<TMotifFinder>::Type
..param.TMotifFinder:A @Class.MotifFinder@ type.
...type:Class.MotifFinder
..returns:The sequence type of $MotifFinder$, i.e. $TValue$ for $MotifFinder<TValue, TSpec>$.
*/

template<typename TValue, typename TSpec>
struct Value< MotifFinder<TValue, TSpec> >
{
	typedef TValue Type;
};
template<typename TValue, typename TSpec>
struct Value< MotifFinder<TValue, TSpec> const>
{
	typedef TValue const Type;
};


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

/**
.Function.findMotif:
..summary:Represents the main function which is used to start the search for noticeable motif patterns.
..cat:Motif Finding
..signature:findMotif(finder,dataset,seq_model)
..param.finder:The @Class.MotifFinder@ object.
...type:Class.MotifFinder
..param.dataset:The dataset object representing the input sequences.
...type:Class.StringSet
..param.seq_model:The seq_model object.
...type:Tag.OOPS
...type:Tag.OMOPS
...type:Tag.ZOOPS
...type:Tag.TCM
...remarks:The sequence models rely on different assumptions about the distribution of motif occurrences
           across the sample sequences. 
..remarks:The PROJECTION algorithm is able to run in OOPS, ZOOPS and TCM mode.
..remarks:The ePatternBranching algorithm is able to run in OOPS and OMOPS mode.
..remarks:The PMS1 and PMSP algorithm is able to run in OOPS, OMOPS, ZOOPS and TCM mode.
*/

/**
.Function.factorial:
..summary:Calculates the factorial value of any value.
..cat:Motif Finding
..signature:factorial(value)
..param.value:The value object.
*/

template<typename TType>
TType factorial(TType n)
{
	TType result = 0;

	if(n==0)
	{
		result = 1;
	}
	else
	{
		result = n*factorial(n-1);
	}
   
	return result;
};

//////////////////////////////////////////////////////////////////////////////

/**
.Function.binomialCoefficient:
..summary:Calculates the binomial coefficient C(n,k).
..cat:Motif Finding
..signature:binomialCoefficient(n,k)
*/

template<typename TType>
TType binomialCoefficient(TType n, TType k)
{
	TType result = 1;
	for(TType i=(n-k+1); i<=n; ++i)
	{
		result*=i;
	}
	result = result/factorial(k);
	
	return result;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.hammingDistance:
..summary:Determines the Hamming distance between two sequences.
..cat:Motif Finding
..signature:hammingDistance(start1,end1,start2)
..param.start1:An iterator pointing to the beginning of the first sequence which is either
              a @Shortcut.DnaString@ or a @Shortcut.Peptide@. 
...type:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
..param.end1:An iterator pointing to the end of the first sequence which is either
            a @Shortcut.DnaString@ or a @Shortcut.Peptide@. 
...type:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
..param.start2:An iterator pointing to the beginning of the second sequence which is either
              a @Shortcut.DnaString@ or a @Shortcut.Peptide@. 
...type:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
*/

template<typename TStringIterator>
size_t hammingDistance(TStringIterator start1, TStringIterator end1, TStringIterator start2)
{
	unsigned int num_of_mismatches = 0;
	while(start1!=end1)
	{
		if(*start1!=*start2)
		{
			++num_of_mismatches;
		}
		++start1;
		++start2;
	}

	return num_of_mismatches;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.inverseHash:
..summary:Determines the corresponding sequence given the hash value.
..cat:Motif Finding
..signature:inverseHash<TValue>(hash_value,alphabet_size,seq_size)
..param.hash_value:The hash_value object.
..param.alphabet_size:The alphabet_size object which is four for nucleotide sequences and twenty for amino acid sequences.
..param.seq_size:The seq_size object representing the size of the corresponding sequence.
..remarks:String<Dna> pattern = inverseHash<Dna>(hash_value,alphabet_size, seq_size)
*/

template<typename TValue, typename TType>
String<TValue>
inverseHash(TType const & hash_value, 
			typename Size<TValue>::Type const & alp_size, 
			typename Size< String<TValue> >::Type const & seq_size)
{
	typedef String<TValue> TString;
	TString seq;
	resize(seq, seq_size);

	TType hash_val = hash_value;
	typedef typename Position<TString>::Type TPos;
	for(TPos i=0; i<seq_size; ++i)
	{
		int letter = hash_val%alp_size;
		seq[i] = static_cast<TValue>(letter);
		hash_val = (hash_val-letter)/alp_size;
	}

	std::reverse(begin(seq), end(seq));
	return seq;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.displayResult:
..summary:Displays all found motif candidates. In the case of the Projection Motif Finder
          the function displays the consensus pattern of the found motif candidate.
..cat:Motif Finding
..signature:displayResult(motif_finder)
..param.motif_finder:The @Class.MotifFinder@ object.
...type:Class.MotifFinder
..param.dataset:The dataset object representing the input sequences.
...type:Class.StringSet
*/

template<typename TValue, typename TAlgorithm>
void
displayResult(MotifFinder<TValue, TAlgorithm> & finder)
{
	typedef String<TValue> TString;
	typedef String<TString> TStrings;

	if(length(finder.set_of_motifs)!=0)
	{
		unsigned int counter = 0;
		typename Iterator<TStrings>::Type iter = begin(finder.set_of_motifs);
		for(; !atEnd(iter, finder.set_of_motifs); goNext(iter))
		{
			std::cout << "[" << counter << "]: " << *iter << "\n";
			++counter;
		}
		std::cout << "\n";
	}
	else
	{
		std::cout << "NO MOTIF HAS BEEN FOUND!!!\n";
	}
}

/////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...