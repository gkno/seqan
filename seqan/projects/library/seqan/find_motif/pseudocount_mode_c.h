#ifndef SEQAN_HEADER_PSEUDOCOUNT_MODE_C_H
#define SEQAN_HEADER_PSEUDOCOUNT_MODE_C_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// CMode
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.CMode:
..summary: Represents the C ("constant") computation scheme for handling "zero" probabilities.
..general:Class.Pseudocount
..cat:Motif Finding
..signature:Pseudocount<TValue, CMode>
..param.TValue:The type of sequence which is considered.
...type:Spec.Dna
...type:Spec.AminoAcid
..remarks:The pseudocount is identical for each residue (pseudocount = epsilon/alphabet_size).
*/

///.Class.Pseudocount.param.TSpec.type:Spec.CMode

struct _CMode;
typedef Tag<_CMode> CMode;


template<typename TValue>
class Pseudocount<TValue, CMode>
{

//_____________________________________________________________________________________________

public:
	double pseudocount;
	double epsilon;

//_____________________________________________________________________________________________

	Pseudocount():
		pseudocount(0),
		epsilon(0)
	{
	}
	Pseudocount(double epsilon_):
		epsilon(epsilon_),
		pseudocount(0)
	{
		_computePseudocount();
	}
	Pseudocount(Pseudocount const & other_):
		pseudocount(other_.pseudocount),
		epsilon(other_.epsilon)
	{
	}
	~Pseudocount()
	{
	}

	Pseudocount const &
	operator = (Pseudocount const & other_)
	{
		pseudocount = other_.pseudocount;
		epsilon = other_.epsilon;

		return *this;
	}

//_____________________________________________________________________________________________

private:
	inline void
	_computePseudocount() 
	{
		// alphabet_size = ValueSize<TValue>::VALUE
		pseudocount = 
			static_cast<double>(epsilon/ValueSize<TValue>::VALUE);
	}

//_____________________________________________________________________________________________
	
};

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

// Function.normalize (s. profile.h)

template<typename TProfile, typename TValue>
void 
normalize(TProfile & profile, Pseudocount<TValue, CMode> & mode)
{
	typedef typename Value<TProfile>::Type TFreqDist;
	typedef typename FrequencyType<TFreqDist>::Type TFrequencyType;

	typename Size<TProfile>::Type profile_size = length(profile);
	for(typename Position<TProfile>::Type i=0; 
		i<profile_size; 
		++i)
	{
		typename Iterator<TFreqDist>::Type fd_begin = begin(profile[i]);
		typename Iterator<TFreqDist>::Type fd_end = end(profile[i]);
		if(std::find(fd_begin, fd_end, 0)!= fd_end)
		{
			// N:=row sum
			TFrequencyType N = sum(profile[i]);

			// add pseudocounts
			for(typename Position<TFreqDist>::Type j=0; 
				j<length(profile[i]); 
				++j)
			{
				profile[i][j] = 
					static_cast<TFrequencyType>(profile[i][j]+mode.pseudocount)/
					static_cast<TFrequencyType>(N+mode.epsilon);
			}
		}
		else
		{
			normalize(profile[i]);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...