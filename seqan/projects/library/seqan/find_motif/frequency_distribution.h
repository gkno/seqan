#ifndef SEQAN_HEADER_FREQUENCY_DISTRIBUTION_H
#define SEQAN_HEADER_FREQUENCY_DISTRIBUTION_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

/**
.Class.FrequencyDistribution:
..summary:Holds a collection of objects of a specific type, where each object represents
          the frequency (absolute or relative probability) of a particular residue which is a member
		  of a fixed sequence alphabet.
..cat:Motif Finding
..signature:MotifFinder<TValue[, TSpec]>
..param.TValue:The type of sequence which is considered.
...type:@Spec.Dna@ or @Spec.AminoAcid@
..param.TSpec:The type of probability distribution. 
...default: $double$
...type:$double$ for relative probabilities, $int$ for absolute probabilities
...remarks: It is preferable to use $double$.
..remarks:The number of objects in "FrequencyDistribution" equals the size of the sequence alphabet.
*/

template<typename TValue, typename TSpec = double>
class FrequencyDistribution
{
//____________________________________________________________________________________________

	enum { SIZE = ValueSize<TValue>::VALUE }; 

//____________________________________________________________________________________________

public:
	TSpec frequency_list[SIZE];

	// constructor & destructor
	FrequencyDistribution()
	{
		std::fill(&frequency_list[0], &frequency_list[SIZE], 0);
	}
	FrequencyDistribution(TValue const & letter_)
	{
		convertResidueToFrequencyDist(*this, letter_);
	}
	FrequencyDistribution(FrequencyDistribution const & other_)
	{
		std::copy(&other_.frequency_list[0], &other_.frequency_list[SIZE], this->frequency_list); 
	}
	~FrequencyDistribution()
	{
	}

	// overloading operators
	FrequencyDistribution & 
	operator = (FrequencyDistribution const & other_)
	{
		if(this!=&other_)
		{
			std::copy(&other_.frequency_list[0], &other_.frequency_list[SIZE], this->frequency_list);
		}
		return *this;
	}

	FrequencyDistribution &
	operator += (FrequencyDistribution const & other_)
	{
		for(size_t i=0; i<SIZE; ++i)
		{
			frequency_list[i]+=other_.frequency_list[i];
		}

		return *this;
	}

	friend FrequencyDistribution 
	operator + (FrequencyDistribution const & lhs_, FrequencyDistribution const & rhs_)
	{
		FrequencyDistribution ret(lhs_);
		ret+=rhs_;

		return ret;
	}

	FrequencyDistribution &
	operator -= (FrequencyDistribution const & other_)
	{
		for(size_t i=0; i<SIZE; ++i)
		{
			frequency_list[i]-=other_.frequency_list[i];
		}

		return *this;
	}

	friend FrequencyDistribution 
	operator - (FrequencyDistribution const & lhs_, FrequencyDistribution const & rhs_)
	{
		FrequencyDistribution ret(lhs_);
		ret-=rhs_;

		return ret;
	}

	template<typename TType>
	FrequencyDistribution &
	operator *= (TType value_)
	{
		for(size_t i=0; i<SIZE; ++i)
		{
			frequency_list[i]*= static_cast<TSpec>(value_);
		}

		return *this;
	}

	template<typename TType>
	friend FrequencyDistribution 
	operator * (FrequencyDistribution const & fd_, TType value_)
	{
		FrequencyDistribution ret(fd_);
		ret*=value_;

		return ret;
	}

	template<typename TPos>
	inline TSpec &
	operator [] (TPos index_)
	{
		return frequency_list[index_];
	}

	template<typename TPos>
	inline TSpec const & 
	operator [] (TPos const index_) const
	{
		return frequency_list[index_];
	}

    friend inline std::ostream & 
	operator << (std::ostream & ostr, FrequencyDistribution & fd_) 
	{ 
		size_t i;
		for(i=0; i<SIZE; ++i)
		{	
			ostr.width(15);
			ostr << std::left << fd_.frequency_list[i];
		}
		
		return ostr;  
	} 

//____________________________________________________________________________________________

};

//////////////////////////////////////////////////////////////////////////////
//Metafunctions
//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.Value:
..summary:Returns the sequence type of a @Class.FrequencyDistribution@ type.
..cat:Motif Finding
..signature:Value<TFrequencyDistribution>::Type
..param.TFrequencyDistribution:A @Class.FrequencyDistribution@ type.
...type:Class.FrequencyDistribution
..returns:The sequence type of $TFrequencyDistribution$, i.e. $TValue$ for $FrequencyDistribution<TValue, TSpec>$.
*/

template<typename TValue, typename TSpec>
struct Value< FrequencyDistribution<TValue, TSpec> >
{
	typedef TValue Type;
};
template<typename TValue, typename TSpec>
struct Value< FrequencyDistribution<TValue, TSpec> const>
{
	typedef TValue const Type;
};

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.FrequencyType:
..summary:Returns the type of probability distribution of a @Class.FrequencyDistribution@ type.
..cat:Motif Finding
..signature:FrequencyType<TFrequencyDistribution>::Type
..param.TFrequencyDistribution:A @Class.FrequencyDistribution@ type.
...type:Class.FrequencyDistribution
..returns:The probability type of $TFrequencyDistribution$, i.e. $TSpec$ for $FrequencyDistribution<TValue, TSpec>$.
*/

template<typename TFrequencyDistribution>
struct FrequencyType;

template<typename TValue, typename TSpec>
struct FrequencyType< FrequencyDistribution<TValue, TSpec> >
{
	typedef TSpec Type;
};
template<typename TValue, typename TSpec>
struct FrequencyType< FrequencyDistribution<TValue, TSpec> const>
{
	typedef TSpec const Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Size.param.T.type:Class.FrequencyDistribution
template<typename TValue, typename TSpec>
struct Size< FrequencyDistribution<TValue, TSpec> >
{
	typedef typename Size<TValue>::Type Type;
};
template<typename TValue, typename TSpec>
struct Size<FrequencyDistribution<TValue, TSpec> const>
{
	typedef typename Size<TValue const>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Position.param.T.type:Class.FrequencyDistribution

template<typename TValue, typename TSpec>
struct Position< FrequencyDistribution<TValue, TSpec> >
{
	typedef typename Position<TValue>::Type Type;
};
template<typename TValue, typename TSpec>
struct Position< FrequencyDistribution<TValue, TSpec> const>
{
	typedef typename Position<TValue const>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Iterator.param.T.type:Class.FrequencyDistribution

/*template <typename TValue, typename TSpec>
struct Iterator< FrequencyDistribution<TValue, TSpec> >
{
	typedef typename Iterator< FrequencyDistribution<TValue, TSpec>, Standard >::Type Type;
};*/

template<typename TValue, typename TSpec>
struct Iterator< FrequencyDistribution<TValue, TSpec>, Standard >
{
	typedef TSpec * Type;
};
template<typename TValue, typename TSpec>
struct Iterator< FrequencyDistribution<TValue, TSpec> const, Standard >
{
	typedef TSpec const * Type;
};

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

/**
.Function.absFreqOfLettersInSeq:
..summary:Counts the number of times each residue of a fixed sequence alphabet occurs in a given sequence.
..cat:Motif Finding
..signature:absFreqOfLettersInSeq(frequencies,start,end)
..param.frequencies:The @Class.FrequencyDistribution@ object which holds the calculated frequencies.
...type:Class.FrequencyDistribution
..param.start:An iterator pointing to the beginning of a given sequence which is either
              a string of @Spec.Dna@ or a string of @Spec.AminoAcid@. 
...type:Concept.Iterator Iterator
....remarks:Standard conform iterator
...type:Shortcut.DnaIterator
....remarks:Iterator for @Shortcut.DnaString@ (a string of @Spec.Dna@).
....see:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
....remarks:Iterator for @Shortcut.Peptide@ (a string of @Spec.AminoAcid@).
....see:Shortcut.PeptideIterator
..param.end:An iterator pointing to the end of a given sequence which is either
            a string of @Spec.Dna@ or a string of @Spec.AminoAcid@.  
...type:Concept.Iterator Iterator
....remarks:Standard conform iterator
...type:Shortcut.DnaIterator
....remarks:Iterator for @Shortcut.DnaString@ (a string of @Spec.Dna@).
....see:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
....remarks:Iterator for @Shortcut.Peptide@ (a string of @Spec.AminoAcid@).
....see:Shortcut.PeptideIterator
*/

template<typename TValue, typename TSpec, typename TSeqIter> 
void 
absFreqOfLettersInSeq(FrequencyDistribution<TValue, TSpec> & fd,
					  TSeqIter & seq_start,
					  TSeqIter & seq_end) 
{	
	while(seq_start!=seq_end)
	{
		++fd[(int)*seq_start];
		++seq_start;
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.absFreqOfLettersInSetOfSeqs:
..summary:Counts the number of times each residue of a fixed sequence alphabet occurs in a given set of sequences.
..cat:Motif Finding
..signature:absFreqOfLettersInSetOfSeqs(frequencies,start,end)
..param.frequencies:The @Class.FrequencyDistribution@ object which holds the calculated frequencies.
...type:Class.FrequencyDistribution
..param.start:An iterator pointing to the first sequence of a given set of sequences which is considered. 
...type:Concept.Iterator Iterator
....remarks:Standard conform iterator
..param.end:An iterator pointing to the last sequence of a given set of sequences which is considered. 
...type:Concept.Iterator Iterator
....remarks:Standard conform iterator
..remarks.text:This function is similar to @Function.absFreqOfLettersInSeq@ except that the function is performed
               on a set of sequences.
*/

template<typename TValue, typename TSpec, typename TIter>
void
absFreqOfLettersInSetOfSeqs(FrequencyDistribution<TValue, TSpec> & fd,
							TIter & seq_start,
							TIter & seq_end)
{
	while(seq_start!=seq_end)
	{
		absFreqOfLettersInSeq(fd, begin(*seq_start), end(*seq_start));
		++seq_start;
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.addValue:
..summary:Adds a value of a specific type to each element of a given @Class.FrequencyDistribution@ object.
..cat:Motif Finding
..signature:addValue(frequencies,value)
..param.frequencies:The @Class.FrequencyDistribution@ object which holds the calculated frequencies.
...type:Class.FrequencyDistribution
..param.value:The value object which is added to each element of a @Class.FrequencyDistribution@ object.
...remarks:The value object should be identical in type to the elements of the @Class.FrequencyDistribution@ object.
*/

template<typename TValue, typename TSpec, typename TType>
void 
addValue(FrequencyDistribution<TValue, TSpec> & fd, TType const & val)
{
	typedef typename Position< FrequencyDistribution<TValue, TSpec> >::Type TPos;
	for(TPos i=0; i<length(fd); ++i)
	{
		fd[i]+= static_cast<TSpec>(val);
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.backgroundFrequency:
..summary:Determines the background letter frequencies in a given dataset
..cat:Motif Finding
..signature:backgroundFrequency(frequencies,start,end)
..param.frequencies:The @Class.FrequencyDistribution@ object which holds the calculated frequencies.
...type:Class.FrequencyDistribution
..param.start:An iterator pointing to the first sequence of a given dataset (set of sequences) which is considered. 
...type:Concept.Iterator Iterator
....remarks:Standard conform iterator
..param.end:An iterator pointing to the last sequence of a given dataset (set of sequences) which is considered. 
...type:Concept.Iterator Iterator
....remarks:Standard conform iterator
*/

template<typename TValue, typename TSpec,typename TDatasetIter> 
void 
backgroundFrequency(FrequencyDistribution<TValue, TSpec> & fd,
					TDatasetIter dataset_start,
					TDatasetIter dataset_end)
{
	absFreqOfLettersInSetOfSeqs(fd, dataset_start, dataset_end);

	// check for zero entries
	if(std::find(begin(fd), end(fd), static_cast<TSpec>(0))!= end(fd))
	{
		// add pseudocounts
		double epsilon = 0.1;
		seqan::Pseudocount<TValue, CMode> p(epsilon);
		addValue(fd, p.pseudocount);
	}
	normalize(fd);
}

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
/// .Function.begin.param.object.type:Class.FrequencyDistribution

template<typename TValue, typename TSpec>
inline typename Iterator< FrequencyDistribution<TValue, TSpec>, Standard >::Type
begin(FrequencyDistribution<TValue, TSpec> & me, Standard)
{
	return me.frequency_list;
}
template<typename TValue, typename TSpec>
inline typename Iterator< FrequencyDistribution<TValue, TSpec> const, Standard >::Type
begin(FrequencyDistribution<TValue, TSpec> const & me, Standard)
{
	return me.frequency_list;
}

/*
template<typename TValue, typename TSpec>
inline typename Iterator< FrequencyDistribution<TValue, TSpec>, Standard >::Type
begin(FrequencyDistribution<TValue, TSpec> & me)
{
	return me.frequency_list;
}
template<typename TValue, typename TSpec>
inline typename Iterator< FrequencyDistribution<TValue, TSpec> const, Standard >::Type
begin(FrequencyDistribution<TValue, TSpec> const & me)
{
	return me.frequency_list;
}
*/

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
/// .Function.end.param.object.type:Class.FrequencyDistribution

template<typename TValue, typename TSpec>
void 
clear(FrequencyDistribution<TValue, TSpec> & fd) 
{
	std::fill(begin(fd), end(fd), static_cast<TSpec>(0));
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.convertResidueToFrequencyDist:
..summary:Coverts a residue to a frequency distribution (profile).
..cat:Motif Finding
..signature:convertResidueToFrequencyDist(frequencies,residue)
..param.frequencies:The @Class.FrequencyDistribution@ object representing the profile for a specific residue.
...type:Class.FrequencyDistribution
..param.residue:The residue object which is considered.
..remarks:This function is used to convert a sequence pattern into a profile (see Function.convertPatternToProfile).
*/

template<typename TValue, typename TSpec>
void 
convertResidueToFrequencyDist(FrequencyDistribution<TValue, TSpec> & fd, TValue const & residue)
{
	typedef typename Position< FrequencyDistribution<TValue, TSpec> >::Type TPos;
	TSpec probability = 
		static_cast<TSpec>(0.5/(ValueSize<TValue>::VALUE-1));

	for(TPos i=0; i<length(fd); ++i)
	{
		if(i==residue)
		{
			fd[i] = 0.5;
		}
		else
		{
			fd[i] = probability;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
/// .Function.end.param.object.type:Class.FrequencyDistribution

template<typename TValue, typename TSpec>
inline typename Iterator< FrequencyDistribution<TValue, TSpec>, Standard >::Type
end(FrequencyDistribution<TValue, TSpec> & me)
{
	return begin(me)+length(me);
}
template<typename TValue, typename TSpec>
inline typename Iterator< FrequencyDistribution<TValue, TSpec> const, Standard >::Type
end(FrequencyDistribution<TValue, TSpec> const & me)
{
	return begin(me)+length(me);
}

//////////////////////////////////////////////////////////////////////////////
/// .Function.length.param.object.type:Class.FrequencyDistribution

template<typename TValue, typename TSpec>
inline typename Size< FrequencyDistribution<TValue, TSpec> >::Type
length(FrequencyDistribution<TValue, TSpec> & me)
{
	return sizeof(me.frequency_list)/sizeof(me.frequency_list[0]);
}

template<typename TValue, typename TSpec>
inline typename Size< FrequencyDistribution<TValue, TSpec> >::Type
length(FrequencyDistribution<TValue, TSpec> const & me)
{
	return sizeof(me.frequency_list)/sizeof(me.frequency_list[0]);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.logarithmize:
..summary:Logarithmizes each element of a given @Class.FrequencyDistribution@ object.
..cat:Motif Finding
..signature:logarithmize(frequencies)
..param.frequencies:The @Class.FrequencyDistribution@ object.
...type:Class.FrequencyDistribution
*/

template<typename TValue, typename TSpec>
void 
logarithmize(FrequencyDistribution<TValue, TSpec> & fd) 
{
	typedef FrequencyDistribution<TValue, TSpec> TFrequencyDistribution;
	typedef typename Position<TFrequencyDistribution>::Type TPos;
	
	for(TPos i=0; i<length(fd); ++i)
	{
		fd[i] = static_cast<TSpec>(log(fd[i]));
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.normalize:
..summary:Determines the normalized frequencies of a given @Class.FrequencyDistribution@ object.
..cat:Motif Finding
..signature:normalize(frequencies)
..param.frequencies:The @Class.FrequencyDistribution@ object.
...type:Class.FrequencyDistribution
*/

template<typename TValue, typename TSpec>
void 
normalize(FrequencyDistribution<TValue, TSpec> & fd) //2)
{
	typedef FrequencyDistribution<TValue, TSpec> TFrequencyDistribution;
	typedef typename Position<TFrequencyDistribution>::Type TPos;
	
	TSpec amount = sum(fd);
	for(TPos i=0; i<length(fd); ++i)
	{
		fd[i] = static_cast<TSpec>(fd[i]/amount);
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.posOfMax:
..summary:Determines the residue position in a given @Class.FrequencyDistribution@ object with the maximum frequency.
..cat:Motif Finding
..signature:posOfMax(frequencies)
..param.frequencies:The @Class.FrequencyDistribution@ object.
...type:Class.FrequencyDistribution
*/

template<typename TValue, typename TSpec>
typename Position< FrequencyDistribution<TValue, TSpec> >::Type
posOfMax(FrequencyDistribution<TValue, TSpec> & me)
{
	typedef FrequencyDistribution<TValue, TSpec> TFrequencyDistribution; 
	typedef typename Position<TFrequencyDistribution>::Type TPos;

	TPos position = 0;
	TSpec max_value = static_cast<TSpec>(0);
	for(TPos i=0; i<length(me); ++i)
	{
		if(me[i]>max_value)
		{
			position = i;
			max_value = me[i];
		}
	}

	return position;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.setValue:
..summary:Sets a given value to the element at the specified position of a @Class.FrequencyDistribution@ object.
..cat:Motif Finding
..signature:setValue(frequencies,position,value)
..param.frequencies:The @Class.FrequencyDistribution@ object.
...type:Class.FrequencyDistribution
..param.position:Position of a given @Class.FrequencyDistribution@ object which is considered.
...type:Metafunction.Position
..param.value:The value object which represents a frequency value.
...remarks:The value object should be identical in type to the elements of the @Class.FrequencyDistribution@ object.
*/

template<typename TValue, typename TSpec, typename TType>
inline void
setValue(FrequencyDistribution<TValue, TSpec> & me, 
		 typename Position< FrequencyDistribution<TValue, TSpec> >::Type pos,
		 TType prob)
{
	me[pos] = static_cast<TSpec>(prob);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.sum:
..summary:Determines the sum of all frequencies in a given @Class.FrequencyDistribution@ object.
..cat:Motif Finding
..signature:sum(frequencies)
..param.frequencies:The @Class.FrequencyDistribution@ object.
...type:Class.FrequencyDistribution
*/

template<typename TValue, typename TSpec>
TSpec 
sum(FrequencyDistribution<TValue, TSpec> & me)
{
	TSpec amount = 
		std::accumulate(begin(me), end(me), static_cast<TSpec>(0));

	return amount;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.value:
..summary:Returns the frequency value at a specific position of a @Class.FrequencyDistribution@ object.
..cat:Motif Finding
..signature:value(frequencies,position)
..param.frequencies:The @Class.FrequencyDistribution@ object.
...type:Class.FrequencyDistribution
..param.position:Position of the specific element being searched in a given @Class.FrequencyDistribution@ object.
*/

template<typename TValue, typename TSpec, typename TPos>
inline TSpec &
value(FrequencyDistribution<TValue, TSpec> & me, TPos pos)
{
	return me[pos];
}
template<typename TValue, typename TSpec, typename TPos>
inline TSpec const &
value(FrequencyDistribution<TValue, TSpec> const & me, TPos pos)
{
	return me[pos];
}
								
//////////////////////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...