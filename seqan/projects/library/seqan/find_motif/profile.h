#ifndef SEQAN_HEADER_PROFILE_H
#define SEQAN_HEADER_PROFILE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
/**
.Shortcut.Profile:
..cat:Strings
..summary:A string of @Class.FrequencyDistribution@.
..signature:Profile
..shortcutfor:Spec.Alloc String
...signature:String<FrequencyDistribution, Alloc<> >
*/

//typedef String<FrequencyDistribution, Alloc<void> > Profile;

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

/**
.Function.convertPatternToProfile:
..summary:Converts a pattern into a profile which consists of a set of frequency distributions.
..cat:Motif Finding
..signature:convertPatternToProfile(profile,begin,end)
..param.profile:The  @Shortcut.Profile@ object which is a set of @Class.FrequencyDistribution|frequency distributions@.
...type:Shortcut.Profile
..param.begin:An iterator pointing to the beginning of a given sequence pattern which is either
              a @Shortcut.DnaString@ or a @Shortcut.Peptide@.
...type:Concept.Iterator
...type:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
..param.end:An iterator pointing to the end of a given sequence pattern which is either
            a @Shortcut.DnaString@ or a @Shortcut.Peptide@.
...type:Concept.Iterator
...type:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
..remarks:The number of @Class.FrequencyDistribution@ objects which together form the @Shortcut.Profile@ 
          equals the length of the given sequence.
..remarks:e.g.:$profile[0]$ represents the frequency distribution for the first residue of
          the given sequence.
..see:Function.convertResidueToFrequencyDist
*/

template<typename TProfile, typename TIterStr>
void 
convertPatternToProfile(TProfile & profile,
						TIterStr str_begin,
						TIterStr str_end)
{
	typedef typename Position<TProfile>::Type TPos;
	unsigned int str_size = str_end-str_begin;
	resize(profile, str_size);
	TPos pos = 0;
	while(str_begin!=str_end)
	{
		convertResidueToFrequencyDist(profile[pos], *str_begin);
		++str_begin;
		++pos;
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.convertSetOfPatternsToProfile:
..summary:Converts a set of sequence patterns into a profile.
..cat:Motif Finding
..signature:convertSetOfPatternsToProfile(profile,l_mers,pseudocount_mode)
..param.profile:The  @Shortcut.Profile@ object which is a set of @Class.FrequencyDistribution|frequency distributions@.
...type:Shortcut.Profile
..param.l_mers:The set of sequence patterns.
...type:Class.StringSet
..param.pseudocount_mode:The @Class.Pseudocount@ object for determining the pseudocount method.
...type:Class.Pseudocount
..remarks:This function is used, for example, in the refinement step of the PROJECTION algorithm to convert
          the collection of l-mers inside the corresponding buckets into a profile. 
*/

template<typename TProfile, typename TStrings, typename TPseudocountMode>
void
convertSetOfPatternsToProfile(TProfile & profile,
					 TStrings & l_mers, 
					 TPseudocountMode & pseudocount)
{
	typedef typename Value<TStrings>::Type TString;
	typedef typename Value<TProfile>::Type TFrequencyDistribution;
	typedef typename Position<TString>::Type TPos;

	typename Size<TString>::Type l = length(l_mers[0]);
	resize(profile, l);
	for(TPos i=0; i<l; ++i)
	{
		TFrequencyDistribution fd;
		typename Iterator<TStrings>::Type l_mers_iter = begin(l_mers);
		typename Iterator<TStrings>::Type l_mers_end = end(l_mers);
		while(l_mers_iter!=l_mers_end)
		{
			++fd[(int) *(begin(*l_mers_iter)+i)];
			++l_mers_iter;
		}
		profile[i] = fd;
	}

	normalize(profile, pseudocount);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.normalize:
..summary:Determines the normalized frequencies.
..cat:Motif Finding
..signature:normalize(container)
..signature:normalize(profile,pseudocount_mode)
..param.container:The @Class.FrequencyDistribution@ or @Shortcut.Profile@ object.
...type:Class.FrequencyDistribution
...type:Shortcut.Profile
..param.profile:The @Shortcut.Profile@ object which is a set of @Class.FrequencyDistribution|frequency distributions@.
...type:Shortcut.Profile
..param.pseudocount_mode:The @Class.Pseudocount@ object for determining the pseudocount method.
...type:Class.Pseudocount
..remarks:If necessary, pseudocounts are first added to the frequency values before normalizing them 
          when the parameter is of type @Shortcut.Profile@.
*/

template<typename TProfile>
void 
normalize(TProfile & profile)
{
	typename Iterator<TProfile>::Type iter = begin(profile);
	typename Iterator<TProfile>::Type iter_end = end(profile);
	while(iter!=iter_end)
	{
		normalize(*iter);
		++iter;
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.completeProfile:
..summary:Concatenates the background frequency with the profile for the motif component.
..cat:Motif Finding
..signature:completeProfile(profile,background_distribution)
..param.profile:The  @Shortcut.Profile@ object which is a set of @Class.FrequencyDistribution|frequency distributions@.
...type:Shortcut.Profile
..param.background_distribution:The @Class.FrequencyDistribution@ object which represents the backround distribution.
...type:Class.FrequencyDistribution
..remarks:The first row of the final @Shortcut.Profile|probability matrix@ represents the @Class.FrequencyDistribution|background distribution@.
*/

template<typename TProfile>
void 
completeProfile(TProfile & profile,
				typename Value<TProfile>::Type & background_distribution)
{
	TProfile copy(profile);
	resize(profile, length(copy)+1);
	profile[0] = background_distribution;

	typename Iterator<TProfile>::Type iter = begin(copy);
	int counter = 1;
	for(; !atEnd(iter, copy); goNext(iter))
	{
		profile[counter] = *iter;
		++counter;
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.display:
..summary:Displays a given set of strings.
..cat:Motif Finding
..signature:display(strings)
..param.strings:The set of strings.
...type:Class.StringSet
...type:Shortcut.Profile
..remarks:This function can also be used to display a @Shortcut.Profile|probability matrix@ 
          which is a set of @Class.FrequencyDistribution|frequency distributions@.
*/


template<typename TStrings>
void 
display(TStrings & strings)
{
	if(length(strings)!=0)
	{
		typename Iterator<TStrings>::Type iter = begin(strings);
		int counter = 0;
		for(; !atEnd(iter, strings); goNext(iter))
		{
			std::cout << "[" << counter << "]: " << *iter << "\n";
			++counter;
		}
		std::cout << "\n";
	}
	else
	{
		std::cout << "EMPTY STRINGS !!!\n";
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////

} // namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...