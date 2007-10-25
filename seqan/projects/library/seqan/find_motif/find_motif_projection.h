#ifndef SEQAN_HEADER_FIND_MOTIF_PROJECTION_H
#define SEQAN_HEADER_FIND_MOTIF_PROJECTION_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Projection
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Projection:
..summary: Represents the PROJECTION algorithm.
..general:Class.MotifFinder
..cat:Motif Finding
..signature:MotifFinder<TValue, Projection>
..param.TValue:The type of sequences to be analyzed.
...type:Spec.Dna
...type:Spec.AminoAcid
*/

///.Class.MotifFinder.param.TSpec.type:Spec.Projection

struct _Projection;
typedef Tag<_Projection> Projection;

//////////////////////////////////////////////////////////////////////////////
// MotifFinder - Projection Spec
//
// t:=dataset size (number of sequences)
// l:=motif size
// m:=number of possible l-mers
// d:=number of substitutions
// k:=projection size
// s:=bucket size
// tr:=number of independent trials
//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
class MotifFinder<TValue, Projection>
{
//____________________________________________________________________________________________

	enum { ALPHABET_SIZE = ValueSize<TValue>::VALUE };
	typedef String<TValue> TString;
	typedef String<TString> TStrings;
	typedef typename Size<TStrings>::Type TSize1;
	typedef typename Size<TString>::Type TSize2;

//____________________________________________________________________________________________

public:
	TSize1 dataset_size;
	TSize2 motif_size;
	unsigned int total_num_of_l_mers;
	unsigned int num_of_substitutions;
	unsigned int projection_size;
	unsigned int bucket_threshold;
	unsigned int num_of_trials;
	TString consensus_pattern;

//____________________________________________________________________________________________

	MotifFinder():
		dataset_size(0),
		motif_size(0),
		total_num_of_l_mers(0),
		num_of_substitutions(0),
		projection_size(0),
		bucket_threshold(0),
		num_of_trials(0)
	{
	}
	MotifFinder(TSize1 const & t_, 
				TSize2 const & l_, 
				unsigned int const & m_total_,
				unsigned int const & d_,
				unsigned int const & k_,
				unsigned int const & s_,
				unsigned int const & tr_):
		dataset_size(t_),
		motif_size(l_),
        total_num_of_l_mers(m_total_),
		num_of_substitutions(d_),
		projection_size(k_),
		bucket_threshold(s_),
		num_of_trials(tr_)
	{
	}
	MotifFinder(TSize1 const & t_, 
				TSize2 const & l_, 
				unsigned int const & m_total_,
				unsigned int const & d_):
		dataset_size(t_),
		motif_size(l_),
		total_num_of_l_mers(m_total_),
		num_of_substitutions(d_),
		projection_size(0),
		bucket_threshold(0),
		num_of_trials(0)
	{
		projection_size = 
			_computeProjectionSize<unsigned int>(ALPHABET_SIZE,
			                                     motif_size,
								                 num_of_substitutions,
								                 total_num_of_l_mers);
		
		bucket_threshold = 
			_computeBucketThreshold<unsigned int>(ALPHABET_SIZE,
			                                      motif_size,
								                  num_of_substitutions,
								                  total_num_of_l_mers,
												  projection_size);

		double prob_q = static_cast<double>(0.95);
		num_of_trials = _computeNumOfTrials(dataset_size,
								            motif_size,
								            num_of_substitutions,
								            projection_size,
								            bucket_threshold,
								            prob_q);
	}
	MotifFinder(MotifFinder const & other_):
		dataset_size(other_.dataset_size),
		motif_size(other_.motif_size),
		total_num_of_l_mers(other_.total_num_of_l_mers),
		num_of_substitutions(other_.num_of_substitutions),
		projection_size(other_.projection_size),
		bucket_threshold(other_.bucket_threshold),
		num_of_trials(other_.num_of_trials)
	{
	}
	~MotifFinder()
	{
	}

	MotifFinder const &
	operator = (MotifFinder const & other_)
	{
		if(this!=&other_)
		{
			this->dataset_size = other_.dataset_size;
			this->motif_size = other_.motif_size;
			this->total_num_of_l_mers = other_.total_num_of_l_mers;
			this->num_of_substitutions = other_.number_of_substitutions;
			this->projection_size = other_.projection_size;
			this->bucket_threshold = other_.bucket_threshold;
			this->num_of_trials = other_.num_of_trials;
		}

		return *this;
	}

//____________________________________________________________________________________________

}; // class MotifFinder<TValue, Projection>


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

/*
.Function._computeProjectionSize:
..summary:Computes the projection size (k).
..cat:Motif Finding
..signature:_computeProjectionSize(alp_size,l,d,m)
..param.alp_size:The size of the sequence alphabet.
...remarks:The alp_size object is four for nucleotide sequences and twenty for amino acid sequences.
..param.l:The size of the motif.
..param.d:The number of substitutions.
..param.m:The total number of possible l-mers of a given dataset.
*/

template<typename TType> 
TType
_computeProjectionSize(TType const & alp_size,
					   TType const & l, 
					   TType const & d,
					   TType const & m)
{
	TType result = static_cast<TType>(0);
	double numerator = log(static_cast<double>(m));
	double denominator = log(static_cast<double>(alp_size));

	//if((numerator/denominator)<=l-d+1)
	if(m<=exp(log(static_cast<double>(alp_size))*(l-d+1)) || (numerator/denominator)>l-d+1)
	{
		result = l-d-1;
	}
	else
	{
		result = floor(numerator/denominator)+1;
	}

	return result;
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._computeBucketThreshold:
..summary:Computes the bucket threshold size (s).
..cat:Motif Finding
..signature:_computeBucketThreshold(alp_size,l,d,m,k)
..param.alp_size:The size of the sequence alphabet.
...remarks:The alp_size object is four for nucleotide sequences and twenty for amino acid sequences.
..param.l:The size of the motif.
..param.d:The number of substitutions.
..param.m:The total number of possible l-mers of a given dataset.
..param.k:The projection size.
*/

template<typename TType> 
TType
_computeBucketThreshold(TType const & alp_size,
					    TType const & l, 
					    TType const & d,
					    TType const & m, 
						TType const & k)
{
	TType result = 3; // or 4

	if(m<=exp(log(static_cast<double>(alp_size))*(l-d+1)))
	{
		result = 
			ceil(static_cast<double>(m/pow(static_cast<double>(alp_size), static_cast<double>(k))*2));
	}

	return result;
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._computeNumOfTrials:
..summary:Computes the number of independent trials (tr).
..cat:Motif Finding
..signature:_computeNumOfTrials(t,l,d,k,s,prob_q)
..param.t:The number of input sequences.
..param.l:The size of the motif.
..param.d:The number of substitutions.
..param.k:The projection size.
..param.s:The bucket threshold size.
..param.prob_q:
...remarks:The prob_q object represents the probability that the planted bucket contains s or 
           more planted motif instances in at least one of the tr trials. Normally, we use 
		   prob_q=0.95.	
...type:$double$
..remarks:tr>= log(1-q)/log(B), where p is the probability that each motif occurence hashes 
          to the planted bucket and B is the probability that fewer than s planted occurences hash
          to the planted buckes in a given trial
*/

template<typename TType> 
TType
_computeNumOfTrials(TType const & t,
					TType const & l,
					TType const & d, 
					TType const & k, 
					TType const & s,
					double const & prob_q)
{
	double prob_p =
		static_cast<double>(binomialCoefficient( (l-d), k ))
	   /static_cast<double>(binomialCoefficient(l,k));
	
	double prob_B = static_cast<double>(0);
	for(unsigned int i=0; i<s; ++i)
	{
		prob_B+=
			static_cast<double>(binomialCoefficient(t,i))
		   *pow(prob_p, static_cast<double>(i))
		   *pow(static_cast<double>(1)-prob_p, static_cast<double>(t-i));
	}

	double numerator = log(static_cast<double>(1)-prob_q);
	double denominator = log(static_cast<double>(prob_B));
	TType result = 
		static_cast<TType>(ceil(static_cast<double>(numerator/denominator)));

	if(result<1)
	{
		result = 1;
	}
	
	return result;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStrings, typename TModel>
inline void
findMotif(MotifFinder<typename Value< typename Value<TStrings>::Type >::Type, Projection> & finder, 
		  TStrings & dataset,
		  TModel & seq_model)
{
	projection(finder.consensus_pattern,
		       dataset,
			   finder.motif_size,
			   finder.projection_size, 
			   finder.bucket_threshold, 
			   finder.num_of_trials,
			   seq_model);
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function.projection:
..summary:Represents the PROJECTION algorithm.
..cat:Motif Finding
..signature:projection(result,dataset,l,k,s,tr,model_type)
..param.result:The result object.
...remarks:The result object holds the consensus pattern of the found motif candidate.
...type:Class.String
..param.dataset:The dataset object representing the input sequences.
...type:Class.String
...signature:String<TString>
...param.TString:A @Class.String@ type
....type:Class.String
..param.l:The size of the motif.
..param.k:The projection size.
..param.s:The bucket threshold size.
..param.tr:The number of independent trials.
..param.model_type:The model_type object.
...type:Tag.OOPS
...type:Tag.ZOOPS
...type:Tag.TCM
..remarks:The PROJECTION algorithm which consists of two steps, the filtering and the refinement step,
          is able to run in OOPS, ZOOPS and TCM mode.
..remarks:The algorithm uses the EM procedure during the refinement phase which was introduced by Bailey and Elkan.
*/

template<typename TStrings, typename TType, typename TModel>
void 
projection(typename Value<TStrings>::Type & result,
		   TStrings & dataset, 
		   typename Size<typename Value<TStrings>::Type>::Type const & l,
		   TType const & k, 
		   TType const & s, 
		   TType const & tr,
		   TModel const & seq_model)
{
	typedef typename Value<TStrings>::Type TString;
	typedef typename Value<TString>::Type TValue;
	typedef typename Position<TString>::Type TPos;
	typedef String<int> TArray; //TBucket = TArray
	typedef String<TArray> TBuckets;

	// dataset information
	typename Size<TStrings>::Type t = length(dataset);

	// count_ar:=array of votes for each h(k-mer)
	typename Size<TArray>::Type ar_size = 
		pow(static_cast<double>(ValueSize<TValue>::VALUE), static_cast<int>(k));

	// array of collection of l-mers
	TBuckets bucket_ar;
	
	double maximum_score = static_cast<double>(-1)*DBL_MAX;
	for(unsigned int trial=0; trial<tr; ++trial)
	{
		TArray count_ar;
		resize(count_ar, ar_size);

		// ----------------------------------------------------------------------------
		// STEP 1:
		// filtering phase (:=key random projection phase)
		// ----------------------------------------------------------------------------
		// choose randomly k different positions
		std::set<int> positions;
		choosePositions(positions,l,k);

		// create shape
		String<char> shape_str;
		createShapeStr(shape_str, positions, l);
		//String<char> shape_str = "x_x";
		Shape<TValue, GappedShape3> shape;
		stringToShape(shape_str, shape);

		unsigned int num_of_relevant_buckets = 0;
		_filteringStep(bucket_ar,count_ar,num_of_relevant_buckets,dataset,shape,l,s);

		// ----------------------------------------------------------------------------
		// STEP 2:
		// checking phase (:= local search-based refinement procedure)
		// ----------------------------------------------------------------------------
		TPos i = 0;
		TPos j = 0;
		while(j<num_of_relevant_buckets & i<ar_size)
		{
			typename Size<TArray>::Type bucket_size = length(bucket_ar[i]);
			if(bucket_size>=s)
			{
				TStrings l_mers;
				resize(l_mers, bucket_size);
				typename Iterator<TArray>::Type bucket_iter = begin(bucket_ar[i]);
				int bucket_element = 0;
				for(; !atEnd(bucket_iter, bucket_ar[i]); goNext(bucket_iter))
				{
					TString l_mer = 
						inverseHash<TValue>(*bucket_iter, ValueSize<TValue>::VALUE, l);
					l_mers[bucket_element] = l_mer;
					++bucket_element;
				}

				TString consensus_pat;
				double score = _refinementStep(consensus_pat, l_mers, dataset, t, l, seq_model);
				if(score>maximum_score)
				{
					maximum_score = score;
					result = consensus_pat;
				}
				++j;
			}
			++i;
		}
		clear(bucket_ar);
	} 
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._filteringStep:
..summary:Given a position set with k different positions we compute a projection value
          for each l-mer in the input sequences and store the specific l-mer in the appropriate
		  bucket which is labeled with the specific projection value.
..cat:Motif Finding
..signature:_filteringStep(buckets,count_ar,num_of_relevant_buckets,dataset,shape,l,s)
..param.buckets:The projection buckets.
...type:Class.String
...signature:String< String<int> >
...remarks:For each hashed l-mers we store their corresponding integer value in the appropriate hash bucket.
..param.count_ar:The integer array which holds the number of entries (l-mers) contained in each projection bucket. 
...type:Class.String
...signature;String<int>
..param.num_of_relevant_buckets:The number of buckets having at least s entries.
..param.dataset:The dataset object representing the input sequences.
...type:Class.String
...signature:String<TString>
...param.TString:A @Class.String@ type
....type:Class.String
..param.shape:The @Class.Shape@ object.
...type:Class.Shape
...signature:Shape<TValue, GappedShape3>
...remarks:Is used for projecting each l-mer.
..param.l:The size of the motif.
..param.s:The bucket threshold size.
*/

template<typename TBucketAr, typename TArray, typename TType, typename TStrings, typename TShape>
void 
_filteringStep(TBucketAr & buckets, 
			   TArray & count_ar,
			   TType & num_of_relevant_buckets,
			   TStrings & dataset,
			   TShape & shape,
			   typename Size< typename Value<TStrings>::Type >::Type const & l,
			   TType const & s)
{
	typedef typename Value<TStrings>::Type TString;
	typedef typename Value<TString>::Type TValue;
	typedef typename Position<TString>::Type TPos;
	typedef typename Value<TBucketAr>::Type TBucket;
	typename Iterator<TStrings>::Type ds_iter = begin(dataset);
	typename Size<TArray>::Type ar_size = length(count_ar);
	Shape<TValue, SimpleShape> shape2(l, ValueSize<TValue>::VALUE); //to compute hash value of l-mer x
 
	// initialize pointer by setting it to null 
	// (=std::fill(begin(count_ar),end(count_ar),0))
	typename Iterator<TArray>::Type count_ar_iter = begin(count_ar) ;
	typename Iterator<TArray>::Type count_ar_end = end(count_ar);
	while(count_ar_iter!=count_ar_end)
	{
		*count_ar_iter = 0;
		++count_ar_iter;
	}

	// go over input sequences & increment corresponding counter in count_ar
	// fill l_mer_index with entries
	int y = 0; //hash-value of created k-mer
	int x = 0; //hash-value of l-mer
	for(; !atEnd(ds_iter, dataset); goNext(ds_iter))
	{
		typename Size<TString>::Type seq_length = length(*ds_iter);
		typename Iterator<TString>::Type seq_iter = begin(*ds_iter);
		typename Iterator<TString>::Type seq_end = begin(*ds_iter)+(seq_length-l+1);
		while( seq_iter!=seq_end )
		{
			y = hash(shape, seq_iter);
	    	++count_ar[y];
			++seq_iter;
		}
	}

	num_of_relevant_buckets = 
		std::count_if(begin(count_ar),
		              end(count_ar),
					  bind2nd(std::greater_equal<int>(),static_cast<int>(s)));
	
	resize(buckets, ar_size);
	TPos i = 0; //for count_ar
	TPos j = 0; //for buckets 
	while(j<num_of_relevant_buckets & i<ar_size)
	{
		if(count_ar[i]>=s)
		{
			TBucket bucket;
			resize(bucket, count_ar[i]);
			buckets[i] = bucket;
			++j;
		}
		++i;
	}

	// go over input sequences a second time and 
	// set the right sequence positions to the corresponding pos[i]
	ds_iter = begin(dataset);
	y = 0;
	for(; !atEnd(ds_iter, dataset); goNext(ds_iter))
	{
		typename Size<TString>::Type seq_length = length(*ds_iter);
		typename Iterator<TString>::Type seq_iter = begin(*ds_iter);
		typename Iterator<TString>::Type seq_end = begin(*ds_iter)+(seq_length-l+1);
		while(seq_iter!=seq_end)
		{
			y = hash(shape, seq_iter);
			typename Size<TBucket>::Type bucket_size = length(buckets[y]);
			if(bucket_size>=s)
			{
				x = hash(shape2, seq_iter);
				*(begin(buckets[y])+(bucket_size-count_ar[y])) = x;
				--count_ar[y];
			}
			++seq_iter;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._refinementStep:
..summary:Refines the collection of l-mers in each relevant bucket which contains at least s l-mers.
..cat:Motif Finding
..signature:_refinementStep(consensus_seq,l_mers,dataset,t,l,model_type)
..param.consensus_seq:The consensus pattern of the resulted set of l-mers after the refinement step.
..param.l_mers:The collection of l-mers inside a relevant bucket.
..param.dataset:The dataset object representing the input sequences.
...type:Class.String
...signature:String<TString>
...param.TString:A @Class.String@ type
....type:Class.String
..param.t:The number of input sequences.
..param.l:The size of the motif.
..param.model_type:The model_type object.
...type:Tag.OOPS
...type:Tag.ZOOPS
...type:Tag.TCM
*/

//////////////////////////////////////////////////////////////////////////////
//	OOPS model
//////////////////////////////////////////////////////////////////////////////

template<typename TStrings>
double 
_refinementStep(typename Value<TStrings>::Type & consensus_seq,
			     TStrings const & l_mers,
			     TStrings & dataset,
				 typename Size<TStrings>::Type const & t,
			     typename Size<typename Value<TStrings>::Type>::Type const & l,
				 OOPS const & oops)
{
	typedef typename Value<TStrings>::Type TString;
	typedef typename Value<TString>::Type TValue;

	// compute background frequency
	FrequencyDistribution<TValue> background;
	backgroundFrequency(background, begin(dataset), end(dataset));

	// step1: initial guess (profile) from bucket
	double epsilon = 0.1;
	Pseudocount<TValue, CMode> pseudocount_mode_c(epsilon);
	String< FrequencyDistribution<TValue> > profile;
	convertSetOfPatternsToProfile(profile, l_mers, pseudocount_mode_c);
	completeProfile(profile, background);

	// step2: refinement of initial profile with em: 5 trials
	double likelihood_score = 0;
	unsigned int num_of_trials = 5;
	unsigned int trial_nr = 0;
	while(trial_nr<num_of_trials)
	{
		//likelihood_score = em(profile, dataset, oops);
		likelihood_score = em(profile, begin(dataset), t, l, oops);
		++trial_nr;
	}

	// step3: form a guess for the planted motif by selecting from each input sequence
	//		  the l-mer x with the largest likelihood ratio
	TStrings set_of_l_mers;
	_getLMersWithTheLargestLikelihoodRatio(set_of_l_mers, begin(dataset), end(dataset), t, profile, l);

	// step4: form profile of l-mers inside the multiset &
	//		  score multiset	
	String< FrequencyDistribution<TValue> > new_profile;
	convertSetOfPatternsToProfile(new_profile, set_of_l_mers, pseudocount_mode_c);
	completeProfile(new_profile, background);

	// score each refined motif (multiset) by its likelihood ratio score
	double score = _computeLikelihoodRatioOfLMers(set_of_l_mers, new_profile);
	determineConsensusSeq(consensus_seq, new_profile, l);

	return score;
}

//////////////////////////////////////////////////////////////////////////////
//	ZOOPS model
//////////////////////////////////////////////////////////////////////////////

template<typename TStrings>
double 
_refinementStep(typename Value<TStrings>::Type & consensus_seq,
			     TStrings const & l_mers,
			     TStrings & dataset,
				 typename Size<TStrings>::Type const & t,
			     typename Size<typename Value<TStrings>::Type>::Type const & l,
				 ZOOPS const & zoops)
{
	typedef typename Value<TStrings>::Type TString;
	typedef typename Value<TString>::Type TValue;

	// compute background frequency
	FrequencyDistribution<TValue> background;
	backgroundFrequency(background, begin(dataset), end(dataset));

	// step1: initial guess (profile) from bucket
	double epsilon = 0.1;
	Pseudocount<TValue, CMode> pseudocount_mode_c(epsilon);
	String< FrequencyDistribution<TValue> > profile;
	convertSetOfPatternsToProfile(profile, l_mers, pseudocount_mode_c);
	completeProfile(profile, background);

	// step2: refinement of initial profile with em: 5 trials
	double gamma = static_cast<double>(1)/sqrt(static_cast<double>(t));
	double likelihood_score = 0;
	unsigned int num_of_trials = 5;//5;
	unsigned int trial_nr = 0;
	while(trial_nr<num_of_trials)
	{
		likelihood_score = em(profile, begin(dataset), t, l, gamma, zoops);
		++trial_nr;
	}

	// step3: form a guess for the planted motif by selecting from each input sequence
	//		  the l-mer x with the largest likelihood ratio
	TStrings set_of_l_mers;
	_getLMersWithTheLargestLikelihoodRatio(set_of_l_mers, begin(dataset), end(dataset), t, profile, l);

	// step4: form profile of l-mers inside the multiset &
	//		  score multiset	
	String< FrequencyDistribution<TValue> > new_profile;
	convertSetOfPatternsToProfile(new_profile, set_of_l_mers, pseudocount_mode_c);
	completeProfile(new_profile, background);

	// score each refined motif (multiset) by its likelihood ratio score
	double score = _computeLikelihoodRatioOfLMers(set_of_l_mers, new_profile);
	determineConsensusSeq(consensus_seq, new_profile, l);

	return score;
}

//////////////////////////////////////////////////////////////////////////////
//	TCM model
//////////////////////////////////////////////////////////////////////////////

template<typename TStrings, typename TType>
double 
_refinementStep(typename Value<TStrings>::Type & consensus_seq,
				TStrings const & l_mers,
				TStrings & dataset,
				TType const & t,
			    TType const & l,
				TCM const & tcm)
{
	typedef typename Value<TStrings>::Type TString;
	typedef typename Value<TString>::Type TValue;
	TType avg_m = 0;
	for(TType i=0; i<t; ++i)
	{
		avg_m += length(dataset[i])-l+1;
	}
	avg_m = avg_m/t;

	// compute background frequency
	FrequencyDistribution<TValue> background;
	backgroundFrequency(background, begin(dataset), end(dataset));

	// step1: initial guess (profile) from bucket
	double epsilon = 0.1;
	Pseudocount<TValue, CMode> pseudocount_mode_c(epsilon);
	String< FrequencyDistribution<TValue> > profile;
	convertSetOfPatternsToProfile(profile, l_mers, pseudocount_mode_c);
	completeProfile(profile, background);

	// step2: refinement of initial profile with em: 5 trials
	double lambda = 
		static_cast<double>(1)/(sqrt(static_cast<double>(t))*static_cast<double>(avg_m));
	double likelihood_score = 0;
	unsigned int num_of_trials = 5;
	unsigned int trial_nr = 0;
	while(trial_nr<num_of_trials)
	{
		likelihood_score = em(profile, begin(dataset), t, l, lambda, tcm);
		++trial_nr;
	}

	// step3: form a guess for the planted motif by selecting from each input sequence
	//		  the l-mer x with the largest likelihood ratio
	TStrings set_of_l_mers;
	_getLMersWithTheLargestLikelihoodRatio(set_of_l_mers,begin(dataset),end(dataset),t,profile,l);

	// step4: form profile of l-mers inside the multiset &
	//		  score multiset	
	String< FrequencyDistribution<TValue> > new_profile;
	convertSetOfPatternsToProfile(new_profile, set_of_l_mers, pseudocount_mode_c);
	completeProfile(new_profile, background);

	// score each refined motif (multiset) by its likelihood ratio score
	double score = _computeLikelihoodRatioOfLMers(set_of_l_mers, new_profile);
	determineConsensusSeq(consensus_seq, new_profile, l);

	return score;
}

//////////////////////////////////////////////////////////////////////////////
//Subfunctions
//////////////////////////////////////////////////////////////////////////////

/*
.Function.choosePositions:
..summary:Chooses randomly k different positions from {0,1,...,(l-1)}
..cat:Motif Finding
..signature:choosePositions(positions,l,k)
..param.positions:The set of k chosen positions.
...type:$set<int>$
..param.l:The size of the motif.
..param.k:The projection size.
*/

template<typename TAssociativeContainer, typename TType>
void
choosePositions(TAssociativeContainer & positions, TType const & l, TType const & k)
{
	//srand((unsigned) time(NULL));

	while(positions.size()<k)
	{
		int position = rand() % l;
		positions.insert(position);
	}
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function.createShapeStr:
..summary:Given a set of k positions createShapeStr creates a shape string of length l
          to mark the k chosen positions which are represented by 'x'.
..signature:createShapeStr(shape_str,positions,l)
..param.shape_str:The shape_str object.
...type:Class.String
...signature:String<char>
...remarks:The shape_str object uses two different characters, '-' and 'x' for marking the
           position which is chosen.
..param.positions:The set of k chosen positions.
...type:$set<int>$
..param.l:The size of the motif.
*/

template<typename TString, typename TAssociativeContainer>
void
createShapeStr(TString & shape_str, 
			   TAssociativeContainer & positions,
			   typename Size<TString>::Type const & l)
{
	resize(shape_str, l);
	std::fill(begin(shape_str), end(shape_str), '_');

	TAssociativeContainer::iterator iter_pos = positions.begin();
	TAssociativeContainer::iterator iter_end = positions.end();
	while(iter_pos!=iter_end)
	{
		shape_str[*iter_pos] = 'x';
		++iter_pos;
	}
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._getLMersWithTheLargestLikelihoodRatio:
..summary:Forms a guess for the planted motif by selecting from each input sequence
		  the l-mer x with the largest likelihood ratio.
..signature:_getLMersWithTheLargestLikelihoodRatio(l_mers,dataset_start,dataset_end,t,profile,l)
..param.l_mers:The collection of t l-mers.
..param.dataset_start:An iterator pointing to the first input sequence of a given dataset.
..param.dataset_end:An iterator pointing to the last input sequence of a given dataset.
..param.t:The number of input sequences.
..param.profile:The profile object which is a set of frequency distributions.
...type:Class.String
....signature:String<TFrequencyDistribution>
..param.l:The size of the motif.
*/

template<typename TStrings, typename TIter, typename TType, typename TProfile>
void
_getLMersWithTheLargestLikelihoodRatio(TStrings & l_mers,
									   TIter dataset_start,
									   TIter dataset_end,
									   TType const & t,
									   TProfile const & profile,
								       TType const & l)
{
	typedef typename Value<TStrings>::Type TString;
	resize(l_mers, t);
	TType pos = static_cast<TType>(0);
	while(dataset_start!=dataset_end)
	{
		TType m = length(*dataset_start)-l+1;
		double maximum =
			_computeLikelihoodRatioOfLMer(begin(*dataset_start), begin(*dataset_start)+l, profile);
		TString l_mer = seqan::infix(*dataset_start, 0, l);

		for(TType i=1; i<m; ++i)
		{
			double likelihood_ratio = 
				_computeLikelihoodRatioOfLMer(begin(*dataset_start)+i, begin(*dataset_start)+i+l, profile);
			if(likelihood_ratio>maximum)
			{
				l_mer = seqan::infix(*dataset_start, i, i+l);
				maximum = likelihood_ratio;
			}
		}
		l_mers[pos] = l_mer;
		++pos;
		++dataset_start;
	}
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._computeLikelihoodRatioOfLMer:
..summary:Computes the likelihood ratio of a given l-mer.
..signature:_computeLikelihoodRatioOfLMer(l_mer_begin,l_mer_end,profile)
..param.l_mer_begin:An iterator pointing to the beginning of a given l-mer pattern.
...type:Concept.Iterator Iterator
....remarks:Standard conform iterator
...type:Shortcut.DnaIterator
....remarks:Iterator for @Shortcut.DnaString@ (a string of @Spec.Dna@).
....see:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
....remarks:Iterator for @Shortcut.Peptide@ (a string of @Spec.AminoAcid@).
....see:Shortcut.PeptideIterator
..param.l_mer_end:An iterator pointing to the end of a given l-mer pattern.
...type:Concept.Iterator Iterator
....remarks:Standard conform iterator
...type:Shortcut.DnaIterator
....remarks:Iterator for @Shortcut.DnaString@ (a string of @Spec.Dna@).
....see:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
....remarks:Iterator for @Shortcut.Peptide@ (a string of @Spec.AminoAcid@).
....see:Shortcut.PeptideIterator
..param.profile:The profile object which is a set of frequency distributions.
...type:Class.String
....signature:String<TFrequencyDistribution>
..remarks:Computes the sum of log probabilites instead of the product of probabilites
*/

template<typename TStrIter, typename TProfile>
double
_computeLikelihoodRatioOfLMer(TStrIter l_mer_begin, 
							  TStrIter l_mer_end,
							  TProfile const & profile)
{
	double result = static_cast<double>(0);
	typedef typename Position<TProfile>::Type TPos;
	TProfile log_profile = profile;
	for(TPos i=0; i<length(log_profile); ++i)
	{
		logarithmize(log_profile[i]);
	}

	double motif_component = static_cast<double>(0);
	double backgr_component = static_cast<double>(0);
	unsigned int pos = 0;
	while(l_mer_begin!=l_mer_end)
	{
		motif_component += log_profile[pos+1][(int)*l_mer_begin];
		backgr_component += log_profile[0][(int)*l_mer_begin];
		++l_mer_begin;
		++pos;
	}
	result = motif_component-backgr_component;

	return exp(result);
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function._computeLikelihoodRatioOfLMers:
..summary:Computes the likelihood ratio of a given set of l-mers.
..signature:_computeLikelihoodRatioOfLMers(l_mers,profile)
..param.l_mers:The collection of l-mers.
..param.profile:The profile object which is a set of frequency distributions.
...type:Class.String
....signature:String<TFrequencyDistribution>
*/

template<typename TStrings, typename TProfile>
double
_computeLikelihoodRatioOfLMers(TStrings const & l_mers, 
							   TProfile const & profile)
{
	double score = static_cast<double>(1);
	typename Size<TStrings>::Type num_of_l_mers = length(l_mers);
	for(typename Position<TStrings>::Type i=0; i<num_of_l_mers; ++i)
	{
		score *= _computeLikelihoodRatioOfLMer(begin(l_mers[i]), end(l_mers[i]), profile);
	}

	return score;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.determineConsensusSeq:
..summary:Determines the consensus pattern of a given profile.
..cat:Motif Finding
..signature:determineConsensusSeq(consensus_seq,profile,l)
..param.consensus_seq:The consensus pattern.
...type:Class.String
...type:Shortcut.DnaString
...type:Shortcut.Peptide
..param.profile:The  @Shortcut.Profile@ object which is a set of @Class.FrequencyDistribution|frequency distributions@.
...type:Shortcut.Profile
..param.l:The size of the motif.
*/

template<typename TString, typename TProfile>
void
determineConsensusSeq(TString & consensus_seq,
					  TProfile & profile,
					  typename Size<TString>::Type const & l)
{
	typedef Value<TString>::Type TValue;
	typename Position<TString>::Type i;

	resize(consensus_seq, l);
	for(i=1; i<=l; ++i)
	{
		consensus_seq[i-1] = 
			static_cast<TValue>(posOfMax(profile[i]));
	}
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function.displayResult:
..summary:Displays the consensus pattern of the found motif candidate.
..cat:Motif Finding
..signature:displayResult(motif_finder)
..param.motif_finder:The @Class.MotifFinder@ object.
...type:Class.MotifFinder
..param.dataset:The dataset object representing the input sequences.
...type:Class.String
...signature:String<TString>
...param.TString:A @Class.String@ type
....type:Class.String
*/

template<typename TValue>
void
displayResult(MotifFinder<TValue, Projection> & projection)
{
	if(length(projection.consensus_pattern)!=0)
	{
		std::cout << projection.consensus_pattern << "\n";
	}
	else
	{
		std::cout << "NO MOTIF HAS BEEN FOUND!!!\n";
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...