#ifndef SEQAN_HEADER_FIND_MOTIF_EPatternBranching_H
#define SEQAN_HEADER_FIND_MOTIF_EPatternBranching_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// EPatternBranching
//////////////////////////////////////////////////////////////////////////////

/**
.Spec.EPatternBranching:
..summary: Represents the ePatternBranching algorithm.
..general:Class.MotifFinder
..cat:Motif Finding
..signature:MotifFinder<TValue, EPatternBranching>
..param.TValue:The type of sequences to be analyzed.
...type:Spec.Dna
...type:Spec.AminoAcid
*/

///.Class.MotifFinder.param.TSpec.type:Spec.EPatternBranching

struct _EPatternBranching;
typedef Tag<_EPatternBranching> EPatternBranching;

//////////////////////////////////////////////////////////////////////////////
// MotifFinder - EPatternBranching Spec
//
// t:=dataset size (number of sequences)
// n:=average sequence size
// l:=motif size
// d:=number of substitutions
// h:=size of the neighborhood considering at first
//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
class MotifFinder<TValue, EPatternBranching>
{
//_________________________________________________________________________________
	
	typedef String<TValue> TString;
	typedef String<TString> TStrings;
	typedef typename Size<TStrings>::Type TSize1;
	typedef typename Size<TString>::Type TSize2;
	typedef String<int> TIntAr;

//_________________________________________________________________________________

public:
	TSize1 dataset_size;
	TSize2 motif_size;
	unsigned int num_of_substitutions;
	bool has_exact_substitutions;
	unsigned int neighborhood_size;
	TStrings set_of_motifs; // result set

	MotifFinder()
	{
	}
	MotifFinder(TSize1 const & t_, 
		        TSize2 const & l_, 
				unsigned int const & d_, 
				bool const & is_exact_, 
				unsigned int const & h_):
		dataset_size(t_),
	    motif_size(l_),
		num_of_substitutions(d_),
		has_exact_substitutions(is_exact_),
		neighborhood_size(h_)
	{
	}
	MotifFinder(TSize1 const & t_, 
		        TSize2 const & l_, 
				unsigned int const & d_, 
				bool const & is_exact_, 
				TIntAr & n_ar_):
		dataset_size(t_),
	    motif_size(l_),
		num_of_substitutions(d_),
		has_exact_substitutions(is_exact_),
		neighborhood_size(0)
	{
		neighborhood_size = 
			computeH(dataset_size, motif_size, num_of_substitutions, is_exact_, n_ar_); 
	}
	MotifFinder(MotifFinder const & other_):
	    dataset_size(other_.dataset_size),
		motif_size(other_.motif_size),
		num_of_substitutions(other_.num_of_substitutions),
		has_exact_substitutions(other_.has_exact_substitutions)
		neighborhood_size(other_.neighborhood_size)
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
			dataset_size = other_.dataset_size;
			motif_size = other_.motif_size;
			num_of_substitutions = other_.num_of_substitutions;
			has_exact_substitutions = other_.has_exact_substitutions;
			neighborhood_size = other_.neighborhood_size;
		}

		return *this;
	}

//_________________________________________________________________________________

}; // class MotifFinder<TValue, EPatternBranching>

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

/*
.Function.computeH:
..summary:Computes the size of the neighborhood considering at first.
..cat:Motif Finding
..signature:computeH(t,l,d,is_exact,n_ar)
..param.t:The number of input sequences.
..param.l:The size of the motif.
..param.d:The number of substitutions.
..param.is_exact:The size of Hamming distance
...type:$bool$
..param.n_ar:The array with the sequence lengths of each input sequence.
...remarks:If all input sequences have the same sequence lengths n_ar consists of t equal values.
*/

template<typename TType, typename TIntAr>
TType
computeH(TType const & t, TType const & l, TType const & d, bool const & is_exact, TIntAr & n_ar) 
{
	TType d_bar = d-1;
	double sum = static_cast<double>(0); //probability p_d
	double prob_P = static_cast<double>(1);

	if(is_exact)
	{
		sum = static_cast<double>(binomialCoefficient(l,d_bar))
			*pow(0.75, static_cast<double>(d_bar))
			*pow(0.25, static_cast<double>(l-d_bar));
	}
	else
	{
		for(unsigned int i=0; i<=d_bar; ++i)
		{
			sum+= 
				static_cast<double>(binomialCoefficient(l,i))
			   *pow(0.75, static_cast<double>(i))
			   *pow(0.25, static_cast<double>(l-i));
		}
	}

	do{
		++d_bar;
		if(is_exact)
		{
			sum = static_cast<double>(0);
			sum = static_cast<double>(binomialCoefficient(l,d_bar))
				*pow(0.75, static_cast<double>(d_bar))
		        *pow(0.25, static_cast<double>(l-d_bar));
		}
		else
		{
			sum += static_cast<double>(binomialCoefficient(l,d_bar))
				  *pow(0.75, static_cast<double>(d_bar))
		          *pow(0.25, static_cast<double>(l-d_bar));
		}
		prob_P = static_cast<double>(1);
		for(unsigned int i=0; i<t; ++i)
		{
			TType n = n_ar[i];
			prob_P*=
				static_cast<double>(1)-pow(static_cast<double>(1)-sum, static_cast<double>(n-l+1));
		}
	} while( (prob_P<0.95) & (d_bar<2*d-1) );

	TType h = 2*d-d_bar-1;
	if(h>=d)
	{
		h=0;
	}
	
	return h;
}

/////////////////////////////////////////////////////////////////////////

template<typename TStrings, typename TModel>
inline void
findMotif(MotifFinder<typename Value<typename Value<TStrings>::Type>::Type,EPatternBranching> & epb2, 
		  TStrings & dataset, 
		  TModel & seq_model)
{
	ePatternBranching(epb2.set_of_motifs,
		               dataset,
					   epb2.motif_size,
					   epb2.num_of_substitutions,
					   epb2.has_exact_substitutions,
					   epb2.neighborhood_size,
					   seq_model);
}

//////////////////////////////////////////////////////////////////////////////

/*
.Function.ePatternBranching:
..summary:Represents the ePatternBranching algorithm.
..cat:Motif Finding
..signature:ePatternBranching(result_set,dataset,l,d,is_exact,h,seq_model)
..param.result_set:The result_set object.
..param.dataset:The dataset object representing the input sequences.
...type:Class.String
...signature:String<TString>
...param.TString:A @Class.String@ type
....type:Class.String
..param.l:The size of the motif.
..param.d:The number of substitutions.
..param.is_exact:The size of Hamming distance.
...type:$bool$
..param.h:The size of the neighborhood considering at first.
..param.seq_model:The seq_model object.
...type:Tag.OOPS
...type:Tag.OMOPS
..remarks:The ePatternBranching algorithm is able to run in OOPS and OMOPS mode.
*/

//////////////////////////////////////////////////////////////////////////////
//	original version = OMOPS model
//////////////////////////////////////////////////////////////////////////////

//considering only the first input sequence
template<typename TStrings, typename TType>
void 
ePatternBranching(TStrings & result_set,
				   TStrings & dataset,
				   TType const & l,
				   TType const & d,
				   bool const & is_exact,
				   TType & h,
				   OMOPS const & omops)
{
	typedef typename Value<TStrings>::Type TString;
	typedef typename Value<TString>::Type TValue;
	typedef typename String<int> TIntAr;
	typename Size<TStrings>::Type t = length(dataset);
	typename Iterator<TStrings>::Type ds_iter = begin(dataset);
	Shape<TValue, SimpleShape> shape(l, ValueSize<TValue>::VALUE);
	std::set<int> result;

	//we only consider l-mers from the first input sequence as starting points
	typename Size<TString>::Type seq_len = length(*ds_iter);
	typename Iterator<TString>::Type seq_iter = begin(*ds_iter);
	typename Iterator<TString>::Type seq_end = begin(*ds_iter)+(seq_len-l+1);
	while(seq_iter!=seq_end)
	{	
		//build set of h-neighborhood of l-mer x
		std::vector<int> V;
		createDVariants(V, seq_iter, l, h, is_exact, shape);
		std::vector<int>::iterator V_iter = V.begin();
		std::vector<int>::iterator V_end = V.end();
		while(V_iter!=V_end)
		{
			TIntAr Chi;
			resize(Chi, 1);
			Chi[0] = *V_iter;
			for(unsigned int j=h; j<=(d-1); ++j)
			{
				std::set<int> good_neighbors;
				typename Iterator<TIntAr>::Type Chi_iter = begin(Chi);
				typename Iterator<TIntAr>::Type Chi_end = end(Chi);
				while(Chi_iter!=Chi_end)
				{
					std::set<int> neighbors;
					bestNeighbors(neighbors,*Chi_iter,j,l,d,dataset);
					std::set_union(good_neighbors.begin(),good_neighbors.end(),
								   neighbors.begin(),neighbors.end(),
								   std::inserter(good_neighbors, good_neighbors.end()));
					++Chi_iter;
				}
				clear(Chi);
				resize(Chi, good_neighbors.size());
				std::copy(good_neighbors.begin(), good_neighbors.end(), begin(Chi));
				good_neighbors.clear();
			}
			// filtering relevant motif candidates
			// check if pattern occurs at least once in each sequence
			typename Iterator<TIntAr>::Type ar_iter = begin(Chi);
			typename Iterator<TIntAr>::Type ar_end = end(Chi);
			while(ar_iter!=ar_end)
			{
				TString pattern =
					inverseHash<TValue>(*ar_iter, ValueSize<TValue>::VALUE, l);
				ds_iter = begin(dataset);
				typename Iterator<TStrings>::Type ds_end = end(dataset);
				bool hasOMOPS = true;
				while((ds_iter!=ds_end) & hasOMOPS)
				{
					if(!(hasAtLeastOneOccurrence(begin(pattern),begin(*ds_iter),end(*ds_iter),l,d,is_exact)))
					{
						hasOMOPS = false;
					}
					++ds_iter;
				}
				if(hasOMOPS)
				{
					result.insert(*ar_iter);
				}
				++ar_iter;
			}
			++V_iter;
		}
		++seq_iter;
	}
	resize(result_set, result.size());
	std::set<int>::iterator set_iter = result.begin();
	std::set<int>::iterator set_end = result.end();
	typename Position<TStrings>::Type pos = 0;
	while(set_iter!=set_end)
	{
		TString l_mer = 
			inverseHash<TValue>(*set_iter, ValueSize<TValue>::VALUE, l);
		result_set[pos] = l_mer;
		++set_iter;
		++pos;
	}
}

//////////////////////////////////////////////////////////////////////////////
//	OOPS model
//////////////////////////////////////////////////////////////////////////////

//we only consider l-mers from the first input sequence as starting points
template<typename TStrings, typename TType>
void 
ePatternBranching(TStrings & result_set,
				  TStrings & dataset,
				  TType const & l,
				  TType const & d, 
				  bool const & is_exact,
				  TType & h,
				  OOPS const & oops)
{
	typedef typename Value<TStrings>::Type TString;
	typedef typename Value<TString>::Type TValue;
	typedef typename String<int> TIntAr;
	typename Size<TStrings>::Type t = length(dataset);
	typename Iterator<TStrings>::Type ds_iter = begin(dataset);
	Shape<TValue, SimpleShape> shape(l, ValueSize<TValue>::VALUE);
	std::set<int> result;

	//we only consider l-mers from the first input sequence as starting points
	typename Size<TString>::Type seq_len = length(*ds_iter);
	typename Iterator<TString>::Type seq_iter = begin(*ds_iter);
	typename Iterator<TString>::Type seq_end = begin(*ds_iter)+(seq_len-l+1);
	while(seq_iter!=seq_end)
	{
		//build set of h-neighborhood of l-mer x
		std::vector<int> V;
		createDVariants(V, seq_iter, l, h, is_exact, shape);
		std::vector<int>::iterator V_iter = V.begin();
		std::vector<int>::iterator V_end = V.end();
		while(V_iter!=V_end)
		{
			TIntAr Chi;
			resize(Chi, 1);
			Chi[0] = *V_iter;
			for(unsigned int j=h; j<=(d-1); ++j)
			{
				std::set<int> good_neighbors;
				typename Iterator<TIntAr>::Type Chi_iter = begin(Chi);
				typename Iterator<TIntAr>::Type Chi_end = end(Chi);
				while(Chi_iter!=Chi_end)
				{
					std::set<int> neighbors;
					bestNeighbors(neighbors,*Chi_iter,j,l,d,dataset);
					std::set_union(good_neighbors.begin(),good_neighbors.end(),
								   neighbors.begin(),neighbors.end(),
								   std::inserter(good_neighbors, good_neighbors.end()));
					++Chi_iter;
				}
				clear(Chi);
				resize(Chi, good_neighbors.size());
				std::copy(good_neighbors.begin(), good_neighbors.end(), begin(Chi));
				good_neighbors.clear();
			}
			// filtering relevant motif candidates
			typename Iterator<TIntAr>::Type ar_iter = begin(Chi);
			typename Iterator<TIntAr>::Type ar_end = end(Chi);
			while(ar_iter!=ar_end)
			{
				TString pattern =
					inverseHash<TValue>(*ar_iter, ValueSize<TValue>::VALUE, l);
				ds_iter = begin(dataset);
				typename Iterator<TStrings>::Type ds_end = end(dataset);
				bool hasExactOOPS = true;
				while((ds_iter!=ds_end) & hasExactOOPS)
				{
					if(!(hasExactOneOccurrence(begin(pattern),begin(*ds_iter),end(*ds_iter),l,d,is_exact)))
					{
						hasExactOOPS = false;
					}
					++ds_iter;
				}
				if(hasExactOOPS)
				{
					result.insert(*ar_iter);
				}
				++ar_iter;
			}
			++V_iter;
		}
		V.clear();
		++seq_iter;
	}
	resize(result_set, result.size());
	std::set<int>::iterator set_iter = result.begin();
	std::set<int>::iterator set_end = result.end();
	typename Position<TStrings>::Type pos = 0;
	while(set_iter!=set_end)
	{
		TString l_mer = 
			inverseHash<TValue>(*set_iter, ValueSize<TValue>::VALUE, l);
		result_set[pos] = l_mer;
		++set_iter;
		++pos;
	}
}

//////////////////////////////////////////////////////////////////////////////
//Subfunctions
//////////////////////////////////////////////////////////////////////////////

/*
.Function.bestNeighbors:
..summary:Represents the second version of GoodNeighbors described in the paper.
..cat:Motif Finding
..signature:bestNeighbors(neighbors,l_mer,j,l,d,dataset)
..param.neighbors:The set of best neighbors being in the Hamming distance 1-vicinity of l_mer.
..param.l_mer:The respective pattern represented by an integer value (hash value) 
              which is mapped to its best 1-neighborhood.
...type:Class.String
..param.j:The current iteration phase.
..param.l:The size of the motif.
..param.d:The number of substitutions.
..param.dataset:The dataset object representing the input sequences.
...type:Class.String
...signature:String<TString>
...param.TString:A @Class.String@ type
....type:Class.String
..remarks:The function uses equations formulated in theorem 1.
*/

template<typename TIntSet, typename TType, typename TStrings>
void 
bestNeighbors(TIntSet & neighbors, 
			   int & l_mer, 
			   TType const & j, 
			   TType const & l,
			   TType const & d,
			   TStrings & dataset) 
{
	typedef typename Value<TStrings>::Type TString;
	typedef typename Value<TString>::Type TValue;
	typedef String<int> TIntArray;
	typedef typename Position<TString>::Type TPos;
	Shape<TValue, SimpleShape> shape(l, ValueSize<TValue>::VALUE);
	typename Size<TStrings>::Type t = length(dataset);
	int seq_nr = 0;
	int seq_pos = 0;
	int beta = static_cast<int>((t-1)*(2*d-j-1));
	TString l_mer_x = 
		inverseHash<TValue>(l_mer,ValueSize<TValue>::VALUE,l);

	//create distance matrix for x & compute maximum distance of x from dataset X
	int max_dist_of_x = 0;
	int ** hd_mat = new int*[t];
	typename Iterator<TStrings>::Type ds_iter = begin(dataset);
	for(; !atEnd(ds_iter, dataset); goNext(ds_iter))
	{
		int min_dist = INT_MAX;
		typename Size<TString>::Type seq_len = length(*ds_iter);
		typename Iterator<TString>::Type seq_iter = begin(*ds_iter);
		typename Iterator<TString>::Type seq_end = begin(*ds_iter)+(seq_len-l+1);
		hd_mat[seq_nr] = new int[seq_len-l+1];
		seq_pos = 0;
		while(seq_iter!=seq_end)
		{
			typename Iterator<TString>::Type x_iter = begin(l_mer_x);
			typename Iterator<TString>::Type x_end = end(l_mer_x);
			int dist = hammingDistance(x_iter,x_end,seq_iter);
			hd_mat[seq_nr][seq_pos] = dist;
			if(dist<min_dist)
			{
				min_dist = dist;
			}
			++seq_iter;
			++seq_pos;
		}

		if(min_dist>max_dist_of_x)
		{
			max_dist_of_x = min_dist;
		}
		++seq_nr;
	}

	// compute delta = max_dist_of_x-d
	int uppercase_delta = max_dist_of_x-d;

	//build set of 1-neighborhood of l-mer x
	std::vector<int> V; 
	createDVariants(V, begin(l_mer_x), l, (unsigned int)1, false, shape);
	std::vector<int>::iterator V_iter = V.begin();
	std::vector<int>::iterator V_end = V.end();
	while(V_iter!=V_end)
	{
		if(l_mer==*V_iter)
		{
			neighbors.insert(*V_iter);
		}
		else
		{
			ds_iter = begin(dataset);
			++ds_iter; // begin with the second sequence
			TString l_mer_y = inverseHash<TValue>(*V_iter, ValueSize<TValue>::VALUE, l);
			int sum_of_dist_y = 0;
			int max_dist_of_y = 0;

			// find out the mismatch position of both l-mers x and y 
			TPos mismatch_pos = 0;
			TValue letter_in_x; //base in x at the mismatch position
			TValue letter_in_y; //base in y at the mismatch position
			for(TPos i=0; i<l; ++i)
			{
				if(l_mer_x[i]!=l_mer_y[i])
				{
					mismatch_pos = i;
					letter_in_x = l_mer_x[mismatch_pos];
					letter_in_y = l_mer_y[mismatch_pos];
					i = l;
				}
			}
			seq_nr = 0;
			for(; !atEnd(ds_iter, dataset); goNext(ds_iter))
			{
				typename Size<TString>::Type seq_len = length(*ds_iter);
				typename Iterator<TString>::Type seq_iter = begin(*ds_iter);
				int m = seq_len-l+1; //#l-mers
				typename Iterator<TString>::Type seq_end = begin(*ds_iter)+m;

				int dist_x = *std::min_element(&hd_mat[seq_nr][0],&hd_mat[seq_nr][m]);
				// num1:=#l-mers having distance dist_x to l_mer x
				unsigned int num1 = 
					std::count(&hd_mat[seq_nr][0], 
							   &hd_mat[seq_nr][m], dist_x);
				// num2:=#l-mers having distance (dist_x+1) to l_mer x
				unsigned int num2 = 
					std::count(&hd_mat[seq_nr][0], 
							   &hd_mat[seq_nr][m], 
							   dist_x+1);
				unsigned int counter1 = 0;
				unsigned int counter2 = 0;
				unsigned int counter3 = 0;
				seq_pos = 0;
				while(seq_iter!=seq_end & counter1<1)
				{
					if( (hd_mat[seq_nr][seq_pos]==dist_x) &
						(letter_in_y==*(seq_iter+mismatch_pos)) )
					{
						++counter1;
					}
					else if( (hd_mat[seq_nr][seq_pos]==dist_x) &
							 (letter_in_x==*(seq_iter+mismatch_pos)) )
					{
						++counter2;
					}
					else if( (hd_mat[seq_nr][seq_pos]==dist_x+1)& 
							 (letter_in_y!=*(seq_iter+mismatch_pos)) )
					{
						++counter3;
					}
					++seq_iter;
					++seq_pos;
				}
				int dist_y = dist_x;
				if(counter1>0)
				{
					dist_y-=1;
				}
				else if( (counter2==num1) & (counter3==num2) )
				{
					dist_y+=1;
				}

				sum_of_dist_y+=dist_y;
				if(dist_y>max_dist_of_y)
				{
					max_dist_of_y = dist_y;
				}
				++seq_nr;
			}
			int lowercase_delta = max_dist_of_y-max_dist_of_x;
			if( (uppercase_delta==d-j) & (lowercase_delta==-1) & (sum_of_dist_y<=beta) )
			{
				neighbors.insert(*V_iter);
			}
			else if( (uppercase_delta==d-j-1) & (lowercase_delta<=0) & (sum_of_dist_y<=beta) )
			{
				neighbors.insert(*V_iter);
			}
			else if( (uppercase_delta<d-j-1) & (sum_of_dist_y<=beta))
			{
				neighbors.insert(*V_iter);
			}
		}
		++V_iter;
	}
	// delete hd_mat
	for(unsigned int i=0; i<t; ++i)
	{
		delete[] hd_mat[i];
	}
	delete[] hd_mat;
}

/////////////////////////////////////////////////////////////////////////

/*
.Function.hasAtLeastOneOccurrence:
..summary:Checks if a given l-mer occurs at least once in a given sequence.
..cat:Motif Finding
..signature:hasAtLeastOneOccurrence(l_mer_begin,seq_begin,seq_end,l,d,is_exact)
..param.l_mer_begin:An iterator pointing to the beginning of a given l-mer pattern.
...type:Concept.Iterator Iterator
....remarks:Standard conform iterator
...type:Shortcut.DnaIterator
....remarks:Iterator for @Shortcut.DnaString@ (a string of @Spec.Dna@).
....see:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
....remarks:Iterator for @Shortcut.Peptide@ (a string of @Spec.AminoAcid@).
....see:Shortcut.PeptideIterator
..param.seq_begin:An iterator pointing to the beginning of a given sequence which is either
              a string of @Spec.Dna@ or a string of @Spec.AminoAcid@. 
...type:Concept.Iterator Iterator
....remarks:Standard conform iterator
...type:Shortcut.DnaIterator
....remarks:Iterator for @Shortcut.DnaString@ (a string of @Spec.Dna@).
....see:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
....remarks:Iterator for @Shortcut.Peptide@ (a string of @Spec.AminoAcid@).
....see:Shortcut.PeptideIterator
..param.seq_end:An iterator pointing to the end of a given sequence pattern which is either
            a string of @Spec.Dna@ or a string of @Spec.AminoAcid@.  
...type:Concept.Iterator Iterator
....remarks:Standard conform iterator
...type:Shortcut.DnaIterator
....remarks:Iterator for @Shortcut.DnaString@ (a string of @Spec.Dna@).
....see:Shortcut.DnaIterator
...type:Shortcut.PeptideIterator
....remarks:Iterator for @Shortcut.Peptide@ (a string of @Spec.AminoAcid@).
....see:Shortcut.PeptideIterator
..param.l:The size of the motif.
..param.d:The number of substitutions.
..param.is_exact:The size of Hamming distance
...type:$bool$
*/

template<typename TStringIter, typename TType>
bool
hasAtLeastOneOccurrence(TStringIter l_mer_begin,
		 TStringIter seq_begin,
		 TStringIter seq_end,
		 TType const & l,
		 TType const & d,
		 bool const & is_exact)
{
	bool result = false;
	TType counter = 0;
	while( (seq_begin!=(seq_end-l+1)) &  counter<1)
	{
		if(is_exact)
		{
			counter += (hammingDistance(seq_begin, seq_begin+l, l_mer_begin)==d) ? 1 : 0;
		}
		else
		{
			counter += (hammingDistance(seq_begin, seq_begin+l, l_mer_begin)<=d) ? 1 : 0;
		}		
		++seq_begin;
	}

	if(counter==1)
	{
		result = true;
	}
	return result;
}

/////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...