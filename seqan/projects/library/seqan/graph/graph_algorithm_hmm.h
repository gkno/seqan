#ifndef SEQAN_HEADER_GRAPH_ALGORITHM_HMM_H
#define SEQAN_HEADER_GRAPH_ALGORITHM_HMM_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////


template<typename TAlphabet, typename TProb = double, typename TSpec = Default>
class HiddenMarkovModel;

//////////////////////////////////////////////////////////////////////////////


template<typename TAlphabet, typename TProb, typename TSpec>
class HiddenMarkovModel 
{
	public:
		typedef String<TProb> TEmissionProb;
		typedef String<TProb> TTransitionProb;
		typedef String<TProb> TInitialProb;
		
		TTransitionProb data_state;
		TEmissionProb data_emission;
		TInitialProb data_initial;

	public:
		HiddenMarkovModel()
		{
			SEQAN_CHECKPOINT
		}

		~HiddenMarkovModel() 
		{
			SEQAN_CHECKPOINT
		}

		HiddenMarkovModel(HiddenMarkovModel const & _other)
		{
			SEQAN_CHECKPOINT
			data_state = _other.data_state;
			data_emission = _other.data_emission;
			data_initial = _other.data_initial;
		}

		HiddenMarkovModel const& 
		operator = (HiddenMarkovModel const& _other) 
		{
			SEQAN_CHECKPOINT
			if (this == &_other) return *this;
			data_state = _other.data_state;
			data_emission = _other.data_emission;
			data_initial = _other.data_initial;
			return *this;
		}
};


//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TProb, typename TSpec, typename TSize, typename TTransitionProb, typename TEmissionProb, typename TInitialProb>
inline void
initializeModel(HiddenMarkovModel<TAlphabet, TProb, TSpec>& hmm,
				TSize const nStates,
				TSize const alphSize,
				TTransitionProb const& trans,
				TEmissionProb const& emis,
				TInitialProb const& init)
{
	SEQAN_CHECKPOINT
		
	typedef HiddenMarkovModel<TAlphabet, TProb, TSpec> THMM;
	typedef typename Iterator<String<TProb> >::Type TStringIter;

	// Resize all tables
	resize(hmm.data_state, nStates * nStates);
	resize(hmm.data_emission, nStates * alphSize);
	resize(hmm.data_initial, nStates);

	// Assign values
	TStringIter it = begin(hmm.data_initial);
	for(TSize i=0;i<nStates;++i) {
		value(it) = getValue(init, i);
		goNext(it);
	}
	it = begin(hmm.data_state);
	for(TSize i=0;i<nStates;++i) {
		for(TSize j=0;j<nStates;++j) {
			value(it) = getValue(trans, i*nStates + j);
			goNext(it);
		}
	}
	it = begin(hmm.data_emission);
	for(TSize i=0;i<nStates;++i) {
		for(TSize j=0;j<alphSize;++j) {
			value(it) = getValue(emis, i*alphSize + j);
			goNext(it);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TProb, typename TSpec, typename TSize, typename TSequence>
inline void
generateSequence(HiddenMarkovModel<TAlphabet, TProb, TSpec> const& hmm,
				 TSize const& len,
				 TSequence& seq)
{
	SEQAN_CHECKPOINT		
	typedef HiddenMarkovModel<TAlphabet, TProb, TSpec> THMM;
	typedef typename Iterator<String<TProb> >::Type TStringIter;

	// Initialization
	mtRandInit();
	clear(seq);
	resize(seq, len);
	TSize nStates = length(hmm.data_initial);
	TSize alphSize = length(hmm.data_emission) / nStates;

	// Initial state
	TProb rmd = (TProb) (mtRand() % 101) / 100.0;
	TSize current_state = nStates - 1;
	TProb incSum = 0;
	for(TSize i=0; i<nStates; ++i) {
		incSum += getValue(hmm.data_initial, i);
		if (incSum > rmd) {
			current_state = i;
			break;
		}
	}

	// Generate sequence
	for(TSize i=0; i<len; ++i) {
		// Generate random character
		rmd = (TProb) (mtRand() % 101) / 100.0;
		incSum = 0;
		seq[i] = TAlphabet(alphSize - 1);
		for(TSize j=0; j<alphSize; ++j) {
			incSum += getValue(hmm.data_emission, current_state*alphSize+j);
			if (incSum > rmd) {
				seq[i] = TAlphabet(j);
				break;
			}
		}

		//std::cout << current_state << ": " << (unsigned int) seq[i] << std::endl;

		// Generate random next state
		rmd = (TProb) (mtRand() % 101) / 100.0;
		incSum = 0;
		TSize new_state = nStates - 1;
		for(TSize j=0; j<nStates; ++j) {
			incSum += getValue(hmm.data_state, current_state*nStates+j);
			if (incSum > rmd) {
				new_state = j;
				break;
			}
		}
		current_state = new_state;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TProb, typename TSpec, typename TSequence>
inline TProb
forwardAlgorithm(HiddenMarkovModel<TAlphabet, TProb, TSpec> const& hmm,
				 TSequence const& seq)
{
	SEQAN_CHECKPOINT		
	typedef HiddenMarkovModel<TAlphabet, TProb, TSpec> THMM;
	typedef typename Size<THMM>::Type TSize;
	typedef typename Iterator<String<TProb> >::Type TStringIter;

	// Initialization
	TSize nStates = length(hmm.data_initial);
	TSize len = length(seq);
	TSize alphSize = length(hmm.data_emission) / nStates;
	String<TProb> alpha;
	resize(alpha, nStates * len);

	// Initialization
	for(TSize i=0; i<nStates; ++i) {
		value(alpha, i) = getValue(hmm.data_initial, i) * getValue(hmm.data_emission, i*alphSize + (TSize) seq[0]);
	}

	// Recursion
	for(TSize pos=0; pos<len-1; ++pos) {
		for(TSize j=0; j<nStates; ++j) {
			TProb sum = 0;
			for(TSize i=0; i<nStates; ++i) {
				sum += getValue(alpha, pos * nStates + i) * getValue(hmm.data_state, i*nStates + j);
			}
			value(alpha, (pos+1) * nStates + j) = sum * getValue(hmm.data_emission, j*alphSize + (TSize) seq[pos+1]);
		}
	}

	// Termination
	TProb total = 0;
	for(TSize i=0; i<nStates; ++i) {
		total += getValue(alpha, (len-1) * nStates + i);
	}

	return total;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TAlphabet, typename TProb, typename TSpec, typename TIDString>
inline void
write(TFile &,
	  HiddenMarkovModel<TAlphabet, TProb, TSpec> const& hmm,
	  TIDString const &,
	  Raw)
{
	typedef HiddenMarkovModel<TAlphabet, TProb, TSpec> const THMM;
	typedef typename Size<THMM>::Type TSize;
	typedef typename Iterator<String<TProb> const, Rooted>::Type TStringIter;
	
	TSize nStates = length(hmm.data_initial);
	TSize alphSize = length(hmm.data_emission) / nStates;
	std::cout << "Number of states: " << nStates << std::endl;
	TStringIter it = begin(hmm.data_initial);
	TStringIter itEnd = end(hmm.data_initial);
	std::cout << "Initial state distribution:" << std::endl;
	while(it != itEnd) {
		std::cout << *it;
		++it;
		if (it != itEnd) std::cout << ", ";
	}
	std::cout << std::endl;
	it = begin(hmm.data_state);
	itEnd = end(hmm.data_state);
	std::cout << "State transition probabilities:" << std::endl;
	while(it != itEnd) {
		std::cout << *it;
		++it;
		if (position(it) % nStates == 0) {
			if (it!= itEnd) std::cout << std::endl;
		} else {
			std::cout << ", ";
		}
	}
	std::cout << std::endl;
	it = begin(hmm.data_emission);
	itEnd = end(hmm.data_emission);
	std::cout << "State emission probabilities:" << std::endl;
	while(it != itEnd) {
		std::cout << *it;
		++it;
		if (position(it) % alphSize == 0) {
			if (it!= itEnd) std::cout << std::endl;
		} else {
			std::cout << ", ";
		}
	}
	std::cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStream, typename TAlphabet, typename TProb, typename TSpec>
inline TStream &
operator <<(TStream & target, 
			HiddenMarkovModel<TAlphabet, TProb, TSpec> const& source)
{
	SEQAN_CHECKPOINT
	write(target, source);
	return target;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
