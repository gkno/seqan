#define SEQAN_PROFILE		// enable time measuring
#define FIONA_ALLOWINDELS	// allow for indels (chooses a less compact FragmentStore)
//#define FIONA_MEMOPT		// small suffix array values (<16mio reads of length <256)
#if defined(_OPENMP)        // Only enable parallelism in fiona if OpenMP is enabled.
#define FIONA_PARALLEL		// divide suffix tree into subtrees for each possible 3-gram
							// and use process subtrees in parallel
#endif  // defined(_OPENMP)

// debugging
//#define SEQAN_DEBUG_INDEX
//#define SEQAN_DEBUG
//#define SEQAN_VERBOSE
//#define SEQAN_VVERBOSE

#if defined(_OPENMP)
#include <omp.h>
#define SEQAN_PARALLEL
#if !defined(__MINGW32__) || defined(__amd64__) || defined(__x86_64__) || defined(__ia64__)
// The parallel STL generates warnings in MinGW: "...parallel/compatibility.h:167:42: note: #pragma message: slow __fetch_and_add_64".
// Thus, we do not enable it in this case.
#define _GLIBCXX_PARALLEL
#endif  // #if !defined(__MINGW32__) || defined(__amd64__) || defined(__x86_64__) || defined(__ia64__)
#endif  // #if defined(_OPENMP)

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
//#include <sys/resource.h>


#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <seqan/misc/misc_cmdparser.h>

// TODO(holtgrew): This raises a warning with Boost 1.42. Deactivate warnings, activate again afterwards. The correct #pragma has to be used for each supported compiler.
//#include <boost/math/distributions/normal.hpp>

//#define MEDIAN

using namespace seqan;

#ifdef FIONA_ALLOWINDELS

	// NOTE:
	// Currently we have to change the StringSet spec of the readSeqStore
	// to Owner as ConcatDirect<> (default) is not able to notice if a read
	// changes its size (after correction).
	struct FionaStoreConfig:
		public FragmentStoreConfig<>
	{
		typedef Owner<>	TReadSeqStoreSpec;
	};

	typedef FragmentStore<void, FionaStoreConfig> TFionaFragStore;

#else

	typedef FragmentStore<> TFionaFragStore;

#endif

typedef Index<TFionaFragStore::TReadSeqStore, IndexWotd<> > TFionaIndex; 
typedef Index<TFionaFragStore::TReadSeqStore, IndexQGram< Shape<Dna5, UngappedShape<4> > > > TFionaQgramIndex;


struct FionaPoisson_;
struct FionaExpected_;
struct FionaCount_ ;
struct FionaPoissonSens_ ;

typedef Tag<FionaPoisson_> const FionaPoisson;
typedef Tag<FionaPoissonSens_> const FionaPoissonSens;
typedef Tag<FionaExpected_> const FionaExpected;
typedef Tag<FionaCount_> const FionaCount;


struct FionaCorrectedError
{
	unsigned int   correctReadId;
//	unsigned int   occurrences;
	unsigned short errorPos;
	unsigned short correctPos;
	unsigned short overlap;
	signed   char  indelLength;		// 0..mismatch, <0..deletion, >0..insertion
//	unsigned char  mismatches;
};

struct FionaOptions
{
	__int64 genomeLength;
	double strictness;
	unsigned acceptedMismatches;
	int maxIndelLength;
	int fromLevel;
	int toLevel;
	unsigned cycles;
	double errorrate;
    int debugRead, corrRead;

	FionaOptions()
	{
		genomeLength = 0;
		strictness = 0.0001;
		acceptedMismatches = 1;
		maxIndelLength = 0;
		cycles = 0;
		fromLevel = 0;
		toLevel = 0;
		errorrate = 0;
        debugRead = -1;
        corrRead = -1;
	}
};

namespace seqan  {

/*restriction for certain levels - between max and min, also table with frequency may be to use eventually TODO*/
	struct FionaNodeConstraints 
	{ 
		 unsigned replen_max; 
		 unsigned replen_min;
     std::map<unsigned,double> frequency;
	}; 

	template <>
	struct Cargo<TFionaIndex> { 
		typedef FionaNodeConstraints Type; 
	}; 

#ifdef FIONA_MEMOPT

	typedef Pair<
		unsigned,
		unsigned,
		BitCompressed<24, 8>	// max. 16M reads of length < 256
	> TSAValue;

#else

	typedef Pair<
		unsigned,				// many reads
		unsigned short,			// of arbitrary length
		Compressed
	> TSAValue;

#endif

	template <>
	struct SAValue< TFionaIndex > 
	{
		typedef TSAValue Type;
	};

	template <>
	struct SAValue< TFionaQgramIndex > 
	{
		typedef TSAValue Type;
	};

	// use a mmap string for storing the q-grams
	template <>
	struct Fibre< TFionaQgramIndex, QGramSA > 
	{
		typedef String<TSAValue, MMap<> > Type;
	};

#ifdef FIONA_PARALLEL

	template <>
	struct Fibre< TFionaIndex, QGramSA > 
	{
		typedef Fibre< TFionaQgramIndex, QGramSA >::Type TSA;
		typedef Infix<TSA>::Type Type;
	};

#endif

	/*TODO THIS FONCTION CAN BE CHANGED FOR THE FREQUENCY - here just one experience*/
	/*by the use also the frequency for A,T,G,C*/
	/*higher frequency - high level as min in which we will begin the searching*/

	/*hide the node between certain level*/
	template <typename TSpec>
	inline bool nodeHullPredicate(Iter<TFionaIndex, VSTree<TopDown<TSpec> > > &it)
	{
		return parentRepLength(it) <= cargo(container(it)).replen_max;
	}

	template <typename TSpec>
	inline bool nodePredicate(Iter<TFionaIndex, VSTree<TopDown<TSpec> > > &it)
	{
		return true;
		FionaNodeConstraints &cons = cargo(container(it));
		unsigned repLen = parentRepLength(it);

		/*TODO may utilise >=*/
		return cons.replen_min <= repLen && repLen <= cons.replen_max;
	}

}




/*TODO cumulative poisson distribution determinated just once for a level, not at each node like now*/

/*matching string*/
inline bool strContains(std::string const & inputStr, std::string const & searchStr)
{
	return inputStr.find(searchStr) != std::string::npos;
}


// Expected value - general, use if all reads have the same length
template <typename TExpectedValues, typename TReadSet, typename TGenomeSize>
void expectedValueEqualReadsLength(TExpectedValues & expected, TReadSet const & readSet, TGenomeSize const genomeLength)
{
	//
	//	E(m) = (read_length - suffix_length + 1) * numberReads / genomeLength
	//

	// without reverse complement
	unsigned readCount = length(readSet) / 2;
	unsigned readLength = length(readSet[0]);

	clear(expected);
	for (unsigned suffixLength = 0; suffixLength <= readLength; ++suffixLength)
		append(expected, (readLength - suffixLength + 1) * readCount / (double)genomeLength);
}

// Expected value for set of reads with different length
template <typename TExpectedValues, typename TReadSet, typename TGenomeSize>
void expectedValueTheoretical(TExpectedValues & expected, TReadSet const & readSet, TGenomeSize const genomeLength)
{
	//
	//	E(m) = (read_length - suffix_length + 1) * numberReads / genomeLength
	//
	
	// calculate a read length histogram
	// excluding reverse complements
	String<unsigned> readLengthHistogram;
	for (unsigned i = 0; i < length(readSet); i += 2)
	{
		unsigned readLength = length(readSet[i]);
		if (readLength >= length(readLengthHistogram))
			resize(readLengthHistogram, readLength + 1, 0);
		++readLengthHistogram[readLength];
	}

	// a = read_length - suffix_length + 1
	clear(expected);
	resize(expected, length(readLengthHistogram), 0.0);
	for (unsigned i = 0; i < length(readLengthHistogram); ++i)
	{
		for (unsigned suffixLength = 0; suffixLength <= i; ++suffixLength)
		{
			double a = i - suffixLength + 1;
			double ratio = a * readLengthHistogram[i] / (double)genomeLength;
			expected[suffixLength] += ratio;
		}
	}
}

/* Standard Deviation */
template <typename TDeviationValues, typename TReadSet, typename TGenomeSize>
void standardDeviation(TDeviationValues & deviation, TReadSet const & readSet, TGenomeSize const genomeLength)
{
	// without reverse complement
	unsigned readCount = length(readSet) / 2;
	unsigned readLength = length(readSet[0]);

	//
	//	SD(m)= numberReads*((read_length - suffix_length + 1)/genomeLength 
	//			- (read_length - suffix_length + 1)^2/genomeLength ^2)
	//

	double valueFirst;
	double valueSecond;
	resize(deviation, readLength + 1);
	for (unsigned suffixLength = 0; suffixLength <= readLength; ++suffixLength)
	{
		valueFirst  = (readLength - suffixLength + 1) / (double)genomeLength;
		valueSecond = pow(valueFirst, 2);
		deviation[suffixLength] = sqrt((valueFirst - valueSecond) * readCount);	
	}
}

// factorial
template <typename TValue>
inline double factorial(TValue n)
{
	double fact = 1.0;
	for (TValue i = 2; i <= n; ++i)
		fact *= i;
	return fact;
}

// cumulative poisson distribution
template <typename TValue, typename TMean>
inline double pValue(TValue k, TMean mean)
{
	// return gsl_cdf_poisson_P(k,mean);
	double negExp = exp(-mean);
	double pValue = 0.0;

	for (TValue i = 0; i <= k; ++i)
		pValue += pow(mean, (double)i) * negExp / factorial(i);
	return pValue;
}

/*estimated a median value for a given level*/
template < typename TIndex, class TSpec >
double medianLevel(Iter<TIndex, VSTree<TSpec> > iter){

	double totalOccs = 0.0;
	double sumMedian = 0.0;
	double median = 0.0;
	double mediumTotalOccs = 0.0;

  std::map<unsigned, unsigned> vectorOccurences;

	goBegin(iter);
	for (; !atEnd(iter); ++iter)
	{
		unsigned numOccs = countOccurrences(iter);
		++vectorOccurences[numOccs];
		totalOccs += numOccs;
	}

	mediumTotalOccs = totalOccs / 2.0;

  std::map<unsigned,unsigned>::iterator iterMap;
	for (iterMap = vectorOccurences.begin (); iterMap != vectorOccurences.end (); ++iterMap)
	{
		sumMedian += iterMap->second*iterMap->first;
		if (sumMedian >= mediumTotalOccs)
		{
			median = iterMap->first;
			break;
		}
	}
	return median;
}

template <typename TPercentage, typename TSize>
inline double probabilityOneError(TPercentage percentageErr, TSize repLen)
{
	return 1.0 - pow(1.0 - percentageErr, repLen);
}

template <typename TFragmentStore, typename TCorrection>
inline void _dumpCorrection(
	TFragmentStore &store,
	TCorrection const &correction,
	unsigned errorReadId)
{
	std::cerr << std::endl;
	std::cerr << "error___read   \t" << store.readSeqStore[errorReadId][correction.errorPos] << '\t';
	for (unsigned i = 0; i < correction.correctPos; ++i)
		std::cout << ' ';
    std::cout << store.readSeqStore[errorReadId] << std::endl;
	std::cerr << "correct_read   \t" << store.readSeqStore[correction.correctReadId][correction.correctPos] << '\t';
	for (unsigned i = 0; i < correction.errorPos; ++i)
		std::cout << ' ';
    std::cerr << store.readSeqStore[correction.correctReadId] << std::endl;
	std::cerr << "error___read_id\t" << errorReadId << std::endl;
	std::cerr << "correct_read_id\t" << correction.correctReadId << std::endl;
	std::cerr << "error_pos      \t" << correction.errorPos << std::endl;
	std::cerr << "correct_pos    \t" << correction.correctPos << std::endl;
	std::cerr << "overlap        \t" << correction.overlap << std::endl;
	std::cerr << "indel_len      \t" << (int)correction.indelLength << std::endl;
}

/*change the erroneous nucleotide in all reads identify with errors*/
template <typename TFragmentStore, typename TCorrections>
void applyReadErrorCorrections(
	TFragmentStore &store,
	TCorrections const &corrections,
    FionaOptions const & options)
{
	typedef typename Value<TCorrections>::Type TCorrection;
	int readCount = length(corrections);

	// we have to make a temp-copy in order to use original (not corrected) reads for correction
	StringSet<typename TFragmentStore::TReadSeq> originalReads;
	resize(originalReads, length(store.readSeqStore), Exact());

#ifdef FIONA_PARALLEL
	#pragma omp parallel for schedule(guided)
#endif
	for (int readId = 0; readId < readCount; ++readId)
	{
		TCorrection const &corr = corrections[readId];
		if (corr.overlap != 0 && corr.indelLength <= 0)
			originalReads[corr.correctReadId] = store.readSeqStore[corr.correctReadId];
	}


#ifdef FIONA_PARALLEL
	#pragma omp parallel for
#endif
	for (int readId = 0; readId < readCount; ++readId)
	{
		TCorrection const &corr = corrections[readId];
		if (corr.overlap == 0) continue;

        if (readId == options.debugRead)
        {
            std::cout << "BEFORE:" << std::endl;
            _dumpCorrection(store, corr, readId);
        }

    std::ostringstream m;
		if (strContains(toCString(store.readNameStore[readId]), "corrected"))
			m << "," << corr.errorPos;
		else
			m << " corrected: " << corr.errorPos;

#ifdef SEQAN_VERBOSE
		_dumpCorrection(store, corr, readId);
#endif
		append(store.readNameStore[readId], m.str());

		if (corr.indelLength == 0)
			store.readSeqStore[readId][corr.errorPos] = originalReads[corr.correctReadId][corr.correctPos];
#ifdef FIONA_ALLOWINDELS
		else if (corr.indelLength > 0)
            erase(store.readSeqStore[readId], corr.errorPos, corr.errorPos + corr.indelLength);
		else
			insert(store.readSeqStore[readId], corr.errorPos, infix(originalReads[corr.correctReadId], corr.correctPos, corr.correctPos + -corr.indelLength));
#endif
        if (readId == options.debugRead)
        {
            std::cout << "AFTER:" << std::endl;
            _dumpCorrection(store, corr, readId);
        }
#ifdef SEQAN_VERBOSE
		std::cout << "corrected:";
		for (unsigned i = 0; i < corr.correctPos; ++i)
			std::cout << ' ';
		std::cout << store.readSeqStore[readId] << std::endl;
#endif
	}
}

template <typename TObserved, typename TExpected, typename TStrictness, typename TErrorRate, typename TPrefixLen>
inline bool potentiallyErroneousNode(
	TObserved observed,
	TExpected expected,
	TStrictness strictness,
	TErrorRate,
	TPrefixLen,
	FionaPoisson const)
{
	// compare the cumulative poisson distribution with the p-value (strictness)
//	return pValue(observed, expected) <= strictness;
    double negExp = exp(-expected);
    double pValue = 0.0;
	double pow = 1.0;
	double fact = 1.0;

    for (TObserved i = 0; i <= observed && pValue <= strictness; ++i, fact *= i){
        pValue += pow * negExp / fact;
		pow *= expected ;
	}
	return pValue <= strictness;
}

template <typename TObserved, typename TExpected, typename TStrictness, typename TErrorRate, typename TPrefixLen>
inline bool potentiallyErroneousNode(
	TObserved observed,
	TExpected expected,
	TStrictness,
	TErrorRate,
	TPrefixLen,
	FionaExpected const)
{
	// compare the weight for a node with a cutoff given by the strictness param.
	return observed < expected;
}

template <typename TObserved, typename TExpected, typename TStrictness, typename TErrorRate, typename TPrefixLen>
inline bool potentiallyErroneousNode(
	TObserved observed,
	TExpected,
	TStrictness cutoff,
	TErrorRate,
	TPrefixLen,
	FionaCount const)
{
	// compare the weight for a node with a fixed cutoff
	return observed < cutoff;
}



template <typename TObserved, typename TExpected, typename TStrictness, typename TErrorRate, typename TPrefixLen>
inline bool potentiallyErroneousNode(
	TObserved observed,
	TExpected expected,
	TStrictness falsenegrate,
	TErrorRate errorrate,
    TPrefixLen prefixlen,
	FionaPoissonSens const)
{
	// Poisson based threshold, given a fixed percentage of missed errors (1 - min sensitivity)
	// consider only the cases with one and two errors
	// the average error rate and the expected value allow to compute the expected count for an error.
	double sensitivity = 1 - falsenegrate ;
	double noerrlm2 = pow(1-errorrate, prefixlen - 2.0);
	double noerrlm1 = noerrlm2 * (1-errorrate) ;
	double errexp1err = expected * noerrlm1 * errorrate ;
	double errexp2err = expected * noerrlm2 * errorrate *errorrate ;
	double perr1 = prefixlen * noerrlm1* errorrate  ;
	double perr2 = (prefixlen * (prefixlen -1) / 2 ) * noerrlm2 * errorrate *errorrate ;
	double sc = perr1 + perr2 ;
	perr1 /= (sc) ; 
	perr2 /= (sc) ; 
	double negExp1err = exp(-errexp1err); //
	double negExp2err = exp(-errexp2err) ;
    double probaerror = 0.0;
	double pow1err = 1.0;
	double pow2err = 1.0;
	double fact = 1.0;
	//std::cout << "Level " << prefixlen << " and params for correc\n" ;
	//std::cout << "Expected and comp: " << expected << " - " << errexp1err << " - " << errexp2err << "\n" ;
	//std::cout << "the basic observed count: " << observed << ", perr is " << perr1 << " - " << perr2  << "\n";
	TObserved i = 0 ;
    for (i = 0; i <= observed && probaerror <= sensitivity; ++i, fact *= i){
        probaerror += perr1 * pow1err * negExp1err / fact;
		probaerror += perr2 * pow2err * negExp2err / fact;
		pow1err *= errexp1err ; 
		pow2err *= errexp2err ;
	}
	//std::cout << "Stopped at observed value **" << i-1 << "** for a sens. of " << probaerror << "and a thr at " << sensitivity << "\n" ;
	return probaerror <= sensitivity;
}

struct Overlap
{
    unsigned overlapSumLeft;
    unsigned overlapSumRight;
    unsigned readId;
    unsigned correctPos;
};


/*detect and repare the reads with errors*/
template <typename TTreeIterator, typename TFragmentStore, typename TCorrections, typename TAlgorithm>
void traverseAndSearchCorrections(
	TTreeIterator iter,
	TFragmentStore &store,
	TCorrections & corrections,
	String<double> & expectedTheoretical,
	FionaOptions const & options,
	Tag<TAlgorithm> const alg)
{
	typedef typename Container<TTreeIterator>::Type TFionaIndex;
	typedef typename Fibre<TFionaIndex, FibreText>::Type TReadSet;
	typedef typename Fibre<TFionaIndex, FibreSA>::Type TSA;
	typedef typename Infix<TSA const>::Type TOccs;
	typedef typename Iterator<TOccs, Standard>::Type TOccsIterator;
	typedef typename Value<TReadSet>::Type TRead;
	typedef typename Value<TRead>::Type TValue;
	typedef typename Iterator<TRead, Standard>::Type TReadIterator;
	typedef typename Value<TCorrections>::Type TCorrection;

	unsigned readCount = length(store.readSeqStore) / 2;
	
	bool debug = false;

	String<TOccs> correctCandidates;    
	const TValue unknownChar = unknownValue<TValue>();
    String<Overlap> bestCorrection; // for a given branch store for every indel-size the best correction
    resize(bestCorrection, options.maxIndelLength * 2 + 1);

	for (goBegin(iter); !atEnd(iter); )                     // do a DFS walk
	{
		unsigned commonPrefix = parentRepLength(iter);			// length of parent label

#ifdef FIONA_DEBUG
		{
			TTreeIterator iterSibling(iter);
			goUp(iterSibling);
			if (representative(iterSibling) == "AGGCGCCTCCTCAA")
			{
				std::cout << "GOTIT!!!!" << std::endl;
				std::cout << "GOTIT!!!!" << std::endl;
				std::cout << "GOTIT!!!!" << std::endl;

				{
					std::cerr << std::endl << "ERROR CANDIDATES ARE:" << std::endl;

					std::cerr << "support         \t" << countOccurrences(iter) << std::endl;
					std::cerr << "occurrences     \t";
					for (unsigned i = 0; i < countOccurrences(iter); ++i)
						std::cerr << '\t' << getOccurrences(iter)[i] << '\t' << suffix(indexText(container(iter)), getOccurrences(iter)[i]) << std::endl;
					std::cerr << std::endl;

					std::cerr << std::endl << "SIBLINGS ARE:" << std::endl;

					goDown(iterSibling);
					do
					{
						if (value(iter).range != value(iterSibling).range)
						{
							std::cerr << "support         \t" << countOccurrences(iterSibling) << std::endl;
							std::cerr << "occurrences     \t";
							for (unsigned i = 0; i < countOccurrences(iterSibling); ++i)
								std::cerr << '\t' << getOccurrences(iterSibling)[i] << '\t' << suffix(indexText(container(iterSibling)), getOccurrences(iterSibling)[i]) << std::endl;
							std::cerr << std::endl;
						}
					} while (goRight(iterSibling));					
				}

			}
		}
#endif //FIONA_DEBUG

		SEQAN_ASSERT_LT(commonPrefix + 1, length(expectedTheoretical));
		TValue firstEdgeChar = parentEdgeFirstChar(iter);
		if (firstEdgeChar != unknownChar &&
			!potentiallyErroneousNode(countOccurrences(iter), expectedTheoretical[commonPrefix+1], options.strictness, options.errorrate, commonPrefix + 1 ,alg))
		{
			goNext(iter);
			continue;
		}

		//
		//	get the id and position (suffix begin) for suspected nodes
		//	for which we can find a more optimal correction
		//
		TOccs errorCandidates = getOccurrences(iter);

#ifdef FIONA_DEBUG		
		bool debugReadFound = false;
		for (unsigned i=0; i<length(errorCandidates); ++i)
			if (((errorCandidates[i].i1 < readCount)? errorCandidates[i].i1: errorCandidates[i].i1 - readCount) == (unsigned)options.debugRead)
				debugReadFound = true;

		bool corrReadFound = false;
		if (debugReadFound)
		{
			TTreeIterator iterSibling(iter);
			goUp(iterSibling); goDown(iterSibling);
            do
            {
                if (value(iter).range != value(iterSibling).range &&
                    parentEdgeFirstChar(iterSibling) != unknownChar)
                {
					TOccs corrCandidates = getOccurrences(iter);
					for (unsigned i=0; i<length(corrCandidates); ++i)
						if (((corrCandidates[i].i1 < readCount)? corrCandidates[i].i1: corrCandidates[i].i1 - readCount) == (unsigned)options.corrRead)
							corrReadFound = true;
                }
            } while (goRight(iterSibling));
		}
		if (debugReadFound || corrReadFound)
			std::cerr << std::endl << "I FOUND " << debugReadFound << " " << corrReadFound << std::endl;
#endif //FIONA_DEBUG


		/*copy the iterator for iterate over the siblings*/
//		typename Iterator<TFionaIndex, TopDown<> >::Type iterSibling(container(iter), nodeUp(iter));
		TTreeIterator iterSibling(iter);
		goUp(iterSibling);

		//
		//	potential reads for make the correction,
		//	because at the same level, with the same prefix
		//

		clear(correctCandidates);
		if (!goDown(iterSibling))
			SEQAN_ASSERT_FAIL("going up and down failed!?");

		// pick potentially correct reads
		do
		{
			if (value(iter).range != value(iterSibling).range &&
				parentEdgeFirstChar(iterSibling) != unknownChar &&
				!potentiallyErroneousNode(countOccurrences(iterSibling), expectedTheoretical[commonPrefix + 1], options.strictness, options.errorrate,commonPrefix +1, alg))
			{
				// save the id and position(where the suffix begin in the reads) in the table of IDs correct
				// also the number of occurrences
				appendValue(correctCandidates, getOccurrences(iterSibling));
			}
		} while (goRight(iterSibling));

        // if there aren't any, try the thickest branch
        if (empty(correctCandidates))
        {
            goUp(iterSibling);
            if (!goDown(iterSibling))
                SEQAN_ASSERT_FAIL("going up and down failed!?");

            unsigned thickestBranch = length(errorCandidates);
            do
            {
                if (value(iter).range != value(iterSibling).range &&
                    thickestBranch < countOccurrences(iterSibling) &&
                    parentEdgeFirstChar(iterSibling) != unknownChar)
                {
                    thickestBranch = countOccurrences(iterSibling);
                    clear(correctCandidates);
                    appendValue(correctCandidates, getOccurrences(iterSibling));
                }
            } while (goRight(iterSibling));
        }

		// continue if we haven't found any correct read
		if (empty(correctCandidates))
		{
#ifdef FIONA_DEBUG
			TOccsIterator errorRead = begin(errorCandidates, Standard());
			TOccsIterator errorReadEnd = end(errorCandidates, Standard());
			for (; errorRead != errorReadEnd; ++errorRead)
			{
				unsigned errorReadId = (*errorRead).i1;
				unsigned errorReadIdFwd = (errorReadId < readCount)? errorReadId: errorReadId - readCount;
				debug = (errorReadIdFwd == (unsigned)options.debugRead);

				unsigned positionError = (*errorRead).i2 + commonPrefix;
				if (debug)
				{
					std::cerr << "NO CORRECTING SIBLING FOUND !!!" << std::endl;
					std::cerr << "prefix_length   \t" << commonPrefix << std::endl;
					std::cerr << "error___read_id \t" << errorReadId << std::endl;
					std::cerr << "error_pos       \t" << positionError << std::endl;
					std::cerr << std::endl << "SIBLINGS ARE:" << std::endl;

					goUp(iterSibling); goDown(iterSibling);
					do
					{
						if (value(iter).range != value(iterSibling).range)
						{
							std::cerr << "support         \t" << countOccurrences(iterSibling) << std::endl;
							std::cerr << "occurrences     \t";
							for (unsigned i = 0; i < countOccurrences(iterSibling); ++i)
								std::cerr << getOccurrences(iterSibling)[i];
							std::cerr << std::endl;
						}
					} while (goRight(iterSibling));					
				}
			}
#endif // FIONA_DEBUG

			// don't descent edges beginning with N
			if (firstEdgeChar == unknownChar)
				goNextRight(iter);
			else
				goNext(iter);
			continue;
		}

		// make the comarison between the substrings(the suffix after the position of error
		TOccsIterator errorRead = begin(errorCandidates, Standard());
		TOccsIterator errorReadEnd = end(errorCandidates, Standard());
		for (; errorRead != errorReadEnd; ++errorRead)
		{
			unsigned errorReadId = (*errorRead).i1;
            unsigned errorReadIdFwd = (errorReadId < readCount)? errorReadId: errorReadId - readCount;
			debug = (errorReadIdFwd == (unsigned)options.debugRead);

#ifdef FIONA_DEBUG
			if (debug)
			{
                std::cout << "ErrorBranch: " << length(errorCandidates) << "   CorrectBranches: ";
                for(unsigned i=0;i<length(correctCandidates);++i)
				{
                    std::cout<<length(correctCandidates[i]) << ", ";
//					for(unsigned j=0;j<length(correctCandidates[i]);++j)
//						if ()
				}
                std::cout<<std::endl;
            }
#endif //FIONA_DEBUG
			// is already identify as erroneus
			TCorrection &corr = corrections[errorReadIdFwd];

			// Branch and Bound
			//
			// if already detect as erroneous(optimal)
			// the total length minus the position of error, and minus 1 mismatch

			// the position where the error is
			unsigned positionError = (*errorRead).i2 + commonPrefix;
/*
            // branch and bound has to be fixed to work with overlapSum
            int maxOverlapPossible = commonPrefix + length(store.readSeqStore[errorReadId]) - positionError - 1;
			if (options.maxIndelLength > 0)
				++maxOverlapPossible;

			if (corr.overlap > maxOverlapPossible) continue;
*/
			TReadIterator itEBegin = begin(store.readSeqStore[errorReadId], Standard());
			TReadIterator itEPrefixBegin = itEBegin + (*errorRead).i2;
			TReadIterator itEEnd = end(store.readSeqStore[errorReadId], Standard());

			for (unsigned c = 0; c < length(correctCandidates); ++c)
			{
				TOccsIterator corrRead = begin(correctCandidates[c], Standard());
				TOccsIterator corrReadEnd = end(correctCandidates[c], Standard());
                
                /////////////////////////////////////////////////////////////////////////////////
                // reset overlap sums
                Iterator<String<Overlap>, Standard>::Type corIt = begin(bestCorrection, Standard());
                Iterator<String<Overlap>, Standard>::Type corItEnd = end(bestCorrection, Standard());
                for (; corIt != corItEnd; ++corIt) 
                {
                    (*corIt).overlapSumLeft = 0;
                    (*corIt).overlapSumRight = 0;
                }
                
				for (; corrRead != corrReadEnd; ++corrRead)
				{					
					/////////////////////////////////////////////////////////////////////////////////
					// compare overlap left of the common prefix

					register TReadIterator itCLeft = begin(store.readSeqStore[(*corrRead).i1], Standard());
					register TReadIterator itE = itEBegin;
					unsigned acceptedMismatchesLeft = options.acceptedMismatches;
					
					int delta = (*errorRead).i2 - (*corrRead).i2;
					if (delta > 0)
						itE += delta;
					else
						itCLeft += -delta;

					for (; itE < itEPrefixBegin; ++itE, ++itCLeft)
						if (*itE != *itCLeft)
							if (--acceptedMismatchesLeft == MaxValue<unsigned>::VALUE) break;

					// too many mismatches left of the common prefix?
					if (acceptedMismatchesLeft == MaxValue<unsigned>::VALUE)
						continue;

                    unsigned overlapLeft = positionError;
                    if (delta > 0)
                        overlapLeft -= delta;

					// the position in the read until which there is the same prefix
					unsigned positionCorrect = (*corrRead).i2 + commonPrefix;

					// search the position detect as erroneous at the level of the read
					TReadIterator itCEnd = end(store.readSeqStore[(*corrRead).i1], Standard());

					/////////////////////////////////////////////////////////////////////////////////
					// compare overlap left of the common prefix
					for (int indel = -options.maxIndelLength; indel <= options.maxIndelLength; ++indel)
					{
                        itE = itEBegin + positionError;
						TReadIterator itC = itCLeft + commonPrefix;
						unsigned acceptedMismatches;

						if (indel == 0)
						{
							// mismatch
							++itE;
							++itC;
							acceptedMismatches = acceptedMismatchesLeft;
						}
						else if (indel > 0)
                        {
                            // gap in correct read (must be deleted in err. read)
                            itE += indel;
                            if (itE + indel >= itEEnd) continue;
                            acceptedMismatches = acceptedMismatchesLeft;
                        }
                        else
                        {
                            // gap in erroneous read (must be filled with an insertion in err. read)
                            itC += -indel;
                            if (itC + -indel >= itCEnd) continue;
                            acceptedMismatches = acceptedMismatchesLeft;
                        }
                        

						TReadIterator itEFirst = itE;
						for (; itE < itEEnd && itC < itCEnd; ++itE, ++itC)
							if (*itE != *itC)
								if (--acceptedMismatches == MaxValue<unsigned>::VALUE)
									break;
                        
                        // too many mismatches right of the common prefix?
                        if (acceptedMismatches == MaxValue<unsigned>::VALUE)
                            continue;


                        // correction candidate
                        Overlap &overlap = bestCorrection[((indel >= 0)? indel * 2: -indel * 2 - 1)];
                        overlap.overlapSumLeft += overlapLeft;
                        overlap.overlapSumRight += itE - itEFirst;
                        overlap.readId = (*corrRead).i1;
                        overlap.correctPos = positionCorrect;
					}
				}
			}
			
            Iterator<String<Overlap>, Standard>::Type corIt = begin(bestCorrection, Standard());
            Iterator<String<Overlap>, Standard>::Type corItEnd = end(bestCorrection, Standard());
            Iterator<String<Overlap>, Standard>::Type corBest = corItEnd;
            unsigned maxOverlapSum = 0;
            for (; corIt != corItEnd; ++corIt) 
            {
                unsigned overlapSum = (*corIt).overlapSumLeft + (*corIt).overlapSumRight;
                if (maxOverlapSum < overlapSum)
                {
                    maxOverlapSum = overlapSum;
                    corBest = corIt;
                }
            }
            // we could also put a minimum overlapSum constraint in here
            if (corBest == corItEnd) continue;// || maxOverlapSum < length(store.readSeqStore[errorReadId])) continue;
            
            int i = corBest - begin(bestCorrection, Standard());
            int bestIndelLength = ((i & 1) == 0)? i / 2: -((i + 1) / 2);
            
			if (debug)
			{
				std::cerr << std::endl;
				std::cerr << "error___read   \t" << store.readSeqStore[errorReadId] << '\t' << store.readSeqStore[errorReadId][positionError] << std::endl;
				std::cerr << "correct_read   \t" << store.readSeqStore[(*corBest).readId] << '\t' << store.readSeqStore[(*corBest).readId][(*corBest).correctPos] << std::endl;
				std::cerr << "error___read_id\t" << errorReadId << std::endl;
				std::cerr << "correct_read_id\t" << (*corBest).readId << std::endl;
				std::cerr << "error_pos      \t" << positionError << std::endl;
				std::cerr << "correct_pos    \t" << (*corBest).correctPos << std::endl;
				std::cerr << "overlapSum     \t" << maxOverlapSum << std::endl;
				std::cerr << "indel_len_alt  \t" << i << std::endl;
				std::cerr << "indel_len      \t" << bestIndelLength << std::endl;
                for (int i = 0; i < (int)length(bestCorrection); ++i)
                    std::cerr << (((i & 1) == 0)? i / 2: -((i + 1) / 2)) << '=' << bestCorrection[i].overlapSumRight << "  ";
                std::cerr << std::endl;
			}

            if (errorReadId >= readCount)
            {
                // switch to reverse-complements
                if ((*corBest).readId < readCount)
                    (*corBest).readId += readCount;
                else
                    (*corBest).readId -= readCount;

                // mirror positions
                positionError = length(store.readSeqStore[errorReadId]) - positionError;
                (*corBest).correctPos = length(store.readSeqStore[(*corBest).readId]) - (*corBest).correctPos;

                if (bestIndelLength == 0)
                {
                    --positionError;
                    --(*corBest).correctPos;
                } else if (bestIndelLength > 0)
                    positionError -= bestIndelLength;
                else
                    (*corBest).correctPos += bestIndelLength;
            }

            /*This is maybe not necessary, because we search to maximisate*/
            /*when the position or the nucleotide differ*/

            /*find longer commun part with certain nb accepted mismatch*/
#ifdef FIONA_PARALLEL
            #pragma omp critical
#endif
            {
                if  (corr.overlap < maxOverlapSum ||
                    (corr.overlap == maxOverlapSum &&
                        (corr.correctReadId > (*corBest).readId ||
                        (corr.correctReadId == (*corBest).readId && corr.correctPos > (*corBest).correctPos))))
                {
                    // update best correction
                    corr.correctReadId = (*corBest).readId;
                    corr.errorPos = positionError;
                    corr.correctPos = (*corBest).correctPos;
                    corr.overlap = maxOverlapSum;
                    corr.indelLength = bestIndelLength;
                }
            }
		}
		// don't descent edges beginning with N
		if (firstEdgeChar == unknownChar)
			goNextRight(iter);
		else
			goNext(iter);
	}
}


/*GC-content*/
/*fonction which allow to determine the frequency for each nucleotide*/
template < typename TFionaIndex, class TSpec >
std::map<unsigned,double>
determineFrequency(Iter< TFionaIndex, VSTree<TSpec> > iter)
{
   /*calculate the frequency for each nucleotide*/
   /*'A' = 0, 'C' = 1, 'G' = 2, 'T' = 3*/
  std::map<unsigned, double> frequency;

   goBegin(iter);

   int position =0;
   /*nombre total nucleotides*/
   double total=countOccurrences(iter);
   /*the first is A (alphabetical ordre)*/
   goDown(iter);
   frequency[position] = (countOccurrences(iter)/total);
   while(goRight(iter)){
	   position+=1;
	   frequency[position] = (countOccurrences(iter)/total);
   }

   /*table of frequency for each nucleotide*/
   return frequency;
}


/*construction Suffix Array */
template <typename TFragmentStore, typename TAlgorithm>
unsigned correctReads(
	TFragmentStore & store,
	FionaOptions & options,
	Tag<TAlgorithm> const alg)
{
	/*iterator with restrictions*/
	typedef Iterator<TFionaIndex, TopDown<ParentLinks<Preorder> > >::Type TConstrainedIterator;
	typedef FionaCorrectedError TCorrected;

	// append their reverse complements
	unsigned readCount = length(store.readSeqStore);
	Dna5String tmp;
	if (options.genomeLength != 1)
		for (unsigned i = 0; i < readCount; ++i)
		{
			tmp = store.readSeqStore[i];
			reverseComplement(tmp);
			appendValue(store.readSeqStore, tmp);
		}

	/*table with the theoreticals values*/
	String<double> expectedTheoretical;
	expectedValueTheoretical(expectedTheoretical, store.readSeqStore, options.genomeLength);

	if (IsSameType<TAlgorithm, FionaExpected_>::VALUE)
	{
		String<double> sd;
		standardDeviation(sd, store.readSeqStore, options.genomeLength);

		/*The strictness value allows to estimate the confidential intervall*/
		for (unsigned i = 0; i < length(expectedTheoretical); ++i)
		{
			double expectedTemporary = expectedTheoretical[i] - options.strictness * sd[i];
			
			/*If the connfidential intervall take value less than 1 ??? not sure for that*/
			/*if(expectedTemporary < 1){
				expectedTheoretical[i] = 1.1;
			}else{*/
				expectedTheoretical[i] = expectedTemporary;
			//}
			//if fixed
			//expectedTheoretical[i]=5;
		}
	}
	
	if (IsSameType<TAlgorithm, FionaCount_>::VALUE)
		for (unsigned i = 0; i < length(expectedTheoretical); ++i)
			expectedTheoretical[i] = options.strictness ;
	
	
	String<TCorrected> corrections;
	TCorrected zeroCorr = { 0, 0, 0, 0, 0 };
	resize(corrections, readCount, zeroCorr);
	
	if (IsSameType<TAlgorithm, FionaExpected_>::VALUE)
		std::cout << std::endl << "Method with expected value for each level" << std::endl;
	if (IsSameType<TAlgorithm, FionaPoisson_>::VALUE) 
		std::cout << std::endl << "Method with p-value and Poisson distribution" << std::endl;
	if (IsSameType<TAlgorithm, FionaPoissonSens_>::VALUE)
		std::cout << std::endl << "Method with sensitivity and Poisson distribution" << std::endl;
	if (IsSameType<TAlgorithm, FionaCount_>::VALUE)
		std::cout << std::endl << "Method with fixed count for each level" << std::endl;

	
	std::cout << "Searching..." << std::endl;
	SEQAN_PROTIMESTART(search);

#ifndef FIONA_PARALLEL
	// FIONA NON-PARALLEL SEARCH
		
	// construct suffix array of the set of reads
	std::cout << "Construct suffix array" << std::endl;
	SEQAN_PROTIMESTART(construct);
	TFionaIndex myIndex(store.readSeqStore);
	TConstrainedIterator myConstrainedIterator(myIndex);

	/*calculate the frequency for each nucleotide, didn't use for the moment*/
	/*'A' = 0, 'C' = 1, 'G' = 2, 'T' = 3*/
  std::map<unsigned,double> frequency;
	frequency = determineFrequency(myConstrainedIterator);
	
	std::cout << "Time required for suffix array construction : " << SEQAN_PROTIMEDIFF(construct) << " seconds." << std::endl;

	/*restrictions just for estimate the genome length if there is no data*/

#ifdef MEDIAN
	unsigned level = fromLevel;
	ofstream out("medianLevels.txt"); 
	ofstream median("medianForEachLevel.txt");
#endif

	if (options.genomeLength == 1)
	{
		std::cout << "Generating Hugues' stats file." << std::endl;
        std::ofstream stats("stats.txt");
		Iterator<TFionaIndex, TopDown<ParentLinks<Preorder> > >::Type it(myIndex);
		goBegin(it);
		CharString tmp;
		String<int> freq;
		stats << "branch\tlength\ttreeDep\tletter\treadPos\tfreq" << std::endl;
		while (!atEnd(it))
		{
			unsigned ofs = parentRepLength(it);
			tmp = parentEdgeLabel(it);

			Iterator<TFionaIndex, TopDown<ParentLinks<Preorder> > >::Type it2(it);
			bool branch = goDown(it2) && goRight(it2);

			for (unsigned i = 0; i < length(tmp); ++i)
			{
				clear(freq);
				for (unsigned j = 0; j < countOccurrences(it); ++j)
				{
					unsigned posInRead = getOccurrences(it)[j].i2 + ofs + i;
					if (length(freq) <= posInRead)
						resize(freq, posInRead + 1, 0);
					++freq[posInRead];
				}
				
				for (unsigned k = 0; k < length(freq); ++k)
					if (freq[k] != 0)
					{
						if ((i + 1 == length(tmp)) && branch)
							stats << "1\t";
						else
							stats << "0\t";
						stats << ofs + i + 1 << '\t' << nodeDepth(it) << '\t' << tmp[i];
						stats << '\t' << k << '\t' << freq[k] << '\t';// << representative(it) << '\t' << tmp << std::endl;
						stats << std::endl;					}
			}
			++it;
		}
		stats.close();
		std::cout << "Done." << std::endl;
		exit(0);
	}
			
	if (options.genomeLength == 0)
	{
#ifdef MEDIAN		
		for (; level < toLevel; ++level)
		{
			//int logRation = (log10(static_cast<double>(length(setReads)/2))/(log10(4.0)));
			//int l = logRation + 1;
			//std::cout << l << std::endl;
			
			cargo(myIndex).replen_min = level;
			cargo(myIndex).replen_max = level+2; 	
			cargo(myIndex).frequency = frequency;
			double numOccs = 0.0; 
			
			median << level << " " << medianLevel(myConstrainedIterator) << std::endl;
			goBegin(myConstrainedIterator);
			while (!atEnd(myConstrainedIterator))
			{
				if (parentRepLength(myConstrainedIterator) > level)
				{
					numOccs = countOccurrences(myConstrainedIterator);
					out << level << " " << numOccs << std::endl; 
				}
				++myConstrainedIterator;
			}
		}
#else
		//int logRation = log10(static_cast<double>(readCount)) / log10(4.0);
		//int l = logRation + 1;
		//std::cout << l << std::endl;
		cargo(myIndex).replen_min = options.fromLevel;
		cargo(myIndex).replen_max = options.fromLevel + 2; 	
		cargo(myIndex).frequency = frequency;

		double expectedValueGivenLevel = medianLevel(myConstrainedIterator);

		/*TODO the length for the reads is always the same, but for real cases it is better to change*/
		/*estimate the genome length by the use of expected value*/
		unsigned readLength = length(store.readSeqStore[0]);

		/* a = readLength - path_label + 1 */
		/*here plus 1 also because the level is between fromLevel and toLevel*/
		double a = readLength - options.fromLevel + 2;
		options.genomeLength = static_cast<__int64>(readCount * a / expectedValueGivenLevel);
		std::cout << "The estimated genome length is " << options.genomeLength << std::endl;
#endif
	}

	/*restrictions for the searching levels*/
	cargo(myIndex).replen_min = options.fromLevel;
	cargo(myIndex).replen_max = options.toLevel;
	cargo(myIndex).frequency = frequency;

	/*the core of the correction method*/
	traverseAndSearchCorrections(myConstrainedIterator, store, corrections, expectedTheoretical, options, alg);
	std::cout << "Time for searching between given levels: "<< SEQAN_PROTIMEDIFF(search) << " seconds." << std::endl;

#else 
	// FIONA PARALLEL SEARCH

	// construct q-gram index
	TFionaQgramIndex qgramIndex(store.readSeqStore);
	String<size_t> packages;

	SEQAN_PROTIMESTART(constructQgramExt);
	std::cout << "Construct external q-gram index ... " << std::flush;
	resize(indexSA(qgramIndex), _qgramQGramCount(qgramIndex), Exact());
	resize(indexDir(qgramIndex), _fullDirLength(qgramIndex), Exact());
	createQGramIndexExt(indexSA(qgramIndex), indexDir(qgramIndex), indexText(qgramIndex), indexShape(qgramIndex), stringSetLimits(indexText(qgramIndex)));
	std::cout << "done. (" << SEQAN_PROTIMEDIFF(constructQgramExt) << " seconds)" << std::endl;

	String<Dna5> tuple;
	resize(tuple, length(indexShape(qgramIndex)));
	appendValue(packages, 0);
	for (unsigned i = 0; i < 5; ++i)
	{
		tuple[0] = i;
		for (unsigned j = 0; j < 5; ++j)
		{
			tuple[1] = j;
			for (unsigned k = 0; k < 5; ++k)
			{
				tuple[2] = k;
				size_t border = indexDir(qgramIndex)[hashUpper(indexShape(qgramIndex), begin(tuple,Standard()), 3)];
				if (border != back(packages))
					appendValue(packages, border);
			}
		}
	}
	clear(indexDir(qgramIndex));

	unsigned finished = 0;
	bool inTerm = isatty(fileno(stdout));

	std::cout << "Suffix tree traversal ............. ";
	if (inTerm)	std::cout << "  0%";
	std::cout << std::flush;
#pragma omp parallel for
	for (int i = 1; i < (int)length(packages) - 1; ++i)
	{
		if (packages[i-1] == packages[i]) continue;
		TFionaIndex myIndex(store.readSeqStore);
		indexSA(myIndex) = infix(indexSA(qgramIndex), packages[i-1], packages[i]);
		cargo(myIndex).replen_min = options.fromLevel;
		cargo(myIndex).replen_max = options.toLevel;

		TConstrainedIterator myConstrainedIterator(myIndex);
		traverseAndSearchCorrections(myConstrainedIterator, store, corrections, expectedTheoretical, options, alg);
		mmapAdvise(indexSA(qgramIndex), MMAP_DONTNEED, packages[i-1], packages[i]);
#pragma omp critical
		{
			++finished;
			if (inTerm)
				std::cout << "\b\b\b\b" << std::setw(3) << (100 * finished) / (length(packages) - 1) << '%' << std::flush;
		}
	}
	if (inTerm)
		std::cout << "\b\b\b\b";
	std::cout << "done. (" << SEQAN_PROTIMEDIFF(search) << " seconds)" << std::endl;

#endif

  std::ofstream out("id");
	unsigned totalCorrections = 0;
	for (unsigned i = 0; i < length(corrections); ++i)
	{
		if (corrections[i].overlap == 0) continue;
		unsigned readId = i;
		if (readId >= readCount)
			readId -= readCount;
		out << readId << " 0" << std::endl;
		++totalCorrections;
	}
	applyReadErrorCorrections(store, corrections, options);
	std::cout << "Total corrected reads number is "<< totalCorrections << std::endl;

	// remove reverse complements
	resize(store.readSeqStore, readCount);
	return totalCorrections;
}

int main(int argc, const char* argv[]) 
{
	CommandLineParser parser;
  std::string rev = "$Revision  0000$";
	addVersionLine(parser, "Fiona version 1.1 2010912 [" + rev.substr(11, 4) + "]");

	FionaOptions options;
	enum { Poisson = 0, Expected = 1, PoissonSens = 2, Count = 3} method = Poisson;

	//////////////////////////////////////////////////////////////////////////////
	// Define options
	addTitleLine(parser, "**********************************************");
	addTitleLine(parser, "***   Fiona - SHREC's Significant Other    ***");
	addTitleLine(parser, "*** (c) Copyright 2010 by Gingerbread Men  ***");
	addTitleLine(parser, "**********************************************");
	addUsageLine(parser, "[OPTION]... <INPUT READS> <CORRECTED READS>");
	addOption(parser, CommandLineOption("f",  "expected",          "use expected value correction with given strictness cutoff", OptionType::Double | OptionType::Label));
	addOption(parser, CommandLineOption("c",  "count",			   "use fixed count correction cutoff", OptionType::Double | OptionType::Label)) ;
	addOption(parser, CommandLineOption("e",  "error-rate" ,       "Give expected error-rate of reads, activate sensitivity mode", OptionType::Double | OptionType::Label ));
	addOption(parser, CommandLineOption("p",  "pvalue",            "set p-value for error detection, (default is false discovery rate,\n\t switch to false negative rate in sensitivity mode)", OptionType::Double | OptionType::Label));
	addOption(parser, CommandLineOption("l",  "levels", 2,         "set lower and upper bound for suffix tree DFS", OptionType::Int | OptionType::Label));
	addOption(parser, CommandLineOption("g",  "genome-length",     "set the length of the underlying genome", OptionType::Int | OptionType::Label));
	addOption(parser, CommandLineOption("m",  "mismatches",        "set the number of accepted mismatches per read", OptionType::Int | OptionType::Label, options.acceptedMismatches));
	addOption(parser, CommandLineOption("i",  "iterations",        "set the number of error correction iterations (0 = auto mode)", OptionType::Int | OptionType::Label, options.cycles));
#ifdef FIONA_ALLOWINDELS
	addOption(parser, CommandLineOption("id", "indel-length",      "set the maximum length of an indel", OptionType::Int | OptionType::Label, options.maxIndelLength));
#endif
#ifdef _OPENMP
	addOption(parser, CommandLineOption("nt", "num-threads",       "set the number of threads used", OptionType::Int | OptionType::Label));
#endif
	addOption(parser, CommandLineOption("dr", "debug-read",        "dump information for a read given its id", OptionType::Int | OptionType::Label));
	addOption(parser, CommandLineOption("cr", "corr-read",         "dump information for a correcting read", OptionType::Int | OptionType::Label));
	requiredArguments(parser, 2);

	bool stop = !parse(parser, argc, argv, std::cerr);
	if (stop) return 0;

	if (isSetLong(parser, "expected"))
	{
		if ((isSetLong(parser, "pvalue") || isSetLong(parser, "count")) && (stop = true))
		{
			std::cout << "You have already chosen the p-value or count correction method." << std::endl;
			std::cout << "Please use only one of the options -f or -p, for more details see README.txt" << std::endl;
		}
		else
		{
			method = Expected;
			getOptionValueLong(parser, "expected", options.strictness);
		}
	}
	else if (isSetLong(parser, "count")){
		method = Count; 
		getOptionValueLong(parser, "count", options.strictness) ;
	}
	else if (isSetLong(parser, "pvalue"))
	{
		method = Poisson;
		getOptionValueLong(parser, "pvalue", options.strictness);
	}

	if (isSetLong(parser, "error-rate")){
		method = PoissonSens ;
		getOptionValueLong(parser, "error-rate", options.errorrate) ;
	}
	
	if (isSetLong(parser, "levels"))
	{
		getOptionValueLong(parser, "levels", 0, options.fromLevel);
		getOptionValueLong(parser, "levels", 1, options.toLevel);
		if ((options.fromLevel > options.toLevel) && (stop = true))
			std::cerr << "Lower bound must be less or equal upper bound." << std::endl;
	}

	getOptionValueLong(parser, "genome-length", options.genomeLength);
	getOptionValueLong(parser, "mismatches", options.acceptedMismatches);
	getOptionValueLong(parser, "iterations", options.cycles);
#ifdef FIONA_ALLOWINDELS
	getOptionValueLong(parser, "indel-length", options.maxIndelLength);
#endif

#ifdef _OPENMP
	if ((options.genomeLength < 2) && (stop = true))
		std::cerr << "A genome length must be given. Length estimation does not work in parallel mode." << std::endl;
	if (isSetLong(parser, "num-threads"))
	{
		unsigned numThreads = 1;
		getOptionValueLong(parser, "num-threads", numThreads);
		omp_set_num_threads(numThreads);
	}
#endif

	if (isSetLong(parser, "debug-read"))
	{
		getOptionValueLong(parser, "debug-read", options.debugRead);
	}
	if (isSetLong(parser, "corr-read"))
	{
		getOptionValueLong(parser, "corr-read", options.corrRead);
	}


	// something went wrong
	if (stop)
	{
		std::cerr << "Exiting ..." << std::endl;
		return 1;
	}

	SEQAN_PROTIMESTART(correction);

	// load original set of reads
	TFionaFragStore store;
	if (!loadReads(store, getArgumentValue(parser, 0)))
	{
		std::cerr << "Failed to open reads file " << getArgumentValue(parser, 0) << std::endl;
		std::cerr << "Exiting ..." << std::endl;
		return 1;
	} else
		std::cout << "Loaded " << length(store.readSeqStore) << " reads." << std::endl;

	// initialise the top and down level by using the log4 from the total number of reads
	if (options.fromLevel == 0)
	{
		int logRation = static_cast<int>(log10(static_cast<double>(length(store.readSeqStore))) / log10(4.0));
		options.fromLevel = logRation + 2;
		options.toLevel   = options.fromLevel + 10;
		std::cout << "The estimated top level is " << options.fromLevel << " and the down level is " << options.toLevel << std::endl;
	}

	unsigned firstCycleCorrections = 0;
	bool autoCycles = (options.cycles == 0);
	if (autoCycles) options.cycles = 20;
	
	unsigned cycle;
	for (cycle = 1; cycle <= options.cycles; ++cycle)
	{
		std::cout << std::endl << "Cycle "<< cycle;
		if (!autoCycles) std::cout << " of " << options.cycles;
		std::cout << std::endl;
		unsigned numCorrected;
		if (method == Poisson)
			/*use of p-value like a limit*/
			numCorrected = correctReads(store, options, FionaPoisson());
		else if (method == Count)
			numCorrected = correctReads(store, options, FionaCount()) ;
		else if (method == PoissonSens)
			numCorrected = correctReads(store, options, FionaPoissonSens()) ;
		else {
			/*use an expected value for a certain level*/
			numCorrected = correctReads(store, options, FionaExpected());
		}


//		if (options.acceptedMismatches > 0) --options.acceptedMismatches;

		// TODO maybe to stop if there is not reads corrected in the cycle before
		// if so after each iteration must save the ID for the reads which are corrected
		// thus we can also show the total number of reads that are corrected at the final stage
		if (cycle == 1)
			firstCycleCorrections = numCorrected;
		else
		{
			if (autoCycles && numCorrected <= firstCycleCorrections / 5)
			{
				++cycle;
				break;
			}
		}
	}

	// write in file all input reads with the corrected one
	std::ofstream out(toCString(getArgumentValue(parser, 1)));
	int numCorrected = 0;
	for (unsigned i = 0; i < length(store.readNameStore); ++i)
	{
		// to give the number of reads corrected for several iteration
		if (strContains(toCString(store.readNameStore[i]), "corrected"))
			++numCorrected;

		out << '>' << store.readNameStore[i] << std::endl;
		out << store.readSeqStore[i] << std::endl;
	}

	if (options.cycles > 1)
		std::cout << "Total number reads corrected for " << cycle-1 << " cycles is " << numCorrected << std::endl;

//	struct rusage usage;
//	getrusage(RUSAGE_SELF, &usage);
	std::cout << std::endl;
	std::cout << "Time required for execution: " << SEQAN_PROTIMEDIFF(correction) << " seconds." << std::endl;
//	std::cout << "Peak resident memory usage:  " << usage.ru_maxrss / (1024*1024) << " Mb." << std::endl;

	return 0;
}
