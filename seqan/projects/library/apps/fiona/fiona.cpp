//#define SEQAN_DEBUG_INDEX
//#define SEQAN_DEBUG
//#define SEQAN_VERBOSE
//#define SEQAN_VVERBOSE
#define SEQAN_PROFILE					// enable time measuring
//#define FIONA_ALLOWINDELS

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>

#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <seqan/misc/misc_cmdparser.h>

// TODO(holtgrew): This raises a warning with Boost 1.42. Deactivate warnings, activate again afterwards. The correct #pragma has to be used for each supported compiler.
//#include <boost/math/distributions/normal.hpp>

//#define MEDIAN

using namespace std;
using namespace seqan;

struct FionaPoisson_;
struct FionaExpected_;

typedef Tag<FionaPoisson_> const FionaPoisson;
typedef Tag<FionaExpected_> const FionaExpected;

struct FionaCorrectedError
{
	unsigned int   correctReadId;
	unsigned int   occurrences;
	unsigned short errorPos;
	unsigned short correctPos;
	unsigned short overlap;
	unsigned char  mismatches;
	signed   char  indelLength;		// 0..mismatch, <0..deletion, >0..insertion
	bool reverse;
};

struct FionaOptions
{
	__int64 genomeLength;
	double strictness;
	unsigned acceptedMismatches;
	enum { Poisson = 0, Expected = 1} method;
	int maxIndelLength;
	int fromLevel;
	int toLevel;
	unsigned cycles;

	FionaOptions()
	{
		method = Poisson;
		genomeLength = 0;
		strictness = 0.0001;
		acceptedMismatches = 1;
		maxIndelLength = 0;
		cycles = 3;
		fromLevel = 0;
		toLevel = 0;
	}
};

#ifdef FIONA_ALLOWINDELS

	struct FionaStoreConfig:
		public FragmentStoreConfig<>
	{
		typedef Owner<>	TReadSeqStoreSpec;
	};

	typedef FragmentStore<void, FionaStoreConfig> TFionaFragStore;

#else

	typedef FragmentStore<> TFionaFragStore;

#endif

typedef Index<TFionaFragStore::TReadSeqStore, Index_Wotd<> > TFionaIndex; 

namespace seqan  {

/*restriction for certain levels - between max and min, also table with frequency may be to use eventually TODO*/
	struct FionaNodeConstraints { 
		 int replen_max; 
		 int replen_min;
		 map<unsigned,double> frequency;
		 bool _cachedPred; 
	}; 

	template <>
	struct Cargo<TFionaIndex> { 
		typedef FionaNodeConstraints Type; 
	}; 

	/*TODO THIS FONCTION CAN BE CHANGED FOR THE FREQUENCY - here just one experience*/
	/*by the use also the frequency for A,T,G,C*/
	/*higher frequency - high level as min in which we will begin the searching*/

	/*hide the node between certain level*/
	template <>
	bool nodePredicate(Iter<TFionaIndex, VSTree<TopDown<ParentLinks<Postorder> > > > &it)
	{
		FionaNodeConstraints &cons = cargo(container(it));
		int valueT = parentRepLength(it);
		
		/*necessary if we use the frequency*/
		int level_min = cons.replen_min;
		
		/*'A' = 0, 'C' = 1, 'G' = 2, 'T' = 3*/
		/*
		String<Dna> d = representative(it);
		typedef Iterator<String<Dna> >::Type TIterator;
		for (TIterator ic = begin(d); ic != end(d); ++ic){
			if(value(ic)==0 || value(ic)==3){
				if((cons.frequency[0]+cons.frequency[3])>(0.1+cons.frequency[1]+cons.frequency[2])){
					level_min = cons.replen_min + 3;
					cons.replen_max = cons.replen_max;
				}
			}else{
				if((cons.frequency[1] + cons.frequency[2])>(0.1 + cons.frequency[1] + cons.frequency[2]) ){level_min = cons.replen_min + 3;	
					level_min = cons.replen_min + 3;
				}
			}
			break;
		}*/

		/*TODO may utilise >=*/
		if (valueT > level_min && valueT < cons.replen_max) return cons._cachedPred = true;
		return cons._cachedPred = false;
	}
}




/*TODO cumulative poisson distribution determinated just once for a level, not at each node like now*/

/*matching string*/
inline bool strContains(string const & inputStr, string const & searchStr)
{
	return inputStr.find(searchStr) != string::npos;
}

/*save in table, information about the erroneous reads and the respective correct one */
Pair<unsigned,Pair<vector<unsigned>, Dna> > 
dataErroneousNodes(
	int bestReadId, 
	int positionError, 
	int commonPrefix, 
	int bestOverlap, 
    int bestMismatches, 
	int bestOccurrences, 
	Dna nucleotide)
{													 
	Pair<unsigned,Pair<vector<unsigned>, Dna> > pairCorrect;	

	/*data for the erroneous node and correct one*/
	/*position error, length parent(path-label), nb match after error position, nb mismatch, nb of occurrences */ 	
	vector<unsigned> data;	
									
	data.push_back(positionError);
	data.push_back(commonPrefix);
	data.push_back(bestOverlap);
	data.push_back(bestMismatches);
	data.push_back(bestOccurrences);
 
	Pair<vector<unsigned>, Dna> dataCorrection;
	dataCorrection.i1 = data;
	dataCorrection.i2 = nucleotide; 
	pairCorrect.i1 =  bestReadId;
	pairCorrect.i2 = dataCorrection;
	return pairCorrect;
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
void expectedValueTheorical(TExpectedValues & expected, TReadSet const & readSet, TGenomeSize const genomeLength)
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
			fill(readLengthHistogram, readLength + 1, 0);
		++readLengthHistogram[readLength];
	}

	// a = read_length - suffix_length + 1
	clear(expected);
	fill(expected, length(readLengthHistogram), 0.0);
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

	map<unsigned, unsigned> vectorOccurences;

	goBegin(iter);
	for (; !atEnd(iter); ++iter)
	{
		unsigned numOccs = countOccurrences(iter);
		++vectorOccurences[numOccs];
		totalOccs += numOccs;
	}
	
	mediumTotalOccs = totalOccs / 2.0;
	
	map<unsigned,unsigned>::iterator iterMap;
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
	cout << endl;
	cout << "indelLen: " << (int)correction.indelLength << endl;
	cout << "reversed: " << correction.reverse << endl;
	cout << "corrId:   " << correction.correctReadId << endl;
	cout << "errorId:  " << errorReadId << endl;
	cout << "overlap:  " << correction.overlap << endl;
	cout << "errors:   " << (unsigned)correction.mismatches << endl;
	cout << "occur:    " << correction.occurrences << endl;
	cout << "corrPos:  " << correction.correctPos << endl;
	cout << "errPos:   " << correction.errorPos << endl;
	cout << "corrRead: ";
	for (unsigned i = 0; i < correction.errorPos; ++i)
		cout << ' ';
	cout << store.readSeqStore[correction.correctReadId] << endl;
	cout << "errRead:  ";
	for (unsigned i = 0; i < correction.correctPos; ++i)
		cout << ' ';
	cout << store.readSeqStore[errorReadId] << endl;
}

/*change the erroneous nucleotide in all reads identify with errors*/
template <typename TFragmentStore, typename TCorrections>
void applyReadErrorCorrections(
	TFragmentStore &store,
	TCorrections const &corrections)
{	
	typedef typename Iterator<TCorrections, Standard>::Type TIter;
	
	TIter it = begin(corrections, Standard());
	TIter itEnd = end(corrections, Standard());
	
	for (unsigned readId = 0; it != itEnd; ++it, ++readId)
		if ((*it).overlap != 0)
		{
			ostringstream m;
			if (strContains(toCString(store.readNameStore[readId]), "corrected"))
				m << "," << (*it).errorPos;
			else
				m << " corrected: " << (*it).errorPos;

#ifdef SEQAN_VERBOSE
			_dumpCorrection(store, *it, readId);
#endif
			append(store.readNameStore[readId], m.str());
			
			if (store.readSeqStore[readId] == "AGGGGGAAATATGATCGCGTATGCGAGAGTAGTGCC")
				cout<<"OK"<<endl;
			
			if ((*it).indelLength == 0)
				store.readSeqStore[readId][(*it).errorPos] = store.readSeqStore[(*it).correctReadId][(*it).correctPos];
#ifdef FIONA_ALLOWINDELS
			else if ((*it).indelLength > 0)
				erase(store.readSeqStore[readId], (*it).errorPos, (*it).errorPos + (*it).indelLength);
			else
				insert(store.readSeqStore[readId], (*it).errorPos, infix(store.readSeqStore[(*it).correctReadId], (*it).correctPos, (*it).correctPos + -(*it).indelLength));
#endif		
#ifdef SEQAN_VERBOSE
			cout << "corrected:";
			for (unsigned i = 0; i < (*it).correctPos; ++i)
				cout << ' ';
			cout << store.readSeqStore[readId] << endl;
#endif
		}
}

template <typename TObserved, typename TExpected, typename TStrictness>
inline bool potentiallyErroneousNode(
	TObserved observed, 
	TExpected expected,
	TStrictness strictness,
	FionaPoisson const)
{
	// compare the cumulative poisson distribution with the p-value (strictness)
	return pValue(observed, expected) <= strictness;
}

template <typename TObserved, typename TExpected, typename TStrictness>
inline bool potentiallyErroneousNode(
	TObserved observed, 
	TExpected expected,
	TStrictness,
	FionaExpected const)
{
	// compare the weight for a node with its expected value (dependent on the length)
	return observed < expected;
}

/*detect and repare the reads with errors*/
template <typename TTreeIterator, typename TFragmentStore, typename TAlgorithm>
void traverseAndSearchCorrections(
	TTreeIterator iter,
	TFragmentStore &store,
	FionaOptions const &options,
	Tag<TAlgorithm> const alg)
{
	typedef typename Container<TTreeIterator>::Type TFionaIndex;
	typedef typename Fibre<TFionaIndex, Fibre_Text>::Type TReadSet;
	typedef typename Fibre<TFionaIndex, Fibre_SA>::Type TSA;
	typedef typename Infix<TSA const>::Type TOccs;
	typedef typename Iterator<TOccs, Standard>::Type TOccsIterator;
	typedef typename Value<TReadSet>::Type TRead;
	typedef typename Iterator<TRead, Standard>::Type TReadIterator;
	typedef FionaCorrectedError TCorrected;
	
	if (TYPECMP<TAlgorithm, FionaExpected_>::VALUE)
		cout << endl << "Method with expected value for each level" << endl;
	else
		cout << endl << "Method with p-value and Poisson distribution" << endl;
	cout << "Searching... " << endl;
	
	/*table with the theoricals values*/
	String<double> expectedTheorical;
	expectedValueTheorical(expectedTheorical, store.readSeqStore, options.genomeLength);
	unsigned readCount = length(store.readSeqStore) / 2;
		
	if (TYPECMP<TAlgorithm, FionaExpected_>::VALUE)
	{
		String<double> sd;
		standardDeviation(sd, store.readSeqStore, options.genomeLength);

		/*The strictness value allows to estimate the confidential intervall*/
		for (unsigned i = 0; i < length(expectedTheorical); ++i)
		{
			double expectedTemporary = expectedTheorical[i] - options.strictness * sd[i];
			
			/*If the connfidential intervall take value less than 1 ??? not sure for that*/
			/*if(expectedTemporary < 1){
				expectedTheorical[i] = 1.1;
			}else{*/
				expectedTheorical[i] = expectedTemporary;
			//}
			//if fixed
			//expectedTheorical[i]=5;
		}
	}
	
	String<TCorrected> corrections;
	TCorrected zeroCorr = {0, 0, 0, 0, 0, 0, 0, 0};
	fill(corrections, readCount, zeroCorr);

	String<TOccs> correctCandidates;

	SEQAN_PROTIMESTART(search);
	for (goBegin(iter); !atEnd(iter); goNext(iter))		// do a DFS walk
	{
		int commonPrefix = parentRepLength(iter);			// length of parent label

		SEQAN_ASSERT_LT(commonPrefix + 1, length(expectedTheorical));
		if (!potentiallyErroneousNode(countOccurrences(iter), expectedTheorical[commonPrefix+1], options.strictness, alg))
			continue;
		
		//
		//	get the id and position (suffix begin) for suspected nodes
		//	for which we can find a more optimal correction
		//
		TOccs errorCandidates = getOccurrences(iter);
					
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
				!potentiallyErroneousNode(countOccurrences(iterSibling), expectedTheorical[commonPrefix + 1], options.strictness, alg))
			{
				// save the id and position(where the suffix begin in the reads) in the table of IDs correct
				// also the number of occurrences
				appendValue(correctCandidates, getOccurrences(iterSibling));
			}
		} while (goRight(iterSibling));
		
		// continue if we haven't found any correct read
		if (empty(correctCandidates)) continue;
		
		// make the comarison between the substrings(the suffix after the position of error
		TOccsIterator errorRead = begin(errorCandidates, Standard());
		TOccsIterator errorReadEnd = end(errorCandidates, Standard());
		for (; errorRead != errorReadEnd; ++errorRead)
		{
			unsigned errorReadId = (*errorRead).i1;
			if (errorReadId >= readCount) errorReadId -= readCount;
//			bool dbg = store.readSeqStore[errorReadId] == "AGGGGGAAATATGATCGCGTATGCGAGAGTAGTGCC";

			// is already identify as erroneus
			TCorrected &corr = corrections[errorReadId];

			// Branch and Bound
			//
			// if already detect as erroneous(optimal)
			// the total length minus the position of error, and minus 1 mismatch
			
			// the position where the error is
			unsigned positionError = (*errorRead).i2 + commonPrefix;
			int maxOverlapPossible = commonPrefix + length(store.readSeqStore[errorReadId]) - positionError - 1;
//			if (dbg) std::cout<<"errorReadId:"<<errorReadId<<std::endl;
//			if (dbg) std::cout<<"maxOverlapPossible:"<<maxOverlapPossible<<std::endl;
//			if (dbg) std::cout<<"corr.overlap:"<<corr.overlap<<std::endl;
//			if (dbg) 
//			std::cout<<"positionError:"<<positionError<<std::endl;
			if (corr.overlap > maxOverlapPossible) continue;

			unsigned bestReadId = 0;
			unsigned bestOccurrences = 0;
			unsigned bestCorrectPos = 0;
			unsigned bestOverlap = 0;
			unsigned bestMismatches = 0;
			int bestIndelLength = 0;
			
			for (unsigned c = 0; c < length(correctCandidates); ++c)
			{
				TOccsIterator corrRead = begin(correctCandidates[c], Standard());
				TOccsIterator corrReadEnd = end(correctCandidates[c], Standard());
				for (; corrRead != corrReadEnd; ++corrRead)
				{
					/*if not already correct by the same reads*/
					/*if(mapID!=allErrors.end()){
					
						if(idReadsCorrect[r].i1.i1==correctAndData.i1){	
							continue;
						}
					}*/

					/*the position in the read until which there is the same prefix*/
					unsigned positionCorrect = (*corrRead).i2 + commonPrefix;
										
					/*search the position detect as erroneous at the level of the read*/
					TReadIterator itEEnd = end(store.readSeqStore[(*errorRead).i1], Standard());
					TReadIterator itCEnd = end(store.readSeqStore[(*corrRead).i1], Standard());
					
					for (int indel = -options.maxIndelLength; indel <= options.maxIndelLength; ++indel)
					{
						TReadIterator itE = begin(store.readSeqStore[(*errorRead).i1], Standard());
						TReadIterator itC = begin(store.readSeqStore[(*corrRead).i1], Standard());
						
						if (indel == 0)
						{
							// mismatch
							itE += positionError + 1;
							itC += positionCorrect + 1;
						}
						else if (indel > 0)
						{
							// insert in erroneous read
							itE += positionError + indel;
							itC += positionCorrect;
						} 
						else
						{
							// deletion in erroneous read
							itE += positionError;
							itC += positionCorrect - indel;
						}
						
						unsigned overlap = commonPrefix;
						unsigned mismatches = 0;
						
						for (; itE < itEEnd && itC < itCEnd; ++itE, ++itC, ++overlap)
							if (*itE != *itC)
							{
								if (++mismatches > options.acceptedMismatches)
									break;
							}

						if (itE == itEEnd || itC == itCEnd)
						{
							if (overlap < bestOverlap)
								continue;
							
							mismatches += (indel > 0)? indel: -indel;

							if (overlap == bestOverlap)
							{
								if (mismatches > bestMismatches)
									continue;
								
								if (mismatches == bestMismatches && length(correctCandidates[c]) < bestOccurrences)
									continue;
							}

							bestReadId = (*corrRead).i1;
							bestOccurrences = length(correctCandidates[c]);
							bestCorrectPos = positionCorrect;
							bestOverlap = overlap;
							bestMismatches = mismatches;
							bestIndelLength = indel;
//			if (dbg) std::cout<<"bestReadId:"<<bestReadId<<std::endl;
//			if (dbg) std::cout<<"bestOccurrences:"<<bestOccurrences<<std::endl;
//			if (dbg) std::cout<<"bestCorrectPos:"<<bestCorrectPos<<std::endl;
//			if (dbg) std::cout<<"bestOverlap:"<<bestOverlap<<std::endl;
//			if (dbg)
//			 std::cout<<"bestIndelLength:"<<bestIndelLength<<std::endl;

//					cout <<"corr:"<<store.readSeqStore[(*corrRead).i1]<<endl;
//					cout <<"err: "<<store.readSeqStore[(*errorRead).i1]<<endl;
						}
					}
				}
			}
	
			if (bestOverlap != 0)
			{
				bool rev = false;
			
				if ((*errorRead).i1 >= readCount)
				{
//					cout <<"0corr:"<<store.readSeqStore[bestReadId]<<endl;
//					cout <<"0err: "<<store.readSeqStore[(*errorRead).i1]<<endl;
					// switch to reverse-complements
					if (bestReadId < readCount)
						bestReadId += readCount;
					else
						bestReadId -= readCount;
					
					// mirror positions
					positionError = length(store.readSeqStore[errorReadId]) - positionError - 1;
					bestCorrectPos = length(store.readSeqStore[bestReadId]) - bestCorrectPos - 1;

					if (bestIndelLength > 0)
						positionError -= bestIndelLength;
					else
						bestCorrectPos -= -bestIndelLength;
					rev = true;
		TTreeIterator iterSibling(iter);
		goUp(iterSibling);
//					std::cout << "common:"<<representative(iterSibling)<<endl;
//					cout <<"corr:"<<store.readSeqStore[bestReadId]<<endl;
//					cout <<"err: "<<store.readSeqStore[errorReadId]<<endl;
				}
				
				/*This is maybe not necessary, because we search to maximisate*/
				/*when the position or the nucleotide differ*/
				//if(number.i1[0]!=positionError || number.i2 != nucleotide){
				
				/*find longer commun part with certain nb accepted mismatch*/
				if (corr.overlap > bestOverlap)
					continue;
				
				// if they are equal
				if (corr.overlap == bestOverlap)
				{
					// looking for longer matching part without mismatch
					if (corr.mismatches < bestMismatches)
						continue;
						
					// equal number of match and mismatch, but greater nb occurrences
					if (corr.mismatches == bestMismatches && corr.occurrences > bestOccurrences)
						continue;
				}

				// update best correction
				corr.correctReadId = bestReadId;
				corr.reverse = rev;
				corr.occurrences = bestOccurrences;
				corr.errorPos = positionError;
				corr.correctPos = bestCorrectPos;
				corr.overlap = bestOverlap;
				corr.mismatches = bestMismatches;
				corr.indelLength = bestIndelLength;
//			if (dbg) _dumpCorrection(store, corr, errorReadId);			
			}
		}
	}
	cout << "Time for searching between given levels: "<< SEQAN_PROTIMEDIFF(search) << " seconds." << endl;

	ofstream out("id");
	unsigned totalCorrections = 0;
	for (unsigned i = 0; i < length(corrections); ++i)
	{
		if (corrections[i].overlap == 0) continue;
		unsigned readId = i;
		if (readId >= readCount)
			readId -= readCount;
		out << readId << " 0" << endl;
		++totalCorrections;
	}
	
	applyReadErrorCorrections(store, corrections);
	cout << "Total corrected reads number is "<< totalCorrections << endl;
}


/*GC-content*/
/*fonction which allow to determine the frequency for each nucleotide*/
template < typename TFionaIndex, class TSpec >
map<unsigned,double>
determineFrequency(Iter< TFionaIndex, VSTree<TSpec> > iter)
{
   /*calculate the frequency for each nucleotide*/
   /*'A' = 0, 'C' = 1, 'G' = 2, 'T' = 3*/
   map<unsigned, double> frequency;

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
template <typename TFragmentStore>
void correctReads(
	TFragmentStore & store,
	FionaOptions & options)
{		   
	/*iterator with restrictions*/
	typedef Iterator<TFionaIndex, TopDown<ParentLinks<Postorder> > >::Type TConstrainedIterator; 

	cout << "Construct suffix array" << endl;

	// append their reverse complements
	unsigned readCount = length(store.readSeqStore);
	Dna5String tmp;
	if (options.genomeLength != 1)
		for (unsigned i = 0; i < readCount; ++i)
		{
			tmp = store.readSeqStore[i];
			reverseComplementInPlace(tmp);
			appendValue(store.readSeqStore, tmp);
		}

	SEQAN_PROTIMESTART(construct);

	// construct suffix array of the set of reads
	TFionaIndex myIndex(store.readSeqStore);
	TConstrainedIterator myConstrainedIterator(myIndex); 

	/*calculate the frequency for each nucleotide, didn't use for the moment*/
	/*'A' = 0, 'C' = 1, 'G' = 2, 'T' = 3*/
	map<unsigned,double> frequency;
	frequency = determineFrequency(myConstrainedIterator);
	
	cout << "Time required for suffix array construction : " << SEQAN_PROTIMEDIFF(construct) << " seconds." << endl;

	/*restrictions just for estimate the genome length if there is no data*/

#ifdef MEDIAN
	unsigned level = fromLevel;
	ofstream out("medianLevels.txt"); 
	ofstream median("medianForEachLevel.txt");
#endif

	if (options.genomeLength == 1)
	{
		cout << "Generating Hugues' stats file." << endl;
		ofstream stats("stats.txt");
		Iterator<TFionaIndex, TopDown<ParentLinks<Preorder> > >::Type it(myIndex);
		goBegin(it);
		CharString tmp;
		String<int> freq;
		stats << "branch\tlength\ttreeDep\tletter\treadPos\tfreq" << endl;
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
						fill(freq, posInRead + 1, 0);
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
						stats << '\t' << k << '\t' << freq[k] << '\t';// << representative(it) << '\t' << tmp << endl;
						stats << endl;					}
			}
			++it;
		}
		stats.close();
		cout << "Done." << endl;
		exit(0);
	}
			
	if (options.genomeLength == 0)
	{
#ifdef MEDIAN		
		for (; level < toLevel; ++level)
		{
			//int logRation = (log10(length(setReads)/2))/(log10(4));
			//int l = logRation + 1;
			//cout << l << endl;
			
			cargo(myIndex).replen_min = level;
			cargo(myIndex).replen_max = level+2; 	
			cargo(myIndex).frequency = frequency;
			double numOccs = 0.0; 
			
			median << level << " " << medianLevel(myConstrainedIterator) << endl;
			goBegin(myConstrainedIterator);
			while (!atEnd(myConstrainedIterator))
			{
				if (parentRepLength(myConstrainedIterator) > level)
				{
					numOccs = countOccurrences(myConstrainedIterator);
					out << level << " " << numOccs << endl; 
				}
				++myConstrainedIterator;
			}
		}
#else
		//int logRation = log10(readCount) / log10(4);
		//int l = logRation + 1;
		//cout << l << endl;
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
		options.genomeLength = readCount * a / expectedValueGivenLevel;
		cout << "The estimated genome length is " << options.genomeLength << endl;
#endif
	}
	
	vector<unsigned> idCorrected;
	
	/*restrictions for the searching levels*/
	cargo(myIndex).replen_min = options.fromLevel;
	cargo(myIndex).replen_max = options.toLevel; 	
	cargo(myIndex).frequency = frequency;

	/*the core of the correction method*/
	if (options.method == options.Poisson)
		/*use of p-value like a limit*/
		traverseAndSearchCorrections(myConstrainedIterator, store, options, FionaPoisson());
	else
		/*use an expected value for a certain level*/
		traverseAndSearchCorrections(myConstrainedIterator, store, options, FionaExpected());
	
	// remove reverse complements
	resize(store.readSeqStore, readCount);
}

int main(int argc, const char* argv[]) 
{
	CommandLineParser parser;
	string rev = "$Revision  0000$";
	addVersionLine(parser, "Fiona version 1.1 2010912 [" + rev.substr(11, 4) + "]");

	FionaOptions options;

	//////////////////////////////////////////////////////////////////////////////
	// Define options
	addTitleLine(parser, "**********************************************");
	addTitleLine(parser, "***   Fiona - SHREC's Significant Other    ***");
	addTitleLine(parser, "*** (c) Copyright 2010 by Gingerbread Man  ***");
	addTitleLine(parser, "**********************************************");
	addUsageLine(parser, "[OPTION]... <INPUT READS> <CORRECTED READS>");
	addOption(parser, CommandLineOption("f",  "expected",          "use expected value correction with given strictness", OptionType::Double | OptionType::Label));
	addOption(parser, CommandLineOption("p",  "pvalue",            "use p-value correction with given strictness", OptionType::Double | OptionType::Label));
	addOption(parser, CommandLineOption("l",  "levels", 2,         "set lower and upper bound for suffix tree DFS", OptionType::Int | OptionType::Label));
	addOption(parser, CommandLineOption("g",  "genome-length",     "set the length of the underlying genome", OptionType::Int | OptionType::Label));
	addOption(parser, CommandLineOption("m",  "mismatches",        "set the number of accepted mismatches per read", OptionType::Int | OptionType::Label, options.acceptedMismatches));
	addOption(parser, CommandLineOption("i",  "iterations",        "set the number of error correction iterations", OptionType::Int | OptionType::Label, options.cycles));
#ifdef FIONA_ALLOWINDELS
	addOption(parser, CommandLineOption("id", "indel-length",      "set the maximum length of an indel", OptionType::Int | OptionType::Label, options.maxIndelLength));
#endif
	requiredArguments(parser, 2);

	bool stop = !parse(parser, argc, argv, cerr);
	if (stop) return 0;

	if (isSetLong(parser, "expected"))
	{
		if (isSetLong(parser, "pvalue") && (stop = true))
		{
			cout << "You have already chosen the p-value correction method." << endl;
			cout << "Please use only one of the options -f or -p, for more details see README.txt" << endl;
		}
		else
		{
			options.method = options.Expected;
			getOptionValueLong(parser, "expected", options.strictness);
		}
	} 
	else if (isSetLong(parser, "pvalue"))
	{
		options.method = options.Poisson;
		getOptionValueLong(parser, "pvalue", options.strictness);
	}

	if (isSetLong(parser, "levels"))
	{
		getOptionValueLong(parser, "levels", 0, options.fromLevel);
		getOptionValueLong(parser, "levels", 1, options.toLevel);
		if ((options.fromLevel > options.toLevel) && (stop = true))
			cerr << "Lower bound must be less or equal upper bound." << endl;
	}

	getOptionValueLong(parser, "genome-length", options.genomeLength);
	getOptionValueLong(parser, "mismatches", options.acceptedMismatches);
	getOptionValueLong(parser, "iterations", options.cycles);
#ifdef FIONA_ALLOWINDELS
	getOptionValueLong(parser, "indel-length", options.maxIndelLength);
#endif
	
	// something went wrong
	if (stop)
	{
		cerr << "Exiting ..." << endl;
		return 1;
	}
	
	SEQAN_PROTIMESTART(correction);
	
	// load original set of reads
	TFionaFragStore store;
	if (!loadReads(store, getArgumentValue(parser, 0)))
	{
		cerr << "Failed to open reads file " << getArgumentValue(parser, 0) << endl;
		cerr << "Exiting ..." << endl;
		return 1;
	}
	
	// initialise the top and down level by using the log4 from the total number of reads
	if (options.fromLevel == 0)
	{
		int logRation = log10(length(store.readSeqStore)) / log10(4.0);
		options.fromLevel = logRation + 2;
		options.toLevel   = options.fromLevel + 10;
		cout << "The estimated top level is " << options.fromLevel << " and the down level is " << options.toLevel << endl;
	}
	cout << endl;

	for (unsigned cycle = 1; cycle <= options.cycles; ++cycle)
	{
		cout << "Cycle "<< cycle << " of " << options.cycles << endl;
		correctReads(store, options);
		
		if (options.acceptedMismatches > 0) --options.acceptedMismatches;			

		// TODO maybe to stop if there is not reads corrected in the cycle before
		// if so after each iteration must save the ID for the reads which are corrected
		// thus we can also show the total number of reads that are corrected at the final stage
	}
	
	// write in file all input reads with the corrected one
	ofstream out(toCString(getArgumentValue(parser, 1)));
	int numCorrected = 0;
	for (unsigned i = 0; i < length(store.readNameStore); ++i)
	{
		// to give the number of reads corrected for several iteration 
		if (strContains(toCString(store.readNameStore[i]), "corrected"))
			++numCorrected;

		out << '>' << store.readNameStore[i]<<endl;
		out << store.readSeqStore[i] << endl;
	}

	if (options.cycles > 1)
		cout << "Total number reads corrected for " << options.cycles << " cycles is " << numCorrected << endl; 

	cout << endl << "Time required for execution: " << SEQAN_PROTIMEDIFF(correction) << " seconds." << endl;

	return 0;
}

	
