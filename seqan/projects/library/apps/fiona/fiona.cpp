#define SEQAN_DEBUG_INDEX
#define SEQAN_DEBUG
#define SEQAN_VERBOSE
#define SEQAN_VVERBOSE

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>

#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/store.h>

// TODO(holtgrew): This raises a warning with Boost 1.42. Deactivate warnings, activate again afterwards. The correct #pragma has to be used for each supported compiler.
//#include <boost/math/distributions/normal.hpp>

//#define MEDIAN

using namespace std;
using namespace seqan;

struct FionaPoisson_;
struct FionaExpected_;

typedef Tag<FionaPoisson_> const FionaPoisson;
typedef Tag<FionaExpected_> const FionaExpected;

/*TODO cumulative poisson distribution determinated just once for a level, not at each node like now*/

/*matching string*/
inline bool strContains(string const & inputStr, string const & searchStr)
{
	return inputStr.find(searchStr) != string::npos;
}

/*save in table, information about the erroneous reads and the respective correct one */
Pair<unsigned,Pair<vector<unsigned>, Dna> > 
dataErroneousNodes(
	int idCorr, 
	int positionError, 
	int lenParent, 
	int best, 
    int smallMismatch, 
	int occurBest, 
	Dna nucleotide)
{													 
	Pair<unsigned,Pair<vector<unsigned>, Dna> > pairCorrect;	

	/*data for the erroneous node and correct one*/
	/*position error, length parent(path-label), nb match after error position, nb mismatch, nb of occurrences */ 	
	vector<unsigned> data;	
									
	data.push_back(positionError);
	data.push_back(lenParent);
	data.push_back(best);
	data.push_back(smallMismatch);
	data.push_back(occurBest);
 
	Pair<vector<unsigned>, Dna> dataCorrection;
	dataCorrection.i1 = data;
	dataCorrection.i2 = nucleotide; 
	pairCorrect.i1 =  idCorr;
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
double medianLevel(Iter< TIndex, VSTree<TSpec> > iter){

	double nbOccurence = 0.0; 
	double totalOccurence = 0.0;
	double sommeMedian = 0.0;
	double median = 0.0;
	double mediumTotalOccurence = 0.0;

	map<unsigned,unsigned> vectorOccurences;

	goBegin(iter);
	while (!atEnd(iter))
	{
		nbOccurence = countOccurrences(iter);
		vectorOccurences[nbOccurence]+=1.0;
		totalOccurence+=nbOccurence;
		++iter;
	}
	
	mediumTotalOccurence = totalOccurence/2.0;
	
	map<unsigned,unsigned>::iterator iterMap;
	for (iterMap = vectorOccurences.begin (); iterMap != vectorOccurences.end (); ++iterMap)
	{
		sommeMedian+= iterMap->second*iterMap->first;
		if(sommeMedian >= mediumTotalOccurence){
			median = iterMap->first;
			break;
		}
	}
	return median;	
}

template <typename TPercentage, typename TSize>
inline double probabilityOneError(TPercentage percentageErr, TSize lenPath)
{
	return 1.0 - pow(1.0 - percentageErr, lenPath);
}

/*change the erroneous nucleotide in all reads identify with errors*/
template <typename TFragmentStore>
void correctErroneousReads(
	map<unsigned,Pair<unsigned,Pair<vector<unsigned>, Dna> > > allErrors,
	TFragmentStore &store)
{	
	map<unsigned,Pair<unsigned,Pair<vector<unsigned>, Dna> > >::iterator iR;
	unsigned readCount = length(store.	readSeqStore) / 2;
	for (iR = allErrors.begin (); iR != allErrors.end (); ++iR)
	{
		Pair<unsigned,Pair<vector<unsigned>, Dna> > setR = (iR->second);
		int posErr = (setR.i2).i1[0];
		
		/*information about the read*/
		ostringstream oss;
		ostringstream m;
		unsigned readId = iR->first;
		if (readId > readCount)
			readId -= readCount;
		if (strContains(toCString(store.readNameStore[readId]), "corrected")){
			m << "," << posErr;
			oss << toCString(store.readNameStore[readId]) << m.str();
			String<char> d = oss.str();
			store.readNameStore[readId] = d;
		}else{
			m << " corrected: " << posErr;
			oss << toCString(store.readNameStore[readId])<< m.str();
			String<char> d = oss.str();
			store.readNameStore[readId] = d;
		}

		store.readSeqStore[iR->first][setR.i2.i1[0]] = setR.i2.i2;
	}
}


/*restriction for certain levels - between max and min, also table with frequency may be to use eventually TODO*/
struct TMyConstraints { 
	 unsigned int replen_max; 
	 unsigned int replen_min;
	 map<unsigned,double> frequency;
	 bool _cachedPred; 
}; 
 
namespace seqan  {

	typedef FragmentStore<> TFionaFragStore;
	typedef Index<TFionaFragStore::TReadSeqStore, Index_Wotd<> > TFionaIndex; 
	
	template <>
	struct Cargo<TFionaIndex> { 
		typedef TMyConstraints Type; 
	}; 

	/*TODO THIS FONCTION CAN BE CHANGED FOR THE FREQUENCY - here just one experience*/
	/*by the use also the frequency for A,T,G,C*/
	/*higher frequency - high level as min in which we will begin the searching*/

	/*hide the node between certain level*/
	template <>
	bool nodePredicate(Iter<TFionaIndex, VSTree<TopDown<ParentLinks<Postorder> > > > &it)
	{
		TMyConstraints &cons = cargo(container(it));
		unsigned int valueT = parentRepLength(it);
		
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
template < typename TFionaIndex, class TSpec, typename TFragmentStore, typename TAlgorithm >
void searchNode(
	Iter< TFionaIndex, VSTree<TSpec> > iter,
	TFragmentStore &store,
	double genomeLength, 
	double strictness, 
	unsigned acceptMismatch,
	Tag<TAlgorithm> const alg)
{
	typedef typename Fibre<TFionaIndex, Fibre_Text>::Type TReadSet;
	typedef typename Value<TReadSet>::Type TRead;
	
	if (TYPECMP<TAlgorithm, FionaExpected_>::VALUE)
		cout << "Method with expected value for each level" << endl;
	else
		cout << "Method with p-value and Poisson distribution" << endl;
	cout << "Searching... " << endl;
	
	/*table with the theoricals values*/
	String<double> expectedTheorical;
	expectedValueTheorical(expectedTheorical, store.readSeqStore, genomeLength);
	unsigned readCount = length(store.readSeqStore) / 2;
	
	if (TYPECMP<TAlgorithm, FionaExpected_>::VALUE)
	{
		String<double> sd;
		standardDeviation(sd, store.readSeqStore, genomeLength);

		/*The strictness value allows to estimate the confidential intervall*/
		for (unsigned i = 0; i < length(expectedTheorical); ++i)
		{
			double expectedTemporary = expectedTheorical[i] - strictness * sd[i];
			
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
	
	/*save reads with errors*/
	vector<String<Dna> > readsErrors;

	/*
		map contains : 
			key - ID of read with errors;
			value - pair of position in read and nucleotide for make change with
	*/
	map<unsigned,Pair<unsigned,Dna> > readsErrorsID;

	String<Dna> nodeError ="";

	/*
		map contains : 
			key - ID of read with errors;
			value - pair which contains:
				first value - ID of correct read
				second value - pair with first value a vectors of 5 values : 
					position error, length parent(path-label), nb match after error position, nb mismatch, nb of occurrences 
				and the second value is the correct nucleotide A,T,G or C
	*/
	map<unsigned,Pair<unsigned,Pair<vector<unsigned>, Dna> > > allErrors;

	/* length of parent label*/
	int lenParent = 0;

	/* length of label to root until courant node*/
	int lenPath = 0;

	clock_t start1, end1;

	start1 = clock();
	
	/*the bigining of the tree*/
	goBegin(iter);
	while (!atEnd(iter))
	{

		/*affectation*/
		lenParent = parentRepLength(iter);
		lenPath = repLength(iter);

		SEQAN_ASSERT_LT(lenParent + 1, length(expectedTheorical));
		if (potentiallyErroneousNode(countOccurrences(iter), expectedTheorical[lenParent+1], strictness, alg))
		{
			nodeError = representative(iter);
	
			/*
				save the id and position (suffix begin) for suspected nodes
				for which we can find a more optimal correction
			*/
			vector<Pair<unsigned,unsigned> > idReadsError;

			/*take all reads containing errors at that position in the tree*/
			for (unsigned i=0; i < length(getOccurrences(iter)); ++i){
				idReadsError.push_back(getOccurrences(iter)[i]);
			}
			/*if there is a nodes detect as erroneus*/
			if(idReadsError.size()>0){
				
				/*copy the iterator for iterate over the siblings*/
				typename Iterator< TFionaIndex, TopDown< ParentLinks<Postorder> > >::Type iterParent(iter);
			 
				/*go at the level of the parent node*/
				goUp(iterParent);
					
				/*
					potential reads for make the correction,
					because at the same level, with the same prefix
				*/
				vector<Pair<Pair<unsigned,unsigned>,unsigned> > idReadsCorrect;
		
				if(goDown(iterParent)){
					/*the length of current node*/
					int lenSiblings = lenParent+1;
					/*nb occurrence for this node*/
					int occur = countOccurrences(iterParent);
					if (representative(iterParent) != nodeError && 
						!potentiallyErroneousNode(occur, expectedTheorical[lenSiblings], strictness, alg))
					 {
							/*save the id and position(where the suffix begin in the reads) in the table of IDs correct*/
							/*also the number of occurrences*/
							for (unsigned i=0; i < length(getOccurrences(iterParent)); ++i)
							{
								Pair<Pair<unsigned,unsigned>,unsigned> dataCorrectRead;
								dataCorrectRead.i1 = getOccurrences(iterParent)[i];
								dataCorrectRead.i2 = occur;
								idReadsCorrect.push_back(dataCorrectRead);
							}	
					}
					while(goRight(iterParent)){
						/*the length of current node*/
						int lenSiblings = lenParent+1;

						occur = countOccurrences(iterParent);
						if	(representative(iterParent) != nodeError &&
							!potentiallyErroneousNode(occur, expectedTheorical[lenSiblings], strictness, alg))
						{
							/*save the id and position(where the suffix begin in the reads) in the table of IDs correct*/
							/*also the number of occurrences*/
							for (unsigned i=0; i < length(getOccurrences(iterParent)); ++i)
							{
								Pair<Pair<unsigned,unsigned>,unsigned> dataCorrectRead;
								dataCorrectRead.i1 = getOccurrences(iterParent)[i];
								dataCorrectRead.i2 = occur;
								idReadsCorrect.push_back(dataCorrectRead);
							}
						}
					}
				}

				/*make the comarison between the substrings(the suffix after the position of error*/
				for (unsigned s=0; s<idReadsError.size();s++){
					/*the position where the error is*/
					unsigned positionError = idReadsError[s].i2+lenParent;
					/*is already identify as erroneus*/
					map<unsigned,Pair<unsigned,Pair<vector<unsigned>, Dna> > >::iterator mapID;
					mapID = allErrors.find(idReadsError[s].i1);
					
					Pair<unsigned,Pair<vector<unsigned>, Dna> >  correctAndData;
					/*Branch and Bound*/
					/*if already detect as erroneous(optimal)*/
					if(mapID!=allErrors.end()){
						/*the total length minus the position of error, and minus 1 mismatch*/
						int maxFitPossible = length(store.readSeqStore[idReadsError[s].i1])- positionError - 1;
						int optimalDistance = lenParent + maxFitPossible;
						correctAndData = (mapID->second);
					    
						/*
							Pair contains:
								vector with: position of errors, length parent and best fit;
								nucleotide which is correct
						*/
						Pair<vector<unsigned>, Dna>  data = correctAndData.i2;
						
						if(data.i1[1]+data.i1[2] > optimalDistance){
							continue;
						}
					}else{
						/*if this is the reverse complement*/
						map<unsigned,Pair<unsigned,Pair<vector<unsigned>, Dna> > >::iterator mapIDrevCompl;
					
						if (idReadsError[s].i1 >= readCount)
						{
							mapIDrevCompl = allErrors.find(idReadsError[s].i1-readCount);
							if (mapIDrevCompl!=allErrors.end())
							{
							/*the total length minus the position of error, and minus 1 mismatch*/
								int maxFitPossible = length(store.readSeqStore[idReadsError[s].i1])- positionError - 1;
								int optimalDistance = lenParent + maxFitPossible;
								correctAndData = (mapIDrevCompl->second);
					    
								/*
									Pair contains:
									vector with: position of errors, length parent and best fit;
									nucleotide which is correct
								*/
								Pair<vector<unsigned>, Dna>  data = correctAndData.i2;
						
								if(data.i1[1]+data.i1[2] > optimalDistance){
									continue;
								}
							}
						}
					}
				
					/*iterate over the sequence with error*/
					typedef typename Iterator<TRead, Standard>::Type TIterator;

					unsigned pos_correct = 0;
					unsigned best = 0;
					unsigned smallMismatch = 0;
					unsigned occurBest = 0;
					unsigned idCorr = 0;
					Dna ntCorr; 
					

					for (unsigned r=0; r<idReadsCorrect.size();r++){
						/*if not already correct by the same reads*/
						/*if(mapID!=allErrors.end()){
						
							if(idReadsCorrect[r].i1.i1==correctAndData.i1){	
								continue;
							}
						}*/

						/*the position in the read until which there is the same prefix*/
						unsigned positionCorrect = idReadsCorrect[r].i1.i2+lenParent;
											
						/*count to the position with error*/
						unsigned pos = 0;
						
						/*search the position detect as erroneous at the level of the read*/
						for (TIterator it = begin(store.readSeqStore[idReadsError[s].i1]); it != end(store.readSeqStore[idReadsError[s].i1]); ++it)
						{
							/* when we reache the position containig error*/	
							if (pos == positionError)
							{

								/*go to next nucleotide*/
								++it;
								/*compare the nucleotides after the error position*/
								unsigned nbChar = 0;
								unsigned mismatch = 0;
								/*compare until the end of the erroneous read or the corrected one*/
								while(it!=end(store.readSeqStore[idReadsError[s].i1]) && (positionCorrect+nbChar+1+mismatch)<(length(store.readSeqStore[idReadsCorrect[r].i1.i1])-1))
								{
									++pos;
									Dna currentNtErrRead  = store.readSeqStore[idReadsError[s].i1][pos];
									Dna currentNtCorrRead = store.readSeqStore[idReadsCorrect[r].i1.i1][positionCorrect+nbChar+1+mismatch];
									/*TODO: the user can choice the number of mismatch*/
									if (currentNtErrRead != currentNtCorrRead)
									{
										++mismatch;
										if (mismatch > acceptMismatch)
											break;
									} else
										++nbChar;
									++it;
								}
						

							    unsigned position = positionCorrect+nbChar+1;
								/*
									find the appopriate string if the end of the erroneous reads is reached  
									or at the end of the correct read
								*/
								/*I change just +mismatch*/
								if(it == end(store.readSeqStore[idReadsError[s].i1])|| position==(length(store.readSeqStore[idReadsCorrect[r].i1.i1]))){							
									/*if more matched characters than the other sibling*/
									if(nbChar+mismatch > best+smallMismatch){
										best = nbChar;
										smallMismatch = mismatch;
										idCorr = idReadsCorrect[r].i1.i1;
										ntCorr = store.readSeqStore[idReadsCorrect[r].i1.i1][positionCorrect];
										pos_correct = positionCorrect;
										occurBest = idReadsCorrect[r].i2;
									/*if the same number of characters, but the number of occurences is larger*/
									}else if (nbChar+mismatch == best+smallMismatch){// && occurBest < idReadsCorrect[r].i2){
										if(nbChar > best){
											best = nbChar;
											smallMismatch = mismatch;
											idCorr = idReadsCorrect[r].i1.i1;
											ntCorr = store.readSeqStore[idReadsCorrect[r].i1.i1][positionCorrect];
											pos_correct = positionCorrect;
											occurBest = idReadsCorrect[r].i2;
										/*when they are equal - test the number of ocurrences*/		
										}else if(nbChar == best && occurBest < idReadsCorrect[r].i2){
											best = nbChar;
											smallMismatch = mismatch;
											idCorr = idReadsCorrect[r].i1.i1;
											ntCorr = store.readSeqStore[idReadsCorrect[r].i1.i1][positionCorrect];
											pos_correct = positionCorrect;
											occurBest = idReadsCorrect[r].i2;
										}
									
									}
								}
								break;
							}
							++pos;
						}
					}
			
			
					if(idCorr !=0){
					
						/*the correct nucleotide*/
						Dna nucleotide =  ntCorr;
						
						/*copie the ID*/
						unsigned idReadErrorPos = idReadsError[s].i1;
						
						/*if so: this is the reverse complement*/
						if ((idReadsError[s].i1) >= readCount)
						{
							/*	ID of read give by the user, 
								which is the reverse complement of the identified 
								as erroneous current read
							*/
							idReadErrorPos = idReadsError[s].i1 - readCount;

							/*position for correction in the complement - read user*/
							positionError = length(store.readSeqStore[idReadsError[s].i1]) - positionError - 1;
							FunctorComplement<Dna5> comp;
							nucleotide = comp(nucleotide);
						}
						
						/*search if we have identify the same read already as erroneous*/
						map<unsigned,Pair<unsigned,Pair<vector<unsigned>, Dna> > >::iterator mapID;
						mapID = allErrors.find(idReadErrorPos);
	

						if(mapID != allErrors.end()){
			
							Pair<unsigned,Pair<vector<unsigned>, Dna> > correctID = (mapID->second);
							
							Pair<vector<unsigned>,Dna> number = correctID.i2;
							/*This is maybe not necessary, because we search to maximisate*/
							/*when the position or the nucleotide differ*/
							//if(number.i1[0]!=positionError || number.i2 != nucleotide){
							
							/*find longer commun part with certain nb accepted mismatch*/
							if(number.i1[1]+number.i1[2] + number.i1[3] < lenParent + best + smallMismatch)
							{
								
								allErrors[idReadErrorPos] = dataErroneousNodes(idCorr, positionError, lenParent, best, smallMismatch, occurBest, nucleotide);
							/*when they are equal*/
							} else if(number.i1[1]+number.i1[2] + number.i1[3]==lenParent + best + smallMismatch){
								/*looking for longer matching part without mismatch*/
								if(number.i1[1]+number.i1[2]<lenParent + best){
					
									allErrors[idReadErrorPos] = dataErroneousNodes(idCorr, positionError, lenParent, best, smallMismatch, occurBest, nucleotide);
								/*eqaul number of match and mismatch, but greater nb occurrences*/
								}else if(number.i1[1]+number.i1[2]==lenParent + best &&number.i1[4]<occurBest){

										allErrors[idReadErrorPos] = dataErroneousNodes(idCorr, positionError, lenParent, best, smallMismatch, occurBest, nucleotide);
										
									}
								}
								break;
							//}
						}else{
								allErrors[idReadErrorPos] = dataErroneousNodes(idCorr, positionError, lenParent, best, smallMismatch, occurBest, nucleotide);
						}
					}
				}
		
			}		
		
		}		
		++iter;
	}
	end1 = clock();
	cout << "Time for searching between given levels: "<< (double)(end1-start1)/CLOCKS_PER_SEC << " seconds." << endl;
	ofstream out("id");
	map<unsigned,Pair<unsigned,Pair<vector<unsigned>, Dna> > >::iterator it;

	for (it = allErrors.begin (); it != allErrors.end (); ++it)
	{
		unsigned readId = it->first;
		if (readId > readCount)
			readId -= readCount;
		out << readId << " 0" << endl; 
	}
	
	correctErroneousReads(allErrors, store);
	cout << "Total corrected reads number is "<< allErrors.size() << endl;
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
	double genomeLength,
	bool poisson,
	double strictness,
	int fromLevel,
	int toLevel,
	int acceptMismatch)
{		   
	clock_t start, end1;
	start = clock(); 
	
	cout << "Construct suffix array"<<endl;

	// append their reverse complements
	unsigned readCount = length(store.readSeqStore);
	Dna5String tmp;
	if (genomeLength != 1.0)
		for (unsigned i = 0; i < readCount; ++i)
		{
			tmp = store.readSeqStore[i];
			reverseComplementInPlace(tmp);
			appendValue(store.readSeqStore, tmp);
		}

	// construct suffix array of the set of reads
	TFionaIndex myIndex(store.readSeqStore); 

	/*iterator with restrictions*/
	typedef Iterator<TFionaIndex, TopDown<ParentLinks<Postorder> > >::Type TConstrainedIterator; 
	TConstrainedIterator myConstrainedIterator(myIndex); 

	/*calculate the frequency for each nucleotide, didn't use for the moment*/
	/*'A' = 0, 'C' = 1, 'G' = 2, 'T' = 3*/
	map<unsigned,double> frequency;
	frequency = determineFrequency(myConstrainedIterator);

	/*restrictions just for estimate the genome length if there is no data*/

#ifdef MEDIAN
	unsigned level = fromLevel;
	ofstream out("medianLevels.txt"); 
	ofstream median("medianForEachLevel.txt");
#endif

	if (genomeLength == 1.0)
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
				
				for (int k = 0; k < length(freq); ++k)
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
			
	if (genomeLength == 0.0)
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
			double nbOccurence = 0.0; 
			
			median << level << " " << medianLevel(myConstrainedIterator) << endl;
			goBegin(myConstrainedIterator);
			while (!atEnd(myConstrainedIterator))
			{
				if (parentRepLength(myConstrainedIterator) > level)
				{
					nbOccurence = countOccurrences(myConstrainedIterator);
					out << level << " " << nbOccurence << endl; 
				}
				++myConstrainedIterator;
			}
		}
#else
		//int logRation = log10(readCount) / log10(4);
		//int l = logRation + 1;
		//cout << l << endl;
		cargo(myIndex).replen_min = fromLevel;
		cargo(myIndex).replen_max = fromLevel+2; 	
		cargo(myIndex).frequency = frequency;

		double expectedValueGivenLevel = medianLevel(myConstrainedIterator);

		/*TODO the length for the reads is always the same, but for real cases it is better to change*/
		/*estimate the genome length by the use of expected value*/
		unsigned readLength = length(store.readSeqStore[0]);

		/* a = readLength - path_label + 1 */
		/*here plus 1 also because the level is between fromLevel and toLevel*/
		double a = readLength - fromLevel + 2;
		double mult = readCount*a;
		genomeLength = mult/expectedValueGivenLevel;
		cout << "The estimated genome length is " << genomeLength << endl;
#endif
	}
	end1 = clock();
	cout << "Time required for suffix array construction : "<< (double)(end1-start)/CLOCKS_PER_SEC << " seconds." << "\n\n";
	
	vector<unsigned> idCorrected;
	
	/*restrictions for the searching levels*/
	cargo(myIndex).replen_min = fromLevel;
	cargo(myIndex).replen_max = toLevel; 	
	cargo(myIndex).frequency = frequency;

	/*the core of the correction method*/
	if (poisson)
		/*use of p-value like a limit*/
		searchNode(myConstrainedIterator, store, genomeLength, strictness, acceptMismatch, FionaPoisson());
	else
		/*use an expected value for a certain level*/
		searchNode(myConstrainedIterator, store, genomeLength, strictness, acceptMismatch, FionaExpected());
	
	// remove reverse complements
	resize(store.readSeqStore, readCount);
}

int main(int argc, char* argv[]) 
{ 
	clock_t start, end;
	
	start = clock();
	
	if(argc > 14){
		cout << "Usage: ./Fiona [options: -f $1 -p $1 -i $1 -l $1 $2 -g $1 -m $1] [input_reads_file] [corrected_reads_file]" << endl;
		cout << "See README.txt for details" << endl;
		exit(1);
	}

	bool programmeOprions = true;

	/*the levels between the search will take place*/
	unsigned fromLevel = 0;
	unsigned toLevel = 0;

	/*the methode that will be use*/
	bool poisson=false;

	double genomeLength = 0.0;

	/*p-value*/
	double strictness = 0.0001;
	
	/**number accepted mismatch - default value 1*/
	unsigned acceptMismatch = 1;

	String<char> inFileName = "";
	String<char> outFileName = "";

	/*default value*/
	unsigned numCycles = 3;
	bool expected = false;
	for (unsigned i = 1; i < length(argv); ++i){
		String<char> option = argv[i];
		
		if(strcmp(toCString(option[0]), "-")==0){
			while(programmeOprions){
				option = argv[i];
				/*Method with strictness value will be use - need a strictness to be done*/
				if(strcmp(toCString(option),"-f")==0){
					if(poisson){
						cout << " "<< endl;
						cout << "You have already choiced the method with p-value "<<endl;
						cout << "Use only one of the options -f or -p, for more detail see README.txt "<<endl; 
						cout << "Exit" <<endl;
						exit(2);
					}
					expected = true;
					strictness = atof(argv[i+1]);
					cout << "The correction method will use the strictness value " << strictness <<endl;
					i+=2;
				}else if(strcmp(toCString(option),"-p")==0){
					if(expected){
						cout << " "<< endl;
						cout << "You have already choiced the method with confidential interval "<<endl; 							cout << "Use only one of the options -f or -p, for more detail see README.txt "<<endl;
						cout << "Exit" <<endl;
						exit(3);
					}
					poisson = true;
					strictness = atof(argv[i+1]);
					cout << "The p-value is set to " << strictness << endl;
					i+=2;
				}else if(strcmp(toCString(option),"-i")==0){
					numCycles = atoi(argv[i+1]);
					cout << "The number of cycle is set to " << numCycles << endl;
					i+=2;
				}else if(strcmp(toCString(option),"-l")==0){
					fromLevel = atoi(argv[i+1]);
					toLevel = atoi(argv[i+2]);
					cout << "The searching levels are between " << fromLevel << " and " << toLevel << endl;
					i+=3;
				}else if(strcmp(toCString(option),"-g")==0){
					genomeLength = atof(argv[i+1]);
					cout << "The genome length is set to " << genomeLength << endl;
					i+=2;
				}else if(strcmp(toCString(option),"-m")==0){
					acceptMismatch = atoi(argv[i+1]);
					cout << "The number of accepted mismatches is " << acceptMismatch << endl;
					i+=2;
				}else{
					option = argv[i];
					if(strcmp(toCString(option[0]), "-")!=0){
						programmeOprions = false;
					}else{
						cout << "There is not option " << option << endl;
						exit(4);
					}
				}
			}
					
		}

		/*defalt method use p-value*/
		if(!poisson && !expected){
			poisson = true;
		}
		inFileName = argv[i];
		outFileName = argv[i+1];
		break;		
	}
	
	
	cout << endl;
	
	// load original set of reads
	TFionaFragStore store;
	loadReads(store, inFileName);
	
	// initialise the top and down level by using the log4 from the total number of reads
	if (fromLevel==0)
	{
		int logRation = log10(length(store.readSeqStore)) / log10(4.0);
		fromLevel = logRation + 2;
		toLevel   = fromLevel + 10;
		cout << "The top level is " << fromLevel << " and the down level is " << toLevel <<endl;
	}
	cout << endl;

	for (unsigned cycle = 1; cycle <= numCycles; ++cycle)
	{
		cout << "Cycle "<< cycle <<" from "<< numCycles << endl;
		correctReads(store, genomeLength, poisson, strictness, fromLevel, toLevel, acceptMismatch);
		
		if (acceptMismatch > 0) --acceptMismatch;			

		/*TODO maybe to stop if there is not reads corrected in the cycle before*/
		/*if so after each iteration must save the ID for the reads which are corrected*/
		/*thus we can also show the total number of reads that are corrected at the final stage*/
	}
	
	/*write in file all input reads with the corrected one*/
	ofstream out(toCString(outFileName));
	int numCorrected = 0;
	for (unsigned i = 0; i < length(store.readNameStore); ++i)
	{
		/*to give the number of reads corrected for several iteration*/ 
		if (strContains(toCString(store.readNameStore[i]), "corrected"))
			++numCorrected;

		out << ">" << store.readNameStore[i]<<endl;
		out << store.readSeqStore[i] << endl;
	}

	if (numCycles > 1)
		cout << "Total number reads corrected for " << numCycles << " cycles is " << numCorrected <<endl; 

	end = clock();
	cout << endl << "Time required for execution: "<< (double)(end-start)/CLOCKS_PER_SEC << " seconds." << "\n\n";

	return 0;
}

	
