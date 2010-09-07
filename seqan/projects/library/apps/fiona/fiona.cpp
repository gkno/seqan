#include <iostream>
#include <fstream>
#include <seqan/index.h>
#include <seqan/basic.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sstream>

//#define MEDIAN

using namespace std;
using namespace seqan;

/*TODO cumulative poisson distribution determinated just once for a level, not at each node like now*/

/*matching string*/
bool strContains(const string inputStr, const string searchStr)
{
        size_t contains;
 
        contains = inputStr.find(searchStr);
 
        if(contains != string::npos)
                return true;
        else
                return false;
}

/*save in table, information about the erroneous reads and the respective correct one */
Pair<int,Pair<vector<int>, Dna> > dataErroneousNodes(int idCorr, int positionError, int lenParent, int best, 
						    int smallMismatch, int occurBest, Dna nucleotide){
													 
	Pair<int,Pair<vector<int>, Dna> > pairCorrect;	

	/*data for the erroneous node and correct one*/
	/*position error, length parent(path-label), nb match after error position, nb mismatch, nb of occurrences */ 	
	vector<int> data;	
									
	data.push_back(positionError);
	data.push_back(lenParent);
	data.push_back(best);
	data.push_back(smallMismatch);
	data.push_back(occurBest);
 
	Pair<vector<int>, Dna> dataCorrection;
	dataCorrection.i1 = data;
	dataCorrection.i2 = nucleotide; 
	pairCorrect.i1 =  idCorr;
	pairCorrect.i2 = dataCorrection;
	return pairCorrect;
}


/*Expected value - general, use if all reads have the same length*/
vector <double> expectedValueEqualReadsLength(StringSet<String<Dna>,Owner < > > readsTable, double genomeLength)
{
		double read_length = length(readsTable[0]);

		/*without reverse complement*/
		double numberReads = length(readsTable)/2;
		/* 
			E(m)= (read_length - suffix_length + 1)*numberReads/genomeLength 
		*/
		vector <double> expected;

		for (double suffix_length=0; suffix_length < read_length+1; suffix_length++){
			expected.push_back(((read_length - suffix_length + 1)*numberReads)/genomeLength );		
		}
		return expected;
}

/*Expected value for set of reads with different length*/
vector <double> expectedValueTheorical(StringSet<String<Dna>,Owner < > > readsTable, double genomeLength){

		/* 
			E(m)= (read_length - suffix_length + 1)*numberReads/genomeLength 
		*/
		map <int,double> expected;
		vector<double> d;
		map<int,int> vectorLength;

		/*without reverse complement*/
		for(int i=0;i < length(readsTable);i=i+2){
			vectorLength[length(readsTable[i])]+=1;
		}
	
		/* a = read_length - suffix_length + 1 */
		double a = 0.0;
		map<int,int>::iterator iterMap;
		for (iterMap = vectorLength.begin (); iterMap != vectorLength.end (); ++iterMap){
			for(int suffix_length=0; suffix_length < iterMap->first; suffix_length++){
				a = iterMap->first - suffix_length +1;
				double ration = (a*iterMap->second)/genomeLength;
				expected[suffix_length]+=ration;		
			}	
		}
		map<int,double>::iterator iMap;
		for (iMap = expected.begin (); iMap != expected.end (); ++iMap){
			d.push_back(iMap->second);
		}	
		return d;
}

/* Standard Deviation */
vector <double> standardDeviation(StringSet<String<Dna>,Owner < > > readsTable, double genomeLength)
	{
		double read_length = length(readsTable[0]);

		/*without reverse complement*/
		double numberReads = length(readsTable)/2;

		/* 
			SD(m)= numberReads*((read_length - suffix_length + 1)/genomeLength 
					- (read_length - suffix_length + 1)^2/genomeLength ^2)
		*/

		vector <double> standardDev;
		double valueFirst;
		double valueSecond;
		for (double suffix_length=0; suffix_length < read_length+1; suffix_length++){
			valueFirst  = (read_length - suffix_length + 1);
			valueFirst  = valueFirst/genomeLength;
			valueSecond = pow(valueFirst,2);
			standardDev.push_back(sqrt(numberReads*(valueFirst - valueSecond)));		
		}
		return standardDev;
}

/*factoriel*/
double factoriel(double n)
{
	double factoriel, i;
	factoriel=1;
	i=1;
	while (i<=n)
	{
		factoriel=factoriel*i;
		i=i+1;
	}
	return factoriel;
}

/*cumulative poisson distribution*/
double pValue(const unsigned int k, const double mean){
	//return gsl_cdf_poisson_P(k,mean);
	double pValue = 0.0;
	for(double i=0.0; i <=k; i++){
		double a = exp(-mean);
		double b = pow(mean,i);
		double c = a*b;
		double f = factoriel(i); 
		pValue+= c/f;	
	}
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

	map<int,int> vectorOccurences;

	goBegin(iter);
	while (!atEnd(iter)){
		nbOccurence = countOccurrences(iter);
		vectorOccurences[nbOccurence]+=1.0;
		totalOccurence+=nbOccurence;
		iter++; 
	}
	
	mediumTotalOccurence = totalOccurence/2.0;
	
	map<int,int>::iterator iterMap;
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

double probabilityOneError(double percentageErr, int lenPath){
	double percent = 1 - percentageErr;
	return (1-pow(percent,lenPath));

}

/*change the erroneous nucleotide in all reads identify with errors*/
Pair<StringSet<String<Dna>,Owner < > >, StringSet<String<char> > > correctErroneousRead(map<int,Pair<int,Pair<vector<int>, Dna> > > allErrors, 
	StringSet<String<char> > fastaComment, StringSet<String<Dna>,Owner < > > setReads){
	
	map<int,Pair<int,Pair<vector<int>, Dna> > >::iterator iR;
	int comment = 0;
	for (iR = allErrors.begin (); iR != allErrors.end (); ++iR)
	{
		Pair<int,Pair<vector<int>, Dna> > setR = (iR->second);
		int posErr = (setR.i2).i1[0];
		
		/*information about the read*/
		ostringstream oss;
		ostringstream m;
		if (strContains(toCString(fastaComment[(iR->first)/2]), "corrected")){
			m << "," << posErr;
			oss << toCString(fastaComment[iR->first/2]) << m.str();
			String<char> d = oss.str();
			fastaComment[(iR->first)/2] = d;
		}else{
			m << " corrected: " << posErr;
			oss << toCString(fastaComment[iR->first/2])<< m.str();
			String<char> d = oss.str();
			fastaComment[(iR->first)/2] = d;
		}

		int pos=0;
		typedef Iterator<String<Dna> >::Type TIterator;
		for (TIterator it = begin(setReads[(iR->first)]); it != end(setReads[(iR->first)]); ++it){
			if(pos==(setR.i2).i1[0]){
				assignValue(it,(setR.i2).i2);
				break;
			}
			pos++;
		}

	}
	Pair<StringSet<String<Dna>,Owner < > >, StringSet<String<char> > > fileData;
	fileData.i1 = setReads;
	fileData.i2 = fastaComment;
	return fileData;
}


/*restriction for certain levels - between max and min, also table with frequency may be to use eventually TODO*/
struct TMyConstraints { 
	 unsigned int replen_max; 
	 unsigned int replen_min;
	 map<int,double> frequency;
	 bool _cachedPred; 
}; 
 
namespace seqan  {

	struct TMyIndex; 
	template <typename TText> 
	 
	struct Cargo<Index<TText> > { 
		typedef TMyConstraints Type; 
	}; 

	/*TODO THIS FONCTION CAN BE CHANGED FOR THE FREQUENCY - here just one experience*/
	/*by the use also the frequency for A,T,G,C*/
	/*higher frequency - high level as min in which we will begin the searching*/
	
	/*hide the node between certain level*/
	template <typename TText, typename TSpec>
	bool nodePredicate(Iter<Index<TText >, TSpec> &it) 
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

/*complement of one nucleotide: A==T, G==C*/
Dna complementaire(Dna nucleotide){

	/*'A' = 0, 'C' = 1, 'G' = 2, 'T' = 3*/
	if(nucleotide == 0) return 3;
	if(nucleotide == 1) return 2;
	if(nucleotide == 2) return 1;
	if(nucleotide == 3) return 0;
	else{
		cout << "Accepted Character :A,T,G and C";
		exit(1);
	}
}


/*TODO make just one fonction for p-value and expected by just changing the tests if poisson or the other*/
/*detect and repare the reads with errors*/
template < typename TMyIndex, class TSpec >
Pair<StringSet<String<Dna>,Owner < > >, StringSet<String<char> > >
searchNodePoisson(Iter< TMyIndex, VSTree<TSpec> > iter,
               	   StringSet<String<Dna>,Owner< > > setReads, StringSet<String<char> > fastaComment,
		   double genomeLength, double strictness, int acceptMismatch){
	
	cout << "Method with p-value and Poisson distribution" << endl;
	cout << "Searching... " << endl;
	
	/*table with the theoricals values*/
	vector <double> expectedTheorical = expectedValueTheorical(setReads, genomeLength);
	
	/*save reads with errors*/
	vector<String<Dna> > readsErrors;

	/*
		map contains : 
			key - ID of read with errors;
			value - pair of position in read and nucleotide for make change with
	*/
	map<int,Pair<int,Dna> > readsErrorsID;

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
	map<int,Pair<int,Pair<vector<int>, Dna> > > allErrors;

	/* length of parent label*/
	int lenParent = 0;

	/* length of label to root until courant node*/
	int lenPath = 0;
	
	clock_t start1, end1;

	start1 = clock();

	/*the bigining of the tree*/
	goBegin(iter);

	while (!atEnd(iter)){

		/*affectation*/
		lenParent = parentRepLength(iter);
		lenPath = repLength(iter);
	
	  SEQAN_ASSERT_LT(lenParent + 1, expectedTheorical.size());
		double pValue_Poisson =  pValue(countOccurrences(iter), expectedTheorical[lenParent+1]);
		
		/*compare the cumulative poisson distribution with the p-value(strictness)*/
		if(pValue_Poisson <= strictness){
			nodeError = representative(iter);
			/*
				save the id and position(where the suffix begin in the reads) for the suspected reads belonging to this level
			*/
			vector<Pair<int,int> > idReadsError;

			/*take all reads containing errors at that position in the tree*/
			for(int i=0; i < length(getOccurrences(iter));i++){
				idReadsError.push_back(getOccurrences(iter)[i]);
			}
			/*if there is a nodes detect as erroneus*/
			if(idReadsError.size()>0){
				/*copy the iterator for iterate over the siblings*/
				typename Iterator< TMyIndex, TopDown< ParentLinks<Postorder> > >::Type iterParent(iter);
			 
				/*go at the level of the parent node*/
				goUp(iterParent);
					
				/*
					potential reads for make the correction,
					because at the same level, with the same prefix
				*/
				vector<Pair<Pair<int,int>,int> > idReadsCorrect;
		
				if(goDown(iterParent)){
					/*the length of current node*/
					int lenSiblings = lenParent+1;
					/*nb occurrence for this node*/
					int occur = countOccurrences(iterParent);
					double pValue_Poisson_Sibling =  pValue(occur, expectedTheorical[lenParent+1]);
					/*when different from the erroneous node & the probability to have an error is greater or equal from p-value*/
					if(representative(iterParent) != nodeError && strictness <= pValue_Poisson_Sibling){
							/*save the id and position(where the suffix begin in the reads) in the table of IDs correct*/
							/*also the number of occurrences*/
							for(int i=0; i < length(getOccurrences(iterParent));i++){
								Pair<Pair<int,int>,int> dataCorrectRead;
								dataCorrectRead.i1 = getOccurrences(iterParent)[i];
								dataCorrectRead.i2 = occur;
								idReadsCorrect.push_back(dataCorrectRead);
							}	
					}
					while(goRight(iterParent)){
						/*the length of current node*/
						int lenSiblings = lenParent+1;

						occur = countOccurrences(iterParent);
						pValue_Poisson_Sibling =  pValue(occur, expectedTheorical[lenParent+1]);
						/*when different from the erroneous node & the probability to have an error is greater or equal from p-value*/
						if(representative(iterParent) != nodeError && strictness <= pValue_Poisson_Sibling ){
							/*save the id and position(where the suffix begin in the reads) in the table of IDs correct*/
							/*also the number of occurrences*/
							for(int i=0; i < length(getOccurrences(iterParent));i++){
								Pair<Pair<int,int>,int> dataCorrectRead;
								dataCorrectRead.i1 = getOccurrences(iterParent)[i];
								dataCorrectRead.i2 = occur;
								idReadsCorrect.push_back(dataCorrectRead);
							}
						}
					}
				}

				/*make the comarison between the substrings(the suffix after the position of error*/
				for(int s=0; s<idReadsError.size();s++){
					/*the position where the error is*/
					int positionError = idReadsError[s].i2+lenParent;
					/*is already identify as erroneus*/
					map<int,Pair<int,Pair<vector<int>, Dna> > >::iterator mapID;
					mapID = allErrors.find(idReadsError[s].i1);
					
					Pair<int,Pair<vector<int>, Dna> >  correctAndData;
					/*Branch and Bound*/
					/*if already detect as erroneous(optimal)*/
					if(mapID!=allErrors.end()){
						/*the total length minus the position of error, and minus 1 mismatch*/
						int maxFitPossible = length(setReads[idReadsError[s].i1])- positionError - 1;
						int optimalDistance = lenParent + maxFitPossible;
						correctAndData = (mapID->second);
					    
						/*
							Pair contains:
								vector with: position of errors, length parent and best fit;
								nucleotide which is correct
						*/
						Pair<vector<int>, Dna>  data = correctAndData.i2;
						
						if(data.i1[1]+data.i1[2] > optimalDistance){
							continue;
						}
					}else{
						/*if this is the reverse complement*/
						map<int,Pair<int,Pair<vector<int>, Dna> > >::iterator mapIDrevCompl;
					
						if(idReadsError[s].i1%2!=0){
							mapIDrevCompl = allErrors.find(idReadsError[s].i1-1);
							if(mapIDrevCompl!=allErrors.end()){
							/*the total length minus the position of error, and minus 1 mismatch*/
								int maxFitPossible = length(setReads[idReadsError[s].i1])- positionError - 1;
								int optimalDistance = lenParent + maxFitPossible;
								correctAndData = (mapIDrevCompl->second);
					    
								/*
									Pair contains:
									vector with: position of errors, length parent and best fit;
									nucleotide which is correct
								*/
								Pair<vector<int>, Dna>  data = correctAndData.i2;
						
								if(data.i1[1]+data.i1[2] > optimalDistance){
									continue;
								}
							}
						}
					}
				
					/*iterate over the sequence with error*/
					typedef Iterator<String<Dna> >::Type TIterator;

					int pos_correct = 0;
					int best = 0;
					int smallMismatch = 0;
					int occurBest = 0;
					int idCorr = 0;
					Dna ntCorr; 
					

					for(int r=0; r<idReadsCorrect.size();r++){
						/*if not already correct by the same reads*/
						/*if(mapID!=allErrors.end()){
						
							if(idReadsCorrect[r].i1.i1==correctAndData.i1){	
								continue;
							}
						}*/

						/*the position in the read until which there is the same prefix*/
						int positionCorrect = idReadsCorrect[r].i1.i2+lenParent;
											
						/*count to the position with error*/
						int pos = 0;
						
						/*search the position detect as erroneous at the level of the read*/
						for (TIterator it = begin(setReads[idReadsError[s].i1]); it != end(setReads[idReadsError[s].i1]); ++it){
							/* when we reache the position containig error*/	
							if(pos == positionError){

								/*copy iterator at the position where the error is*/
								TIterator errITer = it;
								/*go to next nucleotide*/
								it++;
								/*compare the nucleotides after the error position*/
								int nbChar = 0;
								int mismatch = 0;
								/*compare until the end of the erroneous read or the corrected one*/
								while(it!=end(setReads[idReadsError[s].i1]) && (positionCorrect+nbChar+1+mismatch)<(length(setReads[idReadsCorrect[r].i1.i1])-1)){
									pos++;
									Dna currentNtErrRead  = setReads[idReadsError[s].i1][pos];
									Dna currentNtCorrRead = setReads[idReadsCorrect[r].i1.i1][positionCorrect+nbChar+1+mismatch];
									/*the user can choice the number of mismatch*/
									if(currentNtErrRead != currentNtCorrRead){
										mismatch++;
										if(mismatch>acceptMismatch){
											break;
										}
									}else{
										nbChar++;
									}
									it++;
								}
						

							    int position = positionCorrect+nbChar+1;
								/*
									find the appopriate string if the end of the erroneous reads is reached  
									or at the end of the correct read
								*/
								/*I change just +mismatch*/
								if(it == end(setReads[idReadsError[s].i1])|| position==(length(setReads[idReadsCorrect[r].i1.i1]))){							
									/*if more matched characters than the other sibling*/
									if(nbChar+mismatch > best+smallMismatch){
										best = nbChar;
										smallMismatch = mismatch;
										idCorr = idReadsCorrect[r].i1.i1;
										ntCorr = setReads[idReadsCorrect[r].i1.i1][positionCorrect];
										pos_correct = positionCorrect;
										occurBest = idReadsCorrect[r].i2;
									/*if the same number of characters*/
									}else if (nbChar+mismatch == best+smallMismatch){
										/*case when the number of mismatch is greater*/			
										if(nbChar > best){
											best = nbChar;
											smallMismatch = mismatch;
											idCorr = idReadsCorrect[r].i1.i1;
											ntCorr = setReads[idReadsCorrect[r].i1.i1][positionCorrect];
											pos_correct = positionCorrect;
											occurBest = idReadsCorrect[r].i2;
										/*when they are equal - test the number of ocurrences*/		
										}else if(nbChar == best && occurBest < idReadsCorrect[r].i2){
											best = nbChar;
											smallMismatch = mismatch;
											idCorr = idReadsCorrect[r].i1.i1;
											ntCorr = setReads[idReadsCorrect[r].i1.i1][positionCorrect];
											pos_correct = positionCorrect;
											occurBest = idReadsCorrect[r].i2;
										}
									
									}
								}
								break;
							}
							pos++;
						}
					}
			
			
					if(idCorr !=0){
					
						/*the correct nucleotide*/
						Dna nucleotide =  ntCorr;
						
						/*copy the ID*/
						int idReadErrorPos = idReadsError[s].i1;
						
						/*if so: this is the reverse complement*/
						if((idReadsError[s].i1)%2!=0){
							/*	ID of read give by the user, 
								which is the reverse complement of the identified 
								as erroneous current read
							*/
							idReadErrorPos = idReadsError[s].i1 - 1;

							/*position for correction in the complement - read user*/
							positionError = length(setReads[idReadsError[s].i1]) - positionError - 1; 
							nucleotide = complementaire(nucleotide);
						}

// debug					if(idReadsError[s].i1==176 || idReadsError[s].i1 ==177){cout << representative(iter);}

						/*search if we have identify the same read already as erroneous*/
						map<int,Pair<int,Pair<vector<int>, Dna> > >::iterator mapID;
						mapID = allErrors.find(idReadErrorPos);

						if(mapID != allErrors.end()){
			
							Pair<int,Pair<vector<int>, Dna> > correctID = (mapID->second);
							
							Pair<vector<int>,Dna> number = correctID.i2;
							/*This is maybe not necessary, because we search to maximisate*/
							/*when the position or the nucleotide differ*/
							//if(number.i1[0]!=positionError || number.i2 != nucleotide){
							
							/*find longer commun part with certain nb accepted mismatch*/
							if(number.i1[1]+number.i1[2] + number.i1[3] < lenParent + best + smallMismatch){
								
								allErrors[idReadErrorPos] = dataErroneousNodes(idCorr, positionError, lenParent, best, smallMismatch, occurBest, nucleotide);
							/*when they are equal*/
							}else if(number.i1[1]+number.i1[2] + number.i1[3]==lenParent + best + smallMismatch){
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
		iter++;
	}

	end1 = clock();
	cout << "Time for searching between given levels: "<< (double)(end1-start1)/CLOCKS_PER_SEC << " seconds." << endl;
	ofstream out("id");
	map<int,Pair<int,Pair<vector<int>, Dna> > >::iterator it;

	for (it = allErrors.begin (); it != allErrors.end (); ++it)
    	{
		out << ((it->first)/2)<< " 0"<<endl; 
	}

	Pair<StringSet<String<Dna>,Owner < > >, StringSet<String<char> > > fileData = correctErroneousRead(allErrors,fastaComment, setReads);

	cout << "Total corrected reads number is "<< allErrors.size()<<"\n\n";
	
	return fileData;
}

/*detect and repare the reads with errors*/
template < typename TMyIndex, class TSpec >
Pair<StringSet<String<Dna>,Owner < > >, StringSet<String<char> > >
searchNodeExpected(Iter< TMyIndex, VSTree<TSpec> > iter,
		   StringSet<String<Dna>,Owner< > > setReads, StringSet<String<char> > fastaComment,
		   double genomeLength, double strictness, int acceptMismatch){

	cout << "Method with expected value for each level" << endl;
	cout << "Searching... " << endl;
	
	/*table with the theoricals values*/
	vector <double> expectedTheorical = expectedValueTheorical(setReads, genomeLength);
	vector <double> sd = standardDeviation(setReads, genomeLength);

	/*The strictness value allows to estimate the confidential intervall*/
	for (int i=0; i < expectedTheorical.size();i++){
		double expectedTemporary = expectedTheorical[i] - strictness*sd[i];
		
		/*If the connfidential intervall take value less than 1 ??? not sure for that*/
		/*if(expectedTemporary < 1){
			expectedTheorical[i] = 1.1;
		}else{*/
			expectedTheorical[i] = expectedTemporary;
		//}
		//if fixed
		//expectedTheorical[i]=5;
	}
	
	/*save reads with errors*/
	vector<String<Dna> > readsErrors;

	/*
		map contains : 
			key - ID of read with errors;
			value - pair of position in read and nucleotide for make change with
	*/
	map<int,Pair<int,Dna> > readsErrorsID;

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
	map<int,Pair<int,Pair<vector<int>, Dna> > > allErrors;

	/* length of parent label*/
	int lenParent = 0;

	/* length of label to root until courant node*/
	int lenPath = 0;

	clock_t start1, end1;

	start1 = clock();
	
	/*the bigining of the tree*/
	goBegin(iter);
	while (!atEnd(iter)){

		/*affectation*/
		lenParent = parentRepLength(iter);
		lenPath = repLength(iter);
	
		/*compare the weight for a node with its expected value(dependent on the length)*/
		if(countOccurrences(iter) < expectedTheorical[lenParent+1]){
			nodeError = representative(iter);
	
			/*
				save the id and position(where the suffix begin in the reads) for suspected nodes
				for which can be find more optimal correction
			*/
			vector<Pair<int,int> > idReadsError;

			/*take all reads containing errors at that position in the tree*/
			for(int i=0; i < length(getOccurrences(iter));i++){
				idReadsError.push_back(getOccurrences(iter)[i]);
			}
			/*if there is a nodes detect as erroneus*/
			if(idReadsError.size()>0){
				
				typename Iterator< TMyIndex, TopDown< ParentLinks<Postorder> > >::Type iterParent(iter);
			 
				/*go at the level of the parent node*/
				goUp(iterParent);
					
				/*
					potential reads for make the correction,
					because at the same level, with the same prefix
				*/
				vector<Pair<Pair<int,int>,int> > idReadsCorrect;
		
				if(goDown(iterParent)){
					/*the length of current node*/
					int lenSiblings = lenParent+1;
					int occur = countOccurrences(iterParent);
					if(representative(iterParent) != nodeError && occur >= expectedTheorical[lenSiblings]){
							/*save the id and position(where the suffix begin in the reads) in the table of IDs correct*/
							/*also the number of occurrences*/
							for(int i=0; i < length(getOccurrences(iterParent));i++){
								Pair<Pair<int,int>,int> dataCorrectRead;
								dataCorrectRead.i1 = getOccurrences(iterParent)[i];
								dataCorrectRead.i2 = occur;
								idReadsCorrect.push_back(dataCorrectRead);
							}	
					}
					while(goRight(iterParent)){
						/*the length of current node*/
						int lenSiblings = lenParent+1;

						occur = countOccurrences(iterParent);
						if(representative(iterParent) != nodeError && occur >= expectedTheorical[lenSiblings]){
							/*save the id and position(where the suffix begin in the reads) in the table of IDs correct*/
							/*also the number of occurrences*/
							for(int i=0; i < length(getOccurrences(iterParent));i++){
								Pair<Pair<int,int>,int> dataCorrectRead;
								dataCorrectRead.i1 = getOccurrences(iterParent)[i];
								dataCorrectRead.i2 = occur;
								idReadsCorrect.push_back(dataCorrectRead);
							}
						}
					}
				}

				/*make the comarison between the substrings(the suffix after the position of error*/
				for(int s=0; s<idReadsError.size();s++){
					/*the position where the error is*/
					int positionError = idReadsError[s].i2+lenParent;
					/*is already identify as erroneus*/
					map<int,Pair<int,Pair<vector<int>, Dna> > >::iterator mapID;
					mapID = allErrors.find(idReadsError[s].i1);
					
					Pair<int,Pair<vector<int>, Dna> >  correctAndData;
					/*Branch and Bound*/
					/*if already detect as erroneous(optimal)*/
					if(mapID!=allErrors.end()){
						/*the total length minus the position of error, and minus 1 mismatch*/
						int maxFitPossible = length(setReads[idReadsError[s].i1])- positionError - 1;
						int optimalDistance = lenParent + maxFitPossible;
						correctAndData = (mapID->second);
					    
						/*
							Pair contains:
								vector with: position of errors, length parent and best fit;
								nucleotide which is correct
						*/
						Pair<vector<int>, Dna>  data = correctAndData.i2;
						
						if(data.i1[1]+data.i1[2] > optimalDistance){
							continue;
						}
					}else{
						/*if this is the reverse complement*/
						map<int,Pair<int,Pair<vector<int>, Dna> > >::iterator mapIDrevCompl;
					
						if(idReadsError[s].i1%2!=0){
							mapIDrevCompl = allErrors.find(idReadsError[s].i1-1);
							if(mapIDrevCompl!=allErrors.end()){
							/*the total length minus the position of error, and minus 1 mismatch*/
								int maxFitPossible = length(setReads[idReadsError[s].i1])- positionError - 1;
								int optimalDistance = lenParent + maxFitPossible;
								correctAndData = (mapIDrevCompl->second);
					    
								/*
									Pair contains:
									vector with: position of errors, length parent and best fit;
									nucleotide which is correct
								*/
								Pair<vector<int>, Dna>  data = correctAndData.i2;
						
								if(data.i1[1]+data.i1[2] > optimalDistance){
									continue;
								}
							}
						}
					}
				
					/*iterate over the sequence with error*/
					typedef Iterator<String<Dna> >::Type TIterator;

					int pos_correct = 0;
					int best = 0;
					int smallMismatch = 0;
					int occurBest = 0;
					int idCorr = 0;
					Dna ntCorr; 
					

					for(int r=0; r<idReadsCorrect.size();r++){
						/*if not already correct by the same reads*/
						/*if(mapID!=allErrors.end()){
						
							if(idReadsCorrect[r].i1.i1==correctAndData.i1){	
								continue;
							}
						}*/

						/*the position in the read until which there is the same prefix*/
						int positionCorrect = idReadsCorrect[r].i1.i2+lenParent;
											
						/*count to the position with error*/
						int pos = 0;
						
						/*search the position detect as erroneous at the level of the read*/
						for (TIterator it = begin(setReads[idReadsError[s].i1]); it != end(setReads[idReadsError[s].i1]); ++it){
							/* when we reache the position containig error*/	
							if(pos == positionError){

								/*copy iterator at the position where the error is*/
								TIterator errITer = it;
								/*go to next nucleotide*/
								it++;
								/*compare the nucleotides after the error position*/
								int nbChar = 0;
								int mismatch = 0;
								/*compare until the end of the erroneous read or the corrected one*/
								while(it!=end(setReads[idReadsError[s].i1]) && (positionCorrect+nbChar+1+mismatch)<(length(setReads[idReadsCorrect[r].i1.i1])-1)){
									pos++;
									Dna currentNtErrRead  = setReads[idReadsError[s].i1][pos];
									Dna currentNtCorrRead = setReads[idReadsCorrect[r].i1.i1][positionCorrect+nbChar+1+mismatch];
									/*TODO: the user can choice the number of mismatch*/
									if(currentNtErrRead != currentNtCorrRead){
										mismatch++;
										if(mismatch>acceptMismatch){
											break;
										}
									}else{
										nbChar++;
									}
									it++;
								}
						

							    int position = positionCorrect+nbChar+1;
								/*
									find the appopriate string if the end of the erroneous reads is reached  
									or at the end of the correct read
								*/
								/*I change just +mismatch*/
								if(it == end(setReads[idReadsError[s].i1])|| position==(length(setReads[idReadsCorrect[r].i1.i1]))){							
									/*if more matched characters than the other sibling*/
									if(nbChar+mismatch > best+smallMismatch){
										best = nbChar;
										smallMismatch = mismatch;
										idCorr = idReadsCorrect[r].i1.i1;
										ntCorr = setReads[idReadsCorrect[r].i1.i1][positionCorrect];
										pos_correct = positionCorrect;
										occurBest = idReadsCorrect[r].i2;
									/*if the same number of characters, but the number of occurences is larger*/
									}else if (nbChar+mismatch == best+smallMismatch){// && occurBest < idReadsCorrect[r].i2){
										if(nbChar > best){
											best = nbChar;
											smallMismatch = mismatch;
											idCorr = idReadsCorrect[r].i1.i1;
											ntCorr = setReads[idReadsCorrect[r].i1.i1][positionCorrect];
											pos_correct = positionCorrect;
											occurBest = idReadsCorrect[r].i2;
										}else if(nbChar == best && occurBest < idReadsCorrect[r].i2){
											best = nbChar;
											smallMismatch = mismatch;
											idCorr = idReadsCorrect[r].i1.i1;
											ntCorr = setReads[idReadsCorrect[r].i1.i1][positionCorrect];
											pos_correct = positionCorrect;
											occurBest = idReadsCorrect[r].i2;
										}
									
									}
								}
								break;
							}
							pos++;
						}
					}
			
			
					if(idCorr !=0){
					
						/*the correct nucleotide*/
						Dna nucleotide =  ntCorr;
						
						/*copie of the ID*/
						int idReadErrorPos = idReadsError[s].i1;
						
						/*if so: this is the reverse complement*/
						if((idReadsError[s].i1)%2!=0){
							/*	ID of read give by the user, 
								which is the reverse complement of the identified 
								as erroneous current read
							*/
							idReadErrorPos = idReadsError[s].i1 - 1;

							/*position for correction in the complement - read user*/
							positionError = length(setReads[idReadsError[s].i1]) - positionError - 1; 
							nucleotide = complementaire(nucleotide);
						}
						
						/*search if there is the same */
						map<int,Pair<int,Pair<vector<int>, Dna> > >::iterator mapID;
						mapID = allErrors.find(idReadErrorPos);
	

						if(mapID != allErrors.end()){
			
							Pair<int,Pair<vector<int>, Dna> > correctID = (mapID->second);
							
							Pair<vector<int>,Dna> number = correctID.i2;
							//if(number.i1[0]!=positionError || number.i2 != nucleotide){
								if(number.i1[1]+number.i1[2] + number.i1[3] < lenParent + best + smallMismatch){
								
									allErrors[idReadErrorPos] = dataErroneousNodes(idCorr, positionError, lenParent, best, smallMismatch, occurBest, nucleotide);
									
								}else if(number.i1[1]+number.i1[2] + number.i1[3]==lenParent + best + smallMismatch){
									if(number.i1[1]+number.i1[2]<lenParent + best){
					
										allErrors[idReadErrorPos] = dataErroneousNodes(idCorr, positionError, lenParent, best, smallMismatch, occurBest, nucleotide);
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
		iter++;
	}
	end1 = clock();
	cout << "Time for searching between given levels: "<< (double)(end1-start1)/CLOCKS_PER_SEC << " seconds." << endl;
	
	Pair<StringSet<String<Dna>,Owner < > >, StringSet<String<char> > > fileData = correctErroneousRead(allErrors,fastaComment, setReads);

	cout << "Total corrected reads number is "<< allErrors.size()<<endl;
	
	return fileData;
}




/*GC-content*/
/*fonction which allow to determine the frequency for each nucleotide*/
template < typename TMyIndex, class TSpec >
map<int,double>
determineFrequency(Iter< TMyIndex, VSTree<TSpec> > iter,
		   StringSet<String<Dna>,Owner< > > setReads){

   /*calculate the frequency for each nucleotide*/
   /*'A' = 0, 'C' = 1, 'G' = 2, 'T' = 3*/
   map<int, double> frequency;

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



/*take the total number of lines in the user file*/
int getNbLine(String<char> fileName){

	ifstream reads;
	reads.open(toCString(fileName),ifstream::in);
	string temp;
	string str2="";
	int count=0;
	while(getline(reads,temp)){
		str2 = temp.substr (0,1); 
		if(str2.compare(">")==0 || temp.empty()){
			continue;
		}else{
			count++;
		}
	}
	/*number of lines*/
	return(count);
}

/*read the file format Fasta, which contain all reads*/
Pair<StringSet<String<Dna>,Owner < > >, StringSet<String<char> > > createStringSetFromFile(String<char> Filename){

	cout <<"Take all reads from file"<<endl;
	
	/*contain all reads - input reads and the respective reverse complement*/
	StringSet<String<Dna> > setReads;
	
	/*sava the fasta commen(intformation) for each read*/
	StringSet<String<char> > fastaComment;
	
	Pair<StringSet<String<Dna>,Owner < > >, StringSet<String<char> > > fileData;

	int i = 0;
	int nbLine = getNbLine(Filename);

	/*because also reverse complement*/
	nbLine = nbLine*2;
	resize(setReads,nbLine);
	
	ifstream reads;
	reads.open(toCString(Filename),ifstream::in);
	string temp;
	string str2="";
	int nbComment = 0;
	resize(fastaComment, nbLine);
	/*initialisation with all reads and resp. all reverse complement*/
	while(getline(reads,temp)){
		str2 = temp.substr (0,1); 
		if(str2.compare(">")==0 || temp.empty()){
			fastaComment[nbComment] = temp;
			nbComment++;
			continue;
		}else{
			/*because a probleme with different file format*/
			string a = temp.substr(length(temp)-1,1);
			if(a.compare("A")!=0 && a.compare("T")!=0 && a.compare("C")!=0 && a.compare("G")!=0){
				temp = temp.substr(0,length(temp)-1);
			}
			/* an input read*/
			setReads[i] = temp;

			/*the reverse complement*/
			Dna5StringReverseComplement revCompl(temp);
			setReads[i+1] = revCompl;

			/*go two indice after*/
			i=i+2;
		}
	}
	
	fileData.i1 = setReads;
	fileData.i2 = fastaComment;
	
	return(fileData);
}


/*construction Suffix Array */
Pair<StringSet<String<Dna>,Owner < > >, StringSet<String<char> > > builtSuffixArray(StringSet<String<Dna>,Owner < > > setReads, 
		   double genomeLength, bool poisson, double strictness, int fromLevel, int toLevel,
		   String<char> fileOutPutReadsErrorsCorrected, StringSet<String<char> > fastaComment, int acceptMismatch){
		   
	clock_t start, end1;
	start = clock(); 
	
	cout << "Construction suffix array"<<endl;

	/* Suffix Array constrution with reads*/
	typedef Index< StringSet<String<Dna> > > TMyIndex; 
	TMyIndex myIndex(setReads); 

	/*iterator with restrictions*/
	typedef Iterator< TMyIndex, TopDown<ParentLinks<Postorder> > >::Type TConstrainedIterator; 
	TConstrainedIterator myConstrainedIterator(myIndex); 

	/*calculate the frequency for each nucleotide, didn't use for the moment*/
	/*'A' = 0, 'C' = 1, 'G' = 2, 'T' = 3*/
	map<int,double> frequency;
	frequency = determineFrequency(myConstrainedIterator,setReads);

	/*restrictions just for estimate the genome length if there is no data*/

#ifdef MEDIAN
	int level = fromLevel;
	ofstream out("medianLevels.txt"); 
	ofstream median("medianForEachLevel.txt");
#endif

	if (genomeLength == 0.0)
	{
#ifdef MEDIAN		
		while(level < toLevel){
		//int logRation = (log10(length(setReads)/2))/(log10(4));
		//int l = logRation + 1;
		//cout << l << endl;
			cargo(myIndex).replen_min = level;
			cargo(myIndex).replen_max = level+2; 	
			cargo(myIndex).frequency = frequency;
			double nbOccurence = 0.0; 
			
			median << level <<" "<< medianLevel(myConstrainedIterator) << endl;
			goBegin(myConstrainedIterator);
			while (!atEnd(myConstrainedIterator)){
				if(parentRepLength(myConstrainedIterator)>level){
					nbOccurence = countOccurrences(myConstrainedIterator);
					out << level<< " " << nbOccurence <<endl; 
				}
				myConstrainedIterator++; 
			}
			level+=1;
		}
#else
		//int logRation = (log10(length(setReads)/2))/(log10(4));
		//int l = logRation + 1;
		//cout << l << endl;
		cargo(myIndex).replen_min = fromLevel;
		cargo(myIndex).replen_max = fromLevel+2; 	
		cargo(myIndex).frequency = frequency;

		double expectedValueGivenLevel = medianLevel(myConstrainedIterator);

		/*TODO the length for the reads is always the same, but for real cases it is better to change*/
		/*estimate the genome length by the use of expected value*/
		double readLength = length(setReads[0]);
		double numberReads = length(setReads)/2;

		/* a = readLength - path_label + 1 */
		/*here plus 1 also because the level is between fromLevel and toLevel*/
		double a = readLength - fromLevel + 2;
		double mult = numberReads*a;
		genomeLength = mult/expectedValueGivenLevel;
		cout << "The estimated genome length is " << genomeLength << endl;
#endif
	}
	end1 = clock();
	cout << "#Time required for suffix array construction : "<< (double)(end1-start)/CLOCKS_PER_SEC << " seconds." << "\n\n";
	
	vector<int> idCorrected;
	
	/*restrictions for the searching levels*/
	cargo(myIndex).replen_min = fromLevel;
	cargo(myIndex).replen_max = toLevel; 	
	cargo(myIndex).frequency = frequency;

	/*the core of the correction method*/
	if(poisson){
		/*use of p-value like a limit*/
		return searchNodePoisson(myConstrainedIterator, setReads, fastaComment, genomeLength, strictness, acceptMismatch);
	}else{
		/*use an expected value for a certain level*/
		return searchNodeExpected(myConstrainedIterator, setReads, fastaComment, genomeLength, strictness, acceptMismatch);
	}
}

int main(int argc, char* argv[]) { 

	clock_t start, end;
	
	start = clock();
	
	if(argc > 14){
		cout << "Usage: ./Fiona [options: -f $1 -p $1 -i $1 -l $1 $2 -g $1 -m $1] [input_reads_file] [corrected_reads_file]" << endl;
		cout << "See README.txt for details" << endl;
		exit(1);
	}

	bool programmeOprions = true;

	/*the levels between the search will take place*/
	int fromLevel = 0;
	int toLevel = 0;

	/*the methode that will be use*/
	bool poisson=false;

	double genomeLength = 0.0;

	/*p-value*/
	double strictness = 0.0001;
	
	/**number accepted mismatch - default value 1*/
	int acceptMismatch = 1;

	String<char> file = "";
	String<char> fileOutPutReadsErrorsCorrected = "";

	/*default value*/
	int nbcycle = 3;
	bool expected = false;
	for(int i = 1; i < length(argv); i++){
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
					nbcycle = atoi(argv[i+1]);
					cout << "The number of cycle is set to " << nbcycle << endl;
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
		file = argv[i];
		fileOutPutReadsErrorsCorrected = argv[i+1];
		break;		
	}
	
	
	cout << " " <<endl;
	/*reading the input file and initialisation an array with all reads and their reverse complement*/
	Pair<StringSet<String<Dna>,Owner < > >, StringSet<String<char> > > fileData = createStringSetFromFile(file);
	
	/*reads*/
	StringSet<String<Dna> > setReads = fileData.i1;
	/*the comments from input file for each reads*/
	StringSet<String<char> > fastaComment = fileData.i2;
	
	/*initialise the top and down level by using the log4 from the total number of reads*/
	if(fromLevel==0){
		
		int logRation = (log10(length(setReads)/2.0))/(log10(4.0));
		fromLevel = logRation + 2;
		toLevel   = fromLevel + 10;
		cout << "The top level is " << fromLevel << " and the down level is " << toLevel <<endl;
	}
	cout << " " <<endl;
	int cmp = 0;

	Pair<StringSet<String<Dna>,Owner < > >, StringSet<String<char> > > fileReadsAndComment;
	while(nbcycle > cmp){
		cmp++;
		cout << "Cycle "<< cmp <<" from "<< nbcycle<<endl;
 
		/*if that is last cycle or if nbcycle equal 1*/
		if(cmp==nbcycle){
			fileReadsAndComment = builtSuffixArray(setReads,genomeLength, poisson, strictness, fromLevel, toLevel, fileOutPutReadsErrorsCorrected, fastaComment, acceptMismatch);
			setReads = fileReadsAndComment.i1;
			fastaComment = fileReadsAndComment.i2;
		}else{
			fileReadsAndComment = builtSuffixArray(setReads,genomeLength, poisson, strictness, fromLevel, toLevel, fileOutPutReadsErrorsCorrected, fastaComment, acceptMismatch);
			setReads = fileReadsAndComment.i1;
			fastaComment = fileReadsAndComment.i2;
			
			/*make a new table with the data, because certain are changed*/
			StringSet<String<Dna> > newSetReads;
			resize(newSetReads,length(setReads));
			
			for(int i=0; i < length(setReads); i++){
				newSetReads[i] = setReads[i];
				Dna5StringReverseComplement revCompl(setReads[i]);
				newSetReads[i+1] = revCompl;
				i+=1;
			}
			setReads = newSetReads;
			if(acceptMismatch>0){
				acceptMismatch--;
			} 
			
		}
		/*TODO maybe to stop if there is not reads corrected in the cycle before*/
		/*if so after each iteration must save the ID for the reads which are corrected*/
		/*thus we can also show the total number of reads that are corrected at the final stage*/
	}
	
	/*write in file all input reads with the corrected one*/
	ofstream out(toCString(fileOutPutReadsErrorsCorrected));
	int comment = 0;
	int nbCorrected = 0;
	for(int i=0; i < length(setReads); i++){
		/*to give the number of reads corrected for several iteration*/ 
		if(strContains(toCString(fastaComment[comment]), "corrected")){
			nbCorrected+=1;
		}
		out << fastaComment[comment]<<endl;
		out << setReads[i]<<endl;
		i+=1;
		comment++;
	}
	end = clock();
	cout << " "<< endl;
	if(nbcycle != 1){
		cout << "Total number reads corrected for " <<nbcycle << " cycles is " << nbCorrected <<endl; 
	}
		

	cout << "Time required for execution: "<< (double)(end-start)/CLOCKS_PER_SEC << " seconds." << "\n\n";

	return 0;
}

	
