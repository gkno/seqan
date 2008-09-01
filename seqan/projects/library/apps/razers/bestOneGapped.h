#ifndef SEQAN_HEADER_BEST_ONE_GAPPED_DP_H
#define SEQAN_HEADER_BEST_ONE_GAPPED_DP_H

#include <bitset>

using namespace std;

namespace SEQAN_NAMESPACE_MAIN
{

// helpers

char best_shape_folder[200];
bool best_shape_helpFolder = false;

unsigned char shape_countBits[256];


//////////////////////////////////////////////////////////////////////////////
// Returns the estimated minimum coverage of a shape with weight q, span s at threshold t
template<typename TValueQ, typename TValueS, typename TValueT>
inline TValueS getMinCov(TValueQ q, TValueS s, TValueT t)
{
	TValueS mincov;
	if(t > s - q + 1){
		mincov = q + 2 * (t - 1) - (t - (s - q + 1));
	}
	else mincov = q + 2 * (t - 1);

	return mincov;
}

//////////////////////////////////////////////////////////////////////////////
// Returns the sum of two probability values in log space
template<typename TLog>
TLog
logAdd (TLog a, TLog b)
{
	if(isinf(a)) return b;
	if(isinf(b)) return a;
	if(isnan(a+log1p(exp(b-a)))) return a;
	return (a + log1p(exp(b-a)));
}


//////////////////////////////////////////////////////////////////////////////
// Returns log probability of q-gram-configuration q ending at position pos in sequence
template<typename TErrDistr>
typename Value<TErrDistr>::Type _getShapeProb(unsigned q, unsigned pos, TErrDistr & logErrorDistr, unsigned span)
{

	typename Value<TErrDistr>::Type prob = 0.0;
	for(unsigned j = 0; j < span; ++j)
	{
		if (q & 1) prob += (logErrorDistr[pos-span+j+1]);
		else prob += (logErrorDistr[(pos+1) + pos-span+j+1]);
		q = q >> 1;
	}
//	if(prob < PREC) return 0.0;
//	else 
	return prob;
	
}

void
initCountBits()
{	// Count 1-bits in a byte
	for(unsigned i = 0; i < 256; ++i) {
		unsigned char bits = 0;
		for(unsigned char b = 0; b < 8; ++b)
			if ((i >> b) & 1) ++bits;
		shape_countBits[i] = bits;
	}
}
//long double PREC =
String<bool> firstTimeK;


template <typename TLossString, typename TErrorDistr, typename TLossRString>
void 
best_shape(unsigned maxN, TLossString &found, TErrorDistr & logErrorDistr, TLossRString & lossRates, unsigned maxQ, unsigned minQ, unsigned minGap, unsigned maxGap, unsigned maxE, unsigned minE, unsigned maxT, bool optionHammingOnly){

	unsigned shape = 0;
	typedef typename Value<TLossString>::Type TLoss;
	typedef typename Value<TErrorDistr>::Type TErrorVal;

//	TLoss lossRate;
//	FILE *datei;

	bool editDistance = false;
//	bool editDistance = true;
	if(length(logErrorDistr)==4*maxN)
		editDistance = true;

	// remembers the loss rate for shape with maximal minimum Coverage for k errors, for each bucket for a weight q
	String<String<TLoss> > loss;

	// remembers the minimum Coverage for shape with maximal minimum Coverage for k errors, for each bucket, for a weight q
	String<String<unsigned> > minCoverage;

	// remembers the current maximal minimum Coverage for k errors, for each bucket, for a weight q
	String<String<unsigned> > currentMinCov;

	// remembers t for shape with maximal minimum Coverage for k errors, for each bucket, for a weight q
	String<String<unsigned> > threshold;
	
	// remembers k for shape with maximal minimum Coverage for each bucket, for a weight q
	String<String<unsigned> > error;

	// remembers the shape with maximal mininimum Coverage for k errors, for each bucket, for a weight q
	String<String<unsigned> > shapes;

	// remembers the index for the shape with maximal minimum Coverage for k errors, for each bucket, for a weight q
	String<String<int> > counter;


	// helpers
	unsigned gminCov;
	unsigned possQ = 0;

	// amount of loss rates
	unsigned lossRatesl =  (unsigned)(lossRates.size()) +1;

	resize(loss, lossRatesl);
	resize(minCoverage, lossRatesl);
	resize(threshold, lossRatesl);
	resize(counter, lossRatesl);
	resize(error, lossRatesl);
	resize(shapes, lossRatesl);
	resize(currentMinCov, lossRatesl);

	fill(firstTimeK,maxE*lossRatesl,true);
//	for(unsigned i = 0; i < maxE; ++i) firstTimeK[i]=true;


	// for each weight q...
	for(unsigned q = minQ; q <= maxQ; ++q){

		// amount of possible shapes with weight q and size maxS		
		unsigned possQgram = ((maxGap-minGap)*(q-1) +1); 
	
		// amount of possible shapes with weight q, size maxS, maxE		
		possQ = (maxGap*(q-1) +1)* maxE;
		
		for(unsigned i = 0; i < lossRatesl; ++i){
			resize(loss[i], possQ);
			resize(minCoverage[i], possQ);
			resize(threshold[i], possQ);
			resize(error[i], possQ);
			resize(currentMinCov[i], possQ);
			resize(counter[i], possQ);
			resize(shapes[i], possQ);
			
		}


		// j = span of shape
		for(unsigned j = q+minGap; j <= q+maxGap /*|| j < q+q-1 */; ++j){

			// set current minimum coverage 0, counter -1 for each e
			for(unsigned a = 0; a < lossRatesl; ++a){
				String<int> underCounter2 = counter[a];
				String<unsigned> underCurrMinCov = currentMinCov[a];
				for(unsigned e = 0; e < maxE; ++e){
				underCurrMinCov[e] = 0;
				underCounter2[e] =  -1;
				}
				counter[a] = underCounter2;
				currentMinCov[a] = underCurrMinCov;
			}
				

			unsigned gaplength = j - q;
			gminCov = 0;
			// k = position of gap
			for(unsigned k = 1; k < (q/2); ++k){
				
				// gap > 0
				if(j > q){
					unsigned last = q - k;
					shape = (1 << k)-1;
					shape =  shape << (last + gaplength);
					shape = shape | ((1 << (last)) -1);
				}
				//  gap = 0
				else{
					shape = (1 << q)-1;
					k = minQ -1;
				}

				if(editDistance)
				{
					std::cout << "Do edit DP...\n";
					CharString shapeString;
					fill(shapeString,j,'0');
					for(unsigned pos = 0; pos < j; ++pos)
						if((shape & (1 << pos)) > 0)
							shapeString[pos] = '1';

					String< State<long double> > states;
					initPatterns(states, shapeString, maxE-1, logErrorDistr, optionHammingOnly, true);
					computeFilteringLoss(found, states, j, maxT, maxE, logErrorDistr, false, true);
				}
				else
				{// compute loss rate for shape
					std::cout << "Do Hamming DP...\n";
					doDP(maxN, found, logErrorDistr, shape, j, maxE, maxT);
				}

				// go through found and find loss rate
				for(unsigned e = minE; e < maxE; ++e) {
					bool highestOptimalFound = false;
					for(unsigned t = maxT-1; t > 0; --t) {
						TErrorVal lossrate = 1.0 - (TErrorVal) _transformBack(found[e*maxT+t]);
						
						typename std::map<TErrorVal, unsigned>::iterator it, itlow, itup, itend, itbegin;
						if(lossrate <= 0.0){
							if(highestOptimalFound) continue;
							else highestOptimalFound = true;
						}
						// points to first item in loss rate file
						itbegin = lossRates.begin();
						if(lossrate < ((*itbegin).first)){ 
							continue;
						}
						itend = lossRates.end();
						// points to last item in loss rate file
						--itend;
						
						if(lossrate > ((*itend).first)){ 
							continue;
						}
						
						TErrorVal helpLoss = (*itend).first; 
						if(lossrate > helpLoss){
							continue;
						}
						
						
						itup = lossRates.upper_bound(lossrate);
						itlow = itup;
						--itlow;

						gminCov = getMinCov(q, j, t);

						unsigned index = unsigned((*itlow).second);		
						String<unsigned> underCurrMin = currentMinCov[index];
						String<int> underCounter = counter[index];
						String<unsigned> underShapes = shapes[index];
						String<long double> underLoss = loss[index];
						String<unsigned> underThreshold = threshold[index];
						String<unsigned> underError = error[index];
						String<unsigned> underMinCoverage = minCoverage[index];
						int currentCounter = underCounter[e];
						//if(gminCov >= underCurrMin[e]){
						//	if(gminCov > underCurrMin[e]){
						//		currentCounter = 0;
						//		underCurrMin[e] = gminCov;
						//		
						//	}
						//	else if(gminCov == underCurrMin[e]){
								++currentCounter;
						//	}
	
							// rememeber the shape (with current maximal minimum coverage) and the appropriate parameters 
							currentMinCov[index] = underCurrMin;
							underLoss[currentCounter + (e * possQgram)] = lossrate;
							underShapes[currentCounter + (e * possQgram)] = shape;
							underThreshold[currentCounter + (e * possQgram)] = t;
							underError[currentCounter + (e * possQgram)] = e;
							underMinCoverage[currentCounter + (e * possQgram)] = gminCov;
							
							underCounter[e] = currentCounter;
							counter[index] = underCounter;
							loss[index] = underLoss;
							shapes[index] = underShapes;
							threshold[index] = underThreshold;
							error[index] = underError;
							minCoverage[index] = underMinCoverage;
						//}
						
					} // t-loop
				} //e-loop
					
			}// k-loop

			for(unsigned x = 0; x < lossRatesl; ++x){
				String<int> underCount = counter[x];
				for(unsigned e = minE; e < maxE; ++e){
					// if a best shape exists at e errors --> write shape in file
					if(underCount[e] >= 0){
					std::cout << "printing...\n";
					print_shape(q, possQgram, maxN, underCount[e], e, x, shapes, loss, error, threshold, minCoverage, maxQ+maxGap, lossRates, maxE, editDistance);
					}
				}
			}

		} // j-loop
				
	}// i-loop

}// best_shape-loop

				

// write best-shapes in data file
template<typename TShapeString, typename TLString, typename TErrString, typename TThreshString, typename TMinCovString, typename TLossRString>
void print_shape(unsigned , unsigned possQgram, unsigned nMax, int counter, unsigned e, unsigned bucket, TShapeString & shapes, TLString & loss, TErrString & error, TThreshString & threshold, TMinCovString & minCoverage, unsigned , TLossRString & lossRates, unsigned maxE, bool editDist){

	typedef typename Value<TErrString>::Type TErrStr;
	typedef typename Value<TShapeString>::Type TShapeStr;
	typedef typename Value<TShapeStr>::Type TShape;
	typedef typename Value<TLString>::Type TLStr;
	typedef typename Value<TThreshString>::Type TThreshStr;
	typedef typename Value<TMinCovString>::Type TMinCovStr;
	typedef typename Value<TLossRString>::Type TLossRStr;
	typedef typename std::map<long double, unsigned>::iterator TBucketIt;

	
	String<String<char *> > lossInString;
	resize(lossInString, int(lossRates.size()));
	
// 	TErrStr& underError = error[bucket];
// 	TShapeStr& underShape = shapes[bucket];
// 	TLStr& underLoss = loss[bucket];
// 	TThreshStr& underThresh = threshold[bucket];
// 	TMinCovStr& underMinCov = minCoverage[bucket];

	// write the best shape with the appropriate parameters in file
	for(int currentCounter = 0; currentCounter <= counter; ++currentCounter){

		TErrStr underError = error[bucket];
		TShapeStr underShape = shapes[bucket];
		TLStr underLoss = loss[bucket];
		TThreshStr underThresh = threshold[bucket];
		TMinCovStr underMinCov = minCoverage[bucket];

		int currIndex = currentCounter + (e * possQgram);
		
		TShape shape = underShape[currIndex];
		CharString shapeString;
		fill(shapeString,50,'0');

		for(int pos = 31; pos >= 0; --pos)
			if((shape & (1 << pos)) == (1 << pos))
				shapeString[pos] = '1';
		unsigned last = 49;
		while(shapeString[last] != '1')
			--last;
		resize(shapeString,last+1);


		TBucketIt it, itup, itlow;
		long double lossrate = underLoss[currIndex]; 
		itup = lossRates.upper_bound(lossrate);
		itlow = itup;
		--itlow;

		// create the whole file name
		stringstream datName;
		if(best_shape_helpFolder) datName << best_shape_folder;
		else datName << "gapped_params";
		datName << "/shapes_" << nMax << "_" << underError[currIndex] << "_";
		if(editDist) datName << "L_";
		else datName <<"H_";
		datName << (*itlow).first << "_" << (*itup).first << ".dat";
		

		// if datName-file doesnt exist, write the title on it
		if(firstTimeK[maxE*bucket + e]==true){
			firstTimeK[maxE*bucket + e] = false;
			ofstream fout(datName.str().c_str(), ios::out);
			fout << "shape\t\tt\t\tloss rate\t\tminCoverage\n\n";
			fout.close();
		}
		// write best shape with its properties on the file
		ofstream fout(datName.str().c_str(), ios::app | ios::out);
		fout << shapeString << "\t\t";
		fout << underThresh[currIndex] << "\t\t";
		fout << underLoss[currIndex] << "\t\t";
		fout << underMinCov[currIndex] << endl; 
		fout.close();

	}
}

		
	
template <typename TLossString, typename TErrorDistr>
void doDP(unsigned maxN, TLossString &found, TErrorDistr & logErrorDistr, unsigned shape, unsigned span, unsigned maxE, unsigned maxT)
{
	typedef TLossString TMatrixCol;
	typedef typename Value<TErrorDistr>::Type TProbValue;
	typedef String<int> TIntCol;

	unsigned QPot = 1 << span;

	//compute span(Q) and possible q-grams
	unsigned possibleQ = 0;
	for (unsigned e = 0; e < maxE; e++){
	possibleQ += binomialCoefficient(span,e);
	}


	// errors in q (with new number)
	TIntCol error_Col;
	resize(error_Col, possibleQ);

	// new index for q
	TIntCol newNr_Col;
	fill(newNr_Col, QPot, -1);

	// new index of shift-q after shift(0)-operation
	TIntCol shift0_Col;
	resize(shift0_Col, possibleQ);
	
	// new index of shift-q after shift(1)-operation
	TIntCol shift1_Col;
	resize(shift1_Col, possibleQ);

	// true if errors of (q & shape) <= errorsPerQgram
	String<bool> matched_Col;
	fill(matched_Col, possibleQ, false);

	// probability for each possible q with < maxE
	TMatrixCol qProb;
	resize(qProb, possibleQ);


	// fill error_Col, newNr_Col
	unsigned nr = 0;
	for(unsigned q = 0; q < QPot; q++){
		unsigned errors = 
					shape_countBits[q & 255] + shape_countBits[(q >> 8) & 255] + 
					shape_countBits[(q >> 16) & 255] + shape_countBits[(q >> 24) & 255];

		if(errors < maxE){
			newNr_Col[q] = nr;
			error_Col[nr] = errors;
			++nr;
		}

	}


	// fill shift0_Col, shift1_Col, matched_Col
	unsigned i = 0;
	for(unsigned q = 0; q < QPot; q++){
		if(newNr_Col[q] >= 0){

			unsigned _q = (q << 1) & (QPot - 1);
			shift0_Col[i] = newNr_Col[(_q|0)];
			shift1_Col[i] = newNr_Col[(_q|1)];
			
			unsigned qshape = q & shape;
			unsigned errors = 
				shape_countBits[qshape & 255] + shape_countBits[(qshape >> 8) & 255] +
				shape_countBits[(qshape >> 16) & 255] + shape_countBits[(qshape >> 24) & 255];

			if(errors <= 0/*errorsPerQGram*/){
				matched_Col[i] = true;
				
			}
			qProb[i] = _getShapeProb(q,maxN-1,logErrorDistr, span);
			++i;
		}
	}

	clear(newNr_Col);
	


	// columns n-1 and n for recursion 
	TMatrixCol col0;
	TMatrixCol col1;
	fill(col0, maxE * possibleQ * maxT, log(0.0));
	resize(col1, maxE * possibleQ * maxT);

	// remembers sum of probabilities of all possible sequences having
	// exactly e errors (recursively updated during DP)
	TMatrixCol count_col;
	fill(count_col, maxE, log(0.0));

	//initialized with first position (can be either 0 or 1)
	count_col[0] = logErrorDistr[maxN];

	if(maxE > 1) count_col[1] = logErrorDistr[0];
	
	// same thing for last column i.e. sequence length maxN
	// (last column gets special treatment)
	TMatrixCol count_colFinal;
	fill(count_colFinal, maxE, log(0.0));

	

	// RECURSION BEGIN
	for(unsigned q = 0; q < possibleQ; ++q) 
	{
		col0[q*maxT] = log(1.0);
		// every bit in q set to 1 marks an error

		// for n==0
		if(matched_Col[q]){
			// we miss no match if q has no or one error
			// --> probability is 1.0
			col0[q*maxT+1] = log(1.0);
			for(unsigned t = 2; t < maxT; ++t)
				col0[q*maxT+t] = log(0.0);	
		}
		else{
			// we miss 1 match for t>0 and q has more than <errorsPerQGram> errors
			// --> probability for finding the q-gram is 0
			for(unsigned t = 1; t < maxT; ++t)
					col0[q*maxT+t] = log(0.0);
				
		}
	}
//	dump(col0,0);

	// iterate over sequence length n
	TMatrixCol *col = &col1;
	TMatrixCol *colPrev = &col0;

	
/*	cout << "2:0";
	dump(col0, 0);
	cout << " :1";
	dump(col0, 1);
*/
	
	// RECURSION
	//
	// found(n,q,t,e) = (1-errorProb[n-Q]) * found(n-1,0|(q>>1),t-delta,e) delta=1/0 <-> q hat 0/>0 fehler
	//               + errorProb[n-Q] * found(n-1,1|(q>>1),t-delta,e-1)
	
	// rekursion (fuer q-gram matches <=1 fehler)
	// found(n,q,t,e) = (1-errorProb[n-Q]) * found(n-1,0|(q>>1),t-delta,e) delta=1/0 <-> q hat <=1/>1 fehler
	//               + errorProb[n-Q] * found(n-1,1|(q>>1),t-delta,e-1)
	
	for(unsigned n = span; n < maxN; ++n)
	{
		if(n>span) // if n-span==0 count_col is already up to date (as this is how it was initialized)
		{	
			for(int e = maxE - 1; e > 0; --e)
			{
				// update with position n-Q: can be either 0 
				count_col[e] = ((logErrorDistr[maxN+n-span]) + count_col[e]);
				// or 1
				count_col[e] = logAdd(count_col[e],logErrorDistr[n-span] + count_col[e-1]);
				
			}
			// for 0 errors, updating with 0 is the only possibility
			count_col[0] += (logErrorDistr[maxN+n-span]);
		}

		for(unsigned e = 0; e < maxE * possibleQ; e += possibleQ)
		{
		
			for(unsigned q = 0; q < possibleQ; ++q)
			{
				
				for(unsigned t = 0; t < maxT; ++t)
				{
		
					unsigned _t = t;
					if (_t > 0 && matched_Col[q])
						--_t;
					
					// col[errors][qgram][threshold]
					// = col[(errors*QPot+q)*maxT+_t]


					long double recovered = ((logErrorDistr[maxN+n-span]) + (*colPrev)[(e+ shift0_Col[q])*maxT+_t]);
					if(e > 0) {
						if (shift1_Col[q] != -1) recovered = logAdd(recovered,(logErrorDistr[n-span] + (*colPrev)[((e-possibleQ)+ shift1_Col[q])*maxT+_t]));
						else recovered = logAdd(recovered, (logErrorDistr[n-span] + log(0.0)));
					}
					//if(recovered < PREC) recovered = 0.0;

					//if this is the last iteration (i.e. last columnof matrix), add log probability for the whole q gram
					if(n==maxN-1) recovered += qProb[q];
				
		
					(*col)[(e+q)*maxT+t] = recovered;
	
				}
				
				// if this is the last iteration, add log probability of whole q-gram to the appropriate entry in count_col
				// k�nnte man auch erst sp�ter machen...
				if(e ==  ((maxE-1) * possibleQ) && n == maxN - 1)
				{
					unsigned errors = error_Col[q];

					long double prob = qProb[q];
					if(errors < maxE) 
						for(unsigned eSum = errors; eSum < maxE; ++eSum)
						{
							count_colFinal[eSum] = logAdd(count_colFinal[eSum], (prob + count_col[eSum-errors]));
							
						}
				}
			}

		}

		TMatrixCol *tmp = col;
		col = colPrev;
		colPrev = tmp;

/*		cout << n+1<<":0";
		dump(*colPrev, 0);
		cout << " :1";
		dump(*colPrev, 1);
		cout << " :2";
		dump(*colPrev, 2);
		*/	
	}


/*	for(unsigned i = 0; i < maxE; ++i)
		cout << count_colFinal[i] << "\t";
	cout << "\n\n";
*/
	// RECURSION END
	for(unsigned eSum = 0; eSum < maxE; ++eSum)
		for(unsigned t = 0; t < maxT; ++t) 
		{
			long double recovered = log(0.0);
			for(unsigned q = 0; q < possibleQ; ++q) 
			{

				unsigned errors = error_Col[q];
				if (errors <= eSum) {
					unsigned e = eSum - errors;
					recovered = logAdd(recovered,(*colPrev)[(e*possibleQ+q)*maxT+t]);
				}
			}
			//divide by probabilitiy of seeing eSum errors
			found[eSum*maxT+t] = recovered - count_colFinal[eSum];
			
		}

		
}

}
#endif
