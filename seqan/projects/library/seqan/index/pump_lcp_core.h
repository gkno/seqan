/*
 *  lcp_core.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_PUMP_LCP_CORE_H
#define SEQAN_HEADER_PUMP_LCP_CORE_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

    struct LcpConfig
    {
        enum { DefaultWindowSize		= 512 * 1024*1024,	// 512 MB
               DefaultAbsoluteSizes     = true };

		unsigned windowSize;
        bool     absoluteSizes;     // when false, sizes are measured in units of TValue
                                    // when true, sizes are measured in bytes

		LcpConfig():
            windowSize(DefaultWindowSize),
            absoluteSizes(DefaultAbsoluteSizes) {}

		template < typename TValue >
        void absolutize(TValue *) {
            if (!absoluteSizes) return;
            windowSize /= sizeof(TValue);
            if (windowSize < 4096) windowSize = 4096;
        }
    };



	template < typename TTextInput, typename TInvertedSAInput, typename TDest >
	static void lcp_process(TTextInput &textIn, TInvertedSAInput &invertedSAIn, 
						    TDest &dest, LcpConfig conf)
	{
		typedef typename Value<TTextInput>::Type			TValue;
		typedef typename Size<TTextInput>::Type				TSize;
		typedef typename BufReadHandler<TTextInput>::Type	TBufReader;

        SEQAN_PROSET(PRODEPTH, 0);
		TSize rest = length(textIn);
        if (rest < 2) {
            resize(dest, 0);
            return;
        }

        // buffer is a window of textIn
        conf.absolutize((TValue*)NULL);
		TBufReader reader(textIn, conf.windowSize);

        Pair<TSize> out;
		TSize windowBegin = 0;
		TSize overlap = 0;
        TSize _pushes = 0;
        //TSize _olaps = 0;
        //char *seenISA = new bool[n];
        //memset(seenISA, 0, length(invertedSAIn));
		
        #ifdef SEQAN_DEBUG_INDEX
            TSize n = rest;
            TSize lcpMax = 0, lcpAvrg = 0, lcpNumer = 0, sigma = 1;	// for lcpMax, lcpMean, |Sigma|
        #endif

        resize(dest, rest);
		beginWrite(dest);
		// read the first block of textIn
		typename TBufReader::Buffer	buffer = reader.first();
//        typename Value<TBufReader>::Type buffer = reader.first();
		while (length(buffer))
        {
            SEQAN_PROADD(PRODEPTH, 1);
            SEQAN_PROMARK("LCP-Durchlauf beginnen");
            beginRead(textIn);
			beginRead(invertedSAIn);

            out.i2 = 0;                                             // out.i2 == h
            TSize windowSize = length(buffer);
			TSize windowEnd = windowBegin + windowSize;
			TSize newOverlap = windowEnd;
            #ifdef SEQAN_DEBUG_INDEX
                printf("  read window[%x,%x) overlay %x rest %x\n", windowBegin, windowEnd, overlap, rest);
            #endif
            rest -= windowSize;

            //unsigned _pops = 0;

            while (!eof(invertedSAIn)) {

				if ((*invertedSAIn).i1 != 0) {

					TSize left = (*invertedSAIn).i2[1] + out.i2;        // left == j + h

					if (overlap <= (*invertedSAIn).i2[1] && left <= windowEnd) {
	//                        printf("* ");
						for(; left < windowBegin; ++left) {
							++textIn;
							++out.i2;
						}

						TSize k = left - windowBegin;                   // k is relative to window
						while (k < windowSize && (!eof(textIn)) && *textIn == buffer[k]) {
							++textIn;
							++k;
						}
						TSize hAdd = k - (left - windowBegin);
						out.i2 += hAdd;                                 // increase h by successful compares
						left += hAdd;

						if (k != windowSize || rest == 0) {
							out.i1 = (*invertedSAIn).i1 - 1;
							push(dest, out);
							//SEQAN_ASSERT(!seenISA[out.i1] && 0 <= out.i1 && out.i1 < n);
							//seen[out.i1] = true;
							++_pushes;

                            #ifdef SEQAN_DEBUG_INDEX
                                if ((lcpNumer += out.i2) > n) {
                                    lcpNumer -= n;
                                    ++lcpAvrg;
                                }
                                if (lcpMax < out.i2) lcpMax = out.i2;
								if (!out.i2) ++sigma;
                            #endif
						}
					}

					if (left >= windowEnd && (*invertedSAIn).i2[1] < newOverlap) {
						#ifdef SEQAN_VERBOSE
							printf("crossing border @ %d len %d\n", (*invertedSAIn).i2[1], out.i2);
						#endif
						newOverlap = (*invertedSAIn).i2[1];    //++_olaps;
					}
				}

				if (out.i2) --out.i2;                               // if (h > 0) h = h - 1
				else ++textIn;

                ++invertedSAIn; //++_pops;
			}

			endRead(invertedSAIn);
			endRead(textIn);

			windowBegin = windowEnd;
			overlap = newOverlap;
			buffer = reader.next();                                 // read the following block

			//printf("pops:%d pushes:%d overlaps:%d\n", _pops, _pushes, _olaps);
		}
		// trailing zero
		push(dest, Pair<TSize>(length(textIn) - 1, 0));

		//printf("pushes:%d length:%d\n", _pushes, length(textIn));
        //SEQAN_ASSERT(_pushes == length(textIn));
        //for (unsigned i = 0; i < n; ++i)
        //    if (!seen[i])
        //        printf("___%d______HUH?\n",i);
        //delete[] seenISA;

        #ifdef SEQAN_DEBUG_INDEX
            ::std::cout << "  n: " << n;
            ::std::cout << "  lcpMax: " << lcpMax;
            ::std::cout << "  lcpAvrg: " << lcpAvrg + lcpNumer / (double)n;
            ::std::cout << "  sigma: " << sigma << ::std::endl;
        #endif

//TODO: uncomment this line
//		reader.end();
		endWrite(dest);
        SEQAN_PROSET(PRODEPTH, 0);
    }

	template < typename TTextInput, typename TInvertedSAInput, typename TDest >
	static inline void lcp_process(TTextInput &textIn, TInvertedSAInput &invertedSAIn, TDest &dest)
	{
		lcp_process(textIn, invertedSAIn, dest, LcpConfig());
	}

//////////////////////////////////////////////////////////////////////////////


	template < typename TTextInput, typename TLimitsString, typename TInvertedSAInput, typename TDest >
	static void lcp_process_multi(
		TTextInput &textIn, TLimitsString const &limits, 
		TInvertedSAInput &invertedSAIn, 
		TDest &dest, LcpConfig conf)
	{
		typedef typename Value<TTextInput>::Type			TValue;
		typedef typename Size<TTextInput>::Type				TSize;
		typedef typename BufReadHandler<TTextInput>::Type	TBufReader;

		typedef typename Value<TInvertedSAInput>::Type::T2::T	TPair;

        SEQAN_PROSET(PRODEPTH, 0);
		TSize rest = length(textIn);
        if (rest < 2) {
            resize(dest, 0);
            return;
        }

        // buffer is a window of textIn
        conf.absolutize((TValue*)NULL);
		TBufReader reader(textIn, conf.windowSize);

        Pair<TSize> out;
		TSize windowBegin = 0;
		TSize overlap = 0;
        TSize _pushes = 0;
        //TSize _olaps = 0;
        //char *seenISA = new bool[n];
        //memset(seenISA, 0, length(invertedSAIn));
		
        #ifdef SEQAN_DEBUG_INDEX
            TSize n = rest;
            TSize lcpMax = 0, lcpAvrg = 0, lcpNumer = 0, sigma = 1;	// for lcpMax, lcpMean, |Sigma|
        #endif

        resize(dest, rest);
		beginWrite(dest);
		// read the first block of textIn
		typename TBufReader::Buffer	buffer = reader.first();
//        typename Value<TBufReader>::Type buffer = reader.first();

		while (length(buffer))
        {
            SEQAN_PROADD(PRODEPTH, 1);
            SEQAN_PROMARK("LCP-Durchlauf beginnen");

            beginRead(textIn);
			beginRead(invertedSAIn);

            out.i2 = 0;                                             // out.i2 == h
            TSize windowSize = length(buffer);
			TSize windowEnd = windowBegin + windowSize;
			TSize newOverlap = windowEnd;
            #ifdef SEQAN_DEBUG_INDEX
                printf("  read window[%x,%x) overlay %x rest %x\n", windowBegin, windowEnd, overlap, rest);
            #endif
            rest -= windowSize;

            //unsigned _pops = 0;

			TSize leftOrig, left;		// begin and begin + seen_lcp of the lower string (global)

            while (!eof(invertedSAIn)) {

				if ((*invertedSAIn).i1 != 0) {

					leftOrig = posGlobalize((*invertedSAIn).i2[1], limits);
					left = leftOrig + out.i2;        // left == j + h

					if (overlap <= leftOrig && left <= windowEnd) {
	//                        printf("* ");
						for(; left < windowBegin; ++left) {
							++textIn;
							++out.i2;
						}

						// end of lower string relative to window
						TSize kEnd = limits[getValueI1((*invertedSAIn).i2[1]) + 1] - windowBegin;
						if (kEnd > windowSize) kEnd = windowSize;

						TSize k = left - windowBegin;						// k is relative to window

						// comparision will break before the end of the upper string
						// -> only lower string needs a clipping

						while (k < kEnd && (!eof(textIn)) && *textIn == buffer[k]) {
							++textIn;
							++k;
						}
						TSize hAdd = k - (left - windowBegin);
						out.i2 += hAdd;                                 // increase h by successful compares
						left += hAdd;

						if (k != windowSize || rest == 0) {
							out.i1 = (*invertedSAIn).i1 - 1;
							push(dest, out);
							//SEQAN_ASSERT(!seenISA[out.i1] && 0 <= out.i1 && out.i1 < n);
							//seen[out.i1] = true;
							++_pushes;

                            #ifdef SEQAN_DEBUG_INDEX
                                if ((lcpNumer += out.i2) > n) {
                                    lcpNumer -= n;
                                    ++lcpAvrg;
                                }
                                if (lcpMax < out.i2) lcpMax = out.i2;
								if (!out.i2) ++sigma;
                            #endif
						}
					}

					if (left >= windowEnd && leftOrig < newOverlap) {
						#ifdef SEQAN_VERBOSE
							printf("crossing border @ %d len %d\n", (*invertedSAIn).i2[1], out.i2);
						#endif
						newOverlap = leftOrig;    //++_olaps;
					}
				}

				if (out.i2) 
					--out.i2;                               // if (h > 0) h = h - 1
				else 
					++textIn; 

                ++invertedSAIn; //++_pops;
			}

			endRead(invertedSAIn);
			endRead(textIn);

			windowBegin = windowEnd;
			overlap = newOverlap;
			buffer = reader.next();                                 // read the following block

			//printf("pops:%d pushes:%d overlaps:%d\n", _pops, _pushes, _olaps);
		}
		// trailing zero
		push(dest, Pair<TSize>(length(textIn) - 1, 0));

		//printf("pushes:%d length:%d\n", _pushes, length(textIn));
        //SEQAN_ASSERT(_pushes == length(textIn));
        //for (unsigned i = 0; i < n; ++i)
        //    if (!seen[i])
        //        printf("___%d______HUH?\n",i);
        //delete[] seenISA;

        #ifdef SEQAN_DEBUG_INDEX
            ::std::cout << "  n: " << n;
            ::std::cout << "  lcpMax: " << lcpMax;
            ::std::cout << "  lcpAvrg: " << lcpAvrg + lcpNumer / (double)n;
            ::std::cout << "  sigma: " << sigma << ::std::endl;
        #endif

//TODO: uncomment this line
//		reader.end();
		endWrite(dest);
        SEQAN_PROSET(PRODEPTH, 0);
    }

	template < typename TTextInput, typename TLimitsString, typename TInvertedSAInput, typename TDest >
	static void lcp_process_multi(
		TTextInput &textIn, TLimitsString const &limits, 
		TInvertedSAInput &invertedSAIn, 
		TDest &dest)
	{
		lcp_process_multi(textIn, limits, invertedSAIn, dest, LcpConfig());
	}



//}

}

#endif
