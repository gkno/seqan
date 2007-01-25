/*
 *  extender7.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_PUMP_EXTENDER7_H
#define SEQAN_HEADER_PUMP_EXTENDER7_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

    template <typename TCompression = void>
    struct Extender7;

    template <typename TTextInput, typename TNameInput, typename TCompression >
    struct Pipe< Bundle2< TTextInput, TNameInput >, Extender7<TCompression> >
    {
        typedef typename Size<Pipe>::Type           SizeType;
        typedef typename Value<TTextInput>::Type    TextInType;
        typedef typename Value<TNameInput>::Type    NameInType;

        typedef Tuple<TextInType, 4, TCompression>  X4Tuple;
        typedef Tuple<TextInType, 5, TCompression>  X5Tuple;
        typedef Tuple<TextInType, 6, TCompression>  X6Tuple;
        typedef Tuple<typename NameInType::T2, 3>   NTuple;
        typedef Triple<SizeType, NTuple, X6Tuple, Compressed>   OutType0;
        typedef Triple<SizeType, NTuple, X6Tuple, Compressed>   OutType3;
        typedef Triple<SizeType, NTuple, X4Tuple, Compressed>   OutType5;
        typedef Triple<SizeType, NTuple, X5Tuple, Compressed>   OutType6;
        typedef Triple<SizeType, NTuple, X6Tuple, Compressed>   OutType124;

        // pipeline interfaces to ease specialization
        typedef Pipe< void, AbstractSource< OutType0, SizeType > > Out0;
        typedef Pipe< void, AbstractSource< OutType3, SizeType > > Out3;
        typedef Pipe< void, AbstractSource< OutType5, SizeType > > Out5;
        typedef Pipe< void, AbstractSource< OutType6, SizeType > > Out6;
        typedef Pipe< void, AbstractSource< OutType124, SizeType > > Out124;

        //Out0 out0;
        //Out3 out3;
        //Out5 out5;
        //Out6 out6;
        //Out124 out124;
        
        Pipe(Bundle2< TTextInput, TNameInput > &_in)
            //out0(_textIn.length() / 7),
            //out3((_textIn.length() + 4) / 7),
            //out5((_textIn.length() + 2) / 7),
            //out6((_textIn.length() + 1) / 7),
            //out124(_nameIn.length())
        {
            addListener(_in); 
        }
    };

	// little copy helper
	// which is needed for bitpacked structs (no =sign possible)
	template <typename Dest, typename Ofs, typename Src>
	static finline Src const __cp_(Dest &dst, Ofs const ofs, Src const src) {
		return dst.i3.assignAt(ofs, src);
	}
        
    template < typename TTextInput, typename TNameInput,
               typename TOut0, typename TOut3, typename TOut5, typename TOut6, typename TOut124 >
    static bool skew7_extend(TTextInput &textIn, TNameInput &nameIn,
                             TOut0 &out0, TOut3 &out3, TOut5 &out5, TOut6 &out6, TOut124 &out124)
    {
        resize(out0, length(textIn) / 7);
        resize(out3, (length(textIn) + 4) / 7);
        resize(out5, (length(textIn) + 2) / 7);
        resize(out6, (length(textIn) + 1) / 7);
        resize(out124, length(nameIn));
        if (!(
            beginRead(textIn) &&
            beginRead(nameIn) &&
            beginWrite(out0) &&
            beginWrite(out3) &&
            beginWrite(out5) &&
            beginWrite(out6) &&
            beginWrite(out124))) return false;

		typename Value<TOut0>::Type  o0;
		typename Value<TOut3>::Type  o3;
		typename Value<TOut5>::Type  o5;
		typename Value<TOut6>::Type  o6;
		typename Value<TOut124>::Type o1, o2, o4;

        typename Size<TTextInput>::Type p = length(textIn);
        unsigned r = (unsigned)(p % 7);


        // BEGIN I: PREFILL

        switch (r) {
        case 6:
/* 6 */                                                                    __cp_(o6,0,    *textIn); ++textIn; o6.i1 = p--;

        case 5:
/* 5 */                                                         __cp_(o5,0,__cp_(o6,1,    *textIn)); ++textIn; o5.i1 = p--;

        case 4:
/* 4 */                                                 o4.i2[0] = o5.i2[0] = o6.i2[0] = (*nameIn).i2; ++nameIn; o4.i1 = p--;
                                                     __cp_(o4,0,__cp_(o5,1,__cp_(o6,2,    *textIn))); ++textIn;
            
        case 3:
/* 3 */                                   __cp_(o3,0,__cp_(o4,1,__cp_(o5,2,__cp_(o6,3,   *textIn)))); ++textIn; o3.i1 = p--;
                
        case 2:
/* 2 */                           o2.i2[0] = o3.i2[0] = o4.i2[1] = o5.i2[1] = o6.i2[1] = (*nameIn).i2; ++nameIn; o2.i1 = p--;
                               __cp_(o2,0,__cp_(o3,1,__cp_(o4,2,__cp_(o5,3,__cp_(o6,4,    *textIn))))); ++textIn;
                
        case 1:
/* 1 */                o1.i2[0] = o2.i2[1] = o3.i2[1] = o4.i2[2] = o5.i2[2] = o6.i2[2] = (*nameIn).i2; ++nameIn; o1.i1 = p--;
                    __cp_(o1,0,__cp_(o2,1,__cp_(o3,2,                                     *textIn))); ++textIn;
            if (r >= 6) push(out6, o6);
            if (r >= 5) push(out5, o5);
            if (r >= 4) push(out124, o4);

        case 0:;
        }

            // BEGIN II: PREFILL AND PUSH FULLY FILLED TRIPLES

        if (!eof(nameIn)) {
/* 0 */  __cp_(o0,0,__cp_(o1,1,__cp_(o2,2,__cp_(o3,3,                                     *textIn)))); ++textIn; o0.i1 = p--;
                
/* 6 */  __cp_(o0,1,__cp_(o1,2,__cp_(o2,3,__cp_(o3,4,                      __cp_(o6,0,    *textIn))))); ++textIn; o6.i1 = p--;
                
/* 5 */  __cp_(o0,2,__cp_(o1,3,__cp_(o2,4,__cp_(o3,5,           __cp_(o5,0,__cp_(o6,1,    *textIn)))))); ++textIn; o5.i1 = p--;

/* 4 */     o0.i2[0] = o1.i2[1] = o2.i2[2] = o3.i2[2] = o4.i2[0] = o5.i2[0] = o6.i2[0] = (*nameIn).i2; ++nameIn; o4.i1 = p--;
         __cp_(o0,3,__cp_(o1,4,                      __cp_(o4,0,__cp_(o5,1,__cp_(o6,2,    *textIn))))); ++textIn;
            if (r >= 3) push(out3, o3);
            if (r >= 2) push(out124, o2);
                
/* 3 */  __cp_(o0,4,__cp_(o1,5,           __cp_(o3,0,__cp_(o4,1,__cp_(o5,2,__cp_(o6,3,    *textIn)))))); ++textIn; o3.i1 = p--;
                
/* 2 */     o0.i2[1] = o1.i2[2] = o2.i2[0] = o3.i2[0] = o4.i2[1] = o5.i2[1] = o6.i2[1] = (*nameIn).i2; ++nameIn; o2.i1 = p--;
         __cp_(o0,5,           __cp_(o2,0,__cp_(o3,1,__cp_(o4,2,__cp_(o5,3,__cp_(o6,4,    *textIn)))))); ++textIn;
            if (r >= 1) push(out124, o1);
                
/* 1 */     o0.i2[2] = o1.i2[0] = o2.i2[1] = o3.i2[1] = o4.i2[2] = o5.i2[2] = o6.i2[2] = (*nameIn).i2; ++nameIn; o1.i1 = p--;
                    __cp_(o1,0,__cp_(o2,1,__cp_(o3,2,                                     *textIn))); ++textIn;
            push(out0, o0);
            push(out6, o6);
            push(out5, o5);
			push(out124, o4);
            r = 7;
        }

        // MAIN LOOP: PUSH FULLY FILLED TRIPLES
        
        while (!eof(nameIn)) {
/* 0 */  __cp_(o0,0,__cp_(o1,1,__cp_(o2,2,__cp_(o3,3,                                     *textIn)))); ++textIn; o0.i1 -= 7;
                
/* 6 */  __cp_(o0,1,__cp_(o1,2,__cp_(o2,3,__cp_(o3,4,                      __cp_(o6,0,    *textIn))))); ++textIn; o6.i1 -= 7;
                
/* 5 */  __cp_(o0,2,__cp_(o1,3,__cp_(o2,4,__cp_(o3,5,           __cp_(o5,0,__cp_(o6,1,    *textIn)))))); ++textIn; o5.i1 -= 7;

/* 4 */     o0.i2[0] = o1.i2[1] = o2.i2[2] = o3.i2[2] = o4.i2[0] = o5.i2[0] = o6.i2[0] = (*nameIn).i2; ++nameIn; o4.i1 -= 7;
         __cp_(o0,3,__cp_(o1,4,                      __cp_(o4,0,__cp_(o5,1,__cp_(o6,2,    *textIn))))); ++textIn;
            push(out3, o3);
            push(out124, o2);
                
/* 3 */  __cp_(o0,4,__cp_(o1,5,           __cp_(o3,0,__cp_(o4,1,__cp_(o5,2,__cp_(o6,3,    *textIn)))))); ++textIn; o3.i1 -= 7;
                
/* 2 */     o0.i2[1] = o1.i2[2] = o2.i2[0] = o3.i2[0] = o4.i2[1] = o5.i2[1] = o6.i2[1] = (*nameIn).i2; ++nameIn; o2.i1 -= 7;
         __cp_(o0,5,           __cp_(o2,0,__cp_(o3,1,__cp_(o4,2,__cp_(o5,3,__cp_(o6,4,    *textIn)))))); ++textIn;
            push(out124, o1);
                
/* 1 */     o0.i2[2] = o1.i2[0] = o2.i2[1] = o3.i2[1] = o4.i2[2] = o5.i2[2] = o6.i2[2] = (*nameIn).i2; ++nameIn; o1.i1 -= 7;
                    __cp_(o1,0,__cp_(o2,1,__cp_(o3,2,                                     *textIn))); ++textIn;
            push(out0, o0);
            push(out6, o6);
            push(out5, o5);
            push(out124, o4);
        }

        // END: FLUSH PARTIALLY FILLED TRIPLES

        {
/* 0 */             __cp_(o1,1,__cp_(o2,2,__cp_(o3,3,   0)));
                
/* 6 */             __cp_(o1,2,__cp_(o2,3,__cp_(o3,4,   0)));

/* 5 */             __cp_(o1,3,__cp_(o2,4,__cp_(o3,5,   0)));

/* 4 */                o1.i2[1] = o2.i2[2] = o3.i2[2] = 0;
                    __cp_(o1,4,                         0);
            if (r >= 3) push(out3, o3);
            if (r >= 2) push(out124, o2);

/* 3 */             __cp_(o1,5,                         0);

/* 2 */                o1.i2[2] =						0;
            if (r >= 1) push(out124, o1);
        
        }
                
        endWrite(out124);
        endWrite(out6);
        endWrite(out5);
        endWrite(out3);
        endWrite(out0);
        endRead(nameIn);
        endRead(textIn);
        return false;
    }
    
//}

}

#endif
