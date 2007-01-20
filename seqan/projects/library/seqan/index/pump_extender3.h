/*
 *  extender3.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_PUMP_EXTENDER3_H
#define SEQAN_HEADER_PUMP_EXTENDER3_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

    struct Extender3;

    template < typename TTextInput, typename TNameInput >
    struct Pipe< Bundle2< TTextInput, TNameInput >, Extender3 >
    {
		enum { maxShift = 2 };
		typedef typename Size<Pipe>::Type           SizeType;
        typedef typename Value<TTextInput>::Type    TextInType;
        typedef typename Value<TNameInput>::Type    NameInType;

        typedef Tuple<TextInType, maxShift> XTuple;
        typedef Tuple<typename NameInType::T2, maxShift> NTuple;
        typedef Triple<SizeType, NTuple, XTuple, Compressed> OutType0;
        typedef Triple<SizeType, NTuple, XTuple, Compressed> OutType12;

        // pipeline interfaces to ease specialization
        typedef Pipe< void, AbstractSource< OutType0, SizeType > > Out0;
        typedef Pipe< void, AbstractSource< OutType12, SizeType > > Out12;
        //Out0 out0;
        //Out12 out12;
        
        Pipe(Bundle2< TTextInput, TNameInput > &_in)
            //out0(_textIn.size() / 3),
            //out12(_nameIn.size())
        {
            addListener(_in); 
        }
    };
        
    template < typename TTextInput, typename TNameInput, typename TOut0, typename TOut12 >
    static bool skew3_extend(TTextInput &textIn, TNameInput &nameIn, TOut0 &out0, TOut12 &out12)
    {
        resize(out0, length(textIn) / 3);
        resize(out12, length(nameIn));
        if (!(
            beginRead(textIn) && 
            beginRead(nameIn) &&
            beginWrite(out0) && 
            beginWrite(out12))) return false;
		
		typename Value<TOut0>::Type  o0;
		typename Value<TOut12>::Type o1, o2;

        unsigned r = (unsigned)(length(textIn) % 3);
        bool filled = (r != 0);
        if (r == 2) {
            o2.i1 = (*nameIn).i1;
            o2.i2[0] = (*nameIn).i2; ++nameIn;
            o2.i3[0] = *textIn; ++textIn;
        }
            
        if (r >= 1) {
            o1.i1 = (*nameIn).i1;
            o1.i2[0] = (*nameIn).i2; ++nameIn;
            o1.i3[0] = *textIn; ++textIn;
        }
        
        if (r == 2) {
            o2.i2[1] = o1.i2[0];
            o2.i3[1] = o1.i3[0];
            push(out12, o2);
        }
        
        while (!eof(nameIn)) {
            o1.i3[1] = o0.i3[0] = *textIn; ++textIn;
            o2.i3[0] = o0.i3[1] = *textIn; ++textIn;
            
            o2.i1 = (*nameIn).i1;
            o0.i2[0] = o2.i2[0] = o1.i2[1] = (*nameIn).i2; ++nameIn;
            o2.i3[1] = *textIn; ++textIn;
            if (filled)
                push(out12, o1);
            else
                filled = true;
            
            o1.i1 = (*nameIn).i1;
            o0.i2[1] = o2.i2[1] = o1.i2[0] = (*nameIn).i2; ++nameIn;
            o1.i3[0] = o2.i3[1];
            push(out12, o2);
            o0.i1 = o2.i1 + 1;
            push(out0, o0);
        }

        o1.i2[1] = 0;
        o1.i3[1] = 0;
        if (filled) push(out12, o1);

        endWrite(out12);
        endWrite(out0);
        endRead(nameIn);
        endRead(textIn);
        return true;
    }
    
//}

}

#endif
