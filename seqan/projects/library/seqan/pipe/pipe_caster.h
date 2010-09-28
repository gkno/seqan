 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id$
 ==========================================================================*/

#ifndef SEQAN_HEADER_PIPE_CASTER_H
#define SEQAN_HEADER_PIPE_CASTER_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

	struct CasterReinterpret;
	struct CasterConvert;
    
    template < typename TValue, typename TSpec = CasterReinterpret >
    struct Caster;

	template < typename TInput, typename TValue, typename TSpec >
    struct Value< Pipe< TInput, Caster<TValue, TSpec> > > {
		typedef TValue Type;
	};


/**
.Spec.Caster:
..cat:Pipelining
..general:Class.Pipe
..summary:Casts the input type in a specific output type.
..signature:Pipe<TInput, Caster<TValue[, TSpec]> >
..param.TInput:The type of the pipeline module this module reads from.
..param.TValue:The new output type.
..param.TSpec:$CasterReinterpret$ (default) or $CasterConvert$.
..remarks: The input stream is casted using $reinterpret_cast<TValue>$.
..include:seqan/pipe.h
*/

    //////////////////////////////////////////////////////////////////////////////
    // caster pipe
    template <typename TInput, typename TValue >
    struct Pipe< TInput, Caster<TValue, CasterReinterpret> >
    {
		TInput      &in;
        
        Pipe(TInput& _in):
            in(_in) {}
        
        inline TValue const & operator*() const {
            return reinterpret_cast<TValue const &>(*in);
        }

        Pipe& operator++() {
            ++in;
            return *this;
        }
	};
    
    template <typename TInput, typename TValue >
    struct Pipe< TInput, Caster<TValue, CasterConvert> >
    {
		TInput      &in;
        
        Pipe(TInput& _in):
            in(_in) {}
        
        inline TValue operator*() const {
            return TValue(*in);
        }

        Pipe& operator++() {
            ++in;
            return *this;
        }
	};
//}

}

#endif
