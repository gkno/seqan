//============================================================================
// SeqAn - The Library for Sequence Analysis
// http://www.seqan.de 
//
//============================================================================
// $Author$
// $Revision$
// $Date$
//============================================================================

#ifndef SEQAN_HEADER_SEQUENCE_H
#define SEQAN_HEADER_SEQUENCE_H

//____________________________________________________________________________
// prerequisites

#include <seqan/basic.h>

//____________________________________________________________________________

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/sequence/sequence_generated_forwards.h>
#endif

#include <seqan/sequence/sequence_interface.h>
#include <seqan/sequence/lexical.h>

//____________________________________________________________________________
// segments (suffix, ...)

#include <seqan/sequence/segment_base.h>
#include <seqan/sequence/segment_infix.h>
#include <seqan/sequence/segment_suffix.h>
#include <seqan/sequence/segment_prefix.h>

//____________________________________________________________________________
// strings

#include <seqan/sequence/string_base.h>
#include <seqan/sequence/string_pointer.h>
#include <seqan/sequence/string_alloc.h>
#include <seqan/sequence/string_array.h>
#include <seqan/sequence/string_cstyle.h>
#include <seqan/sequence/string_stack.h>
#include <seqan/sequence/string_packed.h>
#include <seqan/sequence/string_value_expand.h>

#include <seqan/sequence/std_string.h>

#include <seqan/sequence/sequence_multiple.h>
#include <seqan/sequence/sequence_shortcuts.h>

#endif //#ifndef SEQAN_HEADER_...
