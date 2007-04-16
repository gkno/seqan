#ifndef SEQAN_HEADER_ALIGN_H
#define SEQAN_HEADER_ALIGN_H

//____________________________________________________________________________
// prerequisites

#include <seqan/sequence.h>
#include <seqan/score.h>

//____________________________________________________________________________

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/align/align_generated_forwards.h>
#endif

#include <seqan/align/gaps_base.h>
#include <seqan/align/gaps_iterator_base.h>

#include <seqan/align/gaps_array.h>

#include <seqan/align/align_base.h>
#include <seqan/align/align_cols_base.h>
#include <seqan/align/align_iterator_base.h>

#include <seqan/align/matrix_base.h>
#include <seqan/misc/priority_type_base.h>
#include <seqan/misc/priority_type_heap.h>

#include <seqan/align/align_dynprog.h>
#include <seqan/align/align_local_dynprog.h>

#endif //#ifndef SEQAN_HEADER_...
