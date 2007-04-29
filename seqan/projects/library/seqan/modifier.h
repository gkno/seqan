#ifndef SEQAN_HEADER_MODIFIER_H
#define SEQAN_HEADER_MODIFIER_H

//____________________________________________________________________________
// prerequisites

#include <functional>

//____________________________________________________________________________
// basics

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/modifier/modifier_generated_forwards.h>
#endif

#include <seqan/sequence.h>
#include <seqan/modifier/modifier_iterator.h>
#include <seqan/modifier/modifier_string.h>

//____________________________________________________________________________
// applications

#include <seqan/modifier/modifier_functors.h>
#include <seqan/modifier/modifier_view.h>
#include <seqan/modifier/modifier_reverse.h>
#include <seqan/modifier/modifier_shortcuts.h>


#endif //#ifndef SEQAN_HEADER_...
