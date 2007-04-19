#ifndef SEQAN_HEADER_FIND_H
#define SEQAN_HEADER_FIND_H

//____________________________________________________________________________
// prerequisites

#include <seqan/sequence.h>
#include <seqan/score.h>

//____________________________________________________________________________
// exact pattern matching

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/find/find_generated_forwards.h>
#endif

#include <seqan/find/find_base.h>

#include <seqan/find/find_horspool.h>
#include <seqan/find/find_shiftand.h>
#include <seqan/find/find_shiftor.h>
#include <seqan/find/find_bndm.h>
#include <seqan/find/find_quasar.h>

//____________________________________________________________________________
// approximate pattern matching

#include <seqan/find/find_score.h>
#include <seqan/find/find_myers_ukkonen.h>

#include <seqan/graph.h>
#include <seqan/find/find_bom.h>
#include <seqan/find/find_ahocorasick.h>
#include <seqan/find/find_multiple_shiftand.h>
#include <seqan/find/find_set_horspool.h>


//#include <seqan/find/find_multi.h>

//??? ToDo ???
//#include <seqan/find/find_wumanber.h>

#endif //#ifndef SEQAN_HEADER_...
