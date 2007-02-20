#ifndef SEQAN_HEADER_INDEX_H
#define SEQAN_HEADER_INDEX_H

#include <climits>
#include <functional>
#include <vector>
#include <stack>
#include <queue>
#include <algorithm>
#include <iterator>
#include <utility>
#include <string.h> // memset

//____________________________________________________________________________
// index creation

#include <seqan/sequence.h>
#include <seqan/sequence/sequence_multiple.h>
#include <seqan/pipe.h>
#include <seqan/modifier.h>

#include <seqan/find/find_base.h>

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/index/index_generated_forwards.h>
#endif

#include <seqan/index/index_base.h>

#include <seqan/index/radix.h>
#include <seqan/index/pump_extender3.h>
#include <seqan/index/pipe_merger3.h>
#include <seqan/index/index_skew3.h>

#include <seqan/index/pump_extender7.h>
#include <seqan/index/pipe_merger7.h>
#include <seqan/index/index_skew7.h>

#include <seqan/index/pump_separator7.h>
#include <seqan/index/index_skew7_multi.h>

#include <seqan/index/index_mm.h>
#include <seqan/index/index_sa_btree.h>

#include <seqan/index/pump_lcp_core.h>
#include <seqan/index/index_lcp.h>
#include <seqan/index/index_lcp_tree.h>

#include <seqan/index/index_childtab.h>
#include <seqan/index/index_bwt.h>

#include <seqan/index/shape_base.h>
#include <seqan/index/shape_gapped.h>
#include <seqan/index/index_qgram.h>

//____________________________________________________________________________
// enhanced suffix array (suffix tree)

#include <seqan/index/index_find.h>
#include <seqan/index/index_shims.h>

#include <seqan/index/index_esa_base.h>
#include <seqan/misc/misc_set.h>
#include <seqan/index/index_esa_stree.h>
#include <seqan/index/index_esa_algs.h>



#endif //#ifndef SEQAN_HEADER_...
