#ifndef SEQAN_HEADER_PIPE_H
#define SEQAN_HEADER_PIPE_H

//____________________________________________________________________________
// prerequisites

#include <seqan/file.h>
#include <seqan/basic/basic_volatile_ptr.h>

#include <cstdio>
#include <cassert>
#include <functional>
#include <iterator>
#include <climits>
#include <vector>
#include <queue>

//____________________________________________________________________________
// pipes

#define SEQAN_NAMESPACE_PIPELINING pipe

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/pipe/pipe_generated_forwards.h>
#endif

#include <seqan/pipe/pipe_base.h>
#include <seqan/pipe/pipe_iterator.h>
#include <seqan/pipe/pipe_caster.h>
#include <seqan/pipe/pipe_counter.h>
#include <seqan/pipe/pipe_echoer.h>
#include <seqan/pipe/pipe_edit_environment.h>
#include <seqan/pipe/pipe_filter.h>
#include <seqan/pipe/pipe_joiner.h>
#include <seqan/pipe/pipe_namer.h>
#include <seqan/pipe/pipe_sampler.h>
#include <seqan/pipe/pipe_shifter.h>
#include <seqan/pipe/pipe_source.h>
#include <seqan/pipe/pipe_tupler.h>

//____________________________________________________________________________
// pools

#include <seqan/pipe/pool_base.h>
#include <seqan/pipe/pool_mapper.h>

#include <seqan/misc/priority_type_base.h>
#include <seqan/misc/priority_type_heap.h>
#include <seqan/pipe/pool_sorter.h>

#endif //#ifndef SEQAN_HEADER_...
