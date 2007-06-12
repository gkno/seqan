/*	2006 by Hendrik Woehrle */
#ifndef SEQAN_HEADER_CHAINING
#define SEQAN_HEADER_CHAINING

//#include <iostream>
//#include <sstream>
//#include <iomanip>
//#include <limits>
//#include <typeinfo>
//#include <vector>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/misc/misc_random.h>

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/chaining/chaining_generated_forwards.h>
#endif

//Skip List

#include <seqan/chaining/skip_list_type.h>
#include <seqan/chaining/skip_list_base.h>
#include <seqan/chaining/skip_list_iterator.h>
#include <seqan/chaining/skip_base_element.h>
#include <seqan/chaining/skip_element.h>
#include <seqan/chaining/skip_list_impl.h>
#include <seqan/chaining/skip_list_dynamic.h>

//Range Tree

#include <seqan/chaining/rt_base.h>
#include <seqan/chaining/rt_skip_element.h>
#include <seqan/chaining/rt_skip_base_element.h>
#include <seqan/chaining/rt_sl_impl.h>
#include <seqan/chaining/rt_sl_base.h>
#include <seqan/chaining/rt_base2.h>
#include <seqan/chaining/rt_common_algos.h>
#include <seqan/chaining/rt_impl.h>
#include <seqan/chaining/rt_sl_def_algos.h>
#include <seqan/chaining/rt_sl_compl_algos.h>

#endif // SEQAN_HEADER_SKIP_LIST
