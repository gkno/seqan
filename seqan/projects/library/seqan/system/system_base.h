/*
 *  system_base.h
 *  SeqAn
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_SYSTEM_BASE_H
#define SEQAN_HEADER_SYSTEM_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////

#ifdef SEQAN_DEBUG

#define SEQAN_DO_SYS(_cond) if (!(_cond)) ::SEQAN_NAMESPACE_MAIN::debug::Message< ::SEQAN_NAMESPACE_MAIN::debug::Check >(__FILE__, __LINE__, #_cond " is FALSE");
#define SEQAN_DO_SYS1(_cond) SEQAN_DO_SYS(_cond)
#define SEQAN_DO_SYS2(_cond, _comment) if (!(_cond)) ::SEQAN_NAMESPACE_MAIN::debug::Error< ::SEQAN_NAMESPACE_MAIN::debug::Check >(__FILE__, __LINE__, _comment);

#else

#define SEQAN_DO_SYS(_cond) { (_cond); }
#define SEQAN_DO_SYS1(_cond) SEQAN_DO_SYS(_cond)
#define SEQAN_DO_SYS2(_cond, _comment) SEQAN_DO_SYS(_cond)

#endif

}

#endif
