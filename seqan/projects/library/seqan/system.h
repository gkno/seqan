#ifndef SEQAN_HEADER_SYSTEM_H
#define SEQAN_HEADER_SYSTEM_H

//____________________________________________________________________________
// prerequisites

#include <seqan/platform.h>

#include <cstdio>
#include <ctime>
#include <string>
#include <iostream>

#ifdef PLATFORM_WINDOWS
# include <windows.h>
#else
# include <climits>
# include <pthread.h>
# include <errno.h>
# include <semaphore.h>
# include <aio.h>
#endif

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/system/system_manual_forwards.h>
#include <seqan/system/file_manual_forwards.h>
#endif

//____________________________________________________________________________
// multi-threading

#include <seqan/system/system_base.h>
#include <seqan/system/system_mutex.h>
#include <seqan/system/system_sema.h>
#include <seqan/system/system_event.h>
#include <seqan/system/system_thread.h>

//____________________________________________________________________________
// synchronous and asynchronous files

#include <seqan/system/file_sync.h>
#include <seqan/system/file_async.h>

#endif //#ifndef SEQAN_HEADER_...
