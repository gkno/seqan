#ifndef SEQAN_HEADER_SYSTEM_H
#define SEQAN_HEADER_SYSTEM_H

#include <cstdio>
#include <ctime>
#include <string>
#include <iostream>

#ifdef PLATFORM_WINDOWS
# include <windows.h>
#else
# include <limits>
# include <pthread.h>
# include <errno.h>
# include <semaphore.h>
# include <aio.h>
#endif


#include <seqan/system/system_base.h>
#include <seqan/system/system_mutex.h>
#include <seqan/system/system_sema.h>
#include <seqan/system/system_event.h>
#include <seqan/system/system_thread.h>
//#include <seqan/system/system_ioqueue.h>

#include <seqan/system/system_page.h>
#include <seqan/system/system_page_raid0.h>

#include <seqan/system/system_profile.h>

#endif //#ifndef SEQAN_HEADER_...
