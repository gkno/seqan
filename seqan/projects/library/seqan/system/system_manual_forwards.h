#ifndef SEQAN_HEADER_SYSTEM_MANUAL_FORWARDS_H 
#define SEQAN_HEADER_SYSTEM_MANUAL_FORWARDS_H 

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

//////////////////////////////////////////////////////////////////////////////
// CLASSES
//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN {

//____________________________________________________________________________
// Event

struct Event;       	// "projects/library/seqan/system/system_event.h"(18)

//____________________________________________________________________________
// Mutex

struct Mutex;       	// "projects/library/seqan/system/system_mutex.h"(16)

//____________________________________________________________________________
// Semaphore

struct Semaphore;       	// "projects/library/seqan/system/system_sema.h"(16)

//____________________________________________________________________________
// Thread

template <typename Worker> struct Thread;       	// "projects/library/seqan/system/system_thread.h"(18)


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//____________________________________________________________________________
// close

inline bool close(Event &e);       	// "projects/library/seqan/system/system_event.h"(109)
inline bool close(Mutex &m);       	// "projects/library/seqan/system/system_mutex.h"(79)
template <typename TWorker> inline bool close(Thread<TWorker> &m);       	// "projects/library/seqan/system/system_thread.h"(97)

//____________________________________________________________________________
// kill

template <typename TWorker> inline bool kill(Thread<TWorker> &m);       	// "projects/library/seqan/system/system_thread.h"(102)

//____________________________________________________________________________
// lock

inline bool lock(Mutex &m);       	// "projects/library/seqan/system/system_mutex.h"(83)

//____________________________________________________________________________
// open

inline bool open(Event &e, bool initial);       	// "projects/library/seqan/system/system_event.h"(101)
inline bool open(Event &e);       	// "projects/library/seqan/system/system_event.h"(105)
inline bool open(Mutex &m, bool initial);       	// "projects/library/seqan/system/system_mutex.h"(71)
inline bool open(Mutex &m);       	// "projects/library/seqan/system/system_mutex.h"(75)
template <typename TWorker> inline bool open(Thread<TWorker> &m);       	// "projects/library/seqan/system/system_thread.h"(87)

//____________________________________________________________________________
// run

template <typename TWorker> inline bool run(Thread<TWorker> &m);       	// "projects/library/seqan/system/system_thread.h"(92)

//____________________________________________________________________________
// signal

inline bool signal(Event &e);       	// "projects/library/seqan/system/system_event.h"(131)

//____________________________________________________________________________
// unlock

inline bool unlock(Mutex &m);       	// "projects/library/seqan/system/system_mutex.h"(87)

//____________________________________________________________________________
// waitFor

inline bool waitFor(Event &e);       	// "projects/library/seqan/system/system_event.h"(113)
template <typename TTime > inline bool waitFor(Event &e, TTime timeout_millis);       	// "projects/library/seqan/system/system_event.h"(118)
template <typename TWorker> inline bool waitFor(Thread<TWorker> &m);       	// "projects/library/seqan/system/system_thread.h"(107)

} //namespace SEQAN_NAMESPACE_MAIN


#endif

