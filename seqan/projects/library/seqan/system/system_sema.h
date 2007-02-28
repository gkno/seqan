/*
 *  system_sema.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_HEADER_SYSTEM_SEMAPHORE_H
#define SEQAN_HEADER_SYSTEM_SEMAPHORE_H

namespace SEQAN_NAMESPACE_MAIN
{

#ifdef PLATFORM_WINDOWS

    static SECURITY_ATTRIBUTES SemaphoreDefaultAttributes = {
        sizeof(SECURITY_ATTRIBUTES),
        NULL,
        true
    };

    struct Semaphore
    {
        typedef LONG Type;
        typedef HANDLE Handle;
        enum { MAX_VALUE = MAXLONG };
        
        Handle hSemaphore;

        Semaphore(Type init = 0, Type max = MAX_VALUE) {
            SEQAN_DO_SYS2(hSemaphore = CreateSemaphore(&SemaphoreDefaultAttributes, init, max, NULL), "Could not create Semaphore")
        }

        ~Semaphore() {
            SEQAN_DO_SYS2(CloseHandle(hSemaphore), "Could not destroy Semaphore")
        }

        bool lock(DWORD timeout_millis = INFINITE) {
            return WaitForSingleObject(hSemaphore, timeout_millis) != WAIT_TIMEOUT;
        }

        void unlock() {
            SEQAN_DO_SYS2(ReleaseSemaphore(hSemaphore, 1, NULL), "Could not unlock Semaphore")
        }

    private:

        Semaphore(Semaphore const & origin) {
            // resource copying is not yet supported (performance reason)
            // it needs a reference counting technique
        }

    };

#else

    struct Semaphore
    {
        typedef unsigned int Type;
        typedef sem_t* Handle;
        
        sem_t data, *hSemaphore;

        Semaphore(Type init = 0):
            hSemaphore(&data)
        {
            SEQAN_DO_SYS(!sem_init(hSemaphore, 0, init));
        }

        ~Semaphore() {
            SEQAN_DO_SYS(!sem_destroy(hSemaphore));
        }

        void lock() {
            SEQAN_DO_SYS(!sem_wait(hSemaphore));
        }

        void unlock() {
            SEQAN_DO_SYS(!sem_post(hSemaphore));
        }

    private:

        Semaphore(Semaphore const & origin) {
            // resource copying is not yet supported (performance reason)
            // it needs a reference counting technique
        }

    };


#endif

}

#endif
