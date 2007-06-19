/*
 *  system_thread.h
 *  SeqAn
 *
 *  Created by David Weese on 17.07.05.
 *
 */

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_HEADER_SYSTEM_THREAD_H
#define SEQAN_HEADER_SYSTEM_THREAD_H

namespace SEQAN_NAMESPACE_MAIN
{

#ifdef PLATFORM_WINDOWS

    static SECURITY_ATTRIBUTES ThreadDefaultAttributes = {
        sizeof(SECURITY_ATTRIBUTES),
        NULL,
        true
    };

    template <typename Worker>
    struct Thread
    {
        typedef HANDLE Handle;

        Handle hThread;
        DWORD  hThreadID;
        Worker worker;

        Thread() {}

        template <typename TArg>
        Thread(TArg &arg):
            worker(arg) {}

        ~Thread() {
            if (*this) {
                cancel();
                wait();
            }
        }

        inline bool open(BOOL initital = false) {
            return hThread = CreateThread(
                &ThreadDefaultAttributes     // default security attributes 
                0,                           // use default stack size  
                &_start,                     // thread function 
                this,                        // argument to thread function 
                0,                           // use default creation flags 
                &hThreadID);                 // returns the thread identifier 
        }

        inline bool close() {
            return CloseHandle(hThread) && !(hThread = NULL);
        }

        inline bool cancel(DWORD exitCode = 0) {
            return !TerminateThread(hThread, exitCode);
        }

        inline bool wait(DWORD timeout_millis = INFINITE) {
            return WaitForSingleObject(hMutex, timeout_millis) != WAIT_TIMEOUT;
        }

        inline operator bool() const {
            return hThread != NULL;
        }

    private:

        Thread(Thread const &) {
            // resource copying is not yet supported (performance reason)
            // it needs a reference counting technique
        }

        static DWORD WINAPI _start(LPVOID _this) {
            reinterpret_cast<Thread*>(_this)->worker.run(&reinterpret_cast<Thread*>(_this));
			return 0;	// return value should indicate success/failure
        }
    };
    
#else

    template <typename Worker>
    struct Thread
    {
        typedef pthread_t* Handle;

        pthread_t data, *hThread;
        Worker worker;

        Thread() {}

        template <typename TArg>
        Thread(TArg &arg):
            worker(arg) {}

        ~Thread() {
            if (*this) {
                cancel();
                wait();
            }
        }

        inline bool open()
        {
            if (!pthread_create(&data, NULL, _start, this) && (hThread = &data)) {
                return true;
            } else
                return false;
        }

        inline bool close() {
            return cancel() && wait() && !(hThread == NULL);
        }

        inline bool cancel() {
            return !(pthread_cancel(data));
        }

        inline bool wait() {
            return !(pthread_join(data, NULL));
        }

        inline bool wait(void* &retVal) {
            return !(pthread_join(data, &retVal));
        }

        inline bool detach() {
            return !(pthread_detach(data));
        }

        inline operator bool() const {
            return hThread != NULL;
        }

    private:

        Thread(Thread const &) {
            // resource copying is not yet supported (performance reason)
            // it needs a reference counting technique
        }

        static void* _start(void* _this) {
            reinterpret_cast<Thread*>(_this)->worker.run(&reinterpret_cast<Thread*>(_this));
			return 0;
        }
    };
    
#endif


	//////////////////////////////////////////////////////////////////////////////
	// global thread functions

	template <typename TWorker>
	inline bool open(Thread<TWorker> &m) {
		return m.open();
	}

	template <typename TWorker>
	inline bool run(Thread<TWorker> &m) {
		return m.open();
	}

	template <typename TWorker>
	inline bool close(Thread<TWorker> &m) {
		return m.close();
	}

	template <typename TWorker>
	inline bool kill(Thread<TWorker> &m) {
		return m.close();
	}

	template <typename TWorker>
	inline bool waitFor(Thread<TWorker> &m) {
		return m.wait();
	}

}

#endif
