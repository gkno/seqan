/*
 *  system_ioqueue.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_SYSTEM_IOQUEUE_H
#define SEQAN_HEADER_SYSTEM_IOQUEUE_H

// THIS FILE IS CURRENTLY NOT USED

namespace SEQAN_NAMESPACE_MAIN
{

    //////////////////////////////////////////////////////////////////////////////
    //
    //  Doubly linked list structure.  Can be used as either a list head, or
    //  as link words.
    //

    typedef struct _LIST_ENTRY {
        struct _LIST_ENTRY *Flink;
        struct _LIST_ENTRY *Blink;
    } LIST_ENTRY, *PLIST_ENTRY;

    #ifndef InitializeListHead
    //
    //  VOID
    //  InitializeListHead(
    //      PLIST_ENTRY ListHead
    //      );
    //

    #define InitializeListHead(ListHead) (\
        (ListHead)->Flink = (ListHead)->Blink = (ListHead))

    //
    //  BOOLEAN
    //  IsListEmpty(
    //      PLIST_ENTRY ListHead
    //      );
    //

    #define IsListEmpty(ListHead) \
        ((ListHead)->Flink == (ListHead))

    //
    //  PLIST_ENTRY
    //  RemoveHeadList(
    //      PLIST_ENTRY ListHead
    //      );
    //

    #define ListFirst(ListHead) \
        ((ListHead)->Flink)

    //
    //  PLIST_ENTRY
    //  RemoveHeadList(
    //      PLIST_ENTRY ListHead
    //      );
    //

    #define RemoveHeadList(ListHead) \
        (ListHead)->Flink;\
        {RemoveEntryList((ListHead)->Flink)}

    //
    //  PLIST_ENTRY
    //  RemoveTailList(
    //      PLIST_ENTRY ListHead
    //      );
    //

    #define RemoveTailList(ListHead) \
        (ListHead)->Blink;\
        {RemoveEntryList((ListHead)->Blink)}

    //
    //  VOID
    //  RemoveEntryList(
    //      PLIST_ENTRY Entry
    //      );
    //

    #define RemoveEntryList(Entry) {\
        PLIST_ENTRY _EX_Blink;\
        PLIST_ENTRY _EX_Flink;\
        _EX_Flink = (Entry)->Flink;\
        _EX_Blink = (Entry)->Blink;\
        _EX_Blink->Flink = _EX_Flink;\
        _EX_Flink->Blink = _EX_Blink;\
        }

    //
    //  VOID
    //  InsertTailList(
    //      PLIST_ENTRY ListHead,
    //      PLIST_ENTRY Entry
    //      );
    //

    #define InsertTailList(ListHead,Entry) {\
        PLIST_ENTRY _EX_Blink;\
        PLIST_ENTRY _EX_ListHead;\
        _EX_ListHead = (ListHead);\
        _EX_Blink = _EX_ListHead->Blink;\
        (Entry)->Flink = _EX_ListHead;\
        (Entry)->Blink = _EX_Blink;\
        _EX_Blink->Flink = (Entry);\
        _EX_ListHead->Blink = (Entry);\
        }

    //
    //  VOID
    //  InsertHeadList(
    //      PLIST_ENTRY ListHead,
    //      PLIST_ENTRY Entry
    //      );
    //

    #define InsertHeadList(ListHead,Entry) {\
        PLIST_ENTRY _EX_Flink;\
        PLIST_ENTRY _EX_ListHead;\
        _EX_ListHead = (ListHead);\
        _EX_Flink = _EX_ListHead->Flink;\
        (Entry)->Flink = _EX_Flink;\
        (Entry)->Blink = _EX_ListHead;\
        _EX_Flink->Blink = (Entry);\
        _EX_ListHead->Flink = (Entry);\
        }

    #endif //InitializeListHead


#ifdef PLATFORM_WINDOWS

    struct IOQueue
    {
        typedef HANDLE Handle;
        typedef ULONGLONG FilePtr;

        typedef void aHint;
        typedef void aCallback(aHint*);
		typedef Event aEvent;

        typedef struct _IOREQUEST : public OVERLAPPED, LIST_ENTRY {
            HANDLE      hFile;
            LPVOID      lpBuffer;
            DWORD       dwTransferCount;
            aCallback   *Callback;
            aHint       *Hint;
            bool        bRead;
            IOQueue     *me;
            void print() {
                if (bRead)
                    printf("READ OPERATION (%x)\n", (ULONGLONG)this);
                else
                    printf("WRITE OPERATION (%x)\n", (ULONGLONG)this);
                printf("handle:\t\t%x\n",   (ULONGLONG)hFile);
                printf("buffer:\t\t%x\n",   (ULONGLONG)lpBuffer);
                printf("fileptr:\t%x:%x\n", OffsetHigh, Offset);
                printf("size:\t\t%x\n",     dwTransferCount);
            }
        } IORequest, *PIORequest;

        typedef PIORequest aRequest;

        LIST_ENTRY  IORequestList;
        LIST_ENTRY  IOPendingList;
        LIST_ENTRY  IODoneList;
        
        Mutex       qLock;                   // locks queue pointers
        Semaphore   qSemaphore;              // counts requests open in the queue
        Event       qEmpty;                  // signals emptiness of the queue
        Handle      wThread;
        DWORD       wThreadID;
        DWORD       maxRequests;

        IOQueue(DWORD _maxRequests = 0):
            maxRequests(_maxRequests),
            qSemaphore(0, (_maxRequests)? _maxRequests: Semaphore::MAX_VALUE),
            qEmpty(true)
        {
            InitializeListHead(&IORequestList);
            InitializeListHead(&IOPendingList);
            InitializeListHead(&IODoneList);

            for(DWORD i = 0; i < maxRequests; i++)
                InsertTailList(&IODoneList, newRequest());
            
            wThread = CreateThread(
                NULL,                        // default security attributes 
                0,                           // use default stack size  
                &IOWorker,                   // thread function 
                this,                        // argument to thread function 
                0,                           // use default creation flags 
                &wThreadID);                 // returns the thread identifier 
                        
            // Check the return value for success. 
            
            SEQAN_DO_SYS(wThread != INVALID_HANDLE_VALUE);
        }

        ~IOQueue()
        {
            // wait until IOWorker took his hands off the lists
            qLock.lock();
            CloseHandle(wThread);
            qLock.unlock();

            if (!IsListEmpty(&IORequestList)) {
                printf("ERROR: Queue not empty while destructing ioqueue class.");
			    printList(&IORequestList);
                deleteList(&IORequestList);
            }

            if (!IsListEmpty(&IOPendingList)) {
                printf("ERROR: Requests are pending while destructing ioqueue class.");
			    printList(&IOPendingList);
                deleteList(&IOPendingList);
            }

            deleteList(&IODoneList);
        }

        PIORequest newRequest() {
            PIORequest req = new IORequest;
            memset(req, 0, sizeof(IORequest));
            return req;

        }

		void deleteList(PLIST_ENTRY list) {
            PIORequest p = static_cast<PIORequest>(ListFirst(&IODoneList));
            while (p != &IODoneList) {
                PIORequest q = static_cast<PIORequest>(p->Flink);
                delete p;
                p = q;
            }
        }

        //////////////////////////////////////////////////////////////////////
        // Callback based read/write

        aRequest areadAt(Handle hFile, FilePtr Offset, LPVOID lpBuffer, DWORD dwTransferCount,
                        aCallback *Callback, aHint *Hint, bool bRead = true)
        {
            qLock.lock();

            PIORequest IORequestPacket;
            if (IsListEmpty(&IODoneList))
                if (maxRequests) {
                    flush();
                    if (IsListEmpty(&IODoneList))
                        IORequestPacket = newRequest();
                    else {
                        IORequestPacket = static_cast<PIORequest>(ListFirst(&IODoneList));
                        RemoveEntryList(IORequestPacket);
                    }
                } else
                    IORequestPacket = newRequest();
            else {
                IORequestPacket = static_cast<PIORequest>(ListFirst(&IODoneList));
                RemoveEntryList(IORequestPacket);
            }

            IORequestPacket->Offset = (DWORD)(Offset & MAXDWORD);
            IORequestPacket->OffsetHigh = (DWORD)(Offset >> 32);
            IORequestPacket->hEvent = INVALID_HANDLE_VALUE;

            IORequestPacket->hFile = hFile;
            IORequestPacket->lpBuffer = lpBuffer;
            IORequestPacket->dwTransferCount = dwTransferCount;
            IORequestPacket->Callback = Callback;
            IORequestPacket->Hint = Hint;
            IORequestPacket->bRead = bRead;
            IORequestPacket->me = this;
            
            InsertTailList(&IORequestList, IORequestPacket);

            qEmpty.reset();
            qSemaphore.unlock();
            qLock.unlock();

            return IORequestPacket;
        }

        aRequest awriteAt(LPCVOID lpBuffer, HANDLE hFile, FilePtr Offset, DWORD dwTransferCount,
                         aCallback* Callback, aHint* Hint)
        {
            return areadAt(
                hFile, Offset, const_cast<LPVOID>(lpBuffer), dwTransferCount, 
                Callback, Hint, false);
        }


        //////////////////////////////////////////////////////////////////////
        // Event based read/write

        aRequest areadAt(Handle hFile, FilePtr Offset, LPVOID lpBuffer, DWORD dwTransferCount,
                        Event &event, bool bRead = true)
        {
            event.reset();
            qLock.lock();

            PIORequest IORequestPacket;
            if (IsListEmpty(&IODoneList))
                if (maxRequests) {
                    flush();
                    if (IsListEmpty(&IODoneList))
                        IORequestPacket = newRequest();
                    else {
                        IORequestPacket = static_cast<PIORequest>(ListFirst(&IODoneList));
                        RemoveEntryList(IORequestPacket);
                    }
                } else
                    IORequestPacket = newRequest();
            else {
                IORequestPacket = static_cast<PIORequest>(ListFirst(&IODoneList));
                RemoveEntryList(IORequestPacket);
            }

            IORequestPacket->Offset = (DWORD)(Offset & MAXDWORD);
            IORequestPacket->OffsetHigh = (DWORD)(Offset >> 32);
            IORequestPacket->hEvent = event.hEvent;

            IORequestPacket->hFile = hFile;
            IORequestPacket->lpBuffer = lpBuffer;
            IORequestPacket->dwTransferCount = dwTransferCount;
            IORequestPacket->bRead = bRead;
            IORequestPacket->me = this;

            InsertTailList(&IORequestList, IORequestPacket);

//            IORequestPacket->print();

            qEmpty.reset();
            qSemaphore.unlock();
            qLock.unlock();

            return IORequestPacket;
        }

        aRequest awriteAt(LPVOID lpBuffer, HANDLE hFile, FilePtr Offset, DWORD dwTransferCount,
                         Event &event)
        {
            return areadAt(
                hFile, Offset, lpBuffer, dwTransferCount, 
                event, false);
        }

        // you have to call this method after every event based asynchronous operation
        // to set request state from PENDING to DONE
        void release(PIORequest IORequestPacket) {
            qLock.lock();
            RemoveEntryList(IORequestPacket);
            InsertTailList(&IODoneList, IORequestPacket);
            qLock.unlock();
        }

        void flush() {
            qEmpty.wait();
        }

        void printList(PLIST_ENTRY list) {
			qLock.lock();
			for(PLIST_ENTRY	p = ListFirst(list); p != list; p = p->Flink) {
				PIORequest req = static_cast<PIORequest>(p);
                req->print();
				LPOVERLAPPED o = static_cast<LPOVERLAPPED>(req);
                DWORD trans = -1;
                GetOverlappedResult(req->hFile, req, &trans, false);

				if (HasOverlappedIoCompleted(o))
					printf("completed (%d)\n", trans);
                else
					printf("incomplete (%d)\n", trans);
			}
			qLock.unlock();
		}
        
        void status() {
            printf("REQUESTED\n");
            printList(&IORequestList);
            printf("\nPENDING\n");
            printList(&IOPendingList);
            printf("\nDONE\n");
            printList(&IODoneList);
        }


    private:

        static void CALLBACK IOCompletionRoutine(
            DWORD dwErrorCode,
            DWORD dwNumberOfBytesTransfered,
            LPOVERLAPPED lpOverlapped) 
        { 
            // If an I/O error occurs, display the error and exit. 
            if (dwErrorCode)
            {
                printf("FATAL I/O Callback Error %ld I/O Context %lx.%lx\n",
                    dwErrorCode,
                    reinterpret_cast<ULONGLONG>(lpOverlapped),
                    reinterpret_cast<ULONGLONG>(lpOverlapped->hEvent));
                return;
            }
            
            PIORequest IORequestPacket = static_cast<PIORequest>(lpOverlapped);
            IOQueue &me = *reinterpret_cast<IOQueue*>(IORequestPacket->me);

            aCallback *callback = IORequestPacket->Callback;
            aHint     *hint     = IORequestPacket->Hint;

            // Give the caller a call ...

            if (callback) callback(hint);

            // Now we got what we need and can release the request.

            me.release(IORequestPacket);

        }

        static DWORD WINAPI IOWorker(LPVOID _me)
        {
            IOQueue &me = *reinterpret_cast<IOQueue*>(_me);
            HANDLE HandleVector[2] = {me.qLock.hMutex, me.qSemaphore.hSemaphore};
            DWORD dummy;
         
            for(;;) 
            { 
         
                // Do an alertable wait on the handle vector. Both objects 
                // being signaled at the same time means that there is an 
                // I/O request in the queue and the caller has exclusive 
                // access to the queue. 
         
                DWORD CompletionStatus = WaitForMultipleObjectsEx(2, HandleVector, 
                        TRUE, INFINITE, TRUE); 
         
                // If the wait failed, error out.
         
                if (CompletionStatus == -1) 
                { 
                    printf("FATAL WAIT ERROR %ld\n", GetLastError()); 
                    return 1;
                } 
                // If an I/O completion occurred, wait for another 
                // I/O request or I/O completion. 
         
                if (CompletionStatus != WAIT_IO_COMPLETION) 
                { 
         
                    // The wait was satisfied. Ownership of the I/O 
                    // request queue is exclusive, and there is something in 
                    // the queue. To insert something in the queue, the 
                    // inserter gets the list lock (mutex), inserts an entry, 
                    // signals the list semaphore, and finally releases the 
                    // list lock. 
         
                    PIORequest IORequestPacket =
                        static_cast<PIORequest>(ListFirst(&me.IORequestList));
                    RemoveEntryList(IORequestPacket);
         
                    BOOL IOOperationStatus;

                    if (IORequestPacket->hEvent == INVALID_HANDLE_VALUE) {
                        // Callback based

                        if (IORequestPacket->bRead)
                            IOOperationStatus =
                                ReadFileEx(
                                    IORequestPacket->hFile,
                                    IORequestPacket->lpBuffer,
                                    IORequestPacket->dwTransferCount,
                                    IORequestPacket,
                                    IOCompletionRoutine);
                        else
                            IOOperationStatus =
                                WriteFileEx(
                                    IORequestPacket->hFile,
                                    IORequestPacket->lpBuffer,
                                    IORequestPacket->dwTransferCount,
                                    IORequestPacket,
                                    IOCompletionRoutine);
					} else {
                        // Event based

                        if (IORequestPacket->bRead) 
                            IOOperationStatus = 
                                ReadFile(
                                    IORequestPacket->hFile, 
                                    IORequestPacket->lpBuffer, 
                                    IORequestPacket->dwTransferCount,
                                    &dummy,
                                    IORequestPacket); 
                        else
                            IOOperationStatus = 
                                WriteFile(
                                    IORequestPacket->hFile, 
                                    IORequestPacket->lpBuffer, 
                                    IORequestPacket->dwTransferCount,
                                    &dummy,
                                    IORequestPacket); 
                    }
         
	                // Test to see if the I/O was queued successfully.

					if (IOOperationStatus || (GetLastError() == ERROR_IO_PENDING)) {
						InsertTailList(&me.IOPendingList, IORequestPacket);
						me.qLock.unlock();
					} else {
						printf("FATAL I/O Queuing Error %ld\n", GetLastError());
                        IORequestPacket->print();
						InsertTailList(&me.IODoneList, IORequestPacket);
						me.qLock.unlock();
                        me.status();
//						return 1;
					}

                    // The I/O queued successfully. Go back into the 
                    // alertable wait for I/O completion or for 
                    // more I/O requests. 
         
                } else {
                    me.qLock.lock();
                    if (IsListEmpty(&me.IORequestList) && IsListEmpty(&me.IOPendingList))
                        me.qEmpty.signal();
                    me.qLock.unlock();
                }

            } 

        }

    };

}

#endif
