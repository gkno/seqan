cmake_minimum_required (VERSION 2.6)
project (sandbox_emde_apps_razers3)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(razers3 razers.cpp
			     razers.h
			     job_queue.h
			     outputFormat.h
			     parallel_atomic.h
			     parallel_job_queue.h
			     parallel_misc.h
			     parallel_store.h
			     paramChooser.h
			     profile_timeline.h
			     razers_match_filter.h
			     razers_matepairs.h
			     razers_matepairs_parallel.h
			     razers_parallel.h
			     razers_parallel_reads.h
			     razers_window.h
			     readSimulator.h)