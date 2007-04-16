#ifndef SEQAN_HEADER_FILE_H
#define SEQAN_HEADER_FILE_H

//____________________________________________________________________________
// prerequisites

#include <seqan/sequence.h>

#include <iostream>
#include <climits>
#include <cstdio>
#include <list>
#include <vector>

//____________________________________________________________________________
// file formats

#include <seqan/file/file_forwards.h>

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/file/file_generated_forwards.h>
#endif

#include <seqan/file/chunk_collector.h>
#include <seqan/file/meta.h>

#include <seqan/file/file_format.h>

#include <seqan/file/stream.h>
#include <seqan/file/cstream.h>
#include <seqan/file/stream_algorithms.h>

#include <seqan/file/file_format_raw.h>
#include <seqan/file/file_format_fasta.h>
#include <seqan/file/file_format_fasta_align.h>
#include <seqan/file/file_format_cgviz.h>

#include <seqan/file/file_format_guess.h>

//____________________________________________________________________________
// files

#include <seqan/file/file_base.h>
#include <seqan/file/file_cstyle.h>
#include <seqan/file/file_array.h>

#include <seqan/system.h>	// (include will be removed soon)

#include <seqan/pipe/pipe_base.h>
#include <seqan/pipe/pipe_source.h>

//____________________________________________________________________________
// external string

#include <seqan/file/file_page.h>
#include <seqan/file/file_page_raid0.h>
#include <seqan/file/string_external.h>

#endif //#ifndef SEQAN_HEADER_...
