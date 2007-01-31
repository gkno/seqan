#ifndef SEQAN_HEADER_FILE_H
#define SEQAN_HEADER_FILE_H

// should be set before anything else
#define _FILE_OFFSET_BITS 64

#include <vector>

#include <seqan/sequence.h>

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


#include <seqan/file/file_base.h>
#include <seqan/basic.h>
#include <seqan/system.h>
#include <seqan/file/file_cstyle.h>
#include <seqan/file/file_sync.h>
#include <seqan/file/file_async.h>
#include <seqan/file/file_array.h>

#include <seqan/pipe/pipe_base.h>
#include <seqan/pipe/pipe_source.h>
#include <seqan/sequence/string_external.h>

#endif //#ifndef SEQAN_HEADER_...
