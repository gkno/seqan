#include <iostream>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <time.h>
#define SEQAN_DEBUG
#define SEQAN_TEST

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>  // Header under test.

#include "test_embl.h"
#include "test_file.h"

SEQAN_BEGIN_TESTSUITE(test_file)
{
    SEQAN_CALL_TEST(test_file_stream);
    SEQAN_CALL_TEST(test_file_cstream);

    SEQAN_CALL_TEST(test_file_raw);

    SEQAN_CALL_TEST(test_file_fasta_crlf);
    SEQAN_CALL_TEST(test_file_fasta_lf);
    SEQAN_CALL_TEST(test_file_fasta_cr);
    SEQAN_CALL_TEST(test_file_fasta_write);

    SEQAN_CALL_TEST(test_file_cgviz);
    SEQAN_CALL_TEST(test_file_fasta_align);
    SEQAN_CALL_TEST(test_file_cgviz);

    SEQAN_CALL_TEST(test_file_embl);
    SEQAN_CALL_TEST(test_file_genbank);

//    SEQAN_CALL_TEST(test_file_reader_iterator);
    SEQAN_CALL_TEST(test_file_reader_string);
//    SEQAN_CALL_TEST(test_file_reader_string2_fasta);
//    SEQAN_CALL_TEST(test_file_reader_string2_embl);
//    SEQAN_CALL_TEST(test_file_reader_string2_genbank);
    SEQAN_CALL_TEST(test_file_reader_string3);

	SEQAN_CALL_TEST(test_file_embl_file);
	SEQAN_CALL_TEST(test_file_embl_meta);
    
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/file/stream.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/file/cstream.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/file/meta.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/file/file_base.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/file/file_format_raw.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/file/file_format_fasta.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/file/file_format_fasta_align.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/file/file_format_embl.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/file/file_format_cgviz.h");
    SEQAN_VERIFY_CHECKPOINTS("projects/library/seqan/file/file_format_guess.h");
}
SEQAN_END_TESTSUITE
