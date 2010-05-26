/*
  Usage: seq_sample N IN.fastq OUT.fastq

  Samples N reads from IN.fastq and writes them to OUT.fastq without
  duplicates.
 */

#include <fstream>
#include <iostream>
#include <set>

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/misc/misc_random.h>
#include <seqan/sequence.h>
#include <seqan/store.h>

using namespace seqan;

template < 
	typename TStream, 
	typename TSeq >
void dumpWrapped(
	TStream &out,
	TSeq &seq)
{
	unsigned size = length(seq);
	unsigned i;
	for (i = 0; i + 60 < size; i += 60)
		out << infix(seq, i, i + 60) << std::endl;
	out << infix(seq, i, size) << std::endl;
}


template <typename TStream, typename TStringSetSpec, typename TStringSetSpec2>
void write(TStream & stream,
           StringSet<String<Dna5Q>, TStringSetSpec> const & sequences,
           StringSet<CharString, TStringSetSpec2> const & seqIds,
           Fastq const &) {
    typedef StringSet<String<Dna5Q>, TStringSetSpec> TStringSet;
    typedef typename Position<TStringSet>::Type TPosition;

    CharString qualBuffer;
    for (TPosition i = 0; i < length(sequences); ++i) {
        stream << "@" << seqIds[i] << std::endl;
        dumpWrapped(stream, sequences[i]);
        stream << "+" << seqIds[i] << std::endl;
        resize(qualBuffer, length(sequences[i]), Exact());
        for (TPosition j = 0; j < length(sequences[i]); ++j) {
            qualBuffer[j] = getQualityValue(sequences[i][j]) + '!';
        }
        dumpWrapped(stream, qualBuffer);
    }
}


int main(int argc, char **argv) {
    if (argc != 4) {
        std::cerr << "ERROR: Wrong number of parameters!" << std::endl;
        std::cerr << "USAGE: seq_sample N IN.fastq OUT.fastq" << std::endl;
        return 1;
    }

    // TODO(holtgrew): Actually, we want to seed this!
//     const unsigned SEED = 42;
    mtRandInit();

    unsigned n = atoi(argv[1]);

    // Load reads using the fragment store.
    FragmentStore<> fragStore;
    if (!loadReads(fragStore, argv[2])) return 1;

    // The factor 10 is arbitrarily chosen to make sure the generation
    // is sufficiently fast.
    if (length(fragStore.readSeqStore) < 10 * n) {
        std::cerr << "Number of reads in file should be >= 10*n" << std::endl;
        return 1;
    }

    // Randomly build the resulting string set.
    // Pick ids to select.
    std::set<unsigned> selectedIds;
    while (length(selectedIds) < n) {
        unsigned x = mtRand() % length(fragStore.readSeqStore);
        selectedIds.insert(x);
    }
    // Actually build result.
    StringSet<String<Dna5Q> > resultSeqs;
    reserve(resultSeqs, n);
    StringSet<CharString> resultIds;
    reserve(resultIds, n);
    for (std::set<unsigned>::const_iterator it = selectedIds.begin(); it != selectedIds.end(); ++it) {
        appendValue(resultSeqs, fragStore.readSeqStore[*it]);
        appendValue(resultIds, fragStore.readNameStore[*it]);
    }

    // Write out the result.
    if (CharString("-") == argv[3]) {
        write(std::cout, resultSeqs, resultIds, Fastq());
    } else {
        std::fstream fstrm(argv[3], std::ios_base::out);
        if (not fstrm.is_open()) {
            std::cerr << "Could not open out file \"" << argv[3] << "\"" << std::endl;
            return 1;
        }
        write(fstrm, resultSeqs, resultIds, Fastq());
    }

    return 0;
}
