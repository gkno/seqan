#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/store.h>

using namespace seqan;


template <typename TFragmentStore>
void performEvaluation(TFragmentStore & store) {
    typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignedReadsIter;
    typedef typename TFragmentStore::TContigStore TContigStore;
    typedef typename TFragmentStore::TReadSeqStore TReadSeqStore;
    typedef typename TFragmentStore::TReadSeq TReadSeq;
    typedef typename Value<TContigStore>::Type TContigStoreElement;
    typedef typename Value<TAlignedReadStore>::Type TAlignedReadStoreElement;
    typedef typename TAlignedReadStoreElement::TGapAnchors TReadGapAnchors;
    typedef typename TContigStoreElement::TContigSeq TContigSeq;
    typedef typename TContigStoreElement::TGapAnchors TContigGapAnchors;
    typedef Gaps<TReadSeq, AnchorGaps<TContigGapAnchors> > TReadGaps;
    typedef Gaps<TContigSeq, AnchorGaps<TContigGapAnchors> > TContigGaps;
    typedef typename Iterator<TContigGaps, Standard>::Type TContigGapAnchorsIterator;
    typedef typename Iterator<TReadGaps, Standard>::Type TReadGapAnchorsIterator;

    size_t alignedReadId = 0;
    for (TAlignedReadsIter it = begin(store.alignedReadStore, Standard()); it != end(store.alignedReadStore, Standard()); ++it, ++alignedReadId) {
        // Get contig and read sequences.
        TContigSeq const & contigSeq = store.contigStore[it->contigId].seq;
        TReadSeq readSeq = store.readSeqStore[it->readId];
        // Get gaps for contig and read.
        TContigGaps contigGaps(contigSeq, store.contigStore[it->contigId].gaps);
        TReadGaps readGaps(readSeq, store.alignedReadStore[alignedReadId].gaps);
        // Limit contig gaps to aligned read position.
        setBeginPosition(contigGaps, _min(it->beginPos, it->endPos));
        setEndPosition(contigGaps, _max(it->beginPos, it->endPos));

        // Can reverse complement readSeq in place since gaps store a holder.
        if (it->beginPos > it->endPos) {
//             std::cerr << "complementing..." << std::endl;
            reverseComplementInPlace(readSeq);
        }

//         std::cerr << "Read:      " << readSeq << std::endl;
//         std::cerr << "Qualities: ";
//         for (unsigned i = 0; i < length(readSeq); ++i)
//             std::cerr << getQualityValue(readSeq[i]) << ", ";
//         std::cerr << std::endl;

        // TODO(holtgrew): Switch read and contig gaps it.
        int alignQualScore = 0;
        TContigGapAnchorsIterator contigGapsIt = begin(contigGaps, Standard());
        for (TReadGapAnchorsIterator readGapsIt = begin(readGaps, Standard()); readGapsIt != end(readGaps, Standard()); ++contigGapsIt, ++readGapsIt) {
//             if (isGap(readGapsIt) || isGap(contigGapsIt))
//                 continue;
//             if (convert<Dna5>(*readGapsIt) == convert<Dna5>(*contigGapsIt))
//                 continue;
//             std::cerr << "mismatch: " << convert<Dna5>(*readGapsIt) << " != " << convert<Dna5>(*contigGapsIt) << " score is " << getQualityValue(convert<Dna5Q>(*readGapsIt)) << std::endl;
            alignQualScore += getQualityValue(convert<Dna5Q>(*readGapsIt));

            Dna5Q x = *readGapsIt;
            std::cout << x << std::endl;
        }

        printf("%6lu\t", alignedReadId);
        std::cout << store.readNameStore[it->readId];
        printf("\t%3d\n", alignQualScore);
    }
}


int main(int argc, char **argv) {
    // Check arguments.
    if (argc != 3) {
        std::cerr << "Usage:  sam_info GENOME.FASTA FILE.SAM" << std::endl;
        return 1;
    }

    // Load files.
    typedef FragmentStore<> TFragmentStore;
    TFragmentStore fragments;

    // Load Contigs.
    std::cerr << "Reading FASTA contigs sequence file " << argv[1] << " ..." << std::endl;
    if (!loadContigs(fragments, argv[1])) {
        std::cerr << "Could not read contigs." << std::endl;
        return 1;
    }

    // Load SAM File.
    std::cerr << "Reading SAM file file " << argv[2] << " ..." << std::endl;
    {
        std::fstream fstrm(argv[2], std::ios_base::in | std::ios_base::binary);
        if (not fstrm.is_open()) {
            std::cerr << "Could not open SAM file." << std::endl;
            return 1;
        }
        read(fstrm, fragments, SAM());
    }

    std::cerr << "Evaluating..." << std::endl;
    performEvaluation(fragments);
    
    return 0;
}
