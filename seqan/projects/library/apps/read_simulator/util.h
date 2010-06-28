#ifndef UTIL_H_
#define UTIL_H_

using namespace seqan;

// Helper function to write out a {Dna,Dna5}Q sequence with qualities to a
// FASTQ file.
template <typename TStream, typename TIdStringSet, typename TSeqStringSet>
void write(TStream & stream,
           TIdStringSet & seqIds,
           TSeqStringSet & sequences,
           Fastq const &) {
    typedef TSeqStringSet TStringSet;
    typedef typename Position<TStringSet>::Type TPosition;

    CharString qualBuffer;
    for (TPosition i = 0; i < length(sequences); ++i) {
        stream << "@" << seqIds[i] << std::endl;
        stream << sequences[i] << std::endl;
        stream << "+" << /*seqIds[i] << */std::endl;
        resize(qualBuffer, length(sequences[i]), Exact());
        for (TPosition j = 0; j < length(sequences[i]); ++j)
            qualBuffer[j] = getQualityValue(sequences[i][j]) + '!';
        stream << qualBuffer << std::endl;
    }
}

#endif  // UTIL_H_
