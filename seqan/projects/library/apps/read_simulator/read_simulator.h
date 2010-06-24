/* Globally shared code. */

#ifndef READ_SIMULATOR_H_
#define READ_SIMULATOR_H_

using namespace seqan;

template <
	typename TStream, 
	typename TSeq >
void _dumpWrapped(
	TStream & out,
	TSeq const & seq)
{
	unsigned size = length(seq);
	unsigned i;
	for (i = 0; i + 60 < size; i += 60)
		out << infix(seq, i, i + 60) << std::endl;
	out << infix(seq, i, size) << std::endl;
}


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


enum ErrorType {
    ERROR_TYPE_MATCH    = 0,
    ERROR_TYPE_MISMATCH = 1,
    ERROR_TYPE_INSERT   = 2,
    ERROR_TYPE_DELETE   = 3
};


struct GlobalOptions {
    enum ReadsType {
        READS_TYPE_ILLUMINA
    };

    ReadsType readsType;
};


// Taken from akemde's read simulator.
#ifdef USE_LOGVALUES

	template <typename TValue>
	inline long double
	_transform(TValue a)
	{
		return log(a);
	}

	template <typename TValue>
	inline long double
	_transformBack(TValue a)
	{
		return exp(a);
	}

	//////////////////////////////////////////////////////////////////////////////
	// Returns the sum of two probability values in log space
	template <typename TValue>
	inline void
	_probAdd(TValue &a, TValue b)
	{
		if (isinf(a)) {
			a = b;
			return;
		}
		if (isinf(b)) return;
		if (isnan(a + log(1 + exp(b - a)))) return;
		a += log(1 + exp(b - a));
	}

	template <typename TValue>
	inline TValue
	_probMul(TValue a, TValue b)
	{
		return a + b;
	}

	template <typename TValue>
	inline TValue
	_probDiv(TValue a, TValue b)
	{
		return a - b;
	}

#else  // USE_LOGVALUES

	template <typename TValue>
	inline TValue
	_transform(TValue a)
	{
		return a;
	}

	template <typename TValue>
	inline TValue
	_transformBack(TValue a)
	{
		return a;
	}

	template <typename TValue>
	inline void
	_probAdd(TValue &a, TValue b)
	{
		a += b;
	}

	template <typename TValue>
	inline TValue
	_probMul(TValue a, TValue b)
	{
		return a * b;
	}

	template <typename TValue>
	inline TValue
	_probDiv(TValue a, TValue b)
	{
		return a / b;
	}

#endif  // USE_LOGVALUES

// Write a random DNA sequence of the given length to the file with the given name.
int writeRandomSequence(size_t length, CharString const & fileName) {
    DnaString randomSequence;
    reserve(randomSequence, length);

    for (size_t i = 0; i < length; ++i) {
        Dna c = static_cast<Dna>(mtRand() % ValueSize<Dna>::VALUE);
        appendValue(randomSequence, c);
    }

    std::ofstream file;
    file.open(toCString(fileName), std::ios_base::out | std::ios_base::trunc);
    if (!file.is_open()) {
        std::cerr << "Failed to write random sequence to " << fileName << std::endl;
        return 1;
    }
    write(file, randomSequence, "random_sequence", Fasta());
    file.close();
    return 0;
}

namespace seqan {
struct MyFragmentStoreConfig {};

template<>
struct FragmentStoreConfig<MyFragmentStoreConfig>
{
	typedef String<Dna5Q>	TReadSeq;
	typedef String<Dna5Q>	TContigSeq;
	
	typedef double			TMean;
	typedef double			TStd;
	typedef signed char		TMappingQuality;
		
	typedef void    TReadStoreElementSpec;
	typedef Owner<Default> TReadSeqStoreSpec;
	typedef void    TMatePairStoreElementSpec;
	typedef void    TLibraryStoreElementSpec;
	typedef void    TContigStoreElementSpec;
	typedef void    TContigFileSpec;
	typedef void    TAlignedReadStoreElementSpec;
	typedef Owner<Default>	TAlignedReadTagStoreSpec;
	typedef void    TAnnotationStoreElementSpec;
};
}  // namespace seqan

#endif  // READ_SIMULATOR_H_
