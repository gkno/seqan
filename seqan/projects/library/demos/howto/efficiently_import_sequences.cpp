// FRAGMENT(includes)
#include <seqan/file.h>
#include <iostream>

using namespace seqan;

// FRAGMENT(open_file)
int main (int argc, char const * argv[])
{
	MultiSeqFile multiSeqFile;
	if (argc < 2 || !open(multiSeqFile.concat, argv[1], OPEN_RDONLY))
		return 1;

// FRAGMENT(guess_and_split)
	AutoSeqFormat format;
	guessFormat(multiSeqFile.concat, format);
	split(multiSeqFile, format);

// FRAGMENT(read_sequences)
	unsigned seqCount = length(multiSeqFile);
	StringSet<String<Dna5Q> > seqs;
	StringSet<CharString> seqIDs;

	reserve(seqs, seqCount, Exact());
	reserve(seqIDs, seqCount, Exact());

	String<Dna5Q> seq;
	CharString qual;
	CharString id;

	for(unsigned i = 0; i < seqCount; ++i)
	{
		assignSeq(seq, multiSeqFile[i], format);    // read sequence
		assignQual(qual, multiSeqFile[i], format);  // read ascii quality values
		assignSeqId(id, multiSeqFile[i], format);   // read sequence id

		// convert ascii to values from 0..62
		// store dna and quality together in Dna5Q
		for (unsigned j = 0; j < length(qual) && j < length(seq); ++j)
			assignQualityValue(seq[j], (int)(ordValue(qual[j]) - 33));

		// we use reserve and append, as assign is not supported
		// by StringSet<..., Owner<ConcatDirect<> > >
		appendValue(seqs, seq, Generous());
		appendValue(ids, id, Generous());
	}

// FRAGMENT(output)
	for (unsigned i = 0; i < seqCount && i < 10; ++i)
	{
		std::cout << '>' << ids[i] << std::endl;
		std::cout << seqs[i] << std::endl;
	}

	return 0;
}

