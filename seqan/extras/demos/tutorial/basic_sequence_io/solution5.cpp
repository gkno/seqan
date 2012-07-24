#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main(int argc, char const ** argv)
{
    if (argc < 2)
    {
        std::cerr << "USAGE: basic_seq_io_example FILENAME\n";
        return 1;
    }

    seqan::SequenceStream seqStream(argv[1], seqan::SequenceStream::WRITE);
    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }

    seqan::CharString id = "seq1";
    seqan::Dna5String seq = "CGAT";

    if (writeRecord(seqStream, id, seq) != 0)
    {
        std::cerr << "ERROR: Could not write to file!\n";
        return 1;
    }

    id = "seq1";
    seq = "TTTT";

    if (writeRecord(seqStream, id, seq) != 0)
    {
        std::cerr << "ERROR: Could not write to file!\n";
        return 1;
    }

    return 0;
}
