// FRAGMENT(header)
#include <fstream>
#include <iostream>

#include <seqan/basic.h>
#include <seqan/file.h>

#include <seqan/stream.h>

using namespace seqan;

int main(int argc, char const ** argv)
{
    if (argc != 2)
    {
        std::cerr << "ERROR: Invalid argument count." << std::endl
                  << "USAGE: " << argv[0] << " FILE" << std::endl;
        return 1;
    }
    
// FRAGMENT(read-sequences-single-pass)
    // Open file and create RecordReader.
    std::ifstream fasta(argv[1], std::ios_base::in | std::ios_base::binary);
    if (!fasta.good())
        return 1;
    RecordReader<std::ifstream, SinglePass<> > reader(fasta);
    
    // Define variables for storing the sequences and sequence ids.
    StringSet<CharString> ids;
    StringSet<String<Dna5Q> > seqs;
    if (read2(ids, seqs, reader, Fastq()) != 0)
    {
        std::cerr << "ERROR reading FASTA." << std::endl;
        return 1;
    }
    
    // Write ids, sequences and qualities to output.
    typedef Iterator<StringSet<CharString>, Rooted>::Type TIdIter;
    typedef Iterator<StringSet<String<Dna5Q> >, Standard>::Type TSeqIter;
    TIdIter idIt = begin(ids, Rooted());
    TSeqIter seqIt = begin(seqs, Standard());
    CharString quals;
    for (; !atEnd(idIt); goNext(idIt), ++seqIt)
    {
        resize(quals, length(*seqIt));
        assignQualities(quals, *seqIt);
        std::cout << *idIt << '\t' << *seqIt << '\t' << quals << std::endl;
    }
    

// FRAGMENT(read-sequences-double-pass)

// FRAGMENT(footer)    
    return 0;
}
