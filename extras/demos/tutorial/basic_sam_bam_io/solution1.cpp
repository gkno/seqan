#include <iostream>
#include <seqan/bam_io.h>

int main()
{
    // Open input stream, BamStream can read SAM and BAM files.
    seqan::BamStream bamStreamIn("example.sam");
    if (!isGood(bamStreamIn))
    {
        std::cerr << "ERROR: Could not open example.sam!\n";
        return 1;
    }
    // Open output stream, "-" means stdin on if reading, else stdout.
    seqan::BamStream bamStreamOut("-", seqan::BamStream::WRITE);
    // Copy header.  The header is automatically written out before
    // the first record.
    bamStreamOut.header = bamStreamIn.header;

    seqan::BamAlignmentRecord record;
    while (!atEnd(bamStreamIn))
    {
        if (readRecord(record, bamStreamIn) != 0)
        {
            std::cerr << "ERROR: Could not read record!\n";
            return 1;
        }
        if (writeRecord(bamStreamOut, record) != 0)
        {
            std::cerr << "ERROR: Could not write record!\n";
            return 1;
        }
    }
    
    return 0;
}
