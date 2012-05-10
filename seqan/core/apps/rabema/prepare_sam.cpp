// ---------------------------------------------------------------------------
// Tool to prepare SAM file for Rabema.  Currently, this only consists of
// replacing all sequence and quality fields that store a "*" with the value
// from the primary alignment.
// ---------------------------------------------------------------------------

#include <iostream>

#include <seqan/bam_io.h>
#include <seqan/stream.h>

int fixRecords(seqan::String<seqan::BamAlignmentRecord> & records)
{
    using namespace seqan;

    if (empty(records))
        return 0;  // OK to be empty.

    // Pick indices with sequences for first/second sequence.
    int idxSeqFirst = -1, idxSeqSecond = -1;
    for (unsigned i = 0; i < length(records); ++i)
    {
        if (hasFlagFirst(records[i]) && !empty(records[i].seq))
            idxSeqFirst = i;
        else if (hasFlagLast(records[i]) && !empty(records[i].seq))
            idxSeqSecond = i;
        else if (!hasFlagFirst(records[i]) && !hasFlagLast(records[i]) && !empty(records[i].seq))
            idxSeqFirst = i;
    }

    // Get sequences for first and second.
    Dna5String seqFirst;
    CharString qualFirst;
    if (idxSeqFirst != -1)
    {
        seqFirst = records[idxSeqFirst].seq;
        if (hasFlagRC(records[idxSeqFirst]))
            reverseComplement(seqFirst);
        qualFirst = records[idxSeqFirst].qual;
        if (hasFlagRC(records[idxSeqFirst]))
            reverse(qualFirst);
    }
    Dna5String seqFirstRC = seqFirst;
    reverseComplement(seqFirstRC);
    CharString qualFirstRC = qualFirst;
    reverse(qualFirst);
    Dna5String seqSecond;
    CharString qualSecond;
    if (idxSeqSecond != -1)
    {
        seqSecond = records[idxSeqSecond].seq;
        if (hasFlagRC(records[idxSeqSecond]))
            reverseComplement(seqSecond);
        qualSecond = records[idxSeqSecond].qual;
        if (hasFlagRC(records[idxSeqSecond]))
            reverse(qualSecond);
    }
    Dna5String seqSecondRC = seqSecond;
    reverseComplement(seqSecondRC);
    CharString qualSecondRC = qualSecond;
    reverse(qualSecond);

    // Actually assign sequences into records.
    for (unsigned i = 0; i < length(records); ++i)
    {
        if (hasFlagFirst(records[i]) || (!hasFlagFirst(records[i]) && !hasFlagLast(records[i])))
        {
            if (idxSeqFirst == -1)
            {
                std::cerr << "ERROR: No sequence for first mate of query name " << records[i].qName << ".\n";
                return 1;
            }
            if (hasFlagRC(records[i]))
            {
                if (!empty(records[i].seq) && records[i].seq != seqFirstRC)
                    SEQAN_FAIL("ERROR: Mismatching sequences for query name %s.", toCString(records[i].qName));
                records[i].seq = seqFirstRC;
                records[i].qual = qualFirstRC;
            }
            else
            {
                if (!empty(records[i].seq) && records[i].seq != seqFirst)
                    SEQAN_FAIL("ERROR: Mismatching sequences for query name %s.", toCString(records[i].qName));
                records[i].seq = seqFirst;
                records[i].qual = qualFirst;
            }
        }
        else if (hasFlagLast(records[i]))
        {
            if (idxSeqSecond == -1)
            {
                std::cerr << "ERROR: No sequence for second mate of query name " << records[i].qName << ".\n";
                return 1;
            }
            if (hasFlagRC(records[i]))
            {
                if (!empty(records[i].seq) && records[i].seq != seqSecondRC)
                    SEQAN_FAIL("ERROR: Mismatching sequences for query name %s.", toCString(records[i].qName));
                records[i].seq = seqSecondRC;
                records[i].qual = qualSecondRC;
            }
            else
            {
                if (!empty(records[i].seq) && records[i].seq != seqSecond)
                    SEQAN_FAIL("ERROR: Mismatching sequences for query name %s.", toCString(records[i].qName));
                records[i].seq = seqSecond;
                records[i].qual = qualSecond;
            }
        }
    }
    
    return 0;
}

int main(int argc, char const ** argv)
{
    using namespace seqan;

    // Check arguments.
    if (argc != 2)
    {
        std::cerr << "USAGE: parepare_sam IN.sam > FIXED.sam\n";
        return 1;
    }

    // Open SAM file for reading.
    std::ifstream inSam(argv[1], std::ios::in | std::ios::binary);
    if (!inSam.good())
    {
        std::cerr << "Could not open file " << argv[1] << " for reading.\n";
        return 1;
    }
    RecordReader<std::ifstream, SinglePass<> > reader(inSam);

    // Read header.
    typedef StringSet<CharString>      TNameStore;
    typedef NameStoreCache<TNameStore> TNameStoreCache;
    
    TNameStore refNameStore;
    TNameStoreCache refNameStoreCache(refNameStore);
    BamIOContext<TNameStore> context(refNameStore, refNameStoreCache);

    BamHeader header;
    if (readRecord(header, context, reader, Sam()) != 0)
    {
        std::cerr << "ERROR: Could not read header from SAM file.\n";
        return 1;
    }

    // Check sort order.
    if (getSortOrder(header) != BAM_SORT_QUERYNAME)
    {
        std::cerr << "SAM file not sorted by 'queryname'\n";
        return 1;
    }

    write2(std::cout, header, context, Sam());

    // Read file in chunks, one for each query name.
    String<BamAlignmentRecord> records;
    BamAlignmentRecord record;
    while (!atEnd(reader))
    {
        if (readRecord(record, context, reader, Sam()) != 0)
        {
            std::cerr << "ERROR reading SAM record!\n";
            return 0;
        }

        if (!empty(records) && record.qName != back(records).qName)
        {
            if (fixRecords(records) != 0)
            {
                std::cerr << "Could not fix records!\n";
                return 1;
            }
            for (unsigned i = 0; i < length(records); ++i)
                write2(std::cout, records[i], context, Sam());
            clear(records);
        }

        appendValue(records, record);
    }
    if (fixRecords(records) != 0)
    {
        std::cerr << "Could not fix records!\n";
        return 1;
    }
    for (unsigned i = 0; i < length(records); ++i)
        write2(std::cout, records[i], context, Sam());
    
    return 0;
}