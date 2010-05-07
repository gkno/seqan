#ifndef WIT_BUILDER_WIT_BUILDER_OPTIONS_H_
#define WIT_BUILDER_WIT_BUILDER_OPTIONS_H_

// Container for the program options.
struct Options
{
    // Whether hits should be verified with MyersUkonnen search.
    bool verify;

    // Whether to reverse-complement the genome.  If false, the reads are
    // reverse complemented to find matches on the reverse strand.
//     bool reverseComplementGenome;

    // Whether N should match all characters without penalty.
    bool matchN;
    
    // Maximum errors in percent, relative to the read length.
    double maxError;

    // Distance function to use, also see validDistanceFunction.
    String<char> distanceFunction;

    // Name of reference contigs file name.
    String<char> referenceSeqFilename;

    // Name of SAM file with golden reads.
    String<char> perfectMapFilename;

    // Return true iff distanceFunction is a valid distance function.
    // Can be one of {"hamming", "edit"}.
    bool validDistanceFunction() const
    {
        if (distanceFunction == "hamming") return true;
        if (distanceFunction == "edit") return true;
        return false;
    }
};

#endif  // WIT_BUILDER_WIT_BUILDER_OPTIONS_H_
