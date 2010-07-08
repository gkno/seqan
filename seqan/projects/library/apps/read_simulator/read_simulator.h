/* Globally shared code. */

#ifndef READ_SIMULATOR_H_
#define READ_SIMULATOR_H_

#include <numeric>

#include <seqan/random.h>
#include <seqan/sequence_journaled.h>

#include "store_config.h"
#include "util.h"

using namespace seqan;

// ============================================================================
// Enums, Tags, Classes.
// ============================================================================

// Enum describing the read type to be simulated.
enum ReadsType {
    READS_TYPE_ILLUMINA,
    READS_TYPE_454
};

typedef void Global;

template <typename TTag>
struct Options;

// Program-wide options.
template <>
struct Options<Global>
{
    // true iff help is to be shown.
    bool showHelp;
    // true iff verbosity is enabled.
    bool verbose;
    // true if very verbose is enabled.
    bool veryVerbose;

    // The type of the reads to be simulated.
    ReadsType readsType;

    // Basic Read Simulation Parameters.

    // Seed to use for the random number generator.
    unsigned seed;
    // Number of reads (pairs) to simulate.
    unsigned numReads;
    // true iff a random reference sequence is to be used.
    bool useRandomSequence;
    // Length of random sequence to be simulated.
    unsigned randomSourceLength;
    // true iff only the forward strand is to be simulated.
    bool onlyForward;
    // true iff only the reverse strand is to be simulated.
    bool onlyReverse;
    // The output file.  Defaults to REFERENCE-FILE.reads.fastq,
    // possibly with a ".1" or ".2" before the ".fastq" if mate pairs
    // are simulated.
    CharString outputFile;
    // Path to the SAM file to generate.  Defaults to fastq file name
    // with suffix ".sam"
    CharString samFile;
    // true iff qualities are to be simulated.
    bool simulateQualities;

    // Mate-Pair Related Options.

    // true iff generating mate pairs is enabled.
    bool generateMatePairs;
    // true iff mate pair library sizes are to be uniformly distributed,
    // otherwise standard distribution is used.
    bool libraryLengthIsUniform;
    // Mate-pair library mean length.
    unsigned libraryLengthMean;
    // Mate-pair library length error.  Standard deviation for normally
    // distributed library lengths, interval length around mean for uniform
    // distribution.
    unsigned libraryLengthError;

    // Source Sequence Repeat Parameters.

    // TODO(holtgrew): Mostly interesting for randomly generated source sequences.

    // Haplotype parameters.

    // Number of haplotypes to generated.  All are generated from the input genome.
    unsigned numHaplotypes;
    // SNP rate.
    double haplotypeSnpRate;
    // Indel rate.
    double haplotypeIndelRate;
    // Smallest number of indels.
    unsigned haplotypeIndelRangeMin;
    // Largest number of indels.
    unsigned haplotypeIndelRangeMax;

    Options()
            : showHelp(false),
              verbose(false),
              veryVerbose(false),
              seed(0),
              numReads(1000),
              useRandomSequence(false),
              randomSourceLength(1000*1000),
              onlyForward(false),
              onlyReverse(false),
              outputFile(""),
              samFile(""),
              simulateQualities(false),
              generateMatePairs(false),
              libraryLengthMean(1000),
              libraryLengthError(100),
              numHaplotypes(1),
              haplotypeSnpRate(0.001),
              haplotypeIndelRate(0.001),
              haplotypeIndelRangeMin(1),
              haplotypeIndelRangeMax(6)
    {}
};

// Use this container for model specific parameters generated before the actual
// simulation.  Setup in void simulateReadsSetupModelSpecificData(...).
template <typename TTag>
struct ModelParameters;

// Enum describing the type of an error.
enum ErrorType {
    ERROR_TYPE_MATCH    = 0,
    ERROR_TYPE_MISMATCH = 1,
    ERROR_TYPE_INSERT   = 2,
    ERROR_TYPE_DELETE   = 3
};

template <typename TReadTypeTag>
struct ReadSimulationInstruction;

template <>
struct ReadSimulationInstruction<Global> {
    unsigned haplotype;
    unsigned contigId;
    bool isForward;
    size_t beginPos;
    size_t endPos;
    // Number of characters added/removed to the string by indels.
    unsigned delCount;
    unsigned insCount;
    String<ErrorType> editString;
    String<int> qualities;

    ReadSimulationInstruction() : delCount(0), insCount(0) {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TStream>
TStream & operator<<(TStream & stream, Options<Global> const & options) {
    stream << "global-options {" << std::endl
           << "  seed:                   " << options.seed << std::endl
           << "  numReads:               " << options.numReads << std::endl
           << "  useRandomSequence:      " << (options.useRandomSequence ? "true" : "false") << std::endl
           << "  randomSourceLength:     " << options.randomSourceLength << std::endl
           << "  onlyForward:            " << (options.onlyForward ? "true" : "false") << std::endl
           << "  onlyReverse:            " << (options.onlyReverse ? "true" : "false") << std::endl
           << "  outputFile:             \"" << options.outputFile << "\"" <<std::endl
           << "  samFile:                \"" << options.samFile << "\"" <<std::endl
           << "  simulateQualities:      " << (options.simulateQualities ? "true" : "false") << std::endl
           << "  generateMatePairs:      " << (options.generateMatePairs ? "true" : "false") << std::endl
           << "  libraryLengthMean:      " << options.libraryLengthMean << std::endl
           << "  libraryLengthError:     " << options.libraryLengthError << std::endl
           << "  numHaplotypes:          " << options.numHaplotypes << std::endl
           << "  haplotypeSnpRate:       " << options.haplotypeSnpRate << std::endl
           << "  haplotypeIndelRate:     " << options.haplotypeIndelRate << std::endl
           << "  haplotypeIndelRangeMin: " << options.haplotypeIndelRangeMin << std::endl
           << "  haplotypeIndelRangeMax: " << options.haplotypeIndelRangeMax << std::endl
           << "}" << std::endl;
    return stream;
}

template <typename TStream>
TStream & operator<<(TStream & stream, ReadSimulationInstruction<Global> const & inst) {
    stream << "(haplotype=" << inst.haplotype << ", contigId=" << inst.contigId << ", isForward=" << inst.isForward << ", beginPos=" << inst.beginPos << ", endPos=" << inst.endPos << ", insCount=" << inst.insCount << ", delCount=" << inst.delCount << ", editString=";
    for (unsigned i = 0; i < length(inst.editString); ++i) {
        stream << "MEID"[inst.editString[i]];
    }
    stream << ", qualities=[";
    for (unsigned i = 0; i < length(inst.qualities); ++i) {
        if (i != 0)
            stream << ", ";
        stream << inst.qualities[i];
    }
    stream << "])";
    return stream;
}

void setUpCommandLineParser(CommandLineParser & parser)
{
    addVersionLine(parser, "SeqAn read simulator 0.1");

    addTitleLine(parser, "SeqAn read simulator");
    addUsageLine(parser, "illumina [OPTIONS] SEQUENCE");
    addLine(parser, "");
    addLine(parser, "Use 'random' for the SEQUENCE file name to generate it randomly.");

    addSection(parser, "Main Options");
    
    addOption(parser, CommandLineOption("s",  "seed", "The seed for RNG.  Default: 0.", OptionType::Integer | OptionType::Label));
    addOption(parser, CommandLineOption("N",  "num-reads", "Number of reads (or mate pairs) to simulate.  Default: 1000.", OptionType::Integer));
    addOption(parser, CommandLineOption("sn", "source-length", "Length of random source sequence.  Default: 1,000,000.", OptionType::Integer));
    addOption(parser, CommandLineOption("f",  "forward-only", "Simulate from forward strand only.  Default: false.", OptionType::Bool));
    addOption(parser, CommandLineOption("r",  "reverse-only", "Simulate from reverse strand only.  Default: false.", OptionType::Bool));
    addOption(parser, CommandLineOption("o",  "output-file", "Write results to PARAM.fasta file instead of SEQUENCE.reads.fasta.  Default: \"\".", OptionType::String));
    addOption(parser, CommandLineOption("sq", "simulate-qualities", "Simulate qualities, generate FASTQ instead of FASTA.  Default: false.", OptionType::Bool));
    addOption(parser, CommandLineOption("v", "verbose", "Verbosity mode.  Default: false.", OptionType::Bool));
    addOption(parser, CommandLineOption("vv", "very-verbose", "High verbosity mode, implies verbosity mode.  Default: false.", OptionType::Bool));

    addSection(parser, "Mate-Pair Options");

    addOption(parser, CommandLineOption("ll", "library-length", "Mate-pair library length.  Default: 1000.", OptionType::Integer));
    addOption(parser, CommandLineOption("le", "library-error", "Mate-pair library tolerance.  Default: 100.", OptionType::Integer));
    addOption(parser, CommandLineOption("mp", "mate-pairs", "Enable mate pair simulation.  Default: false.", OptionType::Bool));

    addSection(parser, "Haplotype Options");

    addOption(parser, CommandLineOption("hn", "num-haplotypes", "Number of haplotypes to simulate.  Default: 1.", OptionType::Integer));
    addOption(parser, CommandLineOption("hs", "haplotype-snp-rate", "Haplotype SNP rate.  Default: 0.001.", OptionType::Double));
    addOption(parser, CommandLineOption("hi", "haplotype-indel-rate", "Haplotype indel rate.  Default: 0.001.", OptionType::Double));
    addOption(parser, CommandLineOption("hm", "haplotype-indel-range-min", "Haplotype indel size min.  Default: 1.", OptionType::Integer));
    addOption(parser, CommandLineOption("hM", "haplotype-indel-range-max", "Haplotype indel size max.  Default: 6.", OptionType::Integer));

    // Need reads type and {SEQUENCE.fasta, random}.
    requiredArguments(parser, 2);
}

template <typename TOptions>
int parseCommandLineAndCheck(TOptions & options,
                             CharString & referenceFilename,
                             CommandLineParser & parser,
                             const int argc,
                             const char * argv[]) {
    if (!parse(parser, argc, argv)) {
        if (!isSetShort(parser, "h"))
            shortHelp(parser, std::cerr);
        return 1;
    }
    if (isSetShort(parser, "h")) {
        options.showHelp = true;
        return 0;
    }

    if (isSetLong(parser, "seed"))
        getOptionValueLong(parser, "seed", options.seed);
    if (isSetLong(parser, "num-reads"))
        getOptionValueLong(parser, "num-reads", options.numReads);
    if (isSetLong(parser, "source-length"))
        getOptionValueLong(parser, "source-length", options.randomSourceLength);
    if (isSetLong(parser, "forward-only"))
        options.onlyForward = true;
    if (isSetLong(parser, "reverse-only"))
        options.onlyReverse = true;
    if (isSetLong(parser, "output-file"))
        getOptionValueLong(parser, "output-file", options.outputFile);
    if (isSetLong(parser, "simulate-qualities"))
        options.simulateQualities = true;
    if (isSetLong(parser, "verbose"))
        options.verbose = true;
    if (isSetLong(parser, "very-verbose")) {
        options.verbose = true;
        options.veryVerbose = true;
    }

    if (isSetLong(parser, "library-length-mean"))
        getOptionValueLong(parser, "library-length-mean", options.libraryLengthMean);
    if (isSetLong(parser, "library-error"))
        getOptionValueLong(parser, "library-error", options.libraryLengthError);
    if (isSetLong(parser, "mate-pairs"))
        options.generateMatePairs = true;

    if (isSetLong(parser, "num-haplotypes"))
        getOptionValueLong(parser, "num-haplotypes", options.numHaplotypes);
    if (isSetLong(parser, "haplotype-snp-rate"))
        getOptionValueLong(parser, "haplotype-snp-rate", options.haplotypeSnpRate);
    if (isSetLong(parser, "haplotype-indel-rate"))
        getOptionValueLong(parser, "haplotype-indel-rate", options.haplotypeIndelRate);
    if (isSetLong(parser, "haplotype-indel-range-min"))
        getOptionValueLong(parser, "haplotype-indel-range-min", options.haplotypeIndelRangeMin);
    if (isSetLong(parser, "haplotype-indel-range-max"))
        getOptionValueLong(parser, "haplotype-indel-range-max", options.haplotypeIndelRangeMax);


    // First argument is "illumina", second one name of reference file.
    referenceFilename = getArgumentValue(parser, 1);

    if (referenceFilename == "random")
        options.useRandomSequence = true;

    return parseCommandLineAndCheckModelSpecific(options, parser);
}

template <typename TOptions, typename TReadsTypeTag>
int simulateReads(TOptions options, CharString referenceFilename, TReadsTypeTag const &) {
    // Print options.
    std::cerr << options;
    std::cerr << "reference file: " << referenceFilename << std::endl;

    // Create a RNG with the given seed.
    RNG<MersenneTwister> rng(options.seed);

    // Load or randomly generate the reference sequence.
    FragmentStore<MyFragmentStoreConfig> fragmentStore;
    if (options.useRandomSequence) {
        referenceFilename = "random.fasta";
        std::cerr << "Generating random sequence of length " << options.randomSourceLength
                  << " to file \"" << referenceFilename << "\"." << std::endl;
        int ret = writeRandomSequence(rng, options.randomSourceLength, referenceFilename);
        if (ret != 0)
            return ret;
    }
    // Generate output file name if necessary.
    if (options.outputFile == "") {
        options.outputFile = referenceFilename;
        append(options.outputFile, ".fastq");
    }
    if (options.samFile == "") {
        options.samFile = options.outputFile;
        append(options.samFile, ".sam");
    }
    std::cerr << "Loading reference sequence from \"" << referenceFilename << "\"" << std::endl;
    if (!loadContigs(fragmentStore, referenceFilename)) {
        std::cerr << "Could not load reference sequence from " << referenceFilename << std::endl;
        return 1;
    }

    // Load and/or build the model specific parameters for the simulation.
    ModelParameters<TReadsTypeTag> modelParameters;
    int ret = simulateReadsSetupModelSpecificData(modelParameters, options);
    if (ret != 0)
        return ret;
    
    // Kick off the read generation.
    ret = simulateReadsMain(fragmentStore, rng, options, modelParameters);
    if (ret != 0)
        return ret;

    // Write out results.
    if (options.generateMatePairs) {
        // Build filename with '.1.' infix.
        CharString mateFilename = options.outputFile;
        size_t dotPos = 0;
        for (size_t i = 0; i < length(mateFilename); ++i)
            if (mateFilename[i] == '.')
                dotPos = i;
        infix(mateFilename, dotPos, dotPos + 1) = ".1.";
        // Write out first mates.
        std::cerr << "Writing resulting reads to \"" << mateFilename << "\" mates/1" << std::endl;
        StringSet<String<Dna5Q>, Dependent<> > reads;
        StringSet<CharString, Dependent<> > readNames;
        for (size_t i = 0; i < length(fragmentStore.readNameStore); i += 2) {
            appendValue(readNames, fragmentStore.readNameStore[i]);
            appendValue(reads, fragmentStore.readSeqStore[i]);
        }
        {
            std::fstream fstrm(toCString(mateFilename), std::ios_base::out);
            if (!fstrm.is_open()) {
                std::cerr << "Could not open out file \"" << mateFilename << "\"" << std::endl;
                return 1;
            }
            write(fstrm, readNames, reads, Fastq());
        }
        // Build filename with '.2.' infix.
        infix(mateFilename, dotPos + 1, dotPos + 2) = "2";
        // Write out second mates.
        std::cerr << "Writing resulting reads to \"" << mateFilename << "\" mates/2" << std::endl;
        clear(reads);
        clear(readNames);
        for (size_t i = 1; i < length(fragmentStore.readNameStore); i += 2) {
            appendValue(readNames, fragmentStore.readNameStore[i]);
            appendValue(reads, fragmentStore.readSeqStore[i]);
        }
        {
            std::fstream fstrm(toCString(mateFilename), std::ios_base::out);
            if (!fstrm.is_open()) {
                std::cerr << "Could not open out file \"" << mateFilename << "\"" << std::endl;
                return 1;
            }
            write(fstrm, readNames, reads, Fastq());
        }
    } else {
        std::cerr << "Writing resulting reads to \"" << options.outputFile << "\"" << std::endl;
        std::fstream fstrm(toCString(options.outputFile), std::ios_base::out);
        if (!fstrm.is_open()) {
            std::cerr << "Could not open out file \"" << options.outputFile << "\"" << std::endl;
            return 1;
        }
        write(fstrm, fragmentStore.readNameStore, fragmentStore.readSeqStore, Fastq());
    }
    std::cerr << "Writing SAM file \"" << options.samFile << "\"" << std::endl;
    {
        std::fstream fstrm(toCString(options.samFile), std::ios_base::out);
        if (!fstrm.is_open()) {
            std::cerr << "Could not open SAM file \"" << options.samFile << "\"" << std::endl;
            return 1;
        }
        write(fstrm, fragmentStore, SAM());
    }
    return 0;
}

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
template <typename TRNG>
int writeRandomSequence(TRNG & rng, size_t length, CharString const & fileName) {
    DnaString randomSequence;
    reserve(randomSequence, length);

    for (size_t i = 0; i < length; ++i) {
        Dna c = static_cast<Dna>(pickRandomNumber(rng, PDF<Uniform<unsigned> >(0, ValueSize<Dna>::VALUE - 1)));
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

template <typename TRNG>
void buildHaplotype(StringSet<String<Dna5, Journaled<Alloc<> > > > & haplotype,
                    FragmentStore<MyFragmentStoreConfig> & fragmentStore,
                    TRNG & rng,
                    Options<Global> const & options) {
    resize(haplotype, length(fragmentStore.contigStore), Exact());
    String<Dna5> buffer;
    reserve(buffer, options.haplotypeIndelRangeMax);

    for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i) {
        std::cout << "    contig # " << i+1 << "/" << length(fragmentStore.contigStore) << std::endl;
        clear(haplotype[i]);
        setHost(haplotype[i], fragmentStore.contigStore[i].seq);
        String<Dna5> const & contig = fragmentStore.contigStore[i].seq;
        String<Dna5, Journaled<Alloc<> > > & haplotypeContig = haplotype[i];

        // j is position in original sequence, k is position in haplotype
        for (size_t j = 0, k = 0; j < length(contig);) {
            double x = pickRandomNumber(rng, PDF<Uniform<double> >(0, 1));
            if (x < options.haplotypeSnpRate) {
                // SNP
                Dna5 c = Dna5(pickRandomNumber(rng, PDF<Uniform<int> >(0, 4));
                if (c == contig[j])
                    c = Dna5(ordValue(c) + 1);
                assignValue(haplotypeContig, k, c);
                j += 1;
                k += 1;
            } else if (x < options.haplotypeSnpRate + options.haplotypeIndelRate) {
                // Indel of random length.
                unsigned rangeLen = options.haplotypeIndelRangeMax - options.haplotypeIndelRangeMin;
                unsigned indelLen = options.haplotypeIndelRangeMin + static_cast<unsigned>(pickRandomNumber(rng, PDF<Uniform<double> >(0, 1)) * rangeLen);
                if (pickRandomNumber(rng, PDF<Uniform<double> >(0, 1)) < 0.5) {
                    // Insertion.
                    clear(buffer);
                    for (unsigned ii = 0; ii < indelLen; ++ii)
                        appendValue(buffer, Dna5(5 * pickRandomNumber(rng, PDF<Uniform<double> >(0, 1))));
                    insert(haplotypeContig, k, buffer);
                    k += indelLen;
                } else {
                    // Deletion.
                    erase(haplotypeContig, k, k + indelLen);
                    j += indelLen;
                }
            } else {
                // Match.
                j += 1;
                k += 1;
            }
        }
    }
}

// Build a read simulation instructions for a haplotype.
//
// pick a contig, probability is proportional to the length
// pick a start position, end position = start position + read length
// pick whether to match on the forward or reverse strand
// simulate edit string
// build quality values
// possibly adjust mate if left read has insert at the beginning or right read has insert at the right
template <typename TReadsTag, typename TRNG>
int buildReadSimulationInstruction(
        String<ReadSimulationInstruction<TReadsTag> > & instructions,
        TRNG & rng,
        unsigned const & haplotypeId,
        StringSet<String<Dna5, Journaled<Alloc<> > > > const & haplotype,
        String<double> const & relativeContigLengths,
        ModelParameters<TReadsTag> const & parameters,
        Options<TReadsTag> const & options)
{
    ReadSimulationInstruction<TReadsTag> inst;
    inst.haplotype = haplotypeId;

    // We have to retry simulation if the mate pair did not fit in.
    bool invalid = false;
    do {
		clear(instructions);
        invalid = false;  // By default, we do not want to repeat.
        // Pick contig id, probability is proportional to the length.
        double x = pickRandomNumber(rng, PDF<Uniform<double> >(0, 1));
        for (unsigned i = 0; i < length(relativeContigLengths); ++i) {
            if (x < relativeContigLengths[i]) {
                inst.contigId = i -1;
                break;
            }
        }
        // Pick whether on forward or reverse strand.
        if (options.onlyForward)
            inst.isForward = true;
        else if (options.onlyReverse)
            inst.isForward = false;
        else
            inst.isForward = pickRandomNumber(rng, PDF<Uniform<int> >(0, 1));
        // Pick the length in the haplotype infix the read comes from, possibly randomly.
        unsigned readLength = pickReadLength(rng, options);
        // This cannot work if the haplotype is shorter than the length of the read to simulate.
        if (length(haplotype[inst.contigId]) < readLength) {
            std::cerr << "ERROR: haplotype (== " << length(haplotype[inst.contigId]) << ") < read length!" << std::endl;
            return 1;
        }
        // Pick a start and end position.
        inst.beginPos = pickRandomNumber(rng, PDF<Uniform<size_t> >(0, length(haplotype[inst.contigId]) - readLength - 1));
        inst.endPos = inst.beginPos + readLength;
        // Simulate the read with these parameters.
        buildSimulationInstructions(inst, rng, readLength, haplotype[inst.contigId], parameters, options);
        // Append read to result list.
        appendValue(instructions, inst);

        // Maybe create a mate for this read.
        if (options.generateMatePairs) {
            // Pick a read length, possibly randomly.
            unsigned readLength = pickReadLength(rng, options);
            // Pick a library length, according to the options.
            size_t libraryLength = pickLibraryLength(rng, options);
            // Compute start and end position.
            inst.endPos = inst.beginPos + readLength + libraryLength;
            inst.beginPos = inst.endPos - readLength;
            // Verify that the mate fits right of the originally simulated read.
            size_t contigLength = length(haplotype[inst.contigId]);
            if ((inst.beginPos > contigLength) || (inst.endPos > contigLength)) {
                // Mate did not fit!  Remove previously added read and set
                // invalid to true so we repeat this simulation.
                SEQAN_ASSERT_GT(length(instructions), 0u);
                eraseBack(instructions);
                invalid = true;
                if (options.verbose)
                    std::cerr << "INFO: Mate did not fit! Repeating..." << std::endl;
                continue;
            }
            // Simulate the read with these parameters.
            buildSimulationInstructions(inst, rng, readLength, haplotype[inst.contigId], parameters, options);
            // Append read to result list.
            appendValue(instructions, inst);
        }
    } while (invalid);
	
	if (options.generateMatePairs)
		SEQAN_ASSERT_EQ(length(instructions), 2u);
	else
		SEQAN_ASSERT_EQ(length(instructions), 1u);
	
    return 0;
}

template <typename TRNG>
inline
unsigned pickLibraryLength(TRNG & rng, Options<Global> const & options)
{
    if (options.libraryLengthIsUniform) {
        // Pick uniformly.
        double minLen = options.libraryLengthMean - options.libraryLengthError;
        double maxLen = options.libraryLengthMean + options.libraryLengthError;
        double len = pickRandomNumber(rng, PDF<Uniform<double> >(minLen, maxLen));
        return static_cast<unsigned>(round(len));
    } else {
        // Pick normally distributed.
        double len = pickRandomNumber(rng, PDF<Normal>(options.libraryLengthMean, options.libraryLengthError));
        return static_cast<unsigned>(round(len));
    }
}

/**
..param.fragmentStore:FragmentStore with the contigs and where to write reads to.
..param.rng:Random number generator to use.
..param.options:Options for the simulation.
..param.errorDistribution:Error distribution, indexed by pos * 4 + ERROR_TYPE_{MATCH,MISMATCH,INSERT,DELETE}.
..param.tag:Tag for specifying reads to simulate.
*/
template <typename TRNG, typename TReadsTag, typename TOptions>
int simulateReadsMain(FragmentStore<MyFragmentStoreConfig> & fragmentStore,
                      TRNG & rng,
                      TOptions const & options,
                      ModelParameters<TReadsTag> const & parameters) {
    typedef FragmentStore<MyFragmentStoreConfig> TFragmentStore;
    typedef Value<TFragmentStore::TMatePairStore>::Type TMatePairStoreElement;

    if (options.verbose)
        std::cerr << "Simulating reads..." << std::endl;

    typedef Position<CharString>::Type TPos;

    // First, we randomly pick the haplotype for each read to be simulated.
    String<unsigned> haplotypeIds;
    reserve(haplotypeIds, options.numReads);
    for (size_t i = 0; i < options.numReads; ++i)
        appendValue(haplotypeIds, pickRandomNumber(rng, PDF<Uniform<unsigned> >(0, options.numHaplotypes - 1)));

    // We do not build all haplotypes at once since this could cost a
    // lot of memory.
    //
    // for each haplotype id
    //   simulate haplotype
    //   for each simulation instruction for this haplotype:
    //     build simulated read
    reserve(fragmentStore.readSeqStore, options.numReads, Exact());
    reserve(fragmentStore.readNameStore, options.numReads, Exact());
    char readName[1024];
    char outFileName[151];
    snprintf(outFileName, 150, "%s", toCString(options.outputFile));
    String<bool> flipped;
    for (unsigned haplotypeId = 0; haplotypeId < options.numHaplotypes; ++haplotypeId) {
        std::cerr << "Simulating for haplotype #" << haplotypeId << "..." << std::endl;
        std::cout << "  Building haplotype..." << std::endl;
        StringSet<String<Dna5, Journaled<Alloc<> > > > haplotypeContigs;
        buildHaplotype(haplotypeContigs, fragmentStore, rng, options);

        // Build partial sums over relative contig lengths so we can pick the contigs later on.
        size_t totalLength = 0;
        for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i)
            totalLength += length(fragmentStore.contigStore[i].seq);
        String<double> relativeContigLengths;
        resize(relativeContigLengths, length(fragmentStore.contigStore) + 1, Exact());
        front(relativeContigLengths) = 0.0;
        for (unsigned i = 0; i < length(fragmentStore.contigStore); ++i) {
            double l = static_cast<double>(length(fragmentStore.contigStore[i].seq));
            relativeContigLengths[i + 1] = l / totalLength;
        }
        std::partial_sum(begin(relativeContigLengths), end(relativeContigLengths), begin(relativeContigLengths));
        back(relativeContigLengths) = 1.0;

        // Simulate the reads...
        std::cerr << "  Simulating reads for haplotype #" << haplotypeId << "..." << std::endl;

//         std::cerr << "Journal: " << haplotypeContigs[0]._journalEntries << std::endl;

        for (unsigned j = 0; j < length(haplotypeIds); ++j) {
            if (haplotypeIds[j] != haplotypeId)
                continue;  // Guard against instructions on wrong haplotype.

            // Build simulation instructions.
            String<ReadSimulationInstruction<TReadsTag> > instructions;
            int res = buildReadSimulationInstruction(instructions, rng, haplotypeId, haplotypeContigs, relativeContigLengths, parameters, options);
            if (res != 0)
                return res;

            int previousMateNum = 0;
            for (unsigned k = 0; k < length(instructions); ++k) {
                ReadSimulationInstruction<TReadsTag> const & inst = instructions[k];
                // Apply simulation instructions.
                SEQAN_ASSERT_EQ(length(fragmentStore.readSeqStore), length(fragmentStore.readNameStore));
                // Cut out segment from haplotype.
                String<Dna5Q> read = infix(haplotypeContigs[inst.contigId], inst.beginPos, inst.endPos);
                applySimulationInstructions(read, rng, inst, options);
                // Append read sequence to read seq store and mate pair to read
                // name store.  This also yields the read id.  We will generate
                // and append the read name below, depending on the read id.
                unsigned readId;
                if (options.generateMatePairs)
                    readId = appendRead(fragmentStore, read, length(fragmentStore.matePairStore));
                else
                    readId = appendRead(fragmentStore, read);

                // Get expected begin/end position in the original sequence.  If we decide to flip this read later on, we will modify the WIT store.
                TPos origBeginPos = virtualToHostPosition(haplotypeContigs[inst.contigId], inst.beginPos);
                TPos origEndPos = virtualToHostPosition(haplotypeContigs[inst.contigId], inst.endPos);

                // Print info about read and haplotype.
                if (options.veryVerbose) {
                    std::cout << ",-- Read #" << readId << std::endl
                              << "| original infix:  " << infix(fragmentStore.contigStore[inst.contigId].seq, origBeginPos, origEndPos) << std::endl
                              << "| haplotype infix: " << infix(haplotypeContigs[inst.contigId], inst.beginPos, inst.endPos) << std::endl
                              << "| read:            " << read << std::endl
                              << "`-- " << std::endl;
                }

                // Generate read name.
                if (options.generateMatePairs) {
                    // Generate the mate num \in {1, 2}, randomly but consistent so two entries belonging together have different nums.  This also decides about the flipping.
                    int mateNum = 3 - previousMateNum;
                    if (readId % 2 == 0) {
                        mateNum = pickRandomNumber(rng, PDF<Uniform<int> >(1, 2));
						SEQAN_ASSERT_GEQ(mateNum, 1);
						SEQAN_ASSERT_LEQ(mateNum, 2);
                        previousMateNum = mateNum;
                        appendValue(flipped, mateNum == 2);
                    } else {
						SEQAN_ASSERT_EQ(flipped[readId - 1], mateNum == 1);
                        appendValue(flipped, mateNum == 1);
                    }
                    sprintf(readName, "%s.%09u/%d contig=%s haploid=%u length=%lu orig_begin=%lu orig_end=%lu edit_string=", outFileName, readId / 2, mateNum, toCString(fragmentStore.contigNameStore[inst.contigId]), haplotypeId, length(read), origBeginPos, origEndPos);
                } else {
                    sprintf(readName, "%s.%09u contig=%s haploid=%u length=%lu orig_begin=%lu orig_end=%lu edit_string=", outFileName, readId, toCString(fragmentStore.contigNameStore[inst.contigId]), haplotypeId, length(read), origBeginPos, origEndPos);
                }
                for (unsigned i = 0; i < length(inst.editString); ++i) {
                    char buffer[2] = "*";
                    buffer[0] = "MEID"[static_cast<int>(inst.editString[i])];
                    strcat(readName, buffer);
                }
                appendValue(fragmentStore.readNameStore, readName);

                // Tentatively add matches to aligned read store.  We will
                // maybe flip begin and end position below in the "flipping and
                // reordering" step and convert the matches to a global
                // alignment in the "convertMatchesToGlobalAlignment" call.
                if (options.generateMatePairs)
                    appendAlignedRead(fragmentStore, readId, inst.contigId, origBeginPos, origEndPos, length(fragmentStore.matePairStore));
                else
                    appendAlignedRead(fragmentStore, readId, inst.contigId, origBeginPos, origEndPos);

                // Perform flipping and reordering.
                if (options.generateMatePairs) {
                    if (readId % 2 == 1) {  // Only flip and append mate pair info after simulating second mate.
                        // Append mate pair element to fragment store's mate pair store.
                        TMatePairStoreElement matePair;
                        matePair.readId[0] = readId - 1 + flipped[readId];
                        matePair.readId[1] = readId - flipped[readId];
                        appendValue(fragmentStore.matePairStore, matePair);

                        // The first mate always comes from the forward strand.
                        append(fragmentStore.readNameStore[readId - 1], " strand=forward");
                        // The second read always comes from the reverse strand.
                        reverseComplementInPlace(fragmentStore.readSeqStore[readId]);
                        // Note: readId is also last index of aligned read store because we only have one alignment per read!
                        std::swap(fragmentStore.alignedReadStore[readId].beginPos, fragmentStore.alignedReadStore[readId].endPos);
                        append(fragmentStore.readNameStore[readId], " strand=reverse");

                        // Maybe, we write out the mates in different order, i.e. flip the entries in the stores.
                        if (flipped[readId]) {
                            SEQAN_ASSERT_TRUE(flipped[readId - 1]);
                            std::swap(fragmentStore.readSeqStore[readId - 1], fragmentStore.readSeqStore[readId]);
                            std::swap(fragmentStore.readNameStore[readId - 1], fragmentStore.readNameStore[readId]);
                            std::swap(fragmentStore.alignedReadStore[readId - 1], fragmentStore.alignedReadStore[readId]);
//                             std::cout << "flipped" << std::endl;
//                             std::cout << fragmentStore.readNameStore[readId - 1] << std::endl << fragmentStore.readNameStore[readId] << std::endl;
//                             std::cout << fragmentStore.alignedReadStore[readId - 1] << std::endl << fragmentStore.alignedReadStore[readId] << std::endl;
                        } else {
                            SEQAN_ASSERT_NOT(flipped[readId - 1]);
//                             std::cout << "not flipped" << std::endl;
//                             std::cout << fragmentStore.readNameStore[readId - 1] << std::endl << fragmentStore.readNameStore[readId] << std::endl;
//                             std::cout << fragmentStore.alignedReadStore[readId - 1] << std::endl << fragmentStore.alignedReadStore[readId] << std::endl;
                        }
                    }
                } else {
                    if (pickRandomNumber(rng, PDF<Uniform<double> >(0.0, 1.0)) < 0.5) {
                        reverseComplementInPlace(back(fragmentStore.readSeqStore));
                        append(back(fragmentStore.readNameStore), " strand=reverse");
                        // Note: readId is also last index of aligned read store because we only have one alignment per read!
                        std::swap(fragmentStore.alignedReadStore[readId].beginPos, fragmentStore.alignedReadStore[readId].endPos);
                    } else {
                        append(back(fragmentStore.readNameStore), " strand=forward");
                    }
                }
            }
			if (options.generateMatePairs) {
                // When generating mate pairs, an even number of reads is generated in each step.
 				SEQAN_ASSERT_EQ(length(fragmentStore.alignedReadStore) % 2, 0u);
				SEQAN_ASSERT_EQ(length(fragmentStore.readNameStore) % 2, 0u);
				SEQAN_ASSERT_EQ(length(fragmentStore.readSeqStore) % 2, 0u);
			}
		}
    }

    // Last but not least, convert the matches collected before to a global alignment.
    convertMatchesToGlobalAlignment(fragmentStore, Score<int, EditDistance>());
    
    if (options.verbose)
        std::cerr << "Simulated " << length(fragmentStore.readSeqStore) << " reads" << std::endl;

    return 0;
}

#endif  // READ_SIMULATOR_H_
