/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de 
  ===========================================================================
  Copyright (C) 2010
  
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.
  
  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
  
  ===========================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  ===========================================================================
  Code specific to the Sanger read simulation.
  ===========================================================================
*/

// TODO(holtgrew): Implement linear quality values as in Illumina model.

#ifndef SIMULATE_SANGER_H_
#define SIMULATE_SANGER_H_

// ============================================================================
// Enums, Tags, Classes.
// ============================================================================

struct _SangerReads;
typedef Tag<_SangerReads> SangerReads;

template <>
struct Options<SangerReads> : public Options<Global>
{
    // Read Length Parameters.

    // Iff true, read lengths follow a uniform distribution, otherwise a
    // standard distribution will be used.
    bool readLengthIsUniform;

    // Average read length.
    double readLengthMean;

    // For standard distributed read lengths, this is the standard deviation,
    // for uniform read length the interval around the average to use for
    // picking the read lengths.
    double readLengthError;

    // Base Calling Error Model Parameters.

    // Mismatch probability ramp.
    double probabilityMismatchBegin;
    double probabilityMismatchEnd;

    // Insert probability ramp.
    double probabilityInsertBegin;
    double probabilityInsertEnd;

    // Delete probability ramp.
    double probabilityDeleteBegin;
    double probabilityDeleteEnd;

    Options()
            : readLengthIsUniform(false),
              readLengthMean(400),
              readLengthError(40),
              probabilityMismatchBegin(0.005),
              probabilityMismatchEnd(0.01),
              probabilityInsertBegin(0.0025),
              probabilityInsertEnd(0.005),
              probabilityDeleteBegin(0.0025),
              probabilityDeleteEnd(0.005)
    {}
};

template <>
struct ReadSimulationInstruction<SangerReads> : public ReadSimulationInstruction<Global> {
};

template<>
struct ModelParameters<SangerReads> : public ModelParameters<Global>
{
};

// ============================================================================
// Metafunctions.
// ============================================================================

// ============================================================================
// Functions.
// ============================================================================

template <typename TStream>
TStream & operator<<(TStream & stream, Options<SangerReads> const & options) {
    stream << static_cast<Options<Global> >(options);
    stream << "sanger-options {" << std::endl
           << "  readLengthIsUniform:      " << options.readLengthIsUniform << std::endl
           << "  readLengthMean:           " << options.readLengthMean << std::endl
           << "  readLengthError:          " << options.readLengthError << std::endl
           << "  probabilityMismatchBegin: " << options.probabilityMismatchBegin << std::endl
           << "  probabilityMismatchEnd:   " << options.probabilityMismatchEnd << std::endl
           << "  probabilityInsertBegin:   " << options.probabilityInsertBegin << std::endl
           << "  probabilityInsertEnd:     " << options.probabilityInsertEnd << std::endl
           << "  probabilityDeleteBegin:   " << options.probabilityInsertEnd << std::endl
           << "  probabilityDeleteEnd:     " << options.probabilityDeleteEnd << std::endl
           << "}" << std::endl;
    return stream;
}

void setUpCommandLineParser(CommandLineParser & parser,
                            SangerReads const &)
{
    setUpCommandLineParser(parser);

    addSection(parser, "Sanger Read Length Parameters");

    addOption(parser, CommandLineOption("nu",  "read-length-uniform", "If set, the read lengths are simulated with a uniform distribution, with standard distribution otherwise.  Default: false.", OptionType::Bool));
    addOption(parser, CommandLineOption("nm",  "read-length-mean", "The mean of the read lengths.  Default: 400.", OptionType::Double));
    addOption(parser, CommandLineOption("ne",  "read-length-error", "The standard deviation (for standard distribution) and interval length (for uniform distribution) for the read length.  Default: 40.", OptionType::Double));

    addSection(parser, "Sanger Error Model Parameters");

    addOption(parser, CommandLineOption("pmb",  "probability-mismatch-begin", "Probability for a mismatch at begin of read.  Default: 0.005.", OptionType::Double));
    addOption(parser, CommandLineOption("pme",  "probability-mismatch-begin", "Probability for a mismatch at end of read.  Default: 0.01.", OptionType::Double));
    addOption(parser, CommandLineOption("pib",  "probability-insert-begin", "Probability for a insert at begin of read.  Default: 0.0025.", OptionType::Double));
    addOption(parser, CommandLineOption("pie",  "probability-insert-begin", "Probability for a insert at end of read.  Default: 0.005.", OptionType::Double));
    addOption(parser, CommandLineOption("pdb",  "probability-delete-begin", "Probability for a delete at begin of read.  Default: 0.0025.", OptionType::Double));
    addOption(parser, CommandLineOption("pde",  "probability-delete-begin", "Probability for a delete at end of read.  Default: 0.005.", OptionType::Double));
}

int parseCommandLineAndCheckModelSpecific(Options<SangerReads> & options,
                                          CommandLineParser & parser)
{
    if (isSetLong(parser, "read-length-uniform"))
        options.readLengthIsUniform = true;
    if (isSetLong(parser, "read-length-mean"))
        getOptionValueLong(parser, "read-length-mean", options.readLengthMean);
    if (isSetLong(parser, "read-length-error"))
        getOptionValueLong(parser, "read-length-error", options.readLengthError);

    if (isSetLong(parser, "probability-mismatch-begin"))
        getOptionValueLong(parser, "probability-mismatch-begin", options.probabilityMismatchBegin);
    if (isSetLong(parser, "probability-mismatch-end"))
        getOptionValueLong(parser, "probability-mismatch-end", options.probabilityMismatchEnd);
    if (isSetLong(parser, "probability-insert-begin"))
        getOptionValueLong(parser, "probability-insert-begin", options.probabilityInsertBegin);
    if (isSetLong(parser, "probability-insert-end"))
        getOptionValueLong(parser, "probability-insert-end", options.probabilityInsertEnd);
    if (isSetLong(parser, "probability-delete-begin"))
        getOptionValueLong(parser, "probability-delete-begin", options.probabilityDeleteBegin);
    if (isSetLong(parser, "probability-delete-end"))
        getOptionValueLong(parser, "probability-delete-end", options.probabilityDeleteEnd);

    return 0;
}

// No model specific data for Sanger reads.
int simulateReadsSetupModelSpecificData(ModelParameters<SangerReads> & /*parameters*/,
                                        Options<SangerReads> const & /*options*/)
{
    return 0;
}

// TODO(holtgrew): Same as 454 reads!
template <typename TRNG>
inline
unsigned pickReadLength(TRNG & rng, Options<SangerReads> const & options)
{
    if (options.readLengthIsUniform) {
        // Pick uniformly.
        double minLen = options.readLengthMean - options.readLengthError;
        double maxLen = options.readLengthMean + options.readLengthError;
        double len = pickRandomNumber(rng, PDF<Uniform<double> >(minLen, maxLen));
        return static_cast<unsigned>(round(len));
    } else {
        // Pick normally distributed.
        double len = pickRandomNumber(rng, PDF<Normal>(options.readLengthMean, options.readLengthError));
        return static_cast<unsigned>(round(len));
    }
}


template <typename TRNG, typename TContig>
void buildSimulationInstructions(ReadSimulationInstruction<SangerReads> & inst, TRNG & rng, unsigned readLength, TContig const & contig, ModelParameters<SangerReads> const & /*parameters*/, Options<SangerReads> const & options)
{
    clear(inst.editString);
    reserve(inst.editString, static_cast<size_t>(1.2 * readLength), Generous());
    inst.delCount = 0;
    inst.insCount = 0;

    //
    // Build Edit String.
    //
    for (unsigned i = 0; i < readLength; /*NOP*/) {
        double x = pickRandomNumber(rng, PDF<Uniform<double> >(0, 1));
        double pos = 1.0 * i / (readLength - 1);
        double pMismatch = options.probabilityMismatchBegin + pos * (options.probabilityMismatchEnd - options.probabilityMismatchBegin);
        double pInsert   = options.probabilityInsertBegin + pos * (options.probabilityInsertEnd - options.probabilityInsertBegin);
        double pDelete   = options.probabilityDeleteBegin + pos * (options.probabilityDeleteEnd - options.probabilityDeleteBegin);
        double pMatch    = 1.0 - pMismatch - pInsert - pDelete;
        if (x < pMatch) {
            // match
            i += 1;
            appendValue(inst.editString, ERROR_TYPE_MATCH);
        } else if (x < pMatch + pMismatch) {
            // mismatch
            i += 1;
            appendValue(inst.editString, ERROR_TYPE_MISMATCH);
        } else if (x < pMatch + pMismatch + pInsert) {
            // insert
            if (length(inst.editString) > 0 && back(inst.editString == ERROR_TYPE_DELETE)) {
                inst.delCount -= 1;
                eraseBack(inst.editString);
            } else {
                i += 1;
                inst.insCount += 1;
                appendValue(inst.editString, ERROR_TYPE_INSERT);
            }
        } else {
            // Decrement string size, do not add a delete if string is
            // too short, possibly remove insert from edit string.
            if (length(inst.editString) > 0) {
                if (back(inst.editString == ERROR_TYPE_INSERT)) {
                    i -= 1;
                    inst.insCount -= 1;
                    eraseBack(inst.editString);
                } else {
                    inst.delCount += 1;
                    appendValue(inst.editString, ERROR_TYPE_DELETE);
                }
            }
        }
    }
    SEQAN_ASSERT_EQ(readLength, length(inst.editString) - inst.delCount);

    //
    // Adjust Positions.
    //

    // If the number of deletions does not equal the number of inserts
    // then we have to adjust the read positions.
    if (inst.delCount != inst.insCount) {
        int delta = static_cast<int>(inst.delCount) - static_cast<int>(inst.insCount);
        inst.endPos += delta;
        if (inst.endPos > length(contig)) {
            delta = inst.endPos - length(contig);
            inst.endPos -= delta;
            inst.beginPos -= delta;
        }
        SEQAN_ASSERT_EQ(inst.endPos - inst.beginPos + inst.insCount - inst.delCount,
                        readLength);
    }

    //
    // Quality Simulation.
    //
    // TODO(holtgrew): Need something that makes sense here.
    //

    SEQAN_ASSERT_GT(length(inst.editString), 0u);
    clear(inst.qualities);
    fill(inst.qualities, length(inst.editString), 40, Exact());
}

template <typename TRNG, typename TString>
void applySimulationInstructions(TString & read, TRNG & rng, ReadSimulationInstruction<SangerReads> const & inst, Options<SangerReads> const & /*options*/)
{
    typedef typename Value<TString>::Type TAlphabet;

    SEQAN_ASSERT_EQ(length(inst.qualities), length(inst.editString));
    
    TString tmp;
    reserve(tmp, length(read) + inst.insCount - inst.delCount);
    unsigned j = 0;
    for (unsigned i = 0; i < length(inst.editString); ++i) {
        SEQAN_ASSERT_LEQ(j, i);

        TAlphabet c;
        //int x, xold;
        switch (inst.editString[i]) {
            case ERROR_TYPE_MATCH:
                SEQAN_ASSERT_LT_MSG(j, length(read), "i = %u", i);
                appendValue(tmp, read[j]);
                assignQualityValue(back(tmp), inst.qualities[i]);
                // std::cout << i << " " << getQualityValue(back(tmp)) << " " << inst.qualities[i] << " " << convert<char>(back(tmp)) << " match" << std::endl;
                //std::cout << back(tmp) << " " << read[j] << " " << inst.qualities[i] << std::endl;
                j += 1;
                break;
            case ERROR_TYPE_MISMATCH:
                c = TAlphabet(pickRandomNumber(rng, PDF<Uniform<int> >(0, ValueSize<TAlphabet>::VALUE - 2)));  // -2, N allowed
                //xold = ordValue(c);
                SEQAN_ASSERT_LT_MSG(j, length(read), "i = %u", i);
                if (ordValue(c) >= ordValue(read[j]))
                    c = TAlphabet(ordValue(c) + 1);
                //x = ordValue(c);
                appendValue(tmp, c);
                assignQualityValue(back(tmp), inst.qualities[i]);
                // std::cout << i << " q(q_i)=" << getQualityValue(back(tmp)) << " q(i)=" << inst.qualities[i] << " char=" << convert<char>(back(tmp)) << " c_old=" << xold << " c=" << x << " r_j=" << ordValue(read[j]) << std::endl;
                // std::cout << i << " " << getQualityValue(back(tmp)) << " " << inst.qualities[i] << " " << convert<char>(back(tmp)) << " mismatch" << std::endl;
                //std::cout << "MM " << c << " " << back(tmp) << " " << inst.qualities[i] << std::endl;
                j += 1;
                break;
            case ERROR_TYPE_INSERT:
                appendValue(tmp, TAlphabet(pickRandomNumber(rng, PDF<Uniform<int> >(0, ValueSize<TAlphabet>::VALUE - 1))));  // -1 == N allowed
                assignQualityValue(back(tmp), inst.qualities[i]);
                // std::cout << i << " " << getQualityValue(back(tmp)) << " " << inst.qualities[i] << " " << convert<char>(back(tmp)) << " insertion" << std::endl;
                break;
            case ERROR_TYPE_DELETE:
                j += 1;
                break;
            default:
                SEQAN_ASSERT_FAIL("Invalid error type.");
        }
    }

    move(read, tmp);
}

#endif  // SIMULATE_SANGER_H_
