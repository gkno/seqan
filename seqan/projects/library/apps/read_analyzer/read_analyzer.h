#ifndef READ_ANALYZER_READ_ANALYZER_H_
#define READ_ANALYZER_READ_ANALYZER_H_

using namespace seqan;

/*
.Class.ReadEvaluationResult
..summary:Stores base and quality counts for reads.
*/
struct ReadEvaluationResult
{
    String<size_t> baseCountOverall;  // arr[base]
    String<String<size_t> > baseCountPerPosition;  // arr[base][pos]
    String<String<String<size_t> > > qualityCountsPerPositionAndBase;  // arr[base][quality][pos]

    ReadEvaluationResult() {}
};

/**
.Function.setReadLength
..summary:Set the maximal read lenth for an ReadEvaluationResult object.  This adjusts the internal buffer sizes and has to be called before using the object.
*/
void setReadLength(ReadEvaluationResult & result, unsigned readLength)
{
    clear(result.baseCountOverall);
    fill(result.baseCountOverall, 5, 0);

    clear(result.baseCountPerPosition);
    resize(result.baseCountPerPosition, 5);
    for (unsigned i = 0; i < 5; ++i)
        fill(result.baseCountPerPosition[i], readLength, 0);

    clear(result.qualityCountsPerPositionAndBase);
    resize(result.qualityCountsPerPositionAndBase, 5);
    for (unsigned i = 0; i < 5; ++i) {
        resize(result.qualityCountsPerPositionAndBase[i], 63);
        for (unsigned j = 0; j < 63; ++j) {
            fill(result.qualityCountsPerPositionAndBase[i][j], readLength, 0);
        }
    }
}

/**
.Function.countBaseWithQualityAtPosition
..summary:Tally the given (quality, base, position) combination in the EvaluationResult object.
*/
inline void countBaseWithQualityAtPosition(ReadEvaluationResult & result, Dna5Q base, size_t position)
{
    // Count overall bases.
    result.baseCountOverall[ordValue(base)] += 1;

    // Count overall bases per position.
    result.baseCountPerPosition[ordValue(base)][position] += 1;

    // Count overall quality per base per position.
    int b = ordValue(base);
    int quality = getQualityValue(base);
    result.qualityCountsPerPositionAndBase[b][quality][position] += 1;
}

/**
.Function.performReadEvaluation
..summary:Perform an evaluation of the read base counts and qualities.
*/
template <typename TFragmentStore>
void performReadEvaluation(ReadEvaluationResult & result, TFragmentStore & fragmentStore)
{
    typedef typename TFragmentStore::TReadSeqStore TReadSeqStore;
    typedef typename Iterator<TReadSeqStore, Standard>::Type TReadSeqStoreIterator;

    for (TReadSeqStoreIterator it = begin(fragmentStore.readSeqStore); it != end(fragmentStore.readSeqStore); ++it) {
        for (unsigned i = 0; i < length(*it); ++i)
            countBaseWithQualityAtPosition(result, (*it)[i], i);
    }
}

/**
.Function.printReadEvaluationResults
..summary:Print statistical metrics about the read evaluation results.
*/
void printReadEvaluationResults(ReadEvaluationResult const & result)
{
    std::cout << "#--file:base-frequencies.dat" << std::endl;
    std::cout << "#Overall Base Frequencies" << std::endl;
    size_t sum = 0;
    for (unsigned i = 0; i < 5; ++i)
        sum += result.baseCountOverall[i];
    std::cout << "#base  ratio [%]      count" << std::endl;
    for (unsigned i = 0; i < 5; ++i) {
        std::cout << "    " << Dna5(i) << " ";
        printf("%10.2f %9lu\n", 100.0 * result.baseCountOverall[i] / sum, result.baseCountOverall[i]);
    }

    std::cout << std::endl << std::endl << "#--file:base-frequencies-position.dat" << std::endl;
    std::cout << "#Base Frequencies [%] Per Position" << std::endl;
    std::cout << "#position      A      C      G      T      N" << std::endl;
    for (unsigned i = 0; i < length(result.baseCountPerPosition[0]); ++i) {  // position
        size_t sum = 0;
        for (unsigned j = 0; j < 5; ++j)  // base
            sum += result.baseCountPerPosition[j][i];
        printf("     %4u %6.2f %6.2f %6.2f %6.2f %6.2f\n",
               i,
               100.0 * result.baseCountPerPosition[0][i] / sum,
               100.0 * result.baseCountPerPosition[1][i] / sum,
               100.0 * result.baseCountPerPosition[2][i] / sum,
               100.0 * result.baseCountPerPosition[3][i] / sum,
               100.0 * result.baseCountPerPosition[4][i] / sum);
    }

    std::cout << std::endl << std::endl << "#--file:qualities-base-position.dat" << std::endl;
    std::cout << "#Mean Quality/Std Dev Per Base Per Position" << std::endl;
    std::cout << "#position      A   sd A      C   sd C      G   sd G      T   sd T      N   sd N      *   sd *" << std::endl;
    for (unsigned i = 0; i < length(result.baseCountPerPosition[0]); ++i) {  // position
        printf("     %4u", i);
        for (unsigned j = 0; j < 5; ++j) {  // base
            // Compute mean.
            size_t sum = 0;
            size_t count = 0;
            for (unsigned k = 0; k < 63; ++k) {  // qualities
                count += result.qualityCountsPerPositionAndBase[j][k][i];
                sum += k * result.qualityCountsPerPositionAndBase[j][k][i];
            }
            double mean = 1.0 * sum / count;
            if (count == 0) {
                printf(" %6s %6s", "-", "-");
            } else {
                // Compute standard deviation.
                double devSum = 0;
                for (unsigned k = 0; k < 63; ++k) {  // qualities
                    double x = k - mean;
                    devSum += result.qualityCountsPerPositionAndBase[j][k][i] * x * x;
                }
                double stdDev = sqrt(devSum / count);
                printf(" %6.2f %6.2f", mean, stdDev);
            }
        }
        {
            // Compute mean.
            size_t sum = 0;
            size_t count = 0;
            for (unsigned j = 0; j < 5; ++j) {  // base
                for (unsigned k = 0; k < 63; ++k) {  // qualities
                    count += result.qualityCountsPerPositionAndBase[j][k][i];
                    sum += k * result.qualityCountsPerPositionAndBase[j][k][i];
                }
            }
            double mean = 1.0 * sum / count;
            if (count == 0) {
                printf(" %6s %6s", "-", "-");
            } else {
                // Compute standard deviation.
                double devSum = 0;
                for (unsigned j = 0; j < 5; ++j) {  // base
                    for (unsigned k = 0; k < 63; ++k) {  // qualities
                        double x = k - mean;
                        devSum += result.qualityCountsPerPositionAndBase[j][k][i] * x * x;
                    }
                }
                double stdDev = sqrt(devSum / count);
                printf(" %6.2f %6.2f", mean, stdDev);
            }
        }
        printf("\n");
    }
}

/*
.Class.AlignmentEvaluationResult
..summary:Stores base and quality counts for reads.
*/
struct AlignmentEvaluationResult
{
    // "Mismatch" also stores matches.

    String<size_t> insertCountsPerBase;  // arr[base]
    String<size_t> deleteCountsPerBase;  // arr[base]
    String<size_t> mismatchCountsPerMismatch;  // arr[source base * 5 + target base]

    String<String<size_t> > qualityCountsForInsertPerBase;  // arr[base][quality]
    String<String<size_t> > qualityCountsForMismatchPerBase;  // arr[src*5+tgt][quality]

    String<String<size_t> > insertCountsPerBasePerPosition;  // arr[base][pos]
    String<String<size_t> > deleteCountsPerBasePerPosition;  // arr[base][pos]
    String<String<size_t> > mismatchCountsPerMismatchPerPosition;  // arr[src*5+tgt][pos]

    String<String<String<size_t> > > qualityCountsForInsertPerBasePerPosition;  // arr[base][quality][pos]
    String<String<String<size_t> > > qualityCountsForMismatchPerMismatchPerPosition;  // arr[src*5+tgt][quality][pos]

    AlignmentEvaluationResult() {}
};

/**
.Function.setReadLength
..summary:Set the maximal read lenth for an AlignmentEvaluationResult object.  This adjusts the internal buffer sizes and has to be called before using the object.
*/
void setReadLength(AlignmentEvaluationResult & result, unsigned readLength)
{
    clear(result.insertCountsPerBase);
    fill(result.insertCountsPerBase, 5, 0);
    clear(result.deleteCountsPerBase);
    fill(result.deleteCountsPerBase, 5, 0);
    clear(result.mismatchCountsPerMismatch);
    fill(result.mismatchCountsPerMismatch, 25, 0);

    clear(result.qualityCountsForInsertPerBase);
    resize(result.qualityCountsForInsertPerBase, 5);
    for (int i = 0; i < 5; ++i)
        fill(result.qualityCountsForInsertPerBase[i], 63, 0);
    clear(result.qualityCountsForMismatchPerBase);
    resize(result.qualityCountsForMismatchPerBase, 25);
    for (int i = 0; i < 25; ++i)
        fill(result.qualityCountsForMismatchPerBase[i], 63, 0);
    
    clear(result.insertCountsPerBasePerPosition);
    resize(result.insertCountsPerBasePerPosition, 5);
    for (int i = 0; i < 5; ++i)
        fill(result.insertCountsPerBasePerPosition[i], readLength, 0);
    clear(result.deleteCountsPerBasePerPosition);
    resize(result.deleteCountsPerBasePerPosition, 5);
    for (int i = 0; i < 5; ++i)
        fill(result.deleteCountsPerBasePerPosition[i], readLength, 0);
    clear(result.mismatchCountsPerMismatchPerPosition);
    resize(result.mismatchCountsPerMismatchPerPosition, 25);
    for (int i = 0; i < 25; ++i)
        fill(result.mismatchCountsPerMismatchPerPosition[i], readLength, 0);

    clear(result.qualityCountsForInsertPerBasePerPosition);
    resize(result.qualityCountsForInsertPerBasePerPosition, 5);
    for (int i = 0; i < 5; ++i) {
        resize(result.qualityCountsForInsertPerBasePerPosition[i], 63);
        for (int j = 0; j < 63; ++j)
            fill(result.qualityCountsForInsertPerBasePerPosition[i][j], readLength, 0);
    }
    clear(result.qualityCountsForMismatchPerMismatchPerPosition);
    resize(result.qualityCountsForMismatchPerMismatchPerPosition, 25);
    for (int i = 0; i < 25; ++i) {
        resize(result.qualityCountsForMismatchPerMismatchPerPosition[i], 63);
        for (int j = 0; j < 63; ++j)
            fill(result.qualityCountsForMismatchPerMismatchPerPosition[i][j], readLength, 0);
    }
}

inline void
countInsertAtPositionWithBase(AlignmentEvaluationResult & result,
                              size_t position,
                              Dna5Q readBase)
{
    int b = ordValue(readBase);
    int q = getQualityValue(readBase);

    result.insertCountsPerBase[b] += 1;
    result.qualityCountsForInsertPerBase[b][q] += 1;
    result.insertCountsPerBasePerPosition[b][position] += 1;
    result.qualityCountsForInsertPerBasePerPosition[b][q][position] += 1;
}

inline void
countDeleteAtPositionWithBase(AlignmentEvaluationResult & result,
                              size_t position,
                              Dna5 referenceBase)
{
    int b = ordValue(referenceBase);

    result.deleteCountsPerBase[b] += 1;
    result.deleteCountsPerBasePerPosition[b][position] += 1;
}

inline void
countMismatchAtPositionWithBase(AlignmentEvaluationResult & result,
                                size_t position,
                                Dna5 referenceBase,
                                Dna5Q readBase)
{
    int s = ordValue(referenceBase);
    int t = ordValue(readBase);
    int q = getQualityValue(readBase);

    result.mismatchCountsPerMismatch[s * 5 + t] += 1;
    result.qualityCountsForMismatchPerBase[s * 5 + t][q] += 1;
    result.mismatchCountsPerMismatchPerPosition[s * 5 + t][position] += 1;
    result.qualityCountsForMismatchPerMismatchPerPosition[s * 5 + t][q][position] += 1;
}

template <typename TFragmentStore>
void performAlignmentEvaluation(AlignmentEvaluationResult & result, TFragmentStore & fragmentStore)
{
    // TODO(holtgrew): Weight info about each alignment by 1/c where c is the number of alignment positions for the read.
    typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignedReadsIter;
    typedef typename TFragmentStore::TContigStore TContigStore;
    typedef typename TFragmentStore::TReadSeqStore TReadSeqStore;
    typedef typename TFragmentStore::TReadSeq TReadSeq;
    typedef typename Value<TContigStore>::Type TContigStoreElement;
    typedef typename Value<TAlignedReadStore>::Type TAlignedReadStoreElement;
    typedef typename TAlignedReadStoreElement::TGapAnchors TReadGapAnchors;
    typedef typename TContigStoreElement::TContigSeq TContigSeq;
    typedef typename TContigStoreElement::TGapAnchors TContigGapAnchors;
    typedef Gaps<TReadSeq, AnchorGaps<TContigGapAnchors> > TReadGaps;
    typedef Gaps<TContigSeq, AnchorGaps<TContigGapAnchors> > TContigGaps;
    typedef typename Iterator<TContigGaps, Standard>::Type TContigGapAnchorsIterator;
    typedef typename Iterator<TReadGaps, Standard>::Type TReadGapAnchorsIterator;

    size_t alignedReadId = 0;
    for (TAlignedReadsIter it = begin(fragmentStore.alignedReadStore, Standard()); it != end(fragmentStore.alignedReadStore, Standard()); ++it, ++alignedReadId) {
        // Get contig and read sequences.
        TContigSeq const & contigSeq = fragmentStore.contigStore[it->contigId].seq;
        TReadSeq readSeq = fragmentStore.readSeqStore[it->readId];
        // Get gaps for contig and read.
        TContigGaps contigGaps(contigSeq, fragmentStore.contigStore[it->contigId].gaps);
        TReadGaps readGaps(readSeq, fragmentStore.alignedReadStore[alignedReadId].gaps);
        // Limit contig gaps to aligned read position.
        setBeginPosition(contigGaps, _min(it->beginPos, it->endPos));
        setEndPosition(contigGaps, _max(it->beginPos, it->endPos));
        // Reverse-complement readSeq in-place.
        unsigned readLength = length(readSeq);
        bool flipped = false;
        if (it->beginPos > it->endPos) {
            flipped = true;
            reverseComplementInPlace(readSeq);
        }

        TContigGapAnchorsIterator contigGapsIt = begin(contigGaps);
        unsigned readPos = 0;
        for (TReadGapAnchorsIterator readGapsIt = begin(readGaps); readGapsIt != end(readGaps); ++contigGapsIt, ++readGapsIt, ++readPos) {
            if (isGap(readGapsIt) && isGap(contigGapsIt))
                continue;  // Skip paddings.
            unsigned reportedPos = flipped ? readPos : readLength - readPos;
            if (isGap(readGapsIt)) {
                // Deletion
                countDeleteAtPositionWithBase(result, reportedPos, convert<Dna5Q>(*contigGapsIt));
            } else if (isGap(contigGapsIt)) {
                // Insert
                countInsertAtPositionWithBase(result, reportedPos, convert<Dna5Q>(*contigGapsIt));
            } else {
                // Match / Mismatch.
                countMismatchAtPositionWithBase(result, reportedPos, convert<Dna5Q>(*contigGapsIt), convert<Dna5Q>(*readGapsIt));
            }
        }
    }
}

void printAlignmentEvaluationResults(AlignmentEvaluationResult const & result)
{
    // Print error counts per base.
    std::cout << std::endl << std::endl << "#--file:error-counts-base.dat" << std::endl;
    printf("#base     insert    delete    mismatch     match\n");
    size_t totalInserts = 0;
    size_t totalDeletes = 0;
    size_t totalMismatches = 0;
    size_t totalMatches = 0;
    for (int i = 0; i < 5; ++i) {
        std::cout << "    " << Dna5(i);
        size_t mismatches = 0;
        size_t matches = 0;
        for (int j = 0; j < 5; ++j)
            if (i != j)
                mismatches += result.mismatchCountsPerMismatch[i * 5 + j];
            else
                matches += result.mismatchCountsPerMismatch[i * 5 + j];
        totalMismatches += mismatches;
        totalMatches += matches;
        totalInserts += result.insertCountsPerBase[i];
        totalDeletes += result.deleteCountsPerBase[i];
        printf(" %9lu %9lu    %9lu %9lu\n", result.insertCountsPerBase[i], result.deleteCountsPerBase[i], mismatches, matches);
    }
    printf("    * %9lu %9lu    %9lu %9lu\n", totalInserts, totalDeletes, totalMismatches, totalMatches);

    // Print substitution counts.
    std::cout << std::endl << std::endl << "#--file:substitution-counts.dat" << std::endl;
    std::cout << "#          " << Dna5(0) << "         " << Dna5(1) << "         " << Dna5(2) << "         " << Dna5(3) << "         " << Dna5(4) << "         *" << std::endl;
    String<size_t> sums;
    fill(sums, 6, 0);
    for (int i = 0; i < 5; ++i) {
        size_t sum = 0;
        std::cout << Dna5(i) << " ";
        for (int j = 0; j < 5; ++j) {
            sums[j] += result.mismatchCountsPerMismatch[i * 5 + j];
            sum += result.mismatchCountsPerMismatch[i * 5 + j];
            printf(" %9lu", result.mismatchCountsPerMismatch[i * 5 + j]);
        }
        printf(" %9lu", sum);
        sums[5] += sum;
        std::cout << std::endl;
    }
    std::cout << "* ";
    for (int i = 0; i < 6; ++i)
        printf(" %9lu", sums[i]);
    std::cout << std::endl;

    // Print mean/stddev of qualities per base for inserts.
    std::cout << std::endl << std::endl << "#--file:insert-qualities.dat" << std::endl;
    printf("#base   mean     sd\n");
    double totalMeanSum = 0;
    size_t totalCount = 0;
    for (int i = 0; i < 5; ++i) {
        std::cout << "    " << Dna5(i);
        size_t meanSum = 0;
        size_t count = 0;
        for (size_t j = 0; j < 63; ++j) {
            meanSum += result.qualityCountsForInsertPerBase[i][j] * j;
            count += result.qualityCountsForInsertPerBase[i][j];
        }
        totalMeanSum += meanSum;
        totalCount += count;
        if (count > 0) {
            double mean = 1.0 * meanSum / count;
            double stdDevSum = 0;
            for (size_t j = 0; j < 63; ++j)
                stdDevSum += result.qualityCountsForInsertPerBase[i][j] * (mean - j) * (mean - j);
            double stdDev = ::std::sqrt(stdDevSum / count);
            printf(" %6.2f %6.2f", mean, stdDev);
        } else {
            printf(" %6s %6s", "-", "-");
        }
        std::cout << std::endl;
    }
    std::cout << "    *";
    if (totalCount > 0) {
        double totalMean = 1.0 * totalMeanSum / totalCount;
        double totalStdDevSum = 0;
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 63; ++j)
                totalStdDevSum += result.qualityCountsForInsertPerBase[i][j] * (totalMean - j) * (totalMean - j);
        double totalStdDev = ::std::sqrt(totalStdDevSum / totalCount);
        printf(" %6.2f %6.2f", totalMean, totalStdDev);
    } else {
        printf(" %6s %6s", "-", "-");
    }
    std::cout << std::endl;

    // Print mean/stddev of qualities per base for mismatches.
    std::cout << std::endl << std::endl << "#--file:insert-mismatch.dat" << std::endl;
    printf("#base   mean     sd\n");
    totalMeanSum = 0;
    totalCount = 0;
    for (int i = 0; i < 5; ++i) {
        std::cout << "    " << Dna5(i);
        size_t meanSum = 0;
        size_t count = 0;
        for (size_t j = 0; j < 63; ++j) {
            meanSum += result.qualityCountsForInsertPerBase[i][j] * j;
            count += result.qualityCountsForInsertPerBase[i][j];
        }
        totalMeanSum += meanSum;
        totalCount += count;
        if (count > 0) {
            double mean = 1.0 * meanSum / count;
            double stdDevSum = 0;
            for (size_t j = 0; j < 63; ++j)
                stdDevSum += result.qualityCountsForInsertPerBase[i][j] * (mean - j) * (mean - j);
            double stdDev = ::std::sqrt(stdDevSum / count);
            printf(" %6.2f %6.2f", mean, stdDev);
        } else {
            printf(" %6s %6s", "-", "-");
        }
        std::cout << std::endl;
    }
    std::cout << "    *";
    if (totalCount > 0) {
        double totalMean = 1.0 * totalMeanSum / totalCount;
        double totalStdDevSum = 0;
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 63; ++j)
                totalStdDevSum += result.qualityCountsForInsertPerBase[i][j] * (totalMean - j) * (totalMean - j);
        double totalStdDev = ::std::sqrt(totalStdDevSum / totalCount);
        printf(" %6.2f %6.2f", totalMean, totalStdDev);
    } else {
        printf(" %6s %6s", "-", "-");
    }
    std::cout << std::endl;

    // Qualities per substitution.
    std::cout << std::endl << std::endl << "#--file:qualities-mismatch-base.dat" << std::endl;
    std::cout << "#Mean Quality/Std Dev Per Substitution" << std::endl;
    std::cout << "* means any, + means any but match" << std::endl;
    std::cout << "#       A   sd A      C   sd C      G   sd G      T   sd T      N   sd N      +   sd +      *   sd *" << std::endl;
    for (unsigned i = 0; i < 5; ++i) {  // source base
        std::cout << Dna5(i) << " ";
        for (unsigned j = 0; j < 5; ++j) {  // target base
            // Compute mean.
            size_t sum = 0;
            size_t count = 0;
            for (unsigned k = 0; k < 63; ++k) {  // qualities
                count += result.qualityCountsForMismatchPerBase[i * 5 + j][k];
                sum += result.qualityCountsForMismatchPerBase[i * 5 + j][k] * k;
            }
            double mean = 1.0 * sum / count;
            if (count == 0) {
                printf(" %6s %6s", "-", "-");
            } else {
                // Compute standard deviation.
                double devSum = 0;
                for (unsigned k = 0; k < 63; ++k) {  // qualities
                    double x = k - mean;
                    devSum += result.qualityCountsForMismatchPerBase[i * 5 + j][k] * x * x;
                }
                double stdDev = sqrt(devSum / count);
                printf(" %6.2f %6.2f", mean, stdDev);
            }
        }
        // Source to all.
        {
            // Compute mean.
            size_t sum = 0;
            size_t count = 0;
            size_t sumMismatch = 0;
            size_t countMismatch = 0;
            for (unsigned j = 0; j < 5; ++j) {  // base
                for (unsigned k = 0; k < 63; ++k) {  // qualities
                    if (i != j) {
                        sumMismatch += result.qualityCountsForMismatchPerBase[i * 5 + j][k] * k;
                        countMismatch += result.qualityCountsForMismatchPerBase[i * 5 + j][k];
                    }
                    sum += result.qualityCountsForMismatchPerBase[i * 5 + j][k] * k;
                    count += result.qualityCountsForMismatchPerBase[i * 5 + j][k];
                }
            }
            double mean = 1.0 * sum / count;
            double meanMismatch = 1.0 * sumMismatch / countMismatch;
            if (countMismatch == 0)
                printf(" %6s %6s", "-", "-");
            if (count == 0)
                printf(" %6s %6s", "-", "-");
            if (count != 0 || countMismatch != 0) {
                // Compute standard deviation.
                double devSum = 0;
                double devSumMismatch = 0;
                for (unsigned j = 0; j < 5; ++j) {  // base
                    for (unsigned k = 0; k < 63; ++k) {  // qualities
                        if (i != j) {
                            double x = k - meanMismatch;
                            devSumMismatch += result.qualityCountsForMismatchPerBase[i * 5 + j][k] * x * x;
                        }
                        double x = k - mean;
                        devSum += result.qualityCountsForMismatchPerBase[i * 5 + j][k] * x * x;
                    }
                }
                if (countMismatch != 0) {
                    double stdDevMismatch = sqrt(devSumMismatch / countMismatch);
                    printf(" %6.2f %6.2f", meanMismatch, stdDevMismatch);
                }
                if (count != 0) {
                    double stdDev = sqrt(devSum / count);
                    printf(" %6.2f %6.2f", mean, stdDev);
                }
            }
        }
        printf("\n");
    }
    // Mismatch to target.
    // TODO(holtgrew): Write me!
    std::cout << "+ TODO" << std::endl;
    // All to target.
    // TODO(holtgrew): Write me!
    std::cout << "* TODO" << std::endl;

    // Error probabilities per position.
    std::cout << std::endl << std::endl << "#--file:error-probabilities-position.dat" << std::endl;
    std::cout << "#position insert [%] delete [%] mismatch [%] total [%]" << std::endl;
    for (unsigned i = 0; i < length(result.insertCountsPerBasePerPosition[0]); ++i) {  // position
        size_t inserts = 0;
        size_t deletes = 0;
        size_t mismatches = 0;
        size_t matches = 0;
        for (unsigned b = 0; b < 5; ++b) {
            inserts += result.insertCountsPerBasePerPosition[b][i];
            deletes += result.deleteCountsPerBasePerPosition[b][i];
            for (unsigned a = 0; a < 5; ++a) {
                if (a == b)
                    matches += result.mismatchCountsPerMismatchPerPosition[a * 5 + b][i];
                else
                    mismatches += result.mismatchCountsPerMismatchPerPosition[a * 5 + b][i];
            }
        }
        size_t total = inserts + deletes + mismatches + matches;
        printf("%9u   %8.5f   %8.5f   %8.5f    %8.5f\n", i, 100.0 * inserts / total, 100.0 * deletes / total, 100.0 * mismatches / total, 100.0 * (inserts + deletes + mismatches) / total);
    }
}

#endif  // READ_ANALYZER_READ_ANALYZER_H_
