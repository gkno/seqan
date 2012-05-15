// ==========================================================================
//                      RABEMA Read Alignment Benchmark
// ==========================================================================
// Copyright (C) 2010 Manuel Holtgrewe, FU Berlin
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

// TODO(holtgrew): I think most time is spent reading the GSI file. We should be able to speed this up greatly with a more compact file format.

#ifndef APPS_RABEMA_EVALUATION_H_
#define APPS_RABEMA_EVALUATION_H_

#include <seqan/bam_io.h>
#include <seqan/basic.h>
#include <seqan/find.h>
#include <seqan/misc/misc_interval_tree.h>
#include <seqan/store.h>
#include <seqan/stream.h>

#include "fai_index.h"

#include "rabema.h"

#include "wit_store.h"
#include "find_hamming_simple_ext.h"
#include "find_myers_ukkonen_reads.h"
#include "find_myers_ukkonen_ext.h"
#include "evaluation_options.h"
#include "rabema_stats.h"
#include "ref_id_mapping.h"

// ============================================================================
// Enums, Tags, Classes.
// ============================================================================

// ----------------------------------------------------------------------------
// Helper Class CmpWitRecordLowering.
// ----------------------------------------------------------------------------

// Comparison functor for lexicographically sorting by (readId, contigId, first pos).

struct CmpWitRecordLowering
{
    bool operator()(WitRecord const & lhs, WitRecord const & rhs) const
    {
        return (lhs.readId < rhs.readId) || (lhs.readId == rhs.readId && lhs.contigId < rhs.contigId) ||
                (lhs.readId == rhs.readId && lhs.contigId == rhs.contigId && lhs.firstPos < rhs.firstPos);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function performIntervalLowering()
// ----------------------------------------------------------------------------

// Relabel intervals with the smallest distance of a contained interval.  Filter out intervals with distance > maxError.

void performIntervalLowering(String<WitRecord> & gsiRecords, int maxError)
{
    if (empty(gsiRecords))
        return;

    typedef Iterator<String<WitRecord> >::Type TIterator;
    typedef IntervalAndCargo<unsigned, unsigned> TInterval;

    // Step 1: Adjust distances.
    std::sort(begin(gsiRecords, Standard()), end(gsiRecords, Standard()), CmpWitRecordLowering());

    // Add sentinel interval.
    WitRecord sentinel(back(gsiRecords));
    sentinel.firstPos = MaxValue<size_t>::VALUE;
    sentinel.lastPos = MaxValue<size_t>::VALUE;
    // sentinel.id = MaxValue<size_t>::VALUE;
    appendValue(gsiRecords, sentinel);

    String<TInterval> openIntervals;
    unsigned i = 0;
    for (TIterator it = begin(gsiRecords, Standard()), itEnd = end(gsiRecords, Standard()); it != itEnd; ++it, ++i)
    {
        unsigned count = 0;
        for (unsigned j = 0; j < length(openIntervals); ++j)
        {
            unsigned idx = length(openIntervals) - 1 - j;
            WitRecord const thisIntervalRecord = gsiRecords[cargo(openIntervals[idx])];
            SEQAN_ASSERT_EQ(thisIntervalRecord.readName, it->readName);
            if (thisIntervalRecord.contigId != it->contigId || thisIntervalRecord.lastPos < it->firstPos)
                count += 1;
        }
        resize(openIntervals, length(openIntervals) - count);

        // Perform distance lowering for containing intervals.
        for (unsigned j = 0; j < length(openIntervals); ++j)
        {
            unsigned idx = length(openIntervals) - 1 - j;
            unsigned id = cargo(openIntervals[idx]);
            if (gsiRecords[id].distance <= maxError)
                gsiRecords[id].distance = _min(gsiRecords[id].distance, it->distance);
            else
                break;  // All containing intervals must have a greater distance.
        }

        appendValue(openIntervals, TInterval(it->firstPos, it->lastPos + 1, i));
    }

    // Step 2: Filter out intervals that are contained in intervals of lesser/equal distance.
    String<WitRecord> filteredGsiRecords;
    clear(openIntervals);
    i = 0;
    for (TIterator it = begin(gsiRecords, Standard()), itend = end(gsiRecords, Standard()); it != itend; ++it, ++i)
    {
        // Remove non-overlapping intervals on top of openInterval stack, appending to filtered intervals
        unsigned count = 0;
        for (unsigned j = 0; j < length(openIntervals); ++j) {
            unsigned idx = length(openIntervals) - 1 - j;
            WitRecord const & thisIntervalRecord = gsiRecords[cargo(openIntervals[idx])];
            SEQAN_ASSERT_EQ(thisIntervalRecord.readName, it->readName);
            if (thisIntervalRecord.contigId != it->contigId || thisIntervalRecord.lastPos < it->firstPos)
            {
                count += 1;
                unsigned startDistance = gsiRecords[cargo(openIntervals[idx])].distance;
                if (!empty(filteredGsiRecords)) {
                    if (back(filteredGsiRecords).lastPos >= leftBoundary(openIntervals[idx]))
                    {
                        if (back(filteredGsiRecords).contigId == gsiRecords[cargo(openIntervals[idx])].contigId)
                        {
                            // Assert current containing already written out.
                            SEQAN_ASSERT_GEQ(back(filteredGsiRecords).firstPos, leftBoundary(openIntervals[idx]));
                            SEQAN_ASSERT_LEQ(back(filteredGsiRecords).lastPos + 1, rightBoundary(openIntervals[idx]));
                            // Get start distance.
                            startDistance = back(filteredGsiRecords).distance + 1;
                        }
                    }
                }
                unsigned upperLimit = maxError;
                if ((unsigned)maxError < startDistance)
                    upperLimit = startDistance;
                for (unsigned i = startDistance; i <= upperLimit; ++i)
                {
                    appendValue(filteredGsiRecords, gsiRecords[cargo(openIntervals[idx])]);
                    back(filteredGsiRecords).distance = i;
                }
            }
        }
        resize(openIntervals, length(openIntervals) - count);

        // Add interval to the stack of intervals.
        if (empty(openIntervals) || gsiRecords[cargo(back(openIntervals))].distance > it->distance)
            appendValue(openIntervals, TInterval(it->firstPos, it->lastPos + 1, i));
    }
    move(gsiRecords, filteredGsiRecords);
}

// ----------------------------------------------------------------------------
// Function benchmarkReadResult()
// ----------------------------------------------------------------------------

template <typename TPatternSpec>
int benchmarkReadResult(RabemaStats & result,
                        String<BamAlignmentRecord> const & samRecords,
                        BamIOContext<StringSet<CharString > > const & bamIOContext,
                        String<WitRecord> const & gsiRecords,
                        StringSet<CharString> const & refSeqNames,
                        NameStoreCache<StringSet<CharString> > const & refSeqNamesCache,
                        StringSet<Dna5String> const & refSeqs,
                        RefIdMapping const & refIdMapping,
                        Options<EvaluateResults> const & options,
                        TPatternSpec const & /*tagPattern*/,
                        bool pairedEnd = false,
                        bool second = false)
{
    typedef IntervalAndCargo<unsigned, unsigned> TInterval;

    // Select gold standard intervals (GSI) records.
    //
    // We only select intervals that match the specification of pairedEnd, second.  Also, the intervals must be for an
    // error rate less than or equal to the maximal configured one.  In case of any-best or all-best mode, we only
    // select those of the eligible intervals with the lowest error rate.
    //
    // In case of oracle mode, we ignore the distance of the intervals in the GSI file here but use it later on.
    //
    // Start with picking the smallest distance if *-best mode.
    int smallestDistance = options.oracleWitMode ? maxValue<int>() : options.maxError;
    if (options.oracleWitMode || options.benchmarkCategory == "any-best" || options.benchmarkCategory == "all-best")
        for (unsigned i = 0; i < length(gsiRecords); ++i)
            smallestDistance = std::min(smallestDistance, gsiRecords[i].distance);
    int largestDistance = options.maxError;
    if (options.oracleWitMode && smallestDistance != maxValue<int>())
        largestDistance = smallestDistance;
    // if (options.oracleWitMode)
    //     smallestDistance = 0;
    String<WitRecord> pickedGsiRecords;
    for (unsigned i = 0; i < length(gsiRecords); ++i)
    {
        // Note: In case of oracle mode, we ignore the distance.
        if (!options.oracleWitMode && gsiRecords[i].distance > smallestDistance)
            continue;  // Skip with wrong distance.
        if (!options.oracleWitMode && gsiRecords[i].distance > options.maxError)
            continue;  // Ignore intervals with too high error rate.
        if (!pairedEnd && (gsiRecords[i].flags & WitRecord::FLAG_PAIRED))
            continue;  // Skip paired if non-paired selected.
        if (pairedEnd && !(gsiRecords[i].flags & WitRecord::FLAG_PAIRED))
            continue;  // Skip non-paired if paired selected.
        if (pairedEnd && second && !(gsiRecords[i].flags & WitRecord::FLAG_SECOND_MATE))
            continue;  // Skip if second selected but this interval is not for second.

        appendValue(pickedGsiRecords, gsiRecords[i]);

        // Get index of the sequence from GSI record contig name.
        if (!getIdByName(refSeqNames, back(pickedGsiRecords).contigName, back(pickedGsiRecords).contigId, refSeqNamesCache))
        {
            std::cerr << "ERROR: Could not find reference sequence for name " << back(pickedGsiRecords).contigName << '\n';
            return 1;
        }
    }

    // On these selected GSI records, we now perform interval lowering in case of "any-best" and "all-best".  This means
    // that an interval I with distance k_i containing an interval J with distance k_j < k_i is re-labeled with a
    // distance of k_j for the smallest k_j of all contained intervals J.
    if (!options.oracleWitMode && (options.benchmarkCategory == "any-best" || options.benchmarkCategory == "all-best"))
        performIntervalLowering(pickedGsiRecords, options.maxError);

    // Build string of intervals for each reference sequence from these filtered and lowered records.
    String<String<TInterval> > intervals;
    resize(intervals, length(refSeqs));
    String<unsigned> numIntervalsForErrorRate;
    resize(numIntervalsForErrorRate, options.maxError + 1, 0);
    String<int> intervalDistances;  // Distance of interval i.
    for (unsigned i = 0; i < length(pickedGsiRecords); ++i)
    {
        int distance = pickedGsiRecords[i].distance;
        if (!options.oracleWitMode && options.benchmarkCategory != "all" && distance != options.maxError)
            continue;  // Only search if interval distance matches -- in *-best.
        int originalDistance = pickedGsiRecords[i].originalDistance;

        appendValue(intervals[pickedGsiRecords[i].contigId], TInterval(pickedGsiRecords[i].firstPos, pickedGsiRecords[i].lastPos + 1, length(intervalDistances)));
        appendValue(intervalDistances, originalDistance);
        if (!options.oracleWitMode && options.benchmarkCategory != "any-best")
            numIntervalsForErrorRate[originalDistance] += 1;
    }
    // Marker array that states whether an interval was hit.
    String<bool> intervalHit;
    resize(intervalHit, length(intervalDistances), false);

    // Build interval trees.
    String<IntervalTree<unsigned> > intervalTrees;
    resize(intervalTrees, length(refSeqs));
    for (unsigned i = 0; i < length(intervalTrees); ++i)
        createIntervalTree(intervalTrees[i], intervals[i]);

    // One of the SAM records must have a non-empty SEQ field, extract the read seq.
    Dna5String readSeqL, readSeqR;
    bool seenL = false, seenR = false;
    for (unsigned i = 0; i < length(samRecords); ++i)
    {
        seenL |= hasFlagFirst(samRecords[i]) || (!hasFlagFirst(samRecords[i]) && !hasFlagLast(samRecords[i]));
        seenR |= hasFlagLast(samRecords[i]);
        if ((hasFlagFirst(samRecords[i]) || (!hasFlagFirst(samRecords[i]) && !hasFlagLast(samRecords[i]))) && !empty(samRecords[i].seq))
        {
            readSeqL = samRecords[i].seq;
            if (hasFlagRC(samRecords[i]))
                reverseComplement(readSeqL);
        }
        else if (hasFlagLast(samRecords[i]) && !empty(samRecords[i].seq))
        {
            readSeqR = samRecords[i].seq;
            if (hasFlagRC(samRecords[i]))
                reverseComplement(readSeqR);
        }
        if (!empty(readSeqL) && !empty(readSeqR))
            break;  // Short-circuit and break.
    }
    if (seenL && empty(readSeqL))
    {
        std::cerr << "ERROR: No alignment for query " << front(samRecords).qName << " (left-end)\n";
        return 1;
    }
    if (seenR && empty(readSeqR))
    {
        std::cerr << "ERROR: No alignment for query " << front(samRecords).qName << " (right-end)\n";
        return 1;
    }

    // Try to hit intervals.
    // std::cerr << "NUM SAM RECORDS\t" << length(samRecords) << '\n';
    for (unsigned i = 0; i < length(samRecords); ++i)
    {
        BamAlignmentRecord const & samRecord = samRecords[i];
        int seqId = refIdMapping.map[samRecord.rId];

        // Compute actual alignment score to rule out invalidly reported alignments.
        //
        // If we run in oracle mode then we ignore the actual alignment score and use the best alignment.  We only care
        // about the alignment's end position in this case.
        if (!options.oracleWitMode)
        {
            int bestDistance = minValue<int>();  // Marker for "not set yet".
            // Get best distance from NM tag if set and we are to trust it.
            if (options.trustNM)
            {
                // TODO(holtgrew): Remove const cast once we have const holders!
                BamTagsDict bamTags(const_cast<CharString &>(samRecord.tags));
                unsigned idx = 0;
                if (findTagKey(idx, bamTags, "NM"))
                {
                    if (!extractTagValue(bestDistance, bamTags, idx))
                        bestDistance = minValue<int>();  // Reset to sentinel.
                }
            }
            // Otherwise, perform a realignment.
            Dna5String readSeq;
            Dna5String contigSeq;
            if (bestDistance == minValue<int>())
            {
                if (hasFlagFirst(samRecord) || (!hasFlagFirst(samRecord) && !hasFlagLast(samRecord)))
                    readSeq = readSeqL;
                if (hasFlagLast(samRecord))
                    readSeq = readSeqR;
                int bandwidth = static_cast<int>(ceil(0.01 * options.maxError * length(readSeq)));
                __int32 beginPos = samRecord.pos - bandwidth;
                __int32 endPos = beginPos + getAlignmentLengthInRef(samRecord) - countPaddings(samRecord.cigar) + 2 * bandwidth;
                contigSeq = infix(refSeqs[seqId], beginPos, endPos);
                if (hasFlagRC(samRecord))
                    reverseComplement(contigSeq);
                Finder<Dna5String> finder(contigSeq);
                Pattern<Dna5String, TPatternSpec> pattern(readSeq, -static_cast<int>(length(readSeq)) * 1000);
                _patternMatchNOfPattern(pattern, options.matchN);
                _patternMatchNOfFinder(pattern, options.matchN);
                bool ret = setEndPosition(finder, pattern, length(contigSeq) - bandwidth);
                (void) ret;  // When run without assertions.
                SEQAN_CHECK(ret, "setEndPosition() must not fail!");
                bestDistance = -getScore(pattern);
            }

            // Skip invalid alignments.
            if (bestDistance > static_cast<int>(0.01 * options.maxError * length(readSeq)))
            {
                if (options.showSuperflousIntervals)
                {
                    std::cerr << "SUPERFLOUS/INVALID\t";
                    write2(std::cerr, samRecord, bamIOContext, Sam());
                    std::cerr << "  READ:  \t" << readSeq << '\n'
                              << "  CONTIG:\t" << contigSeq << '\n'
                              << "  DISTANCE:        \t" << bestDistance << '\n'
                              << "  ALLOWED DISTANCE:\t" << static_cast<int>(0.01 * options.maxError * length(readSeq)) << '\n';
                }
                result.invalidAlignments += 1;
                continue;
            }
        }

        // Get sequence id and last position of alignment.  We try to hit the interval with the last position (not
        // C-style end) of the read.
        unsigned lastPos = 0;
        if (!hasFlagRC(samRecord))
            lastPos = samRecord.pos + getAlignmentLengthInRef(samRecord) - countPaddings(samRecord.cigar) - 1;
        else
            lastPos = length(refSeqs[seqId]) - samRecord.pos - 1;
        
        if (options.showTryHitIntervals)
            std::cerr << "TRY HIT\tchr=" << refSeqNames[seqId] << "\tlastPos=" << lastPos << "\tqName=" << samRecord.qName << "\n";
        
        // Try to hit any interval.
        String<unsigned> result;
        findIntervals(intervalTrees[seqId], lastPos, result);
        for (unsigned i = 0; i < length(result); ++i)
            intervalHit[result[i]] = true;
    }

    // Compute number of found intervals.
    unsigned numFound = 0;
    String<unsigned> foundIntervalsForErrorRate;
    if ((int)length(foundIntervalsForErrorRate) <= largestDistance + 1)
        resize(foundIntervalsForErrorRate, largestDistance + 1, 0);
    if (options.oracleWitMode || options.benchmarkCategory == "any-best")
    {
        int bestDistance = maxValue<int>();
        for (unsigned i = 0; i < length(intervalDistances); ++i)
            if (intervalHit[i])
            {
                if (options.showHitIntervals)
                    std::cerr << "HIT\t" << pickedGsiRecords[i] << "\n";
                bestDistance = std::min(bestDistance, intervalDistances[i]);
            }
        if (bestDistance != maxValue<int>())
        {
            numFound += 1;
            foundIntervalsForErrorRate[bestDistance] += 1;
        }
    }
    else  // !options.oracleWitMode && options.benchmarkCategory in ["all-best", "all"]
    {
        for (unsigned i = 0; i < length(intervalDistances); ++i)
        {
            if (options.benchmarkCategory == "all" && intervalDistances[i] != options.maxError)
                continue;  // Only count intervals on our maximal error rate in "all" mode.
            if (intervalHit[i])
            {
                if (options.showHitIntervals)
                    std::cerr << "HIT\t" << pickedGsiRecords[i] << "\n";
                numFound += 1;
                foundIntervalsForErrorRate[intervalDistances[i]] += 1;
            }
            else
            {
                if (options.showMissedIntervals)  // inside braces for consistency with above
                    std::cerr  << "MISSED\t" << pickedGsiRecords[i] << "\n";
            }
        }
        SEQAN_ASSERT_LEQ(numFound, length(intervalDistances));
    }
    
    // Update the resulting RabemaStats.
    updateMaximalErrorRate(result, largestDistance);
    result.totalReads += 1;
    if (options.oracleWitMode || options.benchmarkCategory == "any-best")
    {
        bool found = (numFound > 0u);
        result.intervalsToFind += 1;
        result.intervalsFound += found;
        result.normalizedIntervals += found;
        int d = (smallestDistance == maxValue<int>()) ? 0 : smallestDistance;
        result.intervalsToFindForErrorRate[d] += 1;
        result.intervalsFoundForErrorRate[d] += found;
        result.normalizedIntervalsToFindForErrorRate[d] += 1;
        result.normalizedIntervalsFoundForErrorRate[d] += found;
    }
    else  // all-best or all was selected
    {
        unsigned intervalsToFind = 0;
        unsigned intervalsFound = 0;
        for (unsigned d = 0; d < length(numIntervalsForErrorRate); ++d)
        {
            // In case of "all", we only count the intervals from with maximal error rate.
            if (options.benchmarkCategory != "all" || (int)d == options.maxError)
            {
                intervalsToFind += numIntervalsForErrorRate[d];
                intervalsFound += foundIntervalsForErrorRate[d];;
                result.intervalsToFindForErrorRate[d] += numIntervalsForErrorRate[d];
                result.intervalsFoundForErrorRate[d] += foundIntervalsForErrorRate[d];
            }
        }
        result.intervalsToFind += intervalsToFind;
        result.intervalsFound += intervalsFound;
        result.normalizedIntervals += 1.0 * intervalsFound / intervalsToFind;
        for (unsigned d = 0; d < length(numIntervalsForErrorRate); ++d)
        {
            // In case of "all", we only count the intervals from with maximal error rate.
            if (options.benchmarkCategory != "all" || (int)d == options.maxError)
            {
                if (intervalsToFind > 0u)
                {
                    result.normalizedIntervalsToFindForErrorRate[d] += 1.0 * numIntervalsForErrorRate[d] / intervalsToFind;
                    result.normalizedIntervalsFoundForErrorRate[d] += 1.0 * foundIntervalsForErrorRate[d] / intervalsToFind;
                }
            }
        }
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Function compareAlignedReadsToReference()
// ----------------------------------------------------------------------------

// Stream over both the SAM and GSI file and compare the hits in the SAM file against the intervals in the GSI file.
//
// Both the SAM file and the GSI file have to be sorted by queryname for this to work.

template <typename TSamStream, typename TSamReaderSpec, typename TGsiStream, typename TGsiStreamSpec, typename TPatternSpec>
int
compareAlignedReadsToReference(RabemaStats & result,
                               RecordReader<TSamStream, TSamReaderSpec> & samReader,
                               BamIOContext<StringSet<CharString> > & bamIOContext,
                               StringSet<CharString> const & refNameStore,
                               NameStoreCache<StringSet<CharString> > const & refNameStoreCache,
                               StringSet<Dna5String> const & refSeqs,
                               RecordReader<TGsiStream, TGsiStreamSpec> & gsiReader,
                               Options<EvaluateResults> const & options,
                               TPatternSpec const & tagPattern)
{
    // Read in initial SAM/GSI records.
    BamAlignmentRecord samRecord;
    if (atEnd(samReader) || readRecord(samRecord, bamIOContext, samReader, Sam()) != 0)
    {
        std::cerr << "ERROR: Could not read first SAM record.\n";
        return 1;
    }
    WitRecord gsiRecord;
    if (atEnd(gsiReader) || readRecord(gsiRecord, gsiReader, Gsi()) != 0)
    {
        std::cerr << "ERROR: Could not read first GSI record.\n";
        return 1;
    }

    // Mapping between ref IDs from SAM file and reference sequence (from SAM file to reference sequences).
    RefIdMapping refIdMapping;
    rebuildMapping(refIdMapping, refNameStore, refNameStoreCache, nameStore(bamIOContext));

    // Current SAM and GSI records are stored in these arrays.
    String<BamAlignmentRecord> currentSamRecords;
    String<WitRecord> currentGsiRecords;

    // These flags store whether we processed the last SAM/GSI record.
    bool samDone = false, gsiDone = false;

    // The main loop: We walk over both the SAM and GSI records.
    // unsigned chunkI = 0;
    std::cerr << "Each dot corresponds to 10k processed reads.\n"
              << "\n"
              << "Progress: ";
    unsigned i = 0;
    while (!samDone || !gsiDone)
    {
        if (i > 0u && i % (100*1000) == 0u)
            std::cerr << i / 100 / 1000 << "00k";
        else if (i > 0 && i % (10*1000) == 0u)
            std::cerr << '.';
        ++i;
        
        // We process the record for the next query/read.  Since records for this next query/read might be missing in
        // both files, we need to determine which is the next one.
        CharString currentReadName;
        if (gsiDone)
            currentReadName = samRecord.qName;
        else if (samDone)
            currentReadName = gsiRecord.readName;
        else
            currentReadName = (gsiRecord.readName < samRecord.qName) ? gsiRecord.readName : samRecord.qName;
        // std::cerr << "CURRENT\t" << currentReadName << '\n';

        // These flags determine whether evaluation is run for single-end and/or paired-end reads.
        bool seenSingleEnd = false, seenPairedEnd = false;

        // Read all SAM records with the same query name.
        clear(currentSamRecords);
        // std::cerr << "SAM\t" << samRecord.qName << '\n';
        while (!samDone && samRecord.qName == currentReadName)
        {
            if (!hasFlagUnmapped(samRecord))  // Ignore records with non-aligned reads.
            {
                seenSingleEnd |= !hasFlagMultiple(samRecord);
                seenPairedEnd |= hasFlagMultiple(samRecord);
                appendValue(currentSamRecords, samRecord);
            }
            if (atEnd(samReader))
            {
                // At end of SAM File, do not read next one.
                samDone = true;
                continue;
            }
            if (readRecord(samRecord, bamIOContext, samReader, Sam()) != 0)
            {
                std::cerr << "ERROR: Could not read SAM record.\n";
                return 1;
            }
            if (samRecord.qName < currentReadName)
            {
                std::cerr << "ERROR: Wrong order in SAM file: " << samRecord.qName << " succeeds " << currentReadName << " in file.\n";
                return 1;
            }
            // Rebuild ref ID mapping if we discovered a new reference sequence.
            if (length(nameStore(bamIOContext)) != length(refIdMapping))
                rebuildMapping(refIdMapping, refNameStore, refNameStoreCache, nameStore(bamIOContext));
        }

        // Read in the next block of GSI records.
        // std::cerr << "GSI\t" << gsiRecord.readName << '\n';
        clear(currentGsiRecords);
        while (!gsiDone && gsiRecord.readName == currentReadName)
        {
            seenSingleEnd |= !(gsiRecord.flags & WitRecord::FLAG_PAIRED);
            seenPairedEnd |= (gsiRecord.flags & WitRecord::FLAG_PAIRED);
            appendValue(currentGsiRecords, gsiRecord);
            if (atEnd(gsiReader))
            {
                // At end of GSI File, do not read next one.
                gsiDone = true;
                continue;
            }
            if (readRecord(gsiRecord, gsiReader, Gsi()) != 0)
            {
                std::cerr << "ERROR: Could not read GSI record.\n";
                return 1;
            }
            if (gsiRecord.readName < currentReadName)
            {
                std::cerr << "ERROR: Wrong order in GSI file: " << gsiRecord.readName << " succeeds " << currentReadName << " in file.\n";
                return 1;
            }
        }

        // Now, compare the SAM records against the intervals stored in the GSI records.
        //
        // We collected the records for all queries.  Here, we differentiate between the different cases.
        if (seenSingleEnd)
        {
            benchmarkReadResult(result, currentSamRecords, bamIOContext, currentGsiRecords, refNameStore, refNameStoreCache, refSeqs, refIdMapping, options, tagPattern, /*pairedEnd=*/false);
        }
        if (seenPairedEnd)
        {
            benchmarkReadResult(result, currentSamRecords, bamIOContext, currentGsiRecords, refNameStore, refNameStoreCache, refSeqs, refIdMapping, options, tagPattern, /*pairedEnd=*/true, /*second=*/false);
            benchmarkReadResult(result, currentSamRecords, bamIOContext, currentGsiRecords, refNameStore, refNameStoreCache, refSeqs, refIdMapping, options, tagPattern, /*pairedEnd=*/true, /*second=*/true);
        }
    }
    std::cerr << " DONE\n";

    return 0;
}

// ----------------------------------------------------------------------------
// Function evaluateReadMapperResult()
// ----------------------------------------------------------------------------

// Entry point for the read mapper evaluation subprogram.

int evaluateReadMapperResult(Options<EvaluateResults> const & options)
{
    double startTime = 0;  // For measuring time below.

    typedef StringSet<CharString>      TNameStore;
    typedef NameStoreCache<TNameStore> TNameStoreCache;

    std::cerr << "==============================================================================\n"
              << "                RABEMA - Read Alignment BEnchMArk\n"
              << "==============================================================================\n"
              << "                        Result Comparison\n"
              << "==============================================================================\n"
              << "\n";
    std::cerr << "____OPTIONS___________________________________________________________________\n\n";

    std::cerr << "Max error rate [%]    " << options.maxError << "\n"
              << "Oracle mode           " << (options.oracleWitMode ? (char const *) "yes" : (char const *) "no") << "\n"
              << "Benchmark category    " << options.benchmarkCategory << "\n"
              << "Distance measure      " << options.distanceFunction << "\n"
              << "Match Ns              " << (options.matchN ? (char const *) "yes" : (char const *) "no") << '\n'
              << "GSI File              " << options.witFileName << '\n'
              << "SAM File              " << options.samFileName << '\n'
              << "Reference File        " << options.seqFileName << '\n'
              << "Show\n"
              << "    additional        " << (options.showAdditionalIntervals ? (char const *) "yes" : (char const *) "no") << '\n'
              << "    hit               " << (options.showHitIntervals ? (char const *) "yes" : (char const *) "no") << '\n'
              << "    missed            " << (options.showMissedIntervals ? (char const *) "yes" : (char const *) "no") << '\n'
              << "    superflous        " << (options.showSuperflousIntervals ? (char const *) "yes" : (char const *) "no") << '\n'
              << "    try hit           " << (options.showTryHitIntervals ? (char const *) "yes" : (char const *) "no") << '\n'
              << "\n";

    std::cerr << "____LOADING FILES_____________________________________________________________\n\n";

    // =================================================================
    // Prepare File I/O.
    // =================================================================

    startTime = sysTime();
    // Open reference FAI index.
    std::cerr << "Reference Index           " << options.seqFileName << ".fai ...";
    FaiIndex faiIndex;
    if (load(faiIndex, toCString(options.seqFileName)) != 0)
    {
        std::cerr << " FAILED (not fatal, we can just build it)\n";
        std::cerr << "Building Index            " << options.seqFileName << ".fai ...";
        if (buildIndex(toCString(options.seqFileName), Fai()) != 0)
        {
            std::cerr << "Could not build FAI index.\n";
            return 1;
        }
        std::cerr << " OK\n";
        std::cerr << "Reference Index           " << options.seqFileName << ".fai ...";
        if (load(faiIndex, toCString(options.seqFileName)) != 0)
        {
            std::cerr << "Could not load FAI index we just build.\n";
            return 1;
        }
        std::cerr << " OK (" << length(faiIndex.indexEntryStore) << " seqs)\n";
    }
    else
    {
        std::cerr << " OK (" << length(faiIndex.indexEntryStore) << " seqs)\n";
    }

    std::cerr << "Reference Sequences       " << options.seqFileName << " ...";
    StringSet<Dna5String> refSeqs;
    resize(refSeqs, length(faiIndex.refNameStore));
    for (unsigned i = 0; i < length(faiIndex.refNameStore); ++i)
    {
        reserve(refSeqs[i], sequenceLength(faiIndex, i), Exact());
        if (getSequence(refSeqs[i], faiIndex, i) != 0)
        {
            std::cerr << "ERROR: Could not read sequence " << faiIndex.refNameStore[i] << ".\n";
            return 0;
        }
    }
    std::cerr << " OK\n";

    // Open gold standard intervals (GSI) file and read in header.
    std::cerr << "Gold Standard Intervals   " << options.witFileName << " (header) ...";
    std::ifstream inGsi(toCString(options.witFileName), std::ios::in | std::ios::binary);
    if (!inGsi.is_open())
    {
        std::cerr << "Could not open GSI file.\n";
        return 1;
    }
    RecordReader<std::ifstream, SinglePass<> > gsiReader(inGsi);
    WitHeader gsiHeader;
    if (readRecord(gsiHeader, gsiReader, Gsi()) != 0)
    {
        std::cerr << "Could not read GSI header.\n";
        return 1;
    }
    std::cerr << " OK\n";
    // TODO(holtgrew): Do anything with GSI header?

    // Open SAM file and read in header.
    std::cerr << "Alignments                " << options.samFileName << " (header) ...";
    std::ifstream inSam(toCString(options.samFileName), std::ios_base::in | std::ios_base::binary);
    if (!inSam.is_open())
    {
        std::cerr << "Could not open SAM file." << std::endl;
        return 1;
    }
    RecordReader<std::ifstream, SinglePass<> > samReader(inSam);
    TNameStore refNameStore;
    TNameStoreCache refNameStoreCache(refNameStore);
    BamIOContext<TNameStore> samIOContext(refNameStore, refNameStoreCache);
    BamHeader samHeader;
    if (readRecord(samHeader, samIOContext, samReader, Sam()) != 0)
    {
        std::cerr << "Could not read SAM header.\n";
        return 1;
    }
    // The SAM file has to be sorted by query name.
    //
    // We do not look at the SAM header for read ordering since samtools sort does not update the SAM header.  This
    // means users would have to update it manually after sorting which is painful.  We will check SAM record order on
    // the fly.
    std::cerr << " OK\n";
    
    std::cerr << "\nTook " << sysTime() - startTime << "s\n";

    // =================================================================
    // Stream the SAM hits against gold standard intervals and check.
    // =================================================================

    std::cerr << "\n____COMPARING SAM HITS WITH INTERVALS_________________________________________\n\n";

    startTime = sysTime();
    typedef Position<WitStore::TIntervalStore>::Type TPos;
    // The result will be a list of ids to entries in witStore.
    int res = 0;
    RabemaStats result(options.maxError);
    if (options.distanceFunction == "edit")
        res = compareAlignedReadsToReference(result, samReader, samIOContext, refNameStore, refNameStoreCache, refSeqs, gsiReader, options, MyersUkkonenReads());
    else  // options.distanceFunction == "hamming"
        res = compareAlignedReadsToReference(result, samReader, samIOContext, refNameStore, refNameStoreCache, refSeqs, gsiReader, options, HammingSimple());
    if (res != 0)
        return 1;

    std::cerr << "\nTook " << sysTime() - startTime << " s\n";

    std::cerr << "\n____RESULTING STATISTICS______________________________________________________\n\n";

    std::cerr << "Note that in all-best and any-best mode, we differentiate the intervals we\n"
              << "found and those we have to find by their distance.  This is not possible in\n"
              << "all mode since a multiple lower-error intervals might be contained in an\n"
              << "higher-error interval.\n"
              << '\n'
              << "Alignments will be marked as \"invalid\" if they have a higher error rate than\n"
              << "allowed, i.e. they might become valid when increasing the allowed error rate.\n"
              << '\n'
              << "In \"all\" mode, only intervals from the maximal level are required.  In \"all-best\"\n"
              << "and \"any-best\" mode, the intervals are relabeled with the smallest distance that\n"
              << "a containing interval has.  Contained intervals are then removed.\n\n\n";
    
    write(std::cerr, result, options);

    return 0;
}

#endif  // #ifndef APPS_RABEMA_EVALUATION_H_
