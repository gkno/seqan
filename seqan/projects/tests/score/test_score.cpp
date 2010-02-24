/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007-2010

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
 ============================================================================
 Tests for the SeqAn module score.
 ==========================================================================*/

// TODO(holtgrew): Check all built-in matrices again fixture matrices, there are many more fixture files.
// TODO(holtgrew): Write-out-read-in only necessary once.
// TODO(holtgrew): There is more code in score that should be tested.
// TODO(holtgrew): "testfile.txt" could also come from a macro that gives the default temporary filename.

#include <iostream>
#include <fstream>
#include <sstream>

#define SEQAN_DEBUG
#define SEQAN_TEST

#include <seqan/basic/basic_testing.h>
#include <seqan/sequence.h>
#include <seqan/score.h>

using namespace std;
using namespace seqan;

// Helper function that compares two amino acid matrices for equality.
// TODO(holtgrew): If used somewhere else, put into some place to share.
template <typename TScore1, typename TScore2>
void testCompareAAMatrices(TScore1 const & mat1, TScore2 const & mat2) {
    AminoAcid a, b;
    for (a = 'A'; a <= '*'; ++a) {
        for (b = 'A'; b <= '*'; ++b) {
            SEQAN_ASSERT_EQ(score(mat1, a, b), score(mat2, a, b));
        }
    }
}


// Compare built-in BLOSUM62 matrix with one from a file.
// TODO(holtgrew): Should be split, write-out-read-in not really necessary happens in *_pam, below, too.
SEQAN_DEFINE_TEST(test_score_matrix_score) {
    // TODO(holtgrew): It should be easier to construct these paths.
    // The path to the directory with the test's sources and fixtures.
    String<char> pathToTestSrc = SEQAN_PROGRAM_PATH;
    pathToTestSrc += "/../../../tests/score/";

    // Load fixture BLOSUM62 matrix.
    // TODO(holtgrew): Should be done in a function.
    Score<int, ScoreMatrix<> > sc;
    String<char> meta;
    String<char> pathToBlosum62(pathToTestSrc);
    append(pathToTestSrc, "BLOSUM62");
    // TODO(holtgrew): If the file does not exist, a bus error occurs, code should catch this case and print an error.
    loadScoreMatrix(sc, pathToTestSrc, meta);

    // Compare fixture BLOSUM62 matrix to built-in one.
    {
        Blosum62 blosum62;
        testCompareAAMatrices(blosum62, sc);
    }

    // Perform assertions on scores of the built-in Blosum62 matrix.
    {
        Blosum62 blosum62;
        SEQAN_ASSERT_EQ(scoreGapExtend(blosum62), -1);
        SEQAN_ASSERT_EQ(scoreGapOpen(blosum62), scoreGapExtend(blosum62));
        SEQAN_ASSERT_EQ(scoreGap(blosum62), scoreGapExtend(blosum62));
    }

    // Store and load the fixture BLOSUM62 matrix again.
    {
        FILE *fl = fopen(SEQAN_TEMP_FILENAME("testfile.txt"), "wb");
        write(fl, sc, meta);
        fclose(fl);

        Score<int, ScoreMatrix<> > sc2;
        String<char> meta2;
        loadScoreMatrix(sc2, SEQAN_TEMP_FILENAME("testfile.txt"), meta2);
        testCompareAAMatrices(sc, sc2);
        SEQAN_ASSERT_EQ(meta, meta2);
    }

    // Store and load the built-in matrix again.
    {
        FILE *fl = fopen(SEQAN_TEMP_FILENAME("testfile.txt"), "wb");
        write(fl, Blosum62());
        fclose(fl);

        Score<int, ScoreMatrix<> > sc2;
        loadScoreMatrix(sc2, SEQAN_TEMP_FILENAME("testfile.txt"));
        testCompareAAMatrices(sc2, Blosum62());
    }

    // Test setScore()
    {
        setScore(sc, 'A', '*', 100);
        SEQAN_ASSERT_EQ(score(sc, 'A', '*'), 100);
    }
}


// Write out built-in PamDayhoff matrix, read in again and compare the
// read matrix from the built-in PamDayhoff matrix again.
// TODO(holtgrew): Original comment was: Compare to http://www.bioinformatics.nl/tools/pam.html
SEQAN_DEFINE_TEST(test_score_pam) {
    //      PamJones pam;

    PamDayhoff pam;
    //Score<int, Pam<> > pam(250, -1, 0);
    //Score<int, Pam<AminoAcid, Pam_Data_Dayhoff_MDM78> > pam;

    SEQAN_ASSERT_EQ(getDist(pam), 250);
    SEQAN_ASSERT_EQ(scoreGapExtend(pam), -1);
    SEQAN_ASSERT_EQ(scoreGapOpen(pam), scoreGapExtend(pam));
    SEQAN_ASSERT_EQ(scoreGap(pam), scoreGapExtend(pam));

    write(stdout, pam);

    // Store and load built-in matrix.
    FILE * fl = fopen(SEQAN_TEMP_FILENAME("testfile.txt"), "wb");
    write(fl, pam);
    fclose(fl);

    Score<int, ScoreMatrix<> > sc;
    loadScoreMatrix(sc, SEQAN_TEMP_FILENAME("testfile.txt"));
    testCompareAAMatrices(sc, PamDayhoff());
}


SEQAN_BEGIN_TESTSUITE(test_score)
{
    SEQAN_CALL_TEST(test_score_pam);
    SEQAN_CALL_TEST(test_score_matrix_score);
}
SEQAN_END_TESTSUITE
