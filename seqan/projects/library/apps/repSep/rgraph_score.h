 /*==========================================================================
  $Id$
 ==========================================================================*/

#ifndef REPSEP_HEADER_RGRAPH_SCORE_H
#define REPSEP_HEADER_RGRAPH_SCORE_H

//////////////////////////////////////////////////////////////////////////////

struct GraphScoring{
    typedef double TScoreValue;

    TScoreValue matchScore;
    TScoreValue mismatchScore;
    TScoreValue matePairScore;

    GraphScoring(){
        matchScore = -0.3f;
        mismatchScore = 1.5f;
        matePairScore = -3.5f;    
    }
};

//////////////////////////////////////////////////////////////////////////////

GraphScoring::TScoreValue & matchScore(GraphScoring & me)
{
    return me.matchScore;
}

GraphScoring::TScoreValue const & matchScore(GraphScoring const & me)
{
    return me.matchScore;
}

//////////////////////////////////////////////////////////////////////////////


GraphScoring::TScoreValue & mismatchScore(GraphScoring & me)
{
    return me.mismatchScore;
}


GraphScoring::TScoreValue const & mismatchScore(GraphScoring const & me)
{
    return me.mismatchScore;
}

//////////////////////////////////////////////////////////////////////////////


GraphScoring::TScoreValue & matePairScore(GraphScoring & me)
{
    return me.matePairScore;
}


GraphScoring::TScoreValue const & matePairScore(GraphScoring const & me)
{
    return me.matePairScore;
}

//////////////////////////////////////////////////////////////////////////////

#endif
