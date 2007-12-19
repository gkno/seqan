#ifndef SEQAN_HEADER_SEQUENCE_MODEL_TYPES_H
#define SEQAN_HEADER_SEQUENCE_MODEL_TYPES_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Tags
//////////////////////////////////////////////////////////////////////////////


/**
.Tag.OOPS:
..summary:Represents the One Occurrence Per Sequence model.
..cat:Motif Search
..remarks:The @Tag.OOPS@ model, which was introduced by Lawrence and Reilly permits 
          exactly one motif occurrence in each sequence.
*/

struct OOPS
{
	enum{VALUE=0};
};

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.OMOPS:
..summary:Represents the One or More Occurence Per Sequence model.
..cat:Motif Search
..remarks:The @Tag.OMOPS@ model is comparable with the @Tag.TCM@ model with the one difference
          that zero occurrence in a sample sequence is not permitted.
*/

struct OMOPS
{
	enum{VALUE=1};
};

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.ZOOPS:
..summary:Represents the Zero or One Occurence Per Sequence model.
..cat:Motif Search
..remarks:The @Tag.ZOOPS@ model formulated by Bailey and Elkan permits at most one
          motif occurrence in each sequence.
*/

struct ZOOPS
{
	enum{VALUE=2};
	double threshold;

	ZOOPS():
		threshold((double)0.5)
	{
	}
	ZOOPS(double val):
		threshold(val)
	{
	}
	~ZOOPS()
	{
	}
};

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.TCM:
..summary:Represents the Two-Component-Mixture Sequence model.
..cat:Motif Search
..remarks:The @Tag.TCM@ model formulated by Bailey and Elkan permits any number pf
          non-overlapping motif occurrences per sequence.
*/

struct TCM
{
	enum{VALUE=3};
	double threshold;

	TCM():
		threshold((double)0.5)
	{
	}
	TCM(double val):
		threshold(val)
	{
	}
	~TCM()
	{
	}
};

//////////////////////////////////////////////////////////////////////////////

} // namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...