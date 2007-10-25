#ifndef SEQAN_HEADER_SEQUENCE_MODEL_TYPES_H
#define SEQAN_HEADER_SEQUENCE_MODEL_TYPES_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Tags
//////////////////////////////////////////////////////////////////////////////


// OOPS = One Occurence Per Sequence
struct OOPS
{
	enum{VALUE=0};
};

//////////////////////////////////////////////////////////////////////////////

// OMOPS = One or More Occurence Per Sequence
struct OMOPS
{
	enum{VALUE=1};
};

//////////////////////////////////////////////////////////////////////////////

// ZOOPS = Zero or One Occurence Per Sequence
struct ZOOPS
{
	enum{VALUE=2};
	double threshold;

	ZOOPS():
		threshold(static_cast<double>(0.5))
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

// TCM = Two-Component-Mixture
struct TCM
{
	enum{VALUE=3};
	double threshold;

	TCM():
		threshold(static_cast<double>(0.5))
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