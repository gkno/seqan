#ifndef SEQAN_HEADER_BLAST_BASE_H
#define SEQAN_HEADER_BLAST_BASE_H

//SEQAN_NO_DDDOC: do not generate documentation for this file

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// Blast Report types


//Info types
struct FullInfo;

struct BasicInfo;



/**
.Spec.StoreReport:
..cat:Blast
..general:Class.Blast
..summary:Stores a Blast report.
..signature:BlastReport<TBlastHsp,StoreReport<TSpec> >
..param.TBlastHsp:The type of HSPs to be stored. See @Class.BlastHsp@
...metafunction:Metafunction.Hsp
...default:BlastHsp<BlastN,BasicInfo> 
..param.TSpec:The specializing type.
...default:BasicInfo
...type:Spec.BasicInfo
...type:Spec.FullInfo 
..include:blast.h
*/
//
//...remarks:BasicInfo only stores query name, database name and a String of all hits found. FullInfo also stores the following 
//parameters: lambda, k, h, gapped_lambda, gapped_k, gapped_h, gap_open, gap_extension; String<char> matrix; double min_expect;
//
template<typename TInfoSpec = BasicInfo>
struct StoreReport;		//stores the whole report


/**
.Spec.StreamReport:
..cat:Blast
..general:Class.Blast
..summary:Reads a Blast report from a file stream.
..signature:BlastReport<TBlastHsp,StreamReport<TFile> >
..param.TBlastHsp:The type of HSPs to be stored. See @Class.BlastHsp@
...metafunction:Metafunction.Hsp
...default:BlastHsp<BlastN,BasicInfo> 
..param.TFile:The type of the stream.
...default:std::fstream
..include:blast.h
*/
template<typename TFile = std::fstream>    //works on a stream
struct StreamReport;




//////////////////////////////////////////////////////////////////////////////
//Blast Meta functions

template<typename T>
struct Hit;

template<typename T>
struct Hsp;


//////////////////////////////////////////////////////////////////////////////
// Blast Tag

struct TagBlast_;
typedef Tag<TagBlast_> const Blast;


//////////////////////////////////////////////////////////////////////////////
// Blat Tag
struct TagBlat_;
typedef Tag<TagBlat_> const Blat;

//////////////////////////////////////////////////////////////////////////////





struct TagBlastN_;
struct TagMegaBlast_;
struct TagBlastP_;
struct TagBlastX_;
struct TagTBlastN_;
struct TagTBlastX_;

template<typename TSpec = TagBlastN_>
class NucleotideBlast{
public:
	NucleotideBlast(){}
	~NucleotideBlast(){}
};

template<typename TSpec = TagBlastP_>
class ProteinBlast{
public:
	ProteinBlast(){}
	~ProteinBlast(){}
};

typedef NucleotideBlast<TagBlastN_> BlastN;
typedef NucleotideBlast<TagMegaBlast_> MegaBlast;

typedef ProteinBlast<TagBlastP_> BlastP;
typedef ProteinBlast<TagBlastX_> BlastX;
typedef ProteinBlast<TagTBlastN_> TBlastN;
typedef ProteinBlast<TagTBlastX_> TBlastX;



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
