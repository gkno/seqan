#ifndef STORE_CONFIG_H_
#define STORE_CONFIG_H_

// Fragment Store Configuration.
//
// We change the default fragment store configuration to use normal
// string sets and not concat string sets for the reads since we need
// to be able to swap read sequences.

namespace seqan {
struct MyFragmentStoreConfig {};

template<>
struct FragmentStoreConfig<MyFragmentStoreConfig>
{
	typedef String<Dna5Q>	TReadSeq;
	typedef String<Dna5Q>	TContigSeq;
	
	typedef double			TMean;
	typedef double			TStd;
	typedef signed char		TMappingQuality;
		
	typedef void    TReadStoreElementSpec;
	typedef Owner<Default> TReadSeqStoreSpec;
	typedef void    TMatePairStoreElementSpec;
	typedef void    TLibraryStoreElementSpec;
	typedef void    TContigStoreElementSpec;
	typedef void    TContigFileSpec;
	typedef void    TAlignedReadStoreElementSpec;
	typedef Owner<Default>	TAlignedReadTagStoreSpec;
	typedef void    TAnnotationStoreElementSpec;
};
}  // namespace seqan

#endif  // #ifndef STORE_CONFIG_H_
