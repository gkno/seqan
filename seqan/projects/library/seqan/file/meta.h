#ifndef SEQAN_HEADER_FILE_META_H
#define SEQAN_HEADER_FILE_META_H



namespace SEQAN_NAMESPACE_MAIN
{

	
/*DISABLED
.Type.StringChairPair: $pair$-type for a STL container, with key and content of Type $String<char>$.
..see:Type.MetadataMap
*/
typedef std::pair < String<char>, String<char> > MetadataEntry;	// pair type for metamap

/*DISABLED
.Type.MetadataMap: The type containing meta data of a $EMBL$ or $Genbank$ sequence record.
..remarks.text: Technically a $std::multimap < MetadataEntry, MetadataEntry >$ 
*/
typedef std::multimap < String<char>, String<char> > MetadataMap;	// multimap type for meta data 

/**
.Metafunction.Metadata:
..summary:Metadata storage class.
..signature:Metadata<Value, Format>::Type
..param.Value:The value type of the file.
..param.Format:A file format tag.
...value:Tag.File Format
..returns.param.Type:Datastructure for storing the metadata of a record in files of format $Format$.
..see:Tag.File Format
..see:Class.FileFormat
..see:Class.Metadata
*/

template <typename TValue = char, typename TFormat = Default>
struct Metadata
{
	typedef MetadataMap Type;
};


/* TODO: ALTERNATIVER ANSATZ

template <typename TValue = char, typename TFormat = Default>
struct Metadata
{
	typedef String<TValue> Type;
};


TODO: data structure for more complex parsing

template <typename TValue = char, typename TFormat = Default>
struct Meta
{
	typedef String<TValue> TKey;
	typedef String<TValue> TObject;
	typedef Pair<TKey, TObject> TValue2;
	typedef String<TValue2> Type;
};
*/

//////////////////////////////////////////////////////////////////////////////


/*DISABLED
.Function.clear:
..param.object.type:Class.MetadataMap
*/

inline void
clear(MetadataMap & me)
{
SEQAN_CHECKPOINT
	me.clear(); // clear STL container
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
