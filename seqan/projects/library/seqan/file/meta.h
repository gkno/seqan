#ifndef SEQAN_HEADER_FILE_META_H
#define SEQAN_HEADER_FILE_META_H



namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

// class that stores file record meta data 

//???TODO??? Lade zunaechst Daten in String, parse bei der ersten Abfrage
/**
.Class.Metadata:
..cat:Input/Output
..summary:Class that stores metadata of records in files.
..remarks:(not yet implemented)
*/
template <typename TValue, typename TFormat = void>
class Metadata
{
};


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TFormat>
inline void
clear(Metadata<TValue, TFormat> & me)
{
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
