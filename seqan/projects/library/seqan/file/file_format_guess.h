#ifndef SEQAN_HEADER_FILE_GUESS_H
#define SEQAN_HEADER_FILE_GUESS_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// guessFileFormat
//////////////////////////////////////////////////////////////////////////////

//guessFileFormat braucht auch data, weil die FileFormat-Klasse von TData
//abhaengig, und das ist so, weil sonst die Kombination von Templates mit
//virtuellen Funktionen nicht funktionieren wuerde.
/**
.Function.guessFileFormat:
..cat:Input/Output
..summary:Tries to determine the format of a file.
..signature:guessFileFormat(file, data)
..param.file: An input file.
..param.data: The target container.
...remarks:This container is not modified by this function.
..returns:A file format object instance that represents the determined file format.
...type:$Class.FileFormat$
..remarks:The $data$-argument is used here as a tag to determine the type of the target.
..see:Function.read
..see:Tag.File Format
*/
template <typename TFile, typename TData, typename TMeta>
inline FileFormat<TFile, TData, TMeta> &
guessFileFormat(TFile & file,
				TData & data)
{
SEQAN_CHECKPOINT
	return FileFormat<TFile, TData, TMeta, Fasta>(); //??? TODO
}

//////////////////////////////////////////////////////////////////////////////

/* DOCH NICHT:
template <typename TTarget, typename TSource>
inline void
read(TTarget & target,
	 TSource & source)
{
SEQAN_CHECKPOINT
	read(target, source, guessFileFormat(target, source));
}
*/

//////////////////////////////////////////////////////////////////////////////
} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef SEQAN_HEADER_...
