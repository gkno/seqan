#ifndef SEQAN_HEADER_BASIC_CONVERTER_H
#define SEQAN_HEADER_BASIC_CONVERTER_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
//Convert
//////////////////////////////////////////////////////////////////////////////

//gibt den Typ an, in den TSource konvertiert werden kann (TTarget oder TTarget &)

/**
.Metafunction.Convert:
..summary:Return type of a conversion. 
..signature:Convert<Target, Source>::Type
..param.Target:Type the object should be converted to.
..param.Source:Type of the object that should be converted to $Target$.
..returns.param.Type:Type that is returned by @Function.convert@.
...remarks:This is either $Target$ or $Target &$:
If instances of $Source: /cvsroot/ag-bio/seqan/version7/projects/library/seqan/basic/basic_converter.h,v $ can be re-interpreted as instances of $Target$,
than this metafunction returns a reference, otherwise it returns $Target$, 
that is @Function.convert@ returns a temporary.
..remarks:A constant instance of $Convert$ is (ab)used as tag argument of @Function.convertImpl@.
*/
template <typename TTarget, typename TSource = void>
struct Convert
{
	typedef TTarget Type;
	
};

//////////////////////////////////////////////////////////////////////////////
//convertImpl
//////////////////////////////////////////////////////////////////////////////
/**
.Function.convertImpl:
..hidefromindex
..cat:Alphabets
..summary:Implements @Function.convert@.
..signature:Convert convertImpl(convert, source)
..param.convert:Object that specifies the conversion.
...type:Metafunction.Convert
...remarks:A constant instance of @Metafunction.Convert@ is used to specify the conversion target.
..param.source:An object that should be converted.
..returns:$source$ converted to the type specified by convert.
...metafunction:Metafunction.Convert
..remarks:This function implements @Function.convert@. 
It is recommended to use @Function.convert@ rather than $convertImpl$.
*/
//??? Spezialisiere convertImpl, verwende convert
//??? Konversion eines einzelnen Zeichens in ein einzelnes Zeichen. Konversion von Sequenzen in Sequenzen finden wo anders statt.
//??? Kann entweder kopieren oder re-interpretieren, je nach Convert::Type
template <typename TTarget, typename T, typename TSource>
inline typename Convert<TTarget, TSource>::Type
convertImpl(Convert<TTarget, T> const,
			TSource & source)
{
	return source;
}
template <typename TTarget, typename T, typename TSource>
inline typename Convert<TTarget, TSource const>::Type
convertImpl(Convert<TTarget, T> const,
			TSource const & source)
{
	return source;
}

//////////////////////////////////////////////////////////////////////////////
//convert
//////////////////////////////////////////////////////////////////////////////
/**
.Function.convert:
..cat:Alphabets
..summary:Converts a value into another value.
..signature:Convert convert<Target>(source)
..param.Target:The type $source$ is converted to.
..param.source:An object that is converted to $Target$.
..returns:$source$ converted to $Target$.
...remarks:If $source$ can be re-interpreted as instance of $Target$, then a reference is returned.
Otherwise the function returns a temporary object. 
...metafunction:Metafunction.Convert
..remarks:This function is implemented in @Function.convertImpl@. 
Do not specialize $convert$, specialize @Function.convertImpl@ instead.
..see:Function.convertImpl
*/
template <typename TTarget, typename TSource>
inline typename Convert<TTarget, TSource>::Type
convert(TSource const & source)
{
	return convertImpl(Convert<TTarget, TSource>(), source);
}

//////////////////////////////////////////////////////////////////////////////
}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
