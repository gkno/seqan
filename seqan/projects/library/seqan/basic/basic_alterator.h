#ifndef SEQAN_HEADER_BASIC_ALTERATOR_H
#define SEQAN_HEADER_BASIC_ALTERATOR_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec = void>
class Alterator;
	
//////////////////////////////////////////////////////////////////////////////
// Metafunctions
//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct _HolderSpec;

template <typename THost, typename TSpec>
struct _HolderSpec<Alterator<THost, TSpec> >
{
	typedef Tristate Type;
};

//////////////////////////////////////////////////////////////////////////////
// Alterator
//////////////////////////////////////////////////////////////////////////////
//Class that wrapps another class (the host)
//The alterator can be specialized independent of the host

template <typename THost, typename TSpec>
class Alterator
{
private:
	typedef typename _HolderSpec<Alterator>::Type THolderSpec;

	mutable Holder<THost, THolderSpec> data_host;

//____________________________________________________________________________

public:
	Alterator()
	{
SEQAN_CHECKPOINT
	}
	Alterator(THost & it_):
		data_host(it_)
	{
SEQAN_CHECKPOINT
	}
	Alterator(THost const & it_):
		data_host(it_)
	{
SEQAN_CHECKPOINT
	}
	Alterator(Alterator const & other_):
		data_host(other_.data_host)
	{
SEQAN_CHECKPOINT
	}

	~Alterator()
	{
SEQAN_CHECKPOINT
	}

	Alterator const & 
	operator = (THost & it_)
	{
SEQAN_CHECKPOINT
		data_host = it_;
		return *this;
	}
	Alterator const & 
	operator = (THost const & it_)
	{
SEQAN_CHECKPOINT
		data_host = it_;
		return *this;
	}
	Alterator const & 
	operator = (Alterator const & other_)
	{
SEQAN_CHECKPOINT
		data_host = other_.data_host;
		return *this;
	}

//____________________________________________________________________________

	friend inline Holder<THost, THolderSpec> &
	_dataHost(Alterator & me)
	{
SEQAN_CHECKPOINT
		return me.data_host;
	}
	friend inline Holder<THost, THolderSpec> const &
	_dataHost(Alterator const & me)
	{
SEQAN_CHECKPOINT
		return me.data_host;
	}

//____________________________________________________________________________
/*
	operator THost () const
	{
SEQAN_CHECKPOINT
		return data_host;
	}
*/
//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// METAFUNCTIONS
//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec>
struct Value<Alterator<THost, TSpec> >:
	Value<THost>
{
};
template <typename THost, typename TSpec>
struct Value<Alterator<THost, TSpec> const>:
	Value<THost>
{
};

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec>
struct GetValue<Alterator<THost, TSpec> >:
	GetValue<THost>
{
};
template <typename THost, typename TSpec>
struct GetValue<Alterator<THost, TSpec> const>:
	GetValue<THost>
{
};

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec>
struct Reference<Alterator<THost, TSpec> >:
	Reference<THost>
{
};
template <typename THost, typename TSpec>
struct Reference<Alterator<THost, TSpec> const>:
	Reference<THost>
{
};

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec>
struct Container<Alterator<THost, TSpec> >:
	Container<THost>
{
};
template <typename THost, typename TSpec>
struct Container<Alterator<THost, TSpec> const>:
	Container<THost>
{
};


//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec>
struct Size< Alterator<THost, TSpec> >:
	Size<THost>
{
};
template <typename THost, typename TSpec>
struct Size< Alterator<THost, TSpec> const>:
	Size<THost>
{
};

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec>
struct Difference< Alterator<THost, TSpec> >:
	Difference<THost>
{
};
template <typename THost, typename TSpec>
struct Difference< Alterator<THost, TSpec> const>:
	Difference<THost>
{
};

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec>
struct Host<Alterator<THost, TSpec> >
{
	typedef THost Type;
};
template <typename THost, typename TSpec>
struct Host<Alterator<THost, TSpec> const>
{
	typedef THost Type;
};

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
