namespace seqan{

    struct TagLength{};

    struct TagErase{
        TagErase( size_t p ) : pos( p ) {};
        size_t pos;
    };

    struct TagValue{
        TagValue( size_t p ) : pos( p ) {};
        size_t pos;
    };

    template< typename TSpec >
        seqan::String< char > info( String< TSpec > & ){
        return "Instance of String";
    }

    template< typename TSpec, typename TStringSpec, typename TSloppySpec, typename TValue >
        seqan::String< char > info( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > & ){
        return "Instance of String< Journal >";
    }
    
    template< typename TClass >
    inline bool run( TestFunctor< TClass, TagLength > &, TClass & instance ){
        length( instance );
        return true;
    }

    template< typename TClass >
    inline bool run( TestFunctor< TClass, TagErase > & functor , TClass & instance ){
        erase( instance, functor.tag.pos );
        return true;
    }

    template< typename TClass >
    inline bool run( TestFunctor< TClass, TagValue > & functor , TClass & instance ){
        value( instance, functor.tag.pos );
        return true;
    }
    
}
