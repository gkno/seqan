namespace seqan{

    unsigned int randomNumber( int upper_bound )
    {
        int value = ( (float)rand() / RAND_MAX )*( upper_bound + 1 );
        
        if( value > upper_bound ){
            value = 0;
        }
        
        return (unsigned int)floor(value);
    }

    template< typename TSpec >
    void generate_random_string( unsigned int length, seqan::String< char, TSpec > & the_string, seqan::String< char > const & alphabet_string ){
        
        seqan::resize( the_string, length );
        typename seqan::Iterator< seqan::String< char, TSpec > >::Type it = seqan::begin( the_string );
        typename seqan::Iterator< seqan::String< char, TSpec > >::Type the_end = seqan::end( the_string );
        
        while( it != the_end ){
            *it = alphabet_string[ randomNumber( seqan::length( alphabet_string ) - 1 ) ];
            ++it;
        }
        
    }

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
