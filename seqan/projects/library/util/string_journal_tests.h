#include <util/test_foundry.h>

namespace seqan{

    unsigned int randomNumber( int upper_bound )
    {
        int value = ( (float)rand() / RAND_MAX )*( upper_bound + 1 );
        
        if( value > upper_bound ){
            value = 0;
        }
        
        return value;
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
    
/*template< typename TString1, typename TString2 >
    Test< TString1, TString2 > string_performance_matchup( TString1 & string_one, TString2 & string_two, size_t runs, size_t max_length ){
        Test< TString1, TString2 > test_return( TestFunctor< TString1 >(), string_one, TestFunctor< TString2 >(), string_two ); //generic test to store results
        return test_return;
    }*/

    struct TagLength{};

    struct TagErase{
        TagErase( size_t p ) : pos( p ) {};
        size_t pos;
    };
    
    struct TagEraseMany{
        TagEraseMany( String< size_t > p ) : positions( p ) {
            i = 0;
        };
        String< size_t > positions;
        
        size_t next_pos(){
            return positions[i++]; //TODO: beware of pos range, maybe use circular string kind of thing
        }
        
        size_t i;
    };

    struct TagValueMany{
        TagValueMany( String< size_t > p ) : positions( p ) {
            i = 0;
        };
        String< size_t > positions;
        
        size_t next_pos(){
            return positions[i++]; //TODO: beware of pos range, maybe use circular string kind of thing
        }
        
        size_t i;
    };

    struct TagValue{
        TagValue( size_t p ) : pos( p ) {};
        size_t pos;
    };
    
    template< typename TIndex >
    struct TagSync{
        TagSync( TIndex & i ) : index( i ) {};
        TIndex & index;
    };
    
    struct TagBuildSkew7{ };

    template< typename TSpec >
        seqan::String< char > info( String< TSpec > & ){
        return "Instance of String";
    }

    template< typename TConfig >
        seqan::String< char > info( String< typename TConfig::Type, Journal< TConfig > > & ){
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
    inline bool run( TestFunctor< TClass, TagEraseMany > & functor , TClass & instance ){
        erase( instance, functor.tag.next_pos() );
        return true;
    }

    template< typename TClass >
    inline bool run( TestFunctor< TClass, TagValue > & functor , TClass & instance ){
        value( instance, functor.tag.pos );
        return true;
    }
    
    template< typename TClass >
    inline bool run( TestFunctor< TClass, TagValueMany > & functor , TClass & instance ){
        value( instance, functor.tag.next_pos() );
        return true;
    }
    
    template< typename TClass, typename TIndex >
    inline bool run( TestFunctor< TClass, TagSync<TIndex> > & functor, TClass & instance ){
        synchronize_index( functor.tag.index, instance );
        return true;
    }
    
    template< typename TClass >
    inline bool run( TestFunctor< TClass, TagBuildSkew7 > &, TClass & instance ){
        Index< TClass > index( instance );
        indexCreate( instance, ESA_SA(), Skew7() );
        indexCreate( instance, ESA_LCP() );
        return true;
    }
    
    template< typename TString >
    class JournalTest{

    public:
    JournalTest( TString& string ): the_holder( string ) {
      srand( time(NULL) );
    };

    void random_insertions( unsigned int number, unsigned int length, seqan::String< char > const & alphabet_string ){
      if( number > 0 ){
         seqan::String< char > insertion_string;
         seqan::resize( insertion_string, length );
         for( unsigned int i = 0; i < length; ++i ){
            insertion_string[i] = alphabet_string[ randomNumber( seqan::length( alphabet_string ) - 1 ) ];
         }
         seqan::insert( randomNumber( seqan::length( seqan::value( the_holder ) ) - 1 ), seqan::value( the_holder ), insertion_string ); //insert teh shit into teh string
         random_insertions( --number, length, alphabet_string );
      }
      return; //TODO: best practice, is that useful / necessary
    }

    void random_deletions( unsigned int number, int length = 10 ){
      if( number > 0 ){
         int position = randomNumber( seqan::length( seqan::value( the_holder ) ) - length );
         seqan::erase( seqan::value( the_holder ), position, position + length ); //erase some stuff from teh string
         random_deletions( --number, length );
      }
      return; //TODO: best practice, is that useful / necessary
    }

    private:
    seqan::Holder< TString > the_holder;
    };
    
}
