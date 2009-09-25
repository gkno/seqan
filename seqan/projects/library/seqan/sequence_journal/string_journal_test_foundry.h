#include <cstdlib>

namespace seqan{

   struct TagGeneric{};
   struct TagLength{};
   struct TagErase{
      TagErase( size_t p ) : pos( p ) {};
   size_t pos;
   };

   template< typename TClass, typename TTag >
   struct TestFunctor{
      TestFunctor(){};
      TestFunctor( TTag t ) : tag(t) {};
      TTag tag;
   };

   template< typename TClass, typename TTag >
   seqan::String<char> info( TestFunctor< TClass, TTag > & ){
      return "Generic Functor";
   }

   template< typename TClass >
   seqan::String< char > info( TClass & ){
      return "Generic Instance";
   }

   template< typename TSpec >
   seqan::String< char > info( String< TSpec > & ){
      return "Instance of String";
   }
   
   template< typename TSpec, typename TStringSpec, typename TSloppySpec, typename TValue >
   seqan::String< char > info( String< TValue, Journal< TValue, TSpec, TStringSpec, TSloppySpec > > & ){
      return "Instance of String< Journal >";
   }

   template< typename TClass1, typename TClass2, typename TTag1 = TagGeneric, typename TTag2 = TagGeneric >
   class Test{
   public:
      Test( TestFunctor< TClass1, TTag1 > functor_a, TClass1 & i1, TestFunctor< TClass2, TTag2 > functor_b, TClass2 & i2 )
       : functor_one( functor_a ), functor_two( functor_b ), instance_one( i1 ), instance_two( i2 ) {};
      TestFunctor< TClass1, TTag1 > functor_one;
      TestFunctor< TClass2, TTag2 > functor_two;
      TClass1 & instance_one;
      TClass2 & instance_two;
      String< Pair< bool, float > > results;
   private:
   };

   template< typename TTest >   
   seqan::String< char > results_formatted( TTest & test ){
      typename Iterator< String< Pair< bool, float > > >::Type it_results = begin( test.results );
      String< char > output = "";
      std::ostringstream helper;
      while( it_results != end( test.results ) ){
         if( (*it_results).i1 ){
            helper << (*it_results).i2;
            output +=  helper.str();
            helper.str("");
         }else{
            output +=  "error";
         }
         
         ++it_results;
         
         output += "\t";
         
         if( (*it_results).i1 ){
            helper << (*it_results).i2;
            output +=  helper.str();
            helper.str("");
         }else{
            output +=  "error";
         }
         
         output += "\n";
         ++it_results;
      }
      return output;
   }

   template< typename TClass, typename TTag >
   inline bool run( TestFunctor< TClass, TTag > & func, TClass & instance ){
      std::cout << "Running " << info( func ) << " with "<< info( instance ) << std::endl;
      return true;
   }
   
   template< typename TClass >
   inline bool run( TestFunctor< TClass, TagLength > &, TClass & instance ){
      length( instance );
//      std::cout << "length" << std::endl;
      return true;
   }
   
   template< typename TClass >
   inline bool run( TestFunctor< TClass, TagErase > & functor , TClass & instance ){
      erase( instance, functor.tag.pos );
//      std::cout << "erase" << std::endl;
      return true;
   }

   template< typename TTest, typename TCount >
   inline void run( TTest & test, TCount count){
      Pair< bool, float > result;
      
      SEQAN_PROTIMESTART(_test);
      
      while( count > 0 ){
         result.i1 = run( test.functor_one, test.instance_one );
         result.i2 = SEQAN_PROTIMEUPDATE(_test);
         
         appendValue( test.results, result, Generous() );
         
         result.i1 = run( test.functor_two, test.instance_two );
         result.i2 = SEQAN_PROTIMEUPDATE(_test);
         
         appendValue( test.results, result, Generous() );
         
         --count;
      }
   }

}



   unsigned int randomNumber( int upper_bound )
   {
      int value = ( (float)rand() / RAND_MAX )*( upper_bound + 1 );
      if( value > upper_bound )
      {
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
         //std::cout << *it << " ";
         ++it;
      }
      //std::cout << std::endl;
      return; //TODO: best practice, is that useful / necessary
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

