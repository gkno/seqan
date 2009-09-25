namespace seqan{

   //Tag to determine the type of functor to be used in the associated test
   struct TagGeneric{};

   //Fucntor that provides the template metaprogramming functionality to choose appropriate run functions for differnet tags
   template< typename TClass, typename TTag >
   struct TestFunctor{
      TestFunctor(){};
      TestFunctor( TTag t ) : tag(t) {};
      TTag tag;
   };

   //Generic info function for debug
   template< typename TClass, typename TTag >
   seqan::String<char> info( TestFunctor< TClass, TTag > & ){
      return "Generic Functor";
   }

   //Generic info function for debug
   template< typename TClass >
   seqan::String< char > info( TClass & ){
      return "Generic Instance";
   }

   //Test class stores two instances and corresponding functors defining which functions will be tested
   template< typename TClass1, typename TClass2, typename TTag1 = TagGeneric, typename TTag2 = TagGeneric >
   class Test{
   public:
      Test( TestFunctor< TClass1, TTag1 > functor_a, TClass1 & i1, TestFunctor< TClass2, TTag2 > functor_b, TClass2 & i2 )
      : functor_one( functor_a ),
        functor_two( functor_b ),
        instance_one( i1 ),
        instance_two( i2 ) {};
        
      TestFunctor< TClass1, TTag1 > functor_one;
      TestFunctor< TClass2, TTag2 > functor_two;
      
      TClass1 & instance_one;
      TClass2 & instance_two;
      
      String< Pair< bool, float > > results;
   
   private:
   };

   //formatted result output . . . early development version
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

   //general run function that returns whether the function call was successful, if that can be determined easily
   template< typename TClass, typename TTag >
   inline bool run( TestFunctor< TClass, TTag > & func, TClass & instance ){
      std::cout << "Running " << info( func ) << " with "<< info( instance ) << std::endl;
      return true;
   }

   //running the test 'count' times storing the results
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

