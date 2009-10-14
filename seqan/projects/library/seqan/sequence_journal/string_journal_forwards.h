#ifndef SEQAN_STRING_JOURNAL_FORWARDS_H
#define SEQAN_STRING_JOURNAL_FORWARDS_H

struct Node;
struct tree_visitor;

namespace seqan{
   
   struct Sloppy;
   struct Strict;

   template< typename TValue, typename TSloppySpec = Sloppy, bool ShiftSpec = false, typename TUnderlyingSpec = Alloc<>, typename TInsertSpec = Alloc<> >
   struct JournalConfig{
   
      typedef TValue Type;
      typedef TSloppySpec TSloppy;
//      static bool const Shift;
      enum{ Shift = 0 };
      typedef TUnderlyingSpec TSpec;
      typedef TInsertSpec TStringSpec;
   
   };
   
//   template< typename TValue, typename TSloppySpec, bool ShiftSpec, typename TUnderlyingSpec, typename TInsertSpec >
//   bool const JournalConfig< TValue, TSloppySpec, ShiftSpec, TUnderlyingSpec, TInsertSpec >::Shift = ShiftSpec;
   
   //TODO: Move SloppySpec into different file?
   template < typename TValue, typename TSloppySpec >
   struct SloppyValue;

   template < typename TValue >
   struct SloppyValue< TValue, Sloppy > {
      typedef TValue & Type;
   };

   template < typename TValue >
   struct SloppyValue< TValue, Strict > {
      typedef TValue const & Type;
   };

   template< typename TValue, typename TSloppySpec >
   struct journal_iterator_proxy{      
      typedef typename SloppyValue< TValue, TSloppySpec >::Type Type;
      journal_iterator_proxy( TValue & value ) : m_value( value ){};
      Type m_value;
   };
   
   template< typename TValue, typename TSloppySpec >
   typename SloppyValue< TValue, TSloppySpec >::Type operator*( journal_iterator_proxy< TValue, TSloppySpec > jtp ){
      return jtp.m_value;
   }

   //template< typename TValue, typename TSpec = Alloc<>, typename TStringSpec = Alloc<>, typename TSloppySpec = Strict >
   template< typename TConfig >
   class Journal;
   
   template < typename TValue, typename TSloppySpec, bool TShiftSpec, typename TSpec, typename TStringSpec >
   class String< TValue, Journal< JournalConfig< TValue, TSloppySpec, TShiftSpec, TSpec, TStringSpec > > >;

   template< typename TConfig >
   class jiter;
}

#endif // ndef(SEQAN_STRING_JOURNAL_FORWARDS_H)
