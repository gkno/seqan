#ifndef SEQAN_STRING_JOURNAL_INTERFACE_H
#define SEQAN_STRING_JOURNAL_INTERFACE_H

namespace seqan{

   template< typename TValue, typename TPosition >
   TValue & value( String< TValue, Journal< JournalConfig< TValue, Sloppy, False > > > & string, TPosition position ) {
      return string.getjournal().get( position );
   }
   
   template< typename TValue, typename TPosition >
   TValue const & value( String< TValue, Journal< JournalConfig< TValue, Strict, False > > > & string, TPosition position ) {
      return string.getjournal().get( position );
   }
   
   template< typename TValue, typename TPosition >
   TValue value( String< TValue, Journal< JournalConfig< TValue, Sloppy, True > > > & string, TPosition position ) {
      return string.getjournal().get( position ) + string.getjournal().shift( string.getjournal().get( position ) );
   }

   template< typename TValue, typename TPosition >
   TValue value( String< TValue, Journal< JournalConfig< TValue, Strict, True > > > & string, TPosition position ) {
      return string.getjournal().get( position ) + string.getjournal().shift( string.getjournal().get( position ) );
   }

   template< typename TValue, typename TPosition >
   TValue & value( String< TValue, Journal< JournalConfig< TValue, Sloppy, False > > > const & string, TPosition position ) {
      return string.getjournal().get( position );
   }
   
   template< typename TValue, typename TPosition >
   TValue const & value( String< TValue, Journal< JournalConfig< TValue, Strict, False > > > const & string, TPosition position ) {
      return string.getjournal().get( position );
   }
   
   template< typename TValue, typename TPosition >
   TValue value( String< TValue, Journal< JournalConfig< TValue, Sloppy, True > > > const & string, TPosition position ) {
      return string.getjournal().get( position ) + string.getjournal().shift( string.getjournal().get( position ) );
   }

   template< typename TValue, typename TPosition >
   TValue value( String< TValue, Journal< JournalConfig< TValue, Strict, True > > > const & string, TPosition position ) {
      return string.getjournal().get( position ) + string.getjournal().shift( string.getjournal().get( position ) );
   }


   template< typename TConfig >
   struct Iterator< String< typename TConfig::Type, Journal< TConfig > >, Standard > {
      typedef jiter< TConfig > Type;
   };

   template< typename TConfig >
   struct Iterator< String< typename TConfig::Type, Journal< TConfig > > const, Standard > {
      typedef jiter< TConfig > Type;
   };

   template< typename TConfig >
   struct Value< jiter< TConfig > > {
      typedef typename TConfig::Type Type;
   };

   template< typename TConfig >
   struct Value< jiter< TConfig > const > {
      typedef typename TConfig::Type Type;
   };

   template< typename TConfig >
   struct Reference< String< typename TConfig::Type, Journal< TConfig > > > {
      typedef typename TConfig::Type const & Type;
   };

   template< typename TConfig >
   struct Reference< String< typename TConfig::Type, Journal< TConfig > > const > {
      typedef typename TConfig::Type const & Type;
   };

   //////////////////////////////////////////////////////////////////////////////

   template< typename TConfig >
   struct IsContiguous< String< typename TConfig::Type, Journal< TConfig > > > {
      enum { VALUE = false };
   };

   template< typename TConfig >
   struct IsContiguous< String< typename TConfig::Type, Journal< TConfig > > const> {
      enum { VALUE = false };
   };

   template< typename TConfig >
   struct DefaultIteratorSpec< String< typename TConfig::Type, Journal< TConfig > > > {
      typedef Standard Type;
   };

   template< typename TConfig >
   struct DefaultIteratorSpec< String< typename TConfig::Type, Journal< TConfig > > const > {
      typedef Standard Type;
   };

   //////////////////////////////////////////////////////////////////////////////
   // suffix array type

   template < typename TConfig >
   struct SAValue< Index< String< typename TConfig::Type, Journal< TConfig > >, typename TConfig::TSpec > >{
      typedef typename Size< String< typename TConfig::Type, Journal< TConfig > > >::Type Type;
   };

    template < typename TConfig >
    struct Fibre< Index< String< typename TConfig::Type, Journal< TConfig > > >, Fibre_SA> {
        typedef String< typename SAValue< Index< String< typename TConfig::Type, Journal< TConfig > >, typename TConfig::TSpec > >::Type, Alloc<> > Type;
    };

    template < typename TConfig >
    struct Fibre< Index< String< typename TConfig::Type, Journal< TConfig > > >, Fibre_LCP> {
        typedef String< typename Size< String< typename TConfig::Type, Journal< TConfig > > >::Type, Alloc<> > Type;
    };

   //////////////////////////////////////////////////////////////////////////////

   template< typename TConfig, typename TString >
   void insert( size_t position, String< typename TConfig::Type, Journal< TConfig > > &journal_string, TString &insert_string ){
      journal_string.insert( position, insert_string );
   }

   template< typename TConfig >
   void insert( size_t position, String< typename TConfig::Type, Journal< TConfig > > &journal_string, typename TConfig::Type value ){
      String< typename TConfig::Type > tmpstr;              // create temporary string
      tmpstr += value;       // save in temporary string
      journal_string.insert( position, begin( tmpstr ), 1 ); // insert one TValue from temporary string :)
   }


   template< typename TConfig, typename TIterator >
   void insert( size_t position, String< typename TConfig::Type, Journal< TConfig > > &journal_string, TIterator it, size_t number ){
      journal_string.insert( position, it, number );
   }

   template< typename TConfig, typename TSource >
   void append( String< typename TConfig::Type, Journal< TConfig > > & target, TSource const& source ){
      target.append( begin( source ), length( source ) );
   }

   template< typename TConfig, typename TSource >
   void append( String< typename TConfig::Type, Journal< TConfig > > & target, TSource & source ){
      target.append( begin( source ), length( source ) );
   }

   template< typename TConfig, typename TString, typename TPosition >
   void replace( TPosition position, String< typename TConfig::Type, Journal< TConfig > > &journal_string, TString &insert_string ){
      journal_string.getjournal().replace( position, begin( insert_string ), length( insert_string ) );
   }

   template< typename TConfig, typename TIterator, typename TPosition >
   void replace( TPosition position, String< typename TConfig::Type, Journal< TConfig > > &journal_string, TIterator it, size_t number ){
      journal_string.getjournal().replace( position, it, number );
   }
   
   template< typename TConfig, typename TPosition >
   void replace( String< typename TConfig::Type, Journal< TConfig > > & journal_string, TPosition position, typename TConfig::Type value ){
      journal_string.getjournal().replace( position, value );
   }

   template< typename TConfig, typename TPosition >
   void erase( String< typename TConfig::Type, Journal< TConfig > > &journal_string, TPosition position ){
      journal_string.del( position, 1 );
   }

   template< typename TConfig, typename TPosition >
   void erase( String< typename TConfig::Type, Journal< TConfig > > &journal_string, TPosition position, TPosition position_end ){
      journal_string.del( position, position_end - position );
   }

   template< typename TConfig >
   void flatten( String< typename TConfig::Type, Journal< TConfig > > &journal_string ){
      journal_string.flatten();
   }

   template< typename TConfig >
   inline typename Size< String< typename TConfig::Type, Journal< TConfig > > >::Type resizeSpace( String< typename TConfig::Type, Journal< TConfig > > &me,
	      		size_t size,
			      size_t pos_begin,
			      size_t pos_end,
               size_t limit = 0 )
   {
      SEQAN_CHECKPOINT
	   return me.resizeSpace( size, pos_begin, pos_end, limit );
   }

   template< typename TConfig, typename TPosition, typename TExpand >
   inline typename Size< String< typename TConfig::Type, Journal< TConfig > > >::Type
   resizeSpace( String< typename TConfig::Type, Journal< TConfig > > & me,
			   typename Size< String< typename TConfig::Type, Journal< TConfig > > >::Type size,
			   TPosition pos_begin,
			   TPosition pos_end,
			   Tag< TExpand > const)
	{
		SEQAN_CHECKPOINT
	   return me.resizeSpace( size, pos_begin, pos_end );
	}

   template< typename TConfig, typename TLength, typename TExpand >
   inline typename Size< String< typename TConfig::Type, Journal< TConfig > > >::Type resize( String< typename TConfig::Type, Journal< TConfig > > & me, TLength new_length, Tag< TExpand > const ){
      std::cout << "Resizing JournalString: " << &me << " to length:" << new_length << std::endl;
      if( length( me ) > (unsigned int)new_length ){
         erase( me, (unsigned int)new_length, length( me ) - (unsigned int)new_length );
      }else{
         return me.getjournal().resize_me( new_length );
      }
      return new_length;
   }

   template< typename TConfig, typename TPos >
   inline typename IndexOperatorValue< typename TConfig::Type, typename TConfig::TSloppy, typename TConfig::TShift >::Type getValue( String< typename TConfig::Type, Journal< TConfig > > & me, TPos pos ){
      return me[pos];
   }

   template< typename TConfig, typename TPos >
   inline typename TConfig::Type const & getValue( String< typename TConfig::Type, Journal< TConfig > > const & me, TPos pos ){
      return me[pos];
   }

   template< typename TConfig, typename TSize, typename TExpand >
   inline typename Size< String< typename TConfig::Type, Journal< TConfig > > >::Type reserve( String< typename TConfig::Type, Journal< TConfig > > & me, TSize new_capacity, Tag< TExpand > const tag){
      return length( me );
   }

   template< typename TConfig, typename TSourceSpec >
   void assign( String< typename TConfig::Type, Journal< TConfig > > &target, String< typename TConfig::Type, TSourceSpec > &source ){
      target.assign_string( source );
   }

   template< typename TConfig >
   void assign( String< typename TConfig::Type, Journal< TConfig > > &target, String< typename TConfig::Type, Journal< TConfig > > &source ){
      target.assign_journal( source );
   }

   template< typename TConfig, typename TTargetSpec >
   void assign( String< typename TConfig::Type, TTargetSpec > &target, String< typename TConfig::Type, Journal< TConfig > > &source ){
      assign( target, source.get_coherent_string() );
   }

   template< typename TConfig, typename TSourceSpec >
   void assign( String< typename TConfig::Type, Journal< TConfig > > &target, String< typename TConfig::Type, TSourceSpec > const &source ){
      target.assign_string( source );
   }

   template< typename TConfig >
   void assign( String< typename TConfig::Type, Journal< TConfig > > & target, String< typename TConfig::Type, Journal< TConfig > > const & source ){
      target.assign_journal( source );
   }

   template< typename TConfig, typename TTargetSpec >
   void assign( String< typename TConfig::Type, TTargetSpec > &target, String< typename TConfig::Type, Journal< TConfig > > const &source ){
      assign( target, source.get_coherent_string() );
   }

   template< typename TConfig >
   typename Iterator< String< typename TConfig::Type, Journal< TConfig > >, Standard >::Type begin( String< typename TConfig::Type, Journal< TConfig > > &me, Standard ){
      return me.it_begin();
   }

   template< typename TConfig >
   typename Iterator< String< typename TConfig::Type, Journal< TConfig > >, Standard >::Type end( String< typename TConfig::Type, Journal< TConfig > > &me, Standard ){
      return me.it_end();
   }

   template< typename TConfig >
   typename Iterator< String< typename TConfig::Type, Journal< TConfig > >, Standard >::Type begin( String< typename TConfig::Type, Journal< TConfig > > const &me, Standard ){
      return me.it_begin();
   }

   template< typename TConfig >
   typename Iterator< String< typename TConfig::Type, Journal< TConfig > >, Standard >::Type end( String< typename TConfig::Type, Journal< TConfig > > const &me, Standard ){
      return me.it_end();
   }

   template< typename TConfig >
   size_t length( String< typename TConfig::Type, Journal< TConfig > > const &me ){
      return me.length();
   }

   template< typename TConfig >
   void const * id( String< typename TConfig::Type, Journal< TConfig > > const &me ){
      return me.getjournal().get_id();
   }

   template< typename TConfig >
   void _setLength( String< typename TConfig::Type, Journal< TConfig > > * me, size_t new_length )
   {
   SEQAN_CHECKPOINT
   	me.del( new_length, length( me ) - new_length );
   }

   template< typename TConfig >
   void _setLength( String< typename TConfig::Type, Journal< TConfig > > & me, size_t new_length )
   {
   SEQAN_CHECKPOINT
   	me.del( new_length, length( me ) - new_length );
   }

}

#endif // ndef(SEQAN_STRING_JOURNAL_INTERFACE_H)
