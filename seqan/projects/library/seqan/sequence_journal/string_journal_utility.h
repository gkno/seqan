namespace seqan{

    template< typename TString > //TODO: rewrite to encompass general Container functionality
    struct suffix_compare_functor{

        suffix_compare_functor( TString & string ) : m_string( string ){}

        template< typename TPos1, typename TPos2 >
        bool operator()( TPos1 & first, TPos2 & second ){
            return suffix( m_string, first ) < suffix( m_string, second );
        }

    private:
        TString & m_string;
    };

    template< typename TIteratorA, typename TIteratorB >
    inline bool _suffix_bigger( TIteratorA & it_a, TIteratorB & it_b ){
        std::cout << "Bigger?" << std::endl;
        while( *it_a == *it_b ){
            std::cout << *it_a << "==" << *it_b << std::endl;
            ++it_a;
            ++it_b;
        }
        std::cout << *it_a << "> " << *it_b << " ?" << std::endl;
        return *it_a > *it_b;
    }

    template< typename TIteratorA, typename TIteratorB >
    inline size_t _lcp_length( TIteratorA it_a, TIteratorB it_b ){
        size_t lcplength = 0;
        //std::cout << "LCP ";
        while( *it_a == *it_b ){
            ++it_a;
            ++it_b;
            ++lcplength;
            //std::cout << '.';
        }
        //std::cout << lcplength << std::endl;
        return lcplength;
    }

   template< typename TString, typename TPos >
   inline int get_shift( TString const & limits, TPos position ){
      size_t low = 0;
      size_t high = length( limits );
      size_t mid = (size_t)floor( high / 2 );
      bool seeking = true;
      while( seeking ){
        if( limits[mid].i1 > position ){
            high = mid;
            mid = low + (size_t)floor( ( high - low ) / 2 );
        }else{
            if( limits[mid + 1].i1 > position || mid + 1 == high ){
                seeking = false;
            }else{
                low = mid;
                mid = low + (size_t)floor( ( high - low ) / 2 );
            }
        }
      }
      return limits[mid].i2;
   }

   template< typename TValue1, typename TValue2, typename TTag >
   void print_pairs( String< Pair< TValue1, TValue2, TTag > > const & string ){
      for( unsigned int i = 0; i < length( string ); ++i ){
         std::cout << "< " << string[i].i1 << ", " << string[i].i2 << " > ";
      }
      std::cout << std::endl;
   }
   
//#define NDEBUG_SYNC    
        template< typename TIndex, typename TValue, typename TSloppySpec, typename TSAInv >
    inline void synchronize_index( TIndex & index, String< TValue, Journal< JournalConfig< TValue, TSloppySpec, False > > > & string, TSAInv & sa_inv ){
        typedef String< TValue, Journal< JournalConfig< TValue, TSloppySpec, False > > > TString;
    
        typedef typename Position< TString >::Type TPos;
        typedef String< typename SAValue< TIndex >::Type > TSA;
        typedef String< TPos > TLCP;

        String< typename SAValue< TIndex >::Type, Journal< JournalConfig< typename SAValue< TIndex >::Type > > > fibre_sa( indexSA( index ) );
        String< typename SAValue< TIndex >::Type, Journal< JournalConfig< typename SAValue< TIndex >::Type > > > fibre_sa_inv( sa_inv );
        String< TPos, Journal< JournalConfig< TPos > > > fibre_lcp( indexLCP( index ) );

        String< TPos > indices; //Stores the positions of indices that need to be inserted

        String< size_t > dels; //Stores the positions of indices that need to be removed                

        typename Iterator< String< Node > >::Type it_nodes = string.getjournal().get_first_node();
        typename Iterator< String< Node > >::Type it_nodes_end = string.getjournal().get_dummy_node();
        
        unsigned int k = 0;
        unsigned int suf_idx = 0;
        int current_shift = 0;
        unsigned int verify_length = 0;
        
        do{
#ifndef NDEBUG_SYNC
            std::cout << "Processing Node:" << std::endl;
            it_nodes->print_info();
#endif            
            if( !it_nodes->op()->insertion() && !it_nodes->op()->deletion() ){
                for( unsigned int j = it_nodes->position; j < it_nodes->position + it_nodes->length; ++j ){
                    fibre_sa[ sa_inv[ j - current_shift ] ] += current_shift; //TODO: modify to save rather than apply
                }
#ifndef NDEBUG_SYNC
                std::cout << "Processing block of length: " << it_nodes->length << std::endl;
#endif
                verify_length = it_nodes->length;
                j_goNext(it_nodes);
                continue; //nothing else to do
            }
            
            if( it_nodes->op()->deletion() ){
                for( unsigned int j = it_nodes->position; j < it_nodes->position - it_nodes->op()->by(); ++j ){
                    appendValue( dels, sa_inv[ j - current_shift ], Generous() );
                }
                erase( fibre_sa_inv, it_nodes->position, it_nodes->position - it_nodes->op()->by() );
            }
            
            if( (*it_nodes).op()->insertion() ){
                for( unsigned int j = it_nodes->position; j < it_nodes->position + it_nodes->length; ++j ){
                    appendValue( indices, j, Generous() ); //Adding inserted indices
                }
                insert( it_nodes->position, fibre_sa_inv, begin( sa_inv ), it_nodes->length ); //TODO: not safe, create insert stirng as (TValue)0
            }

            k = 0;
#ifndef NDEBUG_SYNC
            std::cout << "Verifying up to: " << verify_length << std::endl;
#endif
            while( ++k <= verify_length ){
            	suf_idx = value( sa_inv, it_nodes->position - k - current_shift );
                if( value( fibre_lcp, suf_idx ) <= k - 1 ){
                	if( suf_idx != 0 ){
                	    if( value( fibre_lcp, suf_idx - 1 ) <= k - 1 ){
                	    	break;
                	    }
                	}else{
                	    break;
                	}
                }
                
                appendValue( dels, suf_idx, Generous() );
                appendValue( indices, fibre_sa[suf_idx], Generous() );
            }
            current_shift += it_nodes->op()->by();
            j_goNext(it_nodes);
        }while( it_nodes != it_nodes_end );
        
        std::sort( begin( dels ), end( dels ) );
        dels = prefix( dels, std::unique( begin( dels ), end( dels ) ) );
        suffix_compare_functor< TString > cmp( string );
        std::sort( begin( indices ), end( indices ), cmp );
        indices = prefix( indices, std::unique( begin( indices ), end( indices ) ) );
        
        String< TPos > index_lcp;
        resize( index_lcp, length( indices ) );
        if( length( index_lcp ) != 0 ){
            for( unsigned int i = 0; i < length( index_lcp ) - 1; ++i ){
                index_lcp[i] = lcpLength( suffix( string, indices[i]), suffix( string, indices[i + 1] ) );
            }
            index_lcp[ length( index_lcp ) - 1 ] = 0;
        }
#ifndef NDEBUG_SYNC        
        std::cout << "Status( shifted ):" << std::endl;
        std::cout << "SA+LCP:" << std::endl;
        for( unsigned int i = 0; i < length( fibre_sa ); ++i ){
            std::cout << i << "\t: " << fibre_sa[i] << "\t[ " << fibre_lcp[i] << " ] " << suffix( string, fibre_sa[i] ) << std::endl;
        }
        std::cout << "Indices:" << std::endl;
        for( unsigned int i = 0; i < length( indices ); ++i ){
            std::cout << i << "\t: " << indices[i] << "\t[ " << index_lcp[i] << " ] " << suffix( string, indices[i] ) << std::endl;
        }
        std::cout << "Deletions:" << std::endl;
        for( unsigned int i = 0; i < length( dels ); ++i ){
            std::cout << dels[i] << "\t: " << fibre_sa[dels[i]] << "\t[ " << fibre_lcp[dels[i]] << " ] " << suffix( string, fibre_sa[dels[i]] ) << std::endl;
        }
#endif        
        for( size_t i = 0; i < length( dels ); ++i ){
            if( dels[i] - i != 0 ){
                fibre_lcp[ dels[i] - i - 1 ] = _min( fibre_lcp[ dels[i] - i - 1 ], fibre_lcp[ dels[i] - i ] );
            }
            erase( fibre_lcp, dels[i] - i );
            erase( fibre_sa, dels[i] - i );
        }
        
#ifndef NDEBUG_SYNC
        std::cout << "Status( shifted + deleted ):" << std::endl;
        std::cout << "SA+LCP:" << std::endl;
        for( unsigned int i = 0; i < length( fibre_sa ); ++i ){
            std::cout << i << "\t: " << fibre_sa[i] << "\t[ " << fibre_lcp[i] << " ] " << suffix( string, fibre_sa[i] ) << std::endl;
        }
#endif        
        typename Iterator< String< typename SAValue< TIndex >::Type, Journal< JournalConfig< typename SAValue< TIndex >::Type > > > >::Type it_sa = begin( fibre_sa );
        typename Iterator< String< typename SAValue< TIndex >::Type, Journal< JournalConfig< typename SAValue< TIndex >::Type > > > >::Type it_sa_end = end( fibre_sa );
        typename Iterator< String< TPos > >::Type it_index = begin( indices );
        typename Iterator< String< TPos > >::Type it_index_end = end( indices );
        typename Iterator< String< TPos > >::Type it_index_lcp = begin( index_lcp );
                
        unsigned int insert_pos = 0;
        unsigned int blocklength = 1;
        
        while( it_sa < it_sa_end && it_index < it_index_end ){
#ifndef NDEBUG_SYNC
            std::cout << "Processing suffix " << *it_sa << " and " << *it_index << std::endl;
#endif
            if( suffix( string, *it_sa ) > suffix( string, *it_index ) ){
                std::cout << *it_sa << " > " << *it_index << std::endl;
                blocklength = 1;

                while( ++it_index < end( indices ) && suffix( string, *it_index ) < suffix( string, *it_sa )){
                    ++blocklength;
                }
                
                std::cout << "Inserting Block of Length: " << blocklength << " at Position " << insert_pos << std::endl;
                insert( insert_pos, fibre_sa, it_index - blocklength, blocklength );
                insert( insert_pos, fibre_lcp, it_index_lcp, blocklength );
                
                if( insert_pos + blocklength != length( fibre_lcp ) ){
                    fibre_lcp[ insert_pos + blocklength - 1 ] = lcpLength( suffix( string, fibre_sa[ insert_pos + blocklength - 1 ]), suffix( string, fibre_sa[ insert_pos + blocklength ] ) );
                }
                if( insert_pos != 0 ){
                    fibre_lcp[ insert_pos - 1 ] = lcpLength( suffix( string, fibre_sa[ insert_pos - 1 ] ), suffix( string, fibre_sa[ insert_pos ] ) );
                }
                
                for( unsigned int k = 0; k < blocklength; ++k ){
                    fibre_sa_inv[ fibre_sa[ insert_pos + k ] ] = insert_pos + k;
                }
                
                insert_pos += blocklength;
                
                it_sa = begin( fibre_sa ) + insert_pos;
                it_sa_end = end( fibre_sa );

                if( !(it_index < it_index_end) ){
                    std::cout << "Breaking!" << std::endl;
                    break;
                }
            }
            ++it_sa;
            ++insert_pos;
        }
        
        if( it_index < it_index_end ){
            unsigned int pos = length( fibre_sa );
#ifndef NDEBUG_SYNC
            std::cout << "Appending remaining indices" << std::endl;
#endif
            append( fibre_sa, suffix( indices, it_index ) );
            append( fibre_lcp, suffix( index_lcp, it_index_lcp ) );
            
            for( ; it_index < it_index_end; ++it_index, ++pos ){
                fibre_sa_inv[ *it_index ] = pos;
            }
        }

        current_shift = 0;
        it_nodes = fibre_sa.getjournal().get_first_node();
        it_nodes_end = fibre_sa.getjournal().get_dummy_node();
#ifndef NDEBUG_SYNC
        fibre_sa.getjournal().print_nodes_inorder();
#endif
        do{
            if( !it_nodes->op()->insertion() && !it_nodes->op()->deletion() ){
#ifndef NDEBUG_SYNC
                std::cout << "Shifting by: " << current_shift << std::endl;
#endif
                for( unsigned int j = it_nodes->position; j < it_nodes->position + it_nodes->length; ++j ){
#ifndef NDEBUG_SYNC
                    std::cout << fibre_sa[ j ] << "th Suffix: ";
                    std::cout << fibre_sa_inv[ fibre_sa[ j ] ] << " -> ";
#endif
                    fibre_sa_inv[ fibre_sa[ j ] ] += current_shift; //TODO: modify to save rather than apply
#ifndef NDEBUG_SYNC
                    std::cout << fibre_sa_inv[ fibre_sa[ j ] ] << std::endl;
#endif
                }
                j_goNext(it_nodes);
                continue; //nothing else to do
            }
            
            current_shift += it_nodes->op()->by();
            j_goNext(it_nodes);
        }while( it_nodes != it_nodes_end );

        flatten( fibre_lcp );
        flatten( fibre_sa );
#ifndef NDEBUG_SYNC        
        std::cout << "Status( final ):" << std::endl;
        std::cout << "SA+LCP:" << std::endl;
        for( unsigned int i = 0; i < length( fibre_sa ); ++i ){
            std::cout << i << "\t: " << fibre_sa[i] << "\t[ " << fibre_lcp[i] << " ] " << suffix( string, fibre_sa[i] ) << "\t| " << fibre_sa_inv[i] <<  std::endl;
        }
#endif
    }

//#undef NDEBUG_SYNC
//#define NDEBUG_SYNC    
    template< typename TIndex, typename TValue, typename TSloppySpec, typename TSAInv >
    inline void synchronize_index( TIndex & index, String< TValue, Journal< JournalConfig< TValue, TSloppySpec, True > > > & string, TSAInv & sa_inv ){
    
        typedef String< TValue, Journal< JournalConfig< TValue, TSloppySpec, False > > > TString;
               
        typedef typename Position< TString >::Type TPos;
        typedef String< typename SAValue< TIndex >::Type > TSA;
        typedef String< TPos > TLCP;

        String< typename SAValue< TIndex >::Type, Journal< JournalConfig< typename SAValue< TIndex >::Type > > > fibre_sa( indexSA( index ) );
        String< typename SAValue< TIndex >::Type, Journal< JournalConfig< typename SAValue< TIndex >::Type > > > fibre_sa_inv( sa_inv );
        String< TPos, Journal< JournalConfig< TPos > > > fibre_lcp( indexLCP( index ) );

        String< TPos > indices; //Stores the positions of indices that need to be inserted

        String< size_t > dels; //Stores the positions of indices that need to be removed                

        typename Iterator< String< Node > >::Type it_nodes = string.getjournal().get_first_node();
        typename Iterator< String< Node > >::Type it_nodes_end = string.getjournal().get_dummy_node();
        
        unsigned int k = 0;
        unsigned int suf_idx = 0;
        int current_shift = 0;
        unsigned int verify_length = 0;
        
        do{
#ifndef NDEBUG_SYNC
            std::cout << "Processing Node:" << std::endl;
            it_nodes->print_info();
#endif            
            if( !it_nodes->op()->insertion() && !it_nodes->op()->deletion() ){
                for( unsigned int j = it_nodes->position; j < it_nodes->position + it_nodes->length; ++j ){
                    fibre_sa[ sa_inv[ j - current_shift ] ] += current_shift; //TODO: modify to save rather than apply
                }
#ifndef NDEBUG_SYNC
                std::cout << "Processing block of length: " << it_nodes->length << std::endl;
#endif
                verify_length = it_nodes->length;
                j_goNext(it_nodes);
                continue; //nothing else to do
            }
            
            if( it_nodes->op()->deletion() ){
                for( unsigned int j = it_nodes->position; j < it_nodes->position - it_nodes->op()->by(); ++j ){
                    appendValue( dels, sa_inv[ j - current_shift ], Generous() );
                }
                erase( fibre_sa_inv, it_nodes->position, it_nodes->position - it_nodes->op()->by() );
            }
            
            if( (*it_nodes).op()->insertion() ){
                for( unsigned int j = it_nodes->position; j < it_nodes->position + it_nodes->length; ++j ){
                    appendValue( indices, j, Generous() ); //Adding inserted indices
                }
                insert( it_nodes->position, fibre_sa_inv, begin( sa_inv ), it_nodes->length ); //TODO: not safe, create insert stirng as (TValue)0
            }

            k = 0;
#ifndef NDEBUG_SYNC
            std::cout << "Verifying up to: " << verify_length << std::endl;
#endif
            while( ++k <= verify_length ){
            	suf_idx = value( sa_inv, it_nodes->position - k - current_shift );
                if( value( fibre_lcp, suf_idx ) <= k - 1 ){
                	if( suf_idx != 0 ){
                	    if( value( fibre_lcp, suf_idx - 1 ) <= k - 1 ){
                	    	break;
                	    }
                	}else{
                	    break;
                	}
                }
                
                appendValue( dels, suf_idx, Generous() );
                appendValue( indices, fibre_sa[suf_idx], Generous() );
            }
            current_shift += it_nodes->op()->by();
            j_goNext(it_nodes);
        }while( it_nodes != it_nodes_end );
        
        std::sort( begin( dels ), end( dels ) );
        dels = prefix( dels, std::unique( begin( dels ), end( dels ) ) );
        suffix_compare_functor< TString > cmp( string );
        std::sort( begin( indices ), end( indices ), cmp );
        indices = prefix( indices, std::unique( begin( indices ), end( indices ) ) );
        
        String< TPos > index_lcp;
        resize( index_lcp, length( indices ) );
        if( length( index_lcp ) != 0 ){
            for( unsigned int i = 0; i < length( index_lcp ) - 1; ++i ){
                index_lcp[i] = lcpLength( suffix( string, indices[i]), suffix( string, indices[i + 1] ) );
            }
            index_lcp[ length( index_lcp ) - 1 ] = 0;
        }
#ifndef NDEBUG_SYNC        
        std::cout << "Status( shifted ):" << std::endl;
        std::cout << "SA+LCP:" << std::endl;
        for( unsigned int i = 0; i < length( fibre_sa ); ++i ){
            std::cout << i << "\t: " << fibre_sa[i] << "\t[ " << fibre_lcp[i] << " ] " << suffix( string, fibre_sa[i] ) << std::endl;
        }
        std::cout << "Indices:" << std::endl;
        for( unsigned int i = 0; i < length( indices ); ++i ){
            std::cout << i << "\t: " << indices[i] << "\t[ " << index_lcp[i] << " ] " << suffix( string, indices[i] ) << std::endl;
        }
        std::cout << "Deletions:" << std::endl;
        for( unsigned int i = 0; i < length( dels ); ++i ){
            std::cout << dels[i] << "\t: " << fibre_sa[dels[i]] << "\t[ " << fibre_lcp[dels[i]] << " ] " << suffix( string, fibre_sa[dels[i]] ) << std::endl;
        }
#endif        
        for( size_t i = 0; i < length( dels ); ++i ){
            if( dels[i] - i != 0 ){
                fibre_lcp[ dels[i] - i - 1 ] = _min( fibre_lcp[ dels[i] - i - 1 ], fibre_lcp[ dels[i] - i ] );
            }
            erase( fibre_lcp, dels[i] - i );
            erase( fibre_sa, dels[i] - i );
        }
        
#ifndef NDEBUG_SYNC
        std::cout << "Status( shifted + deleted ):" << std::endl;
        std::cout << "SA+LCP:" << std::endl;
        for( unsigned int i = 0; i < length( fibre_sa ); ++i ){
            std::cout << i << "\t: " << fibre_sa[i] << "\t[ " << fibre_lcp[i] << " ] " << suffix( string, fibre_sa[i] ) << std::endl;
        }
#endif        
        typename Iterator< String< typename SAValue< TIndex >::Type, Journal< JournalConfig< typename SAValue< TIndex >::Type > > > >::Type it_sa = begin( fibre_sa );
        typename Iterator< String< typename SAValue< TIndex >::Type, Journal< JournalConfig< typename SAValue< TIndex >::Type > > > >::Type it_sa_end = end( fibre_sa );
        typename Iterator< String< TPos > >::Type it_index = begin( indices );
        typename Iterator< String< TPos > >::Type it_index_end = end( indices );
        typename Iterator< String< TPos > >::Type it_index_lcp = begin( index_lcp );
                
        unsigned int insert_pos = 0;
        unsigned int blocklength = 1;
        
        while( it_sa < it_sa_end && it_index < it_index_end ){
#ifndef NDEBUG_SYNC
            std::cout << "Processing suffix " << *it_sa << " and " << *it_index << std::endl;
#endif
            if( suffix( string, *it_sa ) > suffix( string, *it_index ) ){
                std::cout << *it_sa << " > " << *it_index << std::endl;
                blocklength = 1;

                while( ++it_index < end( indices ) && suffix( string, *it_index ) < suffix( string, *it_sa )){
                    ++blocklength;
                }
                
                std::cout << "Inserting Block of Length: " << blocklength << " at Position " << insert_pos << std::endl;
                insert( insert_pos, fibre_sa, it_index - blocklength, blocklength );
                insert( insert_pos, fibre_lcp, it_index_lcp, blocklength );
                
                if( insert_pos + blocklength != length( fibre_lcp ) ){
                    fibre_lcp[ insert_pos + blocklength - 1 ] = lcpLength( suffix( string, fibre_sa[ insert_pos + blocklength - 1 ]), suffix( string, fibre_sa[ insert_pos + blocklength ] ) );
                }
                if( insert_pos != 0 ){
                    fibre_lcp[ insert_pos - 1 ] = lcpLength( suffix( string, fibre_sa[ insert_pos - 1 ] ), suffix( string, fibre_sa[ insert_pos ] ) );
                }
                
                for( unsigned int k = 0; k < blocklength; ++k ){
                    fibre_sa_inv[ fibre_sa[ insert_pos + k ] ] = insert_pos + k;
                }
                
                insert_pos += blocklength;
                
                it_sa = begin( fibre_sa ) + insert_pos;
                it_sa_end = end( fibre_sa );

                if( !(it_index < it_index_end) ){
                    std::cout << "Breaking!" << std::endl;
                    break;
                }
            }
            ++it_sa;
            ++insert_pos;
        }
        
        if( it_index < it_index_end ){
            unsigned int pos = length( fibre_sa );
#ifndef NDEBUG_SYNC
            std::cout << "Appending remaining indices" << std::endl;
#endif
            append( fibre_sa, suffix( indices, it_index ) );
            append( fibre_lcp, suffix( index_lcp, it_index_lcp ) );
            
            for( ; it_index < it_index_end; ++it_index, ++pos ){
                fibre_sa_inv[ *it_index ] = pos;
            }
        }

        current_shift = 0;
        it_nodes = fibre_sa.getjournal().get_first_node();
        it_nodes_end = fibre_sa.getjournal().get_dummy_node();
#ifndef NDEBUG_SYNC
        fibre_sa.getjournal().print_nodes_inorder();
#endif
        do{
            if( !it_nodes->op()->insertion() && !it_nodes->op()->deletion() ){
#ifndef NDEBUG_SYNC
                std::cout << "Shifting by: " << current_shift << std::endl;
#endif
                for( unsigned int j = it_nodes->position; j < it_nodes->position + it_nodes->length; ++j ){
#ifndef NDEBUG_SYNC
                    std::cout << fibre_sa[ j ] << "th Suffix: ";
                    std::cout << fibre_sa_inv[ fibre_sa[ j ] ] << " -> ";
#endif
                    fibre_sa_inv[ fibre_sa[ j ] ] += current_shift; //TODO: modify to save rather than apply
#ifndef NDEBUG_SYNC
                    std::cout << fibre_sa_inv[ fibre_sa[ j ] ] << std::endl;
#endif
                }
                j_goNext(it_nodes);
                continue; //nothing else to do
            }
            
            current_shift += it_nodes->op()->by();
            j_goNext(it_nodes);
        }while( it_nodes != it_nodes_end );

        flatten( fibre_lcp );
        flatten( fibre_sa );
#ifndef NDEBUG_SYNC        
        std::cout << "Status( final ):" << std::endl;
        std::cout << "SA+LCP:" << std::endl;
        for( unsigned int i = 0; i < length( fibre_sa ); ++i ){
            std::cout << i << "\t: " << fibre_sa[i] << "\t[ " << fibre_lcp[i] << " ] " << suffix( string, fibre_sa[i] ) << "\t| " << fibre_sa_inv[i] <<  std::endl;
        }
#endif
    }

//#undef NDEBUG_SYNC
}
