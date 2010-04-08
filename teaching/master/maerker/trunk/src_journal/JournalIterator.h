/*
 * JournalIterator.h
 *
 *  Created on: 08.04.2010
 *      Author: Rene-User
 */

#include "JournalString.h"

#ifndef JOURNALITERATOR_H_
#define JOURNALITERATOR_H_

namespace seqan {

	template< typename TValue, typename TStringSpec, typename TJournalSpec >
	struct JournalIterator{

		typedef typename seqan::Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type Type;

		JournalIterator( seqan::String< TValue, Journal< TStringSpec, TJournalSpec > > & str )
			: m_string( str ),
			m_iter_underlying( seqan::begin( underlying( str ) ) ),
			m_iter_insertion( seqan::begin( insertion_string( str ) ) ),
			m_iter_tree( seqan::begin( journal_tree( str ) ) )
		{
			while( isDeletion( cargo(*m_iter_tree).op ) && seqan::goNext( m_iter_tree ) ) { };
			m_remaining_blocksize = cargo(*m_iter_tree).blocksize - 1;
			m_iter_underlying += cargo(*m_iter_tree).physical_position;
			m_iter_insertion += cargo(*m_iter_tree).physical_position;
		}

 		JournalIterator( seqan::String< TValue, Journal< TStringSpec, TJournalSpec > > & str, seqan::TreeIterator< seqan::SeqTree< STNode< JournalNode > > > it_tree )
			: m_string( str ),
			m_iter_underlying( seqan::begin( str.underlying() ) + cargo(*it_tree).physical_position ),
			m_iter_insertion( seqan::begin( str.insertion_string() ) + cargo(*it_tree).physical_position ),
			m_iter_tree( it_tree )
		{
			m_remaining_blocksize = cargo(*m_iter_tree).blocksize - 1;
		}

		inline void operator= ( Type other ){
			m_string = other.string();
			m_iter_underlying = other.iter_underlying();
			m_iter_insertion = other.iter_insertion();
			m_iter_tree = other.iter_tree();
			m_remaining_blocksize = other.remaining_blocksize();
		}

		inline bool operator== ( Type & other ){
			return ( m_iter_underlying == other.iter_underlying() && m_iter_insertion == other.iter_insertion() && m_remaining_blocksize == other.remaining_blocksize() );
			//m_string == other.string() && && m_iter_tree == other.iter_tree() && m_remaining_blocksize == other.remaining_blocksize() );
		}

		inline bool operator!= ( Type & other ){
			return !( operator==(other) );
		}

		inline TValue operator*(){
			if( !isInsertion( cargo(*m_iter_tree).op ) ){
				return *m_iter_underlying;
			}else{
				return *m_iter_insertion;
			}
		}

		inline Type & operator++(){
			if( m_remaining_blocksize > 0 ){
				++m_iter_insertion;
				++m_iter_underlying;
				--m_remaining_blocksize;
			}else{
				while( seqan::goNext( m_iter_tree ) && isDeletion( cargo(*m_iter_tree).op ) ) { }; //TODO: specify behavior for end of string reached -> is reset to tree root regardless of root-node-type atm.
				m_iter_insertion = seqan::begin( m_string.insertion_string() ) + cargo(*m_iter_tree).physical_position;
				m_iter_underlying = seqan::begin( m_string.underlying() ) + cargo(*m_iter_tree).physical_position;
				m_remaining_blocksize = cargo(*m_iter_tree).blocksize - 1;
			}
			return *this;
		}

		inline Type operator++( int ){
			Type temp_it = *this;
			++(*this);
			return temp_it;
		}

		inline Type & operator--(){
			minus_one();
			return *this;
		}

		inline Type operator--( int ){
			Type temp_it = *this;
			minus_one();
			return temp_it;
		}

		int remaining_blocksize() {
			return m_remaining_blocksize;
		}

		typename seqan::Iterator< String< TValue, TStringSpec > >::Type iter_underlying(){
			return m_iter_underlying;
		}

		typename seqan::Iterator< String< TValue > >::Type iter_insertion(){
			return m_iter_insertion;
		}

		seqan::TreeIterator< seqan::SeqTree< STNode< JournalNode > > > iter_tree(){
			return m_iter_tree;
		}

		seqan::String< TValue, Journal< TStringSpec, TJournalSpec > > & string(){
			return m_string;
		}

	private:

		inline void minus_one(){
			if( m_remaining_blocksize < (int)cargo(*m_iter_tree).blocksize - 1 ){
				--m_iter_insertion;
				--m_iter_underlying;
				++m_remaining_blocksize;
			}else{
				while( seqan::goPrev( m_iter_tree ) && isDeletion( cargo(*m_iter_tree).op ) ) { }; //TODO: specify behavior for end of string reached -> is reset to tree root regardless of root-node-type atm.
				m_iter_insertion = seqan::begin( m_string.insertion_string() ) + cargo(*m_iter_tree).physical_position + cargo(*m_iter_tree).blocksize - 1;
				m_iter_underlying = seqan::begin( m_string.underlying() ) + cargo(*m_iter_tree).physical_position + cargo(*m_iter_tree).blocksize - 1;
				m_remaining_blocksize = 0;
			}
		}

		int m_remaining_blocksize;
		typename seqan::Iterator< String< TValue, TStringSpec > >::Type m_iter_underlying;
		typename seqan::Iterator< String< TValue > >::Type m_iter_insertion;
		seqan::TreeIterator< seqan::SeqTree< STNode< JournalNode > > > m_iter_tree;
		seqan::String< TValue, Journal< TStringSpec, TJournalSpec > > & m_string;
	};

	template< typename TValue, typename TStringSpec, typename TJournalSpec >
	struct Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > >{
		typedef JournalIterator< TValue, TStringSpec, TJournalSpec > Type;
	};

	template< typename TValue, typename TStringSpec, typename TJournalSpec >
	inline typename Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type begin( String< TValue, Journal< TStringSpec, TJournalSpec > > & str ){
		typename Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type begin_iter( str );
		return begin_iter;
	}

	template< typename TValue, typename TStringSpec, typename TJournalSpec >
	inline typename Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type end( String< TValue, Journal< TStringSpec, TJournalSpec > > & str ){
		typename Iterator< String< TValue, Journal< TStringSpec, TJournalSpec > > >::Type end_iter( str, seqan::end( journal_tree( str ) ) );
		return end_iter;
	}
};

#endif /* JOURNALITERATOR_H_ */
