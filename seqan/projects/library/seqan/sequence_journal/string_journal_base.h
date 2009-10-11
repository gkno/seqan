#ifndef SEQAN_STRING_JOURNAL_BASE_H
#define SEQAN_STRING_JOURNAL_BASE_H

#define NDEBUG

#ifndef NDEBUG
#  define DEBUG_OUT(args) \
   std::cerr << args << std::endl;
#else
#  define DEBUG_OUT(args)
#endif

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>

#include "string_journal_forwards.h"
#include "string_journal_operation.h"
#include "string_journal_node.h"
#include "journal.h"
#include "string_journal.h"
#include "iterator_journal.h"
#include "string_journal_interface.h"
#include "string_journal_debug.h"
#include "string_journal_utility.h"
#include "string_journal_test_foundry.h"

namespace seqan {
    template< typename TValue, typename TSpec >
    class ShiftString{
    typedef size_t TPos; //TODO: use Position metafunction
    public:
        ShiftString( String< TValue, TSpec > & string, String< Pair< TPos, TValue > > & shifts ) : m_holder( string ), m_shifts( shifts ) {};
        
        inline TValue value( TPos & pos ){
            return value( value( m_holder ), pos ) + get_shift( m_shifts, pos );
        }
    private:
        String< TPos, TValue > m_shifts;
        Holder< String< TValue, TSpec > > m_holder;
    };
    
    template< typename TValue, typename TSpec, typename TPos >
    TValue value( ShiftString< TValue, TSpec > & sstring, TPos & position ){
        return sstring.value( position );
    }
}
#endif // ndef(SEQAN_STRING_JOURNAL_BASE_H)
