// ==========================================================================
//                                  FMIndex
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================

#ifndef FM_INDEX_LF_TABLE_H_
#define FM_INDEX_LF_TABLE_H_

namespace seqan {

// ==========================================================================
// LFTable is an object storing all necessary information for the LF-mapping.
// To be more precise, the occurrence-table data structure as well as the
// prefix-sum table are stored.
// ==========================================================================

template <typename TOccTable, typename TPrefixSumTable>
struct LFTable
{
    TOccTable occTable;
    TPrefixSumTable prefixSumTable;

    LFTable() :
        occTable(),
        prefixSumTable()
    {}

    LFTable(TOccTable const & occTable, TPrefixSumTable const & prefixSumTable) :
        occTable(occTable),
        prefixSumTable(prefixSumTable)
    {}

    inline LFTable & operator=(LFTable const & other)
    {
        occTable = other.occTable;
        prefixSumTable = other.prefixSumTable;
        return *this;
    }

    inline bool operator==(const LFTable & b) const
    {
        return occTable == b.occTable &&
               prefixSumTable == b.prefixSumTable;
    }

};

struct FibreOccTable_;
struct FibrePrefixSumTable_;

typedef Tag<FibreOccTable_> const FibreOccTable;
typedef Tag<FibrePrefixSumTable_> const FMTablePrefixSumTable;

template <typename TOccTable, typename TPrefixSumTable>
struct Fibre<LFTable<TOccTable, TPrefixSumTable>, FibreOccTable>
{
    typedef TOccTable Type;
};

template <typename TOccTable, typename TPrefixSumTable>
struct Fibre<LFTable<TOccTable, TPrefixSumTable> const, FibreOccTable>
{
    typedef TOccTable const Type;
};

template <typename TOccTable, typename TPrefixSumTable>
struct Fibre<LFTable<TOccTable, TPrefixSumTable>, FMTablePrefixSumTable>
{
    typedef TPrefixSumTable Type;
};

template <typename TOccTable, typename TPrefixSumTable>
struct Reference<LFTable<TOccTable, TPrefixSumTable> >
{
    typedef TPrefixSumTable & Type;
};

template <typename TOccTable, typename TPrefixSumTable>
inline typename Fibre<LFTable<TOccTable, TPrefixSumTable>, FMTablePrefixSumTable>::Type const &
getFibre(LFTable<TOccTable, TPrefixSumTable> const & lfTable, FMTablePrefixSumTable)
{
    return lfTable.prefixSumTable;
}

template <typename TOccTable, typename TPrefixSumTable>
inline typename Fibre<LFTable<TOccTable, TPrefixSumTable>, FMTablePrefixSumTable>::Type &
getFibre(LFTable<TOccTable, TPrefixSumTable>&lfTable, FMTablePrefixSumTable)
{
    return lfTable.prefixSumTable;
}

template <typename TOccTable, typename TPrefixSumTable>
inline typename Fibre<LFTable<TOccTable, TPrefixSumTable>, FibreOccTable>::Type &
getFibre(LFTable<TOccTable, TPrefixSumTable>&lfTable, FibreOccTable)
{
    return lfTable.occTable;
}

template <typename TOccTable, typename TPrefixSumTable>
inline typename Fibre<LFTable<TOccTable, TPrefixSumTable>, FibreOccTable>::Type const &
getFibre(LFTable<TOccTable, TPrefixSumTable> const & lfTable, FibreOccTable)
{
    return lfTable.occTable;
}

template <typename TOccTable, typename TPrefixSumTable>
inline void clear(LFTable<TOccTable, TPrefixSumTable> & lfTable)
{
    clear(lfTable.occTable);
    clear(lfTable.prefixSumTable);
}

template <typename TOccTable, typename TPrefixSumTable>
inline bool open(
    LFTable<TOccTable, TPrefixSumTable> & lfTable,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName;    append(name, ".occ");
    if (!open(getFibre(lfTable, FibreOccTable()), toCString(name), openMode))
    {
        return false;
    }
    name = fileName;    append(name, ".psum");  open(getFibre(lfTable, FMTablePrefixSumTable()), toCString(name), openMode);
    return true;

}

template <typename TOccTable, typename TPrefixSumTable>
inline bool open(
    LFTable<TOccTable, TPrefixSumTable> & lfTable,
    const char * fileName)
{
    return open(lfTable, fileName, DefaultOpenMode<LFTable<TOccTable, TPrefixSumTable> >::VALUE);
}

template <typename TOccTable, typename TPrefixSumTable>
inline bool save(
    LFTable<TOccTable, TPrefixSumTable> const & lfTable,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName;    append(name, ".occ");
    if (!save(getFibre(lfTable, FibreOccTable()), toCString(name), openMode))
    {
        return false;
    }
    name = fileName;    append(name, ".psum");  save(getFibre(lfTable, FMTablePrefixSumTable()), toCString(name), openMode);
    return true;
}

template <typename TOccTable, typename TPrefixSumTable>
inline bool save(
    LFTable<TOccTable, TPrefixSumTable> const & lfTable,
    const char * fileName)
{
    return save(lfTable, fileName, DefaultOpenMode<LFTable<TOccTable, TPrefixSumTable> >::VALUE);
}

}
#endif /* FM_INDEX_LF_TABLE_H_ */