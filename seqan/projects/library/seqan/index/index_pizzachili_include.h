 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id$
 ==========================================================================*/

#include <cstdlib>
#include <algorithm>
#include <vector>

#ifndef SEQAN_HEADER_INDEX_PC_SCAFFOLD_H
#define SEQAN_HEADER_INDEX_PC_SCAFFOLD_H

namespace SEQAN_NAMESPACE_MAIN {
    //
    // Adapter interface for the Pizza & Chili compressed indices.
    //

    namespace impl { namespace pizzachili {
        typedef unsigned char uchar;
        typedef unsigned long ulong;
        typedef void* index_t;
        typedef int error_t;

#if 0

        extern "C" {
#   include "pizzachili_interface.h"

            // Undo header file defines:
#   undef uchar
#   undef uint
#   undef ulong
        }

#else

        struct test_index {
            ulong m_len;
            uchar* m_text;

            test_index(uchar* text, ulong length)
                : m_len(length), m_text(new uchar[m_len])
            {
                ::std::copy(text, text + length, m_text);
            }

            ~test_index() {
                delete [] m_text;
            }

        private:
            test_index(test_index const&);
            test_index& operator=(test_index const&);
        };

        char* error_index(int e) {
            return "No error.";
        }

        //
        // Building the index
        //

        int build_index(
            uchar* text,
            ulong length,
            char* /* build_options */,
            void** index
        ) {
            test_index* idx = static_cast<test_index*>(::std::malloc(sizeof(test_index)));
            new (idx) test_index(text, length);
            *index = idx;

            return 0;
        }

        int save_index(void* index, char* filename) {
            return 0;
        }

        int load_index(char* filename, void** index) {
            return 0;
        }

        int free_index(void* index) {
            test_index* idx = static_cast<test_index*>(index);
            idx->~test_index();
            ::std::free(idx);
            return 0;
        }

        int index_size(void* index, ulong* size) {
            test_index* idx = static_cast<test_index*>(index);
            *size = sizeof(test_index) + sizeof(uchar) * idx->m_len;
            return 0;
        }

        //
        // Querying the index
        //

        int count(
            void* index,
            uchar* pattern,
            ulong length,
            ulong* numocc
        ) {
            *numocc = 2;
            return 0;
        }

        int locate(
            void* index,
            uchar* pattern,
            ulong length,
            ulong** occ,
            ulong* numocc
        ) {
            test_index* idx = static_cast<test_index*>(index);
            uchar* start = idx->m_text;
            uchar* end = start + idx->m_len;

            ::std::vector<ulong> result;
            for (
                uchar* i = ::std::search(start, end, pattern, pattern + length);
                i != end;
                i = ::std::search(i + 1, end, pattern, pattern + length)
            )
                result.push_back(i - start);

            *numocc = result.size();
            *occ = static_cast<ulong*>(::std::malloc(sizeof(ulong) * result.size()));
            ::std::copy(result.begin(), result.end(), *occ);
            return 0;
        }

        //
        // Accessing the indexed text.
        //

        int extract(
            void *index,
            ulong from,
            ulong to,
            uchar **snippet,
            ulong *snippet_length
        ) {
            test_index* idx = static_cast<test_index*>(index);
            *snippet_length = to - from + 1;
            *snippet = static_cast<uchar*>(std::malloc(sizeof(uchar) * *snippet_length));
            std::copy(idx->m_text + from, idx->m_text + to + 1, *snippet);
            return 0;
        }

        int display(
            void *index,
            uchar *pattern,
            ulong length,
            ulong numc,
            ulong *numocc,
            uchar **snippet_text,
            ulong **snippet_lengths
        ) {
            return 0;
        }

        int length(
            void *index,
            ulong *length
        ) {
            test_index* idx = static_cast<test_index*>(index);
            *length = idx->m_len;
            return 0;
        }

#endif

    } } // namespace impl::pizzachili
} // namespace SEQAN_NAMESPACE_MAIN

#endif // SEQAN_HEADER_INDEX_PC_SCAFFOLD_H
