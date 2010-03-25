/*==========================================================================
 RazerS - Fast Read Mapping with Controlled Loss Rate
 http://www.seqan.de/projects/razers.html
 
 ============================================================================
 Copyright (C) 2010
 
 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 3 of the License, or (at your options) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 Lesser General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ==========================================================================*/

#ifndef SEQAN_HEADER_RAZERS_PARALLEL_PIPELINE_H
#define SEQAN_HEADER_RAZERS_PARALLEL_PIPELINE_H

#include <iostream>

namespace SEQAN_NAMESPACE_MAIN
{
	// Filter that writes each buffer to a file.
	class FilterFilter: public tbb::filter {
		public:
			FilterFilter();
			/*override*/void* operator()( void* item );
	};
	
	FilterFilter::FilterFilter(): tbb::filter(serial_in_order) {}
	
	void* MyOutputFilter::operator()( void* item ) {
		return NULL;
	}
	
	class MyInputFilter: public tbb::filter {
		public:
			static const size_t n_buffer = 8;
			MyInputFilter( FILE* input_file_ );
		private:
			FILE* input_file;
			size_t next_buffer;
			char last_char_of_previous_buffer;
			MyBuffer buffer[n_buffer];
			/*override*/ void* operator()(void*);
	};
	
	MyInputFilter::MyInputFilter( FILE* input_file_ ) :
		filter(serial_in_order),
		next_buffer(0),
		input_file(input_file_),
		last_char_of_previous_buffer(' ')
	{}
	
	void* MyInputFilter::operator()(void*) {
		MyBuffer& b = buffer[next_buffer];
		next_buffer = (next_buffer+1) % n_buffer;
		size_t n = fread( b.begin(), 1, b.max_size(), input_file );
		if( !n ) {
			// end of file
			return NULL;
		} else {
			b.begin()[-1] = last_char_of_previous_buffer;
			last_char_of_previous_buffer = b.begin()[n-1];
			b.set_end( b.begin()+n );
			return &b;
		}
	}
	
}

#endif