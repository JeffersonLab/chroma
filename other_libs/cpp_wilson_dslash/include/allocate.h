/*
 * allocate.h
 *
 *  Created on: Aug 26, 2016
 *      Author: bjoo
 */

#ifndef CPP_WILSON_DSLASH_INCLUDE_ALLOCATE_H_
#define CPP_WILSON_DSLASH_INCLUDE_ALLOCATE_H_

#include "dslash_config.h"

#ifdef SSE_USE_QDPXX
#include "qdp.h"
#include "qdp_allocator.h"
#else
#include <cstdlib>
#endif

namespace  CPlusPlusWilsonDslash {

inline
void *alloc(size_t numbytes) {
#ifdef SSE_USE_QDPXX
	void* ret_val = (void *)QDP::Allocator::theQDPAllocator::Instance().allocate(numbytes, QDP::Allocator::DEFAULT);
	if( ret_val == nullptr ) {
	   QDP::QDPIO::cout << "cpp_wilson_dslash: alloc failed. Dumping map and exiting" << std::endl;
	   QDP::Allocator::theQDPAllocator::Instance().dump();
	   QDP::QDP_abort(1);
        }
	return ret_val;
#else
	return (void *)std::malloc(numbytes);
#endif
}

inline
void dealloc(void *mem) {
#ifdef SSE_USE_QDPXX
	QDP::Allocator::theQDPAllocator::Instance().free(mem);
#else
	std::free(mem);
#endif
}

}


#endif /* OTHER_LIBS_CPP_WILSON_DSLASH_INCLUDE_ALLOCATE_H_ */
