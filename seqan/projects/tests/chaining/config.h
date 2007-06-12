#ifndef CONFIG_INCLUDED
#define CONFIG_INCLUDED

#pragma warning(disable:4786 4099 4355 985)

#define GCC_VERSION ( __GNUC__ * 10000 + \
                      __GNUC_MINOR__ * 100 + \
                      __GNUC_PATCHLEVEL__ )
#ifndef _WIN32
#if GCC_VERSION < 29700		/* otherwise STL is used? [tst] */
  #include "common/minmax.h"
#endif
#endif


#include <vector>

#ifdef HAVE_RESTRICT
#define RESTRICT restrict
#else
#define RESTRICT
#endif

template<class T>
void delete_vector(std::vector<T>& v) {
    typename std::vector<T>::const_iterator it;
    for(it=v.begin();it!=v.end();++it) delete *it;
    v.resize(0);
}

#endif

// 64-Bit or 32-Bit architecture

#ifdef _64
#define Long int
#else
#define Long long
#endif

