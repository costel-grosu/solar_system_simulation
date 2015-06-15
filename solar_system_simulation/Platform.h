
#if !defined( __Platform_h__ )
#define __Platform_h__

#include "stdafx.h"
#include <stdio.h>
#include <assert.h>
#include <cmath>

#define now __rdtsc

#define FLOAT double
#define V3 dvec3

#ifdef __GNUC__
#define INLINE __attribute__((always_inline))
#define NOINLINE 
#elif defined _MSC_VER
#define NOINLINE __declspec(noinline)
#define INLINE __forceinline
#else
// unsupported compiler
#endif

#endif // __Platform_h__
