#if !defined(_MICROBENCH_H)
#define _MICROBENCH_H

#include <ctype.h>

static inline uint32_t read_time(void)
{
  uint32_t a;
  __asm__ volatile( "rdtsc" :"=a"(a) ::"edx" );
  return a;
}

#endif
