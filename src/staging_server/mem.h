#ifndef __MEM_H__
#define __MEM_H__

#ifdef __MEM_C__
#define EXT
#else
#define EXT extern
#endif

// Return the available memory but max the input argument if input is > 0
uint64_t mem_get_available (uint64_t max_allowed);

#undef EXT
#endif
