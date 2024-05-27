#ifndef PACE2024EXACT_MACROS_H
#define PACE2024EXACT_MACROS_H

namespace CrossGuard {

#define ASSERT_ENABLED 0

#if ASSERT_ENABLED
#define ASSERT(condition) if(!(condition)) {std::cerr << "Error in file " << __FILE__ << " in function " << __FUNCTION__ << " at line " << __LINE__ << "!" << std::endl; abort(); } ((void)0)
#else
#define ASSERT(condition)  ((void)0)
#endif

}

#endif //PACE2024EXACT_MACROS_H
