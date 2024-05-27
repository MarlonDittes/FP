#ifndef PACE2024EXACT_TRANSLATIONTABLE_H
#define PACE2024EXACT_TRANSLATIONTABLE_H

#include <cstddef>
#include <map>

#include "definitions.h"
#include "macros.h"
#include "misc.h"

namespace CrossGuard {

/**
 * When modifying the graph, we will rename vertices. This class stores the
 * translations, so they can be reverted.
 */
    class TranslationTable {
        std::map<u32, u32> translation_A_new_to_old;
        std::map<u32, u32> translation_A_old_to_new;
        std::map<u32, u32> translation_B_new_to_old;
        std::map<u32, u32> translation_B_old_to_new;

    public:
        explicit TranslationTable() = default;

        inline void add_A(u32 o, u32 n) {
            translation_A_new_to_old[n] = o;
            translation_A_old_to_new[o] = n;
        }

        inline void add_B(u32 o, u32 n) {
            translation_B_new_to_old[n] = o;
            translation_B_old_to_new[o] = n;
        }

        inline u32 get_A_old(u32 n) const {
            ASSERT(translation_A_new_to_old.contains(n));

            return translation_A_new_to_old.at(n);
        }

        inline u32 get_A_new(u32 o) const {
            ASSERT(translation_A_old_to_new.contains(o));

            return translation_A_old_to_new.at(o);
        }

        inline u32 get_B_old(u32 n) const {
            ASSERT(translation_B_new_to_old.contains(n));

            return translation_B_new_to_old.at(n);
        }

        inline u32 get_B_new(u32 o) const {
            ASSERT(translation_B_old_to_new.contains(o));

            return translation_B_old_to_new.at(o);
        }

    };

}

#endif //PACE2024EXACT_TRANSLATIONTABLE_H
