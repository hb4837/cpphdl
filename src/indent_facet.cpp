#include "indent_facet.hpp"

indent_facet::result indent_facet::do_out(state_type &need_indentation,
        const intern_type *from, const intern_type *from_end, const intern_type *&from_next,
        extern_type *to, extern_type *to_end, extern_type *&to_next
        ) const {
    result res = std::codecvt_base::noconv;
    for (; (from < from_end) && (to < to_end); ++from, ++to) {
        // 0 indicates that the last character seen was a newline.
        // thus we will print a tab before it. Ignore it the next
        // character is also a newline
        if ((state(need_indentation) == 0) && (*from != '\n')) {
            res = std::codecvt_base::ok;
            state(need_indentation) = 1;
            for (int i = 0; i < m_indentation_level; ++i) {
                *to = ' ';
                ++to;
            }
            if (to == to_end) {
                res = std::codecvt_base::partial;
                break;
            }
        }
        *to = *from; // Copy the next character.

        // If the character copied was a '\n' mark that state
        if (*from == '\n') {
            state(need_indentation) = 0;
        }
    }

    if (from != from_end) {
        res = std::codecvt_base::partial;
    }
    from_next = from;
    to_next = to;

    return res;
};

// static const int index = std::ios_base::xalloc();

namespace im {

    static bool initialized = false;
    
    std::ostream& push(std::ostream& os) {
        if (!initialized) {
            std::ios_base::sync_with_stdio(false);
            im::raii_guard rg(std::cout);
            initialized = true;
        }

        auto ilevel = ++os.iword(index);
        os.imbue(std::locale(os.getloc(), new indent_facet(ilevel)));
        return os;
    }

    std::ostream& pop(std::ostream& os) {
        auto ilevel = (os.iword(index) > 0) ? --os.iword(index) : 0;
        os.imbue(std::locale(os.getloc(), new indent_facet(ilevel)));
        return os;
    }

    /// Clears the ostream indentation set, but NOT the raii_guard.

    std::ostream& clear(std::ostream& os) {
        os.iword(index) = 0;
        os.imbue(std::locale(os.getloc(), new indent_facet(0)));
        return os;
    }

} /* namespace indent_manip */