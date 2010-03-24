/*
  Helpers for the test_modifier tests.
 */

#ifndef TESTS_MODIFIER_HELPERS_H_
#define TESTS_MODIFIER_HELPERS_H_

#include <functional>

// Cesar Chiffre as a test functor.
template <typename TArgChar, typename TResultChar = TArgChar>
struct CaesarChiffre : public std::unary_function<TArgChar, TResultChar> 
{
    TArgChar _delta;

    CaesarChiffre() : _delta(0) {}

    CaesarChiffre(TArgChar delta) {
        if (delta < 0)
            _delta = ('z' - 'a' + 1) - (-delta) % ('z' - 'a' + 1);
        else
            _delta = delta;
    }

    inline TResultChar operator()(TArgChar const & x) const {
        if (('a' <= x) && (x <= 'z'))
            return (x - 'a' + _delta) % ('z' - 'a' + 1) + 'a';
        if (('A' <= x) && (x <= 'Z'))
            return (x - 'A' + _delta) % ('Z' - 'A' + 1) + 'A';
        return x; 
    }
};


#endif  // TESTS_MODIFIER_HELPERS_H_
