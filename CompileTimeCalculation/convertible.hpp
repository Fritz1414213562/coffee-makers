#ifndef COFFEE_MAKERS_COMPILE_TIME_CONVERTIBLE_HPP
#define COFFEE_MAKERS_COMPILE_TIME_CONVERTIBLE_HPP
#include<type_traits>


namespace makers {

template<typename lhs, typename rhs>
using convertible = std::conjunction<std::is_convertible<lhs, rhs>, std::is_convertible<rhs, lhs>>;

}

#endif // COFFEE_MAKERS_COMPILE_TIME_CONVERTIBLE_HPP
