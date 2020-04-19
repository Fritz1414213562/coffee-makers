#ifndef COFFEE_MAKERS_COMPILE_TIME_CONJUNCTION_HPP
#define COFFEE_MAKERS_COMPILE_TIME_CONJUNCTION_HPP

#include<type_traits>


namespace makers {

template<bool... Bs>
using value_conjunction = std::conjunction<std::bool_constant<Bs>...>;

}


#endif /* COFFEE_MAKERS_COMPILE_TIME_CONJUNCTION_HPP */
