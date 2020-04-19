#ifndef COFFEE_MAKERS_COMPILE_TIME_DISJUNCTION_HPP
#define COFFEE_MAKERS_COMPILE_TIME_DISJUNCTION_HPP

#include<type_traits>


namespace makers {

template<bool... Bs>
using value_disjunction = std::disjunction<std::bool_constant<Bs>...>;

}


#endif /* COFFEE_MAKERS_COMPILE_TIME_DISJUNCTION_HPP */
