#ifndef COFFEE_MAKERS_COMPILE_TIME_RETURN_TYPE_TRAITS_HPP
#define COFFEE_MAKERS_COMPILE_TIME_RETURN_TYPE_TRAITS_HPP
#include<coffee-makers/CompileTimeCalculation/conjunction.hpp>
#include<coffee-makers/CompileTimeCalculation/convertible.hpp>
#include<type_traits>


namespace makers {

template<typename l_realT, typename r_realT>
struct longer_prec {
	static_assert(makers::value_conjunction<
		std::is_floating_point<l_realT>::value,
		std::is_floating_point<r_realT>::value>::value,
	"Template paramters are not numerical type");
	
	using type = typename std::conditional<(sizeof(l_realT) > sizeof(r_realT)), l_realT, r_realT>::type;
};


template<typename realT>
struct longer_prec<int, realT> {
	static_assert(std::is_floating_point<realT>::value,
	"Template paramter is not numerical type");

	using type = realT;
};

template<typename realT>
struct longer_prec<realT, int> {
	static_assert(std::is_floating_point<realT>::value,
	"Template paramter is not numerical type");

	using type = realT;
};

template<>
struct longer_prec<int, int> {
	using type = int;
};


template<typename l_realT, typename r_realT>
using longer_prec_t = typename longer_prec<l_realT, r_realT>::type;

}


#endif // COFFEE_MAKERS_COMPILE_TIME_RETURN_TYPE_TRAITS_HPP
