#ifndef COFFEE_MAKERS_UTILITY_INNER_PRODUCT_HPP
#define COFFEE_MAKERS_UTILITY_INNER_PRODUCT_HPP

#include<coffee-makers/Containers/Containers.hpp>
#include<coffee-makers/CompileTimeCalculation/convertible.hpp>
#include<coffee-makers/CompileTimeCalculation/return_type_traits.hpp>

namespace makers {

template<typename l_realT, int l_Dim, typename r_realT, int r_Dim>
inline auto inner_product(const Vector<l_realT, l_Dim>& lhs, const Vector<r_realT, r_Dim>& rhs)
	-> decltype(
		convertible<l_realT, r_realT>::value,
		longer_prec_t<l_realT, r_realT>()) {
	
	if constexpr ((l_Dim > 0) && (r_Dim > 0))
		static_assert(l_Dim == r_Dim,
		"Invalid operation: The vectors of Both sides are different.");
	else if constexpr ((l_Dim < 0) || (r_Dim < 0))
		if (lhs.size() != rhs.size())
			throw std::invalid_argument(
			"Invalid operation: The vectors of Both sides are different.");
	
	using returnT = longer_prec_t<l_realT, r_realT>;

	const Matrix<returnT, 1, 1>& dot_product = lhs.transpose() * rhs;
	return dot_product.at(0);
}

}


#endif // COFFEE_MAKERS_UTILITY_INNER_PRODUCT_HPP
