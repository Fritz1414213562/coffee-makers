#ifndef COFFEE_MAKERS_VECTOR_HPP
#define COFFEE_MAKERS_VECTOR_HPP
#include<coffee-makers/Containers/Matrix.hpp>

namespace makers {

template<typename realT, int Dim>
using Vector = Matrix<realT, Dim, 1>;

//template<typename realT, int r_Dim, int l_Dim>
//inline realT dot(const Vector<realT, r_Dim>& rhs, const Vector<realT, l_Dim>& lhs) {
//
//	if constexpr ((r_Dim > 0) && (l_Dim > 0)) {
//		static_assert(r_Dim == l_Dim,
//		"Invalid operation: The vectors of Both sides are different.");
//	}
//	else if constexpr ((r_Dim < 0) || (l_Dim < 0)) {
//		if (rhs.size() != lhs.size())
//			throw std::invalid_argument(
//			"Invalid operation: The vectors of Both sides are different.");
//	}
//	const Matrix<realT, 1, 1>& dot_product = rhs.transpose() * lhs;
//	return dot_product.at(0);
//}

}

#endif /* COFFEE_MAKERS_VECTOR_HPP */
