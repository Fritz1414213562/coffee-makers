#ifndef COFFEE_MAKERS_FIXED_VECTOR_HPP
#define COFFEE_MAKERS_FIXED_VECTOR_HPP
#include<coffee-makers/Containers/FixedMatrix.hpp>
#include<cmath>


namespace makers {

template<typename realT, std::size_t dim>
using FixedVector = FixedMatrix<realT, dim, 1>;


template<typename Vec>
struct element_type_of {
	using type = typename Vec::element_type;
};


template<typename Vec>
struct column_size_of {
	constexpr static std::size_t value = Vec::Col_CompileTime;
};



template<typename Vec>
inline typename element_type_of<Vec>::type length_sqr(const Vec& input_vector) {
	static_assert(column_size_of<Vec>::value == 1, "'length_sqr' is operated by Matrix<M, N>. (N > 1)");
	typename element_type_of<Vec>::type result = 0.0;
	for (std::size_t idx = 0; idx < input_vector.size(); ++idx) {
		result += input_vector[idx] * input_vector[idx];
	}
	return result;
}

template<typename Vec>
inline typename element_type_of<Vec>::type length(const Vec& input_vector) {
	static_assert(column_size_of<Vec>::value == 1, "'length_sqr' is operated by Matrix<M, N>. (N > 1)");
	typename element_type_of<Vec>::type result = 0.0;
	for (std::size_t idx = 0; idx < input_vector.size(); ++idx) {
		result += input_vector[idx] * input_vector[idx];
	}
	return std::sqrt(result);
}


template<typename realT, std::size_t dim>
inline realT operator*(const FixedMatrix<realT, 1, dim>& lhs, const FixedVector<realT, dim>& rhs) {
	realT result = 0.0;
	for (std::size_t idx = 0; idx < dim; ++idx)
		result += lhs[idx] * rhs[idx];
	
	return result;
}



}

#endif /* COFFEE_MAKERS_FIXED_VECTOR_HPP */
