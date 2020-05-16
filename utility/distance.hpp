#ifndef COFFEE_MAKERS_UTILITY_DISTANCE_HPP
#define COFFEE_MAKERS_UTILITY_DISTANCE_HPP

#include<coffee-makers/Containers/Containers.hpp>
#include<cmath>

namespace makers {

template<typename realT, int Dim>
inline realT distance(const Vector<realT, Dim>& vec) {
	const Matrix<realT, 1, 1>& dot_product = vec.transpose() * vec;
	return std::sqrt(dot_product.at(0));
}

template<typename realT, int Dim>
inline realT square_distance(const Vector<realT, Dim>& vec) {
	const Matrix<realT, 1, 1>& dot_product = vec.transpose() * vec;
	return dot_product.at(0);
}

}

#endif //COFFEE_MAKERS_UTILITY_DISTANCE_HPP
