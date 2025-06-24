#ifndef COFFEE_MAKERS_UTILITY_SQUAREROOT_HPP
#define COFFEE_MAKERS_UTILITY_SQUAREROOT_HPP
#include<coffee-makers/Containers/Containers.hpp>
#include<cmath>


namespace makers {

template<typename realT, int Row, int Col>
inline Matrix<realT, Row, Col> squareroot(const Matrix<realT, Row, Col>& matrix) {

	Matrix<realT, Row, Col> retval(matrix.rows(), matrix.cols());

	for (int idx = 0; idx < matrix.size(); ++idx)
		retval[idx] = std::sqrt(matrix[idx]);
	
	return retval;
}

}

#endif // COFFEE_MAKERS_UTILITY_SQUAREROOT_HPP
