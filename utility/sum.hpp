#ifndef COFFEE_MAKERS_UTILITY_SUM_HPP
#define COFFEE_MAKERS_UTILITY_SUM_HPP

#include<coffee-makers/Containers/Containers.hpp>

namespace makers {

template<typename realT, int Row, int Col>
inline Vector<realT, Col> row_sum(const Matrix<realT, Row, Col>& matrix) {

	Vector<realT, Col> retval(matrix.cols());
	
	for (int i_row = 0; i_row < matrix.rows(); ++i_row)
		retval += matrix.row(i_row);
	return retval;
}


template<typename realT, int Row, int Col>
inline Vector<realT, Row> col_sum(const Matrix<realT, Row, Col>& matrix) {

	Vector<realT, Row> retval(matrix.rows());

	for (int i_col = 0; i_col < matrix.cols(); ++i_col)
		retval += matrix.col(i_col);

	return retval;
}


template<typename realT, int Row, int Col>
inline realT sum(const Matrix<realT, Row, Col>& matrix) {

	realT retval = 0.;

	for (int idx = 0; idx < matrix.size(); ++idx)
		retval += matrix[idx];

	return retval;
}

}


#endif // COFFEE_MAKERS_UTILITY_SUM_HPP
