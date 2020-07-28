#ifndef COFFEE_MAKERS_UTILITY_SORT_HPP
#define COFFEE_MAKERS_UTILITY_SORT_HPP
#include<coffee-makers/Containers/Matrix.hpp>
#include<coffee-makers/Containers/Vector.hpp>
#include<algorithm>


namespace makers {


template<typename realT, int dimN>
Vector<int, dimN> sort_Indices(const Vector<realT, dimN>& vec) {

	Vector<int, dimN> retval(vec.size());

	for (int idx = 0; idx < vec.size(); ++idx)
		retval[idx] = idx;

	std::sort(
		retval._matrix.begin(),
		retval._matrix.end(),
		[&](const int& lhs, const int& rhs) {return (vec._matrix[lhs] < vec._matrix[rhs]);}
	);

	return retval;
}


template<typename realT, int Row, int Col>
Matrix<realT, Row, Col> sort_Column(
	const Matrix<realT, Row, Col>& mat, const Vector<int, Col>& sort_indices) {

	if constexpr (Col < 0)
		if (mat.cols() != sort_indices.size())
			throw std::out_of_range(
			"The column size of matrix is not consistent with that of indices");

	Matrix<realT, Row, Col> retval(mat.rows(), mat.cols());
	for (int i_col = 0; i_col < mat.cols(); ++i_col)
		retval.col_copy(i_col, mat.col(sort_indices[i_col]));
	
	return retval;
}


template<typename realT, int Row, int Col>
Matrix<realT, Row, Col> sort_Row(
	const Matrix<realT, Row, Col>& mat, const Vector<int, Row>& sort_indices) {
	
	if constexpr (Row < 0)
		if (mat.rows() != sort_indices.size())
			throw std::out_of_range(
			"The row size of matrix is not consistent with that of indices.");
	
	Matrix<realT, Row, Col> retval(mat.rows(), mat.cols());

	for (int i_row = 0; i_row < mat.rows(); ++i_row)
		retval.row_copy(i_row, mat.row(sort_indices[i_row]));
	
	return retval;
}


}

#endif // COFFEE_MAKERS_UTILITY_SORT_HPP
