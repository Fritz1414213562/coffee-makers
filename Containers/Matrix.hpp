#ifndef COFFEE_MAKERS_MATRIX_HPP
#define COFFEE_MAKERS_MATRIX_HPP
#include<coffee-makers/CompileTimeCalculation/conjunction.hpp>
#include<coffee-makers/CompileTimeCalculation/disjunction.hpp>
#include<type_traits>
#include<vector>
#include<array>
#include<iostream>
#include<utility>

namespace makers {



enum { Variable = -1 };


template<typename realT, int Row, int Col>
class Matrix {


	static_assert((Row >= -1) && (Col >= -1),
	"Invalid template parameter, smaller than -1");
	static_assert((Row != 0),
	"Invalid template parameter, row size is 0");
	static_assert((Col != 0),
	"Invalid template parameter, column size is 0");

public:

// type
	using scalar_type = realT;
	using container_type =
		typename std::conditional<(Row < 0) || (Col < 0),
			std::vector<realT>, std::array<realT, Row * Col>>::type;


	int rows() const {return _row;}
	int cols() const {return _col;}
	int size() const {return _size;}

	constexpr static int Row_CompileTime = Row;
	constexpr static int Col_CompileTime = Col;
	constexpr static int Size_CompileTime = Row * Col;


public:


/*

	Constructor

*/



	// constructor for a fixed matrix
	Matrix();

	// constructor for a variable matrix
	template<typename Integer, typename std::enable_if<
		value_conjunction<
			std::is_same<Integer, int>::value,
			value_disjunction<(Row < 0), (Col < 0)>::value>::value, std::nullptr_t>::type = nullptr>
	Matrix(const Integer& runtime_row, const Integer& runtime_col);


	// redundant constructor for a fixed matrix
	template<typename Integer, typename std::enable_if<
		value_conjunction<
			std::is_same<Integer, int>::value,
			value_conjunction<(Row > 0), (Col > 0)>::value>::value, std::nullptr_t>::type = nullptr>
	Matrix(const Integer& runtime_row, const Integer& runtime_col);

	template<typename Integer, typename std::enable_if<
		value_conjunction<
			std::is_same<Integer, int>::value, (Row < 0), (Col == 1)>::value,
			std::nullptr_t>::type = nullptr>
	Matrix(const Integer& runtime_dim);


	template<typename Integer, typename std::enable_if<
		value_conjunction<
			std::is_same<Integer, int>::value, (Row > 0), (Col == 1)>::value,
			std::nullptr_t>::type = nullptr>
	Matrix(const Integer& runtime_dim);


	// copy constructor
	Matrix(const Matrix& input_matrix) = default;


	template<typename T, int r_Row, int r_Col>
	Matrix(const Matrix<T, r_Row, r_Col>& input_matrix);

	
	// destructor
	~Matrix() = default;




	Matrix& operator=(const Matrix& rhs_matrix) = default;
	template<typename T, int r_Row, int r_Col>
	Matrix<realT, Row, Col>& operator=(const Matrix<T, r_Row, r_Col>& rhs);

//	template<typename T, int Dim>
//	Matrix<realT, Row, 1>& operator=(const Matrix<T, 1, Dim>& rhs);

	Matrix& operator+=(const Matrix& rhs_matrix);
	template<typename T, int r_Row, int r_Col>
	Matrix<realT, Row, Col>& operator+=(const Matrix<T, r_Row, r_Col>& rhs);

	Matrix& operator-=(const Matrix& rhs_matrix);
	template<typename T, int r_Row, int r_Col>
	Matrix<realT, Row, Col>& operator-=(const Matrix<T, r_Row, r_Col>& rhs);

	Matrix& operator+=(const scalar_type& rhs);
	Matrix& operator-=(const scalar_type& rhs);
	Matrix& operator*=(const scalar_type& rhs);
	Matrix& operator/=(const scalar_type& rhs);

	scalar_type  at(const int& idx, const int& jdx) const;
	scalar_type& at(const int& idx, const int& jdx);
	scalar_type  operator()(const int& idx, const int& jdx) const;
	scalar_type& operator()(const int& idx, const int& jdx);

	scalar_type  at(const int& idx) const {return _matrix.at(idx);}
	scalar_type& at(const int& idx) {return _matrix.at(idx);}
	scalar_type  operator[](const int& idx) const {return _matrix[idx];}
	scalar_type& operator[](const int& idx) {return _matrix[idx];}



	Matrix<realT, Col, 1> row(const int& idx) const;
	Matrix<realT, Row, 1> col(const int& idx) const;
	void row_swap(const int& i_row, const int& j_row);
	void col_swap(const int& i_col, const int& j_col);
	template<typename T>
	void row_copy(const int& i_row, const Matrix<T, Col, 1>& copied_vec);
	template<typename T>
	void col_copy(const int& i_col, const Matrix<T, Row, 1>& copied_vec);


	// 
	Matrix<realT, Col, Row> transpose() const;


//	scalar_type determinant() const;
//
//	Matrix inverse() const;

	scalar_type determinant() const;
	Matrix inverse() const;



private:
	// friend function

	template<typename T, int N>
	friend Matrix<int, N, 1> sort_Indices(const Matrix<T, N, 1>& vec);


private:
	int _row;
	int _col;
	int _size;
	container_type _matrix;


	template<int Dim, int DummyValue = 0>
	struct DeterminantCalculator;

	template<int Dim, int DummyValue = 0>
	struct InverseCalculator;


};









// implement





template<typename realT, int Row, int Col>
Matrix<realT, Row, Col>::Matrix() :
	_row(Row_CompileTime),
	_col(Col_CompileTime),
	_size(Size_CompileTime),
	_matrix{{}} {
	static_assert((Row_CompileTime > 0) && (Col_CompileTime > 0),
	"In the case that 'makers::Matrix' is variable, the default constructor is forbidden.");
	_matrix.fill(0);
}



// constructor for a variable matrix
template<typename realT, int Row, int Col>
template<typename Integer, typename std::enable_if<
	value_conjunction<
		std::is_same<Integer, int>::value,
		value_disjunction<(Row < 0), (Col < 0)>::value>::value, std::nullptr_t>::type>
Matrix<realT, Row, Col>::Matrix(const Integer& runtime_row, const Integer& runtime_col) :
	_row(runtime_row),
	_col(runtime_col),
	_size(runtime_row * runtime_col),
	_matrix(runtime_row * runtime_col, 0) {


	if ((runtime_row <= 0) || (runtime_col <= 0))
		throw std::invalid_argument("Invalid arguments: row size or column size is not positive.");
	if constexpr (Row_CompileTime > 0) {
		if (Row_CompileTime != runtime_row)
			throw std::invalid_argument(
				"The row size at compile time is not consistent with that at run time.");
	}
	if constexpr (Col_CompileTime > 0) {
		if (Col_CompileTime != runtime_col)
			throw std::invalid_argument(
				"The column size at compile time is not consistent with that at run time.");
	}

}



// redundant constructor for a fixed matrix
template<typename realT, int Row, int Col>
template<typename Integer, typename std::enable_if<
	value_conjunction<
		std::is_same<Integer, int>::value,
		value_conjunction<(Row > 0), (Col > 0)>::value>::value, std::nullptr_t>::type>
Matrix<realT, Row, Col>::Matrix(const Integer& runtime_row, const Integer& runtime_col) :
	_row(runtime_row),
	_col(runtime_col),
	_size(runtime_row * runtime_col),
	_matrix{{}} {

//	static_assert((Row_CompileTime < 0) && (Col_CompileTime < 0),
//	"In the case that 'makers::Matrix' is fixed, the constructor (int const&, int const&) is forbidden.");

	if ((runtime_row <= 0) || (runtime_col <= 0))
		throw std::invalid_argument("Invalid arguments: row size or column size is not positive.");
	if (Row_CompileTime != runtime_row)
		throw std::invalid_argument(
			"The row size at compile time is not consistent with that at run time.");
	
	if (Col_CompileTime != runtime_col)
		throw std::invalid_argument(
			"The column size at compile time is not consistent with that at run time.");

	_matrix.fill(0);
}


// for fixed vector
template<typename realT, int Row, int Col>
template<typename Integer, typename std::enable_if<
	value_conjunction<std::is_same<Integer, int>::value,
	(Row > 0), (Col == 1)>::value, std::nullptr_t>::type>
Matrix<realT, Row, Col>::Matrix(const Integer& runtime_dim) :
	_row(runtime_dim),
	_col(Col),
	_size(runtime_dim),
	_matrix{{}} {

//	static_assert((Row_CompileTime < 0) && (Col_CompileTime < 0),
//	"In the case that 'makers::Matrix' is fixed, the constructor (int const&, int const&) is forbidden.");

	if (runtime_dim <= 0)
		throw std::invalid_argument("Invalid arguments: row size or column size is not positive.");
	_matrix.fill(0);
}


// for variable vector
template<typename realT, int Row, int Col>
template<typename Integer, typename std::enable_if<
	value_conjunction<std::is_same<Integer, int>::value,
	(Row < 0), (Col == 1)>::value, std::nullptr_t>::type>
Matrix<realT, Row, Col>::Matrix(const Integer& runtime_dim) :
	_row(runtime_dim),
	_col(Col),
	_size(runtime_dim),
	_matrix(runtime_dim, 0) {

//	static_assert((Row_CompileTime < 0) && (Col_CompileTime < 0),
//	"In the case that 'makers::Matrix' is fixed, the constructor (int const&, int const&) is forbidden.");

	if (runtime_dim <= 0)
		throw std::invalid_argument("Invalid arguments: row size or column size is not positive.");
}


template<typename realT, int Row, int Col>
template<typename T, int r_Row, int r_Col>
Matrix<realT, Row, Col>::Matrix(const Matrix<T, r_Row, r_Col>& input_matrix) {
	static_assert(std::is_convertible<realT, T>::value,
	"This scalar type convertion is invalid.");

//	std::cout << "Copy Constructor" << std::endl;

	if constexpr (Row < 0)
		(*this)._row = input_matrix.rows();

	else if constexpr (value_conjunction<(Row > 0), (r_Row > 0)>::value) {
		static_assert(Row == r_Row,
		"Copying from a different row size matrix");
		(*this)._row = Row;
	}
	else if constexpr (value_conjunction<(Row > 0), (r_Row < 0)>::value) {
		if (Row != input_matrix.rows())
			throw std::invalid_argument(
			"Invalid operation: Copy from a different row size matrix");
		(*this)._row = Row;
	}


	if constexpr (Col < 0)
		(*this)._col = input_matrix.cols();

	else if constexpr (value_conjunction<(Col > 0), (r_Col > 0)>::value) {
		static_assert(Col == r_Col,
		"Copying from a different column size matrix");
		(*this)._col = Col;
	}

	else if constexpr (value_conjunction<(Col > 0), (r_Col < 0)>::value) {
		if (Col != input_matrix.cols())
			throw std::invalid_argument(
			"Invalid operation: Copy from a different column size matrix");
		(*this)._col = Col;
	}

	(*this)._size = (*this)._row * (*this)._col;

	if constexpr ((Row < 0) || (Col < 0)) {
		(*this)._matrix.resize((*this)._size);
	}
	
	for (int idx = 0; idx < (*this)._size; ++idx)
		(*this)[idx] = input_matrix[idx];
}





// copy operator
template<typename realT, int Row, int Col>
template<typename T, int r_Row, int r_Col>
Matrix<realT, Row, Col>&
Matrix<realT, Row, Col>::operator=(const Matrix<T, r_Row, r_Col>& rhs) {

	static_assert(std::is_convertible<realT, T>::value,
	"This scalar type convertion is invalid.");

	if constexpr (Row * r_Row < 0)
		if ((*this)._row != rhs.rows())
			throw std::invalid_argument(
			"Invalid operation: Copy from a different row size matrix");

	if constexpr (Col * r_Col < 0)
		if ((*this)._col != rhs.cols())
			throw std::invalid_argument(
			"Invalid operation: Copy from a different column size matrix");

	if constexpr ((Row > 0) && (r_Row > 0))
		static_assert(Row == r_Row,
		"Copying from a different row size matrix");
	if constexpr ((Col > 0) && (r_Col > 0))
		static_assert(Col == r_Col,
		"Copying from a different column size matrix");
	
	Matrix<realT, Row, Col> retval(rhs.rows(), rhs.cols());
	for (int idx = 0; idx < rhs.size(); ++idx)
		retval[idx] = static_cast<realT>(rhs[idx]);
	
	return retval;

}



//template<typename realT, int Row, int Col>
//float& operator=(const Matrix<realT, 1, 1>& mat) {
//
//	if constexpr ((Row > 0) && (Col > 0)) {
//		static_assert((Row == 1) && (Col == 1), "Invalid operation: Copy from matrix to scalar");
//	}
//	else if constexpr ((Row < 0) || (Col < 0)) {
//		if ((mat.rows() != 1) || (mat.cols() != 1))
//			throw std::invalid_argument(
//			"Invalid operation: Copy from matrix to scalar");
//	}
//	return mat.at(0);
//}



//template<typename realT, int Row, int Col>
//template<typename T, int Dim>
//Matrix<realT, Row, 1>&
//Matrix<realT, Row, Col>::operator=(const Matrix<T, 1, Dim>& rhs) {
//	static_assert((Col == 1) || (Col < 0),
//	"Invalid operation: conversion from 2D to 1D.");
//
//
//}



template<typename realT, int Row, int Col>
Matrix<realT, Row, Col>&
Matrix<realT, Row, Col>::operator+=(const Matrix<realT, Row, Col>& rhs_matrix) {

	if constexpr ((Row > 0) && (Col > 0)) {
		for (int idx = 0; idx < Row; ++idx)
			for (int jdx = 0; jdx < Col; ++jdx)
				(*this)(idx, jdx) += rhs_matrix(idx, jdx);
	}
	else if constexpr ((Row < 0) && (Col > 0)) {
		if ((*this).rows() != rhs_matrix.rows())
			throw std::invalid_argument("Invalid operation between different size matrices");

		for (int idx = 0; idx < (*this).rows(); ++idx)
			for (int jdx = 0; jdx < Col; ++jdx)
				(*this)(idx, jdx) += rhs_matrix(idx, jdx);
	}
	else if constexpr ((Row > 0) && (Col < 0)) {
		if ((*this).cols() != rhs_matrix.cols())
			throw std::invalid_argument("Invalid operation between different size matrices");

		for (int idx = 0; idx < Row; ++idx)
			for (int jdx = 0; jdx < (*this).cols(); ++jdx)
				(*this)(idx, jdx) += rhs_matrix(idx, jdx);
	}
	else if constexpr ((Row < 0) && (Col < 0)) {
		if (((*this).rows() != rhs_matrix.rows()) || ((*this).cols() != rhs_matrix.cols()))
			throw std::invalid_argument("Invalid operation between different size matrices");

		for (int idx = 0; idx < (*this).rows(); ++idx)
			for (int jdx = 0; jdx < (*this).cols(); ++jdx)
				(*this)(idx, jdx) += rhs_matrix(idx, jdx);
	}

	return *this;
}


template<typename realT, int Row, int Col>
template<typename T, int r_Row, int r_Col>
Matrix<realT, Row, Col>&
Matrix<realT, Row, Col>::operator+=(const Matrix<T, r_Row, r_Col>& rhs_matrix) {

	static_assert(std::conjunction<
		std::is_convertible<realT, T>,
		value_disjunction<(Row * r_Row < 0), (Col * r_Col < 0)>>::value,
		"The matrix size is different between left-hand and right-hand sides"
	);

	if constexpr (Row * r_Row < 0) {
		if ((*this).rows() != rhs_matrix.rows())
			throw std::invalid_argument("Invalid operation between different size matrices");

		if constexpr (Col * r_Col > 0) {

			static_assert(Col == r_Col,
			"The Column size is different between left-hand and right-hand sides");

			for (int idx = 0; idx < (*this).rows(); ++idx)
				for (int jdx = 0; jdx < Col; ++jdx)
					(*this)(idx, jdx) += rhs_matrix(idx, jdx);
		}

		else if constexpr (Col * r_Col < 0) {
			if ((*this).cols() != rhs_matrix.cols())
				throw std::invalid_argument("Invalid operation between differene size matrices");

			for (int idx = 0; idx < (*this).rows(); ++idx)
				for (int jdx = 0; jdx < (*this).cols(); ++jdx)
					(*this)(idx, jdx) += rhs_matrix(idx, jdx);
		}
	}

	else if constexpr (Row * r_Row > 0) {
		static_assert(value_conjunction<(Row * r_Row > 0), (Col * r_Col < 0)>::value,
		"Something wrong! First assertion failed.");
		static_assert(Row == r_Row,
		"The Row size is different between left-hand and right-hand sides");

		for (int idx = 0; idx < Row; ++idx)
			for (int jdx = 0; jdx < (*this).cols(); ++jdx)
				(*this)(idx, jdx) += rhs_matrix(idx, jdx);
	}


	return *this;
}


template<typename realT, int Row, int Col>
Matrix<realT, Row, Col>&
Matrix<realT, Row, Col>::operator-=(const Matrix<realT, Row, Col>& rhs_matrix) {

	if constexpr ((Row > 0) && (Col > 0)) {
		for (int idx = 0; idx < Row; ++idx)
			for (int jdx = 0; jdx < Col; ++jdx)
				(*this)(idx, jdx) -= rhs_matrix(idx, jdx);
	}
	else if constexpr ((Row < 0) && (Col > 0)) {
		if ((*this).rows() != rhs_matrix.rows())
			throw std::invalid_argument("Invalid operation between different size matrices");

		for (int idx = 0; idx < (*this).rows(); ++idx)
			for (int jdx = 0; jdx < Col; ++jdx)
				(*this)(idx, jdx) -= rhs_matrix(idx, jdx);
	}
	else if constexpr ((Row > 0) && (Col < 0)) {
		if ((*this).cols() != rhs_matrix.cols())
			throw std::invalid_argument("Invalid operation between different size matrices");

		for (int idx = 0; idx < Row; ++idx)
			for (int jdx = 0; jdx < (*this).cols(); ++jdx)
				(*this)(idx, jdx) -= rhs_matrix(idx, jdx);
	}
	else if constexpr ((Row < 0) && (Col < 0)) {
		if (((*this).rows() != rhs_matrix.rows()) || ((*this).cols() != rhs_matrix.cols()))
			throw std::invalid_argument("Invalid operation between different size matrices");

		for (int idx = 0; idx < (*this).rows(); ++idx)
			for (int jdx = 0; jdx < (*this).cols(); ++jdx)
				(*this)(idx, jdx) -= rhs_matrix(idx, jdx);
	}

	return *this;
}




template<typename realT, int Row, int Col>
template<typename T, int r_Row, int r_Col>
Matrix<realT, Row, Col>&
Matrix<realT, Row, Col>::operator-=(const Matrix<T, r_Row, r_Col>& rhs_matrix) {

	static_assert(std::conjunction<
		std::is_convertible<realT, T>,
		value_disjunction<(Row * r_Row < 0), (Col * r_Col < 0)>>::value,
		"The matrix size is different between left-hand and right-hand sides"
	);

	if constexpr (Row * r_Row < 0) {
		if ((*this).rows() != rhs_matrix.rows())
			throw std::invalid_argument("Invalid operation between different size matrices");

		if constexpr (Col * r_Col > 0) {

			static_assert(Col == r_Col,
			"The Column size is different between left-hand and right-hand sides");

			for (int idx = 0; idx < (*this).rows(); ++idx)
				for (int jdx = 0; jdx < Col; ++jdx)
					(*this)(idx, jdx) -= rhs_matrix(idx, jdx);
		}

		else if constexpr (Col * r_Col < 0) {
			if ((*this).cols() != rhs_matrix.cols())
				throw std::invalid_argument("Invalid operation between differene size matrices");

			for (int idx = 0; idx < (*this).rows(); ++idx)
				for (int jdx = 0; jdx < (*this).cols(); ++jdx)
					(*this)(idx, jdx) -= rhs_matrix(idx, jdx);
		}
	}

	else if constexpr (Row * r_Row > 0) {
		static_assert(value_conjunction<(Row * r_Row > 0), (Col * r_Col < 0)>::value,
		"Something wrong! First assertion failed.");
		static_assert(Row == r_Row,
		"The Row size is different between left-hand and right-hand sides");

		for (int idx = 0; idx < Row; ++idx)
			for (int jdx = 0; jdx < (*this).cols(); ++jdx)
				(*this)(idx, jdx) -= rhs_matrix(idx, jdx);
	}


	return *this;
}




template<typename realT, int Row, int Col>
Matrix<realT, Row, Col>&
Matrix<realT, Row, Col>::operator+=(const scalar_type& rhs) {

	for (int idx = 0; idx < (*this).rows(); ++idx)
		for (int jdx = 0; jdx < (*this).cols(); ++jdx)
			(*this)(idx, jdx) += rhs;

	return *this;
}


template<typename realT, int Row, int Col>
Matrix<realT, Row, Col>&
Matrix<realT, Row, Col>::operator-=(const scalar_type& rhs) {

	for (int idx = 0; idx < (*this).rows(); ++idx)
		for (int jdx = 0; jdx < (*this).cols(); ++jdx)
			(*this)(idx, jdx) -= rhs;

	return *this;
}



template<typename realT, int Row, int Col>
Matrix<realT, Row, Col>&
Matrix<realT, Row, Col>::operator*=(const scalar_type& rhs) {

	for (int idx = 0; idx < (*this).rows(); ++idx)
		for (int jdx = 0; jdx < (*this).cols(); ++jdx)
			(*this)(idx, jdx) *= rhs;

	return *this;
}


template<typename realT, int Row, int Col>
Matrix<realT, Row, Col>&
Matrix<realT, Row, Col>::operator/=(const scalar_type& rhs) {

	for (int idx = 0; idx < (*this).rows(); ++idx)
		for (int jdx = 0; jdx < (*this).cols(); ++jdx)
			(*this)(idx, jdx) /= rhs;

	return *this;
}


template<typename realT, int Row, int Col>
Matrix<realT, Row, Col>
operator+(const Matrix<realT, Row, Col>& lhs, const Matrix<realT, Row, Col>& rhs) {

	if ((lhs.rows() != rhs.rows()) || (lhs.cols() != rhs.cols()))
		throw std::invalid_argument("Invalid operation between different size matrices");
	Matrix<realT, Row, Col> retval(lhs.rows(), lhs.cols());

	for (int idx = 0; idx < lhs.rows(); ++idx)
		for (int jdx = 0; jdx < lhs.cols(); ++jdx)
			retval(idx, jdx) = lhs(idx, jdx) + rhs(idx, jdx);

	return retval;
}



template<typename realT, int Row, int Col>
Matrix<realT, Row, Col>
operator-(const Matrix<realT, Row, Col>& lhs, const Matrix<realT, Row, Col>& rhs) {

	if ((lhs.rows() != rhs.rows()) || (lhs.cols() != rhs.cols()))
		throw std::invalid_argument("Invalid operation between different size matrices");
	Matrix<realT, Row, Col> retval(lhs.rows(), lhs.cols());

	for (int idx = 0; idx < lhs.rows(); ++idx)
		for (int jdx = 0; jdx < lhs.cols(); ++jdx)
			retval(idx, jdx) = lhs(idx, jdx) - rhs(idx, jdx);

	return retval;
}







template<typename realT, int Row, int Col>
Matrix<realT, Row, Col>
operator*(const Matrix<realT, Row, Col>& lhs, const realT& rhs) {

	Matrix<realT, Row, Col> retval(lhs.rows(), lhs.cols());

	for (int idx = 0; idx < lhs.rows(); ++idx)
		for (int jdx = 0; jdx < lhs.cols(); ++jdx)
			retval(idx, jdx) = lhs(idx, jdx) * rhs;

	return retval;
}


template<typename realT, int Row, int Col>
Matrix<realT, Row, Col>
operator*(const realT& lhs, const Matrix<realT, Row, Col>& rhs) {

	Matrix<realT, Row, Col> retval(rhs.rows(), rhs.cols());

	for (int idx = 0; idx < rhs.rows(); ++idx)
		for (int jdx = 0; jdx < rhs.cols(); ++jdx)
			retval(idx, jdx) = lhs * rhs(idx, jdx);

	return retval;
}


template<typename realT, int Row, int Col>
Matrix<realT, Row, Col>
operator/(const Matrix<realT, Row, Col>& lhs, const realT& rhs) {

	Matrix<realT, Row, Col> retval(lhs.rows(), lhs.cols());

	for (int idx = 0; idx < lhs.rows(); ++idx)
		for (int jdx = 0; jdx < lhs.cols(); ++jdx)
			retval(idx, jdx) = lhs(idx, jdx) / rhs;

	return retval;
}


template<typename realT, int LeftRow, int LeftCol, int RightRow, int RightCol>
Matrix<realT, LeftRow, RightCol>
operator*(const Matrix<realT, LeftRow, LeftCol>& lhs, const Matrix<realT, RightRow, RightCol>& rhs) {


	if constexpr ((LeftCol > 0) && (RightRow > 0)) {
		static_assert(LeftCol == RightRow,
		"The column size of left-hand side is not consistent with the right-hand row size");
	}
	else if constexpr ((LeftCol < 0) || (RightRow < 0)) {
		if (lhs.cols() != rhs.rows())
			throw std::invalid_argument(
			"Invalid operation between matrices which have different row size from the other column size");
	}

	
	Matrix<realT, LeftRow, RightCol> retval(lhs.rows(), rhs.cols());

	for (int idx = 0; idx < lhs.rows(); ++idx)
		for (int jdx = 0; jdx < rhs.cols(); ++jdx)
			for (int kdx = 0; kdx < lhs.cols(); ++kdx)
				retval(idx, jdx) += lhs(idx, kdx) * rhs(kdx, jdx);

	return retval;
}


template<typename realT, int Row, int Col>
typename Matrix<realT, Row, Col>::scalar_type
Matrix<realT, Row, Col>::at(const int& idx, const int& jdx) const {
	return this->_matrix.at(idx * (this->_col) + jdx);
}


template<typename realT, int Row, int Col>
typename Matrix<realT, Row, Col>::scalar_type&
Matrix<realT, Row, Col>::at(const int& idx, const int& jdx) {
	return this->_matrix.at(idx * (this->_col) + jdx);
}


template<typename realT, int Row, int Col>
typename Matrix<realT, Row, Col>::scalar_type
Matrix<realT, Row, Col>::operator()(const int& idx, const int& jdx) const {
	return this->_matrix[idx * (this->_col) + jdx];
}


template<typename realT, int Row, int Col>
typename Matrix<realT, Row, Col>::scalar_type&
Matrix<realT, Row, Col>::operator()(const int& idx, const int& jdx) {
	return this->_matrix[idx * (this->_col) + jdx];
}


template<typename realT, int Row, int Col>
Matrix<realT, Col, 1> Matrix<realT, Row, Col>::row(const int& idx) const {

	if (idx >= (*this)._row)
		throw std::out_of_range("The index is out of column range");

	Matrix<realT, Col, 1> retval((*this)._col, 1);
	for (int jdx = 0; jdx < (*this)._col; ++jdx)
		retval[jdx] = (*this)(idx, jdx);
	return retval;
}



template<typename realT, int Row, int Col>
Matrix<realT, Row, 1> Matrix<realT, Row, Col>::col(const int& idx) const {
	if (idx >= (*this)._col)
		throw std::out_of_range("The index in 'col(const int&)' is out of column range");

	Matrix<realT, Row, 1> retval((*this)._row, 1);
	for (int jdx = 0; jdx < (*this)._row; ++jdx)
		retval[jdx] = (*this)(jdx, idx);
	return retval;
}



template<typename realT, int Row, int Col>
void Matrix<realT, Row, Col>::row_swap(const int& i_row, const int& j_row) {
	if ((i_row >= (*this)._row) || (j_row >= (*this)._row))
		throw std::out_of_range("The index in 'row_swap(const int&, const int&' is out of row range");

	for (int idx = 0; idx < (*this)._col; ++idx)
		std::swap((*this).at(i_row, idx), (*this).at(j_row, idx));
}



template<typename realT, int Row, int Col>
void Matrix<realT, Row, Col>::col_swap(const int& i_col, const int& j_col) {
	if ((i_col >= (*this)._col) || (j_col >= (*this)._col))
		throw std::out_of_range("The index in 'col_swap(const int&, const int&' is out of column range");

	for (int idx = 0; idx < (*this)._row; ++idx)
		std::swap((*this).at(idx, i_col), (*this).at(idx, j_col));
}



template<typename realT, int Row, int Col>
template<typename T>
void Matrix<realT, Row, Col>::row_copy(const int& i_row, const Matrix<T, Col, 1>& copied_vec) {
	static_assert(std::is_convertible<realT, T>::value,
	"The scalar type of Matrix is not convertible from that of copied vector.");

	if constexpr (Col < 0)
		if ((*this)._col != copied_vec.size())
			throw std::out_of_range("The column size is not consistent between source and target");

	if (i_row >= (*this)._row)
		throw std::out_of_range("The index of row is out of row range");

	for (int idx = 0; idx < (*this)._col; ++idx)
		(*this).at(i_row, idx) = copied_vec[idx];
}


template<typename realT, int Row, int Col>
template<typename T>
void Matrix<realT, Row, Col>::col_copy(const int& i_col, const Matrix<T, Row, 1>& copied_vec) {
	static_assert(std::is_convertible<realT, T>::value,
	"The scalar type of Matrix is not convertible from that of copied vector.");

	if constexpr (Row < 0)
		if ((*this)._row != copied_vec.size())
			throw std::out_of_range("The row size is not consistent between source and target");

	if (i_col >= (*this)._col)
		throw std::out_of_range("The index of column is out of row range");
	
	for (int idx = 0; idx < (*this)._row; ++idx)
		(*this).at(idx, i_col) = copied_vec[idx];
}



template<typename realT, int Row, int Col>
Matrix<realT, Col, Row> Matrix<realT, Row, Col>::transpose() const {
//Matrix<realT, row, col> Matrix<realT, row, col>::transpose() const {

	Matrix<realT, Col, Row> retval(this->_col, this->_row);

	for (int idx = 0; idx < this->_row; ++idx)
		for (int jdx = 0; jdx < this->_col; ++jdx)
			retval(jdx, idx) = (*this)(idx, jdx);
	
	return retval;
}


template<typename realT, int Row, int Col>
typename Matrix<realT, Row, Col>::scalar_type Matrix<realT, Row, Col>::determinant() const {

	static_assert(value_disjunction<Row == Col, (Row < 0), (Col < 0)>::value,
	"Invalid Matrix Operation: The fixed matrix should be square");
	if constexpr ((Row < 0) || (Col < 0)) {
		if (this->_row != this->_col)
			throw std::runtime_error("The matrix shape should be square for calculating its determinant");
	}

	if constexpr (Row > 0)
		return DeterminantCalculator<Row>::calculate(*this);
	else if constexpr (Row < 0)
		return DeterminantCalculator<Col>::calculate(*this);
}


template<typename realT, int Row, int Col>
Matrix<realT, Row, Col> Matrix<realT, Row, Col>::inverse() const {

	static_assert(value_disjunction<Row == Col, (Row < 0), (Col < 0)>::value,
	"Invalid Matrix Operation: The fixed matrix should be square");
	if constexpr ((Row < 0) || (Col < 0)) {
		if (this->_row != this->_col)
		throw std::runtime_error("The matrix shape should be square for calculating its inverse matrix");
	}

	if constexpr (Row > 0)
		return InverseCalculator<Row>::calculate(*this);
	else if constexpr (Row < 0)
		return InverseCalculator<Col>::calculate(*this);
}


template<typename realT, int Row, int Col>
template<int Dim, int DummyValue>
struct Matrix<realT, Row, Col>::DeterminantCalculator {
	static scalar_type calculate(const Matrix<realT, Row, Col>& mat) {
		Matrix<realT, Row, Col> mat_copy = mat;

		scalar_type retsign = 1e0;
		for (int idx = 0; idx < mat_copy.rows() - 1; ++idx) {
			int max_pivot_idx = idx;
			realT max_pivot = 0;
			for (int jdx = idx; jdx < mat_copy.rows(); ++jdx) {
				const realT& pivot = mat_copy(jdx, idx);
				const realT& sign = (pivot > 0) - (pivot < 0);
				if (max_pivot < (pivot * sign)) {
					max_pivot = (pivot * sign);
					max_pivot_idx = jdx;
				}
			}
			if (max_pivot == 0) break; // because the determinant is always 0
			else if (max_pivot_idx != idx) {
				mat_copy.row_swap(idx, max_pivot_idx);
				retsign *= -1e0;
			}

			for (int jdx = idx + 1; jdx < mat_copy.rows(); ++jdx) {
				const scalar_type& buffer = mat_copy(jdx, idx) / mat_copy(idx, idx);
				for (int kdx = 0; kdx < mat_copy.cols(); ++kdx)
					mat_copy(jdx, kdx) -= mat_copy(idx, kdx) * buffer;
			}
		}

		scalar_type retval = retsign;


		for (int idx = 0; idx < mat_copy.rows(); ++idx)
			retval *= mat_copy(idx, idx);

		return retval;
	}
};


template<typename realT, int Row, int Col>
template<int Dim, int DummyValue>
struct Matrix<realT, Row, Col>::InverseCalculator {
	static Matrix<realT, Row, Col> calculate(const Matrix<realT, Row, Col>& mat) {
		Matrix<realT, Row, Col> mat_copy = mat;
		Matrix<realT, Row, Col> inv_mat(mat.rows(), mat.cols());

		for (int idx = 0; idx < mat.rows(); ++idx)
			inv_mat(idx, idx) = 1e0;

		for (int idx = 0; idx < mat.rows(); ++idx) {
			int max_pivot_idx = idx;
			realT max_pivot = 0;
			for (int jdx = idx; jdx < mat.rows(); ++jdx) {
				const realT& pivot = mat_copy(jdx, idx);
				const realT& sign = (pivot > 0) - (pivot < 0);
				if (max_pivot < (pivot * sign)) {
					max_pivot = (pivot * sign);
					max_pivot_idx = jdx;
				}
			}
			if (max_pivot == 0)
				std::invalid_argument(
				"Invalid operation: There is no inverse matrix");
			else if (max_pivot_idx != idx) {
				mat_copy.row_swap(idx, max_pivot_idx);
				inv_mat.row_swap(idx, max_pivot_idx);
			}
			
			scalar_type inv_diagonal = 1e0 / mat_copy(idx, idx);
			for (int jdx = 0; jdx < mat.cols(); ++jdx) {
				mat_copy(idx, jdx) *= inv_diagonal;
				inv_mat(idx, jdx) *= inv_diagonal;
			}

			for (int jdx = 0; jdx < mat.rows(); ++jdx) {
				if (idx == jdx) continue;
				scalar_type buffer = mat_copy(jdx, idx);
				for (int kdx = 0; kdx < mat.cols(); ++kdx) {
					mat_copy(jdx, kdx) -= mat_copy(idx, kdx) * buffer;
					inv_mat(jdx, kdx) -= inv_mat(idx, kdx) * buffer;
				}
			}
		}

		return inv_mat;
	}
};


template<typename realT, int Row, int Col>
template<int Dim>
struct Matrix<realT, Row, Col>::DeterminantCalculator<2, Dim> {
	static Matrix<realT, Row, Col>::scalar_type calculate(const Matrix<realT, Row, Col>& mat) {
		return (mat(0, 0) * mat(1, 1) - mat(0, 1) * mat(1, 0));
	}
};


template<typename realT, int Row, int Col>
template<int Dim>
struct Matrix<realT, Row, Col>::DeterminantCalculator<3, Dim> {
	static Matrix<realT, Row, Col>::scalar_type calculate(const Matrix<realT, Row, Col>& mat) {

		return
			mat(0, 0) * mat(1, 1) * mat(2, 2) +
			mat(1, 0) * mat(2, 1) * mat(0, 2) +
			mat(2, 0) * mat(0, 1) * mat(1, 2) -
			mat(0, 0) * mat(2, 1) * mat(1, 2) -
			mat(2, 0) * mat(1, 1) * mat(0, 2) -
			mat(1, 0) * mat(0, 1) * mat(2, 2);
	}
};


template<typename realT, int Row, int Col>
template<int Dim>
struct Matrix<realT, Row, Col>::InverseCalculator<2, Dim> {
	static Matrix<realT, 2, 2> calculate(const Matrix<realT, Row, Col>& mat) {
		const auto det_inv = 1e0 / mat.determinant();

		Matrix<realT, 2, 2> inv;
		inv(0, 0) = det_inv * mat(1, 1);
		inv(0, 1) = - det_inv * mat(0, 1);
		inv(1, 0) = - det_inv * mat(1, 0);
		inv(1, 1) = det_inv * mat(0, 0);

		return inv;
	}
};



template<typename realT, int Row, int Col>
template<int Dim>
struct Matrix<realT, Row, Col>::InverseCalculator<3, Dim> {

	static Matrix<realT, 3, 3> calculate(const Matrix<realT, Row, Col>& mat) {
		const auto det_inv = 1e0 / mat.determinant();
	
		Matrix<realT, 3, 3> inv;
		inv(0, 0) = det_inv * (mat(1, 1) * mat(2, 2) - mat(1, 2) * mat(2, 1));
		inv(1, 1) = det_inv * (mat(0, 0) * mat(2, 2) - mat(0, 2) * mat(2, 0));
		inv(2, 2) = det_inv * (mat(0, 0) * mat(1, 1) - mat(0, 1) * mat(1, 0));
	
		inv(0, 1) = det_inv * (mat(0, 2) * mat(2, 1) - mat(0, 1) * mat(2, 2));
		inv(0, 2) = det_inv * (mat(0, 1) * mat(1, 2) - mat(0, 2) * mat(1, 1));
		inv(1, 2) = det_inv * (mat(0, 2) * mat(1, 0) - mat(0, 0) * mat(1, 2));
	
		inv(1, 0) = det_inv * (mat(1, 2) * mat(2, 0) - mat(1, 0) * mat(2, 2));
		inv(2, 0) = det_inv * (mat(1, 0) * mat(2, 1) - mat(2, 0) * mat(1, 1));
		inv(2, 1) = det_inv * (mat(2, 0) * mat(0, 1) - mat(0, 0) * mat(2, 1));
	
		return inv;
	}
};

}


#endif /* COFFEE_MAKERS_MATRIX_HPP */
