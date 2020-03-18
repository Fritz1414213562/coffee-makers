#ifndef COFFEE_MAKERS_VARIABLE_MATRIX_HPP
#define COFFEE_MAKERS_VARIABLE_MATRIX_HPP
#include<vector>
#include<array>
#include<iostream>

namespace makers {


template<typename realT>
class VariableMatrix {

public:
	using element_type = realT;

	std::size_t size() const {return element_size;}
	std::size_t rows() const {return dim_row;}
	std::size_t cols() const {return dim_col;}


public:
	~VariableMatrix() = default;
	
	VariableMatrix(const std::size_t& row, const std::size_t& col) :
		_matrix(row * col, 0),
		element_size(row * col),
		dim_row(row),
		dim_col(col) {}

	template<typename initT, class = typename std::enable_if<
		std::is_convertible<initT, realT>::value>::type>
	VariableMatrix(const std::size_t& row, const std::size_t& col, const initT& ini_val) :
		_matrix(row * col, ini_val),
		element_size(row * col),
		dim_row(row),
		dim_col(col) {}


	VariableMatrix(const VariableMatrix& input_matrix) = default;

	VariableMatrix& operator=(const VariableMatrix& rhs_matrix) = default;

	VariableMatrix& operator+=(const VariableMatrix& rhs_matrix);
	VariableMatrix& operator-=(const VariableMatrix& rhs_matrix);
	VariableMatrix& operator+=(const element_type& rhs);
	VariableMatrix& operator-=(const element_type& rhs);
	VariableMatrix& operator*=(const element_type& rhs);
	VariableMatrix& operator/=(const element_type& rhs);

	element_type at(const std::size_t& idx, const std::size_t& jdx) const;
	element_type& at(const std::size_t& idx, const std::size_t& jdx);
	element_type operator()(const std::size_t& idx, const std::size_t& jdx) const;
	element_type& operator()(const std::size_t& idx, const std::size_t& jdx);

	element_type at(const std::size_t& idx) const {return _matrix.at(idx);}
	element_type& at(const std::size_t& idx) {return _matrix.at(idx);}
	element_type operator[](const std::size_t& idx) const {return _matrix[idx];}
	element_type& operator[](const std::size_t& idx) {return _matrix[idx];}


	// 
	VariableMatrix transpose() const;

	element_type determinant() const;

	VariableMatrix inverse() const;

protected:
	std::vector<element_type> _matrix;
	const std::size_t element_size = 0;
	std::size_t dim_row = 0;
	std::size_t dim_col = 0;

};





template<typename realT>
VariableMatrix<realT>&
VariableMatrix<realT>::operator+=(const VariableMatrix<realT>& rhs_matrix) {

	if ((this->rows() != rhs_matrix.rows()) || (this->cols() != rhs_matrix.cols())) {
		std::cerr << "operation += : The shape of matrix is different" << std::endl;
		std::exit(1);
	}

	for (std::size_t idx = 0; idx < this->rows(); ++idx) {
		for (std::size_t jdx = 0; jdx < this->cols(); ++jdx) {
			(*this)(idx, jdx) += rhs_matrix(idx, jdx);
		}
	}
	return *this;
}


template<typename realT>
VariableMatrix<realT>&
VariableMatrix<realT>::operator-=(const VariableMatrix<realT>& rhs_matrix) {

	if ((this->rows() != rhs_matrix.rows()) || (this->cols() != rhs_matrix.cols())) {
		std::cerr << "operation -= : The shape of matrix is different" << std::endl;
		std::exit(1);
	}

	for (std::size_t idx = 0; idx < this->rows(); ++idx) {
		for (std::size_t jdx = 0; jdx < this->cols(); ++jdx) {
			(*this)(idx, jdx) -= rhs_matrix(idx, jdx);
		}
	}
	return *this;
}


template<typename realT>
VariableMatrix<realT>&
VariableMatrix<realT>::operator+=(const element_type& rhs) {

	for (std::size_t idx = 0; idx < this->rows(); ++idx) {
		for (std::size_t jdx = 0; jdx < this->cols(); ++jdx) {
			(*this)(idx, jdx) += rhs;
		}
	}
	return *this;
}


template<typename realT>
VariableMatrix<realT>&
VariableMatrix<realT>::operator-=(const element_type& rhs) {

	for (std::size_t idx = 0; idx < this->rows(); ++idx) {
		for (std::size_t jdx = 0; jdx < this->cols(); ++jdx) {
			(*this)(idx, jdx) -= rhs;
		}
	}
	return *this;
}



template<typename realT>
VariableMatrix<realT>&
VariableMatrix<realT>::operator*=(const element_type& rhs) {

	for (std::size_t idx = 0; idx < this->rows(); ++idx) {
		for (std::size_t jdx = 0; jdx < this->cols(); ++jdx) {
			(*this)(idx, jdx) *= rhs;
		}
	}
	return *this;
}


template<typename realT>
VariableMatrix<realT>&
VariableMatrix<realT>::operator/=(const element_type& rhs) {

	for (std::size_t idx = 0; idx < this->rows(); ++idx) {
		for (std::size_t jdx = 0; jdx < this->cols(); ++jdx) {
			(*this)(idx, jdx) /= rhs;
		}
	}
	return *this;
}


template<typename realT>
VariableMatrix<realT>
operator+(const VariableMatrix<realT>& lhs, const VariableMatrix<realT>& rhs) {
	if ((lhs.size() != rhs.size()) || (lhs.rows() != rhs.rows()) || (lhs.cols() != rhs.cols())) {
		std::cerr << "The shapes of matrices are different." << std::endl;
		std::exit(1);
	}

	VariableMatrix<realT> result(lhs.rows(), lhs.cols(), 0);

	for (std::size_t idx = 0; idx < lhs.rows(); ++idx) {
		for (std::size_t jdx = 0; jdx < lhs.cols(); ++jdx) {
			result(idx, jdx) = lhs(idx, jdx) + rhs(idx, jdx);
		}
	}

	return result;
}


template<typename realT>
VariableMatrix<realT>
operator-(const VariableMatrix<realT>& lhs, const VariableMatrix<realT>& rhs) {
	if ((lhs.size() != rhs.size()) || (lhs.rows() != rhs.rows()) || (lhs.cols() != rhs.cols())) {
		std::cerr << "The shapes of matrices are different." << std::endl;
		std::exit(1);
	}

	VariableMatrix<realT> result(lhs.rows(), lhs.cols(), 0);

	for (std::size_t idx = 0; idx < lhs.rows(); ++idx) {
		for (std::size_t jdx = 0; jdx < lhs.cols(); ++jdx) {
			result(idx, jdx) = lhs(idx, jdx) - rhs(idx, jdx);
		}
	}

	return result;
}


template<typename realT>
VariableMatrix<realT>
operator*(const VariableMatrix<realT>& lhs, const realT& rhs) {
	VariableMatrix<realT> result(lhs.rows(), lhs.cols());

	for (std::size_t idx = 0; idx < lhs.rows(); ++idx) {
		for (std::size_t jdx = 0; jdx < lhs.cols(); ++jdx) {
			result(idx, jdx) = lhs(idx, jdx) * rhs;
		}
	}

	return result;
}


template<typename realT>
VariableMatrix<realT>
operator*(const realT& lhs, const VariableMatrix<realT>& rhs) {
	VariableMatrix<realT> result(rhs.rows(), rhs.cols());

	for (std::size_t idx = 0; idx < lhs.rows(); ++idx) {
		for (std::size_t jdx = 0; jdx < lhs.cols(); ++jdx) {
			result(idx, jdx) = lhs * rhs(idx, jdx);
		}
	}

	return result;
}


template<typename realT>
VariableMatrix<realT>
operator/(const VariableMatrix<realT>& lhs, const realT& rhs) {
	VariableMatrix<realT> result(lhs.rows(), lhs.cols());

	for (std::size_t idx = 0; idx < lhs.rows(); ++idx) {
		for (std::size_t jdx = 0; jdx < lhs.cols(); ++jdx) {
			result(idx, jdx) = lhs(idx, jdx) / rhs;
		}
	}

	return result;
}


template<typename realT>
VariableMatrix<realT>
operator*(const VariableMatrix<realT>& lhs, const VariableMatrix<realT>& rhs) {
	if (lhs.cols() != rhs.rows()) {
		std::cerr << "Invalid Matrix Operation!: The column size of matrix is not consistent with the row size of right-handed side matrix." << std::endl;
		std::exit(1);
	}

	VariableMatrix<realT> result(lhs.rows(), rhs.cols(), 0);

	for (std::size_t idx = 0; idx < lhs.rows(); ++idx) {
		for (std::size_t jdx = 0; jdx < rhs.cols(); ++jdx) {
			for (std::size_t kdx = 0; kdx < lhs.cols(); ++kdx) {
				result(idx, jdx) += lhs(idx, kdx) * rhs(kdx, jdx);
			}
		}
	}

	return result;
}


template<typename realT>
typename VariableMatrix<realT>::element_type
VariableMatrix<realT>::at(const std::size_t& idx, const std::size_t& jdx) const {
	return this->_matrix.at(idx * (this->cols()) + jdx);
}


template<typename realT>
typename VariableMatrix<realT>::element_type&
VariableMatrix<realT>::at(const std::size_t& idx, const std::size_t& jdx) {
	return this->_matrix.at(idx * (this->cols()) + jdx);
}


template<typename realT>
typename VariableMatrix<realT>::element_type
VariableMatrix<realT>::operator()(const std::size_t& idx, const std::size_t& jdx) const {
	return this->_matrix[idx * (this->cols()) + jdx];
}


template<typename realT>
typename VariableMatrix<realT>::element_type&
VariableMatrix<realT>::operator()(const std::size_t& idx, const std::size_t& jdx) {
	return this->_matrix[idx * (this->cols()) + jdx];
}


template<typename realT>
VariableMatrix<realT> VariableMatrix<realT>::transpose() const {
	VariableMatrix<realT> result(this->cols(), this->rows());
	for (std::size_t idx = 0; idx < this->rows(); ++idx)
		for (std::size_t jdx = 0; jdx < this->cols(); ++jdx)
			result(jdx, idx) = (*this)(idx, jdx);
	
	return result;
}


template<typename realT>
typename VariableMatrix<realT>::element_type
VariableMatrix<realT>::determinant() const {
	if ((this->cols() != this->rows()) || (this->size() == 0)) {
		std::cerr << "The determinant of non-square matrix or empty matrix is undefined." << std::endl;
		std::exit(1);
	}
	VariableMatrix<realT> buffer_matrix = *this;

	for (std::size_t idx = 0; idx < buffer_matrix.rows(); ++idx) {
		for (std::size_t jdx = 0; jdx < buffer_matrix.cols(); ++jdx) {
			if (idx >= jdx) continue;
			else {
				const element_type& buf = buffer_matrix(jdx, idx) / buffer_matrix(idx, idx);
				for (std::size_t kdx = 0; kdx < buffer_matrix.cols(); ++kdx)
					buffer_matrix(jdx, kdx) -= buffer_matrix(idx, kdx) * buf;
			}
		}
	}

	element_type result = 1e0;

	for (std::size_t idx = 0; idx < buffer_matrix.rows(); ++idx)
		result *= buffer_matrix(idx, idx);
	
	return result;
}


template<typename realT>
inline realT determinant(const VariableMatrix<realT>& mat) {
	if ((mat.cols() != mat.rows()) || (mat.size() == 0)) {
		std::cerr << "The determinant of non-square matrix or empty matrix is undefined." << std::endl;
		std::exit(1);
	}
	VariableMatrix<realT> buffer_matrix = mat;

	for (std::size_t idx = 0; idx < buffer_matrix.rows(); ++idx) {
		for (std::size_t jdx = 0; jdx < buffer_matrix.cols(); ++jdx) {
			if (idx >= jdx) continue;
			else {
				const realT& buf = buffer_matrix(jdx, idx) / buffer_matrix(idx, idx);
				for (std::size_t kdx = 0; kdx < buffer_matrix.cols(); ++kdx)
					buffer_matrix(jdx, kdx) -= buffer_matrix(idx, kdx) * buf;
			}
		}
	}

	realT result = 1e0;

	for (std::size_t idx = 0; idx < buffer_matrix.rows(); ++idx)
		result *= buffer_matrix(idx, idx);
	
	return result;
}



template<typename realT>
inline realT determinant_of2d(const VariableMatrix<realT>& matrix_2d) {
	if ((matrix_2d.cols() != matrix_2d.rows()) || (matrix_2d.rows() == 2)) {
		std::cerr << "The determinant of non-square matrix is undefined," << std::endl;
		std::cerr << "or the matrix shape is not {2, 2}" << std::endl;
		std::exit(1);
	}

	return matrix_2d(0, 0) * matrix_2d(1, 1) - matrix_2d(1, 0) * matrix_2d(0, 1);
}


template<typename realT>
inline realT determinant_of3d(const VariableMatrix<realT>& matrix_3d) {
	if ((matrix_3d.cols() != matrix_3d.rows()) || (matrix_3d.rows() == 3)) {
		std::cerr << "The determinant of non-square matrix is undefined," << std::endl;
		std::cerr << "or the matrix shape is not {3, 3}" << std::endl;
		std::exit(1);
	}

	return 
		(matrix_3d(0, 0)) * (matrix_3d(1, 1)) * (matrix_3d(2, 2)) +
		(matrix_3d(1, 0)) * (matrix_3d(2, 1)) * (matrix_3d(0, 2)) +
		(matrix_3d(2, 0)) * (matrix_3d(0, 1)) * (matrix_3d(1, 2)) -
		(matrix_3d(0, 0)) * (matrix_3d(2, 1)) * (matrix_3d(1, 2)) -
		(matrix_3d(2, 0)) * (matrix_3d(1, 1)) * (matrix_3d(0, 2)) -
		(matrix_3d(1, 0)) * (matrix_3d(0, 1)) * (matrix_3d(2, 2));
}




template<typename realT>
VariableMatrix<realT> VariableMatrix<realT>::inverse() const {

	if ((this->rows() != this->cols()) || (this->determinant() == 0)) {
		std::cerr << "The inverse operation to non-regular matrix is undefined." << std::endl;
		std::exit(1);
	}

	VariableMatrix<realT> buffer_matrix = *this;
	VariableMatrix<realT> inverse_matrix;

	for (std::size_t idx = 0; idx < this->rows(); ++idx)
		inverse_matrix(idx, idx) = 1e0;

	for (std::size_t idx = 0; idx < this->rows(); ++idx) {
		element_type inv_mat_i = 1e0 / buffer_matrix(idx, idx);
		for (std::size_t jdx = 0; jdx < this->cols(); ++jdx) {
			buffer_matrix(idx, jdx) *= inv_mat_i;
			inverse_matrix(idx, jdx) *= inv_mat_i;
		}

		for (std::size_t jdx = 0; jdx < this->rows(); ++jdx) {
			if (idx != jdx) {
				element_type buffer = buffer_matrix(jdx, idx);
				for (std::size_t kdx = 0; kdx < this->cols(); ++kdx) {
					buffer_matrix(jdx, kdx) -= buffer_matrix(idx, kdx) * buffer;
					inverse_matrix(jdx, kdx) -= inverse_matrix(idx, kdx) * buffer;
				}
			}
		}
	}

	return inverse_matrix;
}


template<typename realT>
inline VariableMatrix<realT> inverse(const VariableMatrix<realT>& mat) {

	if ((mat.rows() != mat.cols()) || (mat.determinant() == 0)) {
		std::cerr << "The inverse operation to non-regular matrix is undefined." << std::endl;
		std::exit(1);
	}

	VariableMatrix<realT> buffer_matrix = mat;
	VariableMatrix<realT> inverse_matrix;

	for (std::size_t idx = 0; idx < mat.rows(); ++idx)
		inverse_matrix(idx, idx) = 1e0;

	for (std::size_t idx = 0; idx < mat.rows(); ++idx) {
		realT inv_mat_i = 1e0 / buffer_matrix(idx, idx);
		for (std::size_t jdx = 0; jdx < mat.cols(); ++jdx) {
			buffer_matrix(idx, jdx) *= inv_mat_i;
			inverse_matrix(idx, jdx) *= inv_mat_i;
		}

		for (std::size_t jdx = 0; jdx < mat.rows(); ++jdx) {
			if (idx != jdx) {
				realT buffer = buffer_matrix(jdx, idx);
				for (std::size_t kdx = 0; kdx < mat.cols(); ++kdx) {
					buffer_matrix(jdx, kdx) -= buffer_matrix(idx, kdx) * buffer;
					inverse_matrix(jdx, kdx) -= inverse_matrix(idx, kdx) * buffer;
				}
			}
		}
	}

	return inverse_matrix;
}


template<typename realT>
inline VariableMatrix<realT> inverse_of2d(const VariableMatrix<realT>& matrix_2d) {
	const auto det_inv = 1e0 / determinant_of2d(matrix_2d);

	VariableMatrix<realT> inv(2, 2);
	inv(0, 0) = det_inv * matrix_2d(1, 1);
	inv(0, 1) = - det_inv * matrix_2d(0, 1);
	inv(1, 0) = - det_inv * matrix_2d(1, 0);
	inv(1, 1) = det_inv * matrix_2d(0, 0);

	return inv;
}


template<typename realT>
inline VariableMatrix<realT> inverse_of3d(const VariableMatrix<realT>& matrix_3d) {
	const auto det_inv = 1e0 / determinant_of3d(matrix_3d);

	VariableMatrix<realT> inv(3, 3);
	inv(0, 0) = det_inv * (matrix_3d(1, 1) * matrix_3d(2, 2) - matrix_3d(1, 2) * matrix_3d(2, 1));
	inv(1, 1) = det_inv * (matrix_3d(0, 0) * matrix_3d(2, 2) - matrix_3d(0, 2) * matrix_3d(2, 0));
	inv(2, 2) = det_inv * (matrix_3d(0, 0) * matrix_3d(1, 1) - matrix_3d(0, 1) * matrix_3d(1, 0));

	inv(0, 1) = det_inv * (matrix_3d(0, 2) * matrix_3d(2, 1) - matrix_3d(0, 1) * matrix_3d(2, 2));
	inv(0, 2) = det_inv * (matrix_3d(0, 1) * matrix_3d(1, 2) - matrix_3d(0, 2) * matrix_3d(1, 1));
	inv(1, 2) = det_inv * (matrix_3d(0, 2) * matrix_3d(1, 0) - matrix_3d(0, 0) * matrix_3d(1, 2));

	inv(1, 0) = det_inv * (matrix_3d(1, 2) * matrix_3d(2, 0) - matrix_3d(1, 0) * matrix_3d(2, 2));
	inv(2, 0) = det_inv * (matrix_3d(1, 0) * matrix_3d(2, 1) - matrix_3d(2, 0) * matrix_3d(1, 1));
	inv(2, 1) = det_inv * (matrix_3d(2, 0) * matrix_3d(0, 1) - matrix_3d(0, 0) * matrix_3d(2, 1));

	return inv;
}




}


#endif /* COFFEE_MAKERS_VARIABLE_MATRIX_HPP */
