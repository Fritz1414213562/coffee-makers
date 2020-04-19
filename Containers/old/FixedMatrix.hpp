#ifndef COFFEE_MAKERS_FIXED_MATRIX_HPP
#define COFFEE_MAKERS_FIXED_MATRIX_HPP
//#include<coffee-makers/CompileTimeCalculation/conjunction.hpp>
#include<vector>
#include<array>

namespace makers {


template<typename realT, std::size_t row, std::size_t col>
class FixedMatrix {

public:
	using element_type = realT;

	std::size_t size() const {return element_size;}
	std::size_t rows() const {return Row_CompileTime;}
	std::size_t cols() const {return Col_CompileTime;}

	constexpr static std::size_t element_size = row * col;
	constexpr static std::size_t Row_CompileTime = row;
	constexpr static std::size_t Col_CompileTime = col;


public:
	FixedMatrix() : _matrix{{}} {}
	~FixedMatrix() = default;

//	template<typename... Args, class = typename std::enable_if<
//		(sizeof...(Args) == element_size) &&
//		makers::Conjunction<std::is_convertible<Args, realT>...>::value>::type>
	template<typename... Args, class = typename std::enable_if<
		(sizeof...(Args) == element_size) &&
		std::conjunction<std::is_convertible<Args, realT>...>::value>::type>
	FixedMatrix(const Args&... args) : _matrix{{static_cast<realT>(args)...}} {}

	FixedMatrix(const FixedMatrix& input_matrix) = default;
	FixedMatrix& operator=(const FixedMatrix& rhs_matrix) = default;

	FixedMatrix& operator+=(const FixedMatrix& rhs_matrix);
	FixedMatrix& operator-=(const FixedMatrix& rhs_matrix);
	FixedMatrix& operator+=(const element_type& rhs);
	FixedMatrix& operator-=(const element_type& rhs);
	FixedMatrix& operator*=(const element_type& rhs);
	FixedMatrix& operator/=(const element_type& rhs);

	element_type at(const std::size_t& idx, const std::size_t& jdx) const;
	element_type& at(const std::size_t& idx, const std::size_t& jdx);
	element_type operator()(const std::size_t& idx, const std::size_t& jdx) const;
	element_type& operator()(const std::size_t& idx, const std::size_t& jdx);

	element_type at(const std::size_t& idx) const {return _matrix.at(idx);}
	element_type& at(const std::size_t& idx) {return _matrix.at(idx);}
	element_type operator[](const std::size_t& idx) const {return _matrix[idx];}
	element_type& operator[](const std::size_t& idx) {return _matrix[idx];}


	// 
	FixedMatrix<realT, col, row> transpose() const;

	element_type determinant() const;

	FixedMatrix inverse() const;


private:
	std::array<element_type, element_size> _matrix;

};



//template<typename realT>
//class FixedMatrix<realT, 2, 2> {
//	using element_type = realT;
//
//	element_type determinant() const;
//
//	FixedMatrix inverse() const;
//};
//
//
//template<typename realT>
//class FixedMatrix<realT, 3, 3> {
//	using element_type = realT;
//
//	element_type determinant() const;
//
//	FixedMatrix inverse() const;
//};
//




template<typename realT, std::size_t row, std::size_t col>
FixedMatrix<realT, row, col>&
FixedMatrix<realT, row, col>::operator+=(const FixedMatrix<realT, row, col>& rhs_matrix) {

	for (std::size_t idx = 0; idx < row; ++idx) {
		for (std::size_t jdx = 0; jdx < col; ++jdx) {
			(*this)(idx, jdx) += rhs_matrix(idx, jdx);
		}
	}
	return *this;
}


template<typename realT, std::size_t row, std::size_t col>
FixedMatrix<realT, row, col>&
FixedMatrix<realT, row, col>::operator-=(const FixedMatrix<realT, row, col>& rhs_matrix) {

	for (std::size_t idx = 0; idx < row; ++idx) {
		for (std::size_t jdx = 0; jdx < col; ++jdx) {
			(*this)(idx, jdx) -= rhs_matrix(idx, jdx);
		}
	}
	return *this;
}


template<typename realT, std::size_t row, std::size_t col>
FixedMatrix<realT, row, col>&
FixedMatrix<realT, row, col>::operator+=(const element_type& rhs) {

	for (std::size_t idx = 0; idx < row; ++idx) {
		for (std::size_t jdx = 0; jdx < col; ++jdx) {
			(*this)(idx, jdx) += rhs;
		}
	}
	return *this;
}


template<typename realT, std::size_t row, std::size_t col>
FixedMatrix<realT, row, col>&
FixedMatrix<realT, row, col>::operator-=(const element_type& rhs) {

	for (std::size_t idx = 0; idx < row; ++idx) {
		for (std::size_t jdx = 0; jdx < col; ++jdx) {
			(*this)(idx, jdx) -= rhs;
		}
	}
	return *this;
}



template<typename realT, std::size_t row, std::size_t col>
FixedMatrix<realT, row, col>&
FixedMatrix<realT, row, col>::operator*=(const element_type& rhs) {

	for (std::size_t idx = 0; idx < row; ++idx) {
		for (std::size_t jdx = 0; jdx < col; ++jdx) {
			(*this)(idx, jdx) *= rhs;
		}
	}
	return *this;
}


template<typename realT, std::size_t row, std::size_t col>
FixedMatrix<realT, row, col>&
FixedMatrix<realT, row, col>::operator/=(const element_type& rhs) {

	for (std::size_t idx = 0; idx < row; ++idx) {
		for (std::size_t jdx = 0; jdx < col; ++jdx) {
			(*this)(idx, jdx) /= rhs;
		}
	}
	return *this;
}


template<typename realT, std::size_t row, std::size_t col>
FixedMatrix<realT, row, col>
operator+(const FixedMatrix<realT, row, col>& lhs, const FixedMatrix<realT, row, col>& rhs) {
	FixedMatrix<realT, row, col> result;

	for (std::size_t idx = 0; idx < row; ++idx) {
		for (std::size_t jdx = 0; jdx < col; ++jdx) {
			result(idx, jdx) = lhs(idx, jdx) + rhs(idx, jdx);
		}
	}

	return result;
}


template<typename realT, std::size_t row, std::size_t col>
FixedMatrix<realT, row, col>
operator-(const FixedMatrix<realT, row, col>& lhs, const FixedMatrix<realT, row, col>& rhs) {
	FixedMatrix<realT, row, col> result;

	for (std::size_t idx = 0; idx < row; ++idx) {
		for (std::size_t jdx = 0; jdx < col; ++jdx) {
			result(idx, jdx) = lhs(idx, jdx) - rhs(idx, jdx);
		}
	}

	return result;
}


template<typename realT, std::size_t row, std::size_t col>
FixedMatrix<realT, row, col>
operator*(const FixedMatrix<realT, row, col>& lhs, const realT& rhs) {
	FixedMatrix<realT, row, col> result;

	for (std::size_t idx = 0; idx < row; ++idx) {
		for (std::size_t jdx = 0; jdx < col; ++jdx) {
			result(idx, jdx) = lhs(idx, jdx) * rhs;
		}
	}

	return result;
}


template<typename realT, std::size_t row, std::size_t col>
FixedMatrix<realT, row, col>
operator*(const realT& lhs, const FixedMatrix<realT, row, col>& rhs) {
	FixedMatrix<realT, row, col> result;

	for (std::size_t idx = 0; idx < row; ++idx) {
		for (std::size_t jdx = 0; jdx < col; ++jdx) {
			result(idx, jdx) = lhs * rhs(idx, jdx);
		}
	}

	return result;
}


template<typename realT, std::size_t row, std::size_t col>
FixedMatrix<realT, row, col>
operator/(const FixedMatrix<realT, row, col>& lhs, const realT& rhs) {
	FixedMatrix<realT, row, col> result;

	for (std::size_t idx = 0; idx < row; ++idx) {
		for (std::size_t jdx = 0; jdx < col; ++jdx) {
			result(idx, jdx) = lhs(idx, jdx) / rhs;
		}
	}

	return result;
}


template<typename realT, std::size_t L, std::size_t M, std::size_t N>
FixedMatrix<realT, L, N>
operator*(const FixedMatrix<realT, L, M>& lhs, const FixedMatrix<realT, M, N>& rhs) {
	FixedMatrix<realT, L, N> result;

	for (std::size_t idx = 0; idx < L; ++idx) {
		for (std::size_t jdx = 0; jdx < N; ++jdx) {
			for (std::size_t kdx = 0; kdx < M; ++kdx) {
				result(idx, jdx) += lhs(idx, kdx) * rhs(kdx, jdx);
			}
		}
	}

	return result;
}


template<typename realT, std::size_t row, std::size_t col>
typename FixedMatrix<realT, row, col>::element_type
FixedMatrix<realT, row, col>::at(const std::size_t& idx, const std::size_t& jdx) const {
	return this->_matrix.at(idx * col + jdx);
}


template<typename realT, std::size_t row, std::size_t col>
typename FixedMatrix<realT, row, col>::element_type&
FixedMatrix<realT, row, col>::at(const std::size_t& idx, const std::size_t& jdx) {
	return this->_matrix.at(idx * col + jdx);
}


template<typename realT, std::size_t row, std::size_t col>
typename FixedMatrix<realT, row, col>::element_type
FixedMatrix<realT, row, col>::operator()(const std::size_t& idx, const std::size_t& jdx) const {
	return this->_matrix[idx * col + jdx];
}


template<typename realT, std::size_t row, std::size_t col>
typename FixedMatrix<realT, row, col>::element_type&
FixedMatrix<realT, row, col>::operator()(const std::size_t& idx, const std::size_t& jdx) {
	return this->_matrix[idx * col + jdx];
}


template<typename realT, std::size_t row, std::size_t col>
FixedMatrix<realT, col, row> FixedMatrix<realT, row, col>::transpose() const {
//FixedMatrix<realT, row, col> FixedMatrix<realT, row, col>::transpose() const {
	FixedMatrix<realT, col, row> result;
	for (std::size_t idx = 0; idx < row; ++idx)
		for (std::size_t jdx = 0; jdx < col; ++jdx)
			result(jdx, idx) = (*this)(idx, jdx);
	
	return result;
}


template<typename realT, std::size_t row, std::size_t col>
typename FixedMatrix<realT, row, col>::element_type
FixedMatrix<realT, row, col>::determinant() const {
	FixedMatrix<realT, row, col> buffer_matrix = *this;

	for (std::size_t idx = 0; idx < row; ++idx) {
		for (std::size_t jdx = 0; jdx < row; ++jdx) {
			if (idx >= jdx) continue;
			else {
				const element_type& buf = buffer_matrix(jdx, idx) / buffer_matrix(idx, idx);
				for (std::size_t kdx = 0; kdx < col; ++kdx)
					buffer_matrix(jdx, kdx) -= buffer_matrix(idx, kdx) * buf;
			}
		}
	}

	element_type result = 1e0;

	for (std::size_t idx = 0; idx < row; ++idx)
		result *= buffer_matrix(idx, idx);
	
	return result;
}


template<typename realT, std::size_t row, std::size_t col>
inline realT determinant(const FixedMatrix<realT, row, col>& mat) {
	FixedMatrix<realT, row, col> buffer_matrix = mat;

	for (std::size_t idx = 0; idx < row; ++idx) {
		for (std::size_t jdx = 0; jdx < row; ++jdx) {
			if (idx >= jdx) continue;
			else {
				const realT& buf = buffer_matrix(jdx, idx) / buffer_matrix(idx, idx);
				for (std::size_t kdx = 0; kdx < col; ++kdx)
					buffer_matrix(jdx, kdx) -= buffer_matrix(idx, kdx) * buf;
			}
		}
	}

	realT result = 1e0;

	for (std::size_t idx = 0; idx < row; ++idx)
		result *= buffer_matrix(idx, idx);
	
	return result;
}


template<typename realT>
inline realT determinant(const FixedMatrix<realT, 2, 2>& matrix_2d) {
//FixedMatrix<realT>::determinant() const {

	return (matrix_2d(0, 0)) * (matrix_2d(1, 1)) - (matrix_2d(0, 1)) * (matrix_2d(1, 0));
}


template<typename realT>
inline realT determinant(const FixedMatrix<realT, 3, 3>& matrix_3d) {

	return 
		(matrix_3d(0, 0)) * (matrix_3d(1, 1)) * (matrix_3d(2, 2)) +
		(matrix_3d(1, 0)) * (matrix_3d(2, 1)) * (matrix_3d(0, 2)) +
		(matrix_3d(2, 0)) * (matrix_3d(0, 1)) * (matrix_3d(1, 2)) -
		(matrix_3d(0, 0)) * (matrix_3d(2, 1)) * (matrix_3d(1, 2)) -
		(matrix_3d(2, 0)) * (matrix_3d(1, 1)) * (matrix_3d(0, 2)) -
		(matrix_3d(1, 0)) * (matrix_3d(0, 1)) * (matrix_3d(2, 2));
}


template<typename realT, std::size_t row, std::size_t col>
FixedMatrix<realT, row, col> FixedMatrix<realT, row, col>::inverse() const {
	static_assert(row == col, "Invalid Matrix Operation: The matrix is not square");
	FixedMatrix<realT, row, col> buffer_matrix = *this;
	FixedMatrix<realT, row, col> inverse_matrix;

	for (std::size_t idx = 0; idx < row; ++idx)
		inverse_matrix(idx, idx) = 1e0;

	for (std::size_t idx = 0; idx < row; ++idx) {
		realT inv_mat_i = 1e0 / buffer_matrix(idx, idx);
		for (std::size_t jdx = 0; jdx < col; ++jdx) {
			buffer_matrix(idx, jdx) *= inv_mat_i;
			inverse_matrix(idx, jdx) *= inv_mat_i;
		}

		for (std::size_t jdx = 0; jdx < col; ++jdx) {
			if (idx != jdx) {
				realT buffer = buffer_matrix(jdx, idx);
				for (std::size_t kdx = 0; kdx < col; ++kdx) {
					buffer_matrix(jdx, kdx) -= buffer_matrix(idx, kdx) * buffer;
					inverse_matrix(jdx, kdx) -= inverse_matrix(idx, kdx) * buffer;
				}
			}
		}
	}

	return inverse_matrix;
}


template<typename realT, std::size_t row, std::size_t col>
inline FixedMatrix<realT, row, col> inverse(const FixedMatrix<realT, row, col>& mat) {
	static_assert(row == col, "Invalid Matrix Operation: The matrix is not square");
	FixedMatrix<realT, row, col> buffer_matrix = mat;
	FixedMatrix<realT, row, col> inverse_matrix;

	for (std::size_t idx = 0; idx < row; ++idx)
		inverse_matrix(idx, idx) = 1e0;

	for (std::size_t idx = 0; idx < row; ++idx) {
		realT inv_mat_i = 1e0 / buffer_matrix(idx, idx);
		for (std::size_t jdx = 0; jdx < col; ++jdx) {
			buffer_matrix(idx, jdx) *= inv_mat_i;
			inverse_matrix(idx, jdx) *= inv_mat_i;
		}

		for (std::size_t jdx = 0; jdx < col; ++jdx) {
			if (idx != jdx) {
				realT buffer = buffer_matrix(jdx, idx);
				for (std::size_t kdx = 0; kdx < col; ++kdx) {
					buffer_matrix(jdx, kdx) -= buffer_matrix(idx, kdx) * buffer;
					inverse_matrix(jdx, kdx) -= inverse_matrix(idx, kdx) * buffer;
				}
			}
		}
	}

	return inverse_matrix;
}


template<typename realT>
inline FixedMatrix<realT, 2, 2> inverse(const FixedMatrix<realT, 2, 2>& matrix_2d) {
	const auto det_inv = 1e0 / determinant(matrix_2d);

	FixedMatrix<realT, 2, 2> inv;
	inv(0, 0) = det_inv * matrix_2d(1, 1);
	inv(0, 1) = - det_inv * matrix_2d(0, 1);
	inv(1, 0) = - det_inv * matrix_2d(1, 0);
	inv(1, 1) = det_inv * matrix_2d(0, 0);

	return inv;
}


template<typename realT>
inline FixedMatrix<realT, 3, 3> inverse(const FixedMatrix<realT, 3, 3>& matrix_3d) {
	const auto det_inv = 1e0 / determinant(matrix_3d);

	FixedMatrix<realT, 3, 3> inv;
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


#endif /* COFFEE_MAKERS_FIXED_MATRIX_HPP */
