#ifndef COFFEE_MAKERS_VARIABLE_VECTOR_HPP
#define COFFEE_MAKERS_VARIABLE_VECTOR_HPP
//#include<coffee-makers/CompileTimeCalculation/conjunction.hpp>
#include<coffee-makers/Containers/VariableMatrix.hpp>
#include<cmath>


namespace makers {


// forbid upcast
template<typename realT>
class VariableVector : protected VariableMatrix<realT> {


public:

	using element_type = realT;

	std::size_t size() const {return this->element_size;}
	std::size_t rows() const {return this->dim_row;}
	std::size_t cols() const {return this->dim_col;}

public:
	~VariableVector() = default;
	VariableVector(const std::size_t& dim) : VariableMatrix<realT>(dim, 1) {}

	template<typename initT, class = typename std::enable_if<
		std::is_convertible<initT, realT>::value>::type>
	VariableVector(const std::size_t& dim, const initT& ini_val) :
		VariableMatrix<realT>(dim, 1, ini_val) {}
	
	VariableVector(const VariableVector& input_vector) = default;

	VariableVector(const VariableMatrix<realT>& input_matrix) :
		VariableMatrix<realT>(input_matrix.rows(), input_matrix.cols()) {
	//	_matrix(input_matrix.size()),
	//	element_size(input_matrix.size()),
	//	dim_row(input_matrix.size()),
	//	dim_col(1) {
	
		if ((input_matrix.rows() != 1) && (input_matrix.cols() != 1)) {
			std::cerr << "The shape of right-hand side matrix is not vector-shape" << std::endl;
			std::exit(1);
		}
		for (std::size_t idx = 0; idx < input_matrix.size(); ++idx)
			this->_matrix[idx] = input_matrix[idx];
	}


	VariableVector& operator=(const VariableVector& rhs_vector) = default;
	VariableVector& operator=(const VariableMatrix<realT>& rhs_matrix);

	VariableVector& operator+=(const VariableVector& rhs_vector);
	VariableVector& operator-=(const VariableVector& rhs_vector);
	VariableVector& operator+=(const element_type& rhs);
	VariableVector& operator-=(const element_type& rhs);
	VariableVector& operator*=(const element_type& rhs);
	VariableVector& operator/=(const element_type& rhs);

	// operation
	element_type operator()(const std::size_t& idx) const {return this->_matrix[idx];}
	element_type& operator()(const std::size_t& idx) {return this->_matrix[idx];}
	element_type at(const std::size_t& idx) const {return this->_matrix[idx];}
	element_type& at(const std::size_t& idx) {return this->_matrix[idx];}


	// delete the member functions derived from VariableMatrix
	element_type at(const std::size_t& idx, const std::size_t& jdx) const = delete;
	element_type& at(const std::size_t& idx, const std::size_t& jdx) = delete;
	element_type operator()(const std::size_t& idx, const std::size_t& jdx) const = delete;
	element_type& operator()(const std::size_t& idx, const std::size_t& jdx) = delete;

	element_type determinant() const = delete;
	VariableMatrix<realT> inverse() const = delete;


	// defined at VariableMatrix
//	element_type at(const std::size_t& idx) const {return _matrix.at(idx);}
//	element_type& at(const std::size_t& idx) {return _matrix.at(idx);}
//	element_type operator[](const std::size_t& idx) const {return _matrix[idx];}
//	element_type& operator[](const std::size_t& idx) {return _matrix[idx];}

//	VariableMatrix transpose() const;
	
};




//template<typename realT>
//VariableVector<realT>::VariableVector<realT>(const VariableMatrix<realT>& input_matrix) :
//	_matrix(rhs_matrix.size()),
//	element_size(rhs_matrix.size()),
//	dim_row(rhs_matrix.size()),
//	dim_col(1) {
//
//	if ((rhs_matrix.rows() != 1) && (rhs_matrix.cols() != 1)) {
//		std::cerr << "The shape of right-hand side matrix is not vector-shape" << std::endl;
//		std::exit(1);
//	}
//	for (std::size_t idx = 0; idx < rhs_matrix.size(); ++idx)
//		_matrix[idx] = rhs_matrix[idx];
//}


template<typename realT>
VariableVector<realT>&
VariableVector<realT>::operator=(const VariableMatrix<realT>& rhs_matrix) {

	if ((rhs_matrix.size() != this->size()) || ((rhs_matrix.rows() != 1) && (rhs_matrix.cols() != 1))) {
		std::cerr << "operation = : The shape of matrix is not vector-shape" << std::endl;
		std::exit(1);
	}
	for (std::size_t idx = 0; idx < rhs_matrix.size(); ++idx)
		(*this)(idx) = rhs_matrix(idx);
	
	return *this;
}



template<typename realT>
VariableVector<realT>&
VariableVector<realT>::operator+=(const VariableVector<realT>& rhs_vector) {

	if (this->size() != rhs_vector->size()) {
		std::cerr << "operation += : The shape of vector is different" << std::endl;
		std::exit(1);
	}

	for (std::size_t idx = 0; idx < rhs_vector.size(); ++idx)
		(*this)(idx) += rhs_vector(idx);
	
	return *this;
}


template<typename realT>
VariableVector<realT>&
VariableVector<realT>::operator-=(const VariableVector<realT>& rhs_vector) {

	if (this->size() != rhs_vector->size()) {
		std::cerr << "operation -= : The shape of vector is different" << std::endl;
		std::exit(1);
	}

	for (std::size_t idx = 0; idx < rhs_vector.size(); ++idx)
		(*this)(idx) -= rhs_vector(idx);
	
	return *this;
}


template<typename realT>
VariableVector<realT>&
VariableVector<realT>::operator+=(const element_type& rhs) {

	for (std::size_t idx = 0; idx < this->size(); ++idx)
		(*this)(idx) += rhs;
	
	return *this;
}


template<typename realT>
VariableVector<realT>&
VariableVector<realT>::operator-=(const element_type& rhs) {

	for (std::size_t idx = 0; idx < this->size(); ++idx)
		(*this)(idx) -= rhs;
	
	return *this;
}


template<typename realT>
VariableVector<realT>&
VariableVector<realT>::operator*=(const element_type& rhs) {

	for (std::size_t idx = 0; idx < this->size(); ++idx)
		(*this)(idx) *= rhs;
	
	return *this;
}


template<typename realT>
VariableVector<realT>&
VariableVector<realT>::operator/=(const element_type& rhs) {

	for (std::size_t idx = 0; idx < this->size(); ++idx)
		(*this)(idx) /= rhs;
	
	return *this;
}


template<typename realT>
VariableVector<realT>
operator+(const VariableVector<realT>& lhs, const VariableVector<realT>& rhs) {
	if (lhs.size() != rhs.size()) {
		std::cerr << "The vector size is not consistent." << std::endl;
		std::exit(1);
	}

	VariableVector<realT> result(lhs.size());

	for (std::size_t idx = 0; idx < lhs.size(); ++idx)
		result[idx] = lhs(idx) + rhs(idx);
	
	return result;
}


template<typename realT>
VariableVector<realT>
operator-(const VariableVector<realT>& lhs, const VariableVector<realT>& rhs) {
	if (lhs.size() != rhs.size()) {
		std::cerr << "The vector size is not consistent." << std::endl;
		std::exit(1);
	}

	VariableVector<realT> result(lhs.size());

	for (std::size_t idx = 0; idx < lhs.size(); ++idx)
		result[idx] = lhs(idx) - rhs(idx);
	
	return result;
}


template<typename realT>
VariableVector<realT>
operator*(const VariableVector<realT>& lhs, const realT& rhs) {
	VariableVector<realT> result(lhs.size());

	for (std::size_t idx = 0; idx < lhs.size(); ++idx)
		result(idx) = lhs(idx) * rhs;
	
	return result;
}


template<typename realT>
VariableVector<realT>
operator*(const realT& lhs, const VariableVector<realT>& rhs) {
	VariableVector<realT> result(rhs.size());

	for (std::size_t idx = 0; idx < rhs.size(); ++idx)
		result(idx) = lhs * rhs(idx);

	return result;
}


template<typename realT>
VariableVector<realT>
operator/(const VariableVector<realT>& lhs, const realT& rhs) {
	VariableVector<realT> result(lhs.size());

	for (std::size_t idx = 0; idx < lhs.size(); ++idx)
		result(idx) = lhs(idx) / rhs;
	
	return result;
}


template<typename realT>
VariableVector<realT>
operator*(const VariableMatrix<realT>& lhs_matrix, const VariableVector<realT>& rhs_vector) {
	if (lhs_matrix.cols() != rhs_vector.size()) {
		std::cerr << "Error: Invalid Matrix Opration!: The column size of matrix is not consistent with the row size of vector." << std::endl;
		std::exit(1);
	}

	VariableVector<realT> result(lhs_matrix.rows(), 0);

	for (std::size_t idx = 0; idx < lhs_matrix.rows(); ++idx)
		for (std::size_t jdx = 0; jdx < lhs_matrix.cols(); ++jdx)
			result(idx) += lhs_matrix(idx, jdx) * rhs_vector(jdx);
	
	return result;
}


}


#endif /* COFFEE_MAKERS_VARIABLE_VECTOR_HPP */
