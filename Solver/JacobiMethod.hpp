#ifndef COFFEE_MAKERS_EIGEN_SOLVER_JACOBI_METHOD_HPP
#define COFFEE_MAKERS_EIGEN_SOLVER_JACOBI_METHOD_HPP
#include<coffee-makers/Solver/EigenSolver_Base.hpp>
#include<coffee-makers/Containers/FixedMatrix.hpp>
#include<coffee-makers/Containers/FixedVector.hpp>
#include<iostream>
#include<array>
#include<cmath>
#include<utility>


namespace makers {


template<typename MatT>
class JacobiMethod : public EigenSolverBase<JacobiMethod, MatT> {

public:

	using MatrixType = MatT;
	using ScalarType = typename MatT::element_type;

private:

	using scalarT = ScalarType;

	constexpr static std::size_t _row = MatT::Row_CompileTime;
	constexpr static std::size_t _col = MatT::Col_CompileTime;

	constexpr static std::size_t _dim = _row;
	using VecT = FixedVector<scalarT, _dim>;

	constexpr static std::size_t max_loop = 10000;


public:
	JacobiMethod() = default;	
	~JacobiMethod() = default;

	void solve(const MatT& target);

	MatT eigen_vectors() const {return _eigen_vectors;}

	VecT eigen_values() const {return _eigen_values;}


	constexpr static std::size_t Dim = _dim;



private:

	std::pair<std::size_t, std::size_t> max_element(const MatT& target) const;

	scalarT max_relative_diff(const MatT& lhs, const MatT& rhs) const;


private:

	MatT _eigen_vectors;

	VecT _eigen_values;

	bool is_DebugMode = false;

};




template<typename MatT>
void JacobiMethod<MatT>::solve(const MatT& target) {

	if (not this->is_symmetric(target))
		throw std::invalid_argument("Asymmetric matrix");
	
	MatT symm_mat = target;
	MatT Us;

	for (std::size_t idx = 0; idx < _dim; ++idx)
		for (std::size_t jdx = 0; jdx < _dim; ++jdx)
			Us(idx, jdx) = (idx == jdx) ? 1. : 0.;
	
	std::size_t loop = 0;
	for (; loop < max_loop; ++loop) {
		const std::pair<std::size_t, std::size_t>& elem_index = this->max_element(symm_mat);
		if (std::abs(symm_mat(elem_index.first, elem_index.second)) < this->absolute_tolerance()) break;

		const scalarT alpha = (
			symm_mat(elem_index.first, elem_index.first) -
			symm_mat(elem_index.second, elem_index.second)) * 0.5;
		const scalarT beta = -1. * symm_mat(elem_index.first, elem_index.second);
		const scalarT gamma = std::abs(alpha) / std::sqrt(alpha * alpha + beta * beta);
		const scalarT cos_t = std::sqrt(0.5 + gamma * 0.5);
		const scalarT sin_t = std::copysign(std::sqrt(0.5 - gamma * 0.5), alpha * beta);

		MatT U;
		for (std::size_t idx = 0; idx < _dim; ++idx)
			for (std::size_t jdx = 0; jdx < _dim; ++jdx)
				U(idx, jdx) = (idx == jdx) ? 1. : 0.;

		U(elem_index.first, elem_index.first) = cos_t;
		U(elem_index.first, elem_index.second) = sin_t;
		U(elem_index.second, elem_index.first) = -sin_t;
		U(elem_index.second, elem_index.second) = cos_t;

		MatT tmp = U.transpose() * symm_mat * U;
		if (this->max_relative_diff(symm_mat, tmp) < this->relative_tolerance()) break;
		tmp(elem_index.first, elem_index.second) = 0.;
		tmp(elem_index.second, elem_index.first) = 0.;

		symm_mat = tmp;
		tmp = Us * U;
		Us = tmp;
	}

	if (is_DebugMode) {
		std::cout << "Loop: " << loop << std::endl;
	}

	if (loop == max_loop)
		throw std::logic_error("cannot solve with the tolerance");
	
	_eigen_vectors = Us;
	for (std::size_t idx = 0; idx < _dim; ++idx)
		_eigen_values[idx] = symm_mat(idx, idx);
}




template<typename MatT>
std::pair<std::size_t, std::size_t> JacobiMethod<MatT>::max_element(const MatT& target) const {
	scalarT max_elem = std::abs(target(0, 1));
	std::pair<std::size_t, std::size_t> retval = std::make_pair(0, 1);

	for (std::size_t idx = 0; idx < _dim - 1; ++idx)
		for (std::size_t jdx = idx + 1; jdx < _dim; ++jdx)
			if (max_elem < std::abs(target(idx, jdx))) {
				max_elem = std::abs(target(idx, jdx));
				retval = std::make_pair(idx, jdx);
			}
	return retval;
}


template<typename MatT>
typename JacobiMethod<MatT>::scalarT JacobiMethod<MatT>::max_relative_diff(const MatT& lhs, const MatT& rhs) const {
	scalarT retval = 0.0;
	for (std::size_t idx = 0; idx < _dim; ++idx) {
		const scalarT& tmp = std::abs(lhs(idx, idx) / rhs(idx, idx));
		if (retval < tmp) retval = tmp;
	}
	return retval;
}


}



#endif /* COFFEE_MAKERS_EIGEN_SOLVER_JACOBI_METHOD_HPP */
