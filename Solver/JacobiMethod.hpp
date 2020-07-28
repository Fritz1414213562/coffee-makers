#ifndef COFFEE_MAKERS_EIGEN_SOLVER_JACOBI_METHOD_HPP
#define COFFEE_MAKERS_EIGEN_SOLVER_JACOBI_METHOD_HPP
#include<coffee-makers/Solver/EigenSolver_Base.hpp>
#include<coffee-makers/Containers/Containers.hpp>
#include<coffee-makers/CompileTimeCalculation/matrix_traits.hpp>
#include<coffee-makers/utility/sort.hpp>
#include<iostream>
#include<array>
#include<cmath>
#include<utility>


namespace makers {


template<typename MatT>
class JacobiMethod : public EigenSolverBase<JacobiMethod, MatT> {

public:

	using MatrixType = MatT;
	using ScalarType = typename MatT::scalar_type;

private:

	using scalarT = ScalarType;

	constexpr static int _dim = MatrixType::Row_CompileTime;
	using VecT = Vector<scalarT, _dim>;

	constexpr static int max_loop = 10000000;


public:
	JacobiMethod() {
		static_assert(is_fixed<MatrixType>::value,
		"In the case that template Matrix Type is variable, the default constructor of JacobiMethod is forbidden.");
	}
	JacobiMethod(const int& dim_runtime) :
	_eigen_vectors(dim_runtime, dim_runtime),
	_eigen_values(dim_runtime, 1) {

		if (dim_runtime <= 0)
			throw std::invalid_argument("Invalid arguments: dim size is not positive");
	}


	~JacobiMethod() = default;

	void solve(const MatT& target);

	MatT eigen_vectors() const {return _eigen_vectors;}

	VecT eigen_values() const {return _eigen_values;}


	constexpr static int Dim = _dim;



private:

	std::pair<int, int> max_element(const MatT& target) const;

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
	MatT Us(target.rows(), target.cols());

	for (int idx = 0; idx < target.rows(); ++idx)
		for (int jdx = 0; jdx < target.cols(); ++jdx)
			Us(idx, jdx) = (idx == jdx) ? 1. : 0.;
	
	int loop = 0;
	for (; loop < max_loop; ++loop) {
		const std::pair<int, int>& elem_index = this->max_element(symm_mat);
		if (std::abs(symm_mat(elem_index.first, elem_index.second)) < this->absolute_tolerance()) break;

		const scalarT alpha = (
			symm_mat(elem_index.first, elem_index.first) -
			symm_mat(elem_index.second, elem_index.second)) * 0.5;
		const scalarT beta = -1. * symm_mat(elem_index.first, elem_index.second);
		const scalarT gamma = std::abs(alpha) / std::sqrt(alpha * alpha + beta * beta);
		const scalarT cos_t = std::sqrt(0.5 + gamma * 0.5);
		const scalarT sin_t = std::copysign(std::sqrt(0.5 - gamma * 0.5), alpha * beta);

		MatT U(target.rows(), target.cols());
		for (int idx = 0; idx < target.rows(); ++idx)
			for (int jdx = 0; jdx < target.cols(); ++jdx)
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
	
	VecT diagonal_values(target.rows());
	for (int idx = 0; idx < target.rows(); ++idx)
		diagonal_values[idx] = symm_mat(idx, idx);
	
	const makers::Vector<int, _dim>& sort_indices = makers::sort_Indices(diagonal_values);

	_eigen_vectors = makers::sort_Column(Us, sort_indices);
	_eigen_values = makers::sort_Row(diagonal_values, sort_indices);
}




template<typename MatT>
std::pair<int, int> JacobiMethod<MatT>::max_element(const MatT& target) const {
	scalarT max_elem = std::abs(target(0, 1));
	std::pair<int, int> retval = std::make_pair(0, 1);

	for (int idx = 0; idx < target.rows() - 1; ++idx)
		for (int jdx = idx + 1; jdx < target.cols(); ++jdx)
			if (max_elem < std::abs(target(idx, jdx))) {
				max_elem = std::abs(target(idx, jdx));
				retval = std::make_pair(idx, jdx);
			}
	return retval;
}


template<typename MatT>
typename JacobiMethod<MatT>::scalarT JacobiMethod<MatT>::max_relative_diff(const MatT& lhs, const MatT& rhs) const {
	scalarT retval = 0.0;
	for (int idx = 0; idx < lhs.rows(); ++idx) {
		const scalarT& tmp = std::abs(lhs(idx, idx) / rhs(idx, idx));
		if (retval < tmp) retval = tmp;
	}
	return retval;
}


}



#endif /* COFFEE_MAKERS_EIGEN_SOLVER_JACOBI_METHOD_HPP */
