#ifndef COFFEE_MAKERS_SINGULAR_VALUE_SOLVER_HPP
#define COFFEE_MAKERS_SINGULAR_VALUE_SOLVER_HPP
#include<coffee-makers/Containers/Containers.hpp>
#include<coffee-makers/CompileTimeCalculation/matrix_traits.hpp>
#include<coffee-makers/Solver/EigenSolver_Base.hpp>
#include<coffee-makers/Solver/EigenSolver.hpp>
#include<coffee-makers/Solver/JacobiMethod.hpp>
#include<coffee-makers/utility/sqrt.hpp>
#include<cmath>


namespace makers {


template<
	template<typename> typename EigenSolverMethod,
	typename MatT,
	int dimN = std::min(MatT::Row_CompileTime, MatT::Col_CompileTime),
	bool isExtended = std::is_base_of<
		EigenSolverBase<EigenSolverMethod, Matrix<typename MatT::scalar_type, dimN, dimN>>,
		EigenSolverMethod<Matrix<typename MatT::scalar_type, dimN, dimN>>
	>::value>
class SingularValueSolver {
	static_assert(isExtended, "EigenSolverMethod does not extend 'EigenSolver_Base'");
};




template<template<typename> typename EigenSolverMethod, typename MatT, int dimN>
class SingularValueSolver<EigenSolverMethod, MatT, dimN, true> {

	static_assert((dimN == std::min(MatT::Row_CompileTime, MatT::Col_CompileTime)),
	"When you use second argument, the value must be minimal value between row and column size.");


	using scalarT = typename MatT::scalar_type;
	using vectorT = Vector<scalarT, dimN>;
	using TransposedMatT = Matrix<scalarT, MatT::Col_CompileTime, MatT::Row_CompileTime>;

	using left_matrixT = Matrix<scalarT, MatT::Row_CompileTime, MatT::Row_CompileTime>;
	using right_matrixT = Matrix<scalarT, MatT::Col_CompileTime, MatT::Col_CompileTime>;

public:

	SingularValueSolver() :
		rows(MatT::Row_CompileTime),
		cols(MatT::Col_CompileTime),
		dims(dimN) {
		
		static_assert(is_fixed<MatT>::value,
		"In the case that template matrix type is variable, the use of default constructor is forbedden.");
		_eigen_solver4U = std::make_unique<EigenSolverMethod<left_matrixT>>();
		_eigen_solver4V = std::make_unique<EigenSolverMethod<right_matrixT>>();
	}

	SingularValueSolver(const int& row_runtime, const int& col_runtime) :
		rows(row_runtime),
		cols(col_runtime),
		dims(std::min(row_runtime, col_runtime)) {

		if ((row_runtime < 0) || (col_runtime < 0))
			throw std::invalid_argument("Invalid arguments: row / column size is not positive");

	//	if (dims != dimN)
	//		throw std::runtime_error(
	//			"Runtime error: Dimension size of run time is not that of compile time");

		_eigen_solver4U = std::make_unique<EigenSolverMethod<left_matrixT>>(rows);
		_eigen_solver4V = std::make_unique<EigenSolverMethod<right_matrixT>>(cols);
	}

	~SingularValueSolver() = default;


	template<typename T>
	void solve(const T& target);


	vectorT singular_values() const;

	left_matrixT left_singular_vectors() const {
		return _eigen_solver4U->eigen_vectors();
	}

	right_matrixT right_singular_vectors() const {
		return _eigen_solver4V->eigen_vectors();
	}


private:
	std::unique_ptr<EigenSolverMethod<left_matrixT>> _eigen_solver4U;
	std::unique_ptr<EigenSolverMethod<right_matrixT>> _eigen_solver4V;
	MatT _origin_matrix;

	const int rows = -1;
	const int cols = -1;
	const int dims = -1;

};


template<template<typename> typename EigenSolverMethod, typename MatT, int dimN>
template<typename T>
void SingularValueSolver<EigenSolverMethod, MatT, dimN, true>::solve(const T& target) {
	static_assert(std::is_same<T, MatT>::value,
	"The matrix type of a solver is not consistent with that of an argument");

	/*
	if row size is smaller than column one, choose AAt.
	Otherwise, choose AtA
	*/

	// perform eigendecomposition of AAt or AtA.
	_origin_matrix = target;
	_eigen_solver4U->solve(target * target.transpose());
	_eigen_solver4V->solve(target.transpose() * target);

}

template<template<typename> typename EigenSolverMethod, typename MatT, int dimN>
typename SingularValueSolver<EigenSolverMethod, MatT, dimN, true>::vectorT SingularValueSolver<EigenSolverMethod, MatT, dimN, true>::singular_values() const {
	const MatT& sigma = left_singular_vectors().transpose() * _origin_matrix * right_singular_vectors();
	vectorT retval(dims);
	for (int idim = 0; idim < dims; ++idim) {
		retval[dims - 1 - idim] = sigma(sigma.rows() - 1 - idim, sigma.cols() - 1 - idim);
	}
	return retval;

//	if constexpr (not is_fixed<MatT>::value) {
//		if (rows < cols) return squareroot(_eigen_solver4U->eigen_values());
//		else return squareroot(_eigen_solver4V->eigen_values());
//	}
//	else if constexpr (MatT::Row_CompileTime < MatT::Col_CompileTime)
//		return squareroot(_eigen_solver4U->eigen_values());
//	else if constexpr (MatT::Row_CompileTime >= MatT::Col_CompileTime)
//		return squareroot(_eigen_solver4V->eigen_values());
}

}

#endif // COFFEE_MAKERS_SINGULAR_VALUE_SOLVER_HPP
