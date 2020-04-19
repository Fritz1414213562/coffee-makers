#ifndef COFFEE_MAKERS_EIGEN_SOLVER_HPP
#define COFFEE_MAKERS_EIGEN_SOLVER_HPP
#include<coffee-makers/Solver/EigenSolver_Base.hpp>
#include<coffee-makers/Containers/Containers.hpp>
#include<coffee-makers/CompileTimeCalculation/CompileTimeCalculation.hpp>
#include<coffee-makers/CompileTimeCalculation/matrix_traits.hpp>
#include<type_traits>
#include<typeinfo>
#include<vector>
#include<memory>

namespace makers {



template<
	template<typename> typename Method, typename MatT,
		bool isExtended = std::is_base_of<EigenSolverBase<Method, MatT>, Method<MatT>>::value>
class EigenSolver {
	static_assert(isExtended, "EigenSolver_Base is not extended.");
};

//template<typename MatT>
//using is_variable = value_disjunction<
//	(MatT::Row_CompileTime < 0), (MatT::Col_CompileTime < 0)>;
//
template<template<typename> typename Method, typename MatT>
class EigenSolver<Method, MatT, true> {

	using scalarT = typename Method<MatT>::ScalarType;
	using matrixT = MatT;

	static_assert(
	value_disjunction<
		matrixT::Row_CompileTime == matrixT::Col_CompileTime,
		is_variable<matrixT>::value>::value,
	"Eigen value decomposition is undefined for a non-regular matrix.");

	const static int dimN = Method<MatT>::Dim;

	using vectorT = Vector<scalarT, dimN>;



public:
	EigenSolver() {
		_object = std::make_unique<Method<MatT>>();
	}
	EigenSolver(const int& dim_rumtime) {
		_object = std::make_unique<Method<MatT>>(dim_rumtime);
	}
	~EigenSolver() = default;

	template<typename T>
	void solve(const T& target) {
		static_assert(std::is_same<T, matrixT>::value,
		"The Matrix Type of a Solver is not different from that of a argument");

		if constexpr (is_variable<T>::value)
			if (target.rows() != target.cols())
				throw std::invalid_argument(
				"Invalid operation: Eigen value decomposition is undefined for a non-regular matrix");
		_object->solve(target);
	}

	matrixT eigen_vectors() const {
		return _object->eigen_vectors();
	}

	vectorT eigen_values() const {
		return _object->eigen_values();
	}


	void set_DebugMode() {_object->set_DebugMode();}

private:
	std::unique_ptr<Method<MatT>> _object;

};


}


#endif /* COFFEE_MAKERS_EIGEN_SOLVER_HPP */
