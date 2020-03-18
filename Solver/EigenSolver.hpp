#ifndef COFFEE_MAKERS_EIGEN_SOLVER_HPP
#define COFFEE_MAKERS_EIGEN_SOLVER_HPP
#include<coffee-makers/Solver/EigenSolver_Base.hpp>
#include<coffee-makers/Containers/FixedMatrix.hpp>
#include<coffee-makers/Containers/FixedVector.hpp>
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


template<template<typename> typename Method, typename MatT>
class EigenSolver<Method, MatT, true> {


	using scalarT = typename Method<MatT>::ScalarType;
	using matrixT = MatT;

	const static std::size_t dimN = Method<MatT>::Dim;

	using vectorT = FixedVector<scalarT, dimN>;



public:
	EigenSolver() {
		_object = std::make_unique<Method<MatT>>();
	}
	~EigenSolver() = default;

	void solve(const matrixT& target) {
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
