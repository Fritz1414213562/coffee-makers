#ifndef COFFEE_MAKERS_EIGEN_SOLVER_BASE_HPP
#define COFFEE_MAKERS_EIGEN_SOLVER_BASE_HPP
#include<coffee-makers/Containers/FixedMatrix.hpp>
#include<coffee-makers/Containers/FixedVector.hpp>
#include<cmath>


namespace makers {


template<template<typename> typename Method, typename MatT>
class EigenSolverBase {


	using matrixT = MatT;
	using scalarT = typename matrixT::element_type;


	static_assert(matrixT::Row_CompileTime == matrixT::Col_CompileTime, "Eigen value decomposition is undefined for a non-regular matrix.");
	
	constexpr static std::size_t dimN = matrixT::Row_CompileTime;
	
	using vectorT = FixedVector<scalarT, dimN>;
	

public:

	void solve(const matrixT& target) {static_cast<Method<MatT>&>(this)->solve(target);}

	matrixT eigen_vectors() const {return static_cast<Method<MatT>&>(this)->eigen_vectors();}

	vectorT eigen_values() const {return static_cast<Method<MatT>&>(this)->eigen_values();}

	void set_DebugMode() {static_cast<Method<MatT>&>(this)->set_DebugMode();}

protected:

	bool is_symmetric(const matrixT& mat) const {

		bool retval = true;

		for (std::size_t idx = 0; idx < dimN - 1; ++idx)
			for (std::size_t jdx = idx + 1; jdx < dimN; ++jdx)
				retval = retval && (std::abs(mat(idx, jdx) / mat(jdx, idx) - 1.) <= relative_tolerance());
		return retval;
	}


	static scalarT relative_tolerance();

	static scalarT absolute_tolerance();


private:

	template<typename Type, typename DummyType = void>
	struct Inner_RelativeTolerance;

	template<typename Type, typename DummyType = void>
	struct Inner_AbsoluteTolerance;

};



template<template<typename> typename Method, typename MatT>
template<typename Type>
struct EigenSolverBase<Method, MatT>::Inner_RelativeTolerance<int, Type> {


	using type = int;
	constexpr static int value = 0;

};


template<template<typename> typename Method, typename MatT>
template<typename Type>
struct EigenSolverBase<Method, MatT>::Inner_RelativeTolerance<float, Type> {

	using type = float;
	constexpr static float value = 1e-5;

};



template<template<typename> typename Method, typename MatT>
template<typename Type>
struct EigenSolverBase<Method, MatT>::Inner_RelativeTolerance<double, Type> {

	using type = double;
	constexpr static double value = 1e-13;

};



template<template<typename> typename Method, typename MatT >
template<typename Type>
struct EigenSolverBase<Method, MatT>::Inner_AbsoluteTolerance<int, Type> {

	using type = int;
	constexpr static int value = 0;

};



template<template<typename> typename Method, typename MatT>
template<typename Type>
struct EigenSolverBase<Method, MatT>::Inner_AbsoluteTolerance<float, Type> {

	using type = float;
	constexpr static float value = 1e-5;

};



template<template<typename> typename Method, typename MatT>
template<typename Type>
struct EigenSolverBase<Method, MatT>::Inner_AbsoluteTolerance<double, Type> {

	using type = double;
	constexpr static double value = 1e-13;

};



template<template<typename> typename Method, typename MatT>
typename EigenSolverBase<Method, MatT>::scalarT EigenSolverBase<Method, MatT>::relative_tolerance() {return Inner_RelativeTolerance<scalarT>::value;}

template<template<typename> typename Method, typename MatT>
typename EigenSolverBase<Method, MatT>::scalarT EigenSolverBase<Method, MatT>::absolute_tolerance() {return Inner_AbsoluteTolerance<scalarT>::value;} 
}


#endif /* COFFEE_MAKERS_EIGEN_SOLVER_BASE_HPP */
