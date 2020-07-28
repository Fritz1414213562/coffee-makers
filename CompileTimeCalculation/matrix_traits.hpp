#ifndef COFFEE_MAKERS_COMPILE_TIME_MATRIX_TRAITS_HPP
#define COFFEE_MAKERS_COMPILE_TIME_MATRIX_TRAITS_HPP
#include<type_traits>
#include<coffee-makers/CompileTimeCalculation/disjunction.hpp>
#include<coffee-makers/CompileTimeCalculation/conjunction.hpp>


namespace makers {

template<typename MatrixType>
using is_variable = value_disjunction<
	(MatrixType::Row_CompileTime < 0), (MatrixType::Col_CompileTime < 0)>;

template<typename MatrixType>
using is_fixed = value_conjunction<
	(MatrixType::Row_CompileTime > 0), (MatrixType::Col_CompileTime > 0)>;

template<typename MatrixType>
using is_square_matrix = std::bool_constant<
	(MatrixType::Row_CompileTime == MatrixType::Col_CompileTime)>;

}

#endif /* COFFEE_MAKERS_COMPILE_TIME_MATRIX_TRAITS_HPP */
