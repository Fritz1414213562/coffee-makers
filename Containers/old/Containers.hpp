#ifndef COFFEE_MAKERS_CONTAINERS_HPP
#define COFFEE_MAKERS_CONTAINERS_HPP



// include


#include<coffee-makers/Containers/FixedMatrix.hpp>
#include<coffee-makers/Containers/VariableMatrix.hpp>
#include<coffee-makers/Containers/FixedVector.hpp>
#include<coffee-makers/Containers/VariableVector.hpp>
#include<iostream>





// utility




namespace makers {


template<typename inputT, typename outputT, std::size_t row, std::size_t col,
	class = typename std::enable_if<std::is_convertible<inputT, outputT>::value>::type>
inline FixedMatrix<outputT, row, col> fix_Matrix(const VariableMatrix<inputT>& input_matrix) {
	if ((input_matrix.rows() != row) || (input_matrix.cols() != col)) {
		std::cerr << "Error: Invalid Matrix conversion." << std::endl;
		std::cerr << "Input: {" << input_matrix.rows() << ", " << input_matrix.cols() << "} >> ";
		std::cerr << "Output: {" << row << ", " << col << "}" << std::endl;
		std::exit(1);
	}

	FixedMatrix<outputT, row, col> result;

	for (std::size_t idx = 0; idx < row; ++idx) 
		for (std::size_t jdx = 0; jdx < col; ++jdx)
			result(idx, jdx) = input_matrix(idx, jdx);
	
	return result;
}


template<typename inputT, typename outputT, std::size_t dim,
	class = typename std::enable_if<std::is_convertible<inputT, outputT>::value>::type>
inline FixedVector<outputT, dim> fix_Vector(const VariableVector<inputT>& input_vector) {
	if (input_vector.size() != dim) {
		std::cerr << "Error: Invalid Vector conversion." << std::endl;
		std::cerr << "Input: {" << input_vector.size() << ", 1} >> ";
		std::cerr << "Output: {" << dim << ", 1}" << std::endl;
		std::exit(1);
	}
	
	FixedVector<outputT, dim> result;

	for (std::size_t idx = 0; idx < dim; ++idx)
		result(idx) = input_vector(idx);
	
	return result;
}


}


#endif /* COFFEE_MAKERS_CONTAINERS_HPP */
