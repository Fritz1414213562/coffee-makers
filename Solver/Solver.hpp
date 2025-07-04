#ifndef COFFEE_MAKERS_SOLVER_HPP
#define COFFEE_MAKERS_SOLVER_HPP
// interface
#include<coffee-makers/Solver/EigenSolver.hpp>
#include<coffee-makers/Solver/SingularValueSolver.hpp>


// calculation methods
#include<coffee-makers/Solver/JacobiMethod.hpp>







namespace makers {



//// alias

template<typename MatT>
//using JacobiEigenSolver = EigenSolver<JacobiMethod<MatT>>;
using JacobiEigenSolver = EigenSolver<JacobiMethod, MatT>;

template<typename MatT>
using JacobiSVSolver = SingularValueSolver<JacobiMethod, MatT>;


}


#endif /* COFFEE_MAKERS_SOLVER_HPP */
