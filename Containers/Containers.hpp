#ifndef COFFEE_MAKERS_CONTAINERS_HPP
#define COFFEE_MAKERS_CONTAINERS_HPP



// include


#include<coffee-makers/Containers/Matrix.hpp>
#include<coffee-makers/Containers/Vector.hpp>
#include<iostream>


namespace makers {


// utility



// alias

using Mat2i = Matrix<int, 2, 2>;
using Mat3i = Matrix<int, 3, 3>;
using Mat4i = Matrix<int, 4, 4>;
using Mat5i = Matrix<int, 5, 5>;
using MatXi = Matrix<int, Variable, Variable>;

using Mat2f = Matrix<float, 2, 2>;
using Mat3f = Matrix<float, 3, 3>;
using Mat4f = Matrix<float, 4, 4>;
using Mat5f = Matrix<float, 5, 5>;
using MatXf = Matrix<float, Variable, Variable>;

using Mat2d = Matrix<double, 2, 2>;
using Mat3d = Matrix<double, 3, 3>;
using Mat4d = Matrix<double, 4, 4>;
using Mat5d = Matrix<double, 5, 5>;
using MatXd = Matrix<double, Variable, Variable>;

using Vec2i = Vector<int, 2>;
using Vec3i = Vector<int, 3>;
using Vec4i = Vector<int, 4>;
using Vec5i = Vector<int, 5>;
using VecXi = Vector<int, Variable>;

using Vec2f = Vector<float, 2>;
using Vec3f = Vector<float, 3>;
using Vec4f = Vector<float, 4>;
using Vec5f = Vector<float, 5>;
using VecXf = Vector<float, Variable>;

using Vec2d = Vector<double, 2>;
using Vec3d = Vector<double, 3>;
using Vec4d = Vector<double, 4>;
using Vec5d = Vector<double, 5>;
using VecXd = Vector<double, Variable>;

}




#endif /* COFFEE_MAKERS_CONTAINERS_HPP */
