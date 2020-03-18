#ifndef COFFEE_MAKERS_COMPILE_TIME_CONJUNCTION_HPP
#define COFFEE_MAKERS_COMPILE_TIME_CONJUNCTION_HPP

namespace makers {

template<typename...> struct Conjunction : std::true_type {};
template<typename T> struct Conjunction<T> : T {};
template<typename T, typename... Ts>
struct Conjunction<T, Ts...> : std::conditional<static_cast<bool>(T::value), Conjunction<Ts...>, T>::type {};

}


#endif /* COFFEE_MAKERS_COMPILE_TIME_CONJUNCTION_HPP */
