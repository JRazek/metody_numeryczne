#include <concepts>
#include <functional>

namespace utils {

template <std::floating_point T>
using RealFunction = std::function<T(T)>;

template <std::floating_point T>
struct Integral {
  T low_{};
  T high_{};
  RealFunction<T> function_{};
};

}  // namespace utils
