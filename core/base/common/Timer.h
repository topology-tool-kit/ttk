#include <chrono>

namespace ttk {

  class Timer {
    using ClockType = std::chrono::steady_clock;
    using TimeType = std::chrono::time_point<ClockType>;
    using DiffType = std::chrono::duration<double>;

  public:
    Timer() = default;

    inline double getElapsedTime() {
      const auto end = this->getTimeStamp();
      const DiffType diff = end - start_;
      return diff.count();
    }

    inline void reStart() {
      start_ = this->getTimeStamp();
    }

  protected:
    inline TimeType getTimeStamp() {
      return ClockType::now();
    }

    TimeType start_{this->getTimeStamp()};
  };

} // namespace ttk
