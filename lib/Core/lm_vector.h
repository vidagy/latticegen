#ifndef LATTICEGEN_LM_VECTOR_H
#define LATTICEGEN_LM_VECTOR_H

#include <vector>
#include <cmath>
#include "Exceptions.h"

namespace Core
{
  template<typename T>
  class lm_vector
  {
  public:
    lm_vector(const lm_vector<T> &other) = default;

    lm_vector(lm_vector<T> &&other) = default;

    lm_vector<T> &operator=(lm_vector<T> &other) = default;

    lm_vector<T> &operator=(lm_vector<T> &&other) = default;

    explicit lm_vector(unsigned int l_)
      : l_max(l_), data(std::vector<T>((l_max + 1) * (l_max + 1), T())) {}

    lm_vector(const std::vector<T> &data_)
      : l_max((const unsigned int) (lround(sqrt(data_.size())) - 1)), data(data_)
    {
      if ((l_max + 1) * (l_max + 1) != data_.size())
        THROW_INVALID_ARGUMENT("input data size is not square number data_.size() = " + std::to_string(data_.size()));
    }

    T &at(unsigned int l_, int m_)
    {
      if (l_ > l_max)
        THROW_INVALID_ARGUMENT(
          "Vector over-indexing: l_ = " + std::to_string(l_) + " l_max = " + std::to_string(l_max));
      if (abs(m_) > (int) l_)
        THROW_INVALID_ARGUMENT("m_ = " + std::to_string(m_) + " l_ = " + std::to_string(l_));
      return data[l_ * l_ + l_ + m_];
    }

    const T &at(unsigned int l_, int m_) const
    {
      if (l_ > l_max)
        THROW_INVALID_ARGUMENT(
          "Vector over-indexing: l_ = " + std::to_string(l_) + " l_max = " + std::to_string(l_max));
      if (abs(m_) > (int) l_)
        THROW_INVALID_ARGUMENT("m_ = " + std::to_string(m_) + " l_ = " + std::to_string(l_));
      return data[l_ * l_ + l_ + m_];
    }

    const unsigned int l_max;
  private:
    std::vector<T> data;
  };
}

#endif //LATTICEGEN_LM_VECTOR_H
