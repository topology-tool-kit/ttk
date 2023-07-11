#pragma once

#include <Debug.h>

namespace ttk {
  /**
   * @brief Replacement for std::vector<std::vector<SimplexId>>
   *
   * Use this when instead of a std::vector<std::vector<SimplexId>>
   * when the data is set once and not modified afterwards.
   */
  class FlatJaggedArray {
    // flattened sub-vectors data
    std::vector<SimplexId> data_;
    // offset for every sub-vector
    std::vector<SimplexId> offsets_;

  public:
    // ############## //
    // Initialization //
    // ############## //

    /**
     * @brief Set internal data from pre-existing vectors
     */
    inline void setData(std::vector<SimplexId> &&data,
                        std::vector<SimplexId> &&offsets) {
      this->data_ = std::move(data);
      this->offsets_ = std::move(offsets);
    }

    /**
     * @brief Clear the underlying vectors
     */
    inline void clear() {
      this->data_.clear();
      this->offsets_.clear();
    }

    // ############################## //
    // Mimic the vector of vector API //
    // ############################## //

    /**
     * @brief Get the size of a particular sub-vector
     */
    inline SimplexId size(SimplexId id) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(id < 0 || id > (SimplexId)offsets_.size() - 1) {
        return -1;
      }
#endif
      return this->offsets_[id + 1] - this->offsets_[id];
    }

    /**
     * @brief Get the offset of a particular sub-vector
     */
    inline SimplexId offset(SimplexId id) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(id < 0 || id > (SimplexId)offsets_.size()) {
        return -1;
      }
#endif
      return this->offsets_[id];
    }

    struct Slice {
      const size_t len;
      const SimplexId *const ptr;
      inline SimplexId operator[](const size_t id) const {
#ifndef TTK_ENABLE_KAMIKAZE
        if(id >= this->len || ptr == nullptr) {
          return -1;
        }
#endif
        return this->ptr[id];
      }
      inline const SimplexId *data() const {
        return this->ptr;
      }
      inline const SimplexId *begin() const {
        return this->ptr;
      }
      inline const SimplexId *end() const {
        return this->ptr + this->len;
      }
      inline size_t size() const {
        return this->len;
      }
    };

    struct Iterator {
      size_t id;
      const FlatJaggedArray &parent;
      inline void operator++() {
        this->id++;
      }
      inline Slice operator*() const {
        return this->parent[this->id];
      }
      inline bool operator!=(const Iterator &other) const {
        return this->id != other.id;
      }
    };

    inline Iterator begin() const {
      return {0, *this};
    }
    inline Iterator end() const {
      return {this->size(), *this};
    }

    inline Slice operator[](const size_t id) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(id >= this->offsets_.size()) {
        return {0, nullptr};
      }
#endif
      return {static_cast<size_t>(this->size(id)), this->get_ptr(id, 0)};
    }

    /**
     * @brief Returns the data inside the sub-vectors
     */
    inline SimplexId get(SimplexId id, SimplexId local) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(id < 0 || id > (SimplexId)offsets_.size() - 1) {
        return -1;
      }
      if(local < 0 || local >= this->size(id)) {
        return -2;
      }
#endif
      return this->data_[this->offsets_[id] + local];
    }

    /**
     * @brief Returns a const pointer to the data inside the sub-vectors
     */
    inline const SimplexId *get_ptr(SimplexId id, SimplexId local) const {
#ifndef TTK_ENABLE_KAMIKAZE
      if(id < 0 || id > (SimplexId)offsets_.size() - 1) {
        return {};
      }
      if(local < 0 || local >= this->size(id)) {
        return {};
      }
#endif
      return &this->data_[this->offsets_[id] + local];
    }

    /**
     * @brief Returns a const pointer to the offset member
     */
    inline const SimplexId *offset_ptr() const {
      return offsets_.data();
    }

    /**
     * @brief Returns the number of sub-vectors
     */
    inline size_t size() const {
      if(this->empty()) {
        return 0;
      }
      return this->offsets_.size() - 1;
    }

    /**
     * @brief Returns the size of the data_ member
     */
    inline size_t dataSize() const {
      return this->data_.size();
    }

    /**
     * @brief If the underlying buffers are empty
     */
    inline bool empty() const {
      return this->data_.empty() || this->offsets_.empty();
    }

    /**
     * @brief Computes the memory footprint of the array
     */
    inline size_t footprint() const {
      return (this->data_.size() + this->offsets_.size()) * sizeof(SimplexId);
    }

    // #################### //
    // Conversion utilities //
    // #################### //

    /**
     * @brief Fill buffers from a std::vector<std::vector<SimplexId>>
     *
     * Templated to also accept Boost small_vectors.
     */
    template <typename T>
    void fillFrom(const std::vector<T> &src, int threadNumber = 1) {
      this->offsets_.resize(src.size() + 1);
      for(size_t i = 0; i < src.size(); ++i) {
        this->offsets_[i + 1] = this->offsets_[i] + src[i].size();
      }
      this->data_.resize(this->offsets_.back());
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber)
#endif // TTK_ENABLE_OPENMP
      for(size_t i = 0; i < src.size(); ++i) {
        for(size_t j = 0; j < src[i].size(); ++j) {
          this->data_[this->offsets_[i] + j] = src[i][j];
        }
      }
      TTK_FORCE_USE(threadNumber);
    }

    /**
     * @brief Copy buffers to a std::vector<std::vector<SimplexId>>
     */
    void copyTo(std::vector<std::vector<SimplexId>> &dst,
                int threadNumber = 1) const {
      dst.resize(this->size());
      for(size_t i = 0; i < this->size(); ++i) {
        dst[i].resize((*this)[i].size());
      }
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber)
#endif // TTK_ENABLE_OPENMP
      for(size_t i = 0; i < this->size(); ++i) {
        for(size_t j = 0; j < dst[i].size(); ++j) {
          dst[i][j] = (*this)[i][j];
        }
      }
      TTK_FORCE_USE(threadNumber);
    }

    // ##################### //
    // Useless I/O utilities //
    // ##################### //

    /**
     * @brief Write content into a file
     */
    inline void writeToFile(const std::string &fName) const {
      std::ofstream out(fName);
      for(const auto slice : *this) {
        for(const auto el : slice) {
          out << el << " ";
        }
        out << '\n';
      }
    }

    /**
     * @brief Also write std::vector<std::vector<SimplexId>> to disk
     */
    static inline void
      writeToFile(const std::string &fName,
                  const std::vector<std::vector<SimplexId>> &src) {
      std::ofstream out(fName);
      for(const auto &vec : src) {
        for(const auto el : vec) {
          out << el << " ";
        }
        out << '\n';
      }
    }
  };
} // namespace ttk
