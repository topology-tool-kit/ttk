#pragma once

#include <list>
#include <map>

namespace ttk {
  /**
   * @brief LRU cache implementation
   *
   * Adapted from boost/compute/details/lru_cache.hpp
   */
  template <class KeyType, class ValueType>
  class LRUCache {
  public:
    // default cache with at most 8 entries
    LRUCache() : capacity_{8} {
    }

    LRUCache(const std::size_t capacity) : capacity_(capacity) {
    }

    inline bool empty() const {
      return this->map_.empty();
    }

    inline std::size_t capacity() const {
      return this->capacity_;
    }

    inline std::size_t size() const {
      return this->map_.size();
    }

    inline bool contains(const KeyType &key) const {
      return this->map_.find(key) != this->map_.end();
    }

    inline void clear() {
      this->map_.clear();
      this->queue_.clear();
    }

    /**
     * @brief Insert new (key, value) entry
     *
     * Do nothing if key already in use
     */
    inline void insert(const KeyType &key, const ValueType &value) {
      if(this->contains(key)) {
        return; // key already in use
      }

      if(this->size() >= this->capacity_) {
        // cache is full, evict the least recently used entry
        this->map_.erase(this->queue_.back());
        this->queue_.pop_back();
      }

      // insert new entry
      this->queue_.push_front(key);
      this->map_.emplace(key, std::make_pair(value, this->queue_.begin()));
    }

    /**
     * @brief Get value pointer from key
     *
     * @return nullptr if cache miss
     */
    inline ValueType *get(const KeyType &key) {
      const auto i = this->map_.find(key);
      if(i == this->map_.end()) {
        // cache miss, nullptr
        return {};
      }

      // iterator to entry in LRU list
      const auto oldIt = i->second.second;
      if(oldIt != this->queue_.begin()) {
        // update entry, now first in the LRU list
        this->queue_.erase(oldIt);
        this->queue_.push_front(key);

        // update queue position in map
        this->map_[key].second = this->queue_.begin();
      }

      // return value address
      return &i->second.first;
    }

  private:
    using ListType = std::list<KeyType>;
    using MapType
      = std::map<KeyType, std::pair<ValueType, typename ListType::iterator>>;

    MapType map_;
    ListType queue_;
    std::size_t capacity_;
  };

} // namespace ttk
