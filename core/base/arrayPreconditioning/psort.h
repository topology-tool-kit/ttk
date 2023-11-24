/// \ingroup base
/// \author Eve Le Guillou <eve.le-guillou@lip6.fr>
/// \date November 2023
///
/// In this file is implemented a distributed sort based on the algorithm psort
/// by David Cheng et al. The present software has been significantly modified.
/// In particular, the alltoall function has been modified to accept up to
/// INT_MAX*INT_MAX elements in the MPI call. If a wider range of identifiers
/// is needed, please check if alltoallv_c is available in the MPI
/// implementation. The changes to the two functions psort_split and psort_merge
/// are limited. Below can be found the licence of the original code.

/// \b Related \b publication \n
/// "A Novel Parallel Sorting Algorithm for Contemporary Architectures" \n
/// David Cheng, Viral Shah, John Gilbert, Alan Edelman \n
/// Submitted to ALENEX, (2006).

/*
Copyright (c) 2009, David Cheng, Viral B. Shah.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#pragma once
#include <Debug.h>

#ifdef TTK_ENABLE_MPI
namespace p_sort {

  template <typename T, class _Compare>
  class PermCompare {
  private:
    T *weights;
    _Compare comp;

  public:
    PermCompare(T *w, _Compare c) : weights(w), comp(c) {
    }
    bool operator()(int a, int b) {
      return comp(weights[a], weights[b]);
    }
  };

  template <class _Compare, typename _RandomAccessIter>
  void psort_split(_RandomAccessIter first,
                   _RandomAccessIter last,
                   ttk::SimplexId *dist,
                   _Compare comp,
                   std::vector<std::vector<ttk::SimplexId>> &right_ends,
                   MPI_Datatype &MPI_valueType,
                   MPI_Datatype &MPI_distanceType,
                   int nThreads) {

    using dataType =
      typename std::iterator_traits<_RandomAccessIter>::value_type;

    int n_real = ttk::MPIsize_;
    for(int i = 0; i < ttk::MPIsize_; ++i)
      if(dist[i] == 0) {
        n_real = i;
        break;
      }

    std::copy(dist, dist + ttk::MPIsize_, right_ends[ttk::MPIsize_].begin());

    // union of [0, right_end[i+1]) on each processor produces dist[i] total
    // values
    // ttk::SimplexId targets[ttk::MPIsize_-1];
    std::vector<ttk::SimplexId> targets(ttk::MPIsize_ - 1);
    std::partial_sum(dist, dist + (ttk::MPIsize_ - 1), targets.data());

    // keep a list of ranges, trying to "activate" them at each branch
    std::vector<std::pair<_RandomAccessIter, _RandomAccessIter>> d_ranges(
      ttk::MPIsize_ - 1);
    std::vector<std::pair<ttk::SimplexId *, ttk::SimplexId *>> t_ranges(
      ttk::MPIsize_ - 1);
    d_ranges[0] = std::make_pair(first, last);
    t_ranges[0]
      = std::make_pair(targets.data(), targets.data() + ttk::MPIsize_ - 1);

    std::vector<std::vector<ttk::SimplexId>> subdist(
      ttk::MPIsize_ - 1, std::vector<ttk::SimplexId>(ttk::MPIsize_));
    std::copy(dist, dist + ttk::MPIsize_, subdist[0].begin());

    // for each processor, d_ranges - first
    std::vector<std::vector<ttk::SimplexId>> outleft(
      ttk::MPIsize_ - 1, std::vector<ttk::SimplexId>(ttk::MPIsize_, 0));

    for(int n_act = 1; n_act > 0;) {

      for(int k = 0; k < n_act; ++k) {
        assert(subdist[k][ttk::MPIrank_]
               == d_ranges[k].second - d_ranges[k].first);
      }

      //------- generate n_act guesses

      // for the allgather, make a flat array of ttk::MPIsize_ chunks, each with
      // n_act elts
      std::vector<dataType> medians(ttk::MPIsize_ * n_act);
      for(int k = 0; k < n_act; ++k) {
        if(d_ranges[k].first != last) {
          dataType *ptr = &(*d_ranges[k].first);
          ttk::SimplexId index = subdist[k][ttk::MPIrank_] / 2;
          medians[ttk::MPIrank_ * n_act + k] = ptr[index];
        } else
          medians[ttk::MPIrank_ * n_act + k] = *(last - 1);
      }

      MPI_Allgather(MPI_IN_PLACE, n_act, MPI_valueType, &medians[0], n_act,
                    MPI_valueType, ttk::MPIcomm_);

      // compute the weighted median of medians
      std::vector<dataType> queries(n_act);

      std::vector<ttk::SimplexId> ms_perm(n_real);
      for(int k = 0; k < n_act; ++k) {

        for(int i = 0; i < n_real; ++i)
          ms_perm[i] = i * n_act + k;
        TTK_PSORT(nThreads, ms_perm.data(), ms_perm.data() + n_real,
                  PermCompare<dataType, _Compare>(&medians[0], comp));

        ttk::SimplexId mid
          = std::accumulate(
              subdist[k].begin(), subdist[k].end(), (ttk::SimplexId)0)
            / 2;
        ttk::SimplexId query_ind = -1;
        for(int i = 0; i < n_real; ++i) {
          if(subdist[k][ms_perm[i] / n_act] == 0)
            continue;

          mid -= subdist[k][ms_perm[i] / n_act];
          if(mid <= 0) {
            query_ind = ms_perm[i];
            break;
          }
        }

        assert(query_ind >= 0);
        queries[k] = medians[query_ind];
      }
      //------- find min and max ranks of the guesses
      std::vector<ttk::SimplexId> ind_local(2 * n_act);

      for(int k = 0; k < n_act; ++k) {
        std::pair<_RandomAccessIter, _RandomAccessIter> ind_local_p
          = equal_range(
            d_ranges[k].first, d_ranges[k].second, queries[k], comp);

        ind_local[2 * k] = ind_local_p.first - first;
        ind_local[2 * k + 1] = ind_local_p.second - first;
      }

      std::vector<ttk::SimplexId> ind_all(2 * n_act * ttk::MPIsize_);
      MPI_Allgather(ind_local.data(), 2 * n_act, MPI_distanceType,
                    ind_all.data(), 2 * n_act, MPI_distanceType, ttk::MPIcomm_);
      // sum to get the global range of indices
      std::vector<std::pair<ttk::SimplexId, ttk::SimplexId>> ind_global(n_act);
      for(int k = 0; k < n_act; ++k) {
        ind_global[k] = std::make_pair(0, 0);
        for(int i = 0; i < ttk::MPIsize_; ++i) {
          ind_global[k].first += ind_all[2 * (i * n_act + k)];
          ind_global[k].second += ind_all[2 * (i * n_act + k) + 1];
        }
      }

      // state to pass on to next iteration
      std::vector<std::pair<_RandomAccessIter, _RandomAccessIter>> d_ranges_x(
        ttk::MPIsize_ - 1);
      std::vector<std::pair<ttk::SimplexId *, ttk::SimplexId *>> t_ranges_x(
        ttk::MPIsize_ - 1);
      std::vector<std::vector<ttk::SimplexId>> subdist_x(
        ttk::MPIsize_ - 1, std::vector<ttk::SimplexId>(ttk::MPIsize_));
      std::vector<std::vector<ttk::SimplexId>> outleft_x(
        ttk::MPIsize_ - 1, std::vector<ttk::SimplexId>(ttk::MPIsize_, 0));
      int n_act_x = 0;

      for(int k = 0; k < n_act; ++k) {
        ttk::SimplexId *split_low = std::lower_bound(
          t_ranges[k].first, t_ranges[k].second, ind_global[k].first);
        ttk::SimplexId *split_high = std::upper_bound(
          t_ranges[k].first, t_ranges[k].second, ind_global[k].second);

        // iterate over targets we hit
        for(ttk::SimplexId *s = split_low; s != split_high; ++s) {
          assert(*s > 0);
          // a bit sloppy: if more than one target in range, excess won't zero
          // out
          ttk::SimplexId excess = *s - ind_global[k].first;
          // low procs to high take excess for stability
          for(int i = 0; i < ttk::MPIsize_; ++i) {
            ttk::SimplexId amount
              = std::min(ind_all[2 * (i * n_act + k)] + excess,
                         ind_all[2 * (i * n_act + k) + 1]);
            right_ends[(s - targets.data()) + 1][i] = amount;
            excess -= amount - ind_all[2 * (i * n_act + k)];
          }
        }

        if((split_low - t_ranges[k].first) > 0) {
          t_ranges_x[n_act_x] = std::make_pair(t_ranges[k].first, split_low);
          // lop off local_ind_low..end
          d_ranges_x[n_act_x]
            = std::make_pair(d_ranges[k].first, first + ind_local[2 * k]);
          for(int i = 0; i < ttk::MPIsize_; ++i) {
            subdist_x[n_act_x][i]
              = ind_all[2 * (i * n_act + k)] - outleft[k][i];
            outleft_x[n_act_x][i] = outleft[k][i];
          }
          ++n_act_x;
        }

        if((t_ranges[k].second - split_high) > 0) {
          t_ranges_x[n_act_x] = std::make_pair(split_high, t_ranges[k].second);
          // lop off begin..local_ind_high
          d_ranges_x[n_act_x]
            = std::make_pair(first + ind_local[2 * k + 1], d_ranges[k].second);
          for(int i = 0; i < ttk::MPIsize_; ++i) {
            subdist_x[n_act_x][i] = outleft[k][i] + subdist[k][i]
                                    - ind_all[2 * (i * n_act + k) + 1];
            outleft_x[n_act_x][i] = ind_all[2 * (i * n_act + k) + 1];
          }
          ++n_act_x;
        }
      }

      t_ranges = t_ranges_x;
      d_ranges = d_ranges_x;
      subdist = subdist_x;
      outleft = outleft_x;
      n_act = n_act_x;
    }
    TTK_FORCE_USE(nThreads);
  }

  /**
   * This function as been heavily modified compared to the original psort to
   * account for the case where more than INT_MAX elements need to be sent in
   * the alltoall.
   */
  template <typename dataType>
  static void alltoall(std::vector<std::vector<ttk::SimplexId>> &right_ends,
                       std::vector<dataType> &data,
                       std::vector<dataType> &trans_data,
                       ttk::SimplexId *boundaries,
                       MPI_Datatype &MPI_valueType,
                       MPI_Datatype &MPI_distanceType) {
    // MPI_Alltoallv requires integers for displacements and counts.
    // If there is a need to send more, then a second strategy is used
    ttk::SimplexId n_loc_ = data.size();
    // char and not boolean because MPI doesn't know boolean
    char overflowInt{0};
    ttk::SimplexId chunkSize = n_loc_ / INT_MAX + 1;
    if(n_loc_ > INT_MAX) {
      overflowInt = true;
    }
    int n_loc = static_cast<int>(n_loc_);
    // Calculate the counts for redistributing data and detects potential
    // overflows of INT_MAX
    std::vector<ttk::SimplexId> send_counts(ttk::MPIsize_);
    std::vector<ttk::SimplexId> recv_counts(ttk::MPIsize_);
#pragma omp parallel for reduction(+ : overflowInt)
    for(int i = 0; i < ttk::MPIsize_; ++i) {
      ttk::SimplexId scount
        = right_ends[i + 1][ttk::MPIrank_] - right_ends[i][ttk::MPIrank_];
      ttk::SimplexId rcount
        = right_ends[ttk::MPIrank_ + 1][i] - right_ends[ttk::MPIrank_][i];
      if(scount > INT_MAX || rcount > INT_MAX) {
        overflowInt = 1;
      }
      send_counts[i] = scount;
      recv_counts[i] = rcount;
    }
    MPI_Allreduce(
      MPI_IN_PLACE, &overflowInt, 1, MPI_CHAR, MPI_SUM, ttk::MPIcomm_);
    MPI_Allreduce(
      MPI_IN_PLACE, &chunkSize, 1, MPI_distanceType, MPI_MAX, ttk::MPIcomm_);
    if(overflowInt == 0) {
      std::vector<int> send_counts_int(ttk::MPIsize_);
      std::vector<int> send_disps_int(ttk::MPIsize_);
      std::vector<int> recv_counts_int(ttk::MPIsize_);
      std::vector<int> recv_disps_int(ttk::MPIsize_);

      for(int i = 0; i < ttk::MPIsize_; i++) {
        send_counts_int[i] = static_cast<int>(send_counts[i]);
        recv_counts_int[i] = static_cast<int>(recv_counts[i]);
      }

      recv_disps_int[0] = 0;
      std::partial_sum(recv_counts_int.data(),
                       recv_counts_int.data() + ttk::MPIsize_ - 1,
                       recv_disps_int.data() + 1);

      send_disps_int[0] = 0;
      std::partial_sum(send_counts_int.data(),
                       send_counts_int.data() + ttk::MPIsize_ - 1,
                       send_disps_int.data() + 1);
      // Do the transpose
      MPI_Alltoallv(data.data(), send_counts_int.data(), send_disps_int.data(),
                    MPI_valueType, trans_data.data(), recv_counts_int.data(),
                    recv_disps_int.data(), MPI_valueType, ttk::MPIcomm_);

      for(int i = 0; i < ttk::MPIsize_; ++i)
        boundaries[i] = (ttk::SimplexId)recv_disps_int[i];
      boundaries[ttk::MPIsize_] = (ttk::SimplexId)n_loc; // for the merging

      return;
    }

    // if the alltoall needs to send more than INT_MAX elements,
    // then multiple messages, of each at most INT_MAX elements,
    // are sent until all elements are sent.
    std::vector<ttk::SimplexId> send_disps(ttk::MPIsize_);
    std::vector<ttk::SimplexId> recv_disps(ttk::MPIsize_);
    send_disps[0] = 0;
    std::partial_sum(send_counts.data(), send_counts.data() + ttk::MPIsize_ - 1,
                     send_disps.data() + 1);
    recv_disps[0] = 0;
    std::partial_sum(recv_counts.data(), recv_counts.data() + ttk::MPIsize_ - 1,
                     recv_disps.data() + 1);

    int moreToSend{1};
    int count{0};
    std::vector<int> partial_recv_count(ttk::MPIsize_, 0);
    std::vector<int> partial_send_count(ttk::MPIsize_, 0);
    std::vector<int> partial_recv_displs(ttk::MPIsize_, 0);
    std::vector<int> partial_send_displs(ttk::MPIsize_, 0);
    std::vector<dataType> send_buffer_64bits(INT_MAX);
    std::vector<dataType> recv_buffer_64bits(INT_MAX);
    ttk::SimplexId messageSize = std::max(INT_MAX / ttk::MPIsize_ - 1, 1);
    while(moreToSend) {
      // Preparing the send buffer of at most INT_MAX elements
      send_buffer_64bits.resize(0);
      moreToSend = 0;
      for(int i = 0; i < ttk::MPIsize_; i++) {
        if(i > 0) {
          partial_send_displs[i]
            = partial_send_displs[i - 1] + partial_send_count[i - 1];
          partial_recv_displs[i]
            = partial_recv_displs[i - 1] + partial_recv_count[i - 1];
        }
        if(send_counts[i] - count * messageSize > 0) {
          moreToSend = 1;
          if(send_counts[i] - count * messageSize > messageSize) {
            partial_send_count[i] = messageSize;
          } else {
            partial_send_count[i] = send_counts[i] - count * messageSize;
          }
          std::copy(data.begin() + send_disps[i] + count * messageSize,
                    data.begin() + send_disps[i] + count * messageSize
                      + partial_send_count[i],
                    send_buffer_64bits.begin() + partial_send_displs[i]);
        } else {
          partial_send_count[i] = 0;
        }
        if(recv_counts[i] - count * messageSize > 0) {
          if(recv_counts[i] - count * messageSize > messageSize) {
            partial_recv_count[i] = messageSize;
          } else {
            partial_recv_count[i] = recv_counts[i] - count * messageSize;
          }
        } else {
          partial_recv_count[i] = 0;
        }
      }
      MPI_Alltoallv(send_buffer_64bits.data(), partial_send_count.data(),
                    partial_send_displs.data(), MPI_valueType,
                    recv_buffer_64bits.data(), partial_recv_count.data(),
                    partial_recv_displs.data(), MPI_valueType, ttk::MPIcomm_);
      // Receiving the buffer and placing it in the right vector
      for(int i = 0; i < ttk::MPIsize_; i++) {
        if(partial_recv_count[i] > 0) {
          std::copy(recv_buffer_64bits.begin() + partial_recv_displs[i],
                    recv_buffer_64bits.begin() + partial_recv_displs[i]
                      + partial_recv_count[i],
                    trans_data.begin() + recv_disps[i] + count * messageSize);
        }
      }
      count++;
      // Test to see if more messages need to be sent
      MPI_Allreduce(
        MPI_IN_PLACE, &moreToSend, 1, MPI_INTEGER, MPI_SUM, ttk::MPIcomm_);
    }

    for(int i = 0; i < ttk::MPIsize_; ++i)
      boundaries[i] = recv_disps[i];
    boundaries[ttk::MPIsize_] = n_loc_; // for the merging

    return;
  }

  template <class _Compare, typename _RandomAccessIter>
  void psort_merge(_RandomAccessIter in,
                   _RandomAccessIter out,
                   ttk::SimplexId *disps,
                   _Compare comp,
                   _Compare oppositeComp) {

    if(ttk::MPIsize_ == 1) {
      std::copy(in, in + disps[ttk::MPIsize_], out);
      return;
    }

    _RandomAccessIter bufs[2] = {in, out};

    std::vector<ttk::SimplexId> locs(ttk::MPIsize_, 0);

    ttk::SimplexId next = 1;
    while(true) {
      ttk::SimplexId stride = next * 2;
      if(stride >= ttk::MPIsize_)
        break;

      for(ttk::SimplexId i = 0; i + next < ttk::MPIsize_; i += stride) {
        ttk::SimplexId end_ind
          = std::min(i + stride, (ttk::SimplexId)ttk::MPIsize_);

        std::merge(bufs[locs[i]] + disps[i], bufs[locs[i]] + disps[i + next],
                   bufs[locs[i + next]] + disps[i + next],
                   bufs[locs[i + next]] + disps[end_ind],
                   bufs[1 - locs[i]] + disps[i], comp);
        locs[i] = 1 - locs[i];
      }

      next = stride;
    }

    // now have 4 cases for final merge
    if(locs[0] == 0) {
      // 00, 01 => out of place
      std::merge(in, in + disps[next], bufs[locs[next]] + disps[next],
                 bufs[locs[next]] + disps[ttk::MPIsize_], out, comp);
    } else if(locs[next] == 0) {
      // 10 => backwards out of place
      std::merge(
        std::reverse_iterator<_RandomAccessIter>(in + disps[ttk::MPIsize_]),
        std::reverse_iterator<_RandomAccessIter>(in + disps[next]),
        std::reverse_iterator<_RandomAccessIter>(out + disps[next]),
        std::reverse_iterator<_RandomAccessIter>(out),
        std::reverse_iterator<_RandomAccessIter>(out + disps[ttk::MPIsize_]),
        oppositeComp);
    } else {
      // 11 => in-place
      std::inplace_merge(
        out, out + disps[next], out + disps[ttk::MPIsize_], comp);
    }
  }

  /**
   * @brief This function sorts a data vector in a distributed-memory context.
   * It will first sort locally using multi-threading and then operate a global
   * sort using a split and merge strategy. After sorting, each process will
   * have the same number of vertices, but they will be sorted globally (process
   * 0 will have all the smallest ones, etc).
   *
   * @param data The data vector to be sorted in distributed
   * @param comp The comparator to use for the sort
   * @param oppositeComp The opposite of the above comparator
   * (for example, for greater than, its opposite would be less than)
   * @param dist The number of elements to sort on each process
   * @param MPI_valueType The MPI type of the data
   * @param MPI_distanceType The MPI type of the indices
   * @param nThreads The number of threads to use
   */
  template <typename dataType, typename _Compare>
  void parallel_sort(std::vector<dataType> &data,
                     _Compare comp,
                     _Compare oppositeComp,
                     std::vector<ttk::SimplexId> &dist,
                     MPI_Datatype &MPI_valueType,
                     MPI_Datatype &MPI_distanceType,
                     int nThreads) {

    // Sort the data locally
    TTK_PSORT(nThreads, data.begin(), data.end(), comp);

    if(ttk::MPIsize_ == 1)
      return;

    // Find splitters
    std::vector<std::vector<ttk::SimplexId>> right_ends(
      ttk::MPIsize_ + 1, std::vector<ttk::SimplexId>(ttk::MPIsize_, 0));
    psort_split<_Compare>(data.begin(), data.end(), dist.data(), comp,
                          right_ends, MPI_valueType, MPI_distanceType,
                          nThreads);

    // Communicate to destination
    ttk::SimplexId n_loc = data.size();
    std::vector<dataType> trans_data(n_loc);

    std::vector<ttk::SimplexId> boundaries(ttk::MPIsize_ + 1);
    alltoall(right_ends, data, trans_data, boundaries.data(), MPI_valueType,
             MPI_distanceType);

    psort_merge<_Compare>(
      trans_data.data(), data.data(), boundaries.data(), comp, oppositeComp);

    return;
  }

} // namespace p_sort

#endif // TTK_ENABLE_MPI
