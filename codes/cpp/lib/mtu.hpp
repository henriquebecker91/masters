#ifndef HBM_MTU_HPP
#define HBM_MTU_HPP

#include "ukp_common.hpp"
#include "wrapper.hpp"
#include "type_name.hpp"

#include <vector>
#include <chrono>
#include <cassert>

#ifndef HBM_PROFILE_PRECISION
  #define HBM_PROFILE_PRECISION 5
#endif

#ifndef HBM_PRINT_VAR
  #define HBM_PRINT_VAR(var) out << #var ": " << var << std::endl
#endif

namespace hbm {
  template <typename W, typename P, typename I>
  struct mtu1_extra_info_t : extra_info_t {
    std::string algorithm_name{"cpp-mtu1"};
    #ifdef HBM_PROFILE
    double sort_time{0};        ///< Time used sorting items.
    double vector_alloc_time{0};///< Time used allocating vectors for DP.
    double linear_comp_time{0}; ///< Time used by linear time preprocessing.
    double bb_time{0};     ///< Time used creating partial solutions.
    double sol_time{0};    ///< Time used to assemble solution.
    double total_time{0};  ///< Time used by all the algorithm.
    #endif //HBM_PROFILE
    I n{0};      ///< Instance number of items.
    W c{0};      ///< Instance capacity.

    virtual std::string gen_info(void) {
      std::stringstream out("");

      HBM_PRINT_VAR(algorithm_name);
      out <<  "git_head_at_compilation: "
          << HBM_GIT_HEAD_AT_COMPILATION << std::endl;
      out << "type_W: " << hbm::type_name<W>::get() << std::endl;
      out << "type_P: " << hbm::type_name<P>::get() << std::endl;
      out << "type_I: " << hbm::type_name<I>::get() << std::endl;
      HBM_PRINT_VAR(c);
      HBM_PRINT_VAR(n);

      #ifdef HBM_PROFILE
      out << "HBM_PROFILE: defined" << std::endl;
      const double sum_time = sort_time + vector_alloc_time +
        linear_comp_time + bb_time + sol_time;

      std::streamsize old_precision = out.precision(HBM_PROFILE_PRECISION);
      const int two_first_digits_and_period = 3;
      const int percent_size = two_first_digits_and_period+HBM_PROFILE_PRECISION;
      std::ios_base::fmtflags old_flags = out.setf(std::ios::fixed, std:: ios::floatfield);
      char old_fill = out.fill(' ');

      #ifndef HBM_PRINT_TIME
        #define HBM_PRINT_TIME(name, var)\
          out << name << " time: " << var << "s (";\
          out << std::setw(percent_size) << (var/total_time)*100.0;\
          out << "%)" << std::endl
      #endif

      HBM_PRINT_TIME("Sort", sort_time);
      HBM_PRINT_TIME("Vect", vector_alloc_time);
      HBM_PRINT_TIME("O(n)", linear_comp_time);
      HBM_PRINT_TIME("B&B ", bb_time);
      HBM_PRINT_TIME("sol ", sol_time);
      HBM_PRINT_TIME("Sum ", sum_time);
      HBM_PRINT_TIME("All ", total_time);

      out.fill(old_fill);
      out.setf(old_flags);
      out.precision(old_precision);
      #else  //HBM_PROFILE
      out << "HBM_PROFILE: undefined" << std::endl;
      #endif //HBM_PROFILE

      return out.str();
    }
  };

  template <typename W, typename P, typename I>
  struct mtu2_extra_info_t : mtu1_extra_info_t<W, P, I> {
    mtu2_extra_info_t(void) : mtu1_extra_info_t<W, P, I>() {
      mtu1_extra_info_t<W, P, I>::algorithm_name = "cpp-mtu2";
    }
  };

  namespace hbm_mtu_impl {
    using namespace std;
    using namespace std::chrono;

    // "Knapsack Problems" p. 93, equation 3.19 
    template<typename W, typename P>
    inline P u1(
      const item_t<W, P> &bi,
      const item_t<W, P> &bi2,
      const W &c
    ) {
      W c_ = c % bi.w;
      return (c/bi.w)*bi.p + (c_*bi2.p)/bi2.w;
    }

    template<typename W, typename P>
    inline P u3(
      const item_t<W, P> &bi,
      const item_t<W, P> &bi2,
      const item_t<W, P> &bi3,
      const W &c
    ) {
      // The numbers at the right are the equations reference number on
      // "Knapsack Problems" from S. Martello and P. Toth (section 3.6).

      // Auxiliar values
      W c_ = c % bi.w;                              // 3.18
      W c_prime = c_ % bi2.w;                       // 3.22
      P z_prime = (c/bi.w)*bi.p + (c_/bi2.w)*bi2.p; // 3.21

      // u0 computation
      P u0 = z_prime + (c_prime*bi3.p)/bi3.w; // 3.23

      // u1_ computation (all bellow is 3.25)
      // We use the remainder of bi.w (rem) to simulate the
      // rounding-up of a result used in the article.
      W w2_less_c_prime = bi2.w - c_prime;
      W rem = w2_less_c_prime % bi.w;
      W quo = w2_less_c_prime / bi.w;
      W res = rem ? quo + 1 : quo;

      // On the equation immediatly below, at the book, p2/w2 is shown as a
      // fraction multiplying the rest of the equation. So there's no round
      // down in the equation but there's one round down here (assuming that
      // W and P are integers). This is safe to do because this value is
      // subtracted by an already rounded value (integer) and then the result
      // is rounded down (would have thrown away the smaller-than-one part
      // of the value anyway).
      P l = ((c_prime + res*bi.w)*bi2.p)/bi2.w;
      P u1_ = z_prime + (l - res*bi.p);

      // u3 computation
      return max(u1_, u0); // 3.26
    }

    // See comment over hbm::inner_mtu1.
    template<typename W, typename P, typename I>
    void inner_mtu1(
      myvector<W> &w_,
      myvector<P> &p_,
      const I &n,
      const W &c,
      P &z,
      myvector<W> &x
    ) {
      const item_t<W, P>  bi(w_[1], p_[1]),
                          bi2(w_[2], p_[2]),
                          bi3(w_[3], p_[3]);
      // Variables needs to be initialized before we start jumping.
      W y;
      P u;
      I h, i;

      // 1. Initialize
      //step_1: // unused step
      z = 0;
      P z_ = 0;
      W c_ = c;

      w_[n+1] = numeric_limits<W>::max();
      p_[n+1] = 0;
      // Avoid changing w_ and p_ for the rest of the algorithm.
      const auto &w = w_;
      const auto &p = p_;

      vector<W> x_(n+1, 0);
      P upper_u = hbm_mtu_impl::u3(bi, bi2, bi3, c);

      myvector<W> m;
      m.resize(n+1); // This resize don't initialize the m contents with zero

      m[n] = w[n + 1];
      for (I k = n - 1; k > 0; --k) m[k] = min(m[k + 1], w[k + 1]);

      I j = 1;

      // 2. build a new current solution
      step_2:
      while (w[j] > c_)
        if (z >= z_ + (c_*p[j+1])/w[j+1]) goto step_5; else j = j + 1;

      y = c_/w[j];
      u = ((c_ - y*w[j])*p[j+1])/w[j+1];

      if (z >= z_ + y*p[j] + u) goto step_5;
      if (u == 0) goto step_4;

      // 3. save the current solution
      //step_3: // unused step
      c_ = c_ - y*w[j];
      z_ = z_ + y*p[j];
      x_[j] = y;
      j = j + 1;
      if (c_ >= m[j-1]) goto step_2;
      if (z >= z_) goto step_5;
      y = 0;

      // 4. update the best solution so far
      step_4:
      z = z_ + y*p[j];
      for (I k = 1; k < j; ++k) x[k] = x_[k];
      x[j] = y;
      for (I k = j + 1; k <= n; ++k) x[k] = 0;
      if (z == upper_u) goto stop;

      // 5. backtrack
      step_5:
      i = 0; // We use zero as an sentinel value
      for (I k = j - 1; k > 0; --k) if (x_[k] > 0) { i = k; break; }
      if (i == 0) goto stop; // "if no such i then return"
      c_ = c_ + w[i];
      z_ = z_ - p[i];
      x_[i] = x_[i] - 1;
      if (z >= z_ + (c_*p[i+1])/w[i+1]) {
        c_ = c_ + w[i]*x_[i];
        z_ = z_ - p[i]*x_[i];
        x_[i] = 0;
        j = i;
        goto step_5;
      }
      j = i + 1;
      if (c_ - w[i] >= m[i]) goto step_2;
      h = i;

      // 6. try to replace one item of type i with items of type h
      step_6:
      h = h + 1;
      if (z >= z_ + (c_*p[h])/w[h]) goto step_5;
      if (w[h] == w[i]) goto step_6;
      if (w[h] > w[i] ) {
        if (w[h] > c_ || z >= z_ + p[h]) goto step_6;
        z = z_ + p[h];
        for (I k = 1; k <= n; ++k) x[k] = x_[k];
        x[h] = 1;
        if (z == upper_u) goto stop;
        i = h;
        goto step_6;
      } else {
        if (c_ - w[h] < m[h-1]) goto step_6;
        j = h;
        goto step_2;
      }

      stop:
      return;
    }

    // A wrapper around inner_mtu1.
    template<typename W, typename P, typename I>
    void mtu1(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
      // Extra Info Pointer == eip
      mtu1_extra_info_t<W, P, I>* eip = new mtu1_extra_info_t<W, P, I>();
      extra_info_t* upcast_ptr = dynamic_cast<extra_info_t*>(eip);
      sol.extra_info = std::shared_ptr<extra_info_t>(upcast_ptr);

      #ifdef HBM_PROFILE
      // Used to compute the time all the execution algorithm.
      steady_clock::time_point all_mtu1_begin = steady_clock::now();

      // Used by HBM_START_TIMER and HBM_STOP_TIMER if HBM_PROFILE is defined.
      steady_clock::time_point begin;
      #endif

      const W c = eip->c = ukpi.c;
      const I n = eip->n = ukpi.items.size();
      const auto &items = ukpi.items;

      if (!already_sorted) {
        sort_by_eff(ukpi.items);
      }

      // The arrays are of size 'n + 2' because: 1º) notation starts at 1;
      // 2º) inner_mtu1 adds a dummy item at the end.
      myvector<W> w;
      w.resize(n + 2);
      w[0] = 0; // position does not exist, notation begins at 1
      myvector<P> p;
      p.resize(n + 2);
      p[0] = 0; // position does not exist, notation begins at 1
      item_t<W, P> it;
      for (I i = 0, j = 1; i < n; ++i, ++j) {
        it = items[i];
        w[j] = it.w;
        p[j] = it.p;
      }

      myvector<W> x;
      // This resize don't initialize the m contents with zero, inner_mtu1
      // will already clean the array for us.
      x.resize(n+1);

      // Variable z will be overwritten by inner_mtu1 too.
      P z;
      hbm_mtu_impl::inner_mtu1(w, p, n, c, z, x);

      sol.opt = 0;
      sol.y_opt = 0;
      for (I k = 1; k <= n; ++k) {
        if (x[k] > 0) {
          sol.used_items.emplace_back(item_t<W, P>(w[k], p[k]), x[k], k);
          sol.opt += x[k]*p[k];
          sol.y_opt += x[k]*w[k];
        }
      }
      assert(z == sol.opt);
      #ifdef HBM_PROFILE
      eip->total_time = difftime_between_now_and(all_mtu1_begin);
      #endif
    }

    // Sorting method internally used by MTU2. Instead of putting the k most
    // efficient items on the beggining of the vector, it will put all the
    // items with the k-1 greatest efficiencies (i.e. if there's more than one
    // item with the same efficiency, this can be more than k items) and the
    // first item of the k greatest efficiency (if there's k different efficie).
    // The items vector is expected to be indexed base 0 (zero), and to have at
    // least one element.
    // The return is the size of the group made of the items with the k
    // greatest efficiencies (i.e. the number of sorted items).
    template<typename V, typename I>
    inline I mtu2_sort_k_most_efficient(V &items, const I &k) {
      const I n = items.size();
      // As we compute the number of different efficiencies found by comparing
      // an element with the previous one, the first efficiency is never
      // accounted for. We start with one less to account for it.
      I qt_effs_to_find = k - 1;
      I qt_sorted_items = 0;
      // If there's at least k different efficiences on the range, i ends up
      // as the index for the first item with the k-esim efficiency.
      I i = 1;

      // NOTE: The original fortran code make use of a very complex algorithm
      // created by the same author ("A hybrid algorithm for finding the kth
      // smallest of n elements in O(n) time", DOI: 10.1007/BF02288326).
      // As our objective is to sort until the kth smallest element, we use
      // the much simpler and specific code below. This code is optimal if
      // all the items have different efficiencies. If all the items have
      // the same efficiency the loop will execute n/k times; taking in
      // account that k is defined with base on n, this means no more than
      // 100 times.
      while (qt_effs_to_find && qt_sorted_items < n) {
        sort_by_eff(items.begin() + qt_sorted_items, items.end(), k);
        for (I limit = min(n, qt_sorted_items + k); qt_effs_to_find && i < limit; ++i) {
          if (!items[i-1].has_same_eff_that(items[i])) --qt_effs_to_find;
        }
        qt_sorted_items += k;
      }

      // To be completely faithful to MTU2, we need to include the first item
      // with the k-esim greatest efficiency if, and only if, the quantity of
      // items more efficient than the k-esim efficiency is k-1 (i.e. there is
      // no item more efficient than the k-esim efficiency that shares its
      // efficiency with another). Otherwise (at least two items share an
      // efficiency and this efficiency is greater than the k-esim efficiency)
      // we should not include the first item with the k-esim efficiency on the
      // return.
      if (qt_effs_to_find) { // there was less than k different efficiencies
        return i; // i isn't the first item with the k-esim efficiency
      } else {
        // 'i' is the index of the first item with the k-esim efficiency.
        if (i == k) {
          // Between the items more efficient than the k-esim efficiency
          // no two item shared an efficiency.
          return k;
        } else {
          // As we found at least k different efficiencies, i can only be
          // greater than k. And we should not include the first item with
          // the k-esim efficiency.
          return i-1;
        }
      }
    }

    // When the item list is already sorted, we yet need to get to know the
    // index of the k-esim most efficient item. This is used to define the
    // next_slice_size on mtu2.
    // The items vector is expected to be indexed base 0 (zero), and to have at
    // least one element.
    // The return is the size of the group made of the items with the k
    // greatest efficiencies.
    template<typename V, typename I>
    inline I mtu2_get_r_from_sorted(V &items, const I &k) {
      I qt_effs_found = 1;
      I i = 1;

      for (; qt_effs_found < k && i < items.size(); ++i) {
          if (!items[i-1].has_same_eff_that(items[i])) ++qt_effs_found; 
      }

      // See the reasoning for the formulae below on mtu2_sort_k_most_efficient.
      return qt_effs_found < k ? i : (i == k ? k : i-1);
    }

    // MTU2 method. Uses inner_mtu1 internally. Is the same algorithm
    // described on "Knapsack Problems" p. 100 (Martello and Toth)
    // but, as the algorithm described on the book is on a very high-level
    // pseudocode, the translation to optimized C++ can be hard to follow.
    template<typename W, typename P, typename I>
    void mtu2(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
      // Extra Info Pointer == eip
      mtu2_extra_info_t<W, P, I>* eip = new mtu2_extra_info_t<W, P, I>();
      extra_info_t* upcast_ptr = dynamic_cast<extra_info_t*>(eip);
      sol.extra_info = std::shared_ptr<extra_info_t>(upcast_ptr);

      #ifdef HBM_PROFILE
      // Used to compute the time all the execution algorithm.
      steady_clock::time_point all_mtu2_begin = steady_clock::now();

      // Used by HBM_START_TIMER and HBM_STOP_TIMER if HBM_PROFILE is defined.
      steady_clock::time_point begin;
      #endif

      const W c = eip->c = ukpi.c;
      const I n = eip->n = ukpi.items.size();
      auto &items = ukpi.items;

      // This 'if' reduces unnecessary overhead when the instance is a
      // subset-sum instance (all item efficiencies equal, no need to sort),
      // at a low overhead for the other unsorted instances.
      if (is_sorted(items.begin(), items.end())) {
        already_sorted = true;
      }

      P z;
      // v is the size the core problem begins and grows at each iteration.
      // TODO: MAKE THIS A PARAMETER
      if (n <= 100) {
        hbm_mtu_impl::mtu1(ukpi, sol, already_sorted);
        return;
      }
      I v = max(static_cast<I>(100), n/100);

      I next_slice_size;
      if (already_sorted) {
        next_slice_size = mtu2_get_r_from_sorted(items, v);
      } else {
        next_slice_size = mtu2_sort_k_most_efficient(items, v);
      }
      //cout << "next_slice_size: " << next_slice_size << endl;

      // The core will be kept as a profit and a weight array, to avoid having
      // to convert for each call to inner_mtu1. The arrays are of size 'n + 2'
      // because: 1º) notation starts at 1; 2º) inner_mtu1 adds a dummy item at
      // the end.
      myvector<W> core_p;
      core_p.resize(n + 2); // this resize don't initialize any values
      core_p[0] = 0; // position does not exist, notation begins at 1
      myvector<P> core_w;
      core_w.resize(n + 2); // this resize don't initialize any values
      core_w[0] = 0; // position does not exist, notation begins at 1
      item_t<W, P> aux_item;
      for (I i = 0, j = 1; i < next_slice_size; ++i, ++j) {
        aux_item = items[i];
        core_w[j] = aux_item.w;
        core_p[j] = aux_item.p;
      }

      // MTU2 describe some steps on a high level. We tried to implement those
      // on the most efficient fashion possible. One of the high-level steps
      // that are hard to manage in an efficient fashion is keeping non_core
      // updated, as it has elements removed frequently (from front and any
      // position in the middle of it).
      myvector< item_t<W, P> > non_core;
      non_core.resize(n - next_slice_size);
      const size_t qt_bytes_to_copy = non_core.size()*sizeof(item_t<W, P>);
      memcpy(non_core.data(), items.data() + next_slice_size, qt_bytes_to_copy);
      // The vector<bool> isn't ideal here, as it trades memory consumption for
      // time. And we are optimizing this algorithm for time.
      myvector<char> items_to_remove(non_core.size(), 0);
      size_t qt_items_to_remove = 0;
      // This vector is used to assist the removal of items from non_core,
      // between the loops we will copy the elements that aren't marked
      // to be removed from non_core to aux_non_core, and then swap both.
      myvector< item_t<W, P> > aux_non_core;

      // The solution vector. We don't initialize it. The inner_mtu1
      // procedure will alredy do it for us.
      myvector<W> x;
      x.resize(n + 1);

      // Solve for the core problem for the first time.
      hbm_mtu_impl::inner_mtu1(core_w, core_p, next_slice_size, c, z, x);

      P upper_u3 = hbm_mtu_impl::u3(items[0], items[1], items[2], c);

      // non_core_start: When we copy elements from non_core to core we need
      // remove them from non_core, but simply doing so would not be efficient.
      // This way, when we copy elements from non_core to core, we set
      // non_core_start to be the index of the first element that wasn't
      // copied. After we will mark dominated items for removal, and then we
      // can remove both the items that were added to the core and the items
      // that were marked for removal on one sinle pass.
      I non_core_start = 0;

      // non_core_size: to avoid calling non_core.size(). non_core_size
      // is always updated immediatly after non_core changes size.
      size_t non_core_size = non_core.size();

      // k: Index of the next position to be written on core_w/core_p,
      //    takes in consideration that both vectors start at index 1.
      I k = next_slice_size + 1;

      // We stop when our current solution hits the upper bound
      // (z == upper_u3), or when the core problem uses all undominated items
      // (non_core_size <= non_core_start).
      while (z != upper_u3 && non_core_size > non_core_start) {
        // On this for loop we mark for removal the dominated items.
        for (size_t j = non_core_start; j < non_core_size; ++j) {
          const P &pj = non_core[j].p;
          const W &wj = non_core[j].w;
          P u = pj + u1(items[0], items[1], c - wj);
          if (u > z) {
            u = pj + hbm_mtu_impl::u3(items[0], items[1], items[2], c - wj);
          }
          // the 'if' below can't be an 'else' of the 'if' above, don't try!
          // the variable 'u' changes value on the 'if' above
          if (u <= z) {
            items_to_remove[j] = true;
            ++qt_items_to_remove;
          }
        }
        // Ends the loop if all items are marked for removal.
        if (non_core_start + qt_items_to_remove >= non_core_size) break;

        // Here we remove the dominated items from non_core (that were marked
        // for removal on the loop above), and also remove the v first items
        // that were added to core (and that are not part of non_core anymore).
        // Then we update the relevant variables (ex.: non_core_size).
        {
          aux_non_core.resize(non_core_size - non_core_start - qt_items_to_remove);

          for (size_t j = non_core_start, i = 0; j < non_core_size; ++j) {
            if (!items_to_remove[j]) { 
              aux_non_core[i++] = non_core[j];
            }
          }

          non_core.swap(aux_non_core);
          non_core_start = 0;
          non_core_size = non_core.size();
          memset(items_to_remove.data(), 0, non_core_size*sizeof(char));
          qt_items_to_remove = 0;
        }

        // Now we get the v most efficient items on non_core, and add them to
        // the core. We do not remove those items of non_core, but we set
        // non_core_start that will remove them for us on the next loop.
        if (already_sorted) {
          next_slice_size = mtu2_get_r_from_sorted(non_core, v);
        } else {
          next_slice_size = mtu2_sort_k_most_efficient(non_core, v);
        }
        //cout << "next_slice_size: " << next_slice_size << endl;
        for (I i = 0; i < next_slice_size; ++i, ++k) {
          aux_item = non_core[i];
          core_w[k] = aux_item.w;
          core_p[k] = aux_item.p;
        }
        non_core_start = next_slice_size;

        // Finally, we solve the core problem again, using MTU1, and then
        // start the process again.
        // We do not need to clean the x solution array or the z variable,
        // inner_mtu1 already does this for us.
        hbm_mtu_impl::inner_mtu1(core_w, core_p, k-1, c, z, x);
      } 

      // Here we put the solution on the sol (solution_t) structure.
      sol.opt = 0;
      sol.y_opt = 0;
      for (I j = 1; j < k; ++j) {
        if (x[j] > 0) {
          auto item = item_t<W, P>(core_w[j], core_p[j]);
          sol.used_items.emplace_back(item, x[j], j);
          sol.opt += x[j]*core_p[j];
          sol.y_opt += x[j]*core_w[j];
        }
      }

      #ifdef HBM_PROFILE
      eip->total_time = difftime_between_now_and(all_mtu2_begin);
      #endif
    }

    // -------------------- WRAPPERS --------------------
    template<typename W, typename P, typename I>
    struct mtu1_wrap : wrapper_t<W, P, I> {
      virtual void operator()(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) const {
        // Calls the overloaded version with the third argument as a bool
        hbm_mtu_impl::mtu1(ukpi, sol, already_sorted);

        return;
      }

      virtual const std::string& name(void) const {
        static const std::string name = "mtu1";
        return name;
      }
    };

    template<typename W, typename P, typename I = size_t>
    void mtu1(
      instance_t<W, P> &ukpi,
      solution_t<W, P, I> &sol,
      int argc,
      argv_t argv
    ) {
      simple_wrapper(mtu1_wrap<W, P, I>(), ukpi, sol, argc, argv);
    }

    template<typename W, typename P, typename I>
    struct mtu2_wrap : wrapper_t<W, P, I> {
      virtual void operator()(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) const {
        // Calls the overloaded version with the third argument as a bool
        hbm_mtu_impl::mtu2(ukpi, sol, already_sorted);

        return;
      }

      virtual const std::string& name(void) const {
        static const std::string name = "mtu2";
        return name;
      }
    };

    template<typename W, typename P, typename I = size_t>
    void mtu2(
      instance_t<W, P> &ukpi,
      solution_t<W, P, I> &sol,
      int argc,
      argv_t argv
    ) {
      simple_wrapper(mtu2_wrap<W, P, I>(), ukpi, sol, argc, argv);
    }
  }
  // -------------------- EXTERNAL FUNCTIONS --------------------

  /// Gives an upper bound for the optimal profit value of an UKP instance,
  /// based on the most efficient item types of the instance and the knapsack
  /// capacity. The item type efficiency is given by the value of its profit
  /// divided by its weight (p1/w1 >= p2/w2 >= p3/w3 >= pj/wj, j > 3 && j <= n).
  ///
  /// @note Only work for integer weight and profit. Depends on the natural
  ///   rounding-down of integers division.
  /// @param bi Most efficient item type.
  /// @param bi2 Second most efficient item type.
  /// @param bi3 Third most efficient item type.
  /// @param c The knapsack capacity.
  /// @return An upper bound for the optimal profit value of the knapsack
  ///   instance.
  template<typename W, typename P>
  inline P u3(
      const item_t<W, P> &bi,
      const item_t<W, P> &bi2,
      const item_t<W, P> &bi3,
      const W &c
    ) {
    return hbm_mtu_impl::u3(bi, bi2, bi3, c);
  }

  /// Solves an UKP instance by the MTU1 algorithm presented at the book
  /// "Knapsack Problems" (from Martello and Toth) p. 96, and stores the
  /// results at sol.
  ///
  /// @note Ony work for integer weight and profit. Depends on the natural
  ///   rounding-down of integers division.
  /// @param ukpi The UKP instance to be solved.
  /// @param sol The object where the results will be written.
  /// @param already_sorted If the ukpi.items vector needs to be sorted by
  ///   non-decreasing efficiency.
  template<typename W, typename P, typename I>
  void mtu1(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
    hbm_mtu_impl::mtu1(ukpi, sol, already_sorted);
  }

  /// An overloaded function, it's used as argument to test_common functions.
  ///
  /// The only parameter recognized is "--already-sorted". If this parameter is
  /// given the ukpi.items will NOT be sorted by non-increasing efficiency. If
  /// it's ommited the ukpi.items will be sorted by non-increasing efficiency.
  /// If more than one item have the same efficiency, then those items will be
  /// sorted by weight (in relation to each other).
  ///
  /// @see main_take_path
  /// @see mtu2(instance_t<W, P> &, solution_t<W, P, I> &, bool)
  template<typename W, typename P, typename I>
  void mtu1(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
    hbm_mtu_impl::mtu1(ukpi, sol, argc, argv);
  }

  /// Solves an UKP instance by the MTU2 algorithm presented at the book
  /// "Knapsack Problems" (from Martello and Toth) p. 100, and stores the
  /// results at sol.
  ///
  /// @note Ony work for integer weight and profit. Depends on the natural
  ///   rounding-down of integers division.
  /// @param ukpi The UKP instance to be solved.
  /// @param sol The object where the results will be written.
  /// @param already_sorted If the ukpi.items vector needs to be sorted by
  ///   non-decreasing efficiency.
  template<typename W, typename P, typename I>
  void mtu2(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
    hbm_mtu_impl::mtu2(ukpi, sol, already_sorted);
  }

  /// An overloaded function, it's used as argument to test_common functions.
  ///
  /// The only parameter recognized is "--already-sorted". If this parameter is
  /// given the ukpi.items will NOT be sorted by non-increasing efficiency. If
  /// it's ommited the ukpi.items will be sorted by non-increasing efficiency.
  /// If more than one item have the same efficiency, then those items will be
  /// sorted by weight (in relation to each other).
  ///
  /// @see main_take_path
  /// @see mtu2(instance_t<W, P> &, solution_t<W, P, I> &, bool)
  template<typename W, typename P, typename I>
  void mtu2(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
    hbm_mtu_impl::mtu2(ukpi, sol, argc, argv);
  }

  /// The mtu1 code without any wrappers. Used internally by mtu1 and
  /// mtu2. Takes the weights and profits list, the number of items,
  /// the capacity, the variable where to store the optimal value
  /// and the vector where to store the optimal solution found.
  /// Both z and x are cleaned before use, x has to have the correct size
  /// ('n' or greater) but can be unintialized.
  /// Note that this method has index notation beggining on 1 (both
  /// w_ and p_ are expected to begin at 1 and go until n).
  /// It's expected that the size of w and p exceeds n by at least one
  /// position (i.e. w and p have at least size n+2). This is because
  /// inner_mtu1 writes a dummy item at the n+1 position (the data
  /// originally on that position is lost, be warned).
  /// This code assumes that w_ and p_ are sorted by non-increasing
  /// efficiency. Exposed for efficiency reasons.
  template<typename W, typename P, typename I>
  void inner_mtu1(myvector<W> &w_, myvector<P> &p_, const I &n, const W &c, P &z, myvector<W> &x) {
    hbm::hbm_mtu_impl::inner_mtu1(w_, p_, n, c, z, x);
  }
}

#endif //HBM_MTU_HPP

