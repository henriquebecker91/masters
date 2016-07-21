#ifndef HBM_DOMINANCE_HPP
#define HBM_DOMINANCE_HPP

#include <memory> // for shared_ptr
#include <string> // for to_string
#include <list>   // for the no sort dominance

#include "ukp_common.hpp"
#include "wrapper.hpp"
#include "type_name.hpp"

namespace hbm {
  /// Class responsible for pretty printing the result from the wrapped
  /// dominance procedures.
  template <typename W, typename P, typename I>
  struct dom_extra_info_t : extra_info_t {
    std::string info;

    dom_extra_info_t(const std::string &algorithm_name,
                           I num_items, I num_not_sm_dominated_items) {
      info = "algorithm_name: " + algorithm_name + "\n" +
             "type_W: " + hbm::type_name<W>::get() + "\n"
             "type_P: " + hbm::type_name<P>::get() + "\n"
             "type_I: " + hbm::type_name<I>::get() + "\n"
             "Number of items: " + std::to_string(num_items) + "\n" +
             "Number of not dominated items: " +
             std::to_string(num_not_sm_dominated_items) + "\n";
    }

    virtual std::string gen_info(void) {
      return info;
    }
  };
  
  namespace hbm_dominance_impl {
    using namespace std;

    template < typename W, typename P, typename RAI >
    void smdom_sort(RAI begin, RAI end,
                    vector < item_t < W, P > >&undominated) {
      if (begin >= end) return;
      // the vector HAS TO BE sorted by non-increasing profitability first
      // AND non-decreasing weight after for this code work

      undominated.reserve((end-begin) + undominated.size());
      undominated.push_back(*begin);

      // If an item i dominate an item j, then pi/wi >= pj/wj.  (Note that the
      // reverse isn't always true, if pi/wi >= pj/wj then i can dominate j, or
      // not. If the reverse was always true then every UKP problem could be
      // reduced to the best item.) The first statement only means that: if i
      // is the index of an item and the items are ordered by non-increasing
      // efficiency, then it can't be simple or multiple dominated by any item
      // of index j where j < i. This ordering guarantees that the most
      // efficient item is never dominated and that we never need to remove
      // items from our undominated items list.
      for (auto it = begin + 1; it < end; ++it) {
        bool i_is_dominated = false;
        for (auto u : undominated) {
          // The line below only works if we assume that W is an integer
          // number and that when an integer number divides another the
          // result is truncated.
          if (static_cast<P>(it->w/u.w) * u.p >= it->p) {
            i_is_dominated = true;
            break;
          }
        }
        if (!i_is_dominated) undominated.push_back(*it);
      }
    }

    template < typename W, typename P >
    void smdom_sort(vector< item_t <W, P> > &items,
                        vector< item_t <W, P> > &undominated,
                        bool already_sorted) {
      // the vector HAS TO BE sorted by non-increasing profitability first
      // AND non-decreasing weight after for this code work
      if (!already_sorted) sort_by_eff(items);
      auto begin = items.begin(), end = items.end();

      smdom_sort(begin, end, undominated);
    }

    template < typename W, typename P, typename RAI >
    void smdom_no_sort(const RAI begin, const RAI end,
                        vector < item_t<W, P> > &undominated) {
      if (begin >= end) return;
      // ucs -> Undominated CandidateS
      list< item_t<W, P> > ucs;
      ucs.push_back(*begin);

      for (auto it = begin + 1; it < end; ++it) {
        bool i_is_dominated = false;
        // uc -> Undominated Candidate
        // the ucs.end MUST be called at each iteration (we are
        // removing elements from the list in the loop)
        for (auto uc = ucs.begin(); uc != ucs.end(); /*empty, see below*/) {
          // The formulae below only works if we assume that W is an integer
          // number and that when an integer number divides another the
          // result is truncated.
          if (static_cast<P>(it->w/uc->w) * uc->p >= it->p) {
            i_is_dominated = true;
            break;
          } else if (static_cast<P>(uc->w/it->w) * it->p >= uc->p) {
            // The element in the candidate list is dominated
            uc = ucs.erase(uc);
          } else {
            ++uc;
          }
        }
        if (!i_is_dominated) ucs.push_back(*it);
      }

      undominated.reserve(undominated.size() + ucs.size());
      for (auto u : ucs) {
        undominated.push_back(u);
      }
    }

    template<typename W, typename P, typename I>
    void smdom_sort_wrapper(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
      vector< item_t<W, P> > undominated;
      hbm_dominance_impl::smdom_sort(ukpi.items, undominated, already_sorted);

      sol.show_only_extra_info = true;
      dom_extra_info_t<W, P, I>* ptr = new dom_extra_info_t<W, P, I>(
        "smdom_sort", ukpi.items.size(), undominated.size()
      );

      extra_info_t* upcast_ptr = dynamic_cast<extra_info_t*>(ptr);
      sol.extra_info = shared_ptr<extra_info_t>(upcast_ptr);
    }

    template<typename W, typename P, typename I>
    void smdom_no_sort_wrapper(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol) {
      vector< item_t<W, P> > undominated;
      auto begin = ukpi.items.begin(), end = ukpi.items.end();
      hbm_dominance_impl::smdom_no_sort(begin, end, undominated);

      sol.show_only_extra_info = true;
      dom_extra_info_t<W, P, I>* ptr = new dom_extra_info_t<W, P, I>(
        "smdom_no_sort", end - begin, undominated.size()
      );

      extra_info_t* upcast_ptr = dynamic_cast<extra_info_t*>(ptr);
      sol.extra_info = shared_ptr<extra_info_t>(upcast_ptr);
    }

    template<typename W, typename P, typename I>
    struct smdom_sort_wrap : wrapper_t<W, P, I> {
      virtual void operator()(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) const {
        // Calls the overloaded version with the third argument as a bool
        hbm_dominance_impl::smdom_sort_wrapper(ukpi, sol, already_sorted);

        return;
      }

      virtual const std::string& name(void) const {
        static const std::string name = "smdom_sort";
        return name;
      }
    };

    template<typename W, typename P, typename I = size_t>
    void smdom_sort_wrapper(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
      simple_wrapper(smdom_sort_wrap<W, P, I>(), ukpi, sol, argc, argv);
    }

    template<typename W, typename P, typename I = size_t>
    void smdom_no_sort_wrapper(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
      hbm_dominance_impl::smdom_no_sort_wrapper(ukpi, sol);
    }

    template < typename P, typename W, typename I, typename RAI >
    void smdom(const RAI begin, const RAI end,
               vector< item_t < W, P > > &undominated,
               I k,
               bool k_already_sorted) {
      if (!k_already_sorted) sort_by_eff(begin, end, k);
      I size = static_cast<I>(end - begin);
      if (k >= size) {
        hbm_dominance_impl::smdom_sort(begin, end, undominated);
      } else {
        hbm_dominance_impl::smdom_no_sort(begin, end, undominated);
      }
    }

    template < typename P, typename W, typename I >
    void smdom(vector< item_t < W, P > > &items,
               vector< item_t < W, P > > &undominated,
               I k,
               bool k_already_sorted) {
      hbm_dominance_impl::smdom(items.begin(), items.end(), undominated, k, k_already_sorted);
    }
  }

  /// @brief Copy the items that aren't simple or multiple dominated
  ///   from a range to a vector.
  ///
  /// Maybe you want to use shrink_to_fit over undominated after this
  /// procedure. We reserve memory as there wasn't simple or multiple
  /// dominated items in 'items'.
  ///
  /// @param begin Random Access Iterator to item values. The element
  ///   pointed by it is used by the procedure.
  /// @param begin Random Access Iterator to item values. The element
  ///   pointed by it ISN'T used by the procedure. Only the values
  ///   before it are used.
  /// @param undominated The non-dominated items will be added at
  ///   the end of this vector. It doesn't need to be empty, but if
  ///   it isn't, the items already there will be used for the
  ///   dominance computation. To have the guarantee that the result
  ///   is correct all the items in undominated must be at least as 
  ///   efficient as the most efficient item in the 'items' argument,
  ///   and items already in undominated can't dominate other items in
  ///   undominated.
  /// @param already_sorted If is false, the begin~end range is
  ///   modified (will be sorted by efficiency).
  template < typename P, typename W, typename RAI>
  void smdom_sort(RAI begin, RAI end,
                  std::vector < item_t < W, P > >&undominated,
                  bool already_sorted) {
    hbm_dominance_impl::smdom_sort(begin, end, undominated, already_sorted);
  }

  /// @brief Same as the documented overload but don't try to sort the
  ///   source range. Also, the elements from the undominated vector are
  ///   not used on the computation.
  ///
  /// This version has assimptotic complexity of O(N*NND), where
  /// N is the number of items in the range (end-begin) and NND
  /// is the number of items in this range that aren't simple or
  /// multiple dominated. The other version, that sorts the range,
  /// have complexity O(N*NND + N*log(N)), but have a lower constant
  /// factor over N*NND. Because of this, the sorted version is bound
  /// to be better when there are no simple/multiple dominated items
  /// on the range, and the no sort version (this version) is bound to
  /// be better when NND gets close to log(N) (can be even faster than
  /// sorting, if the undominated items are few and the range put them
  /// at beggining).
  ///
  /// @see smdom_sort(RAI begin, RAI end,
  ///   std::vector < item_t < W, P > >&undominated,
  ///   bool already_sorted)
  template < typename P, typename W, typename RAI >
  void smdom_no_sort(const RAI begin, const RAI end,
                     std::vector< item_t < W, P > > &undominated) {
    hbm_dominance_impl::smdom_no_sort(begin, end, undominated);
  }

  /// @brief Copy the items that aren't simple or multiple dominated
  ///   from a range to a vector.
  ///
  /// This version works the following way: If already_sorted is true,
  /// it assumes that the first k elements of the range are ordered by
  /// non-increasing efficiency, and that all elements after the k-1
  /// position are of equal efficiency or less efficient than the
  /// element at the k-1 position. If already_sorted is false, then
  /// the procedure makes the assumption above true. If k is greater
  /// than or equal to the range size, the smdom_sort version is used,
  /// otherwise the smdom_no_sort version is used. Pass k as zero
  /// to denote unambiguously that you want to use the no_sort version
  /// (k_already_sorted can be anything in this case).
  ///
  /// @param begin Random Access Iterator to item values. The element
  ///   pointed by it is used by the procedure.
  /// @param begin Random Access Iterator to item values. The element
  ///   pointed by it ISN'T used by the procedure. Only the values
  ///   before it are used.
  /// @param undominated The non-dominated items will be added at
  ///   the end of this vector. If it is empty, the no problem. If
  ///   it isn't empty read the smdom_sort and smdom_no_sort
  ///   documentation to avoid being surprised by unexpected behavior.
  /// @param sort_k_most_eff A integer between 0 and (end-begin).
  /// @param k_already_sorted If is false, the begin to begin+k range
  ///   is modified (will be sorted by non-increasing efficiency).
  template < typename P, typename W, typename I, typename RAI >
  void smdom(const RAI begin, const RAI end,
             std::vector< item_t < W, P > > &undominated,
             I k,
             bool k_already_sorted) {
    hbm_dominance_impl::smdom(begin, end, undominated, k, k_already_sorted);
  }

  /// Overloaded vector version.
  ///
  /// @see void smdom(const RAI begin, const RAI end, std::vector< item_t < W, P > > &undominated, I k, bool k_already_sorted)
  template < typename P, typename W, typename I >
  void smdom(std::vector< item_t < W, P > > &items,
             std::vector< item_t < W, P > > &undominated,
             I k,
             bool k_already_sorted) {
    hbm_dominance_impl::smdom(items, undominated, k, k_already_sorted);
  }

  /// Execute smdom_sort over ukpi.items. Save result to solution_t.extra_info.
  template<typename W, typename P, typename I>
  void smdom_sort_wrapper(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, bool already_sorted) {
    hbm_dominance_impl::smdom_sort_wrapper(ukpi, sol, already_sorted);
  }

  /// Execute smdom_sort over ukpi.items. Save result to solution_t.extra_info.
  /// If argc < 3, already_sorted = false, else already_sorted = argv[2].
  template<typename W, typename P, typename I = size_t>
  void smdom_sort_wrapper(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
    hbm_dominance_impl::smdom_sort_wrapper(ukpi, sol, argc, argv);
  }

  /// Execute smdom_no_sort over ukpi.items. Save result to
  //solution_t.extra_info.
  template<typename W, typename P, typename I>
  void smdom_no_sort_wrapper(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol) {
    hbm_dominance_impl::smdom_no_sort_wrapper(ukpi, sol);
  }

  /// Execute smdom_no_sort over ukpi.items. Save result to
  //solution_t.extra_info.  / Don't do anything with the parameters.
  template<typename W, typename P, typename I = size_t>
  void smdom_no_sort_wrapper(instance_t<W, P> &ukpi, solution_t<W, P, I> &sol, int argc, argv_t argv) {
    hbm_dominance_impl::smdom_no_sort_wrapper(ukpi, sol, argc, argv);
  }
}

#endif //HBM_DOMINANCE_HPP
