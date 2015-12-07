#ifndef HBM_UKP_COMMON_HPP
#define HBM_UKP_COMMON_HPP

/* includes that don't depend on type or macro definitions */
#include <vector>
#include <iostream>
#include <stdexcept> /* for runtime_error */

/* type and macro definitions */
/* UKP_TYPES BEGIN */
/* Basic types */
namespace hbm {
  #ifdef HBM_INT_EFF
  typedef weight int_eff;
  #endif
  #ifdef HBM_FP_EFF
  typedef float fp_eff;
  #endif
  #ifdef HBM_RATIONAL_EFF
  typedef boost::rational<weight> rational_eff;
  #endif
}
/* UKP_TYPES END */

/* includes that depend on type or macro definitions */
#if defined(HBM_TWO_MULT_COMP) || defined(HBM_INT_EFF)
#include <utility> /* to specialize swap */
#endif

/* Includes for the implementation code inside hbm::hbm_ukp_common_impl
 * namespace. They pollute the global namespace a little, but everything
 * is inside their own namespace (std and boost) so this should
 * be an OK thing to do. */
#include <regex>    /* for regex, regex_match */
#include <string>   /* for stoull */
#include <iostream> /* for cout */
#ifdef HBM_INT_EFF
  #include <boost/sort/spreadsort/spreadsort.hpp> /* for integer_sort */
#else
  #include <algorithm> /* for sort */
#endif
#include "workarounds.hpp" /* for from_string */

namespace hbm {
  #ifdef NOT_DEFINED
  /* The following preprocessor directives code doesn't allow any two of the
   * following macros to be defined at the same time: HBM_TWO_MULT_COMP,
   * HBM_INT_EFF, HBM_FP_EFF and HBM_RATIONAL_EFF. Also it defines the
   * efficiency type based on the macro value*/
  #ifndef HBM_RATIONAL_EFF
    #ifndef HBM_FP_EFF
      #ifndef HBM_INT_EFF
        #ifndef HBM_TWO_MULT_COMP
          #error ONE OF: HBM_RATIONAL_EFF, HBM_FP_EFF, HBM_INT_EFF OR HBM_TWO_MULT_COMP MUST BE \
                 DEFINED, THE RECOMMENDED IS HBM_TWO_MULT_COMP
        #endif
      #else /*HBM_INT_EFF*/
        #ifdef HBM_TWO_MULT_COMP
          #error HBM_INT_EFF AND HBM_TWO_MULT_COMP CANNOT BE DEFINED AT THE SAME TIME
        #endif
        typedef int_eff efficiency;
      #endif /*HBM_INT_EFF*/
    #else /*HBM_FP_EFF*/
      #ifdef HBM_TWO_MULT_COMP
        #error HBM_FP_EFF AND HBM_TWO_MULT_COMP CANNOT BE DEFINED AT THE SAME TIME
      #endif
      #ifdef HBM_INT_EFF
        #error HBM_FP_EFF AND HBM_INT_EFF CANNOT BE DEFINED AT THE SAME TIME
      #endif
      typedef fp_eff efficiency;
    #endif /*HBM_FP_EFF*/
  #else /*HBM_RATIONAL_EFF*/
    #ifdef HBM_TWO_MULT_COMP
      #error HBM_RATIONAL_EFF AND HBM_TWO_MULT_COMP CANNOT BE DEFINED AT THE SAME TIME
    #endif
    #ifdef HBM_INT_EFF
      #error HBM_RATIONAL_EFF AND HBM_INT_EFF CANNOT BE DEFINED AT THE SAME TIME
    #endif
    #ifdef HBM_FP_EFF
      #error HBM_RATIONAL_EFF AND HBM_FP_EFF CANNOT BE DEFINED AT THE SAME TIME
    #endif
    typedef rational_eff efficiency;
  #endif /*HBM_RATIONAL_EFF*/

  template <typename W, typename>
  struct item_and_eff_t {
    /* Same as item_t, but with an efficiency field. The item_t type was
     * originally this type. The only utility of caching the efficiency
     * was to speed-up the sorting phase. The speed-up was, at max, of 2 times
     * (only taking in account the sorting phase time). Also, the efficiency is
     * innacurate and fast (HBM_INT_EFF or HBM_FP_EFF), or exact and slow
     * (HBM_RATIONAL_EFF). It's also hard to plan for any type that can
     * be used for profit and weight (this was used when the type wasn't
     * a template). There are some instances where sorting takes a significant
     * fraction of the algorithm total time, but for this instances a better
     * solution would be removing profit-dominance in linear time, and/or
     * removing part of the simple dominance using a O(CN) simple dominance 
     * detection algorithm variant. This code can be removed safely, and
     * was there only to this explanation appear at a commit with
     * the last version of this type code.
     */
    W w;
    P p;
    #if defined(HBM_RATIONAL_EFF) || defined(HBM_INT_EFF) || defined(HBM_FP_EFF)
    efficiency eff;
    #endif

    inline item_t(void) {}
    #ifdef HBM_TWO_MULT_COMP
    inline item_t(const W &w, const P &p) : w(w), p(p) {}
    #elif defined(HBM_RATIONAL_EFF)
    inline item_t(const W &w, const P &p) : w(w), p(p), eff(p, w) {}
    #elif defined(HBM_INT_EFF)
    inline item_t(const W &w, const P &p) : w(w), p(p) {
      eff = (p << 32) / w;
    }
    #elif defined(HBM_FP_EFF)
    inline item_t(const W &w, const P &p) : w(w), p(p) {
      eff = static_cast<efficiency>(p) / static_cast<efficiency>(w);
    }
    #endif

    inline bool operator==(const item_t &o) const {
      return p == o.p && w == o.w;
    }

    /* Sort by non-increasing eff, if the efficiences are equal
     * sort by non-decreasing weight
     */
    #if defined(HBM_TWO_MULT_COMP)
    inline bool operator<(const item_t &o) const {
      P a = p * static_cast<P>(o.w),
            b = o.p * static_cast<P>(w);
      return a > b || (a == b && w < o.w); 
    }
    #elif defined(HBM_RATIONAL_EFF) || defined(HBM_FP_EFF) || defined(HBM_INT_EFF)
    inline bool operator<(const item_t &o) const {
      return eff > o.eff || (eff == o.eff && w < o.w);
    }
    #endif

    #ifdef HBM_INT_EFF
    inline efficiency operator>>(const int s) const {
      /* NOTE: this operator is needed for boost::sort::spreadsort::integer_sort
       * that is the sort we use when HBM_INT_EFF is defined.
       * As we sort by non-increasing eff, we cannot return
       * (eff >> s), because the integer_sort will compare the value
       * trying to sort them in non-decreasing order, as (a < b) iff (~a >= ~b)
       * this solves the problem. */
      return (~eff) >> s;
    }
    #endif
  };
  #endif //NOT_DEFINED

  template <typename W, typename P>
  struct item_t {
    W w;
    P p;

    inline item_t(void) {}
    inline item_t(const W &w, const P &p) : w(w), p(p) {}

    inline bool operator==(const item_t &o) const {
      return p == o.p && w == o.w;
    }

    /* Sort by non-increasing eff, if the efficiences are equal
     * sort by non-decreasing weight
     */
    inline bool operator<(const item_t &o) const {
      P a = p * static_cast<P>(o.w),
            b = o.p * static_cast<P>(w);
      return a > b || (a == b && w < o.w); 
    }
  };

  template <typename W, typename P>
  struct instance_t {
    W c;
    std::vector< item_t<W, P> > items;
  };

  template <typename W, typename P, typename I>
  struct itemqt_t {
    item_t<W, P> it; /* the item */
    W qt; /* its quantity in the result */
    I ix; /* its index in the ordered items list */

    itemqt_t(const item_t<W, P> &it, const W &qt, const I &ix) : it(it), qt(qt), ix(ix) {}

    void print(std::ostream &cout = std::cout) const {
      cout << "ix: " << ix << " qt: " << qt << " w: " << it.w << " p: " << it.p << std::endl;
    }
  };

  template <typename W, typename P, typename I>
  struct solution_t {
    P opt;
    W y_opt;
    std::vector< itemqt_t<W, P, I> > used_items;
    #ifdef HBM_PROFILE
    /* Time of each phase */
    double sort_time;
    double vector_alloc_time;
    double linear_comp_time;
    double phase1_time;
    double phase2_time;
    double total_time;
    /* Some data about instance */
    I n;
    W c, w_min, w_max;
    /* Some data about structures manipulates by ukp5 */
    W last_dy_non_zero_non_n;
    W qt_non_skipped_ys;
    W qt_gy_zeros;
    W qt_inner_loop_executions;
    std::vector<I> qt_i_in_dy;
    std::vector<P> g;
    std::vector<I> d;
    std::vector<I> non_skipped_d;
    #endif /* HBM_PROFILE */
    #ifdef HBM_CHECK_PERIODICITY
    W last_y_value_outer_loop;
    #endif /* HBM_CHECK_PERIODICITY */
  };

  struct ukp_read_error : std::runtime_error {
    explicit ukp_read_error (const std::string &s) noexcept : std::runtime_error(s) {};
    explicit ukp_read_error (const char* s) noexcept : runtime_error(s) {};
  };

  namespace hbm_ukp_common_impl {
    using namespace std;
    using namespace std::regex_constants;

    const string bs("[[:blank:]]*");
    const string bp("[[:blank:]]+");
    /* If you decided to use the ukp format, but instead of "normal" numbers
     * you decided to use another base or encode they differently (as 
     * infinite precision rationals maybe?) you need to change the regex
     * bellow to something that matches your number encoding format. */
    const string nb("([-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?)");

    const string comm_s("[[:blank:]]*(#.*)?");
    const regex  comm  (comm_s, extended | nosubs | optimize);

    void warn_about_premature_end(size_t n, size_t i) {
      cerr << "WARNING: read_ukp_instance: the n value is " << n <<
              " but 'end data' was found after only " << i << " items"
              << endl;

      return;
    }

    void warn_about_no_end(size_t n) {
      cerr  << "WARNING: read_ukp_instance: the n value is " << n <<
            " but no 'end data' was found after this number of items, instead" <<
            " more items were found, ignoring those (only read the first " <<
            n << " items)." << endl;

      return;
    }

    void io_guard(const istream &in, const string &expected) {
      if (!in.good()) {
        throw ukp_read_error("Expected '" + expected + "' but file ended first.");
      }

      return;
    }

    void get_noncomment_line(istream &in, string &line, const string &exp) {
      do getline(in, line); while (in.good() && regex_match(line, comm));
      io_guard(in, exp);

      return;
    }

    void get_noncomm_line_in_data(istream &in, string &line, const string &exp) {
        static const string warn_blank_in_data =
          "WARNING: read_ukp_instance: found a blank line inside data section.";

        getline(in, line);
        while (in.good() && regex_match(line, comm)) {
          cout << warn_blank_in_data << endl;
          getline(in, line);
        }
        io_guard(in, exp);

        return;
    }

    void try_match(const regex &r, const string &l, const string &exp, smatch &m) {
      if (!regex_match(l, m, r)) {
        throw ukp_read_error("Expected '" + exp + "' but found '" + l + "'");
      }
      
      return;
    }

    template<typename W, typename P, typename I = size_t>
    void read_ukp_instance(istream &in, instance_t<W, P> &ukpi) {
      static const string strange_line =
        "error: read_ukp_instance: Found a strange line inside data section: ";

      static const string m_exp     = "n/m: <number>";
      static const string c_exp     = "c: <number>";
      static const string begin_exp = "begin data";
      static const string item_exp  = "<number> <number>";
      static const string end_exp   = "end data";

      static const string mline_s       = bs + "[mn]:" + bs + nb + comm_s;
      static const string cline_s       = bs + "c:" + bs + nb + comm_s;
      static const string begin_data_s  = bs + "begin" + bp + "data" + comm_s;
      static const string item_s        = bs + nb + bp + nb + comm_s;
      static const string end_data_s    = bs + "end" + bp + "data" + comm_s;

      static const regex mline      (mline_s, icase );
      static const regex cline      (cline_s, icase);
      static const regex begin_data (begin_data_s, extended | icase | nosubs);
      static const regex item       (item_s, icase | optimize);
      static const regex end_data   (end_data_s, extended | icase | nosubs);

      smatch what;
      string line;

      get_noncomment_line(in, line, m_exp);
      try_match(mline, line, m_exp, what);

      I n;
      from_string(what[1], n);
      ukpi.items.reserve(n);

      get_noncomment_line(in, line, c_exp);
      try_match(cline, line, c_exp, what);
      
      from_string(what[1], ukpi.c);
      
      get_noncomment_line(in, line, begin_exp);
      try_match(begin_data, line, begin_exp, what);

      for (I i = 0; i < n; ++i) {
        get_noncomm_line_in_data(in, line, item_exp);

        if (regex_match(line, what, item)) {
          W w;
          P p;
          from_string(what[1], w);
          from_string(what[2], p);
          ukpi.items.emplace_back(w, p);
        } else {
          if (regex_match(line, end_data)) {
            warn_about_premature_end(n, i);
            break;
          } else {
            throw ukp_read_error(strange_line + "(" + line + ")");
          }
        }
      }

      get_noncomm_line_in_data(in, line, end_exp);

      if (!regex_match(line, end_data)) {
        if (regex_match(line, item)) {
          warn_about_no_end(n);
        } else {
          throw ukp_read_error(strange_line + "(" + line + ")");
        }
      }
    }

    template<typename W, typename P, typename I = size_t>
    void read_sukp_instance(istream &in, instance_t<W, P> &ukpi) {
      I n;
      in >> n;
      in >> ukpi.c;
      ukpi.items.reserve(n);

      W w;
      P p;
      for (I i = 0; i < n; ++i) {
        in >> w;
        in >> p;
        ukpi.items.emplace_back(w, p);
      }

      return;
    }

    template<typename W, typename P, typename I = size_t>
    void write_sukp_instance(ostream &out, instance_t<W, P> &ukpi) {
      I n = ukpi.items.size();
      out << n << endl;
      out << ukpi.c << endl;

      for (I i = 0; i < n; ++i) {
        item_t<W, P> tmp = ukpi.items[i];
        out << tmp.w << "\t" << tmp.p << endl;
      }

      return;
    }

    template<typename W, typename P>
    void sort_by_eff(vector< item_t<W, P> > &items) {
      #ifdef HBM_INT_EFF
      boost::sort::spreadsort::integer_sort(items.begin(), items.end());
      #else
      std::sort(items.begin(), items.end());
      #endif
      return;
    }
  }

  template <typename W, typename P, typename I = size_t>
  void read_ukp_instance(std::istream &in, instance_t<W, P> &ukpi) {
    hbm_ukp_common_impl::read_ukp_instance(in, ukpi);
  }

  template<typename W, typename P, typename I = size_t>
  void read_sukp_instance(std::istream &in, instance_t<W, P> &ukpi) {
    hbm_ukp_common_impl::read_sukp_instance(in, ukpi);
  }

  template<typename W, typename P, typename I = size_t>
  void write_sukp_instance(std::ostream &out, instance_t<W, P> &ukpi) {
    hbm_ukp_common_impl::write_sukp_instance(out, ukpi);
  }

  template<typename W, typename P>
  void sort_by_eff(std::vector< item_t<W, P> > &items) {
    hbm_ukp_common_impl::sort_by_eff(items);
  }
}

/* Use an optimized swap for the item class (improves sorting time),
 * but only if this is possible (all members have the operator ^= defined) */
#if (defined(HBM_TWO_MULT_COMP) || defined(HBM_INT_EFF)) && !defined(HBM_NO_XOR_SWAP)
#define HBM_XORSWAP(a, b) ((a)^=(b),(b)^=(a),(a)^=(b))
namespace std {
  template<typename W, typename P>
  inline void swap(hbm::item_t<W, P>& a, hbm::item_t<W, P>& b) noexcept
  {
    HBM_XORSWAP(a.w, b.w);
    HBM_XORSWAP(a.p, b.p);
    #ifdef HBM_INT_EFF
    HBM_XORSWAP(a.eff, b.eff);
    #endif
  }
}
#endif //XOR_SWAP SPECIALIZATION

#endif //HBM_UKP_COMMON_HPP
