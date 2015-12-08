#ifndef HBM_UKP_COMMON_HPP
#define HBM_UKP_COMMON_HPP

#include <vector>    // for vector
#include <iostream>  // for istream and ostream
#include <stdexcept> // for runtime_error
#include <utility>   // to specialize swap

// Includes for the implementation code inside hbm::hbm_ukp_common_impl
// namespace. They pollute the global namespace a little, but everything
// is inside their own namespace (std and hbm) so this should
// be an OK thing to do.
#include <regex>            // for regex, regex_match
#include <iostream>         // for cout
#include <algorithm>        // for sort
#include "workarounds.hpp"  // for from_string

/// Namespace that encloses everything about Henrique Becker Master's.
namespace hbm {
  /// The item on an UKP (unbounded knapsack problem).
  ///
  /// @tparam W The weight type. See README for assumptions about this type.
  /// @tparam P The profit type. See README for assumptions about this type.
  template <typename W, typename P>
  struct item_t {
    W w; ///< The public weight field.
    P p; ///< The public profit field.

    /// An empty constructor. Don't initialize anything.
    ///
    /// Useful if you will set the field values after.
    inline item_t(void) {}

    /// The parametrized constructor.
    ///
    /// Simply uses the w and p arguments to initialize the fields w and p.
    ///
    /// @param w Initial weight value.
    /// @param p Initial profit value.
    inline item_t(const W &w, const P &p) : w(w), p(p) {}

    /// Equality test. Depend on W and P defining equality.
    ///
    /// @param o The item at the right of the == operator.
    ///
    /// @return True if the items are equal, False otherwise.
    inline bool operator==(const item_t &o) const {
      return p == o.p && w == o.w;
    }

    /// The standard ordering for items.
    ///
    /// The condition x \< y is true iff the item x is more efficient
    /// than y (x.p/x.w \> y.p/y.w) or, if the efficiencies are equal,
    /// the condition x.w \< y.w is true. 
    ///
    /// @param o The item at the right of the < operator.
    ///
    /// @return True if the criterion described above is met, False otherwise.
    inline bool operator<(const item_t &o) const {
      P a = p * static_cast<P>(o.w),
            b = o.p * static_cast<P>(w);
      return a > b || (a == b && w < o.w); 
    }
  };

  /// Type that represents an instance of the UKP problem.
  template <typename W, typename P>
  struct instance_t {
    /// The capacity of the instance.
    W c;
    /// The items of the instance. Can be in any order.
    std::vector< item_t<W, P> > items;
  };

  /// Auxiliar type used by solution_t.
  ///
  /// This type represents the quantity of an item in a solution.
  /// It aggregates extra info about the item and how to pretty print
  /// this info.
  ///
  /// @tparam I The index type. See README for assumptions about this type.
  template <typename W, typename P, typename I>
  struct itemqt_t {
    item_t<W, P> it; /// The item in question.
    W qt; /// The item quantity (in the optimal result, normally).
    I ix; /// The item index (in the order used by the algorithm, normally).

    /// The parametrized constructor.
    ///
    /// Basically takes every argument and use it to initialize the
    /// matching field.
    ///
    /// @param it The initial value of it.
    /// @param qt The initial value of qt.
    /// @param ix The initial value of ix.
    itemqt_t(const item_t<W, P> &it, const W &qt, const I &ix) : it(it), qt(qt), ix(ix) {}

    /// Write human-readable object representation to a stream.
    ///
    /// @param cout An ostream where the object representation will be
    ///   outputed.
    void print(std::ostream &cout = std::cout) const {
      // The cout below don't have std:: as it refers to the parameter.
      cout << "ix: " << ix << " qt: " << qt << " w: " << it.w << " p: " << it.p << std::endl;
    }
  };

  /// Type that represents a solution of an UKP problem
  ///   (usually the optimal solution).
  ///
  /// @attention If HBM_PROFILE is NOT defined the only fields that exist are
  /// opt, y_opt, used_items and last_y_value_outer_loop (this last only
  /// exists if HBM_CHECK_PERIODICITY is defined).
  template <typename W, typename P, typename I>
  struct solution_t {
    P opt;   ///< The solution value (solution items profit sum).
    W y_opt; ///< The solution weight (solution items weight sum).

    /// The quantities of each item used in the solution.
    std::vector< itemqt_t<W, P, I> > used_items;

    #ifdef HBM_CHECK_PERIODICITY
    /// @attention The last_y_value_outer_loop field exists only if 
    /// HBM_CHECK_PERIODICITY is defined.

    /// The last capacity computed before detecting periodicity and stoping.
    W last_y_value_outer_loop;
    #endif // HBM_CHECK_PERIODICITY

    #ifdef HBM_PROFILE
    // Time of each phase
    double sort_time;        ///< Time used sorting items.
    double vector_alloc_time;///< Time used allocating vectors for DP.
    double linear_comp_time; ///< Time used by linear time preprocessing.
    double phase1_time;      ///< Time used by ukp5 phase1 (find optimal).
    double phase2_time;      ///< Time used by ukp5 phase2 (assemble solution).
    double total_time;       ///< Time used by all previous steps.
    // Some data about instance
    I n;   ///< Instance number of items.
    W c,   ///< Instance capacity.
    w_min, ///< Instance smallest item weight.
    w_max; ///< Instance biggest item weight.
    // Some data about structures manipulated by ukp5
    std::vector<P> g; ///< The vector of size c+w_max+1 and profit values.
    std::vector<I> d; ///< The vector of size c+w_max+1 and item index values.
    // Some statistics
    /// Number of items sized vector, with the quantity of each value i in dy.
    std::vector<W> qt_i_in_dy;
    /// Same as d, but without the positions skipped.
    std::vector<I> non_skipped_d;
    /// Last position of the d vector that wasn't zero or n (number of items).
    W last_dy_non_zero_non_n; 
    /// Quantity of g positions that weren't skipped by ukp5.
    W qt_non_skipped_ys;
    /// Quantity of zeros in g.
    W qt_gy_zeros;
    /// How many times ukp5 phase 1 inner loop executed. Sum of non_skipped_d.
    W qt_inner_loop_executions;
    #endif // HBM_PROFILE
  };

  /// Exception type for UKP instance read errors.
  struct ukp_read_error : std::runtime_error {
    explicit ukp_read_error (const std::string &s) noexcept : std::runtime_error(s) {};
    explicit ukp_read_error (const char* s) noexcept : runtime_error(s) {};
  };

  /// Inner ukp_common implementations. Do not depend.
  namespace hbm_ukp_common_impl {
    using namespace std;
    using namespace std::regex_constants;

    const string bs("[[:blank:]]*");
    const string bp("[[:blank:]]+");
    // If you decided to use the ukp format, but instead of "normal" numbers
    // you decided to use another base or encode they differently (as 
    // infinite precision rationals maybe?) you need to change the regex
    // bellow to something that matches your number encoding format.
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
    void write_ukp_instance(ostream &out, instance_t<W, P> &ukpi) {
      I n = ukpi.items.size();
      out << "n: " << n << endl;
      out << "c: " << ukpi.c << endl;

      out << "begin data" << endl;
      for (I i = 0; i < n; ++i) {
        item_t<W, P> tmp = ukpi.items[i];
        out << tmp.w << "\t" << tmp.p << endl;
      }
      out << "end data" << endl;

      return;
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
      std::sort(items.begin(), items.end());
      return;
    }
  }

  /// Read an instance in the ukp format from the stream to instance_t.
  ///
  /** A description of the ukp format follows:
  \verbatim
  n: <number of itens>
  c: <knapsack capacity>
  begin data
  <weight of first item> <profit of first item>
  <weight of second item> <profit of second item>
  ...                     ...
  <weight of the n item> <profit of the n item>
  end data
  \endverbatim */
  ///
  /// @param in The istream with the instance formatted as described.
  /// @param ukpi The object that will receive the instance.
  /// @exception ukp_read_error If the instance format is wrong.
  template <typename W, typename P, typename I = size_t>
  void read_ukp_instance(std::istream &in, instance_t<W, P> &ukpi) {
    hbm_ukp_common_impl::read_ukp_instance(in, ukpi);
  }

  /// Write an instance_t to a stream in the ukp format.
  ///
  /// See read_ukp_instance for info about this format.
  /// @see read_ukp_instance
  ///
  /// @param out The output stream where the instance will be written.
  /// @param ukpi The instance that will be written.
  template<typename W, typename P, typename I = size_t>
  void write_ukp_instance(std::ostream &out, instance_t<W, P> &ukpi) {
    hbm_ukp_common_impl::write_ukp_instance(out, ukpi);
  }

  /// Read an instance in the sukp format from the stream to instance_t.
  ///
  /** A description of the ukp format follows:
  \verbatim
  <number of itens>
  <knapsack capacity>
  <weight of first item> <profit of first item>
  <weight of second item> <profit of second item>
  ...                     ...
  <weight of the n item> <profit of the n item>
  \endverbatim */
  ///
  /// @param in The istream with the instance formatted as described.
  /// @param ukpi The object that will receive the instance.
  template<typename W, typename P, typename I = size_t>
  void read_sukp_instance(std::istream &in, instance_t<W, P> &ukpi) {
    hbm_ukp_common_impl::read_sukp_instance(in, ukpi);
  }

  /// Write an instance_t to a stream in the sukp format.
  ///
  /// See read_sukp_instance for info about this format.
  /// @see read_sukp_instance
  ///
  /// @param out The output stream where the instance will be written.
  /// @param ukpi The instance that will be written.
  template<typename W, typename P, typename I = size_t>
  void write_sukp_instance(std::ostream &out, instance_t<W, P> &ukpi) {
    hbm_ukp_common_impl::write_sukp_instance(out, ukpi);
  }

  /// @brief Sort the items by non-increasing efficiency and,
  ///   if tied, by non-decreasing weight. Is not guaranteed to
  ///   be stable.
  ///
  /// @param items An item vector that will be changed by the procedure.
  template<typename W, typename P>
  void sort_by_eff(std::vector< item_t<W, P> > &items) {
    hbm_ukp_common_impl::sort_by_eff(items);
  }
}

// Use an optimized swap for the item class (improves sorting time),
// but only if this is possible (all members have the operator ^= defined).
#ifdef HBM_XOR_SWAP
#define HBM_INNER_XOR_SWAP(a, b) ((a)^=(b),(b)^=(a),(a)^=(b))
namespace std {
  template<typename W, typename P>
  inline void swap(hbm::item_t<W, P>& a, hbm::item_t<W, P>& b) noexcept
  {
    HBM_XORSWAP(a.w, b.w);
    HBM_XORSWAP(a.p, b.p);
  }
}
#endif //XOR_SWAP SPECIALIZATION

#endif //HBM_UKP_COMMON_HPP
