#ifndef HBM_UKP_COMMON_HPP
#define HBM_UKP_COMMON_HPP

// Includes for interface.
#include <vector>    // for vector
#include <iostream>  // for istream and ostream
#include <stdexcept> // for runtime_error
#include <utility>   // to specialize swap

// Includes for implementation.
#include <regex>            // for regex, regex_match
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
    /// We take an stream as a parameter to allow people changing
    /// the precision the numbers should be shown.
    ///
    /// @param out An ostream where the object representation will be
    ///   outputed.
    void print(std::ostream &out = std::cout) const {
      // TODO: create a to_string overloaded funtion for this object
      // and use this function with an stringstream as implementation.
      double e = static_cast<double>(it.p)/static_cast<double>(it.w);

      out << "ix: " << ix << " qt: " << qt << " w: " << it.w << " p: " << it.p << " e: " << e << std::endl;
    }
  };

  /// @brief You should subclass this class if you want to give
  ///   the solution_t more info about your method execution.
  ///
  /// If your method generate extra statistics (time of individual
  /// phases, state of data structures, and other miscellany) but
  /// you don't want to do this inside the method to not affect the
  /// measured time, then you should subclass this class.
  struct extra_info_t {
    /// This method should be redefined to examine your method data
    /// and generate a pretty printing of the relevant stats. Your
    /// method data will be stored in the fields you want in the
    /// subclass.
    virtual std::string gen_info(void) {
      return "";
    }
    /// Your subclass will probably be destroyed by a shared_ptr<extra_info_t>.
    /// So, if you don't redefine the destructor, the base class destructor
    /// will be called instead.
    virtual ~extra_info_t(void) { }
  };

  /// Type that represents a solution of an UKP problem
  ///   (usually the optimal solution).
  template <typename W, typename P, typename I>
  struct solution_t {
    P opt;   ///< The solution value (solution items profit sum).
    W y_opt; ///< The solution weight (solution items weight sum).

    /// The quantities of each item used in the solution.
    std::vector< itemqt_t<W, P, I> > used_items;

    /// @brief If only extra_info should be outputted.
    ///
    /// If this flag is true, then procedures that pretty print solution_t
    /// should ignore all the other fields and print only extra_info.
    /// The default behavior should be print extra_info only after printing
    /// all the other fields.
    bool show_only_extra_info;

    /// This pointer point to any extra data that you may want to 
    /// output. It isn't a string because you may want to compute things
    /// out of the main algorithm, to avoid measuring the time used to
    /// process/format this this extra data.
    std::shared_ptr<extra_info_t> extra_info;

    /// Empty constructor.
    solution_t (void) {
      opt = P(0);
      y_opt = W(0);
      extra_info = std::shared_ptr<extra_info_t>(new extra_info_t());
      show_only_extra_info = false;
    }

    /// "Assemble all parameters in a class." constructor.
    solution_t (P opt, W y_opt,
                const std::vector< itemqt_t<W, P, I> > &used_items,
                bool show_only_extra_info,
                const std::shared_ptr<extra_info_t> &extra_info) : 
                opt(opt), y_opt(y_opt),
                used_items(used_items),
                show_only_extra_info(show_only_extra_info),
                extra_info(extra_info) { }
  };

  /// Exception type for UKP instance read errors.
  struct ukp_read_error : std::runtime_error {
    explicit ukp_read_error (const std::string &s) noexcept : std::runtime_error(s) {};
    explicit ukp_read_error (const char* s) noexcept : runtime_error(s) {};
  };

  /// The argv type, a const pointer to const pointers that point to
  /// const chars.
  typedef const char * const * const argv_t;

  /// Alias for generic ukp solver pointer type. Functions that take an
  /// ukp solving procedure should use this type. Any ukp solving procedure
  /// can easily be wrapped around a procedure with this signature.
  template <typename W, typename P, typename I>
  using ukp_solver_t = void (*)(instance_t<W, P> &,
                                solution_t<W, P, I> &,
                                int argc, argv_t argv);

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

      static const regex mline      (mline_s, icase);
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

    template<typename I, typename RAI>
    void sort_by_eff(RAI begin, RAI end, I sort_k_most_eff) {
      // Without this the most efficient element was put first
      if (sort_k_most_eff == 0) return;

      I size = static_cast<I>(end - begin);
      if (sort_k_most_eff >= size) {
        std::sort(begin, end);
        return;
      }

      RAI middle = begin + sort_k_most_eff;
      std::partial_sort(begin, middle, end);
      return;
    }

    template<typename W, typename P, typename I>
    void sort_by_eff(vector< item_t<W, P> > &items, I sort_k_most_eff) {
      // Without this the most efficient element was put first
      hbm_ukp_common_impl::sort_by_eff(items.begin(), items.end(), sort_k_most_eff);
    }

    template<typename W, typename P, typename I>
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

  /// @brief Sort partially the items by non-increasing efficiency and,
  ///   if tied, by non-decreasing weight.
  ///
  /// Is not guaranteed to be stable. The first sort_k_most_eff items
  /// of the vector will be equal to the sort_k_most_eff items of a
  /// sorted vector. The remaining items (sort_k_most_eff+1 until items.size())
  /// can be in any order (there's no guarantee that they will or will not
  /// be sorted too). If sort_k_most_eff is zero nothing is done, if is
  /// greater than or equal to items.size() all the vector is sorted.
  ///
  /// @param items An item vector that will be changed by the procedure.
  /// @param sort_k_most_eff A integer between 0 and items.size()+1.
  template<typename W, typename P, typename I>
  void sort_by_eff(std::vector< item_t<W, P> > &items, I sort_k_most_eff) {
    hbm_ukp_common_impl::sort_by_eff(items, sort_k_most_eff);
  }

  /// Random Access Iterator overload version.
  template<typename I, typename RAI>
  void sort_by_eff(RAI begin, RAI end, I sort_k_most_eff) {
    hbm_ukp_common_impl::sort_by_eff(begin, end, sort_k_most_eff);
  }

  /// @brief Sort the items by non-increasing efficiency and,
  ///   if tied, by non-decreasing weight.
  ///
  /// Is not guaranteed to be stable.
  ///
  /// @param items An item vector that will be changed by the procedure.
  template<typename W, typename P>
  void sort_by_eff(std::vector< item_t<W, P> > &items) {
    std::sort(items.begin(), items.end());

    return;
  }

  /// @brief Sort the items in the range by non-increasing efficiency
  ///   and if tied, by non-decreasing weight.
  ///
  /// Is not guaranteed to be stable.
  ///
  template<typename RAI>
  void sort_by_eff(RAI begin, RAI end) {
    std::sort(begin, end);
    return;
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
