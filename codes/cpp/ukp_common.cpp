#include <regex> /* for regex, regex_match */
#include <string> /* for stoull */
#include <iostream> /* for cout */

#ifdef HBM_INT_EFF
#include <boost/sort/spreadsort/spreadsort.hpp> /* for integer_sort */
#else
#include <algorithm> /* for sort */
#endif

#include "ukp_common.hpp"
#include "workarounds.hpp"

using namespace std;
using namespace std::regex_constants;
using namespace hbm;

static const string bs("[[:blank:]]*");
static const string bp("[[:blank:]]+");
/* If you decided to use the ukp format, but instead of "normal" numbers
 * you decided to use another base or encode they differently (as 
 * infinite precision rationals maybe?) you need to change the regex
 * bellow to something that matches your number encoding format. */
static const string nb("([-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?)");

static const string comm_s("[[:blank:]]*(#.*)?");
static const regex  comm  (comm_s, extended | nosubs | optimize);

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

void hbm::read_ukp_instance(istream &in, ukp_instance_t &ukpi) {
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

  quantity n;
  from_string(what[1], n);
  ukpi.items.reserve(n);

  get_noncomment_line(in, line, c_exp);
  try_match(cline, line, c_exp, what);
  
  from_string(what[1], ukpi.c);
  
  get_noncomment_line(in, line, begin_exp);
  try_match(begin_data, line, begin_exp, what);

  for (quantity i = 0; i < n; ++i) {
    get_noncomm_line_in_data(in, line, item_exp);

    if (regex_match(line, what, item)) {
      weight w;
      profit p;
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

void hbm::read_sukp_instance(istream &in, ukp_instance_t &ukpi) {
  quantity n;
  in >> n;
  in >> ukpi.c;
  ukpi.items.reserve(n);

  weight w;
  profit p;
  for (quantity i = 0; i < n; ++i) {
    in >> w;
    in >> p;
    ukpi.items.emplace_back(w, p);
  }

  return;
}

void hbm::write_sukp_instance(ostream &out, ukp_instance_t &ukpi) {
  quantity n = ukpi.items.size();
  out << n << endl;
  out << ukpi.c << endl;

  for (quantity i = 0; i < n; ++i) {
    item_t tmp = ukpi.items[i];
    out << tmp.w << "\t" << tmp.p << endl;
  }

  return;
}

void hbm::sort_by_eff(vector<item_t> &items) {
  #ifdef HBM_INT_EFF
  boost::sort::spreadsort::integer_sort(items.begin(), items.end());
  #else
  std::sort(items.begin(), items.end());
  #endif
  return;
}
