#include <regex> /* for regex, regex_match */
#include <string> /* for stoull */
#include <iostream> /* for cout */

#ifdef INT_EFF
#include <boost/sort/spreadsort/spreadsort.hpp> /* for integer_sort */
#else
#include <algorithm> /* for sort */
#endif

#include "ukp_common.hpp"

using namespace std;
using namespace std::regex_constants;

static const string bs("[[:blank:]]*");
static const string bp("[[:blank:]]+");
static const string nb("([1-9][0-9]*)");

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

void read_ukp_instance(istream &in, ukp_instance_t &ukpi) {
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

  static const regex mline      (mline_s, extended | icase );
  static const regex cline      (cline_s, extended | icase);
  static const regex begin_data (begin_data_s, extended | icase | nosubs);
  static const regex item       (item_s, extended | icase | optimize);
  static const regex end_data   (end_data_s, extended | icase | nosubs);

  smatch what;
  string line;

  get_noncomment_line(in, line, m_exp);
  try_match(mline, line, m_exp, what);

  size_t n = stoull(what[1]);
  ukpi.items.reserve(n);

  get_noncomment_line(in, line, c_exp);
  try_match(cline, line, c_exp, what);
  
  ukpi.c = stoull(what[1]);
  
  get_noncomment_line(in, line, begin_exp);
  try_match(begin_data, line, begin_exp, what);

  for (size_t i = 0; i < n; ++i) {
    get_noncomm_line_in_data(in, line, item_exp);

    if (regex_match(line, what, item)) {
      ukpi.items.emplace_back(stoull(what[1]), stoull(what[2]));
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

void read_sukp_instance(istream &in, ukp_instance_t &ukpi) {
  size_t n;
  in >> n;
  in >> ukpi.c;
  ukpi.items.reserve(n);

  for (size_t i = 0; i < n; ++i) {
    size_t w, p;
    in >> w;
    in >> p;
    ukpi.items.emplace_back(w, p);
  }

  return;
}

void write_sukp_instance(ostream &out, ukp_instance_t &ukpi) {
  size_t n = ukpi.items.size();
  out << n << endl;
  out << ukpi.c << endl;

  for (size_t i = 0; i < n; ++i) {
    item_t tmp = ukpi.items[i];
    out << tmp.w << "\t" << tmp.p << endl;
  }

  return;
}

void sort_by_efficiency(vector<item_t> &items) {
  #ifdef INT_EFF
  boost::sort::spreadsort::integer_sort(items.begin(), items.end());
  #else
  std::sort(items.begin(), items.end());
  #endif
  return;
}

