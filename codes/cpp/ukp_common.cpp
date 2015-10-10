#include <boost/xpressive/xpressive.hpp> /* for sregex::*, regex_match */
#include <string> /* for stoll */
#include <iostream> /* for cout */

#ifdef INT_EFF
#include <boost/sort/spreadsort/spreadsort.hpp> /* for integer_sort */
#else
#include <algorithm> /* for sort */
#endif

#include "ukp_common.hpp"

using namespace std;
using namespace boost::xpressive;

void read_ukp_instance(istream &in, ukp_instance_t &ukpi) {
  static sregex number     = sregex::compile("[1-9][0-9]*", regex_constants::icase);
  static sregex comm       = sregex::compile("[:blank:]*(#.*)?", regex_constants::icase);
  static sregex mline      = sregex::compile("[[:space:]]*n:[[:space:]]*([1-9][0-9]*)", regex_constants::icase);
  static sregex cline      = "[:blank:]*c:[:blank:]*(\\d+)" >> comm;
  static sregex begin_data = "[:blank:]*begin[:blank:]+data" >> comm;
  static sregex item       = "[:blank:]*(\\d+)[:blank:]+(\\d+)" >> comm;
  static sregex end_data   = "[:blank:]*end[:blank:]+data" >> comm;

  bool end_data_not_reached = true;
  smatch what;
  string line;

  getline(in, line);
  if (!in.good()) throw ukp_read_error("Couldn't read the first line");

  while (in.good() && regex_match(line, comm)) getline(in, line);
  if (!in.good()) throw ukp_read_error("Last line read before i/o error: " + line);
  if (!regex_match(line, what, mline)) throw ukp_read_error("First nocomment line isn't 'n: <number>' (" + line + ")");

  size_t n = stoll(what[1]);
  ukpi.items.reserve(n);

  while (in.good() && regex_match(line, comm)) getline(in, line);
  if (!in.good()) throw ukp_read_error("Last line read before i/o error: " + line);
  if (!regex_match(line, what, cline)) throw ukp_read_error("First nocomment line after 'n: <number>' isn't 'c: <number>' (" + line + ")");
  
  ukpi.c = stoll(what[1]);
  
  while (in.good() && regex_match(line, comm)) getline(in, line);
  if (!in.good()) throw ukp_read_error("Last line read before i/o error: " + line);
  if (!regex_match(line, begin_data)) throw ukp_read_error("First nocomment line after 'c: <number>' isn't 'begin data' (" + line + ")");

  for (size_t i = 0; i < n; ++i) {
    getline(in, line);
    if (!in.good()) throw ukp_read_error("Last line read before i/o error: " + line);

    if (regex_match(line, what, item)) {
      ukpi.items.emplace_back(stoll(what[1]), stoll(what[2]));
    } else {
      if (regex_match(line, end_data)) {
        cout << "WARNING: read_ukp_instance: the n value is " << n <<
             " but end_data was found after only " << i << " items"
             << endl;
        end_data_not_reached = false;
        break;
      } else if (regex_match(line, comm)) {
        cout << "WARNING: read_ukp_instance: found a blank line inside data section." << endl;
      } else {
        throw ukp_read_error("Strange line in the middle of data section (" + line + ")");
      }
    }
  }

  while (end_data_not_reached) {
    getline(in, line);
    while (in.good() && regex_match(line, comm)) {
      cout << "WARNING: read_ukp_instance: found a blank line inside data section." << endl;
      getline(in, line);
    }
    if (!in.good()) throw ukp_read_error("Last line read before i/o error: " + line);
    if (regex_match(line, end_data)) {
      end_data_not_reached = false;
    } else {
      if (regex_match(line, item)) {
        cout << "WARNING: read_ukp_instance: the n value is " << n <<
             " but no end_data was found after this number of items, instead " <<
             " more items were found, ignoring those (only read the first " <<
             n << " items)." << endl;
      } else {
        throw ukp_read_error("Strange line in the middle of data section (" + line + ")");
      }
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

