#include <cstdlib>
#include <vector>
#include <istream>

using namespace std;

struct item_t {
	int w;
	int p;
};

struct ukp_instance_t {
	int c;
	vector<item_t> items;
};

struct ukp_solution_t {
	vector<int> g;
	vector<int> d;
};

void ukp5(ukp_instance_t &ukpi, ) {
}

void read_sukp_instance(istream in, ukp_instance_t &ukpi) {
	int n << in;
	ukpi.c << in;
	ukpi.items.reserve(n)

	for (int i = 0; i < n; ++n) {
		item_t tmp;
		tmp.w << in;
		tmp.p << in;
		ukpi.items.pushback(tmp);
	}

	return;
}

/*
void write_ukp_instance(ukp_instance_t &ukp_i, ostream out) {
	return;
}*/

int main (void) {

	return EXIT_SUCCESS;
}

