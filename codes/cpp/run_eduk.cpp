//#include "test_common.hpp"
#include <iostream>
#include "eduk.hpp"
#include "eduk.cpp"

int main(int argc, char** argv) {
    ukp_instance_t ukpi;
    ukpi.c = 14;
    ukpi.items.push_back({5, 7});
    ukpi.items.push_back({4, 8});
    ukpi.items.push_back({6, 9});
    ukpi.items.push_back({10, 24});

    cout << "Testing Lazyfy + Consume" << endl;
    auto l_items = Lazyfy(ukpi.items);
    vector<item_t> res;
    consume(&l_items, res);

    if (res == ukpi.items) cout << "Lazyfy . consume == id" << endl;
    else cout << "Something wrong" << endl;
    cout << endl;

    cout << "Testing AddTest" << endl;
    l_items = Lazyfy(ukpi.items);
    res.clear();

    auto addtest = AddTest({5, 1}, &l_items, 10);
    vector<item_t> items2;
    items2.push_back({10, 8});
    items2.push_back({9, 9});
//    items2.push_back({7, 11});
//    items2.push_back({11, 26});

    consume(&addtest, res);

    if (res == items2) cout << "AddTest working" << endl;
    else cout << "Something wrong" << endl;
    cout << endl;

    cout << "Testing Filter" << endl;
    vector<item_t> items3;
    items3.push_back({6, 9});
    items3.push_back({6, 9});
    items3.push_back({6, 6});
    items3.push_back({6, 10});
    items3.push_back({6, 9});
    items3.push_back({6, 11});
    l_items = Lazyfy(items3);
    res.clear();
    items2.clear();

    auto filter = Filter(&l_items);
    items2.push_back({6, 9});
    items2.push_back({6, 10});
    items2.push_back({6, 11});

    consume(&filter, res);

    if (res == items2) cout << "Filter working" << endl;
    else cout << "Something wrong" << endl;
    cout << endl;

    cout << "Testing Merge" << endl;
    items2.clear();
    items3.clear();
    vector<item_t> items4;
    res.clear();

    items4.push_back({10, 24});
    items4.push_back({6, 9});
    items4.push_back({5, 7});
    items4.push_back({4, 8});

    items3.push_back({12, 10});
    items3.push_back({5, 9});
    items3.push_back({3, 5});

    l_items = Lazyfy(items4);
    auto l2_items = Lazyfy(items3);

    auto merge = Merge(&l_items, &l2_items);
    
    items2.push_back({12, 10});
    items2.push_back({10, 24});
    items2.push_back({6, 9});
    items2.push_back({5, 9});
    items2.push_back({4, 8});
    items2.push_back({3, 5});

    consume(&merge, res);

    if (res == items2) cout << "Merge working" << endl;
    else cout << "Something wrong" << endl;
    cout << endl;

    cout << "Testing AddHead" << endl;
    items2.clear();
    items3.clear();
    items4.clear();
    res.clear();

    items4.push_back({10, 24});
    items4.push_back({6, 9});
    items4.push_back({5, 7});
    items4.push_back({4, 8});

    items2.push_back({20, 14});
    items2.push_back({10, 24});
    items2.push_back({6, 9});
    items2.push_back({5, 7});
    items2.push_back({4, 8});

    l_items = Lazyfy(items4);

    auto addhead = AddHead({20, 14}, &l_items);
    
    consume(&addhead, res);

    if (res == items2) cout << "AddHead working" << endl;
    else cout << "Something wrong" << endl;
    cout << endl;


    ukp_solution_t sol;
    eduk(ukpi, sol, false);
    ukp_instance_t debug;
    debug.c = 0;
    debug.items = sol.res;
    write_sukp_instance(cout, debug);
}

