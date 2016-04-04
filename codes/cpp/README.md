This is the README of the HBM library.

## Sample programs

All .cpp files in the same directory that this README are sample programs. Calling the "make" command in this directory will compile them. The test_<method>.out binaries will try to execute the <method> over all PYAsUKP bencmark instances (they should be in ../../data/ukp/), the run_<method>.out takes an instance filename as first parameter and execute the <method> over it (other additional parameters will be repassed to the method).

You need the boost library (http://www.boost.org/) installed to compile run_per.out, run_ukp5.out and test_ukp5.out, or to use ukp5.hpp and periodicity.hpp in your own programs. If you can't/"don't want" to use boost you will need adapt ukp5 to not use periodicity.hpp (or adapt periodicity.hpp to not need boost) and then remove the dump utility from ukp5.hpp too.

This library tries to be header-only. This means you don't need to compile it and link it against your program, you only need to include the right header to your code and LINK THE BOOST LIBRARY (-lboost_filesystem -lboost_system). The code was written in valid C++11.

## Namespace and macros

  * HBM is the acronym for the project name: Henrique Becker Master's.
  * All macros are prefixed with "HBM_".
  * All C++ guard macros are prefixed by HBM_ and followed by the filename with all the letters in caps and all the periods replaced by underscores (i.e.: ukp5.hpp -> HBM_UKP5_HPP). This should suffice to avoid conflicts between the project files, and allow other people to include the files in their own projects without too much risk. All files are inside codes/cpp so the directory name isn't necessary.
  * Everything is inside the hbm namespace (except the macros, obviously).
  * All the implementation details in a file (except the types) are inside a hbm_<filename without extension>_impl namespace inside hbm. Functions with the same name are defined inside hbm namespace but outside the hbm_*_impl namespace and simply call the function of same name inside the implementation namespace. The implementation namespace is used to allow the use of "using XXX" or "using namespace XXX" without polluting the hbm namespace.

### run_ukp5.out parameters explanation

To avoid replicating information, simply execute the command without any parameters, this will give you the usage. If you can't/"don't want" to do this, search for the function 'usage' inside ukp5.hpp.

### Macros explanation

#### HBM_DUMP -- THIS MACRO DOESNT EXIST ANYMORE, BUT THE COMMAND-LINE OPTION CREATE DUMPS DO THIS
A macro you define if you want to know how the arrays used by UKP5
were after the algorithm execution.
This macro can only be defined if HBM_PROFILE macro is defined. This
happens because the info that is dumped will only exist as a field
of the result_t struct if the HBM_PROFILE macro is defined.
If this macro is defined:
* The boost filesystem library is required. The g++ arguments are
  "-lboost_filesystem -lboost_system". The Makefile already include
  these arguments, if you don't plan to use HBM_DUMP you don't need
  the boost.
* The ukp5 procedure will not do any extra effort. All the extra effort
  is already done by HBM_PROFILE.
* The time that is shown in the standard output will not change, as
  writing the dump files is not measured.
* A "g_dump_<filename>.dat" and a "d_dump_<filename>.dat" files will be
  created at the directory where ukp5 was called. They are simple tables
  with a header "y gy" or "y dy", the first column is the capacity value
  from 0 to C, and the second one is: gy = summed profit values (not
  guaranteed to be optimal), dy = item index numbers (between 0 and N).
  You need to understand the inner workings of the ukp5 algorithm
  to know what those numbers mean. The files size vary with C.

#### HBM_CHECK_PERIODICITY
A macro you define if you want the UKP5 algorithm to execute much faster
over some easy instances and a little slower at all other instances. Some
instances have items that are very small compared to the capacity, and
these items have very different efficiencies (profit/weight ratio).
Optimal solutions of instances with those attributes frequently have a
big quantity of an item type (the most efficient one) and almost none
of the other item types. If this macro is defined, the ukp5 will try to
detect when we can stop computing and simply fill the remaining capacity
with the most efficient item. The consequence is that the computation will
stop earlier for instances with the properties described above, and will
last a litlle longer for other instances (because of the overhead of trying
to detect this situation).

#### HBM_INIT_G_BY_CHUNKS
The UKP5 mainly use two arrays/vectors with a size between c and two times c
(c is the knapsack capacity). One of those arrays (the d array) don't need
initialization. The other (the g array) needs to be zeroed before use.
The g initialization can be done at the algorithm startup (all g is
initialized) or while the algorithm runs (chunks of size w_max are initialized
when needed). Sometimes the capacity is very big, but the instance is 'easy'
and ukp5 (with HBM_CHECK_PERIODICITY enabled) will stop after using only a
small fraction of both vectors. In these cases, enabling this flag
(HBM_INIT_G_BY_CHUNKS) will avoid UKP5 time being bounded by c (even in linear
time). Note that this will not avoid the memory consumption being bounded
by c.

#### HBM_PROFILE_PRECISION
A macro that defines how many digits after the period will be used to display
some percentages outputed by HBM_PROFILE. Its value should be an integral
number.

#### HBM_NOT_DEFINED
A macro that MUST NOT be ever defined. Its function is to allow to comment
large portions of code that have /**/ style commentaries inside then, as
C++ doesn't allow nested comments.

#### HBM_XOR_SWAP
Specializes the swap of items to a XOR swap. See the **Assumption Rules**
rule number 3 below.

#### HBM_INNER_XOR_SWAP(a, b) ((a)^=(b),(b)^=(a),(a)^=(b))
A function macro definition that is used in the inner workings
of the library. Don't redefine it, or depend on its existence.

### Assumption Rules
The code makes some assumptions about the profit and weight types,
if you break those the code probably won't compile or work as
intended. The W template parameter refer to item weight and instance
capacity, the P refer to the item profit, and the I refer to the
item indexes.
1. The weight and the profit types behave as primitive numbers in c++
   would behave (define the same operators, you can output them
   with an ostream, cast from and to them with static_cast, etc...).
2. For UKP5 (dynamic programming) the weight must always be an
   integral number, as it is used to index arrays. If you want to use 
   big numbers (that can't be used to index an array/vector) you will
   need to rewrite the code to change the vector to a min heap (with
   the capacity as a key).
3. If profit or weight don't define the ^= operator, the macro HBM_XOR_SWAP
   MUST NOT be defined (it applies ^= over those datatypes).
4. If profit and weight aren't the same type: You should understand that any
   operation between a value of type profit and one of type weight
   (like multiplication, or division) will convert the weight type
   to the profit type before performing the operation (with
   static_cast<P>(W)).
5. The I template type is used to loop between the items and to store the 
   index of items. The assumptions made about it are: I is integral, I can
   be used to index vectors, I is sufficient large to contain the 'n' value
   (number of items).
   The quantity of items in a solution is of type weight, as it is
   the result of a capacity value divided by a item weight (both capacity
   and item weight are of the weight type). And if the item weight is 1
   (one) then the quantity can be equal to the capacity.

