THE MAJORITY OF THE CODE INSIDE THIS FOLDER ISN'T OF AUTORSHIP OF HENRIQUE BECKER (repository owner). CHECK THE FILES INDIVIDUALLY FOR ACCESSING AUTORSHIP.

The following files were downloaded from https://people.sc.fsu.edu/~jburkardt/f77_src/knapsack/knapsack.html:
	knapsack.f
	knapsack_prb.f
	knapsack_prb.sh
	knapsack.sh
The two '.sh' scripts were updated by Henrique Becker on 16/05/2016 (for the work of dhis master's thesis). The updated versions are prefixed with 'hbm_', the non-prefixed versions are the original ones. The fortran codes have the original MTU1 and MTU2 codes from Martello and Toth.

The knapsack.sh/hbm_knapsack.sh script used to compile the knapsack.f file needs a f77split binary (whose source is in the f77split.c file, downloaded from http://people.sc.fsu.edu/~jburkardt/c_src/f77split/f77split.html).

The f77split.c file can be compiled by "gcc f77split.c -o f77split" (gcc 5.3.0 used, no warnings shown).

The fortran files need the gfortran compiler to be compiles. On Arch Linux it's provided by the gcc-fortran package.

The file hbm_knapsack.sh was created by Henrique Becker, at 16/05/2016, to substitute knapsack.sh (knapsack.sh needs a directory structure more complex than the really needed).

