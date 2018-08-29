# 0. Purpose #

Enumerate all the combinatorial types of 2level polytopes in dimension D between 1 and 7


# 1. Dependencies #

The code has three main dependencies
 gcc compiler version >=6.0 or clang
 Boost. For more details, refer to http://www.boost.org 
 Nauty software package, available here http://pallini.di.uniroma1.it. 
Nauty's installation directory has to be passed to the compiler as NAUTYHOME=/path/to/nautydir
Otherwise, the compile will assume that nauty is installed in $(HOME)/.opt/nauty


# 2. Installation #

In order to build the code, just run 'make' in the main folder of this repository.


#3. Usage #

In order to run the code, type a command of the form:
````
./2L_enum D verbose_flag
````
where:
- D is an integer between 1 and 7. For every such D, the algorithm enumerates all combinatorial types of D dimensional 2level polytopes and creates a txt file with all the slack matrices. The algorithm takes a file named (D-1).txt in input, containing the list of all (D-1)dimensional slack matrices of 2level polytopes.
Warning: the computation time for D = 7 is expected to take up to ~61 hours of computation time (on AMD Opteron(TM) 6134 2.3 GHz).
- verbose_flag is an integer between 0 and 3.

verbose_flag = 0 : minimal output (the benchmark times in our paper correspond to this flag)
verbose_flag = 1 : output includes all closed sets + tests for classes of 2L polytopes
verbose_flag = 2 : output includes all closed sets + Hembedding of ground set + tests for classes of 2L polytopes
verbose_flag = 3 : output includes all closed sets, ground sets, slabs_sat_by_point and points_sat_by_slab + tests for classes of 2L polytopes

For instance, starting from the $HOME folder, a correct input is the following:

````
./2L_enum 3 3
````

which passes the values D = 3 and verbose_flag = 3 to the enumeration algorithm.

The user is advised that there might be some deviations in the running times due to the fact that the machine used for the experiments produced in the publication had a different processor. In the case of the current server, the code seems to run slightly faster with respect to the execution times in Table 2. 


# 4. Output of the code #

For every D = 1, ..., 7, the code writes a file (D)d.txt containing the list of all slack matrices of D-dimensional polytopes.

We describe the terminal output of the code for verbose_flag = 3. The other cases are subsets of the case verbose_flag = 3.

 The code takes in input and reads the list of all slack matrices of all (D-1)dimensional 2-level polytopes and it uses nauty to compute the canonical form of each slack matrix.

 At this point, it uses every (D-1)dimensional polytope as base for a D-dimensional polytope.

 It checks the top left (D+1)x(D+1) submatrix corresponds to a simplicial core (actually, since slack matrices are stored as vectors of dynamic_bitsets, the code checks that the top *right* (D+1)x(D+1) submatrix corresponds to a simplicial core. Remember: in dynamic_bitsets, lowest index bits are written to the right).

 It computes the corresponding embedding transformation matrix M_d(0) and its inverse.

 It constructs the Hembedding of facets of the base (used to detect whether points in the ground set have more than two different slacks, being thus incompatible). 
For instance, for the case when the base is the 2simplex, we will have:
````
Constructing Hembedding of facets of the base... 010 100 110 OK
````
which means that the Hembedding of the 2simplex in the hyperplane x_1 = 0 in R^3 is described by the slabs:
````
0 <= x_2 <= 1
0 <= x_3 <= 1
0 <= x_2 + x_3 <= 1
````
It computes the V-embedding and H-embedding of the points of the base P_0 and of the ground set.
It uses the Hembedding of the facets of the base to compute the reduced Hembedding of the ground set.

 It computes the set of all slabs that are valid for the base P_0, as defined in Section 4.1 of the paper.

 It builds the incidence vectors of points vs slabs and of slabs vs points:
   * first for every point x in the reduced ground set X, it prints a dynamic_bitsets (of length #{valid slabs for P_0}) having a 1 in the positions indexed by all slabs E such that 0 <= x(E) <= 1.

   * then, dually, for every slab E valid for P_0 and computed at the previous step, it prints a dynamic_bitsets of length |X| having a 1 in the positions indexed by all points x in the reduced ground set X such that 0 <= x(E) <= 1.

 it constructs the incompatibility matrix, as the matrix I corresponding to the incompatibility relation among pairs of points of the reduced ground set X. It is defined in Section 5.1.

 it runs Ganter's next closure algorithm, starting with A = e_1 (the first canonical unit vector in R^d). For a more detailed account about this algorithm, see: B. Ganter and K. Reuter, Finding all closed sets: A general approach, Order 8 (1991), no. 3, 283–290.

In output we obtain several lines like the following one, for the case when the base is the 2simplex:
````
      | 2level | next_cl      | slackmatrix | 2lev_time
0001  | 1       | 5.295e06  s | 2.0757e05 s | 4.829e05  s 
````

⋅⋅⋅* In the first column, we find a dynamic_bitset which is a characteristic vector for the points (represented as ublas::vector of integers) of the ground set that belong to the closed set A. In the example, "0001" tells us that the first point of the ground set belongs to A. 
In this case, the (Hembedding of the reduced) ground set is 
````
	[3](1,0,0) [3](1,0,1) [3](1,1,1) [3](1,1,0) 
````
Thus, we deduce that only the point `[3](1,0,0)` is in A.

   * The second column, labeled "2level", reports the result of the 2level test for P := conv(P_0 U A):
	'1' means that test is positive
	'0' means that test is negative 
	' ' means that the 2level test is skipped. This happens if the number of vertices of the base P_0 is less than of equal to the number of vertices of any of the facets of the Ddimensional polytope we just constructed, adding the point in the closed set A. In this case, in fact, such polytope will be enumerated further on, when such facet will be assumed as base.

   * the third column "next_cl" reports the time for the computation of the next closed set, performed using Ganter's nextclosure algorithm.

   * the forth column "slackmatrix" reports the time for the construction of the slack matrix of P := conv(P_0 U A) starting from the valid slabs and the point in vert(P_0) U A

   * finally, in the fifth column "2lev_time", we report the time for the 2level test.

 for every base, the algorithm prints the number of closed sets found, the number of 2level polytopes and the total elapsed time for that base. 

 In the end, it prints the total number of closed sets found, the total number of performed 2level tests, the final number of nonisomorphic 2level polytopes and the overall elapsed time.

Moreover, the very last table we compare the numbers for combinatorially inequivalent 2level polytopes and their subclasses:
   * D: dimension in exam
   * Deltaf: 2level polytopes with one simplicial facet
   * CS: centrally symmetric 2level polytopes
   * STAB: stable set polytopes of perfect graphs
   * polar: 2level polytopes whose polar is 2level
   * susp: 2level suspensions, as defined in Section 7
   * (under the "2L" column we print again the number of nonisomorphic 2level polytopes, for comparison)

This data in this last table are reported in Table 1 and Table 3 in the article Enumeration of 2level polytopes.


# 5. Outline of the code #

The sorce code is contained in the file `2L_enum.cpp`.

It is organized in several subfunction. Referring to the output in the previous section, we now describe in which part each of them intervenes.

- `sqsubseteq`, `inc`
are used in the Ganter's algorithm, see: B. Ganter and K. Reuter, Finding all closed sets: A general approach, Order 8 (1991), no. 3, 283–290. 

- `discreteconvexhull_cl`, `incompatibility_cl`
implement the closure operators defined in Sections 4.1, 5.1 respectively

- `hash_matrix`
is used for hashing slack matrices: returns (num_cols  1) * 2^D + (num_rows  1)

- `canonicize`
calls Nauty to obtain the canonical form of the nonincidence graph of slack matrix S

- `checksimplicialcore`, `extractM`, `invertM`
are used in to preprocess the slack matrix of the base polytope P_0

- `construct_slack_matrix`
given the vertices of the base polytope P_0, the points in the closed set A, and the family S of valid slabs for vert(P_0) U A, it constructs the slack matrix for the pair (conv(vert(P_0) U A),S), formally defined in Definition 25.

- `istwolevelpolytope`
checks whether a given 0/1 matrix is the slack matrix of a Ddimensional 2level polytope, by using the list of (D1)dimensional 2level polytopes. This test is purely combinatorial and relies on Lemma 26.

- `is_susp`
tests if the slack matrix in input is a 2level suspension

- `is_polar`
takes a slack matrix in input and tests if the polar of the corresponding polytope is 2level

- `to_list`
takes in input a 0/1 slack matrix S, constructed by the function construct_slack_matrix. Checks if any isomorphic copy of it is already in the list of Ddimensional 2level slack matrices we are incrementally building. If the answer is negative, it adds S to the list. Additionally, if verbose_flag is different from 0, it performs the tests for subclasses of 2level polytopes (Deltaf, CS, STAB, polar, susp  see the previous section of this README file)

- `main`
The main function just uses of all the function above and computes the output described in the previous section of this README file.


# 6. Contacts #

mmacchia@ulb.ac.be