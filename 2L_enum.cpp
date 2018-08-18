//
//  2L_enum.cpp
//
//  Created by Samuel Fiorini and Marco Macchia
//  Purpose: Enumerate all the combinatorial types of 2-level polytopes in dimension D between 3 and 7
//
// REMARKS ON INSTALLATION
// - Install Boost. For more details, refer to http://www.boost.org
// - Compile nauty library, available here http://pallini.di.uniroma1.it
// - The source code here assumes that the header file nauty.h and the following object files are in the source code folder:
//   naurng.o, schreier.o, naugraph.o, nautiL.o, nauty1.o
// - Now it is possible to compile this code with
//      g++ -g 2L_enum.cpp nautyL.o naurng.o nautil.o schreier.o naugraph.o -o 2L_enum

//
// Info on data structures:
//
// - The 0/1-matrices we use are *nonincidence* matrices, which means that 0 indicates
//   incidence while 1 indicates nonincidence. This affects the way sets are manipulated.
//

#include "nauty.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>
#include <algorithm>
#include <chrono>

#define BOOST_UBLAS_NDEBUG 1
//#include <boost/serialization/array_wrapper.hpp>
#include <boost/dynamic_bitset.hpp> // Use Boost's dynamic_bitsets for sets
#include <boost/numeric/ublas/vector.hpp> // Use Boost's basic linear algebra (ublas) vectors
#include <boost/numeric/ublas/matrix.hpp> // Use Boost's basic linear algebra (ublas) matrices
#include <boost/numeric/ublas/io.hpp> // Use Boost's basic linear algebra (ublas) io routines
#include <boost/numeric/ublas/lu.hpp>

// Change the value for D to change the dimension in which the enumeration is performed
// Dimension D, insert a value berween 3 and 7
//#define D 3

// verbose = 0 : minimal output (the benchmarks times correspond to verbose = 0)
// verbose = 1 : output includes all closed sets + tests for classes of 2L polytopes
// verbose = 2 : output includes all closed sets + tests for classes of 2L polytopes
// verbose = 3 : output includes all closed sets, ground sets, slabs_sat_by_point and points_sat_by_slab
//               + tests for classes of 2L polytopes
//#define verbose 0

namespace ublas = boost::numeric::ublas;

typedef boost::dynamic_bitset<>::size_type size_type;
size_type npos = boost::dynamic_bitset<>::npos;

// Canonical forms of all nonincidence graphs of (D-1)-dim 2L-poly
// The key value is the hash value of the matrix (for the moment, computed from number of columns and rows)
// The mapped value is the canonical form of the graph as computed by Nauty, which is a vector of setwords
std::multimap<int,std::vector<setword>> atoms_cg;

// List LD of D-dimensional 2-level polytopes
std::multimap<int,std::vector<setword>> LD;

// Counters for the classes of 2L polytopes
int simplicial_facet = 0; // 2L polytopes with a simplicial facet
int cs = 0; // Centrally-symmetric 2L polytopes
int stab = 0; // Stable set polytopes of a perfect graph
int n_suspensions = 0; // 2L suspensions
int n_polar = 0; // Polar 2L polytopes

// output file
std::ofstream myfile;


// Vector of bitsets, giving for each point the set of slabs that are satisfied by it
std::vector<boost::dynamic_bitset<>> slabs_sat_by_point;

// Vector of bitsets, giving for each slab the set of points (in the ground set) that are satisfied by it
std::vector<boost::dynamic_bitset<>> points_sat_by_slab;

// Check if A \sqsubseteq B
bool sqsubseteq(boost::dynamic_bitset<> A, boost::dynamic_bitset<> B) {
    assert (A.size() == B.size());
    
    // Check if A is a subset of B
    if (A.is_subset_of(B)) {
        
        // find first bit set to 1 in B
        // Skip all bits that are set to 0 in A and to 1 in B, and delete them from B
        // until a bit is found that is 1 in A and 1 in B. At that point, compare A to B
        for (size_type pos = B.find_first(); pos != npos && pos < A.find_first(); pos = B.find_next(pos))
            B.flip(pos);
        
        return (A == B);
    }
    else
        return false;
}

boost::dynamic_bitset<> inc(boost::dynamic_bitset<> A, int i) {
    assert (i < A.size());
    
    for (size_type pos=A.find_first(); pos != npos && pos<i; pos=A.find_next(pos))
        A.flip(pos);
    
    A.set(i);
    
    return A;
}

// Compute the discrete convex hull of a point set A
boost::dynamic_bitset<> discreteconvexhull_cl(const boost::dynamic_bitset<> &A, boost::dynamic_bitset<> &B) {
    assert (A.size() == slabs_sat_by_point.size());
    
    //std::cout << "A.size() = " << A.size() << std::endl;
    //std::cout << "points_sat_by_slab.size() = " << points_sat_by_slab.size() << std::endl;
    //std::cout << "slabs_sat_by_point.size() = " << slabs_sat_by_point.size() << std::endl;
    
    //boost::dynamic_bitset<> B(points_sat_by_slab.size());
    boost::dynamic_bitset<> C(slabs_sat_by_point.size());
    
    size_type pos;
    
    // Set all bits of B to 1
    B.set();
    
    // Intersect all sets of slabs belonging to elements of A
    for (pos = A.find_first(); pos != npos; pos = A.find_next(pos))
        B &= slabs_sat_by_point[pos];
    
    //std::cout << "A = " << A << std::endl;
    //std::cout << "B = " << B << std::endl;
    
    
    // Set all bits of C to 1
    C.set();
    
    // Intersect all sets of points belonging to elements of B
    for (pos = B.find_first(); pos != npos; pos = B.find_next(pos))
        C &= points_sat_by_slab[pos];
    
    //std::cout << "C = " << C << std::endl;
    return C;
}

// incompatibility closure operator
boost::dynamic_bitset<> incompatibility_cl(const std::vector<boost::dynamic_bitset<>> &incompatibility_adjM,boost::dynamic_bitset<> &A) {
    
    if (A.count() > 1) { // if A has less than 2 points, then A = {e_1}
        bool accept = true;
        
        for (size_type i = A.find_next(0); i != npos && accept; i = A.find_next(i)) {
            for (size_type j = A.find_next(0);j != i && accept; j = A.find_next(j))
                accept = !(incompatibility_adjM[i].test(j));
        }
        if (accept)
            return A;
        else {
            boost::dynamic_bitset<> incclA(A.size());
            return incclA.set();
        }
    }
    else
        return A;
}

// Used for hashing slack matrices: returns (num_cols - 1) * 2^D + (num_rows - 1)
// Assert that S has at least one row and at least one column
int hash_matrix(const std::vector<boost::dynamic_bitset<>> &S, int D) {
    // Make sure that the slack matrix has at least one row and column
    assert (S.size() > 0);
    assert (S[0].size() > 0);
    
    int num_rows = (int) S.size();
    int num_cols = (int) S[0].size();
    
    return (((num_cols - 1) << D) + num_rows - 1);
}

// Call Nauty to obtain the canonical form of the nonincidence graph of slack matrix S
std::vector<setword> canonicize(const std::vector<boost::dynamic_bitset<>> &S) {
    std::vector<setword> cg_vec;
    
    // Make sure that the slack matrix has at least one row and column
    assert (S.size() > 0);
    assert (S[0].size() > 0);
    
    // Get the number of rows and columns of current atom
    int num_rows = (int) S.size();
    int num_cols = (int) S[0].size();
    
    // Initializations for Nauty
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    DYNALLSTAT(int,map,map_sz);
    DYNALLSTAT(graph,g,g_sz);
    DYNALLSTAT(graph,cg,cg_sz);
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    
    // Select option for canonical labelling
    options.getcanon = TRUE;
    
    // Select option for custom partition
    options.defaultptn = FALSE;
    
    int m, n;
    
    n = num_rows + num_cols;

    
    m = SETWORDSNEEDED(n);
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
    
    
    // Allocate memory for the graph
    DYNALLOC1(int,lab,lab_sz,n,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
    DYNALLOC1(int,map,map_sz,n,"malloc");
    DYNALLOC2(graph,g,g_sz,n,m,"malloc");
    DYNALLOC2(graph,cg,cg_sz,n,m,"malloc");
    
    // Empty the graph
    EMPTYGRAPH(g,m,n);
    
    // Build the custom partition for the graph
    for (int i = 0; i < num_rows+num_cols; i++) {
        lab[i] = i;
        if (i != num_rows-1 && i != num_rows+num_cols-1)
            ptn[i] = 1;
        else
            ptn[i] = 0;
    }
    
    // Build the edges of the nonincidence graph
    
    // Loop through all the entries of the slack matrix and add an edge when there is a one
    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            if ((S[i])[j] == 1)
                ADDONEEDGE(g,i,j+num_rows,m);
        }
    }
    
    // Obtain canonical graph from Nauty
    densenauty(g,lab,ptn,orbits,&options,&stats,m,n,cg);
    
    // Make std::vector<setword> from canonical graph
    for (size_t k = 0; k < m*(size_t)n; ++k)
        cg_vec.push_back(cg[k]);
    
    // Clean up
    DYNFREE(lab,lab_sz);
    DYNFREE(ptn,ptn_sz);
    DYNFREE(orbits,orbits_sz);
    DYNFREE(map,map_sz);
    DYNFREE(g,g_sz);
    DYNFREE(cg,cg_sz);
    
    return cg_vec;
}

// Checks whether a given 0-1 matrix has a (d+1) x (d+1) lower triangular nonsingular submatrix
// at the top left corner
bool checksimplicialcore(const std::vector<boost::dynamic_bitset<>> &S, int d) {
    bool test = true;
    
    // Check that S has at least d rows
    assert (S.size() >= d);
    
    for (int i=0; i<=d && test; ++i) {
        // Check that the ith row of S has at least d columns
        assert (S[i].size() >= d);
        for (int j=i; j<=d && test; ++j) {
            if (i == j)
                test = ((S[i])[j] == 1);
            else if (i < j)
                test = ((S[i])[j] == 0);
        }
    }
    
    return test;
}

// Extract extended embedding transformation matrix M_d(0) from simplicial core
ublas::matrix<int> extractM(std::vector<boost::dynamic_bitset<>> S, int d) {
    ublas::matrix<int> M(d,d);
    
    // Check that S has at least d rows
    assert (S.size() >= d);
    
    for (int i=0; i<d; ++i) {
        // Check that the ith row of S has at least d columns
        assert (S[i].size() >= d);
        for (int j=0; j<d; ++j) {
            if (i == 0 && j == 0)
                M(i,j) = 1;
            else if (i == 0 || j == 0)
                M(i,j) = 0;
            else
                M(i,j) = (S[i-1])[j-1];
        }
    }
    
    return M;
}


bool invertM(const ublas::matrix<int> &M, ublas::matrix<int> &Minv) {
    // Create a duplicate of the input matrix
    ublas::matrix<size_type> W = M;
    
    assert (W.size1() == W.size2());
    size_type d = W.size1();
    
    // Create the permutation matrix
    ublas::permutation_matrix<size_type> P(d);
    
    Minv.resize(d,d);
    
    // Assign the identity matrix to the inverse
    Minv.assign(ublas::identity_matrix<int>(d));
    
    //  LU factorization and substitution
    size_type res = lu_factorize(W, P);
    if (res != 0)
        return false;
    lu_substitute(W,P,Minv);
    return true;
}

// slack matrix construction
std::vector<boost::dynamic_bitset<>> construct_slack_matrix(const std::vector<ublas::vector<int>> &base_H,const std::vector<ublas::vector<int>> &ground_set_H,const boost::dynamic_bitset<> &A,const boost::dynamic_bitset<> &B,const std::vector<ublas::vector<int>> &slabs, const std::vector<boost::dynamic_bitset<>> &S, std::vector<boost::dynamic_bitset<>> &S_new, int &D) {
    std::vector<boost::dynamic_bitset<>> all_ineqs;

    for (size_type i = B.find_first(); i != npos; i = B.find_next(i)){
        
        // constuct each row of S with column indexed by e_1, the vertices of the base and then the remaining ones
        // in this way the first D+1 columns will be the H-embedding of the canonical affine base of R^D.
        boost::dynamic_bitset<> S_row;
        
        int s = ublas::inner_prod(ground_set_H[0],slabs[i]);
        bool slack1 = (s != 0);
        S_row.push_back(slack1);
        
        for (std::vector<ublas::vector<int>>::const_iterator it=base_H.begin(); it!=base_H.end(); it++) {
            int s = ublas::inner_prod(*it,slabs[i]);
            bool slack1 = (s != 0);
            S_row.push_back(slack1);
        }
        for (size_type j = A.find_next(0); j != npos; j = A.find_next(j)){
            int s = ublas::inner_prod(ground_set_H[j],slabs[i]);
            bool slack1 = (s != 0);
            S_row.push_back(slack1);
        }

        if ((~S_row).count() >= D)
            all_ineqs.push_back(S_row);
        if (S_row.count() >= D)
            all_ineqs.push_back(~S_row);
    }
    
    
    // check maximality of rows
    std::vector<boost::dynamic_bitset<>>::iterator it1,it2;
    for (it1=all_ineqs.begin(); it1!=all_ineqs.end(); ++it1) {
        bool is_maximal = true;
        for (it2=all_ineqs.begin(); it2!=all_ineqs.end() && is_maximal; ++it2) {
            if ((*it2).is_subset_of(*it1) && it2!=it1)
                is_maximal = false;
        }
        if (is_maximal)
            S_new.push_back(*it1);
    }
    
    // rearranging rows of S
    int n_row = 0;
    for (std::vector<boost::dynamic_bitset<>>::iterator it = S_new.begin()+n_row+1; it != S_new.end(); it++) {
        //for (std::vector<boost::dynamic_bitset<>>::iterator print = S.begin(); print != S.end(); print++)
        //    std::cout << *print << std::endl;
        //std::cout << n_row + 1 << std::endl;
        bool accept = true;
        for (size_type j = 0; j < S[0].size() && n_row < D && accept; j++)
            accept = ((*it).test(j+1) == (S[n_row]).test(j));
        
        if (accept) {
            swap(*it, S_new[n_row+1]);
            ++n_row;
        }
    }
    /*
     // debugging
     std::cout << S_new[0].size() << ',' << S_new.size() << ' ';
     std::cout << std::endl;
     for (std::vector<boost::dynamic_bitset<>>::iterator print = S_new.begin(); print != S_new.end(); print++)
     std::cout << *print << std::endl;
     std::cout << std::endl;*/
    return S_new;
}


// Checks whether a given 0-1 matrix is the slack matrix of a D-dimensional 2-level polytope,
// by using the list of (D-1)-dimensional 2-level polytopes.
bool istwolevelpolytope(std::vector<boost::dynamic_bitset<>> S, int &D) {
    
    // First test: check that every column contains at least D zeros
    // by construction, every row of S contains at least D zeros
    bool accept = true;
    for (size_type it_cols = 0; it_cols != S[0].size() && accept; it_cols++) {
        int n_zeros_col = 0;
        for (std::vector<boost::dynamic_bitset<>>::iterator it_rows = S.begin(); it_rows != S.end(); it_rows++) {
            n_zeros_col += !(*it_rows).test(it_cols);
        }
        accept = (n_zeros_col >= D);
    }
    if (!accept)
        return false;
    
    std::vector<boost::dynamic_bitset<>>::iterator it1, it2, it3; // iterators to browse rows of matrix
    std::vector<std::vector<int>> neighbors; // Collect the indices of neighbors of each row
    bool notadjacent;
    
    // Initialize the neighbors data structure with S.size() empty lists.
    std::vector<int> empty_neighborhood;
    for (it1=S.begin(); it1!=S.end(); ++it1)
        neighbors.push_back(empty_neighborhood);
    
    // Dual graph computation
    // std::cout << "Dual graph computation..." << std::endl;
    
    // This loop computes the dual graph of S by going through all unordered pairs
    // of rows of S (indexed by it1 and it2), and checking whether the union does
    // not contain a third row (indexed by it3)
    for (it1=S.begin(); it1!=S.end(); ++it1) {
        int i1 = (int)distance(S.begin(),it1);
        
        // Compute the neighbor of first row, with bigger index
        //std::cout << "New neighbors of " << distance(S.begin(),it1) << ": ";
        for (it2=it1+1; it2!=S.end(); ++it2) {
            int i2 = (int)distance(S.begin(),it2);
            
            // Compute the union of the two rows
            boost::dynamic_bitset<> u = (*it1) | (*it2);
            // Go through all the other rows and check whether union contains some other row.
            // If such a row is found, that means that the two rows are not adjacent.
            notadjacent = false;
            for (it3=S.begin(); it3!=S.end() && !notadjacent; ++it3) {
                if (it3 != it1 && it3 != it2 && (*it3).is_subset_of(u))
                    notadjacent = true;
            }
            // If an adjacency is detected, add the neighbor list of both rows
            if (!notadjacent) {
                //std::cout << i2 << ' ';
                neighbors[i1].push_back(i2);
                neighbors[i2].push_back(i1);
            }
        }
        //std::cout << std::endl;
    }
    //std::cout << std::endl;
    
    // Computation and verification of submatrices
    //std::cout << "Computation of submatrices..." << std::endl;
    
    accept = true;
    
    for (std::vector<std::vector<int>>::iterator it=neighbors.begin(); it!=neighbors.end() && accept; it++)
        accept = ((*it).size() >= D);
    
    if (!accept)
        return false;
    
    bool found = true;
    
    // Once the dual graph is known, go through all the rows and build the corresponding
    // submatrix for each of them. If the input is a slack matrix, this will compute the
    // slack matrix of the corresponding facet
    for (it1=S.begin(); it1!=S.end() && found; ++it1) {
        
        std::vector<boost::dynamic_bitset<>> SS; // submatrix of S corresponding to the row
        int i1 = (int)distance(S.begin(),it1);
        std::vector<size_type> vert_idx; // column indices in S of zeroes
        size_type num_zeroes = 0; // number of columns
        
        // Compute submatrix for current row
        //std::cout << "Submatrix for row #" << i1;
        
        // Count the zeroes in the row and record their positions
        for (size_type i = 0; i < (*it1).size(); ++i) {
            // If a zero is found, record its position
            if (!(*it1).test(i)) {
                vert_idx.push_back(i);
                ++num_zeroes;
            }
        }
        
        // std::cout << " (" << (*it1) << " - " << "#zeroes = " << num_zeroes << ")" << std::endl;
        
        // Go through all neighbors of the row
        for (std::vector<int>::iterator it4=neighbors[i1].begin(); it4!=neighbors[i1].end(); ++it4) {
            boost::dynamic_bitset<> new_row(num_zeroes);
            // Go through all vertices of the row, and collect the bits from the slack matrix
            // to form the slack matrix for the row (which is a submatrix of S)
            for (size_type i = 0; i < num_zeroes; ++i) {
                size_type j = vert_idx[i];
                new_row[i] = (S[*it4])[j];
            }
            // Add new_row to submatrix
            SS.push_back(new_row);
            //std::cout << new_row << std::endl;
        }
        
        // Turn the submatrix into a nonincidence graph, and canonicize it
        std::vector<setword> canonical = canonicize(SS);
        
        found = false;
        
        // Obtain range with atoms_cg corresponding to nonincidence graphs with same hash
        std::pair <std::multimap<int,std::vector<setword>>::iterator, std::multimap<int,std::vector<setword>>::iterator> r;
        
        r = atoms_cg.equal_range(hash_matrix(SS,D-1));
        
        // Browse through all nonincidence graphs that have the same hash to see if one of them
        // is isomorphic to the current nonincidence graph
        for (std::multimap<int,std::vector<setword>>::iterator it5=r.first; it5!=r.second && !found; ++it5)
            found = ((*it5).second == canonical);
        
        //if (found) {
        //    std::cout << "OK: nonincidence graph isomorphic to that of an atom" << std::endl;
        //}
        //else {
        //    std::cout << "ERROR: nonincidence graph not isomorphic to that of an atom" << std::endl;
        //}
        //
        //std::cout << "-" << std::endl;
    }
    
    return found;
}

// Check equality of ublas::vector
bool is_equal(ublas::vector<int> A,ublas::vector<int> B) {
    assert(A.size() == B.size());
    bool flag = true;
    for (size_type it = 0; it != A.size() && flag; it++)
        flag = (A[it] == B[it]);
    return flag;
}

// test if the 2-level polytope having slack-matrix S_new is a suspension
bool is_susp(std::vector<boost::dynamic_bitset<>> S_new) {
    
    // For all rows i of the slack matrix M
    // Partition the columns into F_0 = {j : M(i,j) = 0} and F_1 = {j : M(i,j) = 1}
    // For all translation vectors t such that the vertices of F_1 - t are a subset of those of F_0
    // Check that F_1 - t is a face by testing whether the intersection of all facets of P that contain it is exactly F_1 - t
    
    bool flag = false;
    
    std::vector<ublas::vector<int>> S_cols; // slack embedding of the polytope having S_new as slack matrix
    
    for (size_type col_idx = 0; col_idx != S_new[0].size(); col_idx++) {
        ublas::vector<int> column(S_new.size());
        for (int it0 = 0; it0 != S_new.size(); ++it0)
            column(it0) = S_new[it0][col_idx];
        S_cols.push_back(column);
    }
    
    for (std::vector<boost::dynamic_bitset<>>::iterator it0 = S_new.begin(); it0 != S_new.end(); ++it0) {
        std::vector<size_type> zeros_idx; // column indices in S of zeroes
        std::vector<size_type> ones_idx; // column indices in S of ones
        
        // Count zeroes and ones in the row and record their positions
        for (size_type i = 0; i < (*it0).size(); ++i) {
            if ((*it0).test(i))
                ones_idx.push_back(i);
            else
                zeros_idx.push_back(i);
        }
        
        for (std::vector<size_type>::iterator it1 = ones_idx.begin();it1 != ones_idx.end(); it1++) {
            for (std::vector<size_type>::iterator it2 = zeros_idx.begin();it2 != zeros_idx.end(); it2++) {
                ublas::vector<int> translation_vect = S_cols[*it2] - S_cols[*it1];
                
                bool is_contained = true;
                
                std::vector<size_type> idx_translated_F1;
                
                for (std::vector<size_type>::iterator it3 = ones_idx.begin();it3 != ones_idx.end() && is_contained; it3++) {
                    ublas::vector<int> translated_F1_point = S_cols[*it3] + translation_vect;
                    
                    bool is_found = false;
                    for (int h = 0; h != zeros_idx.size() && !is_found; h++) {
                        is_found = (is_equal(translated_F1_point,S_cols[zeros_idx[h]]));
                        if (is_found)
                            idx_translated_F1.push_back(zeros_idx[h]);
                        
                    }
                    
                    is_contained &= is_found;
                }
                
                
                if (is_contained) {
                    
                    boost::dynamic_bitset<> char_F1(S_new[0].size());
                    char_F1.set();
                    
                    for (int h = 0; h != idx_translated_F1.size(); h++)
                        char_F1.flip(idx_translated_F1[h]);
                    
                    boost::dynamic_bitset<> intersect_rows_containing_F1(S_new[0].size());
                    
                    for (std::vector<boost::dynamic_bitset<>>::iterator it0 = S_new.begin(); it0 != S_new.end(); ++it0) {
                        if ((*it0).is_subset_of(char_F1))
                            intersect_rows_containing_F1 |= *it0;
                    }
                    
                    flag = (intersect_rows_containing_F1 == char_F1);
                }
            }
        }
    }
    return flag;
}

// Test if the polar of the 2-level polytope having slack-matrix S_new is a still 2_level
int is_polar(std::vector<boost::dynamic_bitset<>> S_new,int &D) {
    std::vector<boost::dynamic_bitset<>> S_transpose;
    for (size_type col_idx = 0; col_idx != S_new[0].size(); col_idx++) {
        boost::dynamic_bitset<> column;
        for (size_type row_idx = 0; row_idx != S_new.size(); row_idx++)
            column.push_back(S_new[row_idx][col_idx]);
        S_transpose.push_back(column);
    }
    
    int hash_S = hash_matrix(S_new,D);
    int hash_S_T = hash_matrix(S_transpose,D);
    std::vector<setword> canonical_S_T = canonicize(S_transpose);
    
    
    int amount_polar = 0;
    
    if (hash_S == hash_S_T) {
        if (canonicize(S_new) == canonical_S_T)
            amount_polar = 1; // self-polar
    }
    else {
        
        bool is_isomorphic = false;
        std::pair <std::multimap<int,std::vector<setword>>::iterator, std::multimap<int,std::vector<setword>>::iterator> r;
        
        r = LD.equal_range(hash_S_T);
        
        // Browse through all nonincidence graphs that have the same hash to see if one of them
        // is isomorphic to the current nonincidence graph
        for (std::multimap<int,std::vector<setword>>::iterator it1=r.first; it1!=r.second && !is_isomorphic; ++it1)
            is_isomorphic = ((*it1).second == canonical_S_T);
        
        if (is_isomorphic)
            amount_polar  = 2;
    }
    
    return amount_polar;
}

// check if the slack matrix S is already listed in LD; if not, add it to LD
void to_list(std::vector<boost::dynamic_bitset<>> S_new, int &D, const int &verbose){
        
    int hash_S = hash_matrix(S_new,D);
    std::vector<setword> canonical_S = canonicize(S_new);
    
    bool is_isomorphic = false;
    
    // Obtain range with atoms_cg corresponding to nonincidence graphs with same hash
    std::pair <std::multimap<int,std::vector<setword>>::iterator, std::multimap<int,std::vector<setword>>::iterator> r;
    
    r = LD.equal_range(hash_S);
    
    // Browse through all nonincidence graphs that have the same hash to see if one of them
    // is isomorphic to the current nonincidence graph
    for (std::multimap<int,std::vector<setword>>::iterator it1=r.first; it1!=r.second && !is_isomorphic; ++it1)
        is_isomorphic = ((*it1).second == canonical_S);
    
    if (!is_isomorphic) {
        LD.insert(std::pair<int,std::vector<setword>>(hash_S,canonical_S));
        
        for (std::vector<boost::dynamic_bitset<>>::iterator it2 = S_new.begin(); it2 != S_new.end(); ++it2)
            myfile << *it2 << std::endl;

        myfile << '-' << std::endl;
    
        if (verbose != 0) {
            // check is there exists a simplicial facet
            bool has_simplicial = false;
            for (std::vector<boost::dynamic_bitset<>>::iterator row = S_new.begin(); row != S_new.end() && !has_simplicial; row++)
                has_simplicial = ((~*row).count() == D);
            if (has_simplicial)
                simplicial_facet++;
            
            // check if there exists a simple vertex
            bool STAB = false;
            for (size_type col_idx = 0; col_idx != S_new[0].size() && !STAB; col_idx++) {
                boost::dynamic_bitset<> column;
                for (int it0 = 0; it0 != S_new.size(); ++it0)
                    column.push_back(S_new[it0][col_idx]);
                STAB = ((~column).count() == D);
            }
            
            if (STAB)
                stab++;

            // check if the polytope is centrally symmetric
            bool CS = true;
            for (std::vector<boost::dynamic_bitset<>>::iterator row = S_new.begin(); row != S_new.end() && CS; row++)
                CS = ((*row).count() == S_new[0].size()/2);
            if (CS)
                cs++;

            // tests if the polytope is a suspension
            if (is_susp(S_new))
                n_suspensions++;
            
            // tests if the polytope has a polar that is a polytope
            n_polar+=is_polar(S_new,D);
        }
    }
}


int main(int argc, const char *argv[]) {
    if (argc != 3) {
        std::cout << "\nERROR: Wrong input, not enough arguments" << std::endl;
        std::cout << "Please insert two arguments: dimension D, verbose_flag." << std::endl;
        std::cout << "- D is an integer between 3 and 7." << std::endl;
        std::cout << "- verbose_flag is an integer between 0 and 3." << std::endl;
        std::cout << "verbose_flag = 0 : minimal output (the benchmark times correspond to this flag)" << std::endl;
        std::cout << "verbose_flag = 1 : output includes all closed sets + tests for classes of 2L polytopes" << std::endl;
        std::cout << "verbose_flag = 2 : output includes all closed sets + H-embedding of ground set + tests for classes of 2L polytopes" << std::endl;
        std::cout << "verbose_flag = 3 : output includes all closed sets, ground sets, slabs_sat_by_point and points_sat_by_slab + tests for classes of 2L polytopes" << std::endl << std::endl;

        std::cout << "The input should be of the form ./2L_enum D verbose_flag, e.g.:" << std::endl;
        std::cout << "./2L_enum 3 3" << std::endl;
        std::cout << "which corresponds to D = 3 and verbose_flag = 3." << std::endl << std::endl;
        return 1;
    }

    int D = atoi(argv[1]);
    int verbose = atoi(argv[2]);

    // Warnings and bound on D and verbose
    if ((D != 3) && (D != 4) && (D != 5) && (D != 6) && (D != 7)) {
        std::cout << "\nERROR: Input (dimension D) out of bounds." << std::endl;
        std::cout << "Please insert an integer value between 3 and 7 as dimension D." << std::endl;
        std::cout << "The input should be of the form: ./2L_enum D verbose_flag." << std::endl << std::endl;
        return 1;
    }

    if ((verbose != 0) && (verbose != 1) && (verbose != 2) && (verbose != 3)) {
        std::cout << "\nERROR: Input (verbose_flag) out of bounds." << std::endl;
        std::cout << "Please insert an integer value between 0 and 3 as verbose_flag." << std::endl;
        std::cout << "verbose_flag = 0 : minimal output (the benchmark times correspond to this flag)" << std::endl;
        std::cout << "verbose_flag = 1 : output includes all closed sets + tests for classes of 2L polytopes" << std::endl;
        std::cout << "verbose_flag = 2 : output includes all closed sets + H-embedding of ground set + tests for classes of 2L polytopes" << std::endl;
        std::cout << "verbose_flag = 3 : output includes all closed sets, ground sets, slabs_sat_by_point and points_sat_by_slab + tests for classes of 2L polytopes" << std::endl << std::endl;

        std::cout << "The input should be of the form: ./2L_enum D verbose_flag." << std::endl << std::endl;
        return 1;
    }

    if (D == 7) {
        std::cout << "\nWARNING: the computation time for D = 7 is expected to take up to ~61 hours of computation time (on AMD Opteron(TM) 6134 2.3 GHz)." << std::endl;
        std::cout << "Please insert 'y' to confirm: ";
        char ch = getchar();
        if((ch != 'y') && (ch != 'Y')) {
            ungetc(ch, stdin);
            return 1;
        }
    }



    std::vector<std::vector<boost::dynamic_bitset<>>> atoms; // Vector with all the (D-1)-dim 2L-polytopes
    std::vector<boost::dynamic_bitset<>> S;             // A slack matrix of a (D-1)-dim polytope
    boost::dynamic_bitset<> row;                   // Some row a slack matrix
    std::string line;
    
    
    // counters for number of 2-level polytopes
    int total_2level = 0;
    //int index = 0;


    // counters for timers
    std::chrono::time_point<std::chrono::system_clock> tot_start, tot_end;
    std::chrono::time_point<std::chrono::system_clock> start_per_base, end_per_base;

    std::chrono::time_point<std::chrono::system_clock> start_next_closure, end_next_closure;
    std::chrono::time_point<std::chrono::system_clock> start_slack_matrix, end_slack_matrix;
    std::chrono::time_point<std::chrono::system_clock> start_skip_test, end_skip_test;
    std::chrono::time_point<std::chrono::system_clock> start_2level_test, end_2level_test;
    double tot_next_closure = 0;
    double tot_slack_matrix = 0;
    double tot_2level_test = 0;

    
    tot_start = std::chrono::system_clock::now();
    
    std::cout << "\nWarning: in dynamic_bitsets, lowest index bits are written to the right!" << std::endl;
    std::cout << "Reading all " << (D-1) << "-dimensional 2-level polytopes... ";
    
    if (verbose != 0)
        std::cout << std::endl;
    
    // Open the file that contains (D-1)-dim 2L-polytopes
    std::ifstream inputfile("./" + std::to_string(D-1) + "d.txt");
    if (inputfile.is_open()) {
        while (getline(inputfile,line)) {
            // If we read some data
            if (line.compare(std::string("-")) != 0) {
                // Remove spaces
                // string::iterator end_pos = std::remove(line.begin(),line.end(),' ');
                // line.erase(end_pos, line.end());
                // Turn the std::string into a set
                row = (boost::dynamic_bitset<>)line;
                if (verbose != 0)
                    std::cout << row << " (" << row.size() << ")" << std::endl;
                S.push_back(row);
            }
            // Otherwise this is the end
            else {
                // Add the slack matrix to the list
                atoms.push_back(S);
                // Don't forget to clear S
                S.clear();
                if (verbose != 0)
                    std::cout << "-" << std::endl;
            }
        }
        inputfile.close();
        std::cout << "OK" << std::endl;
    }
    else {
        std::cout << "Unable to open file." << std::endl;
        return 1;
    }
    
    std::cout << "Number of polytopes read = " << atoms.size() << std::endl;
    std::cout << "Computing canonical forms for all nonincidence graphs... ";
    
    // Loop trough all the atoms, compute canonical form and store it
    for (std::vector<std::vector<boost::dynamic_bitset<>>>::iterator it1=atoms.begin(); it1!=atoms.end(); ++it1)
        atoms_cg.insert(std::pair<int,std::vector<setword>>(hash_matrix(*it1,D-1),canonicize(*it1)));
    
    std::cout << "OK" << std::endl;
    
    std::cout << "Processing bases..." << std::endl;
    
    myfile.open(std::to_string(D) + "d.txt");
    
    int N_closed_sets = 0;
    
    // Main loop: loop through all the bases
    for (std::vector<std::vector<boost::dynamic_bitset<>>>::iterator it1=atoms.begin(); it1!=atoms.end(); ++it1) {
        start_per_base = std::chrono::system_clock::now();

        std::vector<boost::dynamic_bitset<>> S = *it1;
        
        int n_base = (int)distance(atoms.begin(),it1) +1;
        std::cout << "\nBase #" << n_base << ':' << std::endl;
        
        // Check if there is a simplicial core on the "top left" of the matrix
        std::cout << "Simplicial core? ";
        if (checksimplicialcore(S,D-1))
            std::cout << "OK" << std::endl;
        
        // Extract (extended) embedding transformation matrix M_d(0)
        ublas::matrix<int> M = extractM(S,D);
        if ((verbose==1) || (verbose==2) || (verbose == 3))
            std::cout << "M_d(0) = " << M << std::endl;
        
        // Compute inverse of matrix M_d(0)
        ublas::matrix<int> Minv;
        invertM(M,Minv);
        if ((verbose==1) || (verbose==2) || (verbose == 3))
            std::cout << "M_d(0)^{-1} = " << Minv << std::endl;
        
        std::cout << "Constructing H-embedding of facets of the base... ";
        // computing the facets of the base using the slack matrix S
        std::vector<boost::dynamic_bitset<>> facets_base(S.size());
        for (std::vector<boost::dynamic_bitset<>>::const_iterator it2 = S.begin(); it2 != S.end(); it2++) {
            boost::dynamic_bitset<> E(D);
            if ((*it2).test(D-1)) {
                for (size_type i = 0; i < D-1; i++)
                    E[i+1] = !(*it2)[i];
            }
            else {
                for (size_type i = 0; i < D-1; i++)
                    E[i+1] = (*it2)[i];
            }
            
            bool found = false;
            for (std::vector<boost::dynamic_bitset<>>::iterator it3 = facets_base.begin(); it3 != facets_base.end() && !found; it3++)
                found = ((*it3) == E);
            if (!found) {
                facets_base.push_back(E);
                if (verbose == 3)
                    std::cout << E << ' ';
            }
        }
        std::cout << "OK" << std::endl;
        
        // Create the set Vert(P_0) (in V-embedding)
        std::cout << "Building V-embedding of base... ";
        
        std::vector<ublas::vector<int>> base_V;
        
        // Loop through all vertices
        for (int j=0; j<S[0].size(); j++) {
            ublas::vector<int> point(D);
            
            // Create a point whose first coordinate is 0, and the others are the D-1 first bits of the jth column of the slack matrix S
            point[0] = 0;
            for (int i=0; i<D-1; i++)
                point[i+1] = (S[i])[j];
            
            // Add point to the V-embedding of the ground set
            base_V.push_back(point);
            
            // Print point - for debugging
            if (verbose==3)
                std::cout << point << ' ';
        }
        std::cout << "OK" << std::endl;
        
        
        std::cout << "Building V-embedding of the ground set... ";
        
        std::vector<ublas::vector<int>> ground_set_V;
        
        ublas::vector<int> count(D);
        // initialize count to 0
        for (int i=0; i<D; i++)
            count(i) = 0;
        
        bool carry;
        
        
        // REDUCED GROUND SET
        // Create the set of points in {1} x {-1,0,1}^{D-1} lex greater than e_1
        // Its size is 1 + (3^{D-1} - 1)/2
        
        ublas::vector<int> point(D);
        point(0) = 1;
        
        for (int i=1; i<D; i++)
            point(i) = 0;
        
        ground_set_V.push_back(point);
        if (verbose == 3)
            std::cout << point << ' ';
        
        for (int i=D-1; i>0 ; i--) {
            point(i) = 1;
            
            count.clear();
            
            while (count(D-1) == 0) {
                int j;
                // Extract a vector in {-1,0,1}^{D-i-1} to fill the vector
                for (j=i+1; j<D; j++)
                    point(j) = count(j-1) - 1;
                
                ground_set_V.push_back(point);
                if (verbose == 2)
                    std::cout << point << ' ';
                
                // Increase counter, by performing mod-3 computation
                j = i;
                do {
                    carry = (count(j) == 2);
                    count(j) = (count(j) + 1) % 3;
                    j++;
                } while (carry && j < D);
            }
        }
        std::cout << "OK" << std::endl;
        
        // Create Vert(P_0) (H-embedding this time), the set of fixed points
        std::cout << "Building H-embedding of base... ";
        
        std::vector<ublas::vector<int>> base_H;
        
        for (std::vector<ublas::vector<int>>::iterator it2=base_V.begin(); it2!=base_V.end(); it2++) {
            ublas::vector<int> point(D);
            point = prod(Minv,*it2);
            base_H.push_back(point);
            if (verbose == 3)
                std::cout << point << ' ';
        }
        std::cout << "OK" << std::endl;
        
        // Create ground set
        std::cout << "Building H-embedding of the reduced ground set... ";
        
        std::vector<ublas::vector<int>> ground_set_H;
        
        for (std::vector<ublas::vector<int>>::iterator it2=ground_set_V.begin(); it2!=ground_set_V.end(); it2++) {
            ublas::vector<int> point(D);
            point = prod(Minv,*it2);
            
            
            // Facet reduction of the ground set:
            // we can throw away all the points x of the ground set where we do not have x(E) in {-1,0,1}
            // for x(E) >= 0, x(E) <= 1 facet of the base, E subset of {2,...,d}
            bool accept = true;
            for (std::vector<boost::dynamic_bitset<>>::iterator it4 = facets_base.begin(); it4 != facets_base.end() && accept; it4++) {
                int xE = 0;
                for (size_type k=(*it4).find_first(); k!= npos; k=(*it4).find_next(k))
                    xE += point[k];
                accept = ((xE == -1) || (xE == 0) || (xE == 1));
            }
            
            if (accept) {
                ground_set_H.push_back(point);
                if (verbose == 3)
                    std::cout << point << ' ';
            }
        }
        std::cout << "OK" << std::endl;
        std::cout << "Size of the ground set = " << (int)ground_set_V.size() << std::endl;
        std::cout << "Size of the reduced ground set = " << (int)ground_set_H.size() << std::endl;
    
        // Compute the inequalities x(E) <= 1, x(E) >= 0. What about embedding everything in dim D+1?
        count.resize(D+1);
        count.clear();
        
        // Skip the all-0 vector
        count(0) = 1;
        
        std::cout << "Building slabs... ";
        std::vector<ublas::vector<int>> slabs;
        
        while (count(D) == 0) {
            int i;
            ublas::vector<int> normal_vector(D);
            
            for (i=0; i<D; i++) {
                normal_vector(i) = count(i);
            }
            
            bool check = true;
            
            for (std::vector<ublas::vector<int>>::iterator it2=base_H.begin(); it2!=base_H.end() && check; it2++) {
                int s = ublas::inner_prod(normal_vector,*it2);
                check = ((s == 0) || (s == 1));
            }
            
            // Add normal vector of slab to the list if it contains all points of the base
            if (check) {
                slabs.push_back(normal_vector);
                
                // Print normal vector - for debugging
                if (verbose == 3)
                    std::cout << normal_vector << ' ';
            }
            
            // Increase counter, by performing mod-2 computations
            i = 0;
            do {
                carry = (count(i) == 1);
                count(i) = (count(i) + 1) % 2;
                i++;
            } while (carry && i <= D);
        }
        std::cout << "OK" << std::endl;
        
        // Reinitialize global data structures
        slabs_sat_by_point.clear();
        points_sat_by_slab.clear();
        
        std::cout << "Building incidences between points and slabs... ";
        
        // Check points versus slabs (for each point, list the slabs containing it)
        for (std::vector<ublas::vector<int>>::iterator it2=ground_set_H.begin(); it2!=ground_set_H.end(); it2++) {
            boost::dynamic_bitset<> slabs_sat(slabs.size());
            for (std::vector<ublas::vector<int>>::iterator it3=slabs.begin(); it3!=slabs.end(); it3++) {
                int s = ublas::inner_prod(*it2,*it3);
                bool sat = ((s == 0) || (s == 1));
                slabs_sat[distance(slabs.begin(),it3)] = sat;
            }
            slabs_sat_by_point.push_back(slabs_sat);
            if (verbose == 3)
                std::cout << slabs_sat << ' ';
        }
        
        if (verbose == 3)
            std::cout << "| ";
        
        // Once more, exchanging roles (for each slab, list the points satisfying it)
        for (std::vector<ublas::vector<int>>::iterator it2=slabs.begin(); it2!=slabs.end(); it2++) {
            boost::dynamic_bitset<> points_sat(ground_set_H.size());
            for (std::vector<ublas::vector<int>>::iterator it3=ground_set_H.begin(); it3!=ground_set_H.end(); it3++) {
                int s = ublas::inner_prod(*it2,*it3);
                bool sat = ((s == 0) || (s == 1));
                points_sat[distance(ground_set_H.begin(),it3)] = sat;
            }
            points_sat_by_slab.push_back(points_sat);
            if (verbose == 3)
                std::cout << points_sat << ' ';
        }
        
        std::cout << "OK" << std::endl;
        
        
        std::cout << "Constructing the incompatibility matrix... ";

        std::vector<boost::dynamic_bitset<>> incompatibility_adjM;
        
        for (size_type i = 0; i != ground_set_H.size(); i++) {
            boost::dynamic_bitset<> empty_row;
            incompatibility_adjM.push_back(empty_row);
            for (size_type j = 0; j != i; j++) {
                bool found = false;
                for (std::vector<boost::dynamic_bitset<>>::iterator it4 = facets_base.begin(); it4 !=facets_base.end() && !found;it4++) {
                    int s_i = 0;
                    int s_j = 0;
                    for (size_type k=(*it4).find_first(); k!= npos; k=(*it4).find_next(k)) {
                        s_i += ground_set_H[i][k];
                        s_j += ground_set_H[j][k];
                    }
                    found = (s_i*s_j == -1);
                }
                incompatibility_adjM[i].push_back(found);
            }
        }
        
        /*
        std::cout << std::endl;
        for (std::vector<boost::dynamic_bitset<>>::iterator p = incompatibility_adjM.begin(); p != incompatibility_adjM.end(); p++)
            std::cout << *p << std::endl;
        */
        std::cout << "OK" << std::endl;
    
        std::cout << "Lauching Ganter's next-closure algorithm and checking 2-levelness... ";
        
        if (verbose != 0)
            std::cout << std::endl;
        else if (verbose == 0)
            std::cout << std::endl;
        int c = 0; // counter to print dots in case verbose = 0, every dot corresponds to 100 closed sets
        
        //std::vector<boost::dynamic_bitset<>> closed_sets;
        int N_closed_sets_current_base = 0;
        int N_2level_current_base = 0;

        
        int n = (int)ground_set_H.size();
        
        if (verbose != 0)
            std::cout << std::string(n, ' ') << "  | 2-level | next_cl      | slack-matrix | 2-lev_time" << std::endl;
        
        boost::dynamic_bitset<> A(n);
        
        while (!A.all()) {
            // compute the complement of A
            boost::dynamic_bitset<> I, CI;
            size_type i;
            
            
            while (!A.all()) {
                
                //inizialization of the set of all slabs belonging to elements of A
                boost::dynamic_bitset<> B(points_sat_by_slab.size());
                
                
                i = 0;
                if (verbose != 0)
                    start_next_closure = std::chrono::system_clock::now();
                do {
                    while(A[i])
                        ++i;
                    I = inc(A,(int)i);
                    //CI = discreteconvexhull_cl(I,B);
                    boost::dynamic_bitset<> dchcl = discreteconvexhull_cl(I,B);
                    CI = incompatibility_cl(incompatibility_adjM,dchcl);
                    ++i;
                } while (!sqsubseteq(I,CI));
                if (verbose != 0)
                    end_next_closure = std::chrono::system_clock::now();
                
                A = CI;
                
                if (verbose != 0)
                    std::cout << A << ' ';
                else if (verbose == 0) {
                    c++;
                    if (c % 100 == 0)
                        std::cout << '.';
                    if ((c % 10000 == 0) && (c >= 10000))
                       std::cout << c/1000 << 'k' << std::endl;
                }
                //closed_sets.push_back(A);
                N_closed_sets_current_base++;
                
                if (verbose != 0)
                    start_slack_matrix = std::chrono::system_clock::now();
                // construct the slack matrix S with embedding transformation matrix in top right position
                std::vector<boost::dynamic_bitset<>> S_new;
                construct_slack_matrix(base_H,ground_set_H,A,B,slabs,S,S_new,D);
                
                if (verbose != 0)
                    end_slack_matrix = std::chrono::system_clock::now();

                
                bool accept = true;
                
                start_skip_test = std::chrono::system_clock::now();
                for (std::vector<boost::dynamic_bitset<>>::iterator it_row = S_new.begin(); it_row != S_new.end() && accept; it_row++)
                    accept = ((~(*it_row)).count() <= S[0].size());
                end_skip_test = std::chrono::system_clock::now();
                
                bool istwolevel = 0;
                
                if (accept) {
                    if (verbose != 0)
                        start_2level_test = std::chrono::system_clock::now();
                    
                    istwolevel = istwolevelpolytope(S_new,D);

                    if (verbose != 0)
                        end_2level_test = std::chrono::system_clock::now();
                }
                
            
                
                
                
                if (verbose != 0) {
                    std::chrono::duration<double> time_next_closure = end_next_closure-start_next_closure;
                    tot_next_closure += time_next_closure.count();
                    
                    std::chrono::duration<double> time_slack_matrix = end_slack_matrix-start_slack_matrix;
                    tot_slack_matrix += time_slack_matrix.count();

                    std::cout << " | ";
                    if (accept)
                        std::cout << istwolevel;
                    else
                        std::cout << '-';
                    std::cout << "       | ";
                    std::cout << std::left << std::setw(11) << time_next_closure.count() << "s | ";
                    std::cout << std::left << std::setw(11) << time_slack_matrix.count() << "s | ";
                    if (accept) {
                        std::chrono::duration<double> time_2level_test = end_2level_test-start_2level_test;
                        tot_2level_test += time_2level_test.count();
                        std::cout << std::left << std::setw(11) << time_2level_test.count() << "s ";
                    }
                    else {
                        std::chrono::duration<double> time_skip_test = end_skip_test-start_skip_test;
                        tot_2level_test += time_skip_test.count();
                        //std::cout << left << setw(8) << time_skip_test.count() << "s";
                    }

                    std::cout << std::endl;
                }
            
                if (istwolevel) {
                    to_list(S_new,D,verbose);                    
                    N_2level_current_base++;
                    total_2level++;
                }
                S_new.clear();
            }
            incompatibility_adjM.clear();
            slabs_sat_by_point.clear();
            points_sat_by_slab.clear();
        }
        if (verbose != 0) {
            std::cout << std::string(n+11, ' ');
            std::cout << " | ";
            std::cout << std::left << std::setw(11) << tot_next_closure << "s | ";
            std::cout << std::left << std::setw(11) << tot_slack_matrix << "s | ";
            std::cout << std::left << std::setw(11) << tot_2level_test << "s";
            std::cout << std::endl;
            
            tot_next_closure = 0;
            tot_slack_matrix = 0;
            tot_2level_test = 0;
        }
        
        if (verbose == 0)
            std::cout << ' ';
        std::cout << "OK" << std::endl;
        
        std::cout << "#closed sets found for base #" << n_base << " = " << N_closed_sets_current_base << std::endl;
        std::cout << "#2L polytopes for base #" << n_base << " = " << N_2level_current_base << std::endl;
        //N_2level_current_base=0;

        N_closed_sets += N_closed_sets_current_base;
        end_per_base = std::chrono::system_clock::now();
        std::chrono::duration<double> time_per_base = end_per_base-start_per_base;
        std::cout << "Elapsed time for base #"<< n_base << " = " << time_per_base.count() << "s" << std::endl;
        
    }
    
    tot_end = std::chrono::system_clock::now();
    std::chrono::duration<double> tot_elapsed_seconds = tot_end-tot_start;
    
    std::cout << std::endl;
    std::cout << "Total #closed sets found = "<< N_closed_sets << std::endl;
    std::cout << "#2-level polytopes tested = "<< total_2level << std::endl;
    std::cout << "#non-isomorphic 2-level polytopes found = "<< LD.size() << std::endl;
    std::cout << "Total elapsed time = "<< tot_elapsed_seconds.count() << "s" << std::endl;
    
    if (verbose != 0) {
        std::cout << std::endl;
        std::cout << " D | Delta-f |      CS |    STAB |   polar |    susp |      2L" << std::endl;
        std::cout << ' ' << D << " | ";
        std::cout << std::right << std::setw(7) << simplicial_facet << " | ";
        std::cout << std::right << std::setw(7) << cs << " | ";
        std::cout << std::right << std::setw(7) << stab << " | ";
        std::cout << std::right << std::setw(7) << n_polar << " | ";
        std::cout << std::right << std::setw(7) << n_suspensions << " | ";
        std::cout << std::right << std::setw(7) << LD.size() << std::endl;
    }

    std::cout << std::endl;
    
    
    myfile.close();
    

    
}

