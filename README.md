# Positroids
Python code for computation with positroids, totally nonnegative matrices, Wilson loop diagrams, and H-primes.

Uses Python 2.7 and Sympy 1+. Not tested on Python 3.
Note: Sage uses Python 2.6 and is not compatible. To get the data into Sage (e.g. for homology computations), you need to export the data into a text file and then import that into Sage. See 2_6_data.py for examples.

For the code used to obtain the results on 2,6 Wilson loop diagrams (arXiv:1803.00958), see 2_6_data.py.


### To do
Right now this is just a collection of functions, written as and when I needed that particular thing for a computation. As I learn about OOP, I'm planning to rewrite this so that each positroid cell is an instance of a Positroid class and holds all the attributes of that cell.


### To use
Import positroids.py in your Python shell for interactive computation with the functions listed below.

```python

from positroids import *
from sympy import *



########
# Syntax
########

# WLD: (n,[[vertices of prop 1],[vertices of prop 2],...])
# Le diagram: Matrix([[row 1],[row 2]....])
# GN: [[I1],[I2],...]
# Permutations: written in truncated 2 line notation, i.e. [3,2,4,1] is the permutation that maps 1 -> 3, 2 -> 2, 3 -> 4, 4 -> 1.

### Example:
# A 6 vertex WLD with propagators (1,3) and (3,5):
W = (6, [[1, 2, 3, 4], [3, 4, 5, 6]])

# Its Le diagram is
Matrix([[0,1,1,1],[1,1,1,x]]) 
# (the x denotes a missing square, i.e. the diagram is not rectangular)

# Its Grassmann necklace is
[[1, 3], [2, 3], [3, 4], [4, 5], [5, 1], [6, 1]]

# Its decorated permutation is
[2,4,5,1,6,3]



########
# A few miscellaneous functions
########


Aq2(k,n) 
#	This counts how many Le diagrams of each dimension there are in Gr(k,n)
#	The coefficient of q^i is how many cells of dimension i there are.

allLe(k,n,dim)
#	This generates all Le diagrams in Gr(k,n)
#	If dim = False (this is the default) then the output is a list
#	If dim = True, output is a dictionary, key i produces list of all i-dim Le diagrams

LetoLatex(M)
#	Turns a Le diagram M into a string containing the tex code for that diagram
# (current formatting is for ytableau package)


allWLD(k,n):
#	Generates all WLD of size (k,n)
#	The syntax for a WLD is (n,[P1,P2,...,Pk]) where each P is a list P=[i,i+1,j,j+1] listing the propagator's support


##########
# Now lots of functions for converting between the various objects:
##########

WLDtoGr(W):
#	turns a WLD into an elt of Gr(k,n)

WLDtoGN(W):
#	turns a WLD into a Grassmann necklace

WLDtoDP(W):
#	turns a WLD into a decorated permutation

WLDtoLe(W):
#	turns a WLD into its Le diagram

GrtoGN(M):
#	turns a TNN element of Gr(k,n) into its Grassmann necklace

GrtoDP(M):
#	turns a TNN element of Gr(k,n) into its decorated permutation

GrtoLe(M):
#	turns a TNN element of Gr(k,n) into its Le diagram

GNtoDP(I):
#	turns a Grassmann necklace into a decorated permutation




########
# Some functions for computing reduced words and boundaries
########

LetoRP(L,n)
#	L is the Le diagram, n is whatever dimension you're in, i.e. the n in Gr(k,n)
#	This turns a Le diagram into a reduced permutation, needed for boundaries
#	Output is a pair of permutations (V,W) corresponding to the pair of permutations defining the cell


isbound(V,W,V1,W1,k):
#	Used to check whether one cell is on the boundary of another
#	Specify permutations attached to two Le diagrams, i.e. LetoRP(L1,n) = (V,W) and LetoRP(L2,n) = (V1,W1)
#	Also specify k, i.e the k in whichever Gr(k,n) you're working in
#	This function returns True iff L2 lies on the boundary of L1, i.e. (V1,W1) lies on the boundary of (V,W)




###########
# Example of computing boundaries in Gr(k,n)
###########


### first set your values of k,n and define the WLD

k,n = 2,6
W = ## define your WLD here ## 
V1 = WLDtoLe(W)



V1RP = LetoRP(V1,n)

AL = allLe(k,n,False)
ALRP = [LetoRP(M,n) for M in AL]

## lists are ordered, so ALRP[i] is the reduced permutation of the Le diagram AL[i]

# produce all bounds of V1
boundsV1 = [AL[i] for i in range(len(AL)) if isbound(V1RP[0],V1RP[1],ALRP[i][0],ALRP[i][1],k)]


# pick out just the d-dim ones
d = 5
[M for M in boundsV1 if flatten(M).count(1) == d]


```
