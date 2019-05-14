from sympy import *
from sympy import pprint
import sys
sys.displayhook = pprint
from sympy.combinatorics.subsets import ksubsets
var('x')
var('q')


# first, a few auxilary functions
Intersection = lambda I,J:sorted([i for i in I if i in J])
SetDifference = lambda A,B:[l for l in A if l not in B]

# sometimes we want to replace a matrix with a list of lists
# don't even ask, I got this off stackoverflow
flatten = lambda M:[j for i in M.tolist() for j in i]

# generates all k-subsets of a set I, usually range(1,n+1)
Subsets = lambda I,k:[list(i) for i in list(ksubsets(I,k))] 




# divmod(n,s) returns (q,r) where n = q*s + r
# in other words, divmod(n,s) is the same as (n//s,n%s)
# in particular, if M has n cols then [divmod(i,n) for i in range(len(M))] is a list of the coords in M 


def qbinom(j):
	if j <=0:
		return 0
	else:
		s = 0
		for i in range(j):
			s+= q**i
		return s

#this counts cells
def Aq2(k,n):
	total=0
	for i in range(k):
		total+= (-1)**i * binomial(n,i)*(q**(k*i)*qbinom(k-i)**i*qbinom(k-i+1)**(n-i) - q**((k+1)*i)*qbinom(k-i-1)**i*qbinom(k-i)**(n-i))
	total = total*(q**(-(k**2)))
	return total.expand()


def isYoung(M):
	m,n=M.rows,M.cols
	x_index=[i for i in range(len(M)) if M[i]==x]
	for i in x_index:
		if flatten(M[i//n,i%n:]).count(x) != len(M[i//n,i%n:]) or flatten(M[i//n:,i%n]).count(x) != len(M[i//n:,i%n]):
			return False
	return True

def isLe(M):
	m,n=M.rows,M.cols
	zero_index=[i for i in range(len(M)) if M[i]==0]
	for i in zero_index:
		if flatten(M[i//n,:i%n]).count(1) > 0 and flatten(M[:i//n,i%n]).count(1) >0:
			return False
	return True


def applyYoungShape(M,I):
	m,n=M.rows,M.cols
	a,b=0,n-1
	M1=M[:,:]
	for n in range(0,max(I)):
		if n in I:
			a+=1
		else:
			M1[a:,b]=Matrix(m-a,1,lambda a,b: x)
			b+=(-1)
	return M1

def verifyAllLe(L):
	m,n=L[0].rows,L[0].cols
	P = Poly(Aq2(m,m+n),q)
	C=P.coeffs()[::-1]
	dim_vect = [0]*(m*n+1)
	for M in L:
		dim_vect[flatten(M).count(1)]+=1
	if dim_vect==C:
		return True
	else:
		return False,C,dim_vect




def allLe(k,n,dim=False):
	L=[ones(k,n-k)]
	SC=Subsets(range(n),k)
	SC.pop(0)
	for I in SC:
		L+=[applyYoungShape(ones(k,n-k),I)]
	L1=[]
	for base in L:
		L2=[base[:,:]]
		for i in range(k*(n-k)):
			for j in range(len(L2)):
				if L2[j][divmod(i,n-k)]==x:
					pass
				else:
					M1=L2[j][:,:]
					M1[divmod(i,n-k)]=0
					if isLe(M1):
						L2+=[M1]
		L1+=L2
	if verifyAllLe(L1)==True and dim==False:
		print True
		return L1
	elif verifyAllLe(L1)==True and dim==True:
		D={i:[] for i in range(k*(n-k)+1)}
		for M in L1:
			D[flatten(M).count(1)]+=[M]
		print True
		return D
	else:
		print False
		return L1


# turns a Le diagram from Python into the corresponding LaTeX code
# now using ytableau instead of genyoungtabtikz

# \ytableaushort{+00,+++,0}

def LetoLatex(M):
	n=M.cols
	latex_str = '\\ytableaushort{'
	for j in range(len(M)):
		if M[j]==x:
			pass
		elif M[j]==1:	
			latex_str+='+'
		else:
			latex_str+='0'
		if j%n==n-1:
			latex_str+=','
	while latex_str[-1]==',':
		latex_str = latex_str[:-1]
	latex_str+='}'
	if latex_str == "\ytableaushort{}":
		latex_str = "$\\emptyset$"
	return latex_str



def GNtoLatex(L):
	latex=''
	for I in L:
		for j in I:
			latex+=str(j)
		latex+=', '
	return latex[:-2]




# syntax for new variables: var('a:c1:4') etc

# format for WLD: W = (n,[[prop1],[prop2],...])
# e.g. W = (6, [[1,2,3,4],[3,4,5,6]])


# compose permutations right to left
# this is not necc the RW coming from the Le diagram
def RW(P):
	P1=P[:]
	T=[]
	for i in range(1,len(P)):
		j=i
		while P1[j] < P1[j-1]:
			temp=P1[j]
			P1[j]=P1[j-1]
			P1[j-1]=temp
			T+=[j]
			j+=-1
			if j==0:
				break
	assert P1 == range(1,len(P)+1)
	return T[::-1]




def length(P):
	return len(RW(P))

# maybe we need to specify k here?
def LetoRW(L):
	k,m=L.rows,L.cols
	W,V=[],[]
	for r in range(k-1,-1,-1):
		for c in range(m-1,-1,-1):
			if L[r,c]!=x:
				W+=[c+k-r]
				if L[r,c]==0:
					V+=[c+k-r]
	return V,W

# precompose P with elem perm s_i 
def Si(i,P):
	temp=P[i]
	P[i]=P[i-1]
	P[i-1]=temp
	return P

def RWtoRP(V,n):
	P=range(1,n+1)
	for k in V:
		P=Si(k,P)
	return P

def LetoRP(L,n):
	V,W = LetoRW(L)
	return RWtoRP(V,n),RWtoRP(W,n)








# just the 1-order for now
# returns true if A <= B
def Gale(A,B):
	assert len(A) == len(B)
	n=len(A)
	A1,B1=sorted(A),sorted(B)
	for i in range(n):
		if A1[i] > B1[i]:
			return False
	return True



# returns True if P <= Q in the strong Bruhat order
def Bruhat(P,Q,k):
	assert len(P)==len(Q)
	n=len(P)
	for i in range(k):
		if Gale(P[:i],Q[:i])==False:
			return False
	for i in range(n-1,k-1,-1):
		if Gale(Q[i:],P[i:])==False:
			return False
	return True

# checks whether parabolic Z with V1 <= V2Z has a chance of existing
def checkZ(V1,V2,k):
	VS = sorted(V1[:k]) + sorted(V1[k:])
	return Bruhat(VS,V2,k)



# find z in the parabolic so that V1 <= V2Z
def findZ(V1,V2,k):
	assert len(V1)==len(V2)
	n = len(V1)
	if Bruhat(V1,V2,k)==True:
		return range(1,n+1)
	Z,VZ=range(n),range(n)
	for i in range(0,k):
		ZC = sorted([j for j in range(n) if j not in Z[:i]])
		r=0
		VZ[i]=V2[ZC[r]]
		while Gale(V1[:i+1],VZ[:i+1])==False:
			r+=1
			VZ[i]=V2[ZC[r]]
		Z[i]=ZC[r]
	for i in range(n,k,-1):
		i
		ZC = sorted([j for j in range(n) if j not in Z[:k]+Z[i:]])
		ZC
		r=len(ZC)-1
		VZ[i-1]=V2[ZC[r]]
		Gale(VZ[i-1:],V1[i-1:])
		while Gale(VZ[i-1:],V1[i-1:])==False:
			r+=-1
			VZ[i-1]=V2[ZC[r]]
			Gale(VZ[i-1:],V1[i-1:])
		Z[i-1]=ZC[r]
	assert sorted(Z)==range(n)
	assert Bruhat(V1,[V2[r] for r in Z],k)
	return [z+1 for z in Z]

# having constructed Z above, we need to check that W2Z <= W1 as well
def verifyZ(W2,W1,Z,k):
	WZ = [W2[z-1] for z in Z]
	return Bruhat(WZ,W1,k)


# returns true if V1,W1 is a bound of V,W
# this happens if and only if there exists a parabolic z with V <= V1z and W1z <= W
def isbound(V,W,V1,W1,k):
	if checkZ(V,V1,k)==False:
		return False
	else:
		Z=findZ(V,V1,k)
		if verifyZ(W1,W,Z,k)==False:
			return False
		else:
			return True







# turns a WLD into an element of the Grassmannian
# just puts an indeterminate p_i into all the nonzero entries
# we'll worry about minus signs later if it becomes important

def WLDtoGr(W):
	M=zeros(len(W[1]),W[0])
	n = 4*len(W[1])
	indets = [var('p'+str(i)) for i in range(n)]
	for i in range(len(W[1])):
		for j in range(len(W[1][i])):
			M[i,W[1][i][j]-1] = indets[4*i+j]
	return M

# grassmannian matrix to grassmann necklace
def GrtoGN(M):
	n=M.cols
	I=[0]*n
	for i in range(n):
		Mc=M[:,i:].row_join(M[:,:i])
		Mrref=Mc.rref(simplify=True)
		I[i]=[(a + i)%n +1 for a in Mrref[1]]
	return I

# grassmann necklace to decorated permutation 
# Curently doesn't record the decoration of fixed points, because WLD don't admit coloops
def GNtoDP(I):
	n=len(I)
	P=[0]*n
	for i in range(len(I)):
		J=SetDifference(I[(i+1)%n],I[i])
		if J==[]:
			P[i]=i+1
		else:
			P[i]=J[0]
	assert sorted(P)==range(1,n+1)
	return P



# here are the missing links, for convenience.

# there aren't any going in the opposite direction yet, e.g. GN to WLD. I guess this is because
# we don't have criteria for a GN to give an admissible WLD?

def WLDtoGN(W):
	M=WLDtoGr(W)
	return GrtoGN(M)

def WLDtoDP(W):
	I=WLDtoGN(W)
	return GNtoDP(I)

def GrtoDP(M):
	I=GrtoGN(M)
	return GNtoDP(I)



# now for some derivation deleting
# this one runs the Deleting Derivations algorithm on a matrix (i.e. an actual matrix, not Gr)

def DDMatrix(M):
	M1=M[:,:]
	for i in range(M1.rows-1,0,-1):
		for j in range(M1.cols-1,0,-1):
			if M1[i,j]==0:
				pass
			else:
				for i1 in range(i-1,-1,-1):
					for j1 in range(j-1,-1,-1):
						M1[i1,j1]=(M1[i1,j1]-(M1[i,j1]*M1[i1,j])/M1[i,j]).expand().factor().cancel()
						# sympy can be really precious about identifying when stuff is zero
						# hence the repeated simplifying. 
	for i in range(len(M1)):
		M1[i]=M1[i].simplify()
	return M1


# DD on an element of Gr(k,n)
def DDGr(M):
	# need to row reduce, map over to a matrix, then DDMatrix
	# sympy has built-in rref: returns a tuple (RREF,[pivots])
	# need to use sipmlify=True to guarantee it identifies zeros correctly
	Mrref=M.rref(simplify=True)
	M1=[flatten(Mrref[0][:,i]) for i in range(Mrref[0].cols) if i not in Mrref[1]] # removes the pivot rows
	# note: this is a list of lists, not a matrix right now
	# each list is a column, so we just need to reorder the outer list
	M1=M1[::-1]
	M2=Matrix(M1).transpose()
	for i in range(M2.rows):
		for j in range(M2.cols):
			M2[i,j]*=(-1)**(M2.rows-i+1) 
	return DDMatrix(M2),Mrref[1]


def GrtoLe(M):
	# input a TNN element of Gr(k,n), returns its Le diagram
	# convention: use 0,1 instead of 0,+ in the Le diagram because Python gets confused otherwise
	# use x to indicate the shape of the diagram, i.e. a box with an x is deleted from the final shape	
	Mnew,SC = DDGr(M)
	for i in range(Mnew.rows):
		for j in range(Mnew.cols):
			if Mnew[i,j]!=0:
				Mnew[i,j]=1
	# we need to distinguish between 0 squares in the diagram and the shape of the diagram
	# the shape is given by the Schubert cell, which is SC above
	x = Symbol('x')
	i,j=0,Mnew.cols-1
	for n in range(0,max(SC)):
		if n in SC:
			i+=1
		else:
			Mnew[i:,j]=Matrix(Mnew.rows-i,1,lambda i,j: x)
			j+=(-1)
	return Mnew




# shortcut if you want to go straight from the WLD
# recall WLD notation: W = (n, [a,b,c,d], [e,f,g,h],...) where letters are the vertices supporting the propagator
def WLDtoLe(W):
	M=WLDtoGr(W)
	return GrtoLe(M)



# now we want to generate all wilson loop diagrams of a given size
# these correspond to partial triangulations of the ngon

# these are all possible propagator lines on an ngon
chords = lambda n:[list(i) for i in list(ksubsets(range(1,n+1),2)) if i[1]-i[0] not in [1,n-1]]

# tests a list of chords for noncrossing
def noncrossing(I):
	for a in range(len(I)):
		for b in range(a+1,len(I)):
			A=min(I[a],I[b])
			B=max(I[a],I[b])
			if A[0] < B[0] < A[1] < B[1]:
				return False
			else:
				pass
	return True



def allWLD(k,n):
	I=chords(n)
	indices = [list(i) for i in list(ksubsets(range(len(I)),k))]
	allposs = [[I[a] for a in b] for b in indices]
	WLDlist=[]
	for I in allposs:
		if noncrossing(I) == True:
			I=[sorted(i+[a%n+1 for a in i]) for i in I]
			WLDlist+=[(n,I)]
		else:
			pass
	return WLDlist

