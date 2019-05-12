import os

# We're working in Gr(2,6) throughout
k,n = 2,6

# Generate a list of all 2,6 Le diagrams, sorted by dimension
L = allLe(k,n,True)

# Get the reduced permutation associated to each Le diagram
LRP = {}
for i in range(len(L)):
	LRP[i] = [LetoRP(M,n) for M in L[i]]



# We would like the 6-dim diagrams in sensible WLD order, so we define them individually
# Diagram naming matches that given in the paper

V1 = Matrix([[0,1,0,1],[1,1,1,1]])
V2 = Matrix([[1,0,1,1],[0,1,1,1]])
V3 = Matrix([[0,1,1,1],[1,1,1,x]])
V4 = Matrix([[1,1,1,0],[0,1,1,1]])
V5 = Matrix([[1,1,0,1],[1,1,1,x]])
V6 = Matrix([[1,0,1,0],[1,1,1,1]])

P1 = Matrix([[1,1,0,1],[0,1,1,1]])
P2 = Matrix([[1,0,1,1],[1,1,1,x]])
P3 = Matrix([[0,1,1,0],[1,1,1,1]])

E1 = Matrix([[1,1,1,x],[1,1,1,x]])
E2 = Matrix([[1,1,1,0],[1,1,1,x]])
E3 = Matrix([[1,1,1,0],[1,1,1,0]])
E4 = Matrix([[1,1,0,1],[1,1,0,1]])
E5 = Matrix([[1,0,1,1],[1,0,1,1]])
E6 = Matrix([[0,1,1,1],[0,1,1,1]])

# Put the WLD into a few lists so we can loop over them easily
VWLD = [V1,V2,V3,V4,V5,V6]
PWLD = [P1,P2,P3]
EWLD = [E1,E2,E3,E4,E5,E6]

# Collect the WLD and nonWLD Le diagrams
WLe = VWLD + PWLD + EWLD
NWLe = [M for M in L[n] if M not in WLe]

# Get the reduced permutations of each list (this is where positroids as objects would be so much nicer!)
WRP = [LetoRP(M,n) for M in WLe]
NWRP = [LetoRP(M,n) for M in NWLe]


L[6] = WLe + NWLe
LRP[6] = [LetoRP(M,n) for M in L[6]]


#####

# Now the computations!
# First, we run through the nonWLD diagrams and pick out all the cells on their boundaries
# Everything here is formatted so you can just include the output file in your tex code.


boundsNW = []
for r in range(len(NWLe)):
	bounds = {n: [NWLe[r]]}
	for i in range(n-1,-1,-1):
		res = []
		for j in range(len(LRP[i])):
			if isbound(NWRP[r][0],NWRP[r][1],LRP[i][j][0],LRP[i][j][1],k):
				res +=[L[i][j]]
		bounds[i] = res[:]
	boundsNW+=[bounds]



writepath = 'NWboundsdata.txt' # go ahead and modify this for wherever you want to save stuff

f = open(writepath,"w")
f.close()


for a in range(len(boundsNW)):
	res = "{\\bf Dimension "+str(n) +"} \\newline" + LetoLatex(boundsNW[a][n][0])
	for b in range(len(boundsNW[a])-2,-1,-1):
		res += "\n \n {\\bf Dimension "+str(b) + "} \\newline \n"
		for c in range(len(boundsNW[a][b])):
			res +=  LetoLatex(boundsNW[a][b][c]) + " \\quad "
		res += " \n \n "
	res +=" \\rule{\\textwidth}{.4pt}\\newpage \n \n \n "
	f = open(writepath,"a")
	f.write(res)
	f.close()




##### now the WLD ones



boundsW = []
for r in range(len(WLe)):
	bounds = {n: [WLe[r]]}
	for i in range(n-1,-1,-1):
		res = []
		for j in range(len(LRP[i])):
			if isbound(WRP[r][0],WRP[r][1],LRP[i][j][0],LRP[i][j][1],k):
				res +=[L[i][j]]
		bounds[i] = res[:]
	boundsW+=[bounds]



writepath = 'Wboundsdata.txt'

f = open(writepath,"w")
f.close()


for a in range(len(boundsW)):
	res = ""
	if a in range(6):
		res += "{\Large $V_{" + str(a +1)  + "}$} \\qquad $\\VV{" + str(a+1) + "}$  \\qquad "
	elif a in range(6,9):
		res += "{\Large $P_{" + str(a -6 +1)  + "}$} \\qquad $\\PP{" + str(a-6+1) + "}$  \\qquad "
	else:
		res += "{\Large $E_{" + str(a -9 +1)  + "}$} \\qquad $\\EE{" + str(a-9+1) + "}$  \\qquad "
	res += LetoLatex(boundsW[a][n][0])
	for b in range(len(boundsW[a])-2,-1,-1):
		res += "\n \n {\\bf Dimension "+str(b) + "} \\newline \n"
		for c in range(len(boundsW[a][b])):
			res +=  LetoLatex(boundsW[a][b][c]) + " \\quad "
		res += " \n \n "
	res +=" \\rule{\\textwidth}{.4pt}\\newpage \n \n \n "
	f = open(writepath,"a")
	f.write(res)
	f.close()







# let's generate all of the top-dim boundaries of every cell


writepath = 'topbounds.txt'

f = open(writepath,"w")
f.close()



for a in range(n,0,-1):
	res = "\\newpage \\section*{Dimension "+str(a)+"}"
	res += "\\begin{longtable}{| p{.20\\textwidth} | p{.80\\textwidth} |} \\hline \n"
	for b in range(len(LRP[a])):
		res+= "\\T" + LetoLatex(L[a][b]) + " & \n"
		for c in range(len(LRP[a-1])):
			if isbound(LRP[a][b][0],LRP[a][b][1],LRP[a-1][c][0],LRP[a-1][c][1],k):
				res += LetoLatex(L[a-1][c]) + "\\quad "
		res += " \\B \\\\  \\hline \n \n "
	res+="\\end{longtable}"
	f = open(writepath,"a")
	f.write(res)
	f.close()






#### Now to compute homology!


#**************************************
# First, some extra functions
#Generate a list of all unique Le diagrams associated to WLDs
#i.e. duplicates from exact subdiagrams removed

def Topchain(k,n):
	imagelist = []
	for i in range(len(allWLD(k,n))):
		imagelist.append(WLDtoLe(allWLD(k,n)[i]))
	generators =underlyingset(imagelist)
	return generators


# use following code as a cheat if we want to eliminate V diagrams from computation

allLelist = allLe(2,6,True)

topdim = allLelist[6]
# topdimRP = [LetoRP(L,6) for L in topdim]

# RP5 = [LetoRP(L,6) for L in allLelist[5]]

# # there's no system to this, I sorted them by hand
# V = [1,7,14,9,16,4]
# E = [20,17,13,12,10,6]
# P = [3,8,15]
# N = [0,2,5,11,18,19]

# def Topchain(k,n):
# 	return [topdim[i] for i in E + P]



#**************************************
#Create chain complex

# compute RP of a list of Le diagrams
def chainrp(chain,n):
	pairlist = []
	for i in chain:
		pairlist.append(LetoRP(i,n))
	#pairlist = list(map(lambda x: LetoRP(x,n), chain))
	return pairlist

#print "chainrp = ", chainrp(Topchain(k,n))
#print "allLerp = ", chainrp(allLe(k,n, dim=False))
#print chainrp(allLe(k,n, dim=False))[0][0]


def IDSharedBounds(prev,L,k,n):
	RPl = LetoRP(L,n)
	for M in prev:
		if isbound(M[0],M[1],RPl[0],RPl[1],k):
			return True
	return False


# starts with Topchain(k,n) and collects all lower-dim Le diagrams in the closure
# i.e. if L is in dimension i, all i-1 diagrams on the boundary of L are added.

def Chaincomplex(k,n,top,d):
	complexarray= {} # dictionaries are better than lists for this
	allLelist = allLe(k,n,True)
	complexarray[d]= top
	#allLeRP = {i:[LetoRP(M,n) for M in allLelist[i]] for i in range(d)}
	for i in range(d-1,-1,-1): 
		prev = [LetoRP(L,n) for L in complexarray[i+1]]
		complexarray[i] = [L for L in allLelist[i] if IDSharedBounds(prev,L,k,n)==True]
	return complexarray



# N = [0,2,5,11,18,19]



#******************************************
#Create differential maps in Z/2



def DiffMatrices(k,n,C,d):
	complexarray = C
	DiffMatrixlist = {}
	for i in range(d,0,-1):
		DM = zeros(len(complexarray[i-1]),len(complexarray[i]))
		for c in range(0, len(complexarray[i])):
			for r in range(0, len(complexarray[i-1])):
				source = complexarray[i][c]
				target = complexarray[i-1][r]
				if isbound(LetoRP(source,n)[0],LetoRP(source,n)[1],LetoRP(target,n)[0],LetoRP(target,n)[1],k):
					DM[r,c] = 1
		DiffMatrixlist[i] = DM[:,:]
	return DiffMatrixlist









# def DiffMatricesConcise(k,n):
# 	d = 3*k
# 	complexarray = Chaincomplex(k,n)
# 	DiffMatrixlist = [None]*(d+2)
# 	DiffMatrixlist[d+1] = zeros(len(complexarray[d]),1)
# 	DiffMatrixlist[0] = zeros(1, len(complexarray[0]))
# 	for i in range(0,d):
# 		DiffMatrixlist[d-i]= Matrix(len(complexarray[d-i-1]),len(complexarray[d-i]), lambda r,c: int(isbound(LetoRP(complexarray[d-i][c],n)[0],LetoRP(complexarray[d-i][c],n)[1],LetoRP(complexarray[d-i-1][r],n)[0],LetoRP(complexarray[d-i-1][r],n)[1],k)))
# 	return DiffMatrixlist




#***************************************
# Now we need to export this in a format that Sage can use

import os
cwd = os.getcwd()

# change this line to tell Python where to save the file
writepath = cwd + '/complex.txt'


# where C is the diffmatrices
def ExportDiffMatrices(C,writepath):
	f = open(writepath,"w")
	f.close()
	res = "{"
	for d in C.keys():
		R = [flatten(C[d][i,:]) for i in range(C[d].rows)]
		res+=str(d) + ":["
		for j in range(len(R)):
			res += str(R[j])+ ","
		res = res[:-1]+"],"
	res = res[:-1] + "}"
	f=open(writepath,"a")
	f.write(res)
	f.close()



# now upload this textfile to Sage, then uncomment and run the following in Sage:

# with open(DATA+'complex2.txt') as f:
#     C = eval(f.read())
   
# C1 = {i:matrix(GF(2),C[i]) for i in C.keys()}
# C2 = ChainComplex(C1,degree_of_differential=-1)
# ascii_art(C2)
# C2.homology(base_ring = GF(2))

#this will give you a list of your diff matrices, viewed as matrices over Z/2Z, in the list C1
# 






