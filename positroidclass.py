# main stuff to do:

# printing for positroid, permutation classes
# finish setting up the positroid init so we can use any rep of the positroid
# finish LetoGN
# need DPtoBP (BP = bruhat pair)
# some nice functions to convert to LaTeX



class Positroid:
	def __init__(self,positroid):
		if type(positroid) == Matrix and IsLe(positroid):
			self.le = positroid
			self.gn = LetoGN(self.le) 
			self.dp = GNtoDP(self.gn)
			self.bp = DPtoBP(self.dp)
		elif type(positroid) == list and IsGN(positroid):
			self.gn = positroid
			self.le = GNtoLe(self.gn)
			self.dp = GNtoDP(self.gn)
			self.bp = DPtoBP(self.dp)




# grassmann necklace to decorated permutation 
# we use the convention that if I_{i+1} = (I_i \ \{i\}) \cup \{j\} then \pi(i) = j
	def GNtoDP(I):
		n=len(I)-1
		P=['dp']+[0]*n
		decorations = []
		for i in range(1,len(I)):
			J=SetDifference(I[(i%n)+1],I[i])
			if J==[]:
				# fixed point, need to know if it's a loop (not in any I) or a coloop (in all I)
				P[i]=i
				if i in I[i]:
					decorations+=[(i,1)]
				else:
					decorations+=[(i,-1)]
			else:
				P[i]=J[0]
		return Dperm(P,decorations)

	# checks whether a matrix represents a Le diagram 
	# (non-rectangular Le diagrams use 'x' to indicate missing squares)
	def isLe(M):
		m,n=M.rows,M.cols
		# pick out all the squares in M which contain a zero
		zero_index=[i for i in range(len(M)) if M[i]==0]
		for i in zero_index:
			# check the Le condition on each zero square
			if flatten(M[i//n,:i%n]).count(1) > 0 and flatten(M[:i//n,i%n]).count(1) >0:
				return False
		return True



	def isGN(I):
		if type(I[0])==str:
			I1 = I[:]
		else:
			I1 = ['']+I[:]
		n = len(I1) - 1
		# I is GN iff I_{i+1} \supseteq I_i \ \{i\} for all i in [n]
		for i in range(1,len(I1)):
			if i in I1[i]:
				if len([j for j in I1[i] if j!=i and j in I1[i%n+1]])<len(I1[i])-1:
					return False
				else:
					pass
			else:
				if len(Intersection(I1[i],I1[i%n+1])) < len(I1[i]):
					return False
				else:
					pass
		return True


# to get the GN from the Le diagram, we need I1 (the shape of the Young diagram)
# and then apply Suho Oh's algorithm to extract the remaining terms

	def getYoungShape(M):
		m,n = M.rows,M.cols
		I1,count = [],1
		rowpos,colpos = 0,n-1
		while rowpos <= m-1 and colpos >= 0:
			if M[rowpos,colpos]=='x':
				colpos+=(-1)
			else:
				rowpos+=1
				I1+=[c]
			c+=1
		if colpos < 0 and c <= m+n:
			I1+=range(c,m+n+1)
		return I1



# the chain rooted at square p is the unique path NW through the diagram starting at p,
# where at each step we take a minimal step strictly NW

# p is a tuple, (row,col)
	def getChain(M,p):
		isblack = lambda A:flatten(A).count(1)==0 
		if M[p] == x:
			return []
		elif M[p] == 0:
			M1=M[0:p[0]+1,0:p[1]+1]
			if isblack(M1) == True:
				return []
			else:
				boxes = [divmod(i,M1.cols) for i in range(len(M1))][::-1]
				for r in boxes:
					if M[r]==1:
						chain=[r]
						break
		else:
			chain = [p]
		while isblack(M[0:chain[-1][0],0:chain[-1][1]]) == False:
			M1=M[0:chain[-1][0],0:chain[-1][1]]
			boxes = [divmod(i,M1.cols) for i in range(len(M1))][::-1]
			for r in boxes:
				if M[r]==1:
					chain+=[r]
					break
		return chain


	# the GN is the chains rooted at the squares on the SE boundary of the diagram
	def LetoGN(M):
		m,n = M.rows,M.cols
		GN = ['gn',getYoungShape(M)]
		# we need a list of squares on the SE border of M
		# NEED TO FIX THIS FOR NON-RECTANGULAR DIAGRAMS
		SE = [(n-1) + i*n for i in range(0,m-1)] + range((m-1)*n,m*n)[::-1]
		for i in SE:
			GN+=[getChain(M,divmod(i,n))]
		return GN









class Permutation:
	def __init__(self,perm):
		# a str in the 0th position allows us to index more naturally
		if type(perm[0])== str:
			self.perm = perm
			self.dim = len(perm) - 1
		else:
			self.perm = ['p']+list(perm) 
			self.dim = len(perm)
	# we compose permutations right to left	
	def __mul__(self,other):
		assert self.len == other.len, "Permutations must have same length"
		return Permutation([self.perm[other.perm[i]] for i in range(1,self.dim+1)])
# every permutation can be expressed as a reduced word of elementary transpositions
# this works for any perm, and does not necessarily produce the RW coming from the Le diagram
	def computeRW(self):
		P1=self.perm[1:]
		T=[]
		for i in range(1,len(P1)):
			j=i
			while P1[j] < P1[j-1]:
				temp=P1[j]
				P1[j]=P1[j-1]
				P1[j-1]=temp
				T+=[j]
				j+=-1
				if j==0:
					break
		self.RW = ['rw']+T[::-1]
		self.RWlen = len(T)
		return self.RW











# subclass for decorated permutations
class DPerm(Permutation):
	def __init__(self,perm,decorations = []):
		Permutation.__init__(self,perm)
		self.perm[0] = 'dp' 
		# expect decorations to be of the form [(i,\pm1),(j,\pm1),...]
		self.decorations = decorations



# class RWPerm(Permutation):
# 	def RW(P):
# 		P1=P[1:]
# 		T=[]
# 		for i in range(1,len(P)):
# 			j=i
# 			while P1[j] < P1[j-1]:
# 				temp=P1[j]
# 				P1[j]=P1[j-1]
# 				P1[j-1]=temp
# 				T+=[j]
# 				j+=-1
# 				if j==0:
# 					break
# 		return T[::-1]
# 	def __init__(self,perm):
# 		Permutation.__init__(self,perm)
# 		self.RWexpr = RW(perm)
# 		#self.length = len(self.RWexpr)

