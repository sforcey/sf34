def nestohedron_poset(S,B):
    '''Given a building set B that is a subset of 2^S, construct a poset
    whose minimal elements are singletons of S, while maximal
    elements are those of B which are not singletons
    '''
    elts = B
    fcn = lambda p,q: p.issubset(q) and p.cardinality()==1
    P=Poset((elts,fcn),facade=true)
    return P

def rank_2_bundled_posets_max(dim,max):
    '''
    This function constructs all isomorphism classes of rank 2 posets whose 
    associahedra are #dim-dimensional, and have exactly #max bundles of maximal
    elements. The output is a list, whose entries are lists of the form 
    [elements of poset (as a Set), poset].
    '''
    S = range(1,dim+2)
    singletons=[]
    for x in S:
        singletons.append(Set([x]))
    powerset = Subsets(S)
    possible_boundaries=[x for x in powerset if x.cardinality()>0]
    exp=len(possible_boundaries)
    print powerset.cardinality(),'sets in powerset and',\
    exp,'boundaries possible!'
    Buildings=Subsets(possible_boundaries,max)
    print 'Will check',Buildings.cardinality(),'posets...'
    Buildings2=[singletons+B.list() for B in Buildings]
    for B in Buildings2:
        for delta in B:
            if B.count(delta)>1:
                B.remove(delta)
                B.append( Set([ delta[0],Set([]) ]) )
    for B in Buildings2:
        for delta in B:
            if B.count(delta)>1:
                B.remove(delta)
                B.append( Set([ delta[0],Set([]) ]) )            
    unique = [[Buildings2[0], nestohedron_poset(S,Buildings2[0])]]
    count=1
    for B in Buildings2:
        BP = nestohedron_poset(S,B)
        flag=true
        for BU in unique:
            if BP.is_isomorphic(BU[1]):
                flag=false
                break
        if flag==true:
            unique.append([B,BP])
            if len(unique)%24==0:
                print 'Max =',max,'Checked',count,'of',binomial(exp,max),\
                'and found',len(unique),'unique posets so far...'
        count+=1
    print len(unique),'posets found!'
    return unique


def duplicate_elements(S,B,delta):
    '''
    Given an subset B of the powerset 2^S and a list delta of elements in B that
    are to be "duplicated", this returns the list which functions as the
    concatenation B+delta. 
    '''
    start=S.cardinality()+1
    add=[]
    for i in range(0,len(delta)):
        add.append( delta[i].union( Set([start+i]) ) )
    return B+add

def SubMultisets(L,k):
    '''
    Given a list L, return the Cartesian product L^k as a list, whose entries
    are are elements of L^k in list format.
    '''
    temp=[[x] for x in L]
    temp2=[]
    if k==1:
        return temp
    for i in range(1,k):
        for x in temp:
            for y in L:
                temp2.append(x+[y])
        temp=temp2
        temp2=[]
    return temp
    
def add_dimension(S,B,dim_plus):
    '''
    Given the ground Set S, a particular B in 2^(2^S), and a positive integer
    dim_plus, returns a list of lists [B', PB'] where B' is the set of elements
    of a poset that, when the #dim_plus maximal non-trivial bundles (with 
    respect to the poset order) are collapsed, the result is just B, while BP'
    is the poset.
    '''
    duplicatable=[]
    for b in B:
        if b.cardinality()>1:
            duplicatable.append(b)

    duplicates=SubMultisets(duplicatable,dim_plus)
    start_set=duplicate_elements(S,B,duplicates[0])
    start_poset=nestohedron_poset(S,start_set)
    newly_bundled=[[start_set,start_poset]]
    for delta in duplicates:
        B2=duplicate_elements(S,B,delta)
        BP=nestohedron_poset(S,B2)
        for BU in newly_bundled:
            flag=true
            if BP.is_isomorphic(BU[1]):
                flag=false
                break
        if flag==true:
            newly_bundled.append([B2,BP])    
    return newly_bundled

def add_dimension_to_list(S,L,dim_plus):
	unique=[]
	count=len(unique)
	for B in L:
		temp=add_dimension(S,B[0],dim_plus)
		for x in temp:
			count+=1
			unique.append(x)
			if len(unique)%25==0:
				print 'Given',len(L),'found',len(unique),'by poset #',\
				count,'of',len(L)*(len(L[0][0])-S.cardinality())+1
	print 'Found',len(unique),'total, given',len(L),'to start'
	return unique

def rank_2_associahedra_sets(dim,min,Lstart):
    '''This command will give you a list of lists [E,P] where E is the list of
    elements of a poset P, with P ordered by subset containment, such that the
    dimension of the poset associahedron KP is of dimension dim, where P also 
    has exactly min minimal elements. Lstart should be a list of consecutive
    nonegative integers, whose maximum value must be at most 2^(2^dim)). This
    list is the cardinality of maximal elements of the posets, before 
    'duplication'.

    The output of this command is a list of lists, where the list of index i is
    a list of posets with i maximal elements before duplication, although the 
    output will be post duplication.
    '''
    S=Set( range(1,min+1) )
    start=Lstart[0]
    end=Lstart[-1]
    prebundled = [ [] for i in range(0,start) ]
    bundled = [ [] for i in range(0,start) ]
    for k in Lstart:
	prebundled.append(rank_2_bundled_posets_max(min-1,k))
	#bundled=[[] for i in range(0,start)]
    for k in Lstart:
	if k == 0:
	   bundled.append([])
	   continue
	L=prebundled[k]
	print 'Starting with',len(L),'may find up to',(dim+1-min)*len(L)*k+1
	time bundled.append(add_dimension_to_list(S,L,dim+1-min))
    return bundled