#New Program
'''Here we utilize the fact that, given a maximal tubing T of a poset P,
the union of tubes in T excludes exactly one element in a nontrivial bundle.
We use this fact to compute all maximal tubings of P. We know that every tubing
is a subset of some maximal tubing. Thus, each tubing is an element of the power
set of some maximal tubing.'''

'''This function takes a poset, and returns the list
[maximal_tubings,list_of_tubes,tube_bundle_dict]. Here, maximal_tubings is the
complete list of maximal tubings. list_of_tubes is the complete list of tubes,
and tube_bundle_dict is a dictionary of lists, indexed by tubes of the poset.
The entry for a given tube is the output from add_tube_bundle_info: a list whose
0th entry is the list of nontrivial bundles that intersect the tube nontrivially
and whose 1st entry is the union of entries in the above list. 
'''
def compute_maximal_tubings(P,show=false):
	P=Poset(P,facade=true)
	[boundaries,bundles,nontrivial_bundles]=elt_bundles(P)
	nontrivial_elements=s_o_s_union(Set(nontrivial_bundles))
	dim = nontrivial_elements.cardinality()-len(nontrivial_bundles)
	if dim==0:
		print '0-dimensional'
		return [Set([P])]
	print 'computing vertices of',dim,'dimensional polytope'
	n_t_up={}
	for m in nontrivial_elements:
	   n_t_up[m]=Set(P.order_filter([m])).intersection(nontrivial_elements)

	####### Ready to construct first round of partial tubings! #######

	partial_tubings=[]
	maximal_tubings=[]
	for nt in nontrivial_elements:
		if n_t_up[nt].cardinality()==1:
			rem = P.order_filter([nt])
			t_poset=P.list()
			for x in rem:
				t_poset.remove(x)
			t_poset=P.subposet(t_poset)
			CComp=CC(t_poset,P)
			if CComp.cardinality()<dim:
				partial_tubings.append(CComp)
			else:
				maximal_tubings.append(CComp)
	if show==true:
		print 'First round complete!'

	############ PREP FOR SUBSEQUENT ITERATIONS ############

	list_of_tubes=[]
	tube_bundle_dict={}
	for Tubing in partial_tubings:
		for tube in Tubing:
			list_of_tubes.append(tube)
			tube_bundle_dict[tube]=add_tube_bundle_info(tube,nontrivial_bundles)
	if show==true:
		print 'Prep for phase 2 complete!'

	############ REMAINING: MANY ITERATIONS ############

	count=1
	while len(partial_tubings)>0:
		for Tubing in partial_tubings:
			for tube in minimal_tubes_in_tubing(Tubing):
				new_Tubings = add_CC_tubings(Tubing,tube,P,tube_bundle_dict,n_t_up)
				for new_tubing in new_Tubings:
					if new_tubing.cardinality()<dim:
						partial_tubings.append(new_tubing)
					else:
						maximal_tubings.append(new_tubing)
					for new_tube in new_tubing:
						if list_of_tubes.count(new_tube)==0:
							list_of_tubes.append(new_tube)
							tube_bundle_dict[new_tube] = \
								add_tube_bundle_info(new_tube,nontrivial_bundles)
							
			partial_tubings.remove(Tubing)

		count=count+1
		if show==true:
			print 'Done with round',count,'now',len(partial_tubings),'partial tubings left!'
	return [maximal_tubings,list_of_tubes,tube_bundle_dict]

'''
This function takes a poset, and returns all tubings as a list of lists, so that
list[k] is the list of all k-tubings of the poset. 
'''
def compute_all_tubings(P,show=false):
	start=compute_maximal_tubings(P,show)
	if show==true:
		print 'Vertices acquired, now construcing all tubings'
	dim=start[0][0].cardinality()
	if show==true:
		print 'dimension is',dim
	Tubings=[]
	for i in range(0,dim+1):
		Tubings.append([])
	Tubings[dim]=start[0]
	for max in Tubings[dim]:
		for k in range(0,dim):
			add=Subsets(max,k)
			for sub in add:
				Tubings[k].append(sub)
	for ind in range(0,len(Tubings)):
		Tubings[ind]=Set(Tubings[ind])
	return Tubings

'''
Given a complete list of tubings, return the face-poset of the corresponding
polytope (the lattice of tubings, ordered by reverse inclusion)
'''
def KP_Poset(tubings):
    fcn = lambda p,q: q.issubset(p)
    ret = Poset([s_o_s_union(Set(tubings)),fcn])
    return Poset(ret,facade=true)

'''
This method take a Poset as the input, and outputs a list. list[1] is a 
dictionary of all bundles of P, indexed by list[0], a list of boundaries of
the bundles. list[2] is a list of all nontrivial bundles of P.
'''
def elt_bundles(P):
    bundles={}
    nontrivial_bundles=[]
    for x in P:
        B=Set(P.order_ideal(P.lower_covers(x)))
        if bundles.keys().count(B)==0:
            bundles[B]=Set([x])
        else: 
            bundles[B]=bundles[B]|Set([x])
    for key in bundles.keys():
        if bundles[key].cardinality()>1:
            nontrivial_bundles.append(bundles[key])
    return [Set(bundles.keys()),bundles,nontrivial_bundles]
##RUNTME: O(P)

'''
This function takes a Set Of Sets (sos), and returns the Set of all elements
contained in at least one of the Sets
'''
def s_o_s_union(sos):
	ret = Set([])
	for s in sos:
		ret = ret.union(s)
	return ret
##RUNTIME: O(sos)

'''
This function takes a subset S of elements of a poset P, and returns
the Set of connected component of P that contains S. When used in the main 
functions, S is always a tube, hence connected.
'''
def CC(S,P):
#Pass S, subposet of P, return connected components of S
    ret = []
    have = Set([])
    flag = true
    while flag:
        comp = []
        comp2 = [Set(S.list()).difference(Set(have)).an_element()]
        while len(comp)<len(comp2):
            comp=comp2
            comp2=grow(comp,S)
            have=have.union(Set(comp2))
        flag = have.cardinality()<S.cardinality()
        ret.append(Set(comp))
    return Set(ret)

'''
This function takes a subset s of elements of the poset P, and returns
a set in the connected component of s that contains s. It does this by
adding the upper set (order filer) of s, and then the lower set (order ideal)
of this new set.
'''
def grow(s,P):
    ret = P.order_filter(s)
    ret = P.order_ideal(ret)
    return ret

'''
This function takes a tubing T, a constituent tube t in T, a poset P, the
dictionary of tube bundles, and a dictionary indexed non-trivial elements, whose
value at a non-trivial element nt is the set of non-trivial elements 
greater than or equal to nt in P. The output is the list of all tubings obtained
by adding the connected components of t
'''
def add_CC_tubings(T,t,P,tube_bundle_dict,n_t_up):
	#Need to add CC(t \ x^up) for x maximal among nontrivial t-elements
	max = []
	new_Tubings=[]
	[nt_t_bundles,nt_t_elements] = tube_bundle_dict[t]
	for nt in nt_t_elements:
		if n_t_up[nt].intersection(nt_t_elements).cardinality()==1:
			max.append(nt)
	for m in max:
		old_T = T.list()
		#old_T.remove(t)
		new_T=old_T+CC_of_removal(t,m,P).list()
		new_Tubings.append(Set(new_T))
	return new_Tubings

'''
This function takes a tube of a poset and a list of nontrivial bundles of the
poset. It returns a list whose 0th entry is the a list of non-trivial bundles 
of the tube, and whose 1st entry is the union of said bundles. 
'''
def add_tube_bundle_info(tube,nontrivial_bundles):
	tube_bundles=[]
	#list_of_tubes.append(tube)
	for nt_bundle in nontrivial_bundles:
		inter=nt_bundle.intersection(tube)
		if inter.cardinality()>1:
			tube_bundles.append(inter)	
	return [tube_bundles,s_o_s_union(Set(tube_bundles))]	

'''
This functions takes a tube t (as a Set) of the Poset P, as well as an element
in t (in general, maximal among elements in non-trivial t-bundles). The output
is the Set of connected components of t minus the upper set of nt.
'''
def CC_of_removal(t,nt,P):
	t_poset=P.subposet(t.list())
	up=t_poset.order_filter([nt])
	new_list = t.list()
	for x in up:
		new_list.remove(x)	
	CComp=CC(t_poset.subposet(new_list),t_poset)
	return CComp

'''
This function takes a tubing (a list of poset elements), recreates the
subposet induced by the tubing, and returns a list of minimal elements
'''
def minimal_tubes_in_tubing(T):
	fcn = lambda p,q: p.issubset(q)
	temp=Poset([T.list(),fcn])
	temp=Poset(temp,facade=true)
	return temp.minimal_elements()

##### This function takes a python Graph, and outputs a Poset #####
def graph_to_poset(G):
	V=G.vertices()
	E=G.edges()
	elms=[]
	rels=[]
	for v in V:
		elms.append(Set([v]))
	for e in E:
		eset=Set([e[0],e[1]])
		elms.append(eset)
		for v in V:
			vset=Set([v])
			if vset.issubset(eset):
				rels.append([vset,eset])
	P=Poset((elms,rels))
	return P

