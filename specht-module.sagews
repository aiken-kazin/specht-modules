︠23085164-5457-48d2-9010-e5aeeabaac3bs︠
def StandardTableauxDominancePoset(lam):
    #This poset is used for triangularity in the Specht Module.
    def composition_tuple(s):
        #s a standard tableau of shape lam.
        #last composition is always lam, so no need to include it here...
        Mu=[]
        for i in range(lam.size()):
            Mui = [ len([j for j in r if j<=i+1]) for r in list(s)]
            Mu.append( Composition( Mui ) )
        return tuple(Mu)

    def dominates(s,t):
        #Checks if s dominates t.
        cs=composition_tuple(s)
        ct=composition_tuple(t)
        for i in range(len(cs)):
            cs_partial_sums=cs[i].partial_sums()
            ct_partial_sums=ct[i].partial_sums()
            #print i, 't', cs_partial_sums, 't', ct_partial_sums
            for j in range(len(cs_partial_sums)):
                if not cs_partial_sums[j]>=ct_partial_sums[j]: return False
        return True

    st=StandardTableaux(lam)
    G=DiGraph()
    if st.cardinality()==1:
        G.add_vertex(st[0])
        return Poset(G)
    for s in st:
        for t in st:
            if dominates(s,t): 
                G.add_edge((s,t))
                #print s,t
    return Poset(G)

class SpechtModule(CombinatorialFreeModule):
    def __init__(self, shape, R=QQ):
        self.shape=Partition(shape)
        self.size=shape.size()
        self.S=SymmetricGroup(self.size)
        self.A=self.S.algebra()
        self.st=StandardTableaux(lam)
        L=StandardTableauxDominancePoset(lam).linear_extension()
        #L.reverse()
        self.ordered_st=[t.element for t in L]

        CombinatorialFreeModule.__init__(self, R, Tableaux(shape=shape))

    def __repr__(self):
        return "Tableau module of shape " + str(self.shape)

    def _permutation_action_on_tableau(self, p, t):
        #t a tableau
        #p a permutation
        s=t.to_list()
        n=t.size()
        for i in [1..n]:
            for c in t.cells_containing(i):
                s[c[0]][c[1]]=p(i)
        return self(Tableau(s))

    def _algebra_action_on_tableau(self, x, t):
        #x an element of the group algebra
        #t a tableau
        action=self._permutation_action_on_tableau
        return sum([x.coefficient(p)*action(p,t) for p in x.support()])

    def action(self, x, t):
        #x an element of the Sn group algebra
        #t an element of self.
        x=x.parent().algebra()(x)
        algAct=self._algebra_action_on_tableau
        return sum([t.coefficient(s)*algAct(x, s) for s in t.support()])

    @cached_method
    def polytabloid(self, t):
        cs=t.column_stabilizer()
        ct=sum([a.sign()*self.A(a) for a in cs])
        rs=t.row_stabilizer()
        rt=sum([self.A(a) for a in rs])
        b=rt*ct
        return self._algebra_action_on_tableau(b, t)

    @cached_method
    def spechtBasis(self):
        return [self.polytabloid(Tableau(t)) for t in self.ordered_st]

    def gens(self):
        return self.spechtBasis()

    def matrix_of_permutation_on_basis(self,p):
        m=[]
        for b in self.spechtBasis():
            a=self.action(p,b)
            coeffs=[]
            for t in self.ordered_st:
                #print t, a.coefficient(t)
                c=a.coefficient(t)
                coeffs.append(c)
                if c!=0: a=a-c*self.polytabloid(t)
            assert a==0, "Something's gone terribly wrong."
            m.append(coeffs[:])
        #m.reverse()
        return matrix(m)

    def character(self):
        ch=[]
        for pi in self.S.conjugacy_classes_representatives():
            ch.append( self.matrix_of_permutation_on_basis(pi).trace() )
        return ch

    def fourierTransform(self, f):
        #f a function from the permutation group to the base field.
        mat=self.matrix_of_permutation_on_basis
        return sum([f(p)*mat(p) for p in self.S])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

n=lam.size()
S=SymmetricGroup(n)
x=var('x')

def P_reflect(s):
    assert s in S, "Need permutations. "+str(s) + str(S)
    if s==S.one(): return x
    if s^2==S.one(): return (1-x)/binomial(n,2)
    return 0

convolve = lambda P,Q: lambda s: sum([P(s.inverse()*t)*Q(t) for t in S])
︡a3bc0334-a5c5-46f4-adce-c6c46b7beb2b︡{"stderr":"Error in lines 86-86\nTraceback (most recent call last):\n  File \"/cocalc/lib/python3.10/site-packages/smc_sagews/sage_server.py\", line 1244, in execute\n    exec(\n  File \"\", line 1, in <module>\nNameError: name 'lam' is not defined\n"}︡{"done":true}
︠57f1fcb4-58e8-4c88-8505-e1353d3a6eac︠









