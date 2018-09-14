#!/usr/bin/env python
"""
simp_sec_quant: SIMPLIFICATOR FOR SECOND-QUANTIZATION EXPRESSIONS
The simp_sec_quant.py file contains a set of classes (objects) that allows 
to construct second quantization expressions with a Fock operator
and a fluctuation potential and spin-free molecular orbitals (closed-shell).

The created expressions are similar to the one derived in Chapter 13 
of the Molecular electronic structure theory book by Helgaker.

The expressions can then be simplified to expressions that can be implemented
in standard quantum chemistry packages.

To see examples of usage look into ${comp_chem_py}/tests/expressions.py

TODO: 
   - add triple and quadruple permutations
   - add inner product with doubly excited determinant
"""

__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"

import itertools

iuns_orb = 0
def uns_orb():
    """ return new index for unspecified orbital and increase counter """
    global iuns_orb
    iuns_orb += 1
    return 'u_' + str( iuns_orb )

iocc_orb = 0
def occ_orb():
    """ return new index for occupied orbital and increase counter """
    global iocc_orb 
    iocc_orb += 1
    return 'o_' + str( iocc_orb )

ivir_orb = 0
def vir_orb():
    """ return new index for virtual orbital and increase counter """
    global ivir_orb
    ivir_orb += 1
    return 'v_' + str( ivir_orb )


#------------------------------------------------------------------------------
class expression(object):
    """ an expression is defined as a sum of terms"""

    def __init__(self, terms): 
        # terms must be a list of terms:
        if any(type(t) is not term for t in terms):
            raise Exception("wrong type in terms")
        self.terms = terms

    def simple(self):
        # expression is simplified if all the terms are simplified
        for t in self.terms:
            if not t.simple():
                return False
        return True

    def simplify(self):
        # simplify each term and return new expression
        newterms = []
        for t in self.terms:
            new = t.simplify()
            if type(new) is term:
                newterms.append( new )
            elif type(new) is expression:
                newterms += new.terms 
            else:
                raise Exception("this should not happen")
        return expression( newterms )

    def __mul__(self,other):
        """ an expression times something else gives a term"""
        return term( elmts=[self] ) * other

    def __add__(self,other):
        """ an expression plus something else gives an expression"""
        newexp = self
        newexp += other
        return newexp

    def __iadd__(self,other):
        """ an expression plus something else gives an expression"""
        if type(other) is expression:
            self.terms += other.terms
        elif type(other) is term:
            self.terms.append( other )
        else:
            # try making a term out of the object and add it
            self.terms.append( term( elmts=[other] ) )
        return self

    def output(self):
        """ return printable latex string corresponding to the expression"""
        return ' '.join( [t.output() for t in self.terms] ) 

    def permute(self, index):
        """ perform permutations for the terms in the expression """
        # each term returns a sum of terms, i.e. an expression
        newexp = expression([]) # empty expression
        for t in self.terms:
            newexp += t.permute( index ) 
        return newexp

    def rm_parenthesis(self):
        """ ditribute parenthesis present in individual terms"""
        newexp = expression([]) # empty expression
        for t in self.terms:
            newexp += t.rm_parenthesis()
        return newexp

    def evaluate(self):
        """ use boxes from the bible and logic to evaluate the individual terms"""
        # we are assuming a simplified expression
        if not self.simple():
            raise Exception("this should not happen")

        newexp = expression([]) # empty expression
        for t in self.terms:
            # check if the term vanishes due to rank reduction
            if not t.vanish():
                # evaluate terms that do not vanish
                # i.e. get expression for the commutator on HF
                newexp += t.evaluate()

        return newexp

#------------------------------------------------------------------------------
class term(object):
    """ A term is a list of elmts multiplied together, with a sign in front (+ or -)."""

    def __init__(self, sign='+', elmts=[]): 
        self.sign = sign
        self.elmts = elmts

        # check validity of expression
        if self.sign is not '+' and  self.sign is not '-':
            print "sign is",sign
            raise Exception("wrong sign for term")

        # the allowed elements can be of the following types:
        allowed_types = [int, float, summation, tensor, operator, el_oper, state, commutator, expression]
        if any( type(elmt) not in allowed_types for elmt in self.elmts ):
            raise Exception("forbidden type in term")


    def simple(self):
        """ simple is true for simplified terms"""
        for i, elmt in enumerate(self.elmts):

            # each element in the term need to be simplified
            if type(elmt) not in simple_types:
                if not elmt.simple():
                    # term cannot be simple
                    return False

            # and commutators cannot be followed by elementary operators
            if type(elmt) is el_oper and i > 0:
                if type(self.elmts[i-1]) is commutator:
                    return False

        return True


    def simplify(self):
        """ simplify each element and return the result """
        # for a term to be simplified, each element must be simplified
        # and commutators cannot be followed by elementary operators 

        # 1) simplify each element
        elmts1 = []
        for elmt in self.elmts:
            if type(elmt) in simple_types:
                # just copy the elelement
                elmts1.append(elmt)
            else:
                elmts1.append( elmt.simplify() )

        # build expression containing a single term 
        newter = [term( sign=self.sign, elmts=elmts1)]
        newexp = expression( newter )

        # 2) check what comes after the commutators and remove parenthesis
        modified = True
        while modified:
            modified = False

            # remove all parenthesis from newexp
            clean_exp = newexp.rm_parenthesis()

            newter = []
            for t in clean_exp.terms:

                elmts1 = []
                for i, elmt in enumerate(t.elmts):
                    if type(elmt) is el_oper and i > 0:
                        if type(t.elmts[i-1]) is commutator:
                            # replace the previous entries by the simplified expression
                            elmts1.pop() # delete commutator
                            exp = identity3( t.elmts[i-1], term(sign=t.sign, elmts=t.elmts[i:]) )
                            # add simplified expression to the term
                            elmts1.append( exp )
                            modified = True
                            break
                        else:
                            # just copy current element
                            elmts1.append( elmt )
                    else:
                        # just copy current element
                        elmts1.append( elmt )

                # update list of terms
                newter.append( term(sign=t.sign, elmts=elmts1) )

            # make new expression
            newexp = expression( newter )

        # on top of that we reorganize each term of the new expression
        newter = []
        for t in newexp.terms:
            newter.append( t.reorganize() )
        newexp = expression( newter )

        # now newexp should be a simplified expression without parenthesis
        if not newexp.simple():
            print "expression not simple", newexp.output()
            raise Exception("this should not happen")

        return newexp

    def reorganize(self):
        """ reorganize simplified term"""
        if not self.simple():
            raise Exception("this should not happen")
        # organize term should be as follows:
        #   1- factors
        #   2- summations
        #   3- tensors
        #   4- a single Bra state
        #   5- a series of elementary operator
        #   6- a simplified commutator
        #   7- a single Ket state
        
        numb = []
        summ = []
        tens = []
        matr = [] # elmts between bra-ket
        for elmt in self.elmts:

            if type(elmt) in [int, float]:
                numb.append( elmt )
            elif type(elmt) is summation:
                summ.append( elmt )
            elif type(elmt) is tensor:
                tens.append( elmt )
            elif type(elmt) in matrix_elmts:
                matr.append( elmt )
            else: 
                raise Exception("this should not happen")

        # multiply numbers together and have only one factor
        n = 1
        for i in numb:
            n *= i

        # build new organize term
        if n==1:
            # remove any unity factor
            newelmts = summ + tens + matr
        else:
            newelmts = [n] + summ + tens + matr

        return term(sign=self.sign, elmts=newelmts)


    def rm_parenthesis(self):
        """distribute whatever is inside parenthesis"""
        # the term contains parenthesis if it contains at least one expression
        if not any(type(elmt) is expression for elmt in self.elmts):
            return expression( [self] )

        # if the term contain only one element then this element is an expression
        # which has to be decomposed into different terms
        if len(self.elmts) <= 1:
            return self.elmts[0] # we just return the expression as such

        newterms = []
        # distribute each element with the rest if needed:
        for i, elmt in enumerate(self.elmts):

            if type(elmt) is not expression:
                # nothing to distribute in that elmt
                continue # deal with next elmt

            # (A1 + A2...)*elmt2*... = A1*elmt2*... + A2*elmt2*... + A3*...
            # so the current elmt will generate as many terms as it includes
            for t in elmt.terms:

                if i == 0:
                    # first element

                    newterm = t * term( sign=self.sign, elmts=self.elmts[i+1:])
                elif i == (len(self.elmts) - 1):
                    # last element
                    newterm = term( sign=self.sign, elmts=self.elmts[:i]) * t
                else: 
                    # middle element
                    # self.sign should be used only once!!!
                    newterm = term( sign=self.sign, elmts=self.elmts[:i]) * t * term( sign='+', elmts=self.elmts[i+1:])

                newterms.append( newterm )

        return expression( newterms )


    def output(self):
        """ return printable string for term """
        # start with sign
        string = []
        string.append(self.sign)

        # add other elmts
        for elmt in self.elmts:
            if type(elmt) is int or type(elmt) is float:
                string.append( str(elmt) )
            elif type(elmt) is expression:
                # we need to add parenthesis
                string.append( '( '+ elmt.output() +' )' )
            else:
                string.append( elmt.output() )
        return ' '.join(string)


    def __add__(self, other):
        """ a term plus something else gives an expression"""
        return expression( [self] ) + other


    def __mul__(self, other):
        """ multiplying a term with something else gives a term"""
        newelmts = self.elmts 
        if type(other) is term:
            newelmts += other.elmts
            if self.sign == other.sign:
                sign = '+'
            else:
                sign = '-'
        else:
            newelmts.append( other )
            sign = self.sign
        return term(sign=sign, elmts=newelmts ) 


    def permute(self, index):
        """ perform permutations of the indices  in the term"""
        # index is a list of tuples, e.g. index = [(a, i), (b, j)]
        # in that case we build a dictionary idmap that maps 
        # index a -> b, i -> j, b -> a, j -> i
        # and use it to build a new expression

        dim = len(index)
        # build map of indices as a dictionary
        exp = expression([])

        # loop over all possible permutations
        for perm in itertools.permutations(index, dim):
            idmap = {}
            for idx,p in zip(index, perm):
                idmap[idx[0]] = p[0]
                idmap[idx[1]] = p[1]

            # permute indices of each element
            new_elmts = []
            for elmt in self.elmts:
                if type(elmt) is float or type(elmt) is int:
                    # nothing to do here
                    new_elmts.append( elmt )
                else:
                    new_elmts.append( elmt.permute(idmap) )
             
            exp += term(sign=self.sign, elmts=new_elmts)

        return exp


    def vanish(self):
        """ use rank reduction to remove vanishing commutators"""
        # we are assuming a simplified term
        if not self.simple():
            raise Exception("this should not happen")

        # if one commutator is zero the whole term is zero
        vanish = False
        for elmt in self.elmts:
            if type(elmt) is commutator:
                vanish = elmt.vanish()
                if vanish:
                    break

        return vanish

    def evaluate(self):
        """ use the box from the bible to evaluate a term"""
        # the term to evaluate should consist of a section of the form:
        #  <BRA| el_oper * el_oper ...* [[operator, el_oper], el_oper] |HF>
        # with potential tensors and summations on the left side
        if not self.simple():
            raise Exception("this should not happen")

        # find the commutator evaluate it and return new term
        count = 0
        newelmts = []
        for i, elmt in enumerate(self.elmts):

            if type(elmt) is commutator:
                # this should happen only once if the term has the correct form
                count += 1
                if count > 1:
                    raise Exception("this should not happen")

                # evaluate the sequence [] |HF>
                newelmts.append( evaluate_comm_HF(elmt, self.elmts[i+1]) )
                break

            newelmts.append( elmt )
        
        # simplify the term to get expression
        newexp = term(sign=self.sign, elmts=newelmts).simplify()

        # finally contract the matrix elements using delta functions
        newexp2 = expression([])
        for t in newexp.terms:
            newterm = t.contract()
            if newterm != 0:
                newexp2 += t.contract() 

        return newexp2

    def contract(self):
        # check that matrix element is not zero

        # first get dummy indices matrix element and the rest
        dummy = []
        matrix = []
        rest = []
        for elmt in self.elmts:
            if type(elmt) is summation:
                dummy += elmt.index
                rest.append( elmt )
            elif type(elmt) in matrix_elmts:
                matrix.append( elmt )
            else:
                rest.append( elmt )


        # get delta from matrix element
        nel_oper = 0
        bra_idxo = []
        bra_idxv = []
        ket_idxo = []
        ket_idxv = []
        for elmt in matrix:
            if type(elmt) is state:
                # find rank of bra and ket
                if elmt.ket:
                    ket_rank = elmt.rank
                    ket_idxo += elmt.occ
                    ket_idxv += elmt.vir
                else:
                    bra_rank = elmt.rank
                    bra_idxo += elmt.occ
                    bra_idxv += elmt.vir
            elif type(elmt) is el_oper:
                nel_oper += 1
                #FIXME: assuming first index is virtual and second is occupied!
                ket_idxv.append( elmt.i1 )
                ket_idxo.append( elmt.i2 )
            else:
                raise Exception("this should not happen")

        # check if matrix element is zero
        if bra_rank != ket_rank + nel_oper:
            return 0


        # get deltas from matrix elements
        newterm = term(sign=self.sign, elmts=rest)
        # virtual indices
        for kv, bv in zip(ket_idxv, bra_idxv):
            if kv == bv:
                # nothing to do
                continue
            else:
                if kv in dummy:
                    # all kv must be replaced by bv and kv must be removed from the summation
                    newterm = newterm.delta(kv, bv)
                elif bv in dummy:
                    # all bv must be replaced by kv and bv must be removed from the summation
                    newterm = newterm.delta(bv, kv)
                else:
                    # indices are diferent and none of them is a dummy index
                    # so the whole term is zero
                    return 0

        # occupied indices
        for ko, bo in zip(ket_idxo, bra_idxo):
            if ko == bo:
                # nothing to do
                continue
            else:
                if ko in dummy:
                    # all ko must be replaced by bo and ko must be removed from the summation
                    newterm = newterm.delta(ko, bo)
                elif bo in dummy:
                    # all bo must be replaced by ko and bo must be removed from the summation
                    newterm = newterm.delta(bo, ko)
                else:
                    # indices are diferent and none of them is a dummy index
                    # so the whole term is zero
                    return 0

        #if bra_rank==2:
        #    # we need to make a permutation
        #    # < ^ab_ij | E_ck E_dl |HF> = P^ab_ij delta(abij, ckdl)
        #    a = bra_idxv[0]
        #    b = bra_idxv[1]
        #    i = bra_idxo[0]
        #    j = bra_idxo[1]
        #    return newterm.permute([(a,i), (b,j)])

        if bra_rank > 2:
            raise Exception("cannot do that yet")

        return newterm

    def delta(self, i1, i2):
        """ replace indices i1 by i2 for all elements in the term.
        
        also remove summation over the i1 index"""
        # So this is basically applying a delta(i1, i2) on the term

        newelmts = []
        for elmt in self.elmts:
            if type(elmt) is summation: 
                # remove index from summation
                index = [idx for idx in elmt.index if idx != i1]
                if index:
                    newelmts.append( summation(index) )

            elif type(elmt) is tensor:
                newelmts.append( elmt.replace(i1, i2) )

            else:
                newelmts.append( elmt )

        return term( sign=self.sign, elmts=newelmts )

#------------------------------------------------------------------------------
class commutator(object):
    """ A commutator consist of two elements."""

    def __init__(self,i1,i2):
        allowed_types = [operator, el_oper, commutator, expression, term]
        if type(i1) in allowed_types:
            self.i1 = i1
        else:
            raise Exception("wrong type for first argument")

        if type(i2) in allowed_types:
            self.i2 = i2
        else:
            raise Exception("wrong type for second argument")

    def __add__(self,other):
        """ adding a commutator with something else gives an expression"""
        # make commutator as a term and add the other stuff
        return term( elmts=[self] ) + other

    def __mul__(self,other):
        """ multiplying a commutator with something else gives a term """
        return term( elmts=[self] ) * other

    def output(self):
        """ return printable latex string"""
        return '['+self.i1.output()+','+self.i2.output()+']'

    def simple(self):
        """ is the commutator in a simplified format"""
        # commutators are simplified if the first entry is a simplified
        # commutator or an operator or a single elementary operator and 
        # if the second entry is a single elementary operator!
        simp = False
        if type(self.i2) is el_oper:
            if type(self.i1) is el_oper:
                # type [E,E] => commutator is simplified
                simp = True
            elif type(self.i1) is operator:
                # type [F,E] => commutator is simplified
                simp = True
            elif type(self.i1) is commutator and self.i1.simple():
                # type [[simple],E] => commutator is simplified
                simp = True
        return simp


    def nested(self):
        """ Number of nested commutators if it is simplified """
        nested = None
        if self.simple():
            nested = 1
            part1 = self.i1
            while type(part1) is commutator:
                nested += 1
                part1 = part1.i1
        return nested


    def inner_rank(self):
        """down rank of inner operator in simplified commutators """
        ir = None
        if self.simple():
            part1 = self.i1
            while type(part1) is commutator:
                part1 = part1.i1
            if type(part1) is el_oper:
                # for standard singlet excitation operators the down rank is 0
                ir = 0
            elif type(part1) is operator:
                # down rank given by user
                ir = part1.rank
        return ir


    def simplify(self):
        """ simplifies the commutator which generally returns an expression"""
        if (self.simple()):
            # if commutator is already simplified then return it
            return self

        # simplify first element
        if type(self.i1) in simple_types:
            part1 = self.i1 # already simplified

        elif type(self.i1) is commutator:
            part1 = self.i1.simplify()
            if type(part1) is expression:
                # part1 is now an expression which needs to be distributed in the outer commutator
                # For example [A + B, C] = [A,C] + [B,C]
                # and the new terms need to be simplified
                return commutator(part1, self.i2).distribute()
        else:
            print "test", type(self.i1), self.i1.output()
            print "test", type(self.i2), self.i2.output()
            raise Exception("this should not happen")
     
        # simplify second element
        if type(self.i2) is term:
            # use commutator identity to simplify the expression
            # [A, B1... Bn] = sum_k^n B1...Bk-1 [A,Bk] Bk+1... Bn
            # this return an expression which has to be simplifyed
            return identity1(part1, self.i2).simplify()

        else:
            raise Exception("cannot do that yet")

    def distribute(self):
        """distribute commutator into a simplified expression"""
        # For example [A + B, C] = [A,C] + [B,C]
        # and simplify each term
        if type(self.i1) is not expression:
            raise Exception("I do not think this should happen")

        if type(self.i2) is expression:
            raise Exception("I do not think this should happen")

        # create a new commutator for each term of the expression in self.i1
        newterms = []
        for t in self.i1.terms:

            # check that terms are simplified
            if not t.simple():
                raise Exception("I do not think this should happen")

            # if the term contain more than 1 elements it needs to be decomposed 
            if  len(t.elmts) > 1:
                # for example: [E1 E2, E3] = E1 [E2, E3] + [E1, E3] E2
                newexp = identity2(t, self.i2).simplify()
                newterms += newexp.terms
            else:
                # just build a commutator with the single element of the term
                comm = commutator( t.elmts[0], self.i2).simplify()
                newterms.append( term( elmts=[comm] ) )

        return expression( newterms )


    def vanish(self):
        """ take a single simplified commutator and apply rank reduction"""
        # commutator has the form C = [[[A,E1],E2],E3]
        # C = 0 if the number of nested commutators (here 3) 
        # is larger than twice the down rank of A
        if (self.nested() > 2*self.inner_rank()):
            return True
        else:
            return False
    

    def permute(self, idmap):
        """ permute indices in the commutator"""
        i1 = self.i1.permute(idmap)
        i2 = self.i2.permute(idmap)
        return commutator(i1, i2)


#------------------------------------------------------------------------------
class el_oper(object):
    """ elementary singlet excitation operator E_pq"""

    def __init__(self,i1,i2):
        if type(i1) is not str or type(i2) is not str:
            raise Exception("wrong input for el_oper")
        self.i1 = i1
        self.i2 = i2

    def simple(self):
        # always simple
        return True

    def simplify(self):
        # elementary operators are always simplified 
        return self

    def __add__(self,other):
        """ adding an el. operator to something else gives an expression"""
        return term( elmts=[self] ) + other

    def __mul__(self,other):
        """ an elementary operator time another one gives a term"""
        return term( elmts=[self] ) * other

    def output(self):
        """ return printable latex string"""
        return 'E_{' + self.i1 + self.i2 + '}'

    def permute(self, idmap):
        """ permute indices in the operator """
        newi1, newi2 = permute([self.i1, self.i2], idmap)
        return el_oper(newi1, newi2)


#------------------------------------------------------------------------------
class operator(object):

    def __init__(self,rank=1,name="dummy",symb="D"):
        self.rank = rank
        self.name = name
        self.symb = symb

        # explicit representation: F = sum_pq F_pq E_pq
        index = []
        elmts = []
        for i in range(0,2*rank,2):
            i1 = uns_orb()
            i2 = uns_orb()
            index.append( i1 )
            index.append( i2 )
            elmts.append( el_oper(i1, i2) )

        # declare summation
        self.summ = summation(index)
        # declare tensor F_pq
        self.tensor = tensor(index,self.symb)
        # declare string of elementary operators
        self.opers = term( elmts=elmts )

    def simple(self):
        # always simple
        return True

    def simplify(self):
        # Operators are always simplified 
        return self

    def output(self):
        return self.symb

    def output_exp(self):
        return self.summ.output() + self.tensor.output() + self.opers.output()

    def permute(self, idmap):
        """ nothing to permute here """
        return self


#------------------------------------------------------------------------------
class state(object):
    """ a state can be a bra or a ket vector and has a rank.

    rank = 0 corresponds to HF state
    rank = 1 to a singly excited determinant
    rank = 2 to a doubly excited determinant..."""

    def __init__(self, rank=0, ket=True, vir=[], occ=[]):
        if any(type(i) is not str for i in vir):
            raise Exception("a list of strings is needed here")
        if any(type(i) is not str for i in occ):
            raise Exception("a list of strings is needed here")
        self.rank = rank   # excitation rank, 0 = HF
        self.ket = ket     # ket or bra state?
        self.vir = vir # list of virtual indices for excited determinants
        self.occ = occ # list of occupied indices for excited determinants

    def output(self):
        """ return latex string to print state"""
        if self.rank == 0:
            string = '\\text{HF}'
        else:
            string = '^{'+ ''.join(self.vir) +'}_{'+ ''.join(self.occ) +'}'

        if self.ket:
            return '\ket{'+string+'}'
        else:
            return '\\bar{\\bra{'+string+'}}'

    def permute(self, idmap):
        if self.rank == 0:
            return self
        else:
            index2 = self.vir + self.occ
            newidx = permute(index2, idmap)
            newvir = []
            newocc = []
            rk = self.rank
            for i in range(rk):
                newvir.append(newidx[i])
                newocc.append(newidx[rk + i])

            return state(self.rank, self.ket, newvir, newocc)


#------------------------------------------------------------------------------
class tensor(object):

    def __init__(self, index, symb):
        if any(type(i) is not str for i in index):
            raise Exception("a list of strings is needed here")
        if type(symb) is not str:
            raise Exception("a string is needed here")
        self.index = index # list of indices of tensor
        self.symb = symb   # tensor symbol

    def output(self):
        return self.symb+'_{'+''.join(self.index)+'}'

    def permute(self, idmap):
        newidx = permute(self.index, idmap)
        return tensor(newidx, self.symb)

    def replace(self, i1, i2):
        #FIXME use permute instead!!
        """ replace index i1 by i2 in tensor"""
        newidx = []
        for idx in self.index:
            if idx == i1:
                # replace
                newidx.append( i2 )
            else:
                # copy
                newidx.append( idx )
        return tensor( newidx, self.symb)



#------------------------------------------------------------------------------
class summation(object):

    def __init__(self, index):
        if any(type(idx) is not str for idx in index):
            raise Exception("wrong input for summation")
        self.index = index # list summation indices

    def output(self):
        if not self.index: 
            return ''

        string = '\sum_{'
        for s in self.index:
            string += s
        string += '}'
        return string

    def permute(self, idmap):
        newidx = permute(self.index, idmap)
        return summation( newidx )


# list of elements that are already considered simplified
simple_types = [int, float, summation, tensor, operator, el_oper, state]
matrix_elmts = [state, el_oper, commutator]


#------------------------------------------------------------------------------
def evaluate_comm_HF(comm, hf):
    """ return expression for commutators on HF state"""
    # check validity of input
    if type(comm) is not commutator:
        raise Exception("wrong input")
    if not comm.simple():
        raise Exception("wrong input")
    if type(hf) is not state:
        raise Exception("wrong input")
    if hf.rank != 0 or not hf.ket:
        raise Exception("wrong input")

    if comm.inner_rank() == 1:
        # Type Fock matrix [F,E] |HF>
        return evaluate_rank1(comm, hf)

    elif comm.inner_rank() == 2:
        # Type fluctuation potential matrix [Phi,E] |HF>
        return evaluate_rank2(comm, hf)

    else:
        raise Exception("cannot do that yet")


def evaluate_rank1(comm, hf):
    """ return expression for "rank1" commutators on HF state"""
    # meaning commutator of the form [[F,E],E] |HF>  with F of rank 1
    if comm.inner_rank() != 1:
        raise Exception("wrong input")

    if comm.nested() == 1:
        # [F,E_ai] |HF> = 2 F_ia |HF> + sum_b F_ba E_bi |HF> - sum_j F_ij E_aj |HF>
        a = comm.i2.i1
        i = comm.i2.i2 
        b = vir_orb()
        j = occ_orb()

        # first term: 2 F_ia |HF> 
        tens = tensor([i, a], comm.i1.symb)
        term1 = term( elmts=[2, tens, hf] )

        # second term: sum_b F_ba E_bi |HF>
        tens = tensor([b, a], comm.i1.symb)
        oper = el_oper(b, i)
        term2 = term( elmts=[summation([b]), tens, oper, hf] )

        # third term: - sum_j F_ij E_aj |HF>
        tens = tensor([i, j], comm.i1.symb)
        oper = el_oper(a, j)
        term3 = term( sign='-', elmts=[summation([j]), tens, oper, hf] )
        
        return term1 + term2 + term3

    elif comm.nested() == 2:
        # [[F,E_ai],E_bj] = - P^{ab}_{ij} F_ib E_aj |HF> 
        b = comm.i2.i1
        j = comm.i2.i2
        a = comm.i1.i2.i1
        i = comm.i1.i2.i2

        tens = tensor([i, b], comm.i1.i1.symb)
        oper = el_oper(a, j)
        term1 = term(sign='-', elmts=[tens, oper, hf])

        return term1.permute([(a, i), (b, j)])

    else:
        # this is zero and should have been simplified before
        raise Exception("should not happen")


def evaluate_rank2(comm, hf):
    """ return expression for "rank2" commutators on HF state"""
    # meaning commutator of the form [[Phi,E],E] |HF>  with Phi of rank 2
    if comm.inner_rank() != 2:
        raise Exception("wrong input")

    if comm.nested() == 1:
        # [Phi,E_ai] |HF> = sum_bj L_{bjia} E_{bj} |HF> + ...
        a = comm.i2.i1
        i = comm.i2.i2
        b = vir_orb()
        j = occ_orb()
        c = vir_orb()
        k = occ_orb()

        # first term:  sum_bj L_{bjia} E_{bj} |HF>
        tens = tensor([b, j, i, a], 'L')
        oper = el_oper(b, j)
        term1 = term( elmts=[summation([b, j]), tens, oper, hf] )

        # second term: sum_bjc g_{bjca} E_{bj} E_{ci} |HF> 
        tens = tensor([b, j, c, a], 'g')
        oper1 = el_oper(b, j) 
        oper2 = el_oper(c, i)
        term2 = term( elmts=[summation([b, j, c]), tens, oper1, oper2, hf] )

        # third term: - sum_bjk g_{bjik} E_{bj} E_{ak} |HF> 
        tens = tensor([b, j, i, k], 'g')
        oper1 = el_oper(b, j)
        oper2 = el_oper(a, k)
        term3 = term( sign='-', elmts=[summation([b, j, k]), tens, oper1, oper2, hf] )

        return term1 + term2 + term3

    elif comm.nested() == 2:
        # [[Phi,E_ai],E_bj] |HF> =  ...
        b = comm.i2.i1
        j = comm.i2.i2
        a = comm.i1.i2.i1 
        i = comm.i1.i2.i2 
        c = vir_orb()
        k = occ_orb()
        d = vir_orb()
        l = occ_orb()

        # term 1: 2 L_{iajb} |HF>
        tens = tensor([i, a, j, b], 'L')
        term1 = term( elmts=[2, tens, hf] )

        # term 2: P^{ab}_{ij} ( sum_c L_{cajb} E_ci - sum_k L_{ikjb} E_{ak} ) |HF>
        tens = tensor([c, a, j, b], 'L')
        oper = el_oper(c, i)
        term21 = term( elmts=[summation([c]), tens, oper] )

        tens = tensor([i, k, j, b], 'L')
        oper = el_oper(a, k)
        term22 = term(sign='-', elmts=[summation([k]), tens, oper] )

        exp2 = term21 + term22
        exp2 = exp2.permute([(a, i), (b, j)])

        term2 = exp2 * hf


        # term 3: - sum_ck P^{ab}_{ij} (g_{ibck} E_aj E_ck + g_{ikcb} E_ak E_cj ) |HF>
        tens = tensor([i, b, c, k], 'g')
        oper1 = el_oper(a, j)
        oper2 = el_oper(c, k)
        term31 = term( elmts=[tens, oper1, oper2] )

        tens = tensor([i, k, c, b], 'g')
        oper1 = el_oper(a, k)
        oper2 = el_oper(c, j)
        term32 = term( elmts=[tens, oper1, oper2] )

        exp3 = term31 + term32
        exp3 = exp3.permute([(a, i), (b, j)])

        term3 = term( sign='-', elmts=[summation([c,k]), exp3, hf] )


        # term 4: sum_kl g_{ikjl} E_ak E_bl |HF>
        tens = tensor([i, k, j, l], 'g')
        oper1 = el_oper(a, k)
        oper2 = el_oper(b, l)
        term4 = term( elmts=[summation([k,l]), tens, oper1, oper2, hf] )

        # term 5: sum_cd g_{cadb} E_ci E_dj |HF>
        tens = tensor([c, a, d, b], 'g')
        oper1 = el_oper(c, i)
        oper2 = el_oper(d, j)
        term5 = term( elmts=[summation([c,d]), tens, oper1, oper2, hf] )

        return term1 + term2 + term3 + term4 + term5

    elif comm.nested() == 3:
        # [[[Phi,E_ai],E_bj],E_ck] |HF> =  ...
        c = comm.i2.i1
        k = comm.i2.i2
        b = comm.i1.i2.i1
        j = comm.i1.i2.i2
        a = comm.i1.i1.i2.i1 
        i = comm.i1.i1.i2.i2
        d = vir_orb()
        l = occ_orb()

        # term 1: - L_{jbic} E_ak |HF>
        tens = tensor([j, b, i, c], 'L')
        oper = el_oper(a, k)
        term1 = term(sign='-', elmts=[tens, oper, hf] )

        # term 2: + sum_l g_{iljc} E_al E_bk |HF>
        tens = tensor([i, l, j, c], 'g')
        oper1 = el_oper(a, l)
        oper2 = el_oper(b, k)
        term2 = term(elmts=[summation([l]), tens, oper1, oper2, hf] )

        # term 3: - sum_d g_{ibdc} E_aj E_dk |HF>
        tens = tensor([i, b, d, c], 'g')
        oper1 = el_oper(a, j)
        oper2 = el_oper(d, k)
        term3 = term(sign='-', elmts=[summation([d]), tens, oper1, oper2, hf] )

        # final expression with permutation
        exp = term1 + term2 + term3
        exp = exp.permute([(a, i), (b, j), (c, k)])

        print(exp.output())
        return exp

    elif comm.nested() == 4:
        raise Exception("not yet implemented")

    else:
        # this is zero and should have been simplified before
        raise Exception("should not happen")


#------------------------------------------------------------------------------
def identity1(A, B):
    # assuming B is of type term
    # [A, B1... Bn] = sum_k^n B1...Bk-1 [A,Bk] Bk+1... Bn
    newterms = []
    s = B.sign
    for i, elmt in enumerate(B.elmts):
        if i==0:
            newterm = commutator(A, elmt) * term(sign=s, elmts=B.elmts[i+1:] )
        elif i==len(B.elmts)-1:
            newterm =  term(sign=s, elmts=B.elmts[0:i] ) * commutator(A, elmt)
        else:
            newterm =  term(sign=s, elmts=B.elmts[0:i] ) * commutator(A, elmt) * term(sign=s, elmts=B.elmts[i+1:] )
        newterms.append( newterm )
    return expression( newterms )


#------------------------------------------------------------------------------
def identity2(A, B):
    # assuming A is of type term
    # [A1... An, B] = sum_k^n A1...Ak-1 [Ak,B] Ak+1... An
    newterms = []
    s = A.sign
    for i, elmt in enumerate(A.elmts):
        if i==0:
            newterm = commutator(elmt, B) * term(sign=s, elmts=A.elmts[i+1:] )
        elif i==len(A.elmts)-1:
            newterm =  term(sign=s, elmts=A.elmts[0:i] ) * commutator(elmt, B)
        else:
            newterm =  term(sign=s, elmts=A.elmts[0:i] ) * commutator(elmt, B) * term(sign=s, elmts=A.elmts[i+1:] )
        newterms.append( newterm )
    return expression( newterms )


#------------------------------------------------------------------------------
def identity3(A, B):
    # assuming A is a commutator and B is a term
    # for example we can have: [F,E] E = [[F,E],E] + E [F, E]
    if len(B.elmts) > 1:
        # solve recursively
        term1 = term(sign=B.sign, elmts=[B.elmts[1:]] )
        exp = identity3( commutator(A, B.elmts[0] ), term1 )
        term2 = B.elmts[0] * identity3(A, term1 )
        return exp + term2
    else:
        term1 = term( sign=B.sign, elmts=[commutator(A, B.elmts[0])] )
        term2 = B.elmts[0] * A
        term2.sign = B.sign
        return term1 + term2

            
#------------------------------------------------------------------------------
def new_exc_oper():
    """ return single excitation operator of type E_ai"""
    i1 = vir_orb() 
    i2 = occ_orb() 
    return el_oper(i1, i2)


#------------------------------------------------------------------------------
def bra(rank=0, vir=[], occ=[]):
    """ return a bra state of a given excitation rank"""
    return state(rank=rank, ket=False, vir=vir, occ=occ )

def ket(rank=0, vir=[], occ=[]):
    """ return a ket state of a given excitation rank"""
    return state(rank=rank, ket=True, vir=vir, occ=occ )

#------------------------------------------------------------------------------
def permute(idlist, idmap):
    """replace indices in idlist based on idmap"""
    newlist = []
    for idx in idlist:

        if idx in idmap:
            newlist.append( idmap[idx] )
        else:
            newlist.append( idx )

    return newlist


