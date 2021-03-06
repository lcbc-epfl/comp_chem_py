#!/bin/usr/env python

import sys
import copy

iuns_orb = 0
def uns_orb():
    """ return new index for unspecified orbital and increase counter """
    global iuns_orb
    iuns_orb += 1
    return 'u_{' + str( iuns_orb ) + '}'

iocc_orb = 0
def occ_orb():
    """ return new index for occupied orbital and increase counter """
    global iocc_orb
    iocc_orb += 1
    return 'o_{' + str( iocc_orb ) + '}'

ivir_orb = 0
def vir_orb():
    """ return new index for virtual orbital and increase counter """
    global ivir_orb
    ivir_orb += 1
    return 'v_{' + str( ivir_orb ) + '}'


def new_el_oper():
    """Return a single excitation operator of type E_pq"""
    return el_oper(uns_orb(), uns_orb())

def new_exc_oper():
    """Return a single excitation operator of type E_ai"""
    return el_oper(vir_orb(), occ_orb(), rank=0)


def tensor(rank, fac=1.0):
    """Create a term containing a summation and a string of elementary operators."""
    indices = []
    newterm = term([], fac=fac)
    for i in range(rank):
        newel = new_el_oper()
        indices.append(newel.i1)
        indices.append(newel.i2)
        newterm *= newel

    return summation( indices ) * newterm


class el_oper(object):
    """Elementary singlet excitation operator E_pq."""

    def __init__(self,i1,i2,rank=1):
        assert isinstance(i1, str), "Wrong input for first argument!"
        assert isinstance(i2, str), "Wrong input for second argument!"
        self.i1 = i1
        self.i2 = i2
        self.rank=rank

    def simple(self):
        # always simple
        return True

    def simplify(self):
        # elementary operators are always simplified
        return self

    def __eq__(self,other):
        if isinstance(other, type(self)) \
            and self.i1 == other.i1 \
            and self.i2 == other.i2 \
            and self.rank == other.rank:
            return True
        else:
            return False

    def __add__(self,other):
        """Adding an el. operator to something else gives an expression."""
        return term( [self] ) + other

    def __sub__(self,other):
        """Substracting an el. operator to something else gives an expression."""
        return term( [self] ) - other

    def __mul__(self,other):
        """An elementary operator times something else gives a term."""
        return term( [self] ) * other

    def output(self):
        """Return printable latex string."""
        return 'E_{' + self.i1 + self.i2 + '}'

    #def permute(self, idmap):
    #    """Permute indices in the operator."""
    #    newi1, newi2 = permute([self.i1, self.i2], idmap)
    #    return el_oper(newi1, newi2, rank=self.rank)


class delta(object):
    """Simple kronecker delta operator."""

    def __init__(self,i1,i2):
        assert isinstance(i1, str), "Wrong input for first argument!"
        assert isinstance(i2, str), "Wrong input for second argument!"
        self.i1 = i1
        self.i2 = i2

    def simple(self):
        # always simple
        return True

    def simplify(self):
        # elementary operators are always simplified
        return self

    def __eq__(self,other):
        if isinstance(other, type(self)) and self.i1 == other.i1 and self.i2 == other.i2:
            return True
        else:
            return False

    def __add__(self,other):
        """Adding a detla to something else gives an expression."""
        return term( [self] ) + other

    def __sub__(self,other):
        """Substracting a delta to something else gives an expression."""
        return term( [self] ) - other

    def __mul__(self,other):
        """A delta operator times something else gives a term."""
        return term( [self] ) * other

    def output(self):
        """Return printable latex string."""
        return '\\delta_{' + self.i1 + self.i2 + '}'

    #def permute(self, idmap):
    #    """Permute indices in the operator."""
    #    newi1, newi2 = permute([self.i1, self.i2], idmap)
    #    return delta(newi1, newi2)


class summation(object):

    def __init__(self, index):
        assert all( isinstance(idx, str) for idx in index)
        self.index = index # list summation indices

    def simple(self):
        # always simple
        return True

    def simplify(self):
        # elementary operators are always simplified
        return self

    def output(self):
        if not self.index:
            return ''

        string = '\sum_{'
        for s in self.index:
            string += s
        string += '}'
        return string

    def __mul__(self, other):
        """A summation times something else gives a term."""
        return term( [self] ) * other

    #def permute(self, idmap):
    #    newidx = permute(self.index, idmap)
    #    return summation( newidx )


class commutator(object):
    """A commutator consist of two elements."""

    def __init__(self,i1,i2):
        allowed_types = [el_oper, commutator, term, expression]
        assert type(i1) in allowed_types, "Wrong type for first argument!"
        assert type(i2) in allowed_types, "Wrong type for second argument!"
        self.i1 = i1
        self.i2 = i2

    def __eq__(self,other):
        if isinstance(other, type(self)) and self.i1 == other.i1 and self.i2 == other.i2:
            return True
        else:
            return False

    def __add__(self,other):
        """Adding a commutator with something else gives an expression."""
        # make commutator as a term and add the other stuff
        return term( [self] ) + other

    def __sub__(self,other):
        """Substracting a commutator to something else gives an expression."""
        return term( [self] ) - other

    def __mul__(self,other):
        """Multiplying a commutator with something else gives a term."""
        return term( [self] ) * other

    def output(self):
        """Return printable latex string."""
        return '['+self.i1.output()+','+self.i2.output()+']'

    def simple(self):
        """Check if the commutator is in a simplified format."""
        # commutators are simplified if the first entry is a simplified
        # commutator or a single elementary operator and
        # if the second entry is a single elementary operator!
        simp = False
        if type(self.i2) is el_oper:
            if type(self.i1) is el_oper:
                # type [E,E] => commutator is simplified
                simp = True
            elif type(self.i1) is commutator and self.i1.simple():
                # type [[simple],E] => commutator is simplified
                simp = True
        return simp


    def nested(self):
        """Number of nested commutators if it is simplified."""
        nested = None
        if self.simple():
            nested = 1
            part1 = self.i1
            while type(part1) is commutator:
                nested += 1
                part1 = part1.i1
        return nested


    def inner_rank(self):
        """Down rank of inner operator in simplified commutators."""
        ir = None
        if self.simple():
            part1 = self.i1
            while type(part1) is commutator:
                part1 = part1.i1
            if type(part1) is el_oper:
                # for standard singlet excitation operators the down rank is 0
                ir = part1.rank
        return ir


    def reduce(self):
        """Apply [E_mn,E_pq] = E_mq delta_pn - E_pn delta_mq."""

        # check input
        assert self.simple()
        if self.nested() == 1:

            m = self.i1.i1
            n = self.i1.i2
            p = self.i2.i1
            q = self.i2.i2
            return el_oper(m,q) * delta(p,n) - el_oper(p,n) * delta(m,q)

        else:
            # reduce inner commutator and simplify the outer one
            newexp = commutator(self.i1.reduce(), self.i2).simplify()
            # return reduced expression
            return newexp.reduce()


    def identity1(self):
        """[A, B1... Bn] = sum_k^n B1...Bk-1 [A,Bk] Bk+1... Bn"""
        A = self.i1
        B = self.i2
        # simple checks
        assert isinstance(B, term)

        # isolate factor and delta from term
        fac_term, op_term = B.extract_delta()

        if len(op_term.elmts) == 0:
            # The commutator is zero!
            raise Exception("This should not happen!")
        elif len(op_term.elmts) == 1:
            # nothing to do
            return fac_term * commutator(A,op_term.elmts[0])

        newexp = expression( [] )
        for i, elmt in enumerate(op_term.elmts):
            if i==0:
                newexp += fac_term * commutator(A, elmt) * term( op_term.elmts[i+1:] )

            elif i==len(op_term.elmts)-1:
                newexp += fac_term * term( op_term.elmts[0:i] ) * commutator(A, elmt)

            else:
                newexp += fac_term * term( op_term.elmts[0:i] ) * commutator(A, elmt) * term( op_term.elmts[i+1:] )

        return newexp


    def identity2(self):
        """[A1... An, B] = sum_k^n A1...Ak-1 [Ak,B] Ak+1... An"""
        A = self.i1
        B = self.i2
        # simple checks
        assert isinstance(A, term)

        # isolate factor and delta from term
        fac_term, op_term = A.extract_delta()

        if len(op_term.elmts) == 0:
            # The commutator is zero!
            raise Exception("This should not happen!")
        elif len(op_term.elmts) == 1:
            # nothing to do
            return fac_term * commutator(op_term.elmts[0],B)

        newexp = expression( [] )
        for i, elmt in enumerate(op_term.elmts):
            if i==0:
                newexp += fac_term * commutator(elmt, B) * term( op_term.elmts[i+1:] )

            elif i==len(op_term.elmts)-1:
                newexp +=  fac_term * term( op_term.elmts[0:i] ) * commutator(elmt, B)

            else:
                newexp +=  fac_term * term( op_term.elmts[0:i] ) * commutator(elmt, B) * term( op_term.elmts[i+1:] )

        return newexp


    def identity3(self, B):
        """ Apply: [A] B = [[A],B] + B [A].

        Maybe recursively if B is a chain of el_oper:
            [A] B1 B2 = [[A], B1] B2 + B1 [A] B2
                      = [[[A], B1], B2] + B2 [[A], B1] + B1 ( [[A],B2] + B2 [A] )
        """

        A = self
        assert A.simple(),"identity3: commutator should be simplified!"
        assert isinstance(B, term), "B should always be a term"
        # we assume the term contains only el_oper
        assert all(isinstance(elmt, el_oper) for elmt in B.elmts)

        if len(B.elmts)==0:
            # B is empty, we just return A scale with B.factor
            print "I don't think this should happen..."
            return A * B.fac

        # Get first el_oper in B term
        el_op = B.elmts[0]
        if len(B.elmts)==1:
            # B contains only one el_oper
            # Back to simple case (conserve factor!)
            exp = commutator(A, el_op) * B.fac  + el_op * A * B.fac
            return exp

        if len(B.elmts)>1:
            # B contains more than one el_oper
            # Solve recursively
            newB = term( B.elmts[1:], fac=B.fac)
            exp1 = commutator(A, el_op).identity3( newB )
            exp2 = el_op * A.identity3( newB )
            return exp1 + exp2


    def distribute1(self):
        """Distribute expression inside a commutator.

        [A1+A2+...+An, B] = [A1,B] + [A2,B] +... [An,B]
        """
        A = self.i1
        B = self.i2
        assert isinstance(A, expression)
        # create a new commutator for each term of the expression in A
        newexp = expression( [] )
        for Ai in A.terms:
            newexp += commutator( Ai, B )

        return newexp


    def distribute2(self):
        """Distribute expression inside a commutator.

        [A, B1+B2+...+Bn] = [A,B1] + [A,B2] +... [A,Bn]
        """
        A = self.i1
        B = self.i2
        assert isinstance(B, expression)
        # create a new commutator for each term of the expression in B
        newexp = expression( [] )
        for Bi in B.terms:
            newexp += commutator( A, Bi )

        return newexp


    def vanish(self):
        """ take a single simplified commutator and apply rank reduction"""
        # commutator has the form C = [[[A,E1],E2],E3]
        # C = 0 if the number of nested commutators (here 3)
        # is larger than twice the down rank of A
        if (self.nested() > 2*self.inner_rank()):
            return True
        else:
            return False


    def simplify(self):
        """Simplifies the commutator which generally returns an expression."""

        # check if commutator is already simplified
        if self.simple():
            return self

        # simplify each element of the commutator
        part1 = self.i1.simplify()
        part2 = self.i2.simplify()

        newcom = commutator(part1, part2)

        # check if new commutator is simplified
        if newcom.simple():
            print "Should not happen"
            return newcom

        if isinstance(part1, commutator) or isinstance(part1, el_oper):
            # part1 is a simplified commutator or just an el_oper

            if isinstance(part2, commutator):
                # expand part2: [part1,[E1,E2]] = [part1,E1 E2] - [part1, E2 E1]
                com1 = commutator(part1, part2.i1 * part2.i2)
                com2 = commutator(part1, part2.i2 * part2.i1)
                exp = com1 + term( [com2], fac=-1.0)

                # return simplified new expression
                return exp.simplify()

            elif isinstance(part2, term):
                # use identity1 to simplify commutator
                exp = newcom.identity1()

                # return simplified new expression
                return exp.simplify()

            elif isinstance(part2, expression):
                # distribute expression in the commutator and simplify
                exp = newcom.distribute2()

                # return simplified new expression
                return exp.simplify()
            else:
                raise Exception("This should not happen!")

        elif isinstance(part1, term):
            # use commutator identity2
            exp = newcom.identity2()

            # return simplified new expression
            return exp.simplify()

        elif isinstance(part1, expression):
            # distribute expression and simplify
            exp = newcom.distribute1()

            # return simplified new expression
            return exp.simplify()

        else:
            raise Exception("This should not happen!")


    #def permute(self, idmap):
    #    """ permute indices in the commutator"""
    #    i1 = self.i1.permute(idmap)
    #    i2 = self.i2.permute(idmap)
    #    return commutator(i1, i2)


class expression(object):
    """An expression is defined as a sum of terms."""

    def __init__(self, terms):
        # Make clean list of terms such that it contains only allowed types
        newterms = []
        for t in terms:
            if isinstance(t, expression):
                # an input element is an expression
                # so we add the terms of the expression instead
                newterms += t.terms
            elif isinstance(t, term):
                newterms.append( t )
            else:
                raise Exception("Wrong element in expression!")

        self.terms = newterms

    def simple(self):
        # expression is simplified if all the terms are simplified
        for t in self.terms:
            if not t.simple():
                return False
        return True

    def simplify(self):

        # check if not already simple!
        if self.simple():
            return self

        # simplify each term and return new expression
        newexp = expression([])
        for t in self.terms:
            newexp += t.simplify()

        return newexp

    def reduce(self):

        # reduce each term and return new expression
        newexp = expression([])
        for t in self.terms:
            newexp += t.reduce()

        return newexp

    def __eq__(self,other):
        """Two expressions are equals if each term is present in the other one."""
        if isinstance(other, type(self)) \
           and len(self.terms) == len(other.terms):
            return all(t in other.terms for t in self.terms)
        else:
            return False

    def __mul__(self,other):
        """An expression times something else gives a term."""
        return term( [self] ) * other

    def __sub__(self,other):
        """See __isub__."""
        newexp = copy.deepcopy(self)
        newexp -= other
        return newexp

    def __isub__(self,other):
        """An expression minus something else gives an expression."""
        newterms = self.terms
        if isinstance(other, expression):
            # change the sign of all the terms inside the expression
            tmp = []
            for t in other.terms:
                tmp.append( t * (-1.0) )
            newterms += tmp

        elif isinstance(other, term):
            # change sign of the term
            newterm = other * (-1.0)
            newterms.append( newterm )

        else:
            # try making a term out of the object and substract it
            newterms.append( term( [other] ) * (-1.0))

        return expression( newterms )

    def __add__(self,other):
        """See __iadd__."""
        newexp = copy.deepcopy(self)
        newexp += other
        return newexp

    def __iadd__(self,other):
        """An expression plus something else gives an expression."""
        newterms = self.terms
        if isinstance(other, expression):
            newterms += other.terms

        elif isinstance(other, term):
            newterms.append( other )

        else:
            # try making a term out of the object and add it
            newterms.append( term( [other] ) )

        return expression( newterms )

    def output(self):
        """ return printable latex string corresponding to the expression"""
        return ' '.join( [t.output() for t in self.terms] )

    def rm_parenthesis(self):
        """ ditribute parenthesis present in individual terms"""
        newexp = expression([]) # empty expression
        for t in self.terms:
            newexp += t.rm_parenthesis()
        return newexp

    #def permute(self, index):
    #    """ perform permutations for the terms in the expression """
    #    # each term returns a sum of terms, i.e. an expression
    #    newexp = expression([]) # empty expression
    #    for t in self.terms:
    #        newexp += t.permute( index )
    #    return newexp

    #def evaluate(self):
    #    """ use boxes from the bible and logic to evaluate the individual terms"""
    #    # we are assuming a simplified expression
    #    if not self.simple():
    #        raise Exception("this should not happen")

    #    newexp = expression([]) # empty expression
    #    for t in self.terms:
    #        # check if the term vanishes due to rank reduction
    #        if not t.vanish():
    #            # evaluate terms that do not vanish
    #            # i.e. get expression for the commutator on HF
    #            newexp += t.evaluate()

    #    return newexp


class term(object):
    """A term is a list of elmts multiplied together with factor in front."""

    def __init__(self, elmts, fac=1.0):
        # the allowed elements can be of the following types:
        allowed_types = [summation, delta, el_oper, commutator, expression]

        # Make clean list of elements such that it contains only allowed types
        newelmts = []
        for elmt in elmts:
            if isinstance(elmt, term):
                # an input element is a term
                # so we add the element of the term instead
                newelmts += elmt.elmts
            elif type(elmt) in allowed_types:
                newelmts.append( elmt )
            else:
                raise Exception("Wrong element in term!")

        self.elmts = newelmts
        self.fac = float(fac)


    def simple(self):
        """A term is simplified if:

            1) It does not contain expressions
            2) It contains only simplified commutators and elementary operators
            3) commutators are not followed by elementary operatots

        Simple terms:
            1) E E E
            2) E [E,E]
            3) - 3.0 E E [[E,E],E]

        Non simplified terms:
            1) Contain expression: E [E,E](E+[E,E])
            2) Commutator is not simplified: E [(E+E),E]
            1) Commutator followed by el_oper: [E,E] E
        """
        for i, elmt in enumerate(self.elmts):

            if isinstance(elmt, expression):
                return False

            elif isinstance(elmt, commutator):
                if not elmt.simple():
                    # Commutator not simplified
                    return False

                if i + 1 < len(self.elmts):
                    if not isinstance(self.elmts[i+1], commutator):
                        # Commutator is not the last element
                        # and it is not followed by another commutator
                        return False

        return True


    def simplify(self):

        # check if not already simple!
        if self.simple():
            return self

        # check if it contains expressions
        if any(isinstance(elmt, expression) for elmt in self.elmts):
            # distribute expresions in the term
            exp = self.rm_parenthesis()
            return exp.simplify()

        # check if all elements are simplified
        if not all(elmt.simple() for elmt in self.elmts):
            # if some elements are not simple
            # we simplify them and create a
            # new simplified term

            newterm = term( [], fac=self.fac)
            for elmt in self.elmts:
                newterm *= elmt.simplify()

            return newterm.simplify()


        fac_term, op_term = self.extract_delta()
        # At this stage we should have only terms that contain
        # simplified commutators or elementary operators!
        #
        # First we get a set of index positions for the commutators
        # and check that the term only contains commutators and el_oper
        cpos = []
        for i, elmt in enumerate(op_term.elmts):
            if isinstance(elmt, commutator):
                cpos.append(i)
            elif isinstance(elmt, el_oper):
                continue
            else:
                # The term should only contain commutators and el_oper
                raise Exception("This should not happen!")

        if len(cpos) > 1:
            print "Warning: More than one commutator in the term!"


        # Then we make a list of subterms separated by commutators
        # e.g.:
        #   term1 is E E E
        #   term2 is [E,E] E E
        #   term3 is [E,E]
        #
        # if the first element is a commutator then the first subterm is empty
        subterms = [op_term.elmts[0:cpos[0]]]
        #
        for i, ic in enumerate(cpos):

            if i + 1 >= len(cpos):
                # last commutator
                subterms.append(op_term.elmts[ic:])

            else:
                subterms.append(op_term.elmts[ic:cpos[i+1]])

        # If a subterm contains more than one elements it means its
        # a commutator followed by one ore more elementary operators
        # so we treat it with the identity3 routine.
        for sub in subterms:
            if len(sub) == 0:
                # nothing to do
                continue
            elif len(sub) == 1:
                # just copy the single element to the newter
                fac_term *= sub[0]
            else:
                # identity3 return an expression
                fac_term *= sub[0].identity3( term( sub[1:] ) )

        # now we should be back with a term containing expressions
        # which needs to be distributed...
        return fac_term.simplify()


    def reduce(self):
        assert self.simple()
        newterm = term( [], fac=self.fac)
        for elmt in self.elmts:
            if isinstance(elmt, commutator):
                newterm *= elmt.reduce()
            else:
                newterm *= elmt

        return newterm.simplify()

    def rm_parenthesis(self):
        """distribute whatever is inside parenthesis"""
        # the term contains parenthesis if it contains at least one expression
        if not any(isinstance(elmt, expression) for elmt in self.elmts):
            return expression( [self] )

        newexp = expression([])
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
                    newexp += t * term( self.elmts[i+1:] ) * self.fac

                elif i == (len(self.elmts) - 1):
                    # last element
                    newexp += term( self.elmts[:i]) * t * self.fac

                else:
                    # middle element
                    newexp += term( self.elmts[:i]) * t * term( self.elmts[i+1:] ) * self.fac

        return newexp


    def output(self):
        """ return printable string for term """
        fac_term, op_term = self.extract_delta()
        # take care of sign
        if fac_term.fac > 0:
            string = ['+']
        else:
            string = ['-']

        # take care of factor
        if abs(fac_term.fac) != 1.0:
            string.append(str(abs(fac_term.fac)))

        # put delta in front
        for elmt in fac_term.elmts:
            if isinstance(elmt, expression):
                # we need to add parenthesis
                string.append( '( '+ elmt.output() +' )' )
            else:
                string.append( elmt.output() )

        # add other elmts
        for elmt in op_term.elmts:
            if isinstance(elmt, expression):
                # we need to add parenthesis
                string.append( '( '+ elmt.output() +' )' )
            else:
                string.append( elmt.output() )

        return ' '.join(string)


    def extract_delta(self):
        """Split a term into two: one containing kronecker delta,
        prefactor and summations and one with the rest."""

        fac_term = term([], fac=self.fac)
        op_term = term([])
        for elmt in self.elmts:
            if isinstance(elmt, delta):
                fac_term *= elmt
            elif isinstance(elmt, summation):
                fac_term *= elmt
            else:
                op_term *= elmt

        return fac_term, op_term


    def __eq__(self,other):
        """Two terms are equal if they contains the same elements in the
        same order and if the factors are the same."""
        if isinstance(other, type(self)):
            return (self.elmts == other.elmts) and (self.fac == other.fac)
        else:
            return False

    def __add__(self, other):
        """A term plus something else gives an expression."""
        return expression( [self] ) + other

    def __sub__(self, other):
        """A term minus something else gives an expression."""
        return expression( [self] ) - other

    def __mul__(self, other):
        """See __imul__"""
        newterm = copy.deepcopy(self)
        newterm *= other
        return newterm

    def __imul__(self, other):
        """Multiplying a term with something else gives a term."""
        if isinstance(other, float):
            # just modify factor of the current term
            fac = self.fac * other
            elmts = self.elmts

        elif isinstance(other, term):
            # get factor
            fac = self.fac * other.fac
            # combine elements
            elmts = self.elmts + other.elmts

        else:
            # factor cannot be modified
            fac = self.fac
            elmts = self.elmts
            elmts.append(other)

        return term(elmts, fac=fac)


    #def permute(self, index):
    #    """ perform permutations of the indices  in the term"""
    #    # index is a list of tuples, e.g. index = [(a, i), (b, j)]
    #    # in that case we build a dictionary idmap that maps
    #    # index a -> b, i -> j, b -> a, j -> i
    #    # and use it to build a new expression

    #    dim = len(index)
    #    # build map of indices as a dictionary
    #    exp = expression([])

    #    # loop over all possible permutations
    #    for perm in itertools.permutations(index, dim):
    #        idmap = {}
    #        for idx,p in zip(index, perm):
    #            idmap[idx[0]] = p[0]
    #            idmap[idx[1]] = p[1]

    #        # permute indices of each element
    #        new_elmts = []
    #        for elmt in self.elmts:
    #            if type(elmt) is float or type(elmt) is int:
    #                # nothing to do here
    #                new_elmts.append( elmt )
    #            else:
    #                new_elmts.append( elmt.permute(idmap) )

    #        exp += term(sign=self.sign, elmts=new_elmts)

    #    return exp


    #def vanish(self):
    #    """ use rank reduction to remove vanishing commutators"""
    #    # we are assuming a simplified term
    #    if not self.simple():
    #        raise Exception("this should not happen")

    #    # if one commutator is zero the whole term is zero
    #    vanish = False
    #    for elmt in self.elmts:
    #        if type(elmt) is commutator:
    #            vanish = elmt.vanish()
    #            if vanish:
    #                break

    #    return vanish

    #def evaluate(self):
    #    """ use the box from the bible to evaluate a term"""
    #    # the term to evaluate should consist of a section of the form:
    #    #  <BRA| el_oper * el_oper ...* [[operator, el_oper], el_oper] |HF>
    #    # with potential tensors and summations on the left side
    #    if not self.simple():
    #        raise Exception("this should not happen")

    #    # find the commutator evaluate it and return new term
    #    count = 0
    #    newelmts = []
    #    for i, elmt in enumerate(self.elmts):

    #        if type(elmt) is commutator:
    #            # this should happen only once if the term has the correct form
    #            count += 1
    #            if count > 1:
    #                raise Exception("this should not happen")

    #            # evaluate the sequence [] |HF>
    #            newelmts.append( evaluate_comm_HF(elmt, self.elmts[i+1]) )
    #            break

    #        newelmts.append( elmt )

    #    # simplify the term to get expression
    #    newexp = term(sign=self.sign, elmts=newelmts).simplify()

    #    # finally contract the matrix elements using delta functions
    #    newexp2 = expression([])
    #    for t in newexp.terms:
    #        newterm = t.contract()
    #        if newterm != 0:
    #            newexp2 += t.contract()

    #    return newexp2

    #def contract(self):
    #    # check that matrix element is not zero

    #    # first get dummy indices matrix element and the rest
    #    dummy = []
    #    matrix = []
    #    rest = []
    #    for elmt in self.elmts:
    #        if type(elmt) is summation:
    #            dummy += elmt.index
    #            rest.append( elmt )
    #        elif type(elmt) in matrix_elmts:
    #            matrix.append( elmt )
    #        else:
    #            rest.append( elmt )


    #    # get delta from matrix element
    #    nel_oper = 0
    #    bra_idxo = []
    #    bra_idxv = []
    #    ket_idxo = []
    #    ket_idxv = []
    #    for elmt in matrix:
    #        if type(elmt) is state:
    #            # find rank of bra and ket
    #            if elmt.ket:
    #                ket_rank = elmt.rank
    #                ket_idxo += elmt.occ
    #                ket_idxv += elmt.vir
    #            else:
    #                bra_rank = elmt.rank
    #                bra_idxo += elmt.occ
    #                bra_idxv += elmt.vir
    #        elif type(elmt) is el_oper:
    #            nel_oper += 1
    #            #FIXME: assuming first index is virtual and second is occupied!
    #            ket_idxv.append( elmt.i1 )
    #            ket_idxo.append( elmt.i2 )
    #        else:
    #            raise Exception("this should not happen")

    #    # check if matrix element is zero
    #    if bra_rank != ket_rank + nel_oper:
    #        return 0


    #    # get deltas from matrix elements
    #    newterm = term(sign=self.sign, elmts=rest)
    #    # virtual indices
    #    for kv, bv in zip(ket_idxv, bra_idxv):
    #        if kv == bv:
    #            # nothing to do
    #            continue
    #        else:
    #            if kv in dummy:
    #                # all kv must be replaced by bv and kv must be removed from the summation
    #                newterm = newterm.delta(kv, bv)
    #            elif bv in dummy:
    #                # all bv must be replaced by kv and bv must be removed from the summation
    #                newterm = newterm.delta(bv, kv)
    #            else:
    #                # indices are diferent and none of them is a dummy index
    #                # so the whole term is zero
    #                return 0

    #    # occupied indices
    #    for ko, bo in zip(ket_idxo, bra_idxo):
    #        if ko == bo:
    #            # nothing to do
    #            continue
    #        else:
    #            if ko in dummy:
    #                # all ko must be replaced by bo and ko must be removed from the summation
    #                newterm = newterm.delta(ko, bo)
    #            elif bo in dummy:
    #                # all bo must be replaced by ko and bo must be removed from the summation
    #                newterm = newterm.delta(bo, ko)
    #            else:
    #                # indices are diferent and none of them is a dummy index
    #                # so the whole term is zero
    #                return 0

    #    #if bra_rank==2:
    #    #    # we need to make a permutation
    #    #    # < ^ab_ij | E_ck E_dl |HF> = P^ab_ij delta(abij, ckdl)
    #    #    a = bra_idxv[0]
    #    #    b = bra_idxv[1]
    #    #    i = bra_idxo[0]
    #    #    j = bra_idxo[1]
    #    #    return newterm.permute([(a,i), (b,j)])

    #    if bra_rank > 2:
    #        raise Exception("cannot do that yet")

    #    return newterm

    #def delta(self, i1, i2):
    #    """ replace indices i1 by i2 for all elements in the term.

    #    also remove summation over the i1 index"""
    #    # So this is basically applying a delta(i1, i2) on the term

    #    newelmts = []
    #    for elmt in self.elmts:
    #        if type(elmt) is summation:
    #            # remove index from summation
    #            index = [idx for idx in elmt.index if idx != i1]
    #            if index:
    #                newelmts.append( summation(index) )

    #        elif type(elmt) is tensor:
    #            newelmts.append( elmt.replace(i1, i2) )

    #        else:
    #            newelmts.append( elmt )

    #    return term( sign=self.sign, elmts=newelmts )


