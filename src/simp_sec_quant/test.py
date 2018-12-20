#!/usr/bin/env python

from new import *
import unittest

class test_commutator(unittest.TestCase):

    def setUp(self):
        # Build input
        self.A = new_el_oper()
        self.B1 = new_el_oper()
        self.B2 = new_el_oper()
        self.B3 = new_el_oper()
        self.B = self.B1 * self.B2 * self.B3
        self.comA = commutator(self.A,self.A)

    def test_identity1(self):
        # Build correct answer
        term1 = commutator(self.A,self.B1) * self.B2 * self.B3
        term2 = self.B1 * commutator(self.A,self.B2) * self.B3
        term3 = self.B1 * self.B2 * commutator(self.A,self.B3)
        ref = term1 + term2 + term3

        # get methods response
        exp = commutator(self.A,self.B).identity1()

        # check
        self.assertEqual(exp, ref)

        # check that second argument must be a term
        with self.assertRaises(AssertionError):
            self.comA.identity1()


    def test_identity2(self):

        # Build correct answer
        term1 = commutator(self.B1,self.A) * self.B2 * self.B3
        term2 = self.B1 * commutator(self.B2,self.A) * self.B3
        term3 = self.B1 * self.B2 * commutator(self.B3,self.A)
        ref = term1 + term2 + term3

        # get methods response
        exp = commutator(self.B,self.A).identity2()

        # check
        self.assertEqual(exp, ref)

        # check that first argument must be a term
        with self.assertRaises(AssertionError):
            self.comA.identity2()


    def test_identity3(self):

        # Build correct answer
        term1 = commutator( commutator(self.comA, self.B1), self.B2)
        term2 = self.B2 * commutator(self.comA, self.B1)
        term3 = self.B1 * ( commutator(self.comA, self.B2) + self.B2 * self.comA)
        ref = term1 + term2 + term3

        # get methods response
        exp = self.comA.identity3(self.B1 * self.B2)

        # check
        self.assertEqual(exp, ref)

        with self.assertRaises(AssertionError):
            # check that commutator must be simplified
            commutator(self.A,self.B).identity3(self.B)
            # check that argument must be a term
            self.comA.identity3(self.A)
            # check that argument cannot contain commutator
            self.comA.identity3(self.A * self.comA)

    def test_simplify(self):
        # verify a few very simple commutator identities
        exp = commutator(self.comA, self.B1 * self.B2).simplify()

        term1 = commutator(commutator(self.comA, self.B1), self.B2)
        term2 = self.B1 * commutator(self.comA, self.B2)
        term3 = self.B2 * commutator(self.comA, self.B1)
        ref = term1 + term2 + term3

        # check
        self.assertEqual(exp, ref)

        # new expression
        exp = commutator(commutator(self.A,self.B1 * self.B2), self.B3).simplify()

        term1 = commutator(commutator(commutator(self.A,self.B1),self.B2),self.B3)
        term2 = self.B1 * commutator(commutator(self.A,self.B2),self.B3)
        term3 = commutator(self.B1,self.B3) * commutator(self.A, self.B2)
        term4 = self.B2 * commutator(commutator(self.A,self.B1),self.B3)
        term5 = commutator(self.B2,self.B3) * commutator(self.A, self.B1)

        ref = term1 + term2 + term3 + term4 + term5

        # check
        self.assertEqual(exp, ref)


if __name__=="__main__":
    unittest.main()
