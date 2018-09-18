Partial Density Of States
-------------------------

The DOS is the function,

.. math::
   DOS(E) = \sum_i g(E-\epsilon_i)

where :math:`g` is a gaussian with fixed FWHM and :math:`\epsilon_i` is the energy of the :math:`i`-th molecular
orbital (MO).  The formula for the PDOS of the :math:`\mu`-th atomic orbital (AO) is

.. math::
   PDOS(E)_\mu = \sum_i q_{\mu i} g(E-\epsilon_i)

where :math:`q_{\mu i}` is the Mulliken charge of the :math:`mu`-th AO in the i-th MO.  It is
calculated as

.. math::
   q_{\mu i} = \sum_\nu S_{\mu \nu} P_{\mu \nu}^i,

where 

.. math::
   S_{\mu \nu} = \langle \mu | \nu \rangle

is the AO overlap matrix and 

.. math::
   P_{\mu \nu}^i = C_{\nu i} C_{\mu i}

is the :math:`i`-th MO density matrix calculated from the MO coefficient matrix, :math:`C`.
Alternatively, we propose to calculate PDOS from Löwdin charges to avoid 
possible negative PDOS. Löwdin charges can be calculated as

.. math::
   q^L_{\mu i} = \sum_\nu S^{1/2}_{\mu \nu} P_{\nu \nu}^i S^{1/2}_{\mu \nu}.

Normally we are interested in the PDOS for a group of orbitals (such as all
the d orbitals on a Ruthenium atom).  In that case, the appropriate PDOS is obtained as
a sum over the PDOS of all relevant orbitals,

.. math::
    PDOS(E) = \sum_\mu PDOS(E)_\mu.


