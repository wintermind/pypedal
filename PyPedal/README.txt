PyPedal 2 CHANGELOG
===================

CHANGES in PyPedal 2.0.1 (Vetinari)
===================================
12/30/2010    Wrapped a stray debugging message in pyp_nrm/inbreeding() so that
it will only print when the debugging_messages option is set.

12/30/2010    An incorrect index was being used in pyp_nrm/inbreeding() to
match-up COI with animal IDs when COI are read from an existing NRM. The
procedure was looping from 0 to (number of animals - 1), but 1 was being
subtracted from the index when reaching into the pedigree to assign the COI.
This means that the COI for animal 0 (first in the pedigree) was being
assigned to animal -1 (the last animal in the pedigree), the COI for animal 1
was being assigned to animal 0, etc. This bug affected any inbreeding
calculations based on an attached NRM, regardless of the pedformat.

12/30/2010    Found a bug in pyp_nrm/inbreeding() that affects ASD pedigrees to
which a NRM is attached at load time. When the COI are copied from the NRM to
the results dictionary the wrong COI is being associated with an animal ID.

12/30/2010    Added a new flag, force, to pyp_nrm/inbreeding() to override the
use of an attached NRM for finding COI. This is needed for debugging a possible
problem with correctly mapping IDs in ASD pedigrees when NRM are attached.
