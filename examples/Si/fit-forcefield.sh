SCRIPT=../../scripts/csld_main
# first fit force field coefficients only
$SCRIPT --ldff_step 2 -f csld.in-forcefield --phonon_step -1 
# next read coefficients and fit FCT only
$SCRIPT --ldff_step 2 -f csld.in-forcefield --phonon_step -1 --symC_step 1 --train_step 1 \
  --override '[fitting] solution_known = sol_FF' \
  --override '[fitting] solution_out = sol_LD+FF' \
  --override '[fitting] submodel1 = FCT 0 1 2 3 4'
# now the obtained phonon doesn't make sense, but the FF potential should look okay.

# alternatively, 1 fit orders 0,1,2 for phonons; 2. fit FF (order "5"); 3. fit orders 3 4
# such that the phonon spectrum is correct but the FF looks weird.
# your choice