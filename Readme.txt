Include this directory './mfiles' into Matlab path
as well as Matlab MIRT directories.

I. Overview of the process (4 steps):

	1. Registration   ---  opt_grad_non_rgd_jac.m
        2. Estimation of setup difference --- est_setup.m
	3. Simplified DOV (Proportionality model) --- dovinit_cor.m
	4. LS fitting to BSP motion model --- after_dovinit.m
	5. B-spline based DOV --- opt_4d_cg_cor_nt.m

II. inputs and outputs for each .m file

1. opt_grad_non_rgd_jac.m   
   In: 2 CT volumes     
   Out: a deformation field, spatial knot locations, spatial knot coefficients

2. est_setup.m
   In: a reference volume, three projection views
   out: rotation angles and translations

2. dovinit_cor.m
   In: a reference volume, projection views and the deformation field from registration
   Out: proportionality sequence \alpha, 

3. after_dovinit.m 
   In: \alpha, spatial knot locations, spatial knot coefficients.
   Out: spatial and temporal knot locations of the B-spline motion model, and the corresponding knot coefficients

4. opt_4d_cg_cor_nt.m
   In: a reference volume, projection views, the output from after_dovinit.m  
   Out: knot coefficients


*III. Notes on data processing

     1.  Make sure the projection files are read out in the right order, sometimes the order is messed up by the way the files are numbered.
