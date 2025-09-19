# PPPL-GTS
#
#  20250902 NAG is replaced by netlib quadpack for open source purpose
#  20250902 an alternate is also provided by cmlib
#
   Before you compile GTS, please load the follwoing modules

   on nersc perlmutter:
     PrgEnv-nvidia/8.3.3 
     cray-hdf5-parallel/1.12.2.3

   on pu stellar:
     intel/2021.1.2
     openmpi/intel-2021.1/4.1.0
     fftw/intel-2021.1/openmpi-4.1.0/3.3.9
     hdf5/intel-2021.1/openmpi-4.1.0/1.10.6
     netcdf/intel-2021.1/hdf5-1.10.6/openmpi-4.1.0/4.7.4

   on pppl flux:
     intel/2023.2.0
     openmpi/4.1.6-intel-2023.2.0
     openmpi-4.1.6/intel-2023.2.0/netcdf-c/main-nzscm
     openmpi-4.1.6/intel-2023.2.0/netcdf-fortran/4.6.1-cepgo

   to set up environment. Then install petsc, pspline, sprng packages.

   Then you can compile GTS now. If your system has NAG installed, Please use 
     make _ARCH=[specific setting for your system]
   to compile GTS. example:
     on nersc perlmutter: make NAG=y _ARCH=perlmutter
     on pu stellar:       make NAG=y _ARCH=stellar
     on pppl flux:        make NAG=y _ARCH=flux
   Otherwise, please ignore "NAG=y".


   20220411
https://docs.github.com/en/get-started/getting-started-with-git/managing-remote-repositories#switching-remote-urls-from-https-to-ssh
git remote set-url origin git@github.com:PrincetonUniversity/PPPL-GTS 




 ---------------------------- Below are the OLD records ---------------------------
on SEABORG
module load MASS/3.0
module load sprng64/1.0
module load idl/6.3
module load nag_smp_64/2
module load netcdf/3.6.0
module switch hdf5_par_64/1.6.4 hdf5_par64/1.6.4
module load superlu_mt_64


compile
gmake OPENMP=y 

====================================================================

6/5/2007   revision 48 is the first version of GTS that has the
           particle distribution algorithm running on Seaborg.

====================================================================

7/31/2007  Use propedit to change the svn property on the current directory:

         svn propedit svn:ignore .

This will call vi and open a file that contains the properties for svn:ignore,
which is a list of unrevisioned files that will be ignored by svn. I put
all the *.mod files in there and also the name of the executable "gtc".
Now, when I do an svn status, I don't get a list of the .mod files.

====================================================================

01/18/2013 This is the first svn version
           source code from Weixing Wang on Jan 12, 2013
           It runs on hopper with WORK/JOB_2012/short_shape

====================================================================

01/22/2013 
           GTS_etg_fem with FEM for 2D smoothing dated back to 2008:

             7302 2008-09-09 11:37 fem.F90
            21863 2008-09-09 11:37 triangle.h
           651275 2008-09-09 11:37 triangle.c
            20654 2008-09-09 11:37 therm2d.F90

            28874 2008-09-09 11:37 solver_petsc.F90

            19842 2008-09-09 11:37 smooth.F90

             4668 2008-09-09 11:37 mesh_gen.c
             8538 2008-09-09 11:37 main.F90



====================================================================

Stephane Ethier
2:28 PM (28 minutes ago) 4/3/2013

           to me, Weixing, Edward 
           Hi All,

           I'm not sure if it could help
           but I finally found on an old cvs server the old version of GTC
           that Mark Adams and Yas Nishimura worked on
           to produce a Finite Element version of the solver.
           I haven't tested it so I have no idea if it runs properly.
           It might give you some ideas though.
           You can have a look at the code on portal in the following directory:

                /p/gkp/ethier/GTC/gtc_fem_adams
           cheers,
           Stephane

====================================================================

11/20/2013 
           merge code for WWang rel31-->rel32
           rel31: bug fix for dxdr & dxdz
                  (solver_petsc.F90 L2484 apply MPI_ALLREDUCE on fem%Pijk)

           Status against revision:     31
           ~/gts/trunk/GTS$svn add Er_provider_rescale.F90
           A         Er_provider_rescale.F90
           ~/gts/trunk/GTS$cd ..
           ~/gts/trunk$svn ci GTS -m "merge code from ww"
           Sending        GTS/Er_provider.F90
           Adding         GTS/Er_provider_rescale.F90
           Sending        GTS/Makefile
           Sending        GTS/Theory_jb.f90
           Sending        GTS/charge_ion.F90
           Sending        GTS/field.f90
           Sending        GTS/get_3D_flow.F90
           Sending        GTS/get_Te_fluctuation.F90
           Sending        GTS/get_ZF_shear.F90
           Sending        GTS/get_fluxes.F90
           Sending        GTS/get_fluxes_e.F90
           Sending        GTS/main.F90
           Sending        GTS/module.f90
           Sending        GTS/open_files.f90
           Sending        GTS/poisson.F90
           Sending        GTS/pushe_new.F90
           Sending        GTS/pushi_new.F90
           Sending        GTS/setup_v2.F90
           Sending        GTS/snapshot.f90
           Sending        GTS/solver_petsc.F90
           Sending        GTS/spl_profile.F90
           Sending        GTS/weight_damping_eon_new2.F90
           Sending        GTS/weight_damping_new.F90
           Transmitting file data .......................
           Committed revision 32.
           ~/gts/trunk$cd -
           ~/gts/trunk/GTS$svn update
           At revision 32.

