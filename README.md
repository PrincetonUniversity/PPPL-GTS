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
