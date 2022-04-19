The main code is in the file "wire3.f90" which is a fortran code trying to solve the RG equations for local Sine-Gordon model.
You can input an arbitrary spectral function, the program finds the Matsubara Green's function (via Kramers–Kronig relations).
The RG equations are then solved numerically using the following fortran libraries:"opkda1.f", "opkda2.f "opkdmain.f", and "quadpack.f90".

To run the code:

1. Open terminal and type what follows
2. gfortran -c wire3.f90
3. gfortran wire3.f90 opkda1.f opkda2.f opkdmain.f quapack.f -o a
4. ./a
5. you will get a txt-file including the UV cuoff and renormalized scattering potential
