**wetSAXS**
C++ program for calculating small angle X-ray scattering (SAXS) profiles from atomic coordinates in PDB format.  The program uses the NNLS code from Suvrit Sra, the frsc directory for ensemble fitting.  

# Build Requirements
Boost c++ library 1.86 or greater
SASTools (https://github.com/rambor/SASTools)

The Boost and SASTools must be built and linked with the same compiler, otherwise you will have symbol issues.  On a MAC with default clang installation and an additional GNU gcc compiler, it is easy to mix these up so you need to specify boost to use the same gcc compiler if using gcc for wetSAXS and SASTools

On a MAC, I would build boost and SAStools using a homebrew install of gcc, libaries are built static.
