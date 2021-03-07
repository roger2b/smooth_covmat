# smooth_covmat
Very fast Fortran90 code to make low-resolution smooth covariance matrices from the Planck hitcount maps.

# Usage of code
- define directories: `basedir, mapname, outdat` where high-resolution *Planck* maps are stored (publicly available at: [Planck Archive](http://pla.esac.esa.int))

# compile code
use `exe.make` file on `TELOS`. You need the following installed libraries:
1. fitsio
2. healpix
3. ifort compiler

# code 
- code degrades the map to an intermediate resolution map and computes the smooth low-resolution covariance matrix from there. Therefore, only pixel inside a pre-defined beam width are taken into account. All pixels outside of the beam are set to zero for the computation of pixel i. 

# cite
Please cite de Belsunce et al. (2021 in prep.) if you use the code
  
