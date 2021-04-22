# smooth_covmat
Very fast Fortran90 code to make low-resolution smooth covariance matrices from the Planck hitcount maps.

## Usage of code
- define directories: `basedir, mapname, outdat` where high-resolution *Planck* maps are stored (publicly available at: [Planck Archive](http://pla.esac.esa.int))

## compile code
use `exe.make` file on `TELOS`. You need the following installed libraries:
1. fitsio
2. healpix
3. ifort compiler

## code 
- code degrades the map to an intermediate resolution map and computes the smooth low-resolution covariance matrix from there. Therefore, only pixel inside a pre-defined beam width are taken into account. All pixels outside of the beam are set to zero for the computation of pixel i. 

## speed
This code reduced the number of operation from O(N_pix_lowres ** 2 * N_pix_imres) to O(N_pix_lowres ** 2 * N_pix_imres') where N_pix_imres' are only the pixels inside of the beam, meaning that instead of N_pix_imres=12*128*128 we have approximately a few hundred pixels in N_pix_imres' typically which reduced the number of operations by approx. 10000. This increased computational efficiency comes at basically no loss of information for the user. The differences are below 1.e-5 for the user. Previous computation time was approx. 3 days, now 30 minutes (on a single CPU). 

I can share a parallelised version that distribitues each 'larger' pixel on a seperate CPU and then merges them again afterwards: results are below 1 min. 

## cite
Please cite:

@ARTICLE{2021arXiv210314378D,
       author = {{de Belsunce}, Roger and {Gratton}, Steven and {Coulton}, William and {Efstathiou}, George},
        title = "{Inference of the optical depth to reionization from low multipole temperature and polarisation Planck data}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - Cosmology and Nongalactic Astrophysics},
         year = 2021,
        month = mar,
          eid = {arXiv:2103.14378},
        pages = {arXiv:2103.14378},
archivePrefix = {arXiv},
       eprint = {2103.14378},
 primaryClass = {astro-ph.CO},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2021arXiv210314378D},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}


  
