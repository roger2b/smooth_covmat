ifort -Ofast  -I/usr/local/healpix/include -c "$1".f90
ifort  -o "$1".exe "$1".o -qopenmp -Bstatic -L/usr/local/healpix/lib -lhealpix -lcfitsio 
