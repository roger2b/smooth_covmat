program processncm

  use healpix_types
  use head_fits
  use fitstools
  use alm_tools
  use pix_tools
  use udgrade_nr

  implicit none

  character(len=200) :: infile, outfile, mapname
  character(len=200) :: outdat, ndx_lr,ndx_im,basedir, freq

  integer, parameter :: nside=16
  integer, parameter :: npix=12*nside*nside
  integer, parameter :: nmaps=10
  integer, parameter :: ppmaps=3
  integer, parameter :: np=2*npix
  real(dp), dimension(0:np-1,0:np-1) :: ncm

  integer :: i,j,k,l,m,ik, jk, nlisti, nlistj
  real(dp), dimension(0:10*nside) :: beam

  real(dp) :: x,  zz1,zz2,angij,thsmooth,zdeg_rad,zbounds(2),ell
  integer :: beamtrunc, icase,icase_map, icase_smooth,lmax

  integer, parameter :: hiresnside=2048
  integer, parameter :: hiresnpix=12*hiresnside*hiresnside
  real(dp), dimension(:,:), allocatable :: hiresmap
  real(dp), dimension(:), allocatable :: hiresqqmap,hiresqumap,hiresuumap

  integer, parameter :: imresnside=128
  integer, parameter :: imresnpix=12*imresnside*imresnside
  real(dp), dimension(:), allocatable :: imresqqmap,imresqumap,imresuumap

  real(dp), dimension(0:imresnpix-1) :: im_x, im_y, im_z
  integer,  dimension(0:imresnpix-1) ::  pix_ik, pix_jk
  real(dp), dimension(0:npix-1) :: lo_x, lo_y, lo_z

  complex(dp), allocatable :: almqq(:,:,:),almqu(:,:,:),almuu(:,:,:)
  complex(dp), allocatable :: alm1(:,:,:),alm2(:,:,:),alm3(:,:,:)
  real(dp), dimension(:,:), allocatable :: weight, pixelwindow

  integer, parameter :: reijonmax=3*nside
  real(dp), parameter :: dreijonmax=dble(reijonmax)
!  real(dp), parameter :: pi=4.d0*DATAN(1.d0)
  real(dp) :: sumqq,sumqu,sumuu
  real(dp) :: argik, argjk, fac_lr, fac_im, rad_disc
  real(dp), dimension(3) :: veci, vecj, veck

  integer :: clock_rate, clock_max, t1,t2

  character(len=80), dimension(1:60) :: header

  print*, '---------------------------------'
  print*, 'compute low resolution cov matrix'
  print*, '---------------------------------'
    print*, ' '
  print*,'pi',  pi
! exact computation for icase = 1
! approximate healpix disc method for icase=2
  icase = 2
  zdeg_rad = 2.d+00*pi/360.d+00
  thsmooth = 7.0d+00*zdeg_rad*0.425
!  icase_smooth = 1 for Gaussian smoothing 
!  icase_smooth = 2 for cosine smoothing with characteristic scale Nside_lr
!  icase_smooth = 3 for cosine smoothing and pixelwindow correction
  icase_smooth = 3


  freq='100'
  freq='143'
  freq='30'
  write (ndx_lr,"(I4.4)") nside
  write (ndx_im,"(I4.4)") imresnside

  print*, 'Enter 1 for Planck 2015 (FFP8)'
  print*, 'Enter 2 for SROLL2.0'
  print*, 'Enter 3 for Planck 2018 (FFP10)'
  read*, icase_map

  print*, 'Enter frequency: 100, 143, 217 or 353'
  read(*,'(A)') freq
  print*, freq
  if (icase_map.eq.1) then
  basedir='/data2/rmvd2/dx11/Planck_2015/'
  mapname=TRIM(basedir)//'frequency_maps/HFI_SkyMap_'//TRIM(freq)//'_2048_R2.02_full.fits'
  outdat=TRIM(basedir)//'ncm/DX11d_'//TRIM(freq)//'_ncm_nside'//TRIM(ndx_lr)//'_im'//TRIM(ndx_im)//'_PP_smooth_pixwin.dat'
  else if (icase_map.eq.2) then
  basedir='/data2/rmvd2/dx11/SROLL2.0/'
  mapname=TRIM(basedir)//'SRoll20_SkyMap_'//TRIM(freq)//'psb_full.fits'
  outdat=TRIM(basedir)//'/SROLL20_'//TRIM(freq)//'_ncm_nside'//TRIM(ndx_lr)//'_im'//TRIM(ndx_im)//'_PP_smooth_pixwin.dat'
  else if (icase_map.eq.3) then
  basedir='/data2/rmvd2/dx11/Planck_FFP10/frequency_maps/'
  !mapname=TRIM(basedir)//'HFI_SkyMap_'//TRIM(freq)//'-psb_2048_R3.01_full.fits'
  mapname=TRIM(basedir)//'HFI_SkyMap_'//TRIM(freq)//'_2048_R3.01_full.fits'
  !mapname=TRIM(basedir)//'LFI_SkyMap_030-BPassCorrected_1024_R3.00_full.fits'
  !mapname=TRIM(basedir)//'LFI_SkyMap_030-BPassCorrected_1024_R3.00_full-ringhalf-1.fits'
  !mapname=TRIM(basedir)//'HFI_SkyMap_'//TRIM(freq)//'_2048_R3.01_hm2.fits'
  outdat=TRIM(basedir)//'/FFP10_'//TRIM(freq)//'full_ncm_nside'//TRIM(ndx_lr)//'_im'//TRIM(ndx_im)//'_PP_smooth_pixwin.dat'
  endif

  print*, 'input high resolution map', TRIM(mapname)
  print*, 'degrade to NSIDE=', nside
  print*, 'intermediate resolution map NSIDE=', imresnside
  print*, 'save file in: ',TRIM(outdat)
  beam=0.d0

  !reijo...
  beamtrunc=3
  beam(0:nside)=1.d0
  do i=nside+1,3*nside
       beam(i)=.5d0*(1.d0+dcos((dble(i)-dble(nside))*pi/(2.d0*dble(nside))))
  enddo

  print*, 'allocate arrays for QQ,QU,UU maps'
 ! allocate arrays for QQ,QU,UU maps
  allocate(hiresmap(0:hiresnpix-1,nmaps))

  allocate(hiresqqmap(0:hiresnpix-1))
  allocate(hiresqumap(0:hiresnpix-1))
  allocate(hiresuumap(0:hiresnpix-1))

  allocate(imresqqmap(0:imresnpix-1))
  allocate(imresqumap(0:imresnpix-1))
  allocate(imresuumap(0:imresnpix-1))
  print*, 'sucessfully allocated'

  call input_map(mapname, hiresmap, hiresnpix, nmaps)
  print*, 'read in maps - ok'

  call convert_nest2ring(hiresnside, hiresmap)
  print *,'Loaded hires map in RING- ok'

  do i=0,hiresnpix-1
   hiresqqmap(i)=hiresmap(i,8)
   hiresqumap(i)=hiresmap(i,9)
   hiresuumap(i)=hiresmap(i,10)
  enddo

  fac_lr = (dble(hiresnside)/dble(nside))**2.d0
  fac_im = (dble(hiresnside)/dble(imresnside))**2.d0

! Allocate alm array
  lmax = MAX(4*imresnside,60)

  print*, 'Reijo smoothing and degrading of maps'
  allocate(pixelwindow(0:lmax, ppmaps))
  pixelwindow=0.d+00
  call pixel_window(pixelwindow,imresnside)
  pixelwindow(0,2)=1.d0
  pixelwindow(1,2)=1.d0

  allocate(almqq(1, 0:lmax, 0:lmax))
  allocate(almqu(1, 0:lmax, 0:lmax))
  allocate(almuu(1, 0:lmax, 0:lmax))
  allocate(weight(2*hiresnside, 1))

  zbounds = (/ -1.0_dp, 1.0_dp /)
  weight = 1.0_dp

  print *, 'Doing spherical transform'
  call map2alm(hiresnside, lmax, lmax, hiresqqmap, almqq, zbounds, weight)
  call map2alm(hiresnside, lmax, lmax, hiresqumap, almqu, zbounds, weight)
  call map2alm(hiresnside, lmax, lmax, hiresuumap, almuu, zbounds, weight)

  allocate(alm1(1, 0:lmax, 0:lmax))
  allocate(alm2(1, 0:lmax, 0:lmax))
  allocate(alm3(1, 0:lmax, 0:lmax))

  zz2=1.d0
  do l = 0, lmax
     ell = dfloat(l)
     if(icase_smooth.eq.1) call f_ell(thsmooth, ell, zz2)
     if(icase_smooth.eq.2) call f_ell_reijo(l, imresnside, pi, zz2)
     if(icase_smooth.eq.3) then
        call f_ell_reijo(l, imresnside, pi, zz2)
        zz2 = zz2*pixelwindow(l,2)
     endif
     if(l.lt.0) zz2 = 0.d+00
     print*, l, zz2, zz2*pixelwindow(l,2)
     do m = 0, l
      alm1(1, l, m)  = almqq(1, l, m)*zz2
      alm2(1, l, m)  = almqu(1, l, m)*zz2
      alm3(1, l, m)  = almuu(1, l, m)*zz2
     end do
  end do

  call alm2map(imresnside, lmax, lmax, alm1, imresqqmap)
  call alm2map(imresnside, lmax, lmax, alm2, imresqumap)
  call alm2map(imresnside, lmax, lmax, alm3, imresuumap)

  print *,'Degraded hires map to im res'
  !print*, 'apply factor to degraded map:', fac_im
  !call udgrade_ring(hiresttmap,hiresnside,imresttmap,imresnside)
  !imresttmap = imresttmap*fac_im

  do i=0, npix-1
     call pix2vec_ring(nside, i, veci)
     lo_x(i)=veci(1)
     lo_y(i)=veci(2)
     lo_z(i)=veci(3)
  enddo
  print *,'done lo res angles'

  do i=0, imresnpix-1
     call pix2vec_ring(imresnside, i, veci)
     im_x(i)=veci(1)
     im_y(i)=veci(2)
     im_z(i)=veci(3)
  enddo
  print *,'done im res angles'


  ncm = 0.d0

  if (icase.eq.1) then
  call system_clock ( t1, clock_rate, clock_max )

  do j=0,npix-1
     do i=j,npix-1
        sumqq=0.d0
        sumqu=0.d0
        sumuu=0.d0
        do k=0,imresnpix-1
           argik=lo_x(i)*im_x(k)+lo_y(i)*im_y(k)+lo_z(i)*im_z(k)
           argjk=lo_x(j)*im_x(k)+lo_y(j)*im_y(k)+lo_z(j)*im_z(k)
           zz1=calccorr(argik)
           zz2=calccorr(argjk)
           sumqq=sumqq+zz1*zz2*imresqqmap(k)
           sumqu=sumqu+zz1*zz2*imresqumap(k)
           sumuu=sumuu+zz1*zz2*imresuumap(k)
        enddo
        ncm(i     ,j     )=sumqq
        ncm(i+npix,j     )=sumqu
        ncm(i+npix,j+npix)=sumuu
     enddo
     if (MOD(j,100).eq.0) then
        print*, 'done for row=', j
     endif
  enddo
  call system_clock ( t2, clock_rate, clock_max )
  write ( *, * ) 'Elapsed real time = ', real ( t2 - t1 ) / real (clock_rate),'s'

  else if (icase.eq.2) then
! compute 

  rad_disc = 4.d0*pi/3.d0/dble(nside)
  print*, 'rad disc', rad_disc

  call system_clock ( t1, clock_rate, clock_max )

  do j=0,npix-1
     call query_disc(imresnside,(/lo_x(j),lo_y(j),lo_z(j)/),rad_disc,pix_jk,nlistj)!,inclusive=1)
     do i=j,npix-1
        angij = lo_x(i)*lo_x(j)+lo_y(i)*lo_y(j)+lo_z(i)*lo_z(j)
        if (angij.gt.DCOS(2.d0*rad_disc)) then
        sumqq=0.d0
        sumqu=0.d0
        sumuu=0.d0
        do k=0,nlistj-1
           jk = pix_jk(k)
           argik=lo_x(i)*im_x(jk)+lo_y(i)*im_y(jk)+lo_z(i)*im_z(jk)
           argjk=lo_x(j)*im_x(jk)+lo_y(j)*im_y(jk)+lo_z(j)*im_z(jk)
           zz1=calccorr(argik)
           zz2=calccorr(argjk)
           sumqq=sumqq+zz1*zz2*imresqqmap(jk)
           sumqu=sumqu+zz1*zz2*imresqumap(jk)
           sumuu=sumuu+zz1*zz2*imresuumap(jk)
        enddo
        if (i.eq.j) write(*,"(I6,I6,E12.5,E12.5,E12.5)"), i, j, sumqq, sumqu, sumuu
        ncm(i     ,j     )=sumqq
        ncm(i+npix,j     )=sumqu
        ncm(i+npix,j+npix)=sumuu
        endif
     enddo
     if (MOD(j,100).eq.0) then
        print*, 'done for row=', j
     endif
  enddo

  call system_clock ( t2, clock_rate, clock_max )
  write ( *, * ) 'Elapsed real time = ', real ( t2 - t1 ) / real (clock_rate),'s'
  endif
! mirror upper right quadrant before mirroring entire matrix to make it
! symmetric
    do j=0,npix-1
     do i=j,npix-1
        ncm(j+npix,i) = ncm(i+npix,j)
     enddo
  enddo

    do j=0,np-1
     do i=j,np-1
        ncm(j,i) = ncm(i,j)
     enddo
  enddo

  print*, 'save file'
  open(11,file=outdat,form='UNFORMATTED',status='unknown')
    write(11) ((ncm(i,j), i=0,np-1), j=0, np-1)
  close(11)

!  open(13,file=outdat_format,form='FORMATTED')
!    write(13,'(ES20.10)') ((tt(i,j), i=0,npix-1), j=0, npix-1)
!  close(13)

contains


elemental function calccorr(x)
    real(dp) :: calccorr
    real(dp), intent(in) :: x
    integer :: k
    real(dp) :: y_k, y_kp1, y_kp2, rkp1, rkp2, kd, dkp1, d2kp1

    y_kp2=0.d0
    y_kp1=0.d0
    y_k=0.d0
    kd=dreijonmax
    do k=reijonmax,1,-1
       y_kp2=y_kp1
       y_kp1=y_k
       dkp1=kd+1.d0
       d2kp1=2*kd+1.d0
       rkp1=1.d0/dkp1
       rkp2=1.d0/(dkp1+1.d0)
       y_k=(rkp1*x*y_kp1+beam(k))*d2kp1-rkp2*dkp1*y_kp2
       kd=kd-1.d0
    enddo

    calccorr=-.5d0*y_kp1+x*y_k+beam(0)
    calccorr=calccorr/(4.d0*pi)
  end function calccorr

  subroutine f_ell(theta_s, ell, zz)
  implicit none
  real(kind=dp) :: zz, theta_s, ell
     zz = dexp(-0.5d+00*(theta_s*ell)**2)
  end subroutine f_ell

  subroutine f_ell_reijo(l, Nside, zpi, zz2)
  double precision zpi, zz2
  integer Nside, l
     if(l.le.Nside) then
     zz2 = 1.d+00
     return
     end if
     if(l.gt.Nside.and.l.le.3*Nside) then
  !         if(l.gt.Nside.and.l.le.31) then
     zz2 = 0.5d+00*(1.d+00 +  dcos(dfloat(l-Nside)*zpi/2.d+00/dfloat(Nside)))
     return
     end if
     if(l.gt.3*Nside) then
  !         if(l.gt.31) then
     zz2 = 0.d+00
     return
     end if
  end subroutine f_ell_reijo


end program processncm



