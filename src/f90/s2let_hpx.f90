! ----------------------------------------------------------- !

SUBROUTINE healpix_forward_real_noiter_noweights( alm, map, nside, L )

  use healpix_types
  use alm_tools
  use healpix_modules
  implicit none

  integer(i4b)  :: nside, nlmax, L, el, em
  integer(i8b) :: npix
  real(dp), dimension(0:12*nside*nside-1) :: map
  complex(dp), dimension(0:L*L-1)     :: alm
  complex(dpc), allocatable, dimension(:,:,:)  ::  alm_temp

  nlmax = L - 1
  npix = nside2npix(nside)

  allocate(alm_temp(1:1, 0:nlmax, 0:nlmax))
  alm_temp = 0.0_dpc

  call map2alm(nside, nlmax, nlmax, map(0:npix-1), &
       alm_temp(1:1,0:nlmax, 0:nlmax))

  do el = 0, nlmax
     do em = 0, el
        alm( el * el + el - em ) =  (-1.0)**real(-em) * CONJG( alm_temp(1, el, em) )
        alm( el * el + el + em ) = alm_temp(1, el, em)
     enddo
  enddo

  deallocate(alm_temp)

END SUBROUTINE healpix_forward_real_noiter_noweights

! ----------------------------------------------------------- !

SUBROUTINE healpix_forward_real( alm, map, nside, L )

  use healpix_types
  use alm_tools
  use healpix_modules
  implicit none

  integer(i4b), parameter  :: nbiter = 3

  integer(i4b) 	:: nside, iter, n_plm, nlmax, L, el, em
  integer(i8b) :: npix
  real(dp), dimension(0:12*nside*nside-1) :: map
  complex(dp), dimension(0:L*L-1) 	  :: alm
  complex(dpc), allocatable, dimension(:,:,:)  ::  alm_temp, alm_f
  real(dp), allocatable, dimension(:)  ::  map_temp
  real(dp), dimension(:,:), allocatable :: plm
  real(dp), allocatable, dimension(:,:) :: dw8
  real(dp), dimension(2) :: z
  character(len=FILENAMELEN) :: def_dir, def_file, final_file
  logical :: filefound

  nlmax = L - 1
  npix = nside2npix(nside)
  n_plm =  (nlmax+1) * (nlmax+2) * nside

  allocate(dw8(1:2*nside, 1:1))
  allocate(alm_temp(1:1, 0:nlmax, 0:nlmax))
  allocate(alm_f(1:1, 0:nlmax, 0:nlmax))
  ALLOCATE(map_temp(0:npix-1))
  allocate(plm(0:n_plm-1,1:1))
  dw8 = 1.0_dp
  alm_temp = 0.0_dpc
  alm_f = 0.0_dpc
  map_temp = 0.0_dp
  z = (/-1.d0 , 1.d0/)
  plm = 0.0_dp

  if( .true. ) then
    def_dir = get_healpix_data_dir()
    def_file = trim(get_healpix_ring_weight_file(nside))
    filefound = scan_directories(def_dir, def_file, final_file)
    if( .not. filefound ) then
      PRINT*, "Ring weight file could not be found! Check healpix directories and files."
    else
      call input_map(final_file, dw8, 2*nside, 1, fmissval=0.0_dp)
      dw8 = 1.0_dp + dw8
    endif
  endif

  call plm_gen(nside, nlmax, nlmax, plm)
  map_temp(0:npix-1) = map(0:npix-1)

  do iter = 0, nbiter

    call map2alm(nside, nlmax, nlmax, map_temp(0:npix-1), &
       alm_temp(1:1,0:nlmax,0:nlmax), z, dw8, plm(0:n_plm-1,1))

    alm_f = alm_f + alm_temp
    map_temp = 0.0_dp
    call alm2map(nside, nlmax, nlmax, alm_f, map_temp, plm(:,1))
    map_temp = map - map_temp

  enddo

    do el = 0, nlmax
       do em = 0, el
          alm( el * el + el - em ) = (-1.0)**real(-em) * CONJG( alm_f(1, el, em) )
          alm( el * el + el + em ) = alm_f(1, el, em)
       enddo
    enddo

  deallocate(dw8, alm_temp, alm_f, map_temp, plm)

END SUBROUTINE healpix_forward_real

! ----------------------------------------------------------- !

SUBROUTINE healpix_inverse_real( map, alm, nside, L )

  use healpix_types
  use alm_tools
  use healpix_modules
  implicit none

  integer(i4b) 	:: nside, nlmax, L, el, em
  integer(i8b) :: npix
  real(dp), dimension(0:12*nside*nside-1) :: map
  complex(dp), dimension(0:L*L-1) 	  :: alm
  complex(dpc), allocatable, dimension(:,:,:)  ::  alm_temp, alm_temp2

  nlmax = L - 1
  npix = nside2npix(nside)

  allocate(alm_temp(1:1, 0:nlmax, 0:nlmax))
  alm_temp = 0.0_dpc

  do el = 0, nlmax
     do em = 0, el
        alm_temp(1, el, em) = alm( el * el + el + em )
     enddo
  enddo

  call alm2map(nside, nlmax, nlmax, alm_temp, map)

  !allocate(alm_temp2(1:1, 0:nlmax, 0:nlmax))
  !alm_temp2 = 0.0_dpc
  !call map2alm(nside, nlmax, nlmax, map, alm_temp2)
  !print*, SUM(CDABS(alm_temp2-alm_temp))/real((nlmax+1)*(nlmax+2)/2.0)
  !alm_temp2)

  deallocate(alm_temp)

END SUBROUTINE healpix_inverse_real


! ----------------------------------------------------------- !

SUBROUTINE read_healpix_maps(maps, file, nside, nmaps)

  use healpix_types
  use healpix_modules

  character(len=filenamelen) :: file, ordering, comment
  character(LEN=80), dimension(1:60) :: header
  integer :: n, nside, npix, nmaps
  real(dp), dimension(0:12*nside*nside-1,1:nmaps) :: maps
  real(dp) :: nullval
  logical :: anynull

  npix = nside2npix(nside)
  call read_bintab(file, maps, npix, nmaps, nullval, anynull, header)
  call get_card(header,'ORDERING',ordering,comment)
  if( trim(ordering) .ne. 'ring' .and. trim(ordering) .ne. 'RING' ) then
    do n = 1, nmaps
     call convert_nest2ring(nside, maps(:,n) )
    enddo
  endif

END SUBROUTINE read_healpix_maps


! ----------------------------------------------------------- !

SUBROUTINE read_healpix_map(map, file, nside)

  use healpix_types
  use healpix_modules

  character(len=filenamelen) :: file, ordering, comment
  character(LEN=80), dimension(1:60) :: header
  integer :: nside, npix
  real(dp), dimension(0:12*nside*nside-1,1:1) :: map
  real(dp) :: nullval
  logical :: anynull

  npix = nside2npix(nside)
  call read_bintab(file, map, npix, 1, nullval, anynull, header)
  call get_card(header,'ORDERING',ordering,comment)
  if( trim(ordering) .ne. 'ring' .and. trim(ordering) .ne. 'RING' ) then
     call convert_nest2ring(nside, map)
     !PRINT*, "Input Healpix map must be RING ordered"
     !PRINT*, "End of program..."
     !call exit(0)
  endif

END SUBROUTINE read_healpix_map

! ----------------------------------------------------------- !

SUBROUTINE write_healpix_map(file, map, nside)

  use healpix_types
  use healpix_modules

  character(len=filenamelen) :: file
  CHARACTER(len=80), DIMENSION(1:120):: header
  integer :: nside, npix
  real(dp), dimension(0:12*nside*nside-1,1:1) :: map

  npix = nside2npix(nside)
  header = ''
  call write_minimal_header(header, 'MAP', nside=nside, ordering='ring')
  call write_bintab(map, npix, 1, header, 120, file)

END SUBROUTINE write_healpix_map

! ----------------------------------------------------------- !


! ----------------------------------------------------------- !

SUBROUTINE healpix_forward_spin_real( almE, almB, mapQ, mapU, nside, L, spin )

  use healpix_types
  use alm_tools
  use healpix_modules
  implicit none

  integer(i4b), parameter  :: nbiter = 3

  integer(i4b)  :: nside, iter, n_plm, nlmax, L, el, em, spin
  integer(i8b) :: npix
  real(dp), dimension(0:12*nside*nside-1) :: mapQ, mapU
  complex(dp), dimension(0:L*L-1)     :: almE, almB
  complex(dpc), allocatable, dimension(:,:,:)  ::  alm_temp, alm_f
  real(dp), allocatable, dimension(:,:)  ::  map_temp
  real(dp), dimension(:,:), allocatable :: plm
  real(dp), allocatable, dimension(:,:) :: dw8
  real(dp), dimension(2) :: z
  character(len=FILENAMELEN) :: def_dir, def_file, final_file
  logical :: filefound

  nlmax = L - 1
  npix = nside2npix(nside)

  allocate(dw8(1:2*nside, 1:1))
  allocate(alm_temp(1:2, 0:nlmax, 0:nlmax))
  allocate(alm_f(1:2, 0:nlmax, 0:nlmax))
  ALLOCATE(map_temp(0:npix-1,1:2))
  dw8 = 1.0_dp
  alm_temp = 0.0_dpc
  alm_f = 0.0_dpc
  map_temp = 0.0_dp
  z = (/-1.d0 , 1.d0/)

  if( .true. ) then
    def_dir = get_healpix_data_dir()
    def_file = trim(get_healpix_ring_weight_file(nside))
    filefound = scan_directories(def_dir, def_file, final_file)
    if( .not. filefound ) then
      PRINT*, "Ring weight file could not be found! Check healpix directories and files."
    else
      call input_map(final_file, dw8, 2*nside, 1, fmissval=0.0_dp)
      dw8 = 1.0_dp + dw8
    endif
  endif

  map_temp(0:npix-1,1) = mapQ(0:npix-1)
  map_temp(0:npix-1,2) = mapU(0:npix-1)

  do iter = 0, nbiter

    call map2alm_spin(nside, nlmax, nlmax, spin, map_temp(0:npix-1,1:2), &
       alm_temp(1:2,0:nlmax,0:nlmax), z, dw8)

    alm_f = alm_f + alm_temp
    map_temp = 0.0_dp
    call alm2map_spin(nside, nlmax, nlmax, spin, alm_f, map_temp)
    map_temp(:,1) = mapQ - map_temp(:,1)
    map_temp(:,2) = mapU - map_temp(:,2)

  enddo

    do el = abs(spin), nlmax
       do em = 0, el

          almE( el * el + el + em ) = alm_f(1, el, em)

          if( em .gt. 0 ) almE( el * el + el - em ) = (-1.0)**real(-em) * CONJG( almE( el * el + el + em ) )

          almB( el * el + el + em ) = alm_f(2, el, em)

          if( em .gt. 0 ) almB( el * el + el - em ) = (-1.0)**real(-em) * CONJG( almB( el * el + el + em ) )

       enddo
    enddo

  deallocate(dw8, alm_temp, alm_f, map_temp)

END SUBROUTINE healpix_forward_spin_real

! ----------------------------------------------------------- !

SUBROUTINE healpix_inverse_spin_real( mapQ, mapU, almE, almB, nside, L, spin )

  use healpix_types
  use alm_tools
  use healpix_modules
  implicit none

  integer(i4b)  :: nside, nlmax, L, el, em, spin
  integer(i8b) :: npix
  real(dp), dimension(0:12*nside*nside-1) :: mapQ, mapU
  complex(dp), dimension(0:L*L-1)     :: almE, almB
  complex(dpc), allocatable, dimension(:,:,:)  ::  alm_temp, alm_temp2
  real(dp), allocatable, dimension(:,:)  ::  map_temp

  nlmax = L - 1
  npix = nside2npix(nside)

  ALLOCATE(map_temp(0:npix-1,1:2))

  allocate(alm_temp(1:2, 0:nlmax, 0:nlmax))
  alm_temp = 0.0_dpc

  do el = abs(spin), nlmax
     do em = 0, el
        alm_temp(1, el, em) = almE(el*el+el+em)
        alm_temp(2, el, em) = almB(el*el+el+em)
     enddo
  enddo

  call alm2map_spin(nside, nlmax, nlmax, spin, alm_temp, map_temp)

  mapQ = map_temp(:,1)
  mapU = map_temp(:,2)

  deallocate(alm_temp, map_temp)

END SUBROUTINE healpix_inverse_spin_real

! ----------------------------------------------------------- !
