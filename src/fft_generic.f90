!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2012 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This is the 'generic' implementation of the FFT library

module decomp_2d_fft
  
  use decomp_2d, only : mytype, DECOMP_INFO, &         ! TYPES
                        nrank, DECOMP_2D_COMM_CART_X,& ! VARIABLES common to all "grids"
                        nx_global, ny_global, nz_global, & ! for pointer association and if decomp_2d_fft_init with NO ARGUMENT
                        xstart,xend,xsize,ystart,yend,ysize,zstart,zend,zsize, &  ! for pointer association
                        decomp_info_init, decomp_info_finalize, decomp_2d_abort, & 
                        transpose_y_to_z, transpose_y_to_x, transpose_x_to_y, transpose_z_to_y  
  use glassman
  
  implicit none
  
  private        ! Make everything private unless declared public

  ! engine-specific global variables
  complex(mytype), pointer, dimension(:) :: buf, scratch
  
  ! Derived-type gathering informations for 2decomp_2D_fft for multigrid purpose
  ! Specific variables
  type DECOMP_FFT_MULTIGRID_SPEC
    complex(mytype), allocatable, dimension(:) :: buf, scratch
  end type DECOMP_FFT_MULTIGRID_SPEC
  
  type(DECOMP_FFT_MULTIGRID_SPEC), allocatable, dimension(:), target :: FFT_multigrid_spec
  type(DECOMP_FFT_MULTIGRID_SPEC), allocatable, dimension(:)         :: FFT_multigrid_spec_tmp ! for reallocation purpose (move_alloc)

  ! common code used for all engines, including global variables, 
  ! generic interface definitions and several subroutines  
#include "fft_common.f90"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time initialisations for the FFT engine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_fft_engine(Igrid)

    implicit none
    
    integer, intent(IN) :: Igrid
    integer             :: cbuf_size
    integer             :: Ngrid0

    if (nrank==0) then
       write(*,*) ' '
       write(*,*) '***** Using the generic FFT engine *****'
       write(*,*) ' '
    end if

    ! allocate FFT_multigrid if first call : size = Igrid
    if (initialised == 0) then
      allocate(FFT_multigrid_spec(Igrid))
    end if

    Ngrid0 = size(FFT_multigrid_spec)

    ! reallocate FFT_multigrid_spec if Igrid > size(FFT_multigrid_spec) : size = Igrid
    if (Igrid > Ngrid0) then
       allocate(FFT_multigrid_spec_tmp(Igrid))
       FFT_multigrid_spec_tmp(1:Ngrid0) = FFT_multigrid_spec
       call move_alloc(FROM=FFT_multigrid_spec_tmp,TO=FFT_multigrid_spec)
    end if

    if (.not. allocated(FFT_multigrid_spec(Igrid)%buf)) then ! initialize if not initialized before 
      cbuf_size = max(FFT_multigrid(Igrid)%ph%xsz(1), FFT_multigrid(Igrid)%ph%ysz(2))
      cbuf_size = max(cbuf_size, FFT_multigrid(Igrid)%ph%zsz(3))
      allocate(FFT_multigrid_spec(Igrid)%buf(cbuf_size))
      allocate(FFT_multigrid_spec(Igrid)%scratch(cbuf_size))
      buf  => FFT_multigrid_spec(Igrid)%buf
      scratch => FFT_multigrid_spec(Igrid)%scratch
    else
      call associate_pointers_decomp_2d_fft_spec(Igrid) ! associate pointers if already initialized
    end if
    
    return
  end subroutine init_fft_engine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! associate specific pointers for decomp_2d_fft 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine associate_pointers_decomp_2d_fft_spec(Igrid)

  implicit none
  integer,intent(in) :: Igrid
  
  integer            :: errorcode

  if (Igrid > size(FFT_multigrid)) then
       errorcode = 4
       call decomp_2d_abort(errorcode, &
            'decomp_2d_fft, associate_pointers_decomp_2d_fft_spec : Igrid > size(FFT_multigrid_spec)')
  end if
  
  buf      => FFT_multigrid_spec(Igrid)%buf
  scratch  => FFT_multigrid_spec(Igrid)%scratch

  end subroutine associate_pointers_decomp_2d_fft_spec


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time finalisations for the FFT engine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine finalize_fft_engine

    implicit none
    
    integer :: Igrid
    
    do Igrid = 1, size(FFT_multigrid_spec)

       deallocate(FFT_multigrid_spec(Igrid)%buf,FFT_multigrid_spec(Igrid)%scratch)
    
    end do
    
    deallocate(FFT_multigrid_spec)

    return
  end subroutine finalize_fft_engine


  ! Following routines calculate multiple one-dimensional FFTs to form 
  ! the basis of three-dimensional FFTs.

  ! c2c transform, multiple 1D FFTs in x direction
  subroutine c2c_1m_x(inout, isign, decomp)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: i,j,k
    
    do k=1,decomp%xsz(3)
       do j=1,decomp%xsz(2)
          do i=1,decomp%xsz(1)
             buf(i) = inout(i,j,k)
          end do
          call spcfft(buf,decomp%xsz(1),isign,scratch)
          do i=1,decomp%xsz(1)
             inout(i,j,k) = buf(i)
          end do
       end do
    end do

    return

  end subroutine c2c_1m_x

  ! c2c transform, multiple 1D FFTs in y direction
  subroutine c2c_1m_y(inout, isign, decomp)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: i,j,k

    do k=1,decomp%ysz(3)
       do i=1,decomp%ysz(1)
          do j=1,decomp%ysz(2)
             buf(j) = inout(i,j,k)
          end do
          call spcfft(buf,decomp%ysz(2),isign,scratch)
          do j=1,decomp%ysz(2)
             inout(i,j,k) = buf(j)
          end do
       end do
    end do

    return

  end subroutine c2c_1m_y

  ! c2c transform, multiple 1D FFTs in z direction
  subroutine c2c_1m_z(inout, isign, decomp)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: i,j,k

    do j=1,decomp%zsz(2)
       do i=1,decomp%zsz(1)
          do k=1,decomp%zsz(3)
             buf(k) = inout(i,j,k)
          end do
          call spcfft(buf,decomp%zsz(3),isign,scratch)
          do k=1,decomp%zsz(3)
             inout(i,j,k) = buf(k)
          end do
       end do
    end do

    return

  end subroutine c2c_1m_z

  ! r2c transform, multiple 1D FFTs in x direction
  subroutine r2c_1m_x(input, output)

    implicit none

    real(mytype), dimension(:,:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:,:), intent(OUT) :: output

    integer :: i,j,k, s1,s2,s3, d1

    s1 = size(input,1)
    s2 = size(input,2)
    s3 = size(input,3)
    d1 = size(output,1)

    do k=1,s3
       do j=1,s2
          ! Glassman's FFT is c2c only, 
          ! needing some pre- and post-processing for r2c
          ! pack real input in complex storage
          do i=1,s1
             buf(i) = cmplx(input(i,j,k),0._mytype, kind=mytype)
          end do
          call spcfft(buf,s1,-1,scratch)
          ! note d1 ~ s1/2+1
          ! simply drop the redundant part of the complex output
          do i=1,d1
             output(i,j,k) = buf(i)
          end do
       end do
    end do

    return

  end subroutine r2c_1m_x

  ! r2c transform, multiple 1D FFTs in z direction
  subroutine r2c_1m_z(input, output)

    implicit none

    real(mytype), dimension(:,:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:,:), intent(OUT) :: output

    integer :: i,j,k, s1,s2,s3, d3

    s1 = size(input,1)
    s2 = size(input,2)
    s3 = size(input,3)
    d3 = size(output,3)

    do j=1,s2
       do i=1,s1
          ! Glassman's FFT is c2c only, 
          ! needing some pre- and post-processing for r2c
          ! pack real input in complex storage
          do k=1,s3
             buf(k) = cmplx(input(i,j,k),0._mytype, kind=mytype)
          end do
          call spcfft(buf,s3,-1,scratch)
          ! note d3 ~ s3/2+1
          ! simply drop the redundant part of the complex output
          do k=1,d3
             output(i,j,k) = buf(k)
          end do
       end do
    end do

    return

  end subroutine r2c_1m_z

  ! c2r transform, multiple 1D FFTs in x direction
  subroutine c2r_1m_x(input, output)

    implicit none

    complex(mytype), dimension(:,:,:), intent(IN)  ::  input
    real(mytype), dimension(:,:,:), intent(OUT) :: output

    integer :: i,j,k, d1,d2,d3

    d1 = size(output,1)
    d2 = size(output,2)
    d3 = size(output,3)

    do k=1,d3
       do j=1,d2
          ! Glassman's FFT is c2c only, 
          ! needing some pre- and post-processing for c2r
          do i=1,d1/2+1
             buf(i) = input(i,j,k)
          end do
          ! expanding to a full-size complex array
          ! For odd N, the storage is:
          !  1, 2, ...... N/2+1   integer division rounded down
          !     N, ...... N/2+2   => a(i) is conjugate of a(N+2-i)
          ! For even N, the storage is:
          !  1, 2, ...... N/2  , N/2+1
          !     N, ...... N/2+2  again a(i) conjugate of a(N+2-i)
          do i=d1/2+2,d1
             buf(i) =  conjg(buf(d1+2-i))
          end do
          call spcfft(buf,d1,1,scratch)
          do i=1,d1
             ! simply drop imaginary part
             output(i,j,k) = real(buf(i), kind=mytype)
          end do
       end do
    end do

    return

  end subroutine c2r_1m_x

  ! c2r transform, multiple 1D FFTs in z direction
  subroutine c2r_1m_z(input, output)

    implicit none

    complex(mytype), dimension(:,:,:), intent(IN)  ::  input
    real(mytype), dimension(:,:,:), intent(OUT) :: output

    integer :: i,j,k, d1,d2,d3

    d1 = size(output,1)
    d2 = size(output,2)
    d3 = size(output,3)

    do j=1,d2
       do i=1,d1
          do k=1,d3/2+1
             buf(k) = input(i,j,k)
          end do
          do k=d3/2+2,d3
             buf(k) =  conjg(buf(d3+2-k))
          end do
          call spcfft(buf,d3,1,scratch)
          do k=1,d3
             output(i,j,k) = real(buf(k), kind=mytype)
          end do
       end do
    end do

    return

  end subroutine c2r_1m_z


#include "fft_common_3d.f90"

  
end module decomp_2d_fft
