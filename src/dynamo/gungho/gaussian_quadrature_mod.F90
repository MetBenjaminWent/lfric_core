!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Contains the routines used for Gaussian quadrature.

!> @details This module has a type for the Gaussian quadrature and a static
!> copy of the Gaussian quadrature that is used throughout the model. The first
!> time the Gaussian quadrature is required, it is created and a pointer to it
!> returned. Subsequent times, the pointer to the already created Guassian
!> quadrature is returned.

module gaussian_quadrature_mod
use constants_mod, only: dp, pi, eps
implicit none
private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
type, public :: gaussian_quadrature_type
  private
  !> allocatable arrays which holds the values of the gaussian quadrature
  real(kind=dp), allocatable :: xgp(:), xgp_h(:,:), wgp(:), wgp_h(:)
contains
  !> Function returns a pointer to the Gaussian quadrature. If a Gaussian quadrature
  !> quadrature had not yet been created, it creates one before returning the pointer
  !> to it
  procedure, nopass :: get_instance
  !final     :: final_gauss
  !> Subroutine writes out an answer for a test
  !! @param self the calling gaussian quadrature
  procedure :: test_integrate

  !> Function quassian quadrature integration of a function f 
  !! @param self the calling gp type
  !! @param f real 3D array each of size ngp which holds the sample values of the
  !! function to be integrated
  !! @return real the value of the function thus integrated
  procedure :: integrate

  !> function evaluate 1D basis function of arbitrary order at xk
  !! @param order The order of the basis function
  !! @param ik quadrature point to evaluate basis function at
  !! @param xindex The point at which the function is unity
  !! @param x The grid points
  !! @param bindex The index of the basis function
  procedure :: poly1d

  !> function evaluate 1D basis function of arbitrary order at xk
  !! @param order The order of the basis function
  !! @param ik quadrature point to evaluate basis function at
  !! @param xindex The point at which the function is unity
  !! @param x The grid points
  !! @param bindex The index of the basis function
  procedure :: poly1d_deriv  

end type

!-------------------------------------------------------------------------------
! Module parameters
!-------------------------------------------------------------------------------
!> integer The number of gaussian quadrature points in the vertical
integer, public, parameter      :: ngp_v = 3
!> integer The number of gaussian quadrature points in the horizontal
!! nqp_h=ngp_v*ngp_v for quads. They can be different (triangles or hexes)
!! but there is no setup code for this
integer, public, parameter      :: ngp_h = 9
!> All fields are integrated onto a fixed Guassian quadrature.
!> This is a static copy of that Gaussian quadrature object 
type(gaussian_quadrature_type), target, allocatable, save :: gq

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

function get_instance() result(instance)
  implicit none

  type(gaussian_quadrature_type), pointer :: instance

  if(.not.allocated(gq)) then
    allocate(gq)
    call init_gauss(gq) 
  end if

  instance => gq

  return
end function get_instance
 
subroutine init_gauss(self)
  !-----------------------------------------------------------------------------
  ! Subroutine to compute the Gaussian points (xgp) and (wgp) wgphts 
  !-----------------------------------------------------------------------------
  implicit none

  class(gaussian_quadrature_type) :: self

  integer             :: i, j, m
  real(kind=dp)       :: p1, p2, p3, pp, z, z1
  
  allocate( self%xgp(ngp_v) )
  allocate( self%wgp(ngp_v) ) 
  allocate( self%xgp_h(ngp_h,2) ) 
  allocate( self%wgp_h(ngp_h) ) 

  z1 = 0.0_dp
  m = (ngp_v + 1) / 2

  !Roots are symmetric in the interval - so only need to find half of them  

  do i = 1, m ! Loop over the desired roots 

    z = cos( pi * (i-0.25_dp) / (ngp_v+0.5_dp) )

    !Starting with the above approximation to the ith root, we enter the main
    !loop of refinement by NEWTON'S method   
    do while ( abs(z-z1) > eps )
      p1 = 1.0_dp
      p2 = 0.0_dp

      !Loop up the recurrence relation to get the Legendre polynomial evaluated
      !at z                 
      do j = 1, ngp_v
        p3 = p2
        p2 = p1
        p1 = ((2.0_dp*j-1.0_dp) * z * p2 - (j-1.0_dp)*p3) / j
      end do

      !p1 is now the desired Legendre polynomial. We next compute pp, its
      !derivative, by a standard relation involving also p2, the polynomial of one
      !lower order.      
      pp = ngp_v*(z*p1-p2)/(z*z-1.0_dp)
      z1 = z
      z = z1 - p1/pp             ! Newton's Method  
    end do

    self%xgp(i) =  - z                  ! Roots will be bewteen -1.0 & 1.0 
    self%xgp(ngp_v+1-i) =  + z            ! and symmetric about the origin  
    self%wgp(i) = 2.0/((1.0-z*z)*pp*pp) ! Compute the wgpht and its       
    self%wgp(ngp_v+1-i) = self%wgp(i)          ! symmetric counterpart         

  end do     ! i loop
      
  !Shift quad points from [-1,1] to [0,1]
  do i=1,ngp_v
    self%xgp(i) = 0.5*(self%xgp(i) + 1.0)
  end do

! This is correct for quads (will need modification for hexes/triangles)
  m = 1
  do i=1,ngp_v
    do j=1,ngp_v 
      self%xgp_h(m,1) = self%xgp(i)
      self%xgp_h(m,2) = self%xgp(j)
      self%wgp_h(m) = self%wgp(i)*self%wgp(j)
      
      m = m + 1
    end do
  end do

  return
end subroutine init_gauss

subroutine test_integrate(self)
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------

  use log_mod, only: log_event, log_scratch_space, LOG_LEVEL_INFO

  implicit none

  class(gaussian_quadrature_type) :: self

  integer          :: i,k
  real(kind=dp)    :: func(ngp_v*ngp_v,ngp_v)
  real(kind=dp)    :: answer

  do i=1,ngp_h
    do k=1,ngp_v
      func(i,k) = self%xgp_h(i,1)*self%xgp_h(i,2)*1.0_dp*1.0_dp
    end do
  end do

  answer = self%integrate(func)
  write( log_scratch_space, '(A,F0.0)') 'int(x^2,x=0..1,y=0..1,z=0..1) = ', &
                                        answer
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  return
end subroutine test_integrate
  
function integrate(self,f)
  !-----------------------------------------------------------------------------
  ! Compute 3D Gaussian integration of function f  
  !-----------------------------------------------------------------------------
  implicit none

  class(gaussian_quadrature_type), intent(in) :: self

  real(kind=dp), intent(in) :: f(ngp_v*ngp_v,ngp_v)
  real(kind=dp)             :: integrate

  integer :: i,k

  integrate = 0.0_dp
  do i=1,ngp_h
    do k=1,ngp_v  
      integrate = integrate + self%wgp_h(i)*self%wgp(k)*f(i,k)
    end do
  end do
  
  integrate = 0.5_dp*0.5_dp*0.5_dp*integrate

  return
end function integrate

!----------------------------------------------------------------------------
! Evaluate 1D basis functions of arbitary order at the points of quadrature
!----------------------------------------------------------------------------
function poly1d(self, order, ik, xindex, x , bindex)
  implicit none
  class(gaussian_quadrature_type), intent(in) :: self
  ! order of the basis function 
  integer,                         intent(in) :: order
  ! quadrature point to evaluate basis function at
  integer,                         intent(in) :: ik
  ! Point function is unity at
  real(kind=dp),                   intent(in) :: xindex
  ! Index of basis function
  integer,                         intent(in) :: bindex
  ! grid points
  real(kind=dp),                   intent(in) :: x(order+bindex)

  real(kind=dp) :: poly1d
  ! internal tempories
  ! loop counters
  integer       :: j
  
  poly1d = 1.0 
  
  do j=1,bindex-1
     poly1d = poly1d*(self%xgp(ik)-x(j))/(xindex-x(j))
  end do
  do j=bindex+1,order+1
     poly1d = poly1d*(self%xgp(ik)-x(j))/(xindex-x(j))
  end do
  
  return
end function poly1d

!-----------------------------------------------------------------------------
! evaluate derivative of 1D basis function of arbitrary order at xk
!-----------------------------------------------------------------------------
function poly1d_deriv(self, order,ik,xindex,x,bindex)
  
  implicit none
  class(gaussian_quadrature_type), intent(in) :: self  
  ! Order of basis function
  integer,                         intent(in) :: order
  ! Index of basis function
  integer,                         intent(in) :: bindex
  ! quadrature point to evaluate basis function at
  integer,                         intent(in) :: ik
  ! Point function is unity at
  real(kind=dp),                   intent(in) :: xindex
  ! grid points
  real(kind=dp),                   intent(in) :: x(order+1)

  ! tempories
  real(kind=dp) :: poly1d_deriv
  real(kind=dp) :: denom,t
  integer       :: k,j


  poly1d_deriv = 0.0
  denom = 1.0

  do j=1,bindex-1
    denom = denom * 1.0/(xindex-x(j))  
  end do
  do j=bindex+1,order+1
    denom = denom * 1.0/(xindex-x(j))  
  end do

  do k=1,bindex-1
    t = 1.0
    do j=1,order+1
      if (j .ne. bindex .and. j .ne. k ) then
        t = t * (self%xgp(ik) - x(j))
      end if
    end do
    poly1d_deriv = poly1d_deriv + t*denom
  end do  
  do k=bindex+1,order+1
    t = 1.0
    do j=1,order+1
      if (j .ne. bindex .and. j .ne. k ) then
        t = t * (self%xgp(ik) - x(j))
      end if
    end do
    poly1d_deriv = poly1d_deriv + t*denom
  end do 

  end function poly1d_deriv

end module gaussian_quadrature_mod
