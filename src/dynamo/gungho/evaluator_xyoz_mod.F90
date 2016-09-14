!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Contains evaluator_xyoz_type and quadrature_xyoz_type

!> @details This module contains the xyoz evaluator and quadrature types
!> (evaluator_xyoz_type and quadrature_xyoz_type). The evaluator contains 
!> points and quadrature extends this (inheritance) to include weights. The
!> points and weights points and weights are stored in 2D (x-y - horizontal) 
!> and 1D (z - vertical). A proxy is used to access the data and the 
!> construction of each type is defined below. A type bound procedure 
!> 'compute_evaluate' is also available. This method uses the evaluate_function 
!> defined in objects of class evaluate_function_type (e.g. function space) for 
!> the xyoz data points.
!> 
!> ~Constructors~
!> 
!> evaluator_xyoz_type(np_xy, np_z, points_xy, points_z)
!> provide the number of points (np_xy and np_z) and the points (points_x(2, np_xy) 
!> and points_z(np_z)) in horizontal and vertical directions.
!> The constructor copies the points provided into the evaluator's arrays
!> 
!> quadrature_xyoz_type(np_xy, np_z, rule)
!> provide the number of points (np_xy and np_z) in horizontal and vertical 
!> directions and the quadrature rule (rule) of type quadrature_rule_type. The
!> quadrature rule is used to generate 1D arrays that populate the quadrature points.
!> There is an assumption that the space is symmetrical so the 1D arrays are computted 
!> using np=np_xy**(1./2.) 
!> Currently there are two supported rules: Gaussian and Newton-Cotes

module evaluator_xyoz_mod
use constants_mod,         only: r_def, i_def, PI, EPS
use log_mod,               only: LOG_LEVEL_ERROR, log_event, log_scratch_space
use quadrature_rule_mod,   only: quadrature_rule_type
use evaluator_mod,         only: evaluator_type
use function_space_mod,    only: function_space_type
use evaluate_function_mod, only: evaluate_function_type

implicit none
private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! xyoz evaluator type
!-------------------------------------------------------------------------------

type, public, extends(evaluator_type) :: evaluator_xyoz_type
  private
  !> allocatable arrays which holds the points (x)
  real(kind=r_def), allocatable :: points_xy(:,:), points_z(:)

  !> Number of quadrature points xy direction = np_x*np_y
  integer(kind=i_def) :: np_xy

  !> Number of quadrature points in each direction
  integer(kind=i_def) :: np_z

contains

  ! Get a proxy with public pointers to the data in a evaluator_xyoz type.
  procedure, public :: get_evaluator_proxy

  ! Evaluates the function for given set of 3d points 
  procedure, public :: compute_evaluate

  ! Destroy the evaluator object
  final     :: evaluator_destructor

end type evaluator_xyoz_type

!> Psy layer representation of a evaluator_xyoz type
!>
!> This is an accessor class that allows access to evaluator_xyoz_type 
!> data and information with each element accessed via a public pointer.
!>
type, public :: evaluator_xyoz_proxy_type

  private
  !> allocatable arrays which holds the values of the gaussian quadrature
  real(kind=r_def), pointer, public :: points_xy(:,:) => null()
  real(kind=r_def), pointer, public :: points_z(:)    => null()
  !> Number of points
  integer, public                   :: np_xy, np_z

contains
end type evaluator_xyoz_proxy_type

!-------------------------------------------------------------------------------
! xyoz quadrature type
!-------------------------------------------------------------------------------

type, public, extends(evaluator_xyoz_type) :: quadrature_xyoz_type
  private
  !> allocatable arrays which holds the quadrature weights
  real(kind=r_def), allocatable :: weights_xy(:), weights_z(:)

contains

  ! Get a proxy with public pointers to the data in a quadrature_xyoz type.
  procedure, public :: get_quadrature_proxy

  ! Destroy the quadrature object
  final     :: quadrature_destructor

end type quadrature_xyoz_type

!> Psy layer representation of a quadrature_xyoz type
!>
!> This is an accessor class that allows access to evaluator_xyoz_type 
!> data and information with each element accessed via a public pointer.
!>
type, public :: quadrature_xyoz_proxy_type

  private
  !> allocatable arrays which holds the values of the gaussian quadrature
  real(kind=r_def), pointer, public :: points_xy(:,:) => null()
  real(kind=r_def), pointer, public :: points_z(:)    => null()
  real(kind=r_def), pointer, public :: weights_xy(:)  => null()
  real(kind=r_def), pointer, public :: weights_z(:)   => null()

  !> Number of points
  integer, public          :: np_xy, np_z

contains
end type quadrature_xyoz_proxy_type

!-------------------------------------------------------------------------------
! Module parameters
!-------------------------------------------------------------------------------
interface evaluator_xyoz_type
  module procedure init_evaluator
end interface

interface quadrature_xyoz_type
  module procedure init_quadrature
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!===============================================================================!
!                             xyoz evaluator type                               !
!===============================================================================!

!-------------------------------------------------------------------------------
!> @brief Initialises the xyoz evaluator type
!> @param[in] np_xy integer, The number of points in the horizontal
!> @param[in] np_z integer, The number of points in the vertical
!> @param[in] points_xy array of points to use (np_xy=size(points_xy))
!> @param[in] points_z array of points to use (np_z=size(points_z))
function init_evaluator(np_xy, np_z, points_xy, points_z) result (self)

  implicit none

  type(evaluator_xyoz_type) :: self
  integer(kind=i_def), intent(in)                  :: np_xy, np_z
  real(kind=r_def), dimension(2,np_xy), intent(in) :: points_xy
  real(kind=r_def), dimension(np_z), intent(in)    :: points_z

  ! Allocate space for the points of points weights in the quad type
  allocate( self%points_z(np_z) )
  allocate( self%points_xy(2, np_xy) )

  ! Copy
  self%np_xy = np_xy
  self%np_z = np_z
  self%points_xy = points_xy
  self%points_z = points_z

end function init_evaluator
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> @brief Function to create a proxy with access to the data in the
!>        evaluator_xyoz_type.
!>
!> @return The proxy type with public pointers to the elements of
!> evaluator_xyoz_type.
type(evaluator_xyoz_proxy_type ) function get_evaluator_proxy(self)

  implicit none

  class(evaluator_xyoz_type), target, intent(in)  :: self

  get_evaluator_proxy % points_xy => self % points_xy
  get_evaluator_proxy % points_z  => self % points_z
  get_evaluator_proxy % np_xy     = self % np_xy
  get_evaluator_proxy % np_z      = self % np_z

end function get_evaluator_proxy
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> @brief Evaluates the a given function for on a set of 3d points
!> @param[in] func_to_call enumerator defining the function to call
!> @param[in] ef object containing the function to evaluate
!> @param[in] ndf integer number of dofs
!> @param[out] basis real 3 dimensional array holding the evaluated
!> function
subroutine compute_evaluate(self, func_to_call, ef, ef_dim, ndf, basis)

  implicit none

  class(evaluator_xyoz_type),                                     intent(in)  :: self
  class(evaluate_function_type),                                  intent(in)  :: ef
  integer(kind=i_def),                                            intent(in)  :: func_to_call
  integer(kind=i_def),                                            intent(in)  :: ef_dim
  integer(kind=i_def),                                            intent(in)  :: ndf
  real(kind=r_def), dimension(ef_dim,ndf,self%np_xy,self%np_z),   intent(out) :: basis

  ! local variables - loop counters
  integer(kind=i_def) :: df
  real(kind=r_def)    :: xyz(3)
  integer(kind=i_def) :: qp1
  integer(kind=i_def) :: qp2

  do qp2 = 1, self%np_z
    xyz(3) = self%points_z(qp2)
    do qp1 = 1, self%np_xy
      xyz(1) = self%points_xy(1,qp1)
      xyz(2) = self%points_xy(2,qp1)
      do df = 1, ndf
          basis(:,df,qp1,qp2) = ef%evaluate_function(func_to_call,df,xyz)
      end do
    end do
  end do


end subroutine compute_evaluate
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> @brief Routine to destroy evaluator
subroutine evaluator_destructor(self)
  implicit none
  type(evaluator_xyoz_type) :: self

  if(allocated(self%points_z))  deallocate(self%points_z)
  if(allocated(self%points_xy)) deallocate(self%points_xy)
  
end subroutine evaluator_destructor
!-------------------------------------------------------------------------------

!===============================================================================!
!                             xyoz quadrature type                              !
!===============================================================================!

!-------------------------------------------------------------------------------
!> @brief Initialises the xyoz quadrature type
!> @param[in] np_xy integer, The number of points in the horizontal
!> @param[in] np_z integer, The number of points in the vertical
!> @param[in] rule quadrature_rule_type, quadrature rule
function init_quadrature(np_xy, np_z, rule) result (self)

  implicit none

  type(quadrature_xyoz_type) :: self
  integer, intent(in) :: np_xy, np_z
  class(quadrature_rule_type), intent(in) :: rule
 
  self%np_xy = np_xy
  self%np_z = np_z

  call create_quadrature( self, rule )

end function init_quadrature
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> @brief Distribute quadrature points and weights
!> @param[in] self The calling quadrature_type
!> @param[in] rule quadrature_rule_type quadrature rule to use
!> @todo this code is correct for quads but will need modification for
!>       hexes/triangles)
subroutine create_quadrature(self, rule)

  implicit none

  class(quadrature_xyoz_type) :: self
  class(quadrature_rule_type), intent(in) :: rule

  integer(kind=i_def)    :: i,j,ic,np_1d
  real(kind=r_def), allocatable       :: points_weights(:,:)

  ! Currently assume that the space is symmetric so that can use the 1-D
  ! quadrature rule
  np_1d = int(self%np_xy ** (1.0/2.0))

  ! Allocate space for the points of points weights in the quad type
  allocate( self%points_z(self%np_z) )
  allocate( self%weights_z(self%np_z) )
  allocate( self%points_xy(2,self%np_xy) )
  allocate( self%weights_xy(self%np_xy) )

  ! Initilise all to zero
  self%points_z(:) = 0.0_r_def
  self%weights_z(:) = 0.0_r_def
  self%points_xy(:,:) = 0.0_r_def
  self%weights_xy(:) = 0.0_r_def

  ! Allocate space for the points and weights of the 1D with dimension defined
  ! in quad type
  allocate( points_weights( np_1d,2 ) )

  ! Get a copy of the 1D points and weights
  points_weights = rule % quadrature_rule( np_1d )

  ! Distribute the 1D points and weights
  ! This is correct for quads (will need modification for hexes/triangles)
  self%points_z = points_weights(:,1)
  self%weights_z = points_weights(:,2)

  ic = 1
  do i=1,np_1d
    do j=1,np_1d
      self%points_xy(1,ic) = points_weights(i,1)
      self%points_xy(2,ic) = points_weights(j,1)
      self%weights_xy(ic) = points_weights(i,2)*points_weights(j,2)

      ic = ic + 1
    end do
  end do

  deallocate( points_weights )

  return
end subroutine create_quadrature
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> @brief Function to create a proxy with access to the data in the
!>        quadrature_xyoz_type.
!>
!> @return The proxy type with public pointers to the elements of
!> quadrature_xyoz_type.
type(quadrature_xyoz_proxy_type ) function get_quadrature_proxy(self)

  implicit none

  class(quadrature_xyoz_type), target, intent(in)  :: self

  get_quadrature_proxy % points_xy  => self % points_xy
  get_quadrature_proxy % points_z   => self % points_z
  get_quadrature_proxy % weights_xy => self % weights_xy
  get_quadrature_proxy % weights_z  => self % weights_z
  get_quadrature_proxy % np_xy      = self % np_xy
  get_quadrature_proxy % np_z       = self % np_z

end function get_quadrature_proxy
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> @brief Routine to destroy quadrature
subroutine quadrature_destructor(self)
  implicit none
  type(quadrature_xyoz_type) :: self

  if(allocated(self%points_z))  deallocate(self%points_z)
  if(allocated(self%points_xy)) deallocate(self%points_xy)
  if(allocated(self%weights_z))  deallocate(self%weights_z)
  if(allocated(self%weights_xy)) deallocate(self%weights_xy)
  
end subroutine quadrature_destructor
!-------------------------------------------------------------------------------

end module evaluator_xyoz_mod
