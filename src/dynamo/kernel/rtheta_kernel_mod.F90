!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief The kernel computes the rhs of the thermodynamic equation for the nonlinear equations, 
!>         this constists entirely of the advection term u.grad(theta)
!> @detail Kernel to  compute the rhs of thermodynamic equation for the nonlinear equations, in 
!>         the absense of source terms this is purely an advection term:
!>         rtheta = -u.grad(theta)
module rtheta_kernel_mod
use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, func_type,                     &
                                    GH_FIELD, GH_READ, GH_INC,               &
                                    W0, W2, W3,                              &
                                    GH_BASIS, GH_DIFF_BASIS, GH_ORIENTATION, &
                                    CELLS
use constants_mod,           only : r_def

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: rtheta_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                                  &
       arg_type(GH_FIELD,   GH_INC,  W0),                              &
       arg_type(GH_FIELD,   GH_READ, W0),                              &       
       arg_type(GH_FIELD,   GH_READ, W2),                              &
       arg_type(GH_FIELD,   GH_READ, W3)                               &
       /)
  type(func_type) :: meta_funcs(3) = (/                                &
       func_type(W0, GH_BASIS, GH_DIFF_BASIS),                         &
       func_type(W2, GH_BASIS, GH_ORIENTATION),                        &
       func_type(W3, GH_BASIS)                                         &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass ::rtheta_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface rtheta_kernel_type
   module procedure rtheta_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public rtheta_code
contains

type(rtheta_kernel_type) function rtheta_kernel_constructor() result(self)
  return
end function rtheta_kernel_constructor

!> @brief The subroutine which is called directly by the Psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf_w0 The number of degrees of freedom per cell for w0
!! @param[in] undf_w0  The number of unique degrees of freedom  for w0
!! @param[in] map_w0 Integer array holding the dofmap for the cell at the base of the column for w0
!! @param[in] w0_basis Real 5-dim array holding basis functions evaluated at gaussian quadrature points 
!! @param[inout] r_theta Real array the data 
!! @param[inout] theta Real array the potential temperature
!! @param[in] ndf_w2 The number of degrees of freedom per cell for w2
!! @param[in] undf_w2  The number of unique degrees of freedom  for w2
!! @param[in] map_w2 Integer array holding the dofmap for the cell at the base of the column for w2 
!! @param[in] w2_basis Real 5-dim array holding basis functions evaluated at gaussian quadrature points 
!! @param[in] orientation the orientation arrays for the velocity field
!! @param[in] f the mass flux field
!! @param[in] ndf_w3 The number of degrees of freedom per cell for w3
!! @param[in] undf_w3  The number of unique degrees of freedom  for w3
!! @param[in] map_w3 Integer array holding the dofmap for the cell at the base of the column for w3 
!! @param[in] w3_basis Real 5-dim array holding basis functions evaluated at gaussian quadrature points 
!! @param[in] rho the density field
!! @param[in] nqp_h the number of horizontal quadrature points
!! @param[in] nqp_v the number of vertical quadrature points
!! @param[in] wqp_h the weights of the horizontal quadrature points
!! @param[in] wqp_v the weights of the vertical quadrature points
subroutine rtheta_code(nlayers,                                                &
                       ndf_w0, undf_w0, map_w0, w0_basis, w0_diff_basis,       &
                       r_theta, theta,                                         &
                       ndf_w2, undf_w2, map_w2, w2_basis, orientation, f,      &
                       ndf_w3, undf_w3, map_w3, w3_basis, rho,                 &
                       nqp_h, nqp_v, wqp_h, wqp_v )
                               
  
  !Arguments
  integer, intent(in) :: nlayers, nqp_h, nqp_v
  integer, intent(in) :: ndf_w0, ndf_w2, ndf_w3, undf_w0, undf_w2, undf_w3

  integer, dimension(ndf_w0), intent(in) :: map_w0
  integer, dimension(ndf_w2), intent(in) :: map_w2, orientation
  integer, dimension(ndf_w3), intent(in) :: map_w3

  real(kind=r_def), dimension(1,ndf_w0,nqp_h,nqp_v), intent(in) :: w0_basis  
  real(kind=r_def), dimension(3,ndf_w0,nqp_h,nqp_v), intent(in) :: w0_diff_basis  
  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v), intent(in) :: w2_basis 
  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v), intent(in) :: w3_basis  

  real(kind=r_def), dimension(undf_w0), intent(inout) :: r_theta
  real(kind=r_def), dimension(undf_w0), intent(in)    :: theta
  real(kind=r_def), dimension(undf_w2), intent(in)    :: f
  real(kind=r_def), dimension(undf_w3), intent(in)    :: rho

  real(kind=r_def), dimension(nqp_h), intent(in)      ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

  !Internal variables
  integer               :: df, k 
  integer               :: qp1, qp2
  
  real(kind=r_def), dimension(ndf_w0)          :: rtheta_e, theta_e
  real(kind=r_def), dimension(ndf_w2)          :: f_e
  real(kind=r_def), dimension(ndf_w3)          :: rho_e
  real(kind=r_def) :: f_at_quad(3), grad_theta_at_quad(3) 
  real(kind=r_def) :: rho_at_quad
  real(kind=r_def) :: advective_term             
  
  do k = 0, nlayers-1
  ! Extract element arrays of chi
    do df = 1, ndf_w0
      rtheta_e(df) = 0.0_r_def
      theta_e(df)  = theta(  map_w0(df) + k )
    end do
    do df = 1, ndf_w2
      f_e(df) = f( map_w2(df) + k )*real(orientation(df),r_def)
    end do
    do df = 1, ndf_w3
      rho_e(df) = rho( map_w3(df) + k )
    end do   
  ! compute the RHS integrated over one cell
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        f_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w2
          f_at_quad(:)  = f_at_quad(:)  + f_e(df)*w2_basis(:,df,qp1,qp2)
        end do
        grad_theta_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w0
          grad_theta_at_quad(:) = grad_theta_at_quad(:) &
                                + theta_e(df)*w0_diff_basis(:,df,qp1,qp2)
        end do
        rho_at_quad = 0.0_r_def
        do df = 1, ndf_w3
          rho_at_quad  = rho_at_quad  + rho_e(df)*w3_basis(1,df,qp1,qp2)
        end do

        advective_term = wqp_h(qp1)*wqp_v(qp2) &
                       * dot_product(f_at_quad,grad_theta_at_quad)/rho_at_quad

        do df = 1, ndf_w0
          rtheta_e(df) = rtheta_e(df) - w0_basis(1,df,qp1,qp2)*advective_term
        end do
      end do
    end do
    do df = 1, ndf_w0
      r_theta( map_w0(df) + k ) =  r_theta( map_w0(df) + k ) + rtheta_e(df)
    end do 
  end do
  
end subroutine rtheta_code

end module rtheta_kernel_mod
