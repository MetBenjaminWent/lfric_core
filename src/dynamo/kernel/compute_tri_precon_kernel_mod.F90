!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Compute the tridiagonal vertical only terms for the helmholtz matrix 
!!        for lowest order elements 
!> @details Compute the terms from the helmholtz matrix restricted to the
!!          vertical for lowest order and store them in three fields suitable
!!          for use with a triadiagonal solver

module compute_tri_precon_kernel_mod
use constants_mod,           only: r_def, i_def
use kernel_mod,              only: kernel_type
use argument_mod,            only: arg_type, func_type,                      &
                                   GH_OPERATOR, GH_FIELD, GH_READ, GH_WRITE, &
                                   W0, W3, ANY_SPACE_1,                      &
                                   GH_BASIS, GH_DIFF_BASIS,                  &
                                   CELLS
use planet_config_mod,       only : kappa, cp
use timestepping_config_mod, only : dt, alpha

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
type, public, extends(kernel_type) :: compute_tri_precon_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                                  &
       arg_type(GH_FIELD*3,  GH_WRITE, W3),                            &
       arg_type(GH_FIELD,    GH_READ,  W0),                            &
       arg_type(GH_FIELD,    GH_READ,  W3),                            &
       arg_type(GH_FIELD*3,  GH_READ,  ANY_SPACE_1)                    &
       /)
  type(func_type) :: meta_funcs(1) = (/                                &
       func_type(ANY_SPACE_1, GH_DIFF_BASIS)                           &
       /)
  integer :: iterates_over = CELLS

contains
  procedure, nopass :: compute_tri_precon_code
end type compute_tri_precon_kernel_type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface compute_tri_precon_kernel
   module procedure compute_tri_precon_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public compute_tri_precon_code
contains

type(compute_tri_precon_kernel_type) function compute_tri_precon_constructor() result(self)
  return
end function compute_tri_precon_constructor
  
!> @brief Compute the tridiagonal vertical only terms for the helmholtz matrix 
!!        for lowest order elements 
!! @param[in] nlayers Number of layers.
!! @param[out] tri_0 Diagonal of the tridiagonal matrix
!! @param[out] tri_plus Upper diagonal of the tridiagonal matrix
!! @param[out] tri_mins Lower diagonal of the tridiagonal matrix
!! @param[in]  theta Potential temperature field
!! @param[in]  rho Density field
!! @param[in]  chi1 First coordinate component
!! @param[in]  chi2 Second coordinate component
!! @param[in]  chi3 Third coordinate component
!! @param[in]  ndf_w3 Number of degrees of freedom per cell for the operator space.
!! @param[in]  undf_w3 Number of unique degrees of freedum for the w3 space
!! @param[in]  map_w3 Dofmap for the cell at the base of the column.
!! @param[in]  ndf_w0 Number of degrees of freedom per cell for the operator space.
!! @param[in]  undf_w0 Number of unique degrees of freedum for the w0 space
!! @param[in]  map_w0 Dofmap for the cell at the base of the column.
!! @param[in]  diff_basis_chi Differential basis functions evaluated at W3 nodal points.
subroutine compute_tri_precon_code(nlayers,                         &
                                   tri_0, tri_plus, tri_minus,      &
                                   theta, rho,                      &
                                   chi1, chi2, chi3,                &
                                   ndf_w3, undf_w3, map_w3,         &
                                   ndf_w0, undf_w0, map_w0,         &
                                   ndf_chi, undf_chi, map_chi,      &
                                   diff_basis_chi )

  use calc_exner_pointwise_mod, only: calc_exner_pointwise
  use coordinate_jacobian_mod,  only: coordinate_jacobian

  implicit none
  !Arguments
  integer(kind=i_def), intent(in) :: ndf_w3, ndf_w0, ndf_chi
  integer(kind=i_def), intent(in) :: undf_w3, undf_w0, undf_chi
  integer(kind=i_def), intent(in) :: nlayers

  integer, dimension(ndf_chi), intent(in) :: map_chi
  integer, dimension(ndf_w0),  intent(in) :: map_w0
  integer, dimension(ndf_w3),  intent(in) :: map_w3

  real(kind=r_def), dimension(undf_w3),  intent(out) :: tri_0
  real(kind=r_def), dimension(undf_w3),  intent(out) :: tri_plus
  real(kind=r_def), dimension(undf_w3),  intent(out) :: tri_minus
  real(kind=r_def), dimension(undf_w3),  intent(in)  :: rho
  real(kind=r_def), dimension(undf_w0),  intent(in)  :: theta
  real(kind=r_def), dimension(undf_chi), intent(in)  :: chi1
  real(kind=r_def), dimension(undf_chi), intent(in)  :: chi2
  real(kind=r_def), dimension(undf_chi), intent(in)  :: chi3

  real(kind=r_def), dimension(3,ndf_chi,ndf_w3,1), intent(in) :: diff_basis_chi

  !Internal variables
  integer(kind=i_def)                        :: k, kp, km, df
  real(kind=r_def)                           :: theta_ref, dpdz, theta_m, theta_p, rho_p, rho_m
  real(kind=r_def)                           :: kappa_term
  real(kind=r_def), dimension(0:nlayers-1)   :: exner, dj
  real(kind=r_def), dimension(0:nlayers)     :: HB_inv, dthetadz, Pw, Pt
  real(kind=r_def), dimension(3,3,0:nlayers-1) :: jac
  real(kind=r_def), dimension(3,3)           :: jac_av
  real(kind=r_def), dimension(ndf_chi)       :: chi1_e, chi2_e, chi3_e
  real(kind=r_def)                           :: JTJ

  ! Metric terms: 
  ! J(3)^2 on w points
  ! dj on w and rho points
  ! Currently only for lowest order uniform grid with theta in W0
  kappa_term = (1.0_r_def - kappa)/kappa
  ! Compute layer terms
  do k = 0,nlayers-1 
    theta_ref = 0.125*(theta(map_w0(1) + k) + theta(map_w0(2) + k) &
                     + theta(map_w0(3) + k) + theta(map_w0(4) + k) &
                     + theta(map_w0(5) + k) + theta(map_w0(6) + k) &
                     + theta(map_w0(7) + k) + theta(map_w0(8) + k))
    exner(k) = calc_exner_pointwise(rho(map_w3(1)+k), theta_ref)

    do df = 1,ndf_chi
      chi1_e(df) = chi1(map_chi(df)+k)
      chi2_e(df) = chi2(map_chi(df)+k)
      chi3_e(df) = chi3(map_chi(df)+k)
    end do
    call coordinate_jacobian(ndf_chi, 1, 1, chi1_e, chi2_e, chi3_e,  &
                             diff_basis_chi, jac(:,:,k), dj(k))
  end do  

  ! Compute terms on interfaces
  HB_inv(0)   = 1.0
  dthetadz(0) = 0.0
  do k = 1, nlayers-1
    theta_p = 0.25*(theta(map_w0(1)+k+1) + theta(map_w0(2)+k+1) &
                  + theta(map_w0(3)+k+1) + theta(map_w0(4)+k+1) )
    theta_m = 0.25*(theta(map_w0(1)+k-1) + theta(map_w0(2)+k-1) &
                  + theta(map_w0(3)+k-1) + theta(map_w0(4)+k-1) )
    dthetadz(k) = (theta_p - theta_m)
    dpdz = (exner(k) - exner(k-1))
   
    jac_av = 0.5_r_def*(jac(:,:,k-1) + jac(:,:,k))
    JTJ = jac_av(1,3)**2 + jac_av(2,3)**2 + jac_av(3,3)**2
    HB_inv = 1.0_r_def/max(0.1_r_def,JTJ-cp*alpha*dt*dthetadz(k)*dpdz)
  end do
  HB_inv(nlayers)   = 1.0
  dthetadz(nlayers) = 0.0

  Pw(0) = 0.0_r_def
  Pt(0) = 0.0_r_def
  do k = 1, nlayers-1
    theta_ref = 0.25*(theta(map_w0(1) + k) + theta(map_w0(2) + k) &
                    + theta(map_w0(3) + k) + theta(map_w0(4) + k)  )  
    Pw(k) = -alpha*dt*cp*theta_ref*HB_inv(k)
    Pt(k) = -alpha*dt*dthetadz(k)*Pw(k)
  end do
  Pw(nlayers) = 0.0_r_def
  Pt(nlayers) = 0.0_r_def
 
  do k = 0, nlayers - 1
    kp = min(k+1,nlayers-1)
    km = max(k-1,0)
    theta_m = 0.25*(theta(map_w0(1)+k) + theta(map_w0(2)+k) &
                  + theta(map_w0(3)+k) + theta(map_w0(4)+k) )
    theta_p = 0.25*(theta(map_w0(5)+k) + theta(map_w0(6)+k) &
                  + theta(map_w0(7)+k) + theta(map_w0(8)+k) )

    rho_p = 0.5_r_def*(rho(map_w3(1)+k)*dj(k) + rho(map_w3(1)+kp)*dj(kp))
    rho_m = 0.5_r_def*(rho(map_w3(1)+k)*dj(k) + rho(map_w3(1)+km)*dj(km))

    tri_plus(map_w3(1)+k)  = alpha*dt*Pw(k+1)*rho_p/(rho(map_w3(1)+k)*dj(k)) - 0.5_r_def*Pt(k+1)/theta_p
    tri_minus(map_w3(1)+k) = alpha*dt*Pw(k)  *rho_m/(rho(map_w3(1)+k)*dj(k)) + 0.5_r_def*Pt(k)  /theta_m
    tri_0(map_w3(1)+k) = kappa_term/exner(k) - tri_plus(map_w3(1)+k) - tri_minus(map_w3(1)+k)
  end do

end subroutine compute_tri_precon_code

end module compute_tri_precon_kernel_mod
