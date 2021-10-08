!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Tangent linear for computing vertical fluxes through fitting a
!!        high order 1D upwind reconstruction.
module tl_poly1d_vert_adv_kernel_mod

use argument_mod,         only : arg_type, func_type,   &
                                 GH_FIELD, GH_SCALAR,   &
                                 GH_REAL, GH_INTEGER,   &
                                 GH_READWRITE, GH_READ, &
                                 GH_BASIS, CELL_COLUMN, &
                                 GH_EVALUATOR,          &
                                 ANY_DISCONTINUOUS_SPACE_1
use constants_mod,        only : r_def, i_def
use fs_continuity_mod,    only : W2, Wtheta
use kernel_mod,           only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: tl_poly1d_vert_adv_kernel_type
  private
  type(arg_type) :: meta_args(9) = (/                                            &
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, Wtheta),                    &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W2),                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W2),                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                 &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                 &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                                  &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: tl_poly1d_vert_adv_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: tl_poly1d_vert_adv_code

contains

!> @brief Computes the tangent linear for vertical fluxes for a tracer density.
!! @param[in]  nlayers   Number of layers
!! @param[in,out] advective ACTIVE Advective update to increment
!! @param[in]  wind      ACTIVE Change in wind field
!! @param[in]  tracer    ACTIVE Change in tracer field to advect
!! @param[in]  ls_wind   Lin state for wind field
!! @param[in]  ls_tracer Lin state for tracer field to advect
!! @param[in]  coeff     Array of polynomial coefficients for interpolation
!! @param[in]  ndata     Number of data points per dof location
!! @param[in]  global_order Desired order of polynomial reconstruction
!! @param[in]  logspace  If true (=1), then perform interpolation in log space
!! @param[in]  ndf_wt    Number of degrees of freedom per cell
!! @param[in]  undf_wt   Number of unique degrees of freedom for the tracer
!!                       field
!! @param[in]  map_wt    Cell dofmaps for the tracer space
!! @param[in]  ndf_w2    Number of degrees of freedom per cell
!! @param[in]  undf_w2   Number of unique degrees of freedom for the flux &
!!                       wind fields
!! @param[in]  map_w2    Dofmap for the cell at the base of the column
!! @param[in]  ndf_c     Number of degrees of freedom per cell for the coeff
!!                       space
!! @param[in]  undf_c    Total number of degrees of freedom for the coeff space
!! @param[in]  map_c     Dofmap for the coeff space
subroutine tl_poly1d_vert_adv_code( nlayers,              &
                                    advective,            &
                                    wind,                 &
                                    tracer,               &
                                    ls_wind,              &
                                    ls_tracer,            &
                                    coeff,                &
                                    ndata,                &
                                    global_order,         &
                                    logspace,             &
                                    ndf_wt,               &
                                    undf_wt,              &
                                    map_wt,               &
                                    ndf_w2,               &
                                    undf_w2,              &
                                    map_w2,               &
                                    ndf_c,                &
                                    undf_c,               &
                                    map_c)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                    :: nlayers
  integer(kind=i_def), intent(in)                    :: ndf_wt
  integer(kind=i_def), intent(in)                    :: undf_wt
  integer(kind=i_def), intent(in)                    :: ndf_w2
  integer(kind=i_def), intent(in)                    :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt
  integer(kind=i_def), intent(in)                    :: ndf_c
  integer(kind=i_def), intent(in)                    :: undf_c
  integer(kind=i_def), dimension(ndf_c),  intent(in) :: map_c
  integer(kind=i_def), intent(in)                    :: ndata
  integer(kind=i_def), intent(in)                    :: global_order

  real(kind=r_def), dimension(undf_wt), intent(inout) :: advective
  real(kind=r_def), dimension(undf_w2), intent(in)    :: wind
  real(kind=r_def), dimension(undf_wt), intent(in)    :: tracer
  real(kind=r_def), dimension(undf_w2), intent(in)    :: ls_wind
  real(kind=r_def), dimension(undf_wt), intent(in)    :: ls_tracer
  real(kind=r_def), dimension(undf_c),  intent(in)    :: coeff

  integer(kind=i_def), intent(in) :: logspace

  ! Internal variables
  integer(kind=i_def)            :: k, ij, p, stencil, order, &
                                    m, ijkp, direction, boundary_offset, &
                                    vertical_order, km, idx

  real(kind=r_def)               :: w, tracer_p, tracer_m, ls_tracer_p, ls_tracer_m
  integer(kind=i_def), allocatable, dimension(:,:) :: smap

  ! Ensure that we reduce the order if there are only a few layers
  vertical_order = min(global_order, nlayers-1)

  ! Compute the offset map for all even orders up to order
  allocate( smap(global_order+1,0:global_order) )
  smap(:,:) = 0
  do m = 0,global_order
    do stencil = 1,m+1
      smap(stencil,m) = - m/2 + (stencil-1)
    end do
  end do

  ij = map_wt(1)

  ! Reconstruct upwind tracer at cell centres with stencil
  ! direction determined by w at Wtheta points
  do k = 1, nlayers - 1

    ! Check if this is the upwind cell
    w = sign( 1.0_r_def, ls_wind(map_w2(5) + k ) )
    ! w > 0, direction = 0
    ! w < 0, direction = 1
    direction = int( 0.5_r_def*(abs(w)-w),i_def )

    !---- Linearisation state---------------------------------------
    ! Compute value at W3 point in cell below this cell
    km = k - 1
    order = min(vertical_order, min(2*(km+1), 2*(nlayers-1 - (km-1))))
    ! At the boundaries reduce to centred linear interpolation
    if ( km == 0 .and. direction == 0 ) order = 1
    if (logspace == 1_i_def) then
      ls_tracer_m = 1.0_r_def
      do p = 1,order+1
        ijkp = ij + km + smap(p,order) + direction
        idx = p - 1 + direction*(global_order+1) + km*ndata + map_c(1)
        ls_tracer_m = ls_tracer_m  * abs(ls_tracer( ijkp ))**coeff( idx )
      end do
    else
      ls_tracer_m = 0.0_r_def
      do p = 1,order+1
        ijkp = ij + km + smap(p,order) + direction
        idx = p - 1 + direction*(global_order+1) + km*ndata + map_c(1)
        ls_tracer_m = ls_tracer_m  + ls_tracer( ijkp )*coeff( idx )
      end do
    end if

    ! Compute value at W3 point in this cell
    order = min(vertical_order, min(2*(k+1), 2*(nlayers-1 - (k-1))))
    ! At the boundaries reduce to centred linear interpolation
    if ( k == nlayers - 1 .and. direction == 1 ) order = 1

    ! Offset at the top boundary to ensure we use the correct entries
    boundary_offset = 0
    if ( k == nlayers - 1 .and. direction == 1)  boundary_offset = -1

    if (logspace == 1_i_def) then
      ! Interpolate log(tracer)
      ! I.e. polynomial = exp(c_1*log(tracer_1) + c_2*log(tracer_2) + ...)
      !                 = tracer_1**c_1*tracer_2**c_2...
      ls_tracer_p = 1.0_r_def
      do p = 1,order+1
        ijkp = ij + k + smap(p,order) + direction + boundary_offset
        idx = p - 1 + direction*(global_order+1) + k*ndata + map_c(1)
        ls_tracer_p = ls_tracer_p  * abs(ls_tracer( ijkp ))**coeff( idx )
      end do
    else
      ls_tracer_p = 0.0_r_def
      do p = 1,order+1
        ijkp = ij + k + smap(p,order) + direction + boundary_offset
        idx = p - 1 + direction*(global_order+1) + k*ndata + map_c(1)
        ls_tracer_p = ls_tracer_p  + ls_tracer( ijkp )*coeff( idx )
      end do
    end if

    !--- Perturbation------------------------------------------------
    ! Compute value at W3 point in cell below this cell
    km = k - 1
    order = min(vertical_order, min(2*(km+1), 2*(nlayers-1 - (km-1))))
    ! At the boundaries reduce to centred linear interpolation
    if ( km == 0 .and. direction == 0 ) order = 1

    if (logspace == 1_i_def) then
      tracer_m = 0.0_r_def
      do p = 1,order+1
        ijkp = ij + km + smap(p,order) + direction
        idx = p - 1 + direction*(global_order+1) + km*ndata + map_c(1)
        tracer_m = tracer_m +        &
                   coeff( idx ) *    &
                   tracer( ijkp ) /  &
                   ls_tracer( ijkp )
      end do
      tracer_m = tracer_m * ls_tracer_m
    else
      tracer_m = 0.0_r_def
      do p = 1,order+1
        ijkp = ij + km + smap(p,order) + direction
        idx = p - 1 + direction*(global_order+1) + km*ndata + map_c(1)
        tracer_m = tracer_m  + tracer( ijkp )*coeff( idx )
      end do
    end if

    ! Compute value at W3 point in this cell
    order = min(vertical_order, min(2*(k+1), 2*(nlayers-1 - (k-1))))
    ! At the boundaries reduce to centred linear interpolation
    if ( k == nlayers - 1 .and. direction == 1 ) order = 1

    ! Offset at the top boundary to ensure we use the correct entries
    boundary_offset = 0
    if ( k == nlayers - 1 .and. direction == 1)  boundary_offset = -1

    if (logspace == 1_i_def) then
      ! Interpolate log(tracer)
      ! I.e. polynomial = exp(c_1*log(tracer_1) + c_2*log(tracer_2) + ...)
      !                 = tracer_1**c_1*tracer_2**c_2...
      tracer_p = 0.0_r_def
      do p = 1,order+1
        ijkp = ij + k + smap(p,order) + direction + boundary_offset
        idx = p - 1 + direction*(global_order+1) + k*ndata + map_c(1)
        tracer_p = tracer_p +        &
                   coeff( idx ) *    &
                   tracer( ijkp ) /  &
                   ls_tracer( ijkp )
      end do
      tracer_p = tracer_p * ls_tracer_p
    else
      tracer_p = 0.0_r_def
      do p = 1,order+1
        ijkp = ij + k + smap(p,order) + direction + boundary_offset
        idx = p - 1 + direction*(global_order+1) + k*ndata + map_c(1)
        tracer_p = tracer_p  + tracer( ijkp )*coeff( idx )
      end do
    end if

    ! --- Calculation -------------------------------------------------
    ! Compute advective increment
    advective( map_wt(1) + k ) = advective( map_wt(1) + k )   &
                               + ls_wind( map_w2(5) + k )     &
                               * ( tracer_p - tracer_m )      &
                               + wind( map_w2(5) + k )        &
                               * ( ls_tracer_p - ls_tracer_m )
  end do

  deallocate( smap )
end subroutine tl_poly1d_vert_adv_code

end module tl_poly1d_vert_adv_kernel_mod
