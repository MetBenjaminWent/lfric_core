!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Tangent linear for vertical fluxes with 1D reconstruction.
!> @details The nonlinear model is:
!>          \f F = \rho u \f
!>          The tangent linear model is:
!>          \f F = \rho ls_u + ls_\rho u \f
!>          Near the boundaries the order of reconstruction may be reduced
!>          if there are not enough points to compute desired order.
module tl_poly1d_vert_flux_kernel_mod

use argument_mod,      only : arg_type, func_type,         &
                              reference_element_data_type, &
                              GH_FIELD, GH_SCALAR,         &
                              GH_REAL, GH_INTEGER,         &
                              GH_INC, GH_READ, GH_BASIS,   &
                              CELL_COLUMN, GH_EVALUATOR,   &
                              ANY_DISCONTINUOUS_SPACE_1,   &
                              outward_normals_to_vertical_faces
use constants_mod,     only : r_def, i_def
use fs_continuity_mod, only : W2, W3
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: tl_poly1d_vert_flux_kernel_type
  private
  type(arg_type) :: meta_args(9) = (/                                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_INC,   W2),                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2),                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2),                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                             &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                             &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                              &
       /)
  type(func_type) :: meta_funcs(1) = (/                                 &
       func_type(W2, GH_BASIS)                                          &
       /)
  type(reference_element_data_type) :: meta_reference_element(1) = (/   &
       reference_element_data_type( outward_normals_to_vertical_faces ) &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: tl_poly1d_vert_flux_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public tl_poly1d_vert_flux_code

contains

!> @brief Computes the tangent linear for vertical fluxes.
!! @param[in]  nlayers      Number of layers
!! @param[inout] flux       ACTIVE Change in mass flux field
!! @param[in]  wind         ACTIVE Change in wind field
!! @param[in]  density      ACTIVE Change in tracer density
!! @param[in]  ls_wind      Linearisation state wind field
!! @param[in]  ls_density   Linearisation state tracer density
!! @param[in]  coeff        Array of polynomial coefficients for interpolation
!! @param[in]  ndata        Number of data points per dof location
!! @param[in]  global_order Desired order of polynomial reconstruction
!! @param[in]  logspace     If true (=1), perform interpolation in log space
!! @param[in]  ndf_w2       Number of degrees of freedom per cell
!! @param[in]  undf_w2      Number of unique degrees of freedom for the flux &
!!                          wind fields
!! @param[in]  map_w2       Dofmap for the cell at the base of the column
!! @param[in]  basis_w2     Basis function array evaluated at w2 nodes
!! @param[in]  ndf_w3       Number of degrees of freedom per cell
!! @param[in]  undf_w3      Number of unique degrees of freedom for density
!! @param[in]  map_w3       Cell dofmaps for the density space
!! @param[in]  ndf_c        Number of degrees of freedom per cell for the
!!                          coeff space
!! @param[in]  undf_c       Total number of degrees of freedom for the
!!                          coeff space
!! @param[in]  map_c        Dofmap for the coeff space
!! @param[in]  nfaces_re_v  Number of vertical faces (used by PSyclone to size
!!                          coeff array)
!! @param[in]  outward_normals_to_vertical_faces Vector of normals to the
!!                                               reference element vertical
!!                                               "outward faces"
subroutine tl_poly1d_vert_flux_code( nlayers,                           &
                                     flux,                              &
                                     wind,                              &
                                     density,                           &
                                     ls_wind,                           &
                                     ls_density,                        &
                                     coeff,                             &
                                     ndata,                             &
                                     global_order,                      &
                                     logspace,                          &
                                     ndf_w2,                            &
                                     undf_w2,                           &
                                     map_w2,                            &
                                     basis_w2,                          &
                                     ndf_w3,                            &
                                     undf_w3,                           &
                                     map_w3,                            &
                                     ndf_c,                             &
                                     undf_c,                            &
                                     map_c,                             &
                                     nfaces_re_v,                       &
                                     outward_normals_to_vertical_faces )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                    :: nlayers
  integer(kind=i_def), intent(in)                    :: ndata
  integer(kind=i_def), intent(in)                    :: ndf_w3
  integer(kind=i_def), intent(in)                    :: undf_w3
  integer(kind=i_def), intent(in)                    :: ndf_w2
  integer(kind=i_def), intent(in)                    :: undf_w2
  integer(kind=i_def), intent(in)                    :: ndf_c
  integer(kind=i_def), intent(in)                    :: undf_c
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_c),  intent(in) :: map_c
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
  integer(kind=i_def), intent(in)                    :: global_order, nfaces_re_v
  real(kind=r_def), dimension(undf_w2), intent(out)  :: flux
  real(kind=r_def), dimension(undf_w2), intent(in)   :: wind
  real(kind=r_def), dimension(undf_w3), intent(in)   :: density
  real(kind=r_def), dimension(undf_w2), intent(in)   :: ls_wind
  real(kind=r_def), dimension(undf_w3), intent(in)   :: ls_density
  real(kind=r_def), dimension(undf_c),  intent(in)   :: coeff

  real(kind=r_def), dimension(3,ndf_w2,ndf_w2), intent(in) :: basis_w2

  real(kind=r_def), intent(in)  :: outward_normals_to_vertical_faces(:,:)

  integer(kind=i_def), intent(in) :: logspace

  ! Internal variables
  integer(kind=i_def)            :: k, df, ij, p, stencil, order, id, &
                                    m, ijkp, boundary_offset, vertical_order, &
                                    idx
  integer(kind=i_def)            :: vert_offset
  real(kind=r_def)               :: direction
  real(kind=r_def), dimension(2) :: v_dot_n
  real(kind=r_def)               :: polynomial_density
  real(kind=r_def)               :: ls_polynomial_density

  integer(kind=i_def), allocatable, dimension(:,:) :: smap

  ! Constants

  ! Ensure that we reduce the order if there are only a few layers
  vertical_order = min(global_order, nlayers-1)

  vert_offset = 4_i_def

  ! Compute the offset map for all even orders up to order
  allocate( smap(global_order+1,0:global_order) )
  smap(:,:) = 0
  do m = 0,global_order
    do stencil = 1,m+1
      smap(stencil,m) = - m/2 + (stencil-1)
    end do
  end do

  do id = 1,nfaces_re_v
    df = id + vert_offset
    v_dot_n(id) =  dot_product( basis_w2(:,df,df), &
                                outward_normals_to_vertical_faces(:,id) )
  end do

  ij = map_w3(1)

  ! Vertical flux computation
  do k = 0, nlayers - 1

    order = min( vertical_order, min( 2*(k+1), 2*(nlayers-1 - (k-1)) ) )

    boundary_offset = 0
    if ( order > 0 ) then
      if ( k == 0 )           boundary_offset =  1
      if ( k == nlayers - 1 ) boundary_offset = -1
    end if

    do id = 1,nfaces_re_v

      df = id + vert_offset
      ! Check if this is the upwind cell
      direction = ls_wind( map_w2(df) + k ) * v_dot_n(id)
      if ( direction >= 0.0_r_def ) then

        ! Linearisation state
        if (logspace == 1_i_def) then
          ls_polynomial_density = 1.0_r_def
          do p = 1,order+1
            ijkp = ij + k + smap(p,order) + boundary_offset
            idx = p - 1 + (id-1)*(global_order + 1) + k*ndata + map_c(1)
            ls_polynomial_density = ls_polynomial_density *    &
                                    abs(ls_density( ijkp )) ** &
                                    coeff( idx )
          end do
        else
          ls_polynomial_density = 0.0_r_def
          do p = 1,order+1
            ijkp = ij + k + smap(p,order) + boundary_offset
            idx = p - 1 + (id-1)*(global_order + 1) + k*ndata + map_c(1)
            ls_polynomial_density = ls_polynomial_density + &
                                    ls_density( ijkp ) *    &
                                    coeff( idx )
          end do
        end if

        ! Perturbation
        if (logspace == 1_i_def) then
          polynomial_density = 0.0_r_def
          do p = 1,order+1
            ijkp = ij + k + smap(p,order) + boundary_offset
            idx = p - 1 + (id-1)*(global_order + 1) + k*ndata + map_c(1)
            polynomial_density = polynomial_density + &
                                 coeff( idx ) *       &
                                 density ( ijkp ) /   &
                                 ls_density( ijkp )
          end do
          polynomial_density = polynomial_density * ls_polynomial_density
        else
          polynomial_density = 0.0_r_def
          do p = 1,order+1
            ijkp = ij + k + smap(p,order) + boundary_offset
            idx = p - 1 + (id-1)*(global_order + 1) + k*ndata + map_c(1)
            polynomial_density = polynomial_density + &
                                 density( ijkp ) *    &
                                 coeff( idx )
          end do
        end if

        ! Calculation
        flux( map_w2(df) + k ) = wind(map_w2(df) + k) * ls_polynomial_density &
                               + ls_wind(map_w2(df) + k) * polynomial_density
      end if
    end do
  end do

  ! Boundary conditions
  order = min(vertical_order, 2)
  boundary_offset = 0

  k = 0
  id = 1
  if ( order > 0 ) boundary_offset = 1
  df = id + vert_offset

  ! Linearisation state
  if (logspace == 1_i_def) then
    ls_polynomial_density = 1.0_r_def
    do p = 1,order+1
      ijkp = ij + k + smap(p,order) + boundary_offset
      idx = p - 1 + (id-1)*(global_order + 1) + k*ndata + map_c(1)
      ls_polynomial_density = ls_polynomial_density *    &
                              abs(ls_density( ijkp )) ** &
                              coeff( idx )
    end do
  else
    ls_polynomial_density = 0.0_r_def
    do p = 1,order+1
      ijkp = ij + k + smap(p,order) + boundary_offset
      idx = p - 1 + (id-1)*(global_order + 1) + k*ndata + map_c(1)
      ls_polynomial_density = ls_polynomial_density + &
                              ls_density( ijkp )*     &
                              coeff( idx )
    end do
  end if

  ! Perturbation
  if (logspace == 1_i_def) then
    polynomial_density = 0.0_r_def
    do p = 1,order+1
      ijkp = ij + k + smap(p,order) + boundary_offset
      idx = p - 1 + (id-1)*(global_order + 1) + k*ndata + map_c(1)
      polynomial_density = polynomial_density + &
                           coeff( idx ) *       &
                           density ( ijkp ) /   &
                           ls_density( ijkp )
    end do
    polynomial_density = polynomial_density * ls_polynomial_density
  else
    polynomial_density = 0.0_r_def
    do p = 1,order+1
      ijkp = ij + k + smap(p,order) + boundary_offset
      idx = p - 1 + (id-1)*(global_order + 1) + k*ndata + map_c(1)
      polynomial_density = polynomial_density + &
                           density( ijkp ) *    &
                           coeff( idx )
    end do
  end if

  ! Calculation
  flux( map_w2(df) + k ) = 0.0_r_def
  flux( map_w2(df) + k ) = wind(map_w2(df) + k) * ls_polynomial_density + &
                           ls_wind(map_w2(df) + k) * polynomial_density


  k = nlayers-1
  id = 2
  if ( order > 0 ) boundary_offset =  -1
  df = id + vert_offset

  ! Linearisation state
  if (logspace == 1_i_def) then
    ls_polynomial_density = 1.0_r_def
    do p = 1,order+1
      ijkp = ij + k + smap(p,order) + boundary_offset
      idx = p - 1 + (id-1)*(global_order + 1) + k*ndata + map_c(1)
      ls_polynomial_density = ls_polynomial_density *    &
                              abs(ls_density( ijkp )) ** &
                              coeff( idx )
    end do
  else
    ls_polynomial_density = 0.0_r_def
    do p = 1,order+1
      ijkp = ij + k + smap(p,order) + boundary_offset
      idx = p - 1 + (id-1)*(global_order + 1) + k*ndata + map_c(1)
      ls_polynomial_density = ls_polynomial_density + &
                              ls_density( ijkp ) *    &
                              coeff( idx )
    end do
  end if

  ! Perturbation
  if (logspace == 1_i_def) then
    polynomial_density = 0.0_r_def
    do p = 1,order+1
      ijkp = ij + k + smap(p,order) + boundary_offset
      idx = p - 1 + (id-1)*(global_order + 1) + k*ndata + map_c(1)
      polynomial_density = polynomial_density + &
                           coeff( idx ) *       &
                           density( ijkp ) /    &
                           ls_density( ijkp )
    end do
    polynomial_density = polynomial_density * ls_polynomial_density
  else
    polynomial_density = 0.0_r_def
    do p = 1,order+1
      ijkp = ij + k + smap(p,order) + boundary_offset
      idx = p - 1 + (id-1)*(global_order + 1) + k*ndata + map_c(1)
      polynomial_density = polynomial_density + &
                           density( ijkp ) *    &
                           coeff( idx )
    end do
  end if

  ! Calculation
  flux( map_w2(df) + k ) = 0.0_r_def
  flux( map_w2(df) + k ) = wind(map_w2(df) + k) * ls_polynomial_density + &
                           ls_wind(map_w2(df) + k) * polynomial_density

  deallocate( smap )

end subroutine tl_poly1d_vert_flux_code

end module tl_poly1d_vert_flux_kernel_mod
