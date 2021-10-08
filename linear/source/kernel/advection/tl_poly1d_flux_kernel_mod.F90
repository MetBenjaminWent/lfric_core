!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Tangent linear for horizontal fluxes with 1D reconstruction.
!> @details The nonlinear model is:
!>          \f F = \rho u \f
!>          The tangent linear model is:
!>          \f F = \rho ls_u + ls_\rho u \f
module tl_poly1d_flux_kernel_mod

use argument_mod,      only : arg_type, func_type,         &
                              reference_element_data_type, &
                              GH_FIELD, GH_SCALAR,         &
                              GH_REAL, GH_INTEGER,         &
                              GH_INC, GH_READ,             &
                              STENCIL, CROSS, GH_BASIS,    &
                              CELL_COLUMN, GH_EVALUATOR,   &
                              ANY_DISCONTINUOUS_SPACE_1,   &
                              outward_normals_to_horizontal_faces
use constants_mod,     only : r_def, i_def
use fs_continuity_mod, only : W2, W3
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: tl_poly1d_flux_kernel_type
  private
  type(arg_type) :: meta_args(8) = (/                                        &
       arg_type(GH_FIELD,  GH_REAL,    GH_INC,   W2),                        & ! Flux
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2),                        & ! Wind
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),                        & ! Density
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2),                        & ! LS Wind
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3, STENCIL(CROSS)),        & ! LS Density
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! coeffs
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                             & ! order
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                              & ! ndata
       /)
  type(func_type) :: meta_funcs(1) = (/                                   &
       func_type(W2, GH_BASIS)                                            &
       /)
  type(reference_element_data_type) :: meta_reference_element(1) = (/     &
       reference_element_data_type( outward_normals_to_horizontal_faces ) &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: tl_poly1d_flux_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: tl_poly1d_flux_code

contains

!> @brief Computes the horizontal fluxes for a tracer density.
!! @param[in]  nlayers      Number of layers
!! @param[in,out] flux      ACTIVE Mass flux field to compute
!! @param[in]  wind         ACTIVE Change in Wind field
!! @param[in]  density      ACTIVE Change in Tracer density
!! @param[in]  ls_wind      PASSIVE Lin state for Wind field
!! @param[in]  ls_density   PASSIVE Lin state for Tracer density
!! @param[in]  stencil_size Size of the stencil (number of cells)
!! @param[in]  stencil_map  Dofmaps for the stencil
!! @param[in]  coeff        Array of polynomial coefficients for interpolation
!! @param[in]  order        Desired order of polynomial reconstruction
!! @param[in]  ndata        Number of data points per dof location
!! @param[in]  ndf_w2       Number of degrees of freedom per cell
!! @param[in]  undf_w2      Number of unique degrees of freedom for the flux
!!                          wind fields
!! @param[in]  map_w2       Dofmap for the cell at the base of the column
!! @param[in]  basis_w2     Basis function array evaluated at w2 nodes
!! @param[in]  ndf_w3       Number of degrees of freedom per cell
!! @param[in]  undf_w3      Number of unique degrees of freedom for the density field
!! @param[in]  map_w3       Dofmap for the cell at the base of the column for the density field
!! @param[in]  ndf_c        Number of degrees of freedom per cell for the coeff space
!! @param[in]  undf_c       Total number of degrees of freedom for the coeff space
!! @param[in]  map_c        Dofmap for the coeff space
!! @param[in]  nfaces_re_h  Number of horizontal neighbours
!! @param[in]  outward_normals_to_horizontal_faces Vector of normals to the
!!                                                 reference element horizontal
!!                                                 "outward faces"
subroutine tl_poly1d_flux_code( nlayers,              &
                                flux,                 &
                                wind,                 &
                                density,              &
                                ls_wind,              &
                                ls_density,           &
                                stencil_size,         &
                                stencil_map,          &
                                coeff,                &
                                order,                &
                                ndata,                &
                                ndf_w2,               &
                                undf_w2,              &
                                map_w2,               &
                                basis_w2,             &
                                ndf_w3,               &
                                undf_w3,              &
                                map_w3,               &
                                ndf_c,                &
                                undf_c,               &
                                map_c,                &
                                nfaces_re_h,          &
                                outward_normals_to_horizontal_faces )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                    :: nlayers
  integer(kind=i_def), intent(in)                    :: ndf_w3
  integer(kind=i_def), intent(in)                    :: undf_w3
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
  integer(kind=i_def), intent(in)                    :: ndf_w2
  integer(kind=i_def), intent(in)                    :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), intent(in)                    :: ndf_c
  integer(kind=i_def), intent(in)                    :: undf_c
  integer(kind=i_def), dimension(ndf_c),  intent(in) :: map_c
  integer(kind=i_def), intent(in)                    :: ndata
  integer(kind=i_def), intent(in)                    :: order
  integer(kind=i_def), intent(in)                    :: stencil_size
  integer(kind=i_def), intent(in)                    :: nfaces_re_h

  real(kind=r_def), dimension(undf_w2), intent(inout) :: flux
  real(kind=r_def), dimension(undf_w2), intent(in)    :: wind
  real(kind=r_def), dimension(undf_w3), intent(in)    :: density
  real(kind=r_def), dimension(undf_w2), intent(in)    :: ls_wind
  real(kind=r_def), dimension(undf_w3), intent(in)    :: ls_density
  real(kind=r_def), dimension(undf_c),  intent(in)    :: coeff

  real(kind=r_def), dimension(3,ndf_w2,ndf_w2), intent(in) :: basis_w2

  integer(kind=i_def), dimension(ndf_w3,stencil_size), intent(in) :: stencil_map

  real(kind=r_def), intent(in) :: outward_normals_to_horizontal_faces(3,nfaces_re_h)

  ! Internal variables
  integer(kind=i_def)                      :: k, df, p, face, stencil, &
                                              stencil_depth, depth, face_mod, &
                                              ijkp
  real(kind=r_def)                         :: direction
  real(kind=r_def), dimension(nfaces_re_h) :: v_dot_n
  real(kind=r_def)                         :: polynomial_density
  real(kind=r_def)                         :: ls_polynomial_density

  integer(kind=i_def), dimension(order+1,nfaces_re_h) :: map1d

  ! Compute 1d map from the cross stencil
  ! i.e for order = 2 the stencil map is
  !      | 5 |
  !  | 2 | 1 | 4 |
  !      | 3 |
  ! so map1d is
  ! ( 1, 2, 4 )
  ! ( 1, 3, 5 )
  ! ( 1, 2, 4 )
  ! ( 1, 3, 5 )
  ! First cell is always the centre cell
  stencil_depth = order/2
  map1d(1,:) = 1
  do face = 1,nfaces_re_h
    depth=1
    face_mod = mod(face+1,2) * stencil_depth
    do stencil = 2,stencil_depth+1
      map1d(stencil+depth-1, face) = stencil + face_mod
      map1d(stencil+depth, face) = stencil + order + face_mod
      depth=depth+1
    end do
  end do

  do df = 1,nfaces_re_h
    v_dot_n(df) =  dot_product(basis_w2(:,df,df),outward_normals_to_horizontal_faces(:,df))
  end do

  ! Horizontal flux computation
  do k = 0, nlayers - 1
    do df = 1,nfaces_re_h

      ! Check if this is the upwind cell
      direction = ls_wind(map_w2(df) + k )*v_dot_n(df)
      if ( direction > 0.0_r_def ) then

        ! Linearisation state
        ls_polynomial_density = 0.0_r_def
        do p = 1,order+1
          stencil = map1d(p,df)
          ijkp = p - 1 + (df-1)*(order+1) + k*ndata + map_c(1)
          ls_polynomial_density = ls_polynomial_density &
                                + ls_density( stencil_map(1,stencil) + k )*coeff( ijkp )
        end do

        ! Perturbation
        polynomial_density = 0.0_r_def
        do p = 1,order+1
          stencil = map1d(p,df)
          ijkp = p - 1 + (df-1)*(order+1) + k*ndata + map_c(1)
          polynomial_density = polynomial_density &
                             + density( stencil_map(1,stencil) + k )*coeff( ijkp )
        end do

        flux(map_w2(df) + k ) = wind(map_w2(df) + k) * ls_polynomial_density + &
                                ls_wind(map_w2(df) + k) * polynomial_density
      end if
    end do
  end do

end subroutine tl_poly1d_flux_code

end module tl_poly1d_flux_kernel_mod
