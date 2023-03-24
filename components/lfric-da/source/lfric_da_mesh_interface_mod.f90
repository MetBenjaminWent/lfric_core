!------------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------

!> @brief  Container for methods that retrieve mesh shape data for use by DA

module lfric_da_mesh_interface_mod
  use mesh_mod,                       only: mesh_type
  use local_mesh_mod,                 only: local_mesh_type
  use mesh_collection_mod,            only: mesh_collection
  use constants_mod,                  only: i_def, r_def
  use fs_continuity_mod,              only: Wtheta, W3
  use extrusion_config_mod,           only: stretching_height
  use coord_transform_mod,            only: xyz2llr
  use extrusion_mod,                  only: TWOD

  implicit none

contains

  !> @brief  Gets the number of vertical layers of the given mesh
  !>
  !> @param[in]  mesh_id  The id of the mesh
  !> @return  The number of layers of the mesh
  function get_levels(mesh_id) result(levels)
    implicit none
    integer(i_def), intent(in)  :: mesh_id
    integer(i_def) :: levels

    type(mesh_type), pointer :: mesh

    mesh => mesh_collection%get_mesh(mesh_id)
    levels = mesh%get_nlayers()

  end function get_levels

  !> @brief  Determines the resolution of a cube sphere mesh
  !>
  !> @param[in]   mesh_id        The id of the mesh
  !> @param[out]  is_cubesphere  Returns true if the mesh is determined to be a
  !>                             cubesphere.
  !> @param[out]  grid_size      The horizontal resolution of the cubesphere.
  !>                             If the mesh is not a cubesphere this is instead
  !>                             the total number of cells.
  subroutine get_grid_size(mesh_id, is_cubesphere, grid_size)
    implicit none
    integer(i_def), intent(in)  :: mesh_id
    logical, intent(out)        :: is_cubesphere
    integer(i_def), intent(out) :: grid_size

    type(mesh_type), pointer  :: mesh
    type(local_mesh_type)     :: local_mesh
    integer(i_def)            :: ncells

    mesh => mesh_collection%get_mesh(mesh_id)
    local_mesh = mesh%get_local_mesh()
    ncells = local_mesh%get_ncells_global_mesh()

    ! get_ncells_global_mesh

    ! No metadata exists to describe overall mesh shape
    ! A periodic, spherical mesh with 6n^2 cells is likely to be a cubesphere
    grid_size = nint( sqrt( real(ncells, kind=r_def) / 6 ), kind = i_def )

    is_cubesphere = mesh%is_topology_periodic()               &
              .and. mesh%is_geometry_spherical()              &
              .and. grid_size ** 2 * 6 == ncells

    ! If the mesh is not a cubesphere, return number of cells instead
    if ( .not. is_cubesphere ) then
      grid_size = ncells
    end if

  end subroutine get_grid_size

  !> @brief  Gets the cell-centred lon/lat coordinates for the given mesh
  !> @details  The atlas_field_interface requires a map_horizontal_ptr to
  !>           identify cells in atlas and LFRic. JEDI constructs this based on
  !>           the lonlat coordinates of the cells, retrieved here, and so this
  !>           procedure breaks field encapsulation.
  !>
  !> @param[in]   mesh_id       The id of the mesh
  !> @param[out]  lon_external  Array of longitude points in default order
  !> @param[out]  lat_external  Array of latitude points in default order
  subroutine get_lonlat(mesh_id, lon_external, lat_external)
    implicit none
    integer(i_def), intent(in)              :: mesh_id
    real(r_def), intent(inout), allocatable :: lon_external(:), lat_external(:)

    type(mesh_type), pointer  :: mesh
    integer(i_def)            :: i, ncells
    real(r_def)               :: cell_centre(3), radius

    mesh => mesh_collection%get_mesh(mesh_id)
    mesh => mesh_collection%get_mesh_variant(mesh, TWOD)

    ncells = mesh%get_last_edge_cell()

    allocate(lon_external(ncells))
    allocate(lat_external(ncells))

    do i=1, ncells
      ! Get coords of cell centres
      call mesh%get_cell_centre_coords( i, cell_centre(:) )
      call xyz2llr( cell_centre(1),  &
                    cell_centre(2),  &
                    cell_centre(3),  &
                    lon_external(i), &
                    lat_external(i), &
                    radius )
    end do

  end subroutine get_lonlat

  !> @brief  Get the normalised heights of W3 levels (cell centres)
  !>
  !> @param[in]  mesh_id  The id of the mesh
  !> @return  An array containing the normalised heights
  function get_sigma_w3_levels(mesh_id) result(levels)
    implicit none
    integer(i_def), intent(in) :: mesh_id
    real(r_def), allocatable :: levels(:)

    integer(i_def) :: len
    real(r_def), allocatable :: wtheta_levels(:)
    type(mesh_type), pointer :: mesh

    mesh => mesh_collection%get_mesh(mesh_id)

    allocate( wtheta_levels( mesh%get_nlayers() + 1 ) )
    call mesh%get_eta(wtheta_levels)

    len = size(wtheta_levels) - 1
    levels = ( wtheta_levels(2:len+1) + wtheta_levels(1:len) ) / 2

  end function get_sigma_w3_levels

  !> @brief  Get the normalised heights of Wtheta levels (cell edges)
  !>
  !> @param[in]  mesh_id  The id of the mesh
  !> @return  An array containing the normalised heights
  function get_sigma_wtheta_levels(mesh_id) result(levels)
    implicit none
    integer(i_def), intent(in) :: mesh_id
    real(r_def), allocatable :: levels(:)

    type(mesh_type), pointer :: mesh

    mesh => mesh_collection%get_mesh(mesh_id)

    allocate( levels( mesh%get_nlayers() + 1 ) )
    call mesh%get_eta(levels)

  end function get_sigma_wtheta_levels

  !> @brief Get the physical height above which mesh levels are not affected by
  !>        orography.
  !> @details  Also called constant level height.
  !>           This value is stored only in configuration, not the mesh object.
  !>
  !> @return  stretching_height in physical coordinates
  function get_stretching_height() result(stretching_height_out)
    implicit none
    real(r_def) :: stretching_height_out

    stretching_height_out = stretching_height
  end function get_stretching_height

  !> @brief  Get the physical height of the top of the mesh.
  !> @details  Also called boundary layer height.
  !>
  !> @param[in]  mesh_id  The id of the mesh
  !> @return  domain_top in physical coordinates
  function get_domain_top(mesh_id) result(domain_top)
    implicit none
    integer(i_def)  :: mesh_id
    real(r_def)     :: domain_top

    type(mesh_type), pointer :: mesh

    mesh => mesh_collection%get_mesh(mesh_id)
    domain_top = mesh%get_domain_top()

  end function get_domain_top

end module lfric_da_mesh_interface_mod