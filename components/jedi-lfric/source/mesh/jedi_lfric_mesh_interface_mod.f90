!------------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------

!> @brief  Container for methods that retrieve mesh shape data for use by DA

module jedi_lfric_mesh_interface_mod

  use mesh_mod,             only: mesh_type
  use local_mesh_mod,       only: local_mesh_type
  use mesh_collection_mod,  only: mesh_collection
  use constants_mod,        only: i_def, r_def
  use fs_continuity_mod,    only: Wtheta, W3
  use extrusion_config_mod, only: stretching_height
  use coord_transform_mod,  only: xyz2llr
  use extrusion_mod,        only: TWOD
  use log_mod,              only: log_event, LOG_LEVEL_ERROR, log_scratch_space

  implicit none


  type(mesh_type),       pointer :: mesh       => null()
  type(local_mesh_type), pointer :: local_mesh => null()

  private :: mesh, local_mesh

contains

  !> @brief  Sets the target mesh for functions in this module to the mesh with
  !>         the given name.
  !> @details  The mesh name is searched in the global mesh_collection.
  !>
  !> @param[in]  mesh_name  The name of the desired target mesh.
  subroutine set_target_mesh(mesh_name)

    implicit none

    character(*), intent(in)  :: mesh_name

    mesh => mesh_collection%get_mesh(mesh_name)
    local_mesh => mesh%get_local_mesh()

  end subroutine set_target_mesh

  !> @brief  Gets the number of vertical layers in the target mesh.
  !>
  !> @return  The number of layers in the mesh.
  function get_nlayers() result(layers)

    implicit none

    integer(i_def) :: layers

    layers = mesh%get_nlayers()

  end function get_nlayers

  !> @brief  Gets the number of cells in a horizontal layer of the mesh in the
  !>         local partition, excluding ghost and halo cells.
  !>
  !> @return  The number of cells in a horizontal layer in the mesh.
  function get_layer_ncells() result(horizontal)

    implicit none

    integer(i_def) :: horizontal

    horizontal = mesh%get_last_edge_cell()

  end function get_layer_ncells

  !> @brief  Determines whether the target mesh is a cubesphere.
  !> @details  Cubesphere meshes are spherical with six panels.
  !>
  !> @return  Number of cells along the edge of each cubesphere panel.
  function is_mesh_cubesphere() result(is_cubesphere)

    implicit none

    logical :: is_cubesphere

    is_cubesphere = mesh%is_topology_periodic()  &
              .and. mesh%is_geometry_spherical() &
              .and. local_mesh%get_num_panels_global_mesh() == 6_i_def

  end function is_mesh_cubesphere

  !> @brief  Gets the number of cells along the edge of each panel in a
  !>         cubesphere mesh.
  !> @details  Cubesphere panels are identical squares so the single return
  !>           value is the length of all edges.
  !>
  !> @return  Number of cells along the edge of each cubesphere panel.
  function get_cubesphere_resolution() result(grid_size)

    implicit none

    integer(i_def) :: grid_size

    integer(i_def) :: n_cells, n_panels, cells_per_panel

    n_cells  = local_mesh%get_ncells_global_mesh()
    n_panels = local_mesh%get_num_panels_global_mesh()

    cells_per_panel = n_cells / n_panels
    grid_size = nint( sqrt( real( cells_per_panel, kind=r_def ) ), kind=i_def )

  end function get_cubesphere_resolution

  !> @brief    Gets the cell-centred lon/lat coordinates for the given mesh.
  !> @details  This requires the local_mesh to be defined in lon/lat
  !>           coordinates. A cartesian local_mesh will cause this to throw an
  !>           error.
  !>
  !> @return  lonlat 2xN array of longitude/latitude points in default order.
  function get_lonlat() result(lonlat)

    implicit none

    real(r_def), allocatable :: lonlat(:,:)

    integer(i_def) :: i, ncells

    ncells = local_mesh%get_last_edge_cell()

    allocate( lonlat( 2, ncells ) )

    if ( local_mesh%is_coord_sys_ll() ) then

      do i=1, ncells
        call local_mesh%get_cell_coords( i, lonlat(:,i) )
      end do

    else

      call log_event( "Specified local_mesh does not use lonlat coordinates.", &
                      log_level_error )

    end if

  end function get_lonlat

  !> @brief  Get the normalised heights of W3 levels (cell centres)
  !>
  !> @return  An array containing the normalised heights
  function get_sigma_w3_levels() result(levels)

    implicit none

    real(r_def), allocatable :: levels(:)

    integer(i_def) :: len
    real(r_def), allocatable :: wtheta_levels(:)

    wtheta_levels = get_sigma_wtheta_levels()

    len = size( wtheta_levels ) - 1
    levels = ( wtheta_levels(2:len+1) + wtheta_levels(1:len) ) / 2

  end function get_sigma_w3_levels

  !> @brief  Get the normalised heights of Wtheta levels (cell edges)
  !>
  !> @return  An array containing the normalised heights
  function get_sigma_wtheta_levels() result(levels)

    implicit none

    real(r_def), allocatable :: levels(:)

    allocate( levels( mesh%get_nlayers() + 1 ) )
    call mesh%get_eta( levels )

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
  !> @return  domain_top in physical coordinates
  function get_domain_top() result(domain_top)

    implicit none

    real(r_def) :: domain_top

    domain_top = mesh%get_domain_top()

  end function get_domain_top

end module jedi_lfric_mesh_interface_mod
