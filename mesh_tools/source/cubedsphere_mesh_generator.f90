!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @mainpage Cubedsphere mesh generator
!> @brief   Utility to generate a cubedsphere surface mesh and write to a file
!>          which conforms to the UGRID format convention.
!> @details Usage:
!>
!>          cubedsphere_mesh_generator <filename>
!>          filename - Controlling namelist file
!>
!-----------------------------------------------------------------------------
program cubedsphere_mesh_generator

  use cli_mod,         only: get_initial_filename
  use constants_mod,   only: i_def, imdi
  use cubedsphere_mesh_generator_config_mod,                                  &
                       only: read_cubedsphere_mesh_generator_namelist,        &
                             postprocess_cubedsphere_mesh_generator_namelist, &
                             edge_cells, smooth_passes, nmeshes, mesh_names,  &
                             mesh_filename
  use ESMF
  use gencube_ps_mod,  only: gencube_ps_type
  use io_utility_mod,  only: open_file, close_file
  use log_mod,         only: log_scratch_space, log_event, log_set_level,     &
                             LOG_LEVEL_INFO, LOG_LEVEL_ERROR, LOG_LEVEL_DEBUG
  use ncdf_quad_mod,   only: ncdf_quad_type
  use ugrid_file_mod,  only: ugrid_file_type
  use ugrid_2d_mod,    only: ugrid_2d_type

  implicit none

  type(ESMF_VM)  :: vm
  integer(i_def) :: rc

  character(:), allocatable :: filename
  integer(i_def)            :: namelist_unit

  type(gencube_ps_type),  allocatable :: csgen(:)
  type(ugrid_2d_type),    allocatable :: ugrid_2d(:)
  class(ugrid_file_type), allocatable :: ugrid_file

  integer(i_def) :: fsize
  integer(i_def) :: i

  call log_set_level(LOG_LEVEL_INFO)

  ! Start up ESMF
  CALL ESMF_Initialize( vm=vm, rc=rc,                    &
                        logkindflag=ESMF_LOGKIND_SINGLE, &
                        defaultlogfilename="cubedsphere.log" )
  if (rc /= ESMF_SUCCESS) call log_event( 'Failed to initialise ESMF.', &
                                          LOG_LEVEL_ERROR )

  ! Read mesh generation namelist from file
  call get_initial_filename( filename )
  namelist_unit = open_file( filename )
  call read_cubedsphere_mesh_generator_namelist( namelist_unit, vm, 0 )
  call postprocess_cubedsphere_mesh_generator_namelist()
  call close_file( namelist_unit )
  deallocate( filename )

  ! Check the number of meshes requested?
  if (nmeshes < 1) then
    write(log_scratch_space,'(A,I0,A)') &
        'Invalid number of meshes requested, (',nmeshes,')'
    call log_event(log_scratch_space,LOG_LEVEL_ERROR)
  end if

  ! Check for missing data.
  if (ANY(edge_cells == imdi)) then
    write(log_scratch_space,'(A)') &
        'Missing data in namelist variable, edge_cells'
    call log_event(log_scratch_space,LOG_LEVEL_ERROR)
  end if

  allocate(csgen(nmeshes))
  allocate(ugrid_2d(nmeshes))

  ! Create each mesh in chain by setting details in ugrid_2d object
  ! using mesh_generators
  call log_event( "Generating cubed-sphere mesh(es) with...", LOG_LEVEL_INFO )

  do i=1, nmeshes

    ! Create object which can generate the biperiodic mesh from
    ! specified inputs.
    csgen(i) = gencube_ps_type( mesh_names(i), edge_cells(i), smooth_passes )

    write(log_scratch_space, "(2(A,I0))") &
        "  ndivs: ", edge_cells(i),       &
        ', smoothing passes:', smooth_passes
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )

    ! Mesh for the UGRID conforming NetCDF file is
    ! set by the generator passed to it
    call ugrid_2d(i)%set_by_generator(csgen(i))

  end do
  call log_event( "...generation complete.", LOG_LEVEL_INFO )


  ! Now the write out mesh to the NetCDF file
  do i=1, nmeshes

    if (.not. allocated(ugrid_file)) allocate(ncdf_quad_type::ugrid_file)
    call ugrid_2d(i)%set_file_handler(ugrid_file)

    if ( i==1 ) then
      call ugrid_2d(i)%write_to_file( trim(mesh_filename) )
    else
      call ugrid_2d(i)%append_to_file( trim(mesh_filename) )
    end if

    inquire(file=trim(mesh_filename), size=fsize)
    write( log_scratch_space, '(2(A,I0),A)')                    &
        'Adding ugrid mesh (ndivs:', edge_cells(i),') to ' //   &
        trim(adjustl(mesh_filename)) // ' - ', fsize,           &
        ' bytes written.'

    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    if (allocated(ugrid_file)) deallocate(ugrid_file)

  end do

  call ESMF_Finalize(rc=rc)

  if ( allocated( csgen ) )    deallocate (csgen)
  if ( allocated( ugrid_2d ) ) deallocate (ugrid_2d)

end program cubedsphere_mesh_generator
