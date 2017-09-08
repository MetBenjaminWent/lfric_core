!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @mainpage Biperiodic mesh generator
!>
!> @brief   Utility to generate a biperiodic surface mesh and write to a file
!>          which conforms to the UGRID format convention.
!> @details Usage:
!>
!>          biperiodic_mesh_generator <filename>
!>          filename - Controlling namelist file
!>
!-----------------------------------------------------------------------------
program biperiodic_mesh_generator

  use biperiodic_mesh_generator_config_mod,                                    &
                         only: read_biperiodic_mesh_generator_namelist,        &
                               postprocess_biperiodic_mesh_generator_namelist, &
                               edge_cells_x, edge_cells_y, domain_x, domain_y, &
                               nmeshes, mesh_names, mesh_filename

  use cli_mod,           only: get_initial_filename
  use constants_mod,     only: i_def, str_def
  use ESMF
  use genbiperiodic_mod, only: genbiperiodic_type
  use io_utility_mod,    only: open_file, close_file
  use log_mod,           only: log_scratch_space, log_event, &
                                LOG_LEVEL_INFO, LOG_LEVEL_ERROR
  use ncdf_quad_mod,     only: ncdf_quad_type
  use ugrid_2d_mod,      only: ugrid_2d_type
  use ugrid_file_mod,    only: ugrid_file_type

  implicit none

  type(ESMF_VM)  :: vm
  integer(i_def) :: rc

  character(:), allocatable :: filename
  integer(i_def)            :: namelist_unit

  type(genbiperiodic_type), allocatable :: bpgen(:)
  type(ugrid_2d_type),      allocatable :: ugrid_2d(:)
  class(ugrid_file_type),   allocatable :: ugrid_file

  integer(i_def) :: fsize
  integer(i_def) :: i

  character(str_def) :: rchar
  character(str_def) :: fmt_str

  ! Start up ESMF
  call ESMF_Initialize( vm=vm, rc=rc,                    &
                        logkindflag=ESMF_LOGKIND_SINGLE, &
                        defaultlogfilename="biperiodic.log" )
  if (rc /= ESMF_SUCCESS) call log_event( 'Failed to initialise ESMF.', &
                                          LOG_LEVEL_ERROR )

  ! Read mesh generation namelist from file
  call get_initial_filename( filename )
  namelist_unit = open_file( filename )
  call read_biperiodic_mesh_generator_namelist( namelist_unit, vm, 0 )
  call postprocess_biperiodic_mesh_generator_namelist( )
  call close_file( namelist_unit )
  deallocate( filename )

  ! Create objects to manipulate UGRID conforming NetCDF file
  allocate(bpgen(nmeshes))
  allocate(ugrid_2d(nmeshes))

  ! Create each mesh in chain by setting details in ugrid_2d object
  ! using mesh_generators
  do i=1, nmeshes

    ! Create object which can generate the biperiodic mesh from
    ! specified inputs.
    bpgen(i) = genbiperiodic_type( mesh_names(i),   &
                                   edge_cells_x(i), &
                                   edge_cells_y(i), &
                                   domain_x,        &
                                   domain_y )

    fmt_str ='(A,I0)'
    write(log_scratch_space, fmt_str) &
        "Generating biperiodic mesh: "//trim(mesh_names(i))
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write(log_scratch_space, fmt_str) &
        "  Cells in x:  ", edge_cells_x(i)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write(log_scratch_space, fmt_str) &
        "  Cells in y:  ", edge_cells_y(i)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    fmt_str ='(A)'
    write(rchar, '(F6.1)') domain_x / edge_cells_x(i)
    write(log_scratch_space, fmt_str) &
        "  Cell width:  "//trim(adjustl(rchar))
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    write(rchar, '(F6.1)') domain_x / edge_cells_y(i)
    write(log_scratch_space, fmt_str) &
        "  Cell height: "//trim(adjustl(rchar))
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    ! Mesh for the UGRID conforming NetCDF file is
    ! set by the generator passed to it
    call ugrid_2d(i)%set_by_generator(bpgen(i))

  end do

  call log_event( "...generation complete.", LOG_LEVEL_INFO )


  ! Now the write out mesh to the NetCDF file
  do i=1, nmeshes

    if (.not. allocated(ugrid_file)) allocate(ncdf_quad_type::ugrid_file)

    call ugrid_2d(i)%set_file_handler(ugrid_file)

    if (i==1) then
      call ugrid_2d(i)%write_to_file( trim(mesh_filename) )
    else
      call ugrid_2d(i)%append_to_file( trim(mesh_filename) )
    end if

    inquire(file=mesh_filename, size=fsize)
    write( log_scratch_space, '(A,I0,A)')                          &
        'Writing ugrid mesh to ' // trim(adjustl(mesh_filename))// &
        ' - ', fsize, ' bytes written.'

    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    if (allocated(ugrid_file)) deallocate(ugrid_file)

  end do

  call ESMF_Finalize(rc=rc)

  if ( allocated( bpgen ) )    deallocate (bpgen)
  if ( allocated( ugrid_2d ) ) deallocate (ugrid_2d)

end program biperiodic_mesh_generator
