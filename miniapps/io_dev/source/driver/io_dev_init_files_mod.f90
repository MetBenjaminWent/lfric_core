!-------------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Sets up file configuration fof IO_Dev
!> @details Collects file configuration information and formats it so that it
!>          can be passed to the infrastructure
module io_dev_init_files_mod

  use constants_mod,       only: i_def, i_native, &
                                 str_def, str_max_filename
  use file_mod,            only: file_type
  use lfric_xios_file_mod, only: lfric_xios_file_type
  use driver_io_mod,       only: append_file_to_list
  ! Configuration modules
  use files_config_mod,    only: diag_stem_name,            &
                                 checkpoint_stem_name,      &
                                 start_dump_filename,       &
                                 start_dump_directory,      &
                                 time_varying_input_path,   &
                                 time_data_path
  use io_dev_config_mod,   only: field_initialisation,            &
                                 field_initialisation_start_dump, &
                                 time_variation,                  &
                                 time_variation_ancil
  use io_config_mod,       only: use_xios_io,               &
                                 diagnostic_frequency,      &
                                 checkpoint_write,          &
                                 checkpoint_read,           &
                                 write_diag, write_dump
  use time_config_mod,     only: timestep_start,            &
                                 timestep_end

  implicit none

  private
  public :: init_io_dev_files

  contains

  subroutine init_io_dev_files(files_list)

    implicit none

    class(file_type), allocatable, intent(out) :: files_list(:)

    type(lfric_xios_file_type)      :: tmp_xios_file
    character(len=str_max_filename) :: checkpoint_write_fname, &
                                       checkpoint_read_fname,  &
                                       dump_fname,             &
                                       input_fname
    integer(i_def)                  :: ts_start, ts_end
    integer(i_native)               :: rc

    ! Get time configuration in integer form
    read(timestep_start,*,iostat=rc)  ts_start
    read(timestep_end,*,iostat=rc)  ts_end

    if ( use_xios_io) then

      ! Setup diagnostic output file
      if ( write_diag ) then
        call tmp_xios_file%file_new("io_dev_diag")
        call tmp_xios_file%configure(xios_id="io_dev_diag", freq=diagnostic_frequency)
        call append_file_to_list(tmp_xios_file, files_list)
      end if

      ! Setup dump-writing context information
      if ( write_dump ) then
        ! Create dump filename from base name and end timestep
        write(dump_fname,'(A,A,I6.6)') &
           trim(start_dump_directory)//'/'//trim(start_dump_filename),"_", &
           timestep_end

        ! Setup dump file for end timestep
        call tmp_xios_file%file_new(dump_fname)
        call tmp_xios_file%configure(xios_id="io_dev_dump_out", freq=ts_end)
        call append_file_to_list(tmp_xios_file, files_list)
      end if

      ! Setup dump-reading context information
      if ( field_initialisation == field_initialisation_start_dump ) then
        ! Create dump filename from stem
        write(dump_fname,'(A)') trim(start_dump_directory)//'/'// &
                                trim(start_dump_filename)

        ! Setup dump file
        call tmp_xios_file%file_new(dump_fname)
        call tmp_xios_file%configure(xios_id="io_dev_dump_in")
        call append_file_to_list(tmp_xios_file, files_list)
      end if

      ! Setup checkpoint writing context information
      if ( checkpoint_write ) then
        ! Create checkpoint filename from stem and end timestep
        write(checkpoint_write_fname,'(A,A,I6.6)') &
                             trim(checkpoint_stem_name),"_", ts_end

        call tmp_xios_file%file_new(checkpoint_write_fname)
        call tmp_xios_file%configure( xios_id="io_dev_checkpoint_write",  &
                                 freq=ts_end - ts_start, &
                                 field_group_id="checkpoint_fields" )
        call append_file_to_list(tmp_xios_file, files_list)
      end if

      ! Setup checkpoint reading context information
      if ( checkpoint_read ) then
        ! Create checkpoint filename from stem and (start - 1) timestep
        write(checkpoint_read_fname,'(A,A,I6.6)') &
                     trim(checkpoint_stem_name),"_", (ts_start - 1)

        call tmp_xios_file%file_new(checkpoint_read_fname)
        call tmp_xios_file%configure(xios_id="io_dev_checkpoint_read")
        call append_file_to_list(tmp_xios_file, files_list)
      end if

      ! Setup time-varying input files
      if ( time_variation == time_variation_ancil ) then
        ! Set time-varying input filename from namelist
        write(input_fname,'(A)') trim(start_dump_directory)//'/'// &
                                 trim(time_varying_input_path)
        call tmp_xios_file%file_new(input_fname)
        call tmp_xios_file%configure(xios_id="io_dev_time_varying_input")
        call append_file_to_list(tmp_xios_file, files_list)

        ! Set time data input filename from namelist
        write(input_fname,'(A)') trim(start_dump_directory)//'/'// &
                                 trim(time_data_path)
        call tmp_xios_file%file_new(input_fname)
        call tmp_xios_file%configure(xios_id="io_dev_times")
        call append_file_to_list(tmp_xios_file, files_list)

      end if

    end if ! use_xios_io

  end subroutine init_io_dev_files

end module io_dev_init_files_mod
