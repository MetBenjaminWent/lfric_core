!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!------------------------------------------------------------------------------
!
!>
!> @brief   Holds and manages the multiple global meshes used to setup a model
!>          run.
!>
!> @details A container which holds a collection of global meshes
!>          It will handle the creation and storing of requested global meshes.
!
module global_mesh_collection_mod

  use constants_mod,   only: r_def, i_def, imdi, str_max_filename
  use mesh_mod,        only: mesh_type
  use log_mod,         only: log_event, log_scratch_space,                     &
                             LOG_LEVEL_ERROR, LOG_LEVEL_INFO
  use linked_list_mod, only: linked_list_type, linked_list_item_type
  use global_mesh_mod, only: global_mesh_type
  use partition_mod,   only: partition_type

  implicit none

  private

  type, public :: global_mesh_collection_type
    private
    type(linked_list_type) :: global_mesh_list

  contains
    procedure, public :: add_new_global_meshes ! Will use multimesh ugrid files
    procedure, public :: add_new_global_mesh   ! Deprecated on multimesh ugrid file
    procedure, public :: add_unit_test_global_mesh
    procedure, public :: get_global_mesh

    procedure, public :: clear

    final             :: global_mesh_collection_destructor

  end type global_mesh_collection_type

  interface global_mesh_collection_type
    module procedure global_mesh_collection_constructor
  end interface

  ! Module variable allows access to the single mesh collection
  type(global_mesh_collection_type), public :: global_mesh_collection

contains

!> Constructs the mesh collection object
!> @return self The constructed mesh collection object
function global_mesh_collection_constructor() result(self)

  implicit none
  type(global_mesh_collection_type) :: self

  self%global_mesh_list = linked_list_type()

end function global_mesh_collection_constructor

!========================== Multi-Global_Mesh UGRID file=======================
!
!> @brief Add global mesh objects to the the global_mesh_collection from a
!>        ugrid file which contain one or more global meshes and any
!>        associated intermesh connectivity data. By convention, the 
!>        first global mesh listed in the file is assumed to be the
!>        intended primal mesh for this ugrid file. This primal mesh must
!>        match the primal mesh in the model if a primal mesh already exists
!> @param [in] filename    Filename of ugrid file containing details of global
!>                         mesh objects and associated intermesh connectiviy.
function add_new_global_meshes(self, filename) result (global_mesh_id)

  implicit none

  class(global_mesh_collection_type), intent(inout) :: self
  character(len=str_max_filename),    intent(in)    :: filename

  integer(i_def) :: global_mesh_id

  call log_event( "Add_Mesh: Requires multi-mesh ugrid file functionality" &
                , LOG_LEVEL_INFO )

  ! 1.0 Call routine to read ugrid file and populate global mesh objects
  !     with the data

  ! 2.0 Loop over each global mesh object
  ! 2.1    Calculate a hash and check to see if object already exists
  !        in collection, use existing mesh if present.
  ! 2.2    Check to see if ugrid file contained mapping data to other
  !        meshes in the ugrid file
  ! 2.3    If so obtain hash for target mesh, does the global map in
  !        the collection contain a map to a existing mesh with this hash?
  !        If so do not use this particular mapping data
  !        (Should have warnings if the case of conflicting hashes?)

  return
end function add_new_global_meshes
!
!========================== Multi-Global_Mesh UGRID file=======================

  
!> @brief Add a global mesh objects to the collection from ugrid files
!>        which contain a single global meshes per file. Maps between
!>        each global mesh are created on the fly dependent on the input
!>        order. 
!> @param [in] filenames  File containing details of a single global mesh
!>                        object.
!> @return global_mesh_id Integer id of the global mesh added to collection
function add_new_global_mesh( self, filename) result (global_mesh_id)

  implicit none

  class(global_mesh_collection_type), intent(inout) :: self
  character(len=str_max_filename), intent(in)       :: filename
  
  integer(i_def) :: global_mesh_id

  type (global_mesh_type) :: global_mesh


  ! Need to generate a hash on each global mesh object read in. Though 
  ! Matthew thinks this best wriiten in C

  ! For there is no checking on uniqueness of global meshes read in
  ! Assumption that there is one mesh per file.

  global_mesh = global_mesh_type( trim(filename) )
  global_mesh_id = global_mesh%get_id()
  call self%global_mesh_list%insert_item( global_mesh )

  return
end function add_new_global_mesh



!> @brief Call to a global mesh collection to add a global mesh
!>        object for unit tests.
!> @return global_mesh_id Integer id of the global mesh added to collection
function add_unit_test_global_mesh(self) result(global_mesh_id)

  implicit none

  class(global_mesh_collection_type), intent(inout) :: self
  integer(i_def) :: global_mesh_id

  type (global_mesh_type) :: global_mesh

  global_mesh = global_mesh_type( )
  global_mesh_id = global_mesh%get_id()
  call self%global_mesh_list%insert_item( global_mesh )

  return
end function add_unit_test_global_mesh



!> @brief Call to a global mesh collection to return a global mesh with
!>        a specified global_mesh_id
!> @param [in] global_mesh_id  Integer id of global mesh object required
function get_global_mesh( self, global_mesh_id ) result( global_mesh )

  implicit none

  class(global_mesh_collection_type) :: self
  integer(i_def), intent(in) :: global_mesh_id

  type(global_mesh_type), pointer :: global_mesh

  ! Pointer to linked list - used for looping through the list
  type(linked_list_item_type),pointer :: loop => null()

  ! start at the head of the mesh collection linked list
  loop => self%global_mesh_list%get_head()

  do
    ! If list is empty or we're at the end of list and we didn't find the
    ! mesh_id, return a null pointer
    if ( .not. associated(loop) ) then
      nullify(global_mesh)
      exit
    end if

    ! Otherwise search list for the id we want
    if ( global_mesh_id == loop%payload%get_id() ) then
      ! 'cast' to the global_mesh_type 
      select type(m => loop%payload)
        type is (global_mesh_type)
          global_mesh => m
      end select
      exit
    end if
    loop => loop%next
  end do

end function get_global_mesh


!> Clear all items from the mesh collection linked list
subroutine clear(self)

  implicit none

  class(global_mesh_collection_type), intent(inout) :: self

  call self%global_mesh_list%clear()

end subroutine clear


! Mesh collection destructor
subroutine global_mesh_collection_destructor(self)

  implicit none

  type (global_mesh_collection_type), intent(inout) :: self

  call self%clear()

end subroutine global_mesh_collection_destructor


end module global_mesh_collection_mod
