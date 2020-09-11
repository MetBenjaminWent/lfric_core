!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief   Holds and manages the multiple global meshes used to setup a model
!>          run.
!>
!> @details A container which holds a collection of global meshes.
!>          It will handle the creation and storing of requested global meshes.
!
module global_mesh_collection_mod

  use constants_mod,         only : r_def, i_def, IMDI, i_native, &
                                    str_max_filename, str_def
  use global_mesh_mod,       only : global_mesh_type
  use linked_list_mod,       only : linked_list_type, linked_list_item_type
  use log_mod,               only : log_event, log_scratch_space, &
                                    LOG_LEVEL_ERROR, LOG_LEVEL_TRACE

  use development_config_mod, only : implement_consolidated_multigrid

  implicit none

  private

  type, public :: global_mesh_collection_type

    private

    ! List of the global_mesh_type objects in this collection.
    type(linked_list_type) :: global_mesh_list

    ! Number of panels in the mesh layout. npanels is set to be the same as
    ! the 1st global mesh loaded into the collection. All subsequent global
    ! meshes added should have been specified with the same nume of panels.
    !> @deprecated  Once multiple global meshes and associated mappings
    !>              are available in ugrid files.

    ! At present, global_mesh_type objects which are described from ugrid
    ! files which contain details for only one global mesh and thus contain no
    ! information about mapping between global meshes at different resolutions.
    ! As a consequence of this, mappings are calculated by the
    ! global_mesh_collection as between subsequent global meshes as they are
    ! added to the global mesh_collection. This calculation requires that all
    ! meshes have the same number of panels in the mesh. npanels is set to
    ! be the same as the 1st global mesh loaded into the collection.
    integer(i_def)         :: npanels = IMDI

    ! Pointer to global_mesh_type object in linked list. This global mesh
    ! object will be use as the source mesh when for global mesh map creation
    ! when the next global mesh is added to the collection.
    ! THIS IS TEMPORARY AND SHOULD BE REMOVED WHEN GLOBAL MESH MAPS ARE
    ! READ DIRECTLY FROM THE UGRID MESH FILE
    type(global_mesh_type),  pointer :: source_global_mesh => null()

    character(str_def), allocatable :: name_tags(:)
    integer(i_def),     allocatable :: name_ids(:)

  contains
    procedure, public  :: add_new_global_mesh
    procedure, public  :: add_unit_test_global_mesh
    procedure, public  :: set_next_source_mesh
    procedure, private :: map_global_meshes

    procedure, public  :: n_meshes
    procedure, public  :: get_mesh_names
    procedure, public  :: get_mesh_id

    procedure, public  :: get_mesh_by_name
    procedure, public  :: get_mesh_by_id
    generic,   public  :: get_global_mesh => get_mesh_by_id, &
                                             get_mesh_by_name
    procedure, public  :: check_for

    procedure, public  :: clear
    final              :: global_mesh_collection_destructor

  end type global_mesh_collection_type

  interface global_mesh_collection_type
    module procedure global_mesh_collection_constructor
  end interface

  !> @brief Singleton instance of a global_mesh_collection_type object.
  !>
  type(global_mesh_collection_type), public, allocatable :: &
      global_mesh_collection

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Constructs an empty mesh collection object.
  !>
  !> @return The constructed mesh collection object.
  !>
  function global_mesh_collection_constructor() result(self)

    implicit none

    type(global_mesh_collection_type) :: self

    self%global_mesh_list = linked_list_type()

  end function global_mesh_collection_constructor


  !===========================================================================
  !> @brief Destructor tears down object prior to being freed.
  !>
  subroutine global_mesh_collection_destructor(self)

    ! Object finalizer
    implicit none

    type (global_mesh_collection_type), intent(inout) :: self

    call self%clear()

    return
  end subroutine global_mesh_collection_destructor


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Adds a global mesh object to the collection from a ugrid file
  !>        which contains the requested global mesh topology.
  !>
  !> Maps between each global mesh are created on the fly between subsequent
  !> global meshes in the collection in the order that they are added.
  !>
  !> At present, there is no way to uniquely identify each global mesh that
  !> is read in. Global meshes are currently identified by an integer ID
  !> which is assigned on creation of the global mesh object.
  !>
  !> When multiple global meshes per file are possible there will be a need
  !> to use a hash function to identify each global mesh object read in.
  !>
  !> @param[in] filename  File containing details of global
  !>                      mesh object(s).
  !> @param[in] global_mesh_name The name of the global mesh
  !>                             to read from the file.
  !> @param[in] npanels Number of panels in the passed mesh.
  !>
  !> @return ID of the global mesh added to collection.
  !>
  function add_new_global_mesh( self,             &
                                filename,         &
                                global_mesh_name, &
                                npanels ) result( global_mesh_id )

    implicit none

    class(global_mesh_collection_type), intent(inout) :: self
    character(str_max_filename),        intent(in)    :: filename
    character(str_def),                 intent(in)    :: global_mesh_name
    integer(i_def),                     intent(in)    :: npanels


    integer(i_def)          :: n_global_meshes
    integer(i_def)          :: global_mesh_id
    type (global_mesh_type) :: global_mesh_to_add
    type (global_mesh_type), pointer :: global_mesh_at_tail => null()

    integer(i_def) :: source_global_mesh_id
    integer(i_def) :: target_global_mesh_id
    integer(i_native) :: i

    ! Pointer to linked list - used for looping through the list
    type(linked_list_item_type), pointer :: list_item => null()

    ! 1.0 Perform some checks based on the existings contents
    !     of the collection
    !=================================================================
    n_global_meshes = self%global_mesh_list%get_length()

    ! 1.1 If this is the first global mesh in the collection, it
    !     will be used to set the constraint for the number of
    !     panels used by subsequent meshes added to this collection.
    if (n_global_meshes < 1) then
      self%npanels = npanels
    else
      if (self%npanels /= npanels) then
        write(log_scratch_space,'(A,I0,A)')        &
            'This global mesh collection is '//    &
            'for global meshes of comprising of ', &
            self%npanels, ' panels.'

        call log_event(log_scratch_space,LOG_LEVEL_ERROR)
      end if
    end if

    ! 2.0 Read in the requested mesh topology from file.
    !=================================================================
    global_mesh_to_add = global_mesh_type( trim(filename), &
                                           global_mesh_name )
    global_mesh_id = global_mesh_to_add%get_id()

    if (implement_consolidated_multigrid) then

      ! Protect this until implementation as the current unit-tests
      ! (which will be removed on implementation in gungho/lfric_atm)
      ! use a common mesh tag name (unit-test). they also rely on mesh
      ! maps being created as global mesh maps are added. This will be
      ! removed at implementation.


      ! Check list of tag names to see if mesh is already in collection
      ! As these are assumed to be read in from the same mesh input
      ! file specifed for a given run, if the name is already in the
      ! collection it will be the same mesh.
      if (self%check_for(global_mesh_name)) then
        do i=1, size(self%name_tags)
          if ( trim(self%name_tags(i)) == trim(global_mesh_name) ) then
            return
          end if
        end do
      end if

      ! Read in the named mesh to add from the mesh input file
      global_mesh_to_add = global_mesh_type( filename, &
                                             global_mesh_name )

      global_mesh_id = global_mesh_to_add%get_id()

      ! Now add the requested mesh to the collection
      call self%global_mesh_list%insert_item( global_mesh_to_add )

    else

      ! ==========================================================
      ! ORIGINAL CODE
      ! ==========================================================
      ! This will be removed if the a followup ticket to implement
      ! consolidated multigrid into lfric apps.

      ! If there is at least one other mesh in the collection, attempt
      ! to create a map from it to the global mesh being added
      if (n_global_meshes >= 1) then

        ! In order to map the global meshes the two global
        ! meshes should have the same number of panels
        source_global_mesh_id = self%source_global_mesh%get_id()
        target_global_mesh_id = global_mesh_id

        call self%global_mesh_list%insert_item( global_mesh_to_add )
        call self%map_global_meshes( source_global_mesh_id, &
                                     target_global_mesh_id )

        call self%set_next_source_mesh( target_global_mesh_id )

      else

        call self%global_mesh_list%insert_item( global_mesh_to_add )
        call self%set_next_source_mesh( global_mesh_id )

      end if

    end if ! Implement_consolidated_multigrid

    nullify( global_mesh_at_tail )
    nullify( list_item )

    if ( implement_consolidated_multigrid ) &
        call update_tags( self,             &
                          global_mesh_name, &
                          global_mesh_id )
    return
  end function add_new_global_mesh


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Returns the number of meshes in the collection.
  !>
  !> @detail This function returns the number of unique mesh
  !>         tag names in this collection.
  !>
  !> @return Number of mesh tag names avaible to query.
  !>
  function n_meshes(self) result(number_of_meshes)

    implicit none

    class(global_mesh_collection_type), intent(in) :: self

    integer(i_def) :: number_of_meshes

    if (allocated(self%name_tags)) then
      number_of_meshes = size(self%name_tags)
    else
      number_of_meshes = 0
    end if

    return
  end function n_meshes


  !===========================================================================
  ! Internal routine to generate cell mappings (if possible) between
  ! subsequent global meshes added to the collection. This routine is present
  ! as mapping data is not currently available in the global mesh ugrid
  ! files. Mappings are created on a "per mesh panel" basis and are checked
  ! against the variable npanels.
  !
  !> @deprecated When ugrid files with multiple meshes/mappings are
  !>             available
  !
  subroutine map_global_meshes( self,                   &
                                source_global_mesh_id,  &
                                target_global_mesh_id )

    implicit none

    class(global_mesh_collection_type), intent(inout) :: self
    integer(i_def), intent(in) :: source_global_mesh_id
    integer(i_def), intent(in) :: target_global_mesh_id

    type (global_mesh_type), pointer :: source_mesh => null()
    type (global_mesh_type), pointer :: target_mesh => null()
    type (global_mesh_type), pointer :: coarse_mesh => null()
    type (global_mesh_type), pointer :: fine_mesh   => null()

    integer(i_def) :: i,j,n,count       ! counters
    integer(i_def) :: coarse_panel_start_id
    integer(i_def) :: fine_panel_start_id
    integer(i_def) :: coarse_panel_end_id
    integer(i_def) :: fine_panel_end_id
    integer(i_def) :: coarse_id
    integer(i_def) :: fine_id

    integer(i_def) :: coarse_edge_cells
    integer(i_def) :: fine_edge_cells
    integer(i_def) :: coarse_ncells
    integer(i_def) :: fine_ncells
    integer(i_def) :: source_ncells
    integer(i_def) :: target_ncells

    integer(i_def) :: factor ! ratio of edge_cells

    integer(i_def) :: start_x
    integer(i_def) :: start_y
    integer(i_def) :: end_x
    integer(i_def) :: end_y

    integer(i_def), allocatable :: coarse_panel_ids(:,:)
    integer(i_def), allocatable :: fine_panel_ids(:,:)

    integer(i_def), allocatable :: coarse_to_fine_gid_map(:,:)
    integer(i_def), allocatable :: fine_to_coarse_gid_map(:,:)
    integer(i_def), allocatable :: tmp_panel_ids(:)


    ! Check for identical source and target global meshes
    if (source_global_mesh_id == target_global_mesh_id) then
      write(log_scratch_space,'(A)') &
          "Cannot have identical consective global mesh ids"
      call log_event(log_scratch_space,LOG_LEVEL_TRACE)
      return
    end if

    source_mesh => self%get_global_mesh(source_global_mesh_id)
    target_mesh => self%get_global_mesh(target_global_mesh_id)

    ! Figure out which of the meshes has fewer
    ! cells i.e. the coarse one
    source_ncells = source_mesh%get_ncells()
    target_ncells = target_mesh%get_ncells()

    if ( mod(source_ncells, target_ncells) == 0 .or. &
        mod(target_ncells, source_ncells) == 0 ) then

      if (source_ncells == target_ncells) then
        write(log_scratch_space, '(A)')                           &
            'Meshes have equal number of cells, direct mapping, ' &
            // 'no mapping created'
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
        return

      else if (target_ncells > source_ncells) then
        coarse_mesh   => source_mesh
        fine_mesh     => target_mesh
        coarse_ncells = source_ncells
        fine_ncells   = target_ncells

      else
        coarse_mesh   => target_mesh
        fine_mesh     => source_mesh
        coarse_ncells = target_ncells
        fine_ncells   = source_ncells

      end if

    else
      ! The number of cells in one of these meshes
      ! needs to be a factor of the number cells in the other.
      write(log_scratch_space, '(2(A,I0))')                                   &
          'Unable to generate intermesh connectivity for global meshes: id:', &
          source_global_mesh_id, ';id:', target_global_mesh_id
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
      return

    end if

    coarse_id = coarse_mesh%get_id()
    fine_id   = fine_mesh%get_id()

    ! Now we know which of the source and target are the
    ! coarser and finer meshes, generate a coarse_to_fine
    ! and fine_to_coarse map per panel
    coarse_edge_cells = int(sqrt(real(coarse_ncells/self%npanels, r_def)), &
                            i_def)
    fine_edge_cells   = int(sqrt(real(fine_ncells/self%npanels,   r_def)), &
                            i_def)
    factor            = fine_edge_cells/coarse_edge_cells

    allocate( coarse_panel_ids ( coarse_edge_cells, coarse_edge_cells ))
    allocate( fine_panel_ids   ( fine_edge_cells,   fine_edge_cells   ))
    allocate( coarse_to_fine_gid_map( factor**2, coarse_ncells ) )
    allocate( fine_to_coarse_gid_map( 1, fine_ncells ))

    do n=1, self%npanels

      ! Get the id of the initial cell for the panel
      coarse_panel_start_id = ((n-1)*coarse_ncells / self%npanels) + 1
      fine_panel_start_id   = ((n-1)*fine_ncells   / self%npanels) + 1

      ! =======================================================
      ! Populate arrays with ids and reshape to panel layout
      ! =======================================================

      ! For coarse mesh
      allocate(tmp_panel_ids(coarse_ncells))
      count=1
      coarse_panel_end_id = coarse_panel_start_id + coarse_ncells-1
      do i=coarse_panel_start_id, coarse_panel_end_id
        tmp_panel_ids(count) = i
        count=count+1
      end do

      coarse_panel_ids = reshape( tmp_panel_ids,(/coarse_edge_cells, &
                                                  coarse_edge_cells/) )
      deallocate(tmp_panel_ids)

      ! For fine mesh
      allocate(tmp_panel_ids(fine_ncells))
      count=1
      fine_panel_end_id =  fine_panel_start_id + fine_ncells-1
      do i=fine_panel_start_id, fine_panel_end_id
        tmp_panel_ids(count) = i
        count=count+1
      end do
      fine_panel_ids = reshape(tmp_panel_ids,(/fine_edge_cells, &
                                               fine_edge_cells/))
      deallocate(tmp_panel_ids)

      ! ===============================
      ! Populate intermesh gid maps
      ! ===============================

      ! Coarse to Fine
      do j=1, coarse_edge_cells
        start_y = ((j-1)*factor) + 1
        end_y   = start_y + factor - 1
        do i=1, coarse_edge_cells
          start_x = ((i-1)*factor) + 1
          end_x   = start_x + factor - 1

          coarse_to_fine_gid_map(:,coarse_panel_ids(i,j)) =            &
              reshape( fine_panel_ids( start_x:end_x, start_y:end_y ), &
                      (/factor**2/) )
        end do
      end do

      ! Fine to Coarse
      coarse_panel_end_id = coarse_panel_start_id &
                            + (coarse_ncells/self%npanels) - 1

      do j=coarse_panel_start_id, coarse_panel_end_id
        do i=1, factor**2
          fine_to_coarse_gid_map(1, coarse_to_fine_gid_map(i,j)) = j
        end do
      end do

    end do

    ! Now we have the intermesh gid maps from coarse to fine
    ! and fine to coarse meshes. Create a global_mesh_map for each
    ! direction.
    call coarse_mesh % add_global_mesh_map( fine_mesh, &
                                            coarse_to_fine_gid_map )
    call fine_mesh   % add_global_mesh_map( coarse_mesh, &
                                            fine_to_coarse_gid_map )

    deallocate( coarse_panel_ids )
    deallocate( fine_panel_ids   )
    deallocate( coarse_to_fine_gid_map )
    deallocate( fine_to_coarse_gid_map )

    nullify(source_mesh, target_mesh, coarse_mesh, fine_mesh)

    return
  end subroutine map_global_meshes

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Specifies the global mesh to be used as the source mesh for
  !>        mapping when the next global mesh is added to the collection.
  !>
  !> @deprecated When ugrid files with multiple meshes/mappings are available.
  !>
  subroutine set_next_source_mesh(self, global_mesh_id)

    implicit none

    class(global_mesh_collection_type), intent(inout) :: self

    integer(i_def), intent(in) :: global_mesh_id

    ! Pointer to linked list - used for looping through the list
    class(linked_list_item_type), pointer :: loop => null()

    ! start at the head of the mesh collection linked list
    loop => self%global_mesh_list%get_head()

    if (self%global_mesh_list%item_exists(global_mesh_id)) then
      do
        ! Search list for the id we want
        if ( global_mesh_id == loop%payload%get_id() ) then
          ! 'cast' to the global_mesh_type
          select type(m => loop%payload)
          type is (global_mesh_type)
            self%source_global_mesh => m
          end select
          exit
        end if
        loop => loop%next
      end do
    else
      write(log_scratch_space,'(A)') &
          "Invalid global mesh id: does not exist in collection"
      call log_event(log_scratch_space,LOG_LEVEL_ERROR)
    end if

    nullify(loop)

    return
  end subroutine set_next_source_mesh


  !===========================================================================
  !> @brief Adds hardwired global mesh object to the collection (for unit
  !>        testing).
  !>
  !> @return ID of the global mesh added to collection.
  !>
  function add_unit_test_global_mesh(self) result(global_mesh_id)

    implicit none

    class(global_mesh_collection_type), intent(inout) :: self
    integer(i_def) :: global_mesh_id

    type (global_mesh_type) :: global_mesh

    integer(i_def)     :: n_global_meshes
    character(str_def) :: global_mesh_name

    n_global_meshes = self%global_mesh_list%get_length()

    if (n_global_meshes < 1) then
      self%npanels   = 1
    else
      if (self%npanels /= 1) then
        write(log_scratch_space,'(A,I0,A)')        &
            'This global mesh collection is '//    &
            'for global meshes of comprising of ', &
            self%npanels, ' panels.'

        call log_event(log_scratch_space,LOG_LEVEL_ERROR)
      end if
    end if

    global_mesh      = global_mesh_type()
    global_mesh_name = global_mesh%get_mesh_name()
    global_mesh_id   = global_mesh%get_id()
    call self%global_mesh_list%insert_item( global_mesh )

    if ( implement_consolidated_multigrid ) &
        call update_tags( self,             &
                          global_mesh_name, &
                          global_mesh_id )

    return
  end function add_unit_test_global_mesh


  !===========================================================================
  !> @brief   Requests a global mesh object with the specified mesh name
  !>          from the collection.
  !> @details The name of the global mesh is determined by the name given
  !>          to the mesh topology in the mesh input file.
  !>
  !> @param[in] global_mesh_name Name tag of global mesh object to retrieve.
  !>
  !> @return    global_mesh      Pointer to global mesh object with
  !>                             requested name if present in collection.
  !>                             A null pointer is returned if there is no
  !>                             mesh with the requested name.
  !>
  function get_mesh_by_name( self, global_mesh_name ) result( global_mesh )

    implicit none

    class(global_mesh_collection_type) :: self
    character(str_def), intent(in)     :: global_mesh_name

    type(global_mesh_type), pointer    :: global_mesh

    integer(i_def) :: n_global_meshes, i
    integer(i_def) :: global_mesh_id

    n_global_meshes = size(self%name_tags)

    global_mesh => null()

    do i=1 , n_global_meshes
      if (trim(global_mesh_name) == trim(self%name_tags(i))) then
        global_mesh_id = self%name_ids(i)
        global_mesh => global_mesh_collection % get_global_mesh( global_mesh_id )
        exit
      end if
    end do

    if ( .not. associated(global_mesh) ) then
      write(log_scratch_space,'(A)') &
          trim(global_mesh_name)//' not found in the global mesh collection.'
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if

    return
  end function get_mesh_by_name


  !===========================================================================
  !> @brief   Requests a global mesh object with specified mesh ID from
  !>          the collection.
  !> @details The mesh ID is an internal integer ID assigned when the
  !>          mesh is instantiated. It is dependent on the order in which
  !>          the global mesh objects were created (not to be confused with
  !>          the order in which they were added to the collection).
  !>
  !> @param[in] global_mesh_id   ID of mesh to retrieve.
  !>
  !> @return    global_mesh      Pointer to global mesh object with
  !>                             the requested ID if present in collection.
  !>                             A null pointer is returned if there is no
  !>                             mesh with the requested ID.
  !>
  function get_mesh_by_id( self, global_mesh_id ) result( global_mesh )

    implicit none

    class(global_mesh_collection_type) :: self
    integer(i_def), intent(in)         :: global_mesh_id

    type(global_mesh_type), pointer    :: global_mesh


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

    nullify(loop)

  end function get_mesh_by_id


  !===========================================================================
  !> @brief Returns mesh tag names of mesh objects in the collection.
  !>
  !> @return mesh_names  String array <<allocatable>> of mesh names in
  !>                     collection.
  !>
  function get_mesh_names( self ) result( mesh_names )

    implicit none

    class(global_mesh_collection_type) :: self
    character(str_def), allocatable :: mesh_names(:)

    if (allocated(mesh_names)) deallocate(mesh_names)
    if (allocated(self%name_tags)) &
        allocate(mesh_names, source=self%name_tags)

    return
  end function  get_mesh_names


  !===========================================================================
  !> @brief Returns mesh ID specified mesh in the collection.
  !>
  !> @param[in] mesh_name  Mesh name of object to return ID.
  !>
  !> @return    mesh_id    Integer ID of global mesh object if present
  !>                       in collection.
  !>
  function get_mesh_id( self, mesh_name ) result( mesh_id )

    implicit none

    class(global_mesh_collection_type) :: self
    character(str_def), intent(in) :: mesh_name
    integer(i_def) :: n_meshes,i
    integer(i_def) :: mesh_id

    mesh_id = imdi
    n_meshes = self%global_mesh_list%get_length()

    do i=1 , n_meshes
      if (trim(mesh_name) == trim(self%name_tags(i))) then
        mesh_id = self%name_ids(i)
        exit
      end if
    end do

    return
  end function get_mesh_id


  !===========================================================================
  !> @brief Queries the collection as to the presence of a
  !>        mesh object with the specified mesh name in the collection
  !>
  !> @param[in] global_mesh_name  Mesh name of object to check for.
  !>
  !> @return    logical          .true. if global mesh object present
  !>                             in collection.
  !>
  function check_for(self, global_mesh_name) result(answer)

    implicit none

    class(global_mesh_collection_type), intent(in) :: self
    character(str_def),                 intent(in) :: global_mesh_name

    logical :: answer

    answer = .false.
    if (allocated(self%name_tags)) then
      answer = any(self%name_tags == global_mesh_name)
    end if

    return
  end function check_for


  !===========================================================================
  !> @brief Forced clear of all the global mesh objects in the collection.
  !>
  !> This routine should not need to be called manually except (possibly) in
  !> pfunit tests
  !>
  subroutine clear(self)

    ! Clear all items from the linked list in the collection
    implicit none

    class(global_mesh_collection_type), intent(inout) :: self

    call self%global_mesh_list%clear()

    nullify(self%source_global_mesh)

    if (allocated(self%name_tags)) deallocate(self%name_tags)
    if (allocated(self%name_ids))  deallocate(self%name_ids)

    return
  end subroutine clear

  !===========================================================================
  !> @brief PRIVATE SUBROUTINE: Updates the name_tags and name_ids arrays when
  !>        a mesh is added to the collection
  !>
  !> @param[inout] self    Global mesh collection object to update tags/ids
  !> @param[in] mesh_name  The mesh tag name to be added, must not already exist
  !>                       in name_tags array
  !> @param[in] mesh_id    ID to use to map with the provided mesh name.
  !>
  subroutine update_tags( self, mesh_name, mesh_id )

    implicit none

    class(global_mesh_collection_type), intent(inout) :: self
    character(str_def), intent(in) :: mesh_name
    integer(i_def), intent(in) :: mesh_id

    character(str_def), allocatable :: new_tag_list(:)
    integer(i_def),     allocatable :: new_id_list(:)
    integer(i_def)                  :: i, n_tags

    ! 1. Check that the id provided exists
    if (.not. self%global_mesh_list%item_exists(mesh_id)) then
      write(log_scratch_space,'(A,I0,A)')             &
          'Unable to update tags, global_mesh id: ',  &
           mesh_id,' not found in collection'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    ! 2. Check to see if the mesh name/mesh_id entry exists already
    if ( allocated(self%name_tags) ) then
      do i=1, size(self%name_tags)

        if ( (self%name_tags(i) == mesh_name) .and. &
             (self%name_ids(i)  == mesh_id) ) then
          ! Do nothing
          return

        else if ( (self%name_tags(i) == mesh_name) .and. &
                  (self%name_ids(i)  /= mesh_id) ) then
          ! Tag name name used for different mesh, flag error
          write(log_scratch_space,'(A,I0,A)')      &
              'Global mesh tag name "'//trim(mesh_name)// &
              '" already assigned in collection.'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )

        end if
      end do

      ! Tag name not used and id exists, so append new entry
      n_tags = size(self%name_tags)
      allocate( new_tag_list(n_tags+1) )
      allocate( new_id_list(n_tags+1)  )

      new_tag_list(:n_tags)  = self%name_tags(:)
      new_id_list(:n_tags )  = self%name_ids(:)
      new_tag_list(n_tags+1) = mesh_name
      new_id_list(n_tags+1)  = mesh_id

      call move_alloc( new_tag_list, self%name_tags )
      call move_alloc( new_id_list, self%name_ids )

    else

      ! No tags assigned yet
      allocate(self%name_tags(1))
      allocate(self%name_ids(1))
      self%name_tags(1) = mesh_name
      self%name_ids(1)  = mesh_id

    end if

    return
  end subroutine update_tags

end module global_mesh_collection_mod
