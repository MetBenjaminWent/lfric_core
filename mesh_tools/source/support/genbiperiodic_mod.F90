!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief   Module to define the genbiperiodic_type, a subclass of the
!>          ugrid_generator_type which generates a biperiodic mesh in a
!>          format suitable for storage as a ugrid file.
!> @details Type implements the ugrid_generator_type interface to
!>          construct a biperiodic mesh.  All required connectivity is
!>          calculated and made availabel to the ugrid writer.
!>
!-------------------------------------------------------------------------------
module genbiperiodic_mod
!-------------------------------------------------------------------------------
  use ugrid_generator_mod,   only : ugrid_generator_type
  use constants_mod,         only : r_def, i_def, str_def, str_long
  use log_mod,               only : log_event, LOG_LEVEL_ERROR
  use reference_element_mod, only : W, S, E, N, SWB, SEB, NWB, NEB

  implicit none

  private

  ! Mesh Vertex directions: local aliases for reference_element_mod values
  integer(i_def), parameter :: NW = NWB
  integer(i_def), parameter :: NE = NEB
  integer(i_def), parameter :: SE = SEB
  integer(i_def), parameter :: SW = SWB

  ! Prefix for error messages
  character(len=*),   parameter :: PREFIX = "[Biperiodic Mesh] "
  character(str_def), parameter :: MESH_CLASS = "plane"

  type, extends(ugrid_generator_type), public :: genbiperiodic_type

    private

    character(str_def)          :: mesh_name
    character(str_def)          :: mesh_class
    character(str_long)         :: generator_inputs
    integer(i_def)              :: nx, ny
    real(r_def)                 :: dx, dy
    integer(i_def), allocatable :: cell_next(:,:)     ! (4, nx*ny)
    integer(i_def), allocatable :: verts_on_cell(:,:) ! (4, nx*ny)
    integer(i_def), allocatable :: edges_on_cell(:,:) ! (4, nx*ny)
    integer(i_def), allocatable :: verts_on_edge(:,:) ! (2, nx*ny)
    real(r_def),    allocatable :: vert_coords(:,:)   ! (2, nx*ny)
  contains
    procedure :: calc_adjacency
    procedure :: calc_face_to_vert
    procedure :: calc_edges
    procedure :: calc_coords
    procedure :: generate
    procedure :: get_metadata
    procedure :: get_dimensions
    procedure :: get_coordinates
    procedure :: get_connectivity
    procedure :: write_mesh
  end type genbiperiodic_type

!-------------------------------------------------------------------------------
  interface genbiperiodic_type
    module procedure genbiperiodic_constructor
  end interface genbiperiodic_type

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------
!> @brief   Constructor for genbiperiodic_type
!> @details Accepts mesh dimension and optional domain size arguments
!>          for initialisation.
!>
!> @param[in] mesh_name  Name of this mesh topology
!> @param[in] nx         Number of faces in biperiodic mesh x axis
!> @param[in] ny         Number of faces in biperiodic mesh y axis
!> @param[in] domain_x   Optional: Domain size in x-axis
!> @param[in] domain_y   Optional: Domain size in y-axis
!>
!> @return    self       Instance of genbiperiodic_type
!-------------------------------------------------------------------------------
function genbiperiodic_constructor( mesh_name, nx, ny, domain_x, domain_y ) &
                            result( self )

  implicit none

  character(len=*), intent(in) :: mesh_name
  integer(i_def),   intent(in) :: nx, ny
  real(r_def),      intent(in) :: domain_x, domain_y

  type( genbiperiodic_type ) :: self
  character(str_def)         :: rchar1
  character(str_def)         :: rchar2

  if (nx < 2 .or. ny < 2) then
    call log_event( PREFIX//"Invalid dimension argument.", LOG_LEVEL_ERROR )
  end if

  self%mesh_name  = trim(mesh_name)
  self%mesh_class = trim(MESH_CLASS)
  self%nx = nx
  self%ny = ny

  if (domain_x < 0)                                                       &
      call log_event( PREFIX//" x-domain argument must be non-negative.", &
                      LOG_LEVEL_ERROR )

  if (domain_y < 0)                                                       &
      call log_event( PREFIX//" y-domain argument must be non-negative.", &
                      LOG_LEVEL_ERROR )

  self%dx = domain_x / self%nx
  self%dy = domain_y / self%ny


  write(rchar1,'(F10.2)') domain_x
  write(rchar2,'(F10.2)') domain_y
  write(self%generator_inputs,'(2(A,I0),2(A))') &
      'edge_cells_x=', self%nx,                 &
      ';edge_cells_y=', self%ny,                &
      ';domain_x='//trim(adjustl(rchar1))//     &
      ';domain_y='//trim(adjustl(rchar2))

  return
end function genbiperiodic_constructor
!-------------------------------------------------------------------------------
!> @brief   For each cell, calculates the set of cells to which it is
!>          adjacent.
!> @details Allocates and populates the instance's cell_next(:,:) array
!>          with the id of each cell to which the index cell is adjacent.
!>
!> @param[in]   self       The genbiperiodic_type instance reference.
!> @param[out]  cell_next  A rank 2 (4,ncells)-sized array containing the
!>                         adjacency map.
!-------------------------------------------------------------------------------
subroutine calc_adjacency(self, cell_next)

  implicit none

  class(genbiperiodic_type),   intent(in)  :: self
  integer(i_def), allocatable, intent(out) :: cell_next(:,:)

  integer(i_def) :: nx, ny, ncells
  integer(i_def) :: cell, astat

  nx = self%nx
  ny = self%ny
  ncells = self%nx * self%ny

  allocate(cell_next(4, ncells), stat=astat)

  if (astat /= 0)                                               &
      call log_event( PREFIX//"Failure to allocate cell_next.", &
                      LOG_LEVEL_ERROR )

  do cell=1, ncells
    ! Cell default values
    cell_next(W, cell) = cell-1
    cell_next(S, cell) = cell+nx
    cell_next(E, cell) = cell+1
    cell_next(N, cell) = cell-nx

    ! Top row
    if (cell <= nx) then
      cell_next(N, cell) = (ny-1)*nx + cell
    end if
    ! Bottom row
    if (cell > ncells-nx) then
      cell_next(S, cell) = cell - (ny-1)*nx
    end if
    ! Left edge
    if (mod(cell, nx) == 1) then
      cell_next(W, cell) = cell + nx-1
    end if
    ! Right edge
    if (mod(cell, nx) == 0) then
      cell_next(E, cell) = cell - (nx-1)
    end if
  end do

end subroutine calc_adjacency

!-------------------------------------------------------------------------------
!> @brief   For each cell, calculates the four vertices whih comprise it.
!> @details Allocates and populates the instance's mesh(:,:) array with
!>          the vertices which form each cell.
!>
!> @param[in]   self           The genbiperiodic_type instance reference.
!> @param[out]  verts_on_cell  A rank 2 (4,ncells)-sized integer array of 
!>                             vertices which constitute each cell.
!-------------------------------------------------------------------------------
subroutine calc_face_to_vert(self, verts_on_cell)

  implicit none

  class(genbiperiodic_type),   intent(in)  :: self
  integer(i_def), allocatable, intent(out) :: verts_on_cell(:,:)

  integer(i_def) :: nx, ny, ncells
  integer(i_def) :: y, vert, cell, nxf, astat

  nx = self%nx
  ny = self%ny
  ncells = self%nx * self%ny

  cell = 1
  nxf  = 1

  allocate(verts_on_cell(4, ncells), stat=astat)

  if (astat /= 0)                                          &
      call log_event( PREFIX//"Failure to allocate mesh.", &
                      LOG_LEVEL_ERROR )

  do vert = 1, 4
    verts_on_cell(vert, cell) = nxf
    nxf = nxf+1
  end do

  ! East neighbour
  verts_on_cell(NW , self%cell_next(E, cell)) = verts_on_cell(NE, cell)
  verts_on_cell(SW , self%cell_next(E, cell)) = verts_on_cell(SE, cell)

  ! South neighbour
  verts_on_cell(NW , self%cell_next(S, cell)) = verts_on_cell(SW, cell)
  verts_on_cell(NE , self%cell_next(S, cell)) = verts_on_cell(SE, cell)

  ! First row
  do cell = 2, nx-1
    verts_on_cell(NE, cell) = nxf
    verts_on_cell(SE, cell) = nxf+1
    nxf = nxf + 2

    ! East neighbour
    verts_on_cell(NW , self%cell_next(E, cell)) = verts_on_cell(NE, cell)
    verts_on_cell(SW , self%cell_next(E, cell)) = verts_on_cell(SE, cell)

    ! South neighbour
    verts_on_cell(NW , self%cell_next(S, cell)) = verts_on_cell(SW, cell)
    verts_on_cell(NE , self%cell_next(S, cell)) = verts_on_cell(SE, cell)
  end do

  ! Inner rows
  do y = 1, ny-2
    ! First cell in row
    cell = y*nx+1
    verts_on_cell(SE, cell) = nxf
    verts_on_cell(SW, cell) = nxf+1
    nxf = nxf+2

    ! South neighbour
    verts_on_cell(NW , self%cell_next(S, cell)) = verts_on_cell(SW, cell)
    verts_on_cell(NE , self%cell_next(S, cell)) = verts_on_cell(SE, cell)

    ! Remainder of row, every other cell
    do cell = y*nx+3, (y+1)*nx-1, 2
      verts_on_cell(SW, cell) = nxf
      verts_on_cell(SE, cell) = nxf+1
      nxf = nxf+2

      ! South neighbour
      verts_on_cell(NW, self%cell_next(S, cell)) = verts_on_cell(SW, cell)
      verts_on_cell(NE, self%cell_next(S, cell)) = verts_on_cell(SE, cell)
    end do
    ! Special case at end of row for odd nx
    if (mod(nx, 2) == 1) then
      cell = (y+1)*nx
      verts_on_cell(SW, cell) = nxf
      nxf = nxf+1

      ! South neighbour
      verts_on_cell(NW, self%cell_next(S, cell)) = verts_on_cell(SW, cell)

      ! West neighbour
      verts_on_cell(SE, self%cell_next(W, cell)) = verts_on_cell(SW, cell)
    end if
  end do

  ! Right edge
  do cell = nx, ncells, nx
    ! Copy from left edge
    verts_on_cell(NE, cell) = verts_on_cell(NW, cell-nx+1)
    verts_on_cell(SE, cell) = verts_on_cell(SW, cell-nx+1)
  end do

  ! Inner step 2 cells
  do y = 1, ny-2
    do cell = y*nx+2, (y+1)*nx, 2
      verts_on_cell(NW, cell) = verts_on_cell(SW, self%cell_next(N, cell))
      verts_on_cell(NE, cell) = verts_on_cell(SE, self%cell_next(N, cell))
      verts_on_cell(SE, cell) = verts_on_cell(SW, self%cell_next(E, cell))
      verts_on_cell(SW, cell) = verts_on_cell(SE, self%cell_next(W, cell))
    end do
  end do

  ! Special case at end of row for odd nx
  if (mod(nx, 2) == 1) then
    do cell = 2*nx, ncells, nx
      verts_on_cell(NW, cell) = verts_on_cell(SW, self%cell_next(N, cell))
    end do
  end if

  ! Last row
  do cell = ncells-nx+1, ncells
    ! Copy from first row
    verts_on_cell(SE, cell) = verts_on_cell(NE, cell-(ny-1)*nx)
    verts_on_cell(SW, cell) = verts_on_cell(NW, cell-(ny-1)*nx)

    ! Copy from N
    verts_on_cell(NW, cell) = verts_on_cell(SW, self%cell_next(N, cell))
    verts_on_cell(NE, cell) = verts_on_cell(SE, self%cell_next(N, cell))
  end do

  return
end subroutine calc_face_to_vert

!-------------------------------------------------------------------------------
!> @brief   Calculates the edges which are found on each cell and the
!>          pair of vertices which are found on each edge.
!> @details Allocates and populates both the edges_on_cell and
!>          verts_on_edge arrays for the instance.
!>
!> @param[in]   self           The genbiperiodic_type instance reference.
!> @param[out]  edges_on_cell  A rank-2 (4,ncells)-sized integer array of
!>                             the edges found on each cell.
!> @param[out]  verts_on_edge  A rank-2 (2,2*ncells)-sized integer array
!>                             of the vertices found on each edge.
!-------------------------------------------------------------------------------
subroutine calc_edges(self, edges_on_cell, verts_on_edge)

  implicit none

  class(genbiperiodic_type),   intent(in)  :: self
  integer(i_def), allocatable, intent(out) :: edges_on_cell(:,:)
  integer(i_def), allocatable, intent(out) :: verts_on_edge(:,:)

  integer(i_def) :: nx, ny, ncells
  integer(i_def) :: cell, nxf, astat

  nx = self%nx
  ny = self%ny
  ncells = self%nx * self%ny

  cell = 1
  nxf  = 1

  allocate(edges_on_cell(4, nx*ny), stat=astat)

  if (astat /= 0)                                                   &
      call log_event( PREFIX//"Failure to allocate edges_on_cell.", &
                      LOG_LEVEL_ERROR )

  allocate(verts_on_edge(2, 2*nx*ny), stat=astat)

  if (astat /= 0)                                                   &
      call log_event( PREFIX//"Failure to allocate verts_on_edge.", &
                      LOG_LEVEL_ERROR )

  ! Top row
  do cell = 1, nx
    ! Top edge
    verts_on_edge(1, nxf) = self%verts_on_cell(NW, cell)
    verts_on_edge(2, nxf) = self%verts_on_cell(NE, cell)
    edges_on_cell(N, cell) = nxf
    ! Right edge
    verts_on_edge(1, nxf+1) = self%verts_on_cell(NE, cell)
    verts_on_edge(2, nxf+1) = self%verts_on_cell(SE, cell)
    edges_on_cell(E, cell) = nxf+1
    ! Bottom edge
    verts_on_edge(1, nxf+2) = self%verts_on_cell(SE, cell)
    verts_on_edge(2, nxf+2) = self%verts_on_cell(SW, cell)
    edges_on_cell(S, cell) = nxf+2
    nxf = nxf+3
  end do

  ! Inner rows
  do cell = nx+1, ncells-nx
    ! Right edge
    verts_on_edge(1, nxf) = self%verts_on_cell(NE, cell)
    verts_on_edge(2, nxf) = self%verts_on_cell(SE, cell)
    edges_on_cell(E, cell) = nxf
    ! Bottom edge
    verts_on_edge(1, nxf+1) = self%verts_on_cell(SE, cell)
    verts_on_edge(2, nxf+1) = self%verts_on_cell(SW, cell)
    edges_on_cell(S, cell) = nxf+1
    nxf = nxf+2
  end do

  ! Bottom row
  do cell = ncells-nx+1, ncells
    ! Right edge
    verts_on_edge(1, nxf) = self%verts_on_cell(NE, cell)
    verts_on_edge(2, nxf) = self%verts_on_cell(SE, cell)
    edges_on_cell(E, cell) = nxf
    nxf = nxf+1
  end do

  ! EoC Top row W copy
  do cell = 1, nx
    edges_on_cell(W, cell) = edges_on_cell(E, self%cell_next(W, cell))
  end do

  ! EoC Inner rows N & W copy
  do cell = nx+1, ncells-nx
    edges_on_cell(N, cell) = edges_on_cell(S, self%cell_next(N, cell))
    edges_on_cell(W, cell) = edges_on_cell(E, self%cell_next(W, cell))
  end do

  ! EoC Bottom row N,S,W copy
  do cell = ncells-nx+1, ncells
    edges_on_cell(N, cell) = edges_on_cell(S, self%cell_next(N, cell))
    edges_on_cell(S, cell) = edges_on_cell(N, self%cell_next(S, cell))
    edges_on_cell(W, cell) = edges_on_cell(E, self%cell_next(W, cell))
  end do

  return
end subroutine calc_edges

!-------------------------------------------------------------------------------
!> @brief   Calculates the coordinates of vertices in the mesh.
!> @details Assigns an (x,y) coordinate in units of dx and dy to each mesh
!>          vertex according to its Cartesian position in the mesh.
!>
!> @param[in]   self         The genbiperiodic_type instance reference.
!> @param[out]  vert_coords  A rank 2 (2,ncells)-sized real array of x and y
!>                           coordinates for each vertex.
!-------------------------------------------------------------------------------
subroutine calc_coords(self, vert_coords)

  implicit none

  class(genbiperiodic_type), intent(in)  :: self
  real(r_def), allocatable,  intent(out) :: vert_coords(:,:)

  integer(i_def) :: ncells, nx, ny
  integer(i_def) :: cell, x, y, astat

  nx = self%nx
  ny = self%ny
  ncells = nx*ny

  allocate(vert_coords(2, ncells), stat=astat)

  if (astat /= 0)                                                 &
      call log_event( PREFIX//"Failure to allocate vert_coords.", &
                      LOG_LEVEL_ERROR )

  do cell = 1, ncells
    y = 1+(cell-1)/nx
    x = cell-(y-1)*nx
    ! x, E/W
    vert_coords(1, self%verts_on_cell(SW, cell)) = real(x-1 - nx/2)*self%dx
    ! y  N/S
    vert_coords(2, self%verts_on_cell(SW, cell)) = real(ny/2 - (y-1))*self%dy
  end do

  return
end subroutine calc_coords

!-------------------------------------------------------------------------------
!> @brief Populates the arguments with the dimensions defining
!>        the biperiodic mesh.
!>
!> @param[in]   self                The genbiperiodic_type instance reference.
!> @param[out]  num_nodes           The number of nodes on the mesh.
!> @param[out]  num_edges           The number of edges on the mesh.
!> @param[out]  num_faces           The number of faces on the mesh.
!> @param[out]  num_nodes_per_face  The number of nodes around each face.
!> @param[out]  num_edges_per_face  The number of edges around each face.
!> @param[out]  num_nodes_per_face  The number of nodes around each edge.
!-------------------------------------------------------------------------------
subroutine get_dimensions( self,                                   &
                           num_nodes, num_edges, num_faces,        &
                           num_nodes_per_face, num_edges_per_face, &
                           num_nodes_per_edge )
  implicit none

  class(genbiperiodic_type), intent(in) :: self

  integer(i_def), intent(out) :: num_nodes
  integer(i_def), intent(out) :: num_edges
  integer(i_def), intent(out) :: num_faces
  integer(i_def), intent(out) :: num_nodes_per_face
  integer(i_def), intent(out) :: num_edges_per_face
  integer(i_def), intent(out) :: num_nodes_per_edge

  num_nodes = self%nx*self%ny
  num_edges = 2*self%nx*self%ny
  num_faces = self%nx*self%ny

  num_nodes_per_face = 4
  num_edges_per_face = 4
  num_nodes_per_edge = 2

  return
end subroutine get_dimensions

!-------------------------------------------------------------------------------
!> @brief   Populates the argument array with the coordinates of the
!>          mesh's vertices.
!> @details Exposes the instance's vert_coords array to the caller.
!>
!> @param[in]   self              The genbiperiodic_type instance reference.
!> @param[out]  node_coordinates  The argument to receive the vert_coords data.
!-------------------------------------------------------------------------------
subroutine get_coordinates(self, node_coordinates)

  implicit none

  class(genbiperiodic_type), intent(in)  :: self
  real(r_def),               intent(out) :: node_coordinates(:,:)

  node_coordinates = self%vert_coords

  return
end subroutine get_coordinates

!-------------------------------------------------------------------------------
!> @brief   Populates the argument arrays with the corresponding mesh
!>          connectivity information.
!> @details Implements the connectivity-providing interface required
!>          by the ugrid writer.
!>
!> @param[in]   self
!> @param[out]  face_node_connectivity  Face-node connectivity.
!> @param[out]  edge_node_connectivity  Edge-node connectivity.
!> @param[out]  face_edge_connectivity  Face-edge connectivity.
!> @param[out]  face_face_connectivity  Face-face connectivity.
!-------------------------------------------------------------------------------
subroutine get_connectivity( self,                   &
                             face_node_connectivity, &
                             edge_node_connectivity, &
                             face_edge_connectivity, &
                             face_face_connectivity )

  implicit none

  class(genbiperiodic_type), intent(in) :: self
  integer(i_def), intent(out) :: face_node_connectivity(:,:)
  integer(i_def), intent(out) :: edge_node_connectivity(:,:)
  integer(i_def), intent(out) :: face_edge_connectivity(:,:)
  integer(i_def), intent(out) :: face_face_connectivity(:,:)

  face_node_connectivity = self%verts_on_cell
  edge_node_connectivity = self%verts_on_edge
  face_edge_connectivity = self%edges_on_cell
  face_face_connectivity = self%cell_next

  return
end subroutine get_connectivity

!-------------------------------------------------------------------------------
!> @brief   Generates the mesh and connectivity.
!> @details Calls each of the instance methods which calculate the
!>          specified mesh and populate the arrays.
!>
!> @param[in,out]  self  The genbiperiodic_type instance reference.
!-------------------------------------------------------------------------------
subroutine generate(self)

  implicit none

  class(genbiperiodic_type), intent(inout)         :: self

  call calc_adjacency(self, self%cell_next)
  call calc_face_to_vert(self, self%verts_on_cell)
  call calc_edges(self, self%edges_on_cell, self%verts_on_edge)
  call calc_coords(self, self%vert_coords)

  return
end subroutine generate

!-------------------------------------------------------------------------------
!> @brief Writes out the mesh and connectivity for debugging purposes.
!>
!> @param[in,out]  self  The genbiperiodic_type instance reference.
!-------------------------------------------------------------------------------
subroutine write_mesh(self)

  use iso_fortran_env,     only : stdout => output_unit

  implicit none

  class(genbiperiodic_type), intent(in) :: self

  integer(i_def) :: i, cell, ncells

  ncells = self%nx * self%ny

  write(stdout,*) "cell_next"
  do cell=1, ncells
      write(stdout,"(I3,T8,4I4)") cell, self%cell_next(:,cell)
  end do

  write(stdout,*)
  write(stdout,*) "verts_on_cell"
  do cell=1, ncells
    write(stdout,"(I3,T8,4I4)") cell, self%verts_on_cell(:,cell)
  end do

  write(stdout,*)
  write(stdout,*) "verts_on_edge"
  do i=1, size(self%verts_on_edge, 2)
    write(stdout,"(I3,T8,2I4)") i, self%verts_on_edge(:,i)
  end do

  write(stdout,*)
  write(stdout,*) "edges_on_cell"
  do cell=1, ncells
    write(stdout,"(I3,T8,4I4)") cell, self%edges_on_cell(:,cell)
  end do

  write(*,*)
  write(*,*) "vert_coords"
  do cell=1, ncells
    write(*,"(I3,T8,2F11.4)") cell, self%vert_coords(:,cell)
  end do

  return
end subroutine write_mesh

!-----------------------------------------------------------------------------
!> @brief Returns mesh metadata information.
!>
!> @param[in]   self              The generator strategy object.
!> @param[out]  mesh_name         The name of this mesh instance
!> @param[out]  mesh_class        Primitive shape, i.e. plane
!> @param[out]  generator_inputs  Inputs used to create this mesh from the
!>                                mesh_generator
!-----------------------------------------------------------------------------
subroutine get_metadata( self, mesh_name, mesh_class, generator_inputs )

  implicit none

  class(genbiperiodic_type), intent(in) :: self
  character(str_def),  intent(out) :: mesh_name
  character(str_def),  intent(out) :: mesh_class
  character(str_long), intent(out) :: generator_inputs

  mesh_name        = self%mesh_name
  mesh_class       = self%mesh_class
  generator_inputs = trim(self%generator_inputs)

  return
end subroutine  get_metadata
!-------------------------------------------------------------------------------
end module genbiperiodic_mod
