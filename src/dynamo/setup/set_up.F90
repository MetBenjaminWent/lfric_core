!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief Sets up data required for the the function spaces

!> @details This code generates a mesh and determines the basis functions and
!> dofmaps. This will be replaced with code that reads this in from a mesh
!> generation and paritioning pre-processor stage.  

! There are no tests for this code as this will be replaced.

module set_up_mod

  use constants_mod,              only : dp
  use function_space_mod,         only : function_space_type
  use reference_element_mod,      only : reference_cube

  use mesh_generator_mod,         only : mesh_generator_init,  &
                                   mesh_generator_cubedsphere, &
                                   mesh_generator_biperiodic,  &
                                   mesh_connectivity
  use num_dof_mod,                only : num_dof_init
  use basis_function_mod,         only : get_basis, &
              v0_basis, v1_basis, v2_basis, v3_basis, &
              v0_diff_basis, v1_diff_basis, v2_diff_basis, v3_diff_basis, &
              v0_nodal_coords, v1_nodal_coords, v2_nodal_coords, v3_nodal_coords

  use dofmap_mod,                 only : get_dofmap, &
              v0_dofmap, v1_dofmap, v2_dofmap, v3_dofmap
  use gaussian_quadrature_mod,    only : ngp_h, ngp_v
  implicit none
contains 

!> Generates a mesh and determines the basis functions and dofmaps (this will 
!> be replaced with code that reads the information in)
  subroutine set_up( )

    use log_mod, only : log_event, LOG_LEVEL_INFO
    use mesh_mod, only : num_cells, num_layers, element_order, l_spherical, &
                         v_unique_dofs,v_dof_entity

    implicit none

    real(kind=dp), parameter                 :: delta=1.0_dp
    character(len = 100)                     :: filename

    ! hard-coded these numbers are
    num_cells = 9
    num_layers = 3
    element_order = 0
    l_spherical = .false.
    filename = '../../data/Cubegrid.dat' 
    call log_event( "set_up: generating/reading the mesh", LOG_LEVEL_INFO )

!  ----------------------------------------------------------
!  Mesh generation, really a preprocessor step for reading
! -----------------------------------------------------------

    ! Setup reference cube  
    call reference_cube()
    ! Initialise mesh
    call mesh_generator_init(num_cells,num_layers)
    ! Genereate mesh  
    if ( l_spherical ) then
       call mesh_generator_cubedsphere(filename,num_cells,num_layers,delta)
    else
       call mesh_generator_biperiodic(num_cells,3,3,num_layers,delta,delta,delta)
    end if
    ! Extend connectivity ( cells->faces, cells->edges )  
    call mesh_connectivity(num_cells)    


! -----------------------------------------------------------
! Initialise FE elements on the mesh constructed above
! really another pre-processor step
! ----------------------------------------------------------

    ! initialise numbers of dofs    
    call num_dof_init(num_cells,num_layers,element_order,v_unique_dofs,v_dof_entity)

    call log_event( "set_up: computing basis functions", LOG_LEVEL_INFO )

    ! read the values of the basis functions. 
    call get_basis( k=element_order, &
                    v_unique_dofs=v_unique_dofs,v_dof_entity=v_dof_entity )  

    call log_event( "set_up: computing the dof_map", LOG_LEVEL_INFO )
    ! compute the dof maps for each function space
    call get_dofmap(nlayers=num_layers,v_dof_entity=v_dof_entity, &
                    ncell=num_cells,v_unique_dofs=v_unique_dofs)

    return

  end subroutine set_up

end module set_up_mod
