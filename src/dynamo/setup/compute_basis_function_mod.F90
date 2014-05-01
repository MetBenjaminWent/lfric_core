!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Computes the basis functions on quadrature points for 4 element spaces 
!-------------------------------------------------------------------------------
module compute_basis_function_mod

  use num_dof_mod
  use reference_element_mod

  use constants_mod, only: dp
  use gaussian_quadrature_mod, only: gaussian_quadrature_type, &
                                     ngp ! parameter for how many GQ points
  use function_space_mod, only : function_space_type

  implicit none

  real(kind=dp) :: xgp(ngp)

contains 

  subroutine compute_basis(k,v0,v1,v2,v3, v_unique_dofs,v_dof_entity)
    !-----------------------------------------------------------------------------
    ! Subroutine to compute test/trial functions on quadrature points
    !-----------------------------------------------------------------------------

    ! order of elements
    integer, intent(in) :: k
    type(function_space_type), intent(inout) :: v0, v1, v2, v3
    integer, intent(in) :: v_unique_dofs(4,2), v_dof_entity(4,0:3)

    integer :: i, jx, jy, jz, order, idx, j1, j2
    integer :: j(3), j2l_edge(12,3), j2l_face(6,3), face_idx(6), edge_idx(12,2)
    integer, allocatable :: lx(:), ly(:), lz(:)
    real(kind=dp)    :: fx, fy, fz, gx, gy, gz, dfx, dfy, dfz
    real(kind=dp)    :: x1(k+2), x2(k+1)
    !real(kind=dp)    :: unit_vec_v2(nv2,3), unit_vec_v1(nv1,3)
    real(kind=dp),allocatable    :: unit_vec_v2(:,:), unit_vec_v1(:,:)


    ! Allocate to be larger than should be needed
    allocate ( lx(3*(k+2)**3) )
    allocate ( ly(3*(k+2)**3) )
    allocate ( lz(3*(k+2)**3) )

    allocate(unit_vec_v2(v_unique_dofs(3,2),3))
    allocate(unit_vec_v1(v_unique_dofs(2,2),3))

    ! positional arrays - need two, i.e quadratic and linear for RT1
    do i=1,k+2
      x1(i) = real(i-1)/real(k+1)
    end do
    do i=1,k+1
      x2(i) = real(i-1)/real(k)
    end do
    if ( k == 0 ) x2(1) = 0.5

    ! some look arrays based upon reference cube topology
    face_idx = (/ 1, k+2, k+2, 1, 1, k+2 /)

    edge_idx(:,1) = (/ 1, k+2, k+2, 1, 1, k+2, k+2, 1,   1,   k+2, k+2, 1   /)
    edge_idx(:,2) = (/ 1, 1,   1,   1, 1, 1,   k+2, k+2, k+2, k+2, k+2, k+2 /)

    j2l_face(1,:) = (/ 2, 3, 1 /)
    j2l_face(2,:) = (/ 3, 2, 1 /)
    j2l_face(3,:) = (/ 2, 3, 1 /)
    j2l_face(4,:) = (/ 3, 2, 1 /)
    j2l_face(5,:) = (/ 2, 1, 3 /)
    j2l_face(6,:) = (/ 2, 1, 3 /)

    j2l_edge(1 ,:) = (/ 1, 2, 3 /)
    j2l_edge(2 ,:) = (/ 2, 1, 3 /)
    j2l_edge(3 ,:) = (/ 1, 2, 3 /)
    j2l_edge(4 ,:) = (/ 2, 1, 3 /)
    j2l_edge(5 ,:) = (/ 2, 3, 1 /)
    j2l_edge(6 ,:) = (/ 2, 3, 1 /)
    j2l_edge(7 ,:) = (/ 2, 3, 1 /)
    j2l_edge(8 ,:) = (/ 2, 3, 1 /)
    j2l_edge(9 ,:) = (/ 1, 2, 3 /)
    j2l_edge(10,:) = (/ 2, 1, 3 /)
    j2l_edge(11,:) = (/ 1, 2, 3 /)
    j2l_edge(12,:) = (/ 2, 1, 3 /)

    !-----------------------------------------------------------------------------
    ! Section for test/trial functions of v0 fields
    !-----------------------------------------------------------------------------
    order = k+1

    ! compute indices of functions
    idx = 1

    ! dofs in volume
    do jz=2,k+1
      do jy=2,k+1
        do jx=2,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          idx = idx + 1
        end do
      end do
    end do

    ! dofs on faces
    do i=1,nfaces
      do j1=2,k+1
        do j2=2,k+1 
          j(1) = j1
          j(2) = j2
          j(3) = face_idx(i)
          lx(idx) = j(j2l_face(i,1))
          ly(idx) = j(j2l_face(i,2))
          lz(idx) = j(j2l_face(i,3))
          idx = idx + 1
        end do
      end do
    end do

    ! dofs on edges
    do i=1,nedges
      do j1=2,k+1
        j(1) = j1
        j(2) = edge_idx(i,1)
        j(3) = edge_idx(i,2)
        lx(idx) = j(j2l_edge(i,1))
        ly(idx) = j(j2l_edge(i,2))
        lz(idx) = j(j2l_edge(i,3))
        idx = idx + 1
      end do
    end do

    ! dofs on vertices
    do i=1,nverts
    !  do j1=1,nv0_vert
      do j1=1,v_dof_entity(1,0)
        lx(idx) =  1+(k+1)*int(x_vert(i,1))
        ly(idx) =  1+(k+1)*int(x_vert(i,2))
        lz(idx) =  1+(k+1)*int(x_vert(i,3))
        idx = idx + 1
      end do
    end do
    do i=1,v_unique_dofs(1,2)
    !do i=1,nv0
      do jx=1,ngp
        fx = poly1d(order,xgp(jx),x1(lx(i)),x1,lx(i))
        dfx = poly1d_deriv(order,xgp(jx),x1(lx(i)),x1,lx(i))
        do jy=1,ngp
          fy = poly1d(order,xgp(jy),x1(ly(i)),x1,ly(i))
          dfy = poly1d_deriv(order,xgp(jy),x1(ly(i)),x1,ly(i))
          do jz=1,ngp
            fz = poly1d(order,xgp(jz),x1(lz(i)),x1,lz(i))
            dfz = poly1d_deriv(order,xgp(jz),x1(lz(i)),x1,lz(i))
            call v0%set_basis(fx*fy*fz,i,jx,jy,jz,1)
            call v0%set_diff_basis(dfx*fy*fz,i,jx,jy,jz,1)
            call v0%set_diff_basis(fx*dfy*fz,i,jx,jy,jz,2)
            call v0%set_diff_basis(fx*fy*dfz,i,jx,jy,jz,3)
          end do
        end do
      end do
    end do

    !-----------------------------------------------------------------------------
    ! section for test/trial functions of v1 fields
    !-----------------------------------------------------------------------------
    order = k+1

    !do idx=1,nv1
    do idx = 1,v_unique_dofs(2,2)
      do i=1,3
        unit_vec_v1(idx,i) = 0.0
      end do
    end do

    ! compute indices of functions
    idx = 1

    ! dofs in volume
    ! u components
    do jz=2,k+1
      do jy=2,k+1
        do jx=1,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          unit_vec_v1(idx,1) = 1.0
          idx = idx + 1
        end do
      end do
    end do
    ! v components
    do jz=2,k+1
      do jy=1,k+1
        do jx=2,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          unit_vec_v1(idx,2) = 1.0
          idx = idx + 1
        end do
      end do
    end do
    ! w components
    do jz=1,k+1
      do jy=2,k+1
        do jx=2,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          unit_vec_v1(idx,3) = 1.0
          idx = idx + 1
        end do
      end do
    end do

    ! dofs on faces
    do i=1,nfaces
      do j1=1,k+1
        do j2=2,k+1
          j(1) = j1
          j(2) = j2
          j(3) = face_idx(i)
          lx(idx) = j(j2l_face(i,1))
          ly(idx) = j(j2l_face(i,2))
          lz(idx) = j(j2l_face(i,3))
          unit_vec_v1(idx,:) = tangent_to_edge(edge_on_face(i,1),:)
          idx = idx + 1
        end do
      end do
      do j1=2,k+1
        do j2=1,k+1
          j(1) = j1
          j(2) = j2
          j(3) = face_idx(i)
          lx(idx) = j(j2l_face(i,1))
          ly(idx) = j(j2l_face(i,2))
          lz(idx) = j(j2l_face(i,3))
          unit_vec_v1(idx,:) = tangent_to_edge(edge_on_face(i,2),:)
          idx = idx + 1
        end do
      end do  
    end do

    ! dofs on edges
    do i=1,nedges
      do j1=1,k+1
        j(1) = j1
        j(2) = edge_idx(i,1)
        j(3) = edge_idx(i,2)
        lx(idx) = j(j2l_edge(i,1))
        ly(idx) = j(j2l_edge(i,2))
        lz(idx) = j(j2l_edge(i,3))
        unit_vec_v1(idx,:) = tangent_to_edge(i,:)
        idx = idx + 1
      end do
    end do

    ! this needs correcting
    !do i=1,nv1  
    do i=1,v_unique_dofs(2,2)
      do jx=1,ngp
        fx = poly1d(order,xgp(jx),x1(lx(i)),x1,lx(i))
        dfx = poly1d_deriv(order,xgp(jx),x1(lx(i)),x1,lx(i))
        if (lx(i) <= order) then
          gx = poly1d(order-1,xgp(jx),x2(lx(i)),x2,lx(i))
        else
          gx = 0.0
        end if
        do jy=1,ngp
          fy = poly1d(order,xgp(jy),x1(ly(i)),x1,ly(i))
          dfy = poly1d_deriv(order,xgp(jy),x1(ly(i)),x1,ly(i))
          if (ly(i) <= order) then 
            gy = poly1d(order-1,xgp(jy),x2(ly(i)),x2,ly(i))
          else
            gy = 0.0
          end if
          do jz=1,ngp
            fz = poly1d(order,xgp(jz),x1(lz(i)),x1,lz(i))
            dfz = poly1d_deriv(order,xgp(jz),x1(lz(i)),x1,lz(i))
            if (lz(i) <= order) then     
              gz = poly1d(order-1,xgp(jz),x2(lz(i)),x2,lz(i))
            else
              gz = 0.0
            end if
            
            call v1%set_basis(gx*fy*fz*unit_vec_v1(i,1),i,jx,jy,jz,1)
            call v1%set_basis(fx*gy*fz*unit_vec_v1(i,2),i,jx,jy,jz,2)
            call v1%set_basis(fx*fy*gz*unit_vec_v1(i,3),i,jx,jy,jz,3)

            call v1%set_diff_basis((fx*dfy*gz*unit_vec_v1(i,3) - fx*gy*dfz*unit_vec_v1(i,2)) ,i,jx,jy,jz,1)
            call v1%set_diff_basis((gx*fy*dfz*unit_vec_v1(i,1) - dfx*fy*gz*unit_vec_v1(i,3)) ,i,jx,jy,jz,2)
            call v1%set_diff_basis((dfx*gy*fz*unit_vec_v1(i,2) - gx*dfy*fz*unit_vec_v1(i,1)), i,jx,jy,jz,3)
          end do
        end do
      end do
    end do


    !-----------------------------------------------------------------------------
    ! Section for test/trial functions of v2 fields
    !-----------------------------------------------------------------------------
    order = k + 1

    !do idx=1,nv2
    do idx=1,v_unique_dofs(3,2)
      do i=1,3
        unit_vec_v2(idx,i) = 0.0
      end do
    end do

    idx = 1
    ! dofs in volume
    ! u components
    do jz=1,k+1
      do jy=1,k+1
        do jx=2,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          unit_vec_v2(idx,1) = 1.0
          idx = idx + 1
        end do
      end do
    end do
    ! v components
    do jz=1,k+1
      do jy=2,k+1
        do jx=1,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          unit_vec_v2(idx,2) = 1.0
          idx = idx + 1
        end do
      end do
    end do
    ! w components
    do jz=2,k+1
      do jy=1,k+1
        do jx=1,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          unit_vec_v2(idx,3) = 1.0
          idx = idx + 1
        end do
      end do
    end do

    ! dofs on faces
    do i=1,nfaces
      do j1=1,k+1
        do j2=1,k+1 
          j(1) = j1
          j(2) = j2
          j(3) = face_idx(i)
          lx(idx) = j(j2l_face(i,1))
          ly(idx) = j(j2l_face(i,2))
          lz(idx) = j(j2l_face(i,3))
          unit_vec_v2(idx,:) = normal_to_face(i,:)
          idx = idx + 1
        end do
      end do
    end do

    !do i=1,nv2
    do i=1,v_unique_dofs(3,2)
      do jx=1,ngp
        fx = poly1d(order,xgp(jx),x1(lx(i)),x1,lx(i))
        dfx = poly1d_deriv(order,xgp(jx),x1(lx(i)),x1,lx(i))
        if (lx(i) <= order) then
          gx = poly1d(order-1,xgp(jx),x2(lx(i)),x2,lx(i))
        else
          gx = 0.0
        end if
        do jy=1,ngp
          fy = poly1d(order,xgp(jy),x1(ly(i)),x1,ly(i))
          dfy = poly1d_deriv(order,xgp(jy),x1(ly(i)),x1,ly(i))
          if (ly(i) <= order) then
            gy = poly1d(order-1,xgp(jy),x2(ly(i)),x2,ly(i))
          else
            gy = 0.0
          end if       
          do jz=1,ngp
            fz = poly1d(order,xgp(jz),x1(lz(i)),x1,lz(i))
            dfz = poly1d_deriv(order,xgp(jz),x1(lz(i)),x1,lz(i))
            if (lz(i) <= order) then
              gz = poly1d(order-1,xgp(jz),x2(lz(i)),x2,lz(i))
            else
              gz = 0.0
            end if
            
            call v2%set_basis(fx*gy*gz*unit_vec_v2(i,1),i,jx,jy,jz,1)
            call v2%set_basis(gx*fy*gz*unit_vec_v2(i,2),i,jx,jy,jz,2)
            call v2%set_basis(gx*gy*fz*unit_vec_v2(i,3),i,jx,jy,jz,3)
            
            call v2%set_diff_basis((dfx*gy*gz*unit_vec_v2(i,1) & 
                                      + gx*dfy*gz*unit_vec_v2(i,2) &
                                      + gx*gy*dfz*unit_vec_v2(i,3) ), &
                                      i,jx,jy,jz,1) 
          end do
        end do
      end do
    end do

    !-----------------------------------------------------------------------------
    ! Section for test/trial functions of v3 fields
    !-----------------------------------------------------------------------------
    order = k
    ! compute indices of functions
    idx = 1
    ! dofs in volume
    do jz=1,k+1
      do jy=1,k+1
        do jx=1,k+1
          lx(idx) =  jx
          ly(idx) =  jy
          lz(idx) =  jz
          idx = idx + 1
        end do
      end do
    end do

    !do i=1,nv3
    do i=1,v_unique_dofs(4,2)
      do jx=1,ngp
        gx = poly1d(order,xgp(jx),x2(lx(i)),x2,lx(i))
        do jy=1,ngp
          gy = poly1d(order,xgp(jy),x2(ly(i)),x2,ly(i))
          do jz=1,ngp
            gz = poly1d(order,xgp(jz),x2(lz(i)),x2,lz(i))
            call v3%set_basis(gx*gy*gz,i,jx,jy,jz,1)
          end do
        end do
      end do
    end do

    ! tidy up
    deallocate ( lx )
    deallocate ( ly )
    deallocate ( lz )

    deallocate ( unit_vec_v2, unit_vec_v1)

  end subroutine compute_basis


  function poly1d(order,xk,xl,x,l)
  !-----------------------------------------------------------------------------
  ! evaluate 1D basis function of arbitrary order at xk
  !-----------------------------------------------------------------------------

  implicit none

  real(kind=dp) :: poly1d
  ! Order of basis function
  integer, intent(in) :: order
  ! Index of basis function
  integer, intent(in) :: l
  ! quadrature point to evaluate basis function at
  real(kind=dp), intent(in) :: xk
  ! Point function is unity at
  real(kind=dp), intent(in) :: xl
  ! grid points
  real(kind=dp), intent(in) :: x(order+1)

  integer :: j

  poly1d = 1.0

  do j=1,l-1
    poly1d = poly1d*(xk-x(j))/(xl-x(j))
  end do
  do j=l+1,order+1
    poly1d = poly1d*(xk-x(j))/(xl-x(j))
  end do

  end function poly1d


  function poly1d_deriv(order,xk,xl,x,l)
  !-----------------------------------------------------------------------------
  ! evaluate derivative of 1D basis function of arbitrary order at xk
  !-----------------------------------------------------------------------------
  implicit none

  real(kind=dp) :: poly1d_deriv
  ! Order of basis function
  integer, intent(in) :: order
  ! Index of basis function
  integer, intent(in) :: l
  ! quadrature point to evaluate basis function at
  real(kind=dp), intent(in) :: xk
  ! Point function is unity at
  real(kind=dp), intent(in) :: xl
  ! grid points
  real(kind=dp), intent(in) :: x(order+1)

  real(kind=dp) :: denom,t

  integer :: k,j


  poly1d_deriv = 0.0
  denom = 1.0

  do j=1,l-1
    denom = denom * 1.0/(xl-x(j))  
  end do
  do j=l+1,order+1
    denom = denom * 1.0/(xl-x(j))  
  end do

  do k=1,l-1
    t = 1.0
    do j=1,order+1
      if (j .ne. l .and. j .ne. k ) then
        t = t * (xk - x(j))
      end if
    end do
    poly1d_deriv = poly1d_deriv + t*denom
  end do  
  do k=l+1,order+1
    t = 1.0
    do j=1,order+1
      if (j .ne. l .and. j .ne. k ) then
        t = t * (xk - x(j))
      end if
    end do
    poly1d_deriv = poly1d_deriv + t*denom
  end do 

  end function poly1d_deriv

end module compute_basis_function_mod
