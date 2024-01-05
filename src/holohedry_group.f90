subroutine holohedry_group()
!=========================================================================
! get the operations of holohedry (lattice point group) and then obtain 
! the crystal system based on the number of characteristic rotational axes
!=========================================================================
  use variables, only : a_prim, holo_ele, holo_ord, crys_sys
  use mat_trans, only : mat_inv, mat_det
  use constants, only : tol, tol33, stdout
  implicit none
!-------------------------------------------------------------------------
! local variables
  integer(kind=4) :: ia, ib, ic, i, j, k, nvec, rot2, rot3, rot4, rot6
  real(kind=8) :: la1, la2, la3, la_max, lvec, trace
  real(kind=8) :: vec(3), vec_mat(3,3), vec_list(3,343), rot_mat(3,3,48)
!-------------------------------------------------------------------------
! lengths of basis vectors
  la1 = sqrt(a_prim(1,1)**2+a_prim(2,1)**2+a_prim(3,1)**2)
  la2 = sqrt(a_prim(1,2)**2+a_prim(2,2)**2+a_prim(3,2)**2)
  la3 = sqrt(a_prim(1,3)**2+a_prim(2,3)**2+a_prim(3,3)**2)
  la_max = max(la1,la2,la3)
  !
  nvec = 0
  vec_list = 0.d0   ! (coords,counter)
  !
  ! pick out lattice vectors with lengths less than or equal to the longest basis vector
  do ia = -3, 3
    do ib = -3, 3
      do ic = -3, 3
        vec(:) = ia*a_prim(:,1) + ib*a_prim(:,2) + ic*a_prim(:,3)
        lvec = sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
        if(lvec-la_max < tol) then
          nvec = nvec + 1
          vec_list(:,nvec) = vec(:)
        end if
      end do
    end do
  end do
  !
  ! require vi.vj = ti.tj for every i, j
  do i = 1, nvec
    do j = 1, nvec
      do k = 1, nvec
        vec_mat = reshape((/vec_list(:,i),vec_list(:,j),vec_list(:,k)/),(/3,3/))
        if(all(abs(matmul(transpose(vec_mat),vec_mat)-matmul(transpose(a_prim),a_prim)) < tol33)) then
          ! [x']=[R][x],[e']=[e][M], since [e'][x']=[e][x] -> MR=I
          holo_ord = holo_ord + 1
          rot_mat(:,:,holo_ord) = matmul(mat_inv(vec_mat),a_prim)
          !
        end if
      end do
    end do
  end do
  !
  allocate(holo_ele(3,3,holo_ord))
  holo_ele(:,:,1:holo_ord) = rot_mat(:,:,1:holo_ord)
  !
  ! determine the crystal system based on the determinant and trace
  rot2 = 0   ! 2-fold rotation
  rot3 = 0   ! 3-fold rotation
  rot4 = 0   ! 4-fold rotation
  rot6 = 0   ! 6-fold rotation
  do i = 1, holo_ord
    ! only need considering proper rotations
    if(mat_det(holo_ele(:,:,i)) > 0.d0) then
      trace = holo_ele(1,1,i) + holo_ele(2,2,i) + holo_ele(3,3,i)
      if(abs(trace+1.d0) < tol) rot2 = rot2 + 1
      if(abs(trace) < tol)      rot3 = rot3 + 1
      if(abs(trace-1.d0) < tol) rot4 = rot4 + 1
      if(abs(trace-2.d0) < tol) rot6 = rot6 + 1
    end if
  end do
  !
  ! do not disrupt the order
  if(rot2 == 1) crys_sys = 'Monoclinic'     ! single 2-fold rotation
  if(rot2 == 3) crys_sys = 'Orthorhombic'   ! three 2-fold rotation
  if(rot4 == 2) crys_sys = 'Tetragonal'     ! single 4-fold rotation
  if(rot3 == 2) crys_sys = 'Trigonal'       ! single 3-fold rotation
  if(rot6 == 2) crys_sys = 'Hexagonal'      ! single 6-fold rotation
  if(rot3 == 8) crys_sys = 'Cubic'          ! four 3-fold rotation
  !
  write (stdout,'(a,a)') 'The crystal system of the structure is: ', trim(adjustl(crys_sys))
  !
  return
end subroutine holohedry_group