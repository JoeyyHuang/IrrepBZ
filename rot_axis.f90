subroutine rot_axis(mat,angle,axis)
!=========================================================================
! find the rotational axis (in Crys) based on the given matrix (in Crys) and angle
!=========================================================================
  use variables, only : a_prim
  use mat_trans, only : mat_det, mat_inv
  use constants, only : tpi, err_tol, err_tol3
  implicit none
!-------------------------------------------------------------------------
! arguments
  real(kind=8),intent(in) :: mat(3,3)
  real(kind=8),intent(in) :: angle
  integer(kind=4),intent(out) :: axis(3)
! local variables
  integer(kind=4) :: ia, ib, ic
  integer(kind=4) :: diff(3)
  real(kind=8) :: mat2(3,3), rot(3,3), n(3), vec(3)
!-------------------------------------------------------------------------
  if(mat_det(mat) > 0.d0) mat2 = mat
  if(mat_det(mat) < 0.d0) mat2 = -mat
  ! unit rotational vector in Cart.
  rot = matmul(a_prim,matmul(mat2,mat_inv(a_prim))) ! [e]^-1[R0][e]=[R]  R0: rotation in Cart frame
  !
  if(abs(angle-tpi/2.d0) < err_tol) then     ! pi angle
    if(abs(rot(1,1)+1.d0) > err_tol) then
      n(1) = sqrt(abs(rot(1,1)+1.d0)/2.d0)
      n(2) = rot(1,2)/(2.d0*n(1))
      n(3) = rot(1,3)/(2.d0*n(1))
    else
      n(1) = 0.d0
      if(abs(rot(2,2)+1.d0) > err_tol) then
        n(2) = sqrt(abs(rot(2,2)+1.d0)/2.d0)
        n(3) = rot(2,3)/(2.d0*n(2))
      else
        n(2) = 0.d0
        n(3) = sqrt(abs(rot(3,3)+1.d0)/2.d0)
      end if
    end if
  else   ! non-pi angle
    n(1) = (rot(3,2)-rot(2,3)) / (2.d0*sin(angle))
    n(2) = (rot(1,3)-rot(3,1)) / (2.d0*sin(angle))
    n(3) = (rot(2,1)-rot(1,2)) / (2.d0*sin(angle))
  end if
  !
  ! corresponding axis in crys frame
  axis = 0
  do ia = -3, 3
    do ib = -3, 3
      do ic = -3, 3
        if(abs(ia)+abs(ib)+abs(ic) /= 0) then
          vec(:) = ia*a_prim(:,1) + ib*a_prim(:,2) + ic*a_prim(:,3)
          vec = vec / sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
          if(all(abs(n-vec) < err_tol3)) then
            axis = (/ ia, ib, ic /)
            goto 100
          end if
        end if
      end do
    end do
  end do
  !
100 continue
  ! convert to coprime integers
  diff = 2*(axis/2) - axis
  if(all(abs(dble(diff)) < err_tol3)) then
    axis = axis/2
    return
  end if
  diff = 3*(axis/3) - axis
  if(all(abs(dble(diff)) < err_tol3)) then
    axis = axis/3
    return
  end if
  !
  return
end subroutine rot_axis
  
        