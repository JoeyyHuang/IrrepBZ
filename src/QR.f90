subroutine QR(n,H,Q,R)
!=========================================================================
! Perform QR factorization on the matrix H with order n
!=========================================================================
  implicit none
!-------------------------------------------------------------------------
! arguments
  integer(kind=4),intent(in) :: n
  complex(kind=8),intent(in) :: H(n,n)
  complex(kind=8),intent(out) :: Q(n,n), R(n,n)
! local variables
  integer(kind=4) :: i, j
  complex(kind=8) :: x(n,n)
!-------------------------------------------------------------------------
  ! Gram-Schmidt orthogonalization
  x(:,1) = H(:,1)
  do i = 2, n
    x(:,i) = H(:,i)
    do j = 1, i-1
      ! be careful of the order of two elements in dot_product
      x(:,i) = x(:,i) - dot_product(x(:,j),H(:,i))/dot_product(x(:,j),x(:,j))*x(:,j)
    end do
  end do
  ! Q
  do i = 1, n
    Q(:,i) = x(:,i)/sqrt(dot_product(x(:,i),x(:,i)))
  end do
  ! R
  R = (0.d0,0.d0)
  do i = 1, n 
    R(i,i) = sqrt(dot_product(x(:,i),x(:,i)))
  end do
  do i = 2, n
    do j = 1, i-1
      R(j,i) = dot_product(x(:,j),H(:,i))/dot_product(x(:,j),x(:,j))*sqrt(dot_product(x(:,j),x(:,j)))
    end do
  end do
  !
  return
end subroutine QR