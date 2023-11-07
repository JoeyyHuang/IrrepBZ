subroutine small_rep(k_crys, lit_ord, lit_ele, lit_symb, lit_axis, nirr, lit_irr, lit_cha, spinor)
!=========================================================================
! get the left regular representation of the given little co-group and further 
! construct the irreducible representation (small representation) by 
! numerically diagonalizing a random matrix commute with regular representation
!=========================================================================
  use variables, only : stdout
  use constants, only : im, tpi, err_tol, err_tol3, err_tol22, err_tol33
  implicit none
!-------------------------------------------------------------------------
! arguments
  integer(kind=4),intent(in) :: lit_ord
  integer(kind=4),intent(in) :: lit_axis(3,lit_ord)
  real(kind=8),intent(in) :: k_crys(3), lit_ele(3,4,lit_ord)
  character(len=2),intent(in) :: lit_symb(lit_ord)
  logical,intent(in) :: spinor
  integer(kind=4),intent(out) :: nirr
  complex(kind=8),intent(out) :: lit_cha(lit_ord,lit_ord)   ! #element,#nirrep
  complex(kind=8),intent(out) :: lit_irr(nint(sqrt(dble(lit_ord))),nint(sqrt(dble(lit_ord))),lit_ord,lit_ord)
! local variables
  integer(kind=4) :: i, j, k, l, m, n, p, ia, ib, ic
  real(kind=8) :: eva_tem
  real(kind=8) :: k_crys2(3), fac_spin(lit_ord,lit_ord), g(3,lit_ord), rmat(lit_ord,lit_ord), imat(lit_ord,lit_ord), &
                  err_toln(lit_ord), eva(lit_ord)
  complex(kind=8) :: fac_sys
  complex(kind=8) :: rand_mat(lit_ord,lit_ord), rand_herm(lit_ord,lit_ord), comm_herm(lit_ord,lit_ord), &
                     cha_tem(lit_ord), rep_tem(nint(sqrt(dble(lit_ord))),nint(sqrt(dble(lit_ord))), lit_ord), &
                     lit_reg(lit_ord,lit_ord,lit_ord), evc(lit_ord,lit_ord), evc_tem(lit_ord), &
                     Q(lit_ord,lit_ord), R(lit_ord,lit_ord)
  logical :: is_diag
!-------------------------------------------------------------------------
  if(spinor) call SU2(lit_ord, lit_symb, lit_axis, lit_ele, fac_spin)
  !
  ! left regular projective representation
  lit_reg = (0.d0,0.d0)
  do i = 1, lit_ord
    ! k*alpha_i = k + G_i
    k_crys2 = matmul(k_crys,lit_ele(:,1:3,i))
    do ia = -3, 3
      do ib = -3, 3
        do ic = -3, 3
          if(all(abs(k_crys2-k_crys-(/ia,ib,ic/)) < err_tol3)) g(:,i) = dble((/ia,ib,ic/))
        end do
      end do
    end do
    !
    do j = 1, lit_ord
      fac_sys = exp(-im*tpi*dot_product(g(:,i),lit_ele(:,4,j)))   ! factor system e^(-i*g_i*tau_j)
      do k = 1, lit_ord
        if(all(abs(matmul(lit_ele(:,1:3,i),lit_ele(:,1:3,j))-lit_ele(:,1:3,k)) < err_tol33)) then
          if(spinor) lit_reg(k,j,i) = fac_sys*fac_spin(i,j)
          if(.NOT. spinor) lit_reg(k,j,i) = fac_sys  ! here regular rep is already unitary
          exit
        end if
      end do !k
    end do ! j
  end do ! i
  !
  ! small representation
  !
  ! generate random Hermite matrix and diagonalize
  do while(.true.)
    call random_number(rmat)
    call random_number(imat)
    rand_mat = rmat + im*imat
    rand_herm = rand_mat + conjg(transpose(rand_mat))
    !
    ! construct Hermite matrix commute with reg_rep
    comm_herm = (0.d0,0.d0)
    do i = 1, lit_ord
      comm_herm(:,:) = comm_herm(:,:)+matmul(lit_reg(:,:,i),matmul(rand_herm,conjg(transpose(lit_reg(:,:,i)))))
    end do
    !
    ! diagonalization with QR iteration
    evc = (0.d0,0.d0) ! eigenvectors
    do i = 1, lit_ord
      evc(i,i) = (1.d0,0.d0)
    end do
    n = 3000 ! loops
    do i = 1, n
      call QR(lit_ord,comm_herm,Q,R)
      comm_herm = matmul(R,Q)
      evc = matmul(evc,Q)
    end do
    ! convergence
    is_diag = .true.
    do i = 1, lit_ord
      do j = 1, lit_ord
        if(i /= j .and. abs(comm_herm(i,j)) > err_tol) is_diag = .false.
      end do
    end do
    if(is_diag) then
      write (stdout,*) 'Diagonalization succeeded'
      exit
    else
      write (stdout,*) 'Diagonalization failed'
      write (stdout,*) 'Now it is retrying with other random Hermite matrix...' 
      write (stdout,*) 'if still fail, try to increase the QR iterations'
    end if
    write (stdout,*)
  end do ! random matrix
  !
  ! eigenvalues
  do i = 1, lit_ord
    eva(i) = real(comm_herm(i,i))
  end do
  ! sort the eigenvalues in ascending order
  do i = 1, lit_ord-1
    do j = i+1, lit_ord
      if(eva(i) > eva(j)) then
        eva_tem = eva(i); eva(i) = eva(j); eva(j) = eva_tem
        evc_tem(:) = evc(:,i); evc(:,i) = evc(:,j); evc(:,j) = evc_tem(:)
      end if
    end do
  end do
  !
  ! get the irreps
  i = 1
  nirr = 1
  err_toln = err_tol
  lit_irr = (0.d0,0.d0)
  lit_cha = (0.d0,0.d0)
  do while(.true.)
    eva_tem = eva(i)
    do p = i, lit_ord
      if(abs(eva_tem-eva(p)) < err_tol) j = p
    end do
    !
    cha_tem = (0.d0,0.d0)
    do k = 1, lit_ord
      rep_tem(1:j-i+1,1:j-i+1,k) = matmul(conjg(transpose(evc(:,i:j))),matmul(lit_reg(:,:,k),evc(:,i:j)))
      do l = 1, j-i+1
        cha_tem(k) = cha_tem(k) + rep_tem(l,l,k)*exp(-im*tpi*dot_product(k_crys,lit_ele(:,4,k))) ! trace include translation
      end do
    end do
    !
    ! ignore the duplicated irreps
    do m = 1, nirr
      if(all(abs(cha_tem(:)-lit_cha(:,m)) < err_toln(:))) goto 100  ! abs(complex) = real
    end do
    lit_irr(1:j-i+1,1:j-i+1,:,nirr) = rep_tem(1:j-i+1,1:j-i+1,:)  ! irreducible projective representations of litle co-group
    lit_cha(:,nirr) = cha_tem(:)
    nirr = nirr + 1
100 continue
    !
    if(j == lit_ord) exit
    i = j + 1
  end do
  nirr = nirr - 1
  !
  ! sort the irreps
  do i = 1, nirr-1
    do j = i+1, nirr
      if(abs(lit_cha(1,i))-abs(lit_cha(1,j)) > err_tol .or. &
        (abs(lit_cha(1,i)-lit_cha(1,j)) < err_tol .and. abs(sum(lit_cha(:,i)))-abs(sum(lit_cha(:,j))) < -err_tol)) then
        cha_tem = lit_cha(:,i); rep_tem = lit_irr(:,:,:,i)
        lit_cha(:,i) = lit_cha(:,j); lit_irr(:,:,:,i) = lit_irr(:,:,:,j)
        lit_cha(:,j) = cha_tem; lit_irr(:,:,:,j) = rep_tem
      end if
    end do
  end do
  !
  return
end subroutine small_rep