subroutine small_rep_t(k_crys, lit_ord, lit_ord_t, lit_ele_t, lit_symb_t, lit_axis_t, nirr_t, lit_irr_t, lit_cha_t, spinor)
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
  integer(kind=4),intent(in) :: lit_ord, lit_ord_t
  integer(kind=4),intent(in) :: lit_axis_t(3,lit_ord_t)
  real(kind=8),intent(in) :: k_crys(3), lit_ele_t(3,4,lit_ord_t)
  character(len=2),intent(in) :: lit_symb_t(lit_ord_t)
  logical,intent(in) :: spinor
  integer(kind=4),intent(out) :: nirr_t
  complex(kind=8),intent(out) :: lit_cha_t(lit_ord_t,lit_ord_t)   ! #element,#nirrep
  complex(kind=8),intent(out) :: lit_irr_t(nint(sqrt(dble(lit_ord_t))),nint(sqrt(dble(lit_ord_t))),lit_ord_t,lit_ord_t)
! local variables
  integer(kind=4) :: i, j, k, l, m, n, p, ia, ib, ic
  real(kind=8) :: eva_tem
  real(kind=8) :: k_crys2(3), fac_spin(lit_ord_t,lit_ord_t), g(3,lit_ord_t), eva(lit_ord_t), &
                  rmat(lit_ord_t,lit_ord_t), imat(lit_ord_t,lit_ord_t), err_toln(lit_ord_t)
  complex(kind=8) :: fac_sys
  complex(kind=8) :: rand_mat(lit_ord_t,lit_ord_t), rand_herm(lit_ord_t,lit_ord_t), comm_herm(lit_ord_t,lit_ord_t), &
                     cha_tem(lit_ord_t), rep_tem(nint(sqrt(dble(lit_ord_t))),nint(sqrt(dble(lit_ord_t))), lit_ord_t), &
                     lit_reg(lit_ord_t,lit_ord_t,lit_ord_t), evc(lit_ord_t,lit_ord_t), &
                     evc_tem(lit_ord_t), Q(lit_ord_t,lit_ord_t), R(lit_ord_t,lit_ord_t)
  logical :: is_diag
!-------------------------------------------------------------------------
  if(spinor) call SU2(lit_ord_t, lit_symb_t, lit_axis_t, lit_ele_t, fac_spin)
  !
  ! left projective regular representation
  lit_reg = (0.d0,0.d0)
  do i = 1, lit_ord_t
    k_crys2 = matmul(k_crys,lit_ele_t(:,1:3,i))
    do ia = -3, 3
      do ib = -3, 3
        do ic = -3, 3
          if(i <= lit_ord) then  ! k*alpha_i = k + G_i
            if(all(abs(k_crys2-k_crys-(/ia,ib,ic/)) < err_tol3)) g(:,i) = dble((/ia,ib,ic/))
          else  ! k*alpha_i = -k + G_i
            if(all(abs(k_crys2+k_crys-(/ia,ib,ic/)) < err_tol3)) g(:,i) = -dble((/ia,ib,ic/))
          end if
        end do
      end do
    end do
    !
    do j = 1, lit_ord_t
      fac_sys = exp(-im*tpi*dot_product(g(:,i),lit_ele_t(:,4,j)))   ! factor system e^(-i*g_i*tau_j)
      if((i <= lit_ord .and. j <= lit_ord) .or. (i > lit_ord .and. j > lit_ord)) then
        do k = 1, lit_ord
          if(all(abs(matmul(lit_ele_t(:,1:3,i),lit_ele_t(:,1:3,j))-lit_ele_t(:,1:3,k)) < err_tol33)) then
            if(spinor) lit_reg(k,j,i) = fac_sys*fac_spin(i,j)
            if(.NOT. spinor) lit_reg(k,j,i) = fac_sys  ! here regular rep is already unitary
            exit
          end if
        end do !k
      else if((i <= lit_ord .and. j > lit_ord) .or. (i > lit_ord .and. j <= lit_ord)) then
        do k = lit_ord+1, lit_ord_t
          if(all(abs(matmul(lit_ele_t(:,1:3,i),lit_ele_t(:,1:3,j))-lit_ele_t(:,1:3,k)) < err_tol33)) then
            if(spinor) lit_reg(k,j,i) = fac_sys*fac_spin(i,j)
            if(.NOT. spinor) lit_reg(k,j,i) = fac_sys  ! here regular rep is already unitary
            exit
          end if
        end do !k
      end if
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
    do i = 1, lit_ord_t
      comm_herm(:,:) = comm_herm(:,:)+matmul(lit_reg(:,:,i),matmul(rand_herm,conjg(transpose(lit_reg(:,:,i)))))
    end do
    !
    ! diagonalization with QR iteration
    evc = (0.d0,0.d0) ! eigenvectors
    do i = 1, lit_ord_t
      evc(i,i) = (1.d0,0.d0)
    end do
    n = 3000 ! loops
    do i = 1, n
      call QR(lit_ord_t,comm_herm,Q,R)
      comm_herm = matmul(R,Q)
      evc = matmul(evc,Q)
    end do
    ! convergence
    is_diag = .true.
    do i = 1, lit_ord_t
      do j = 1, lit_ord_t
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
  do i = 1, lit_ord_t
    eva(i) = real(comm_herm(i,i))
  end do
  ! sort the eigenvalues in ascending order
  do i = 1, lit_ord_t-1
    do j = i+1, lit_ord_t
      if(eva(i) > eva(j)) then
        eva_tem = eva(i); eva(i) = eva(j); eva(j) = eva_tem
        evc_tem(:) = evc(:,i); evc(:,i) = evc(:,j); evc(:,j) = evc_tem(:)
      end if
    end do
  end do
  !
  ! get the irreps
  i = 1
  nirr_t = 1
  err_toln = err_tol
  lit_irr_t = (0.d0,0.d0)
  lit_cha_t = (0.d0,0.d0)
  do while(.true.)
    eva_tem = eva(i)
    do p = i, lit_ord_t
      if(abs(eva_tem-eva(p)) < err_tol) j = p
    end do
    !
    cha_tem = (0.d0,0.d0)
    do k = 1, lit_ord_t
      rep_tem(1:j-i+1,1:j-i+1,k) = matmul(conjg(transpose(evc(:,i:j))),matmul(lit_reg(:,:,k),evc(:,i:j)))
      do l = 1, j-i+1
        cha_tem(k) = cha_tem(k) + rep_tem(l,l,k)*exp(-im*tpi*dot_product(k_crys,lit_ele_t(:,4,k))) ! trace including translation
      end do
    end do
    !
    ! ignore the duplicated irreps
    do m = 1, nirr_t
      if(all(abs(cha_tem(:)-lit_cha_t(:,m)) < err_toln(:))) goto 100  ! abs(complex) = real
    end do
    lit_irr_t(1:j-i+1,1:j-i+1,:,nirr_t) = rep_tem(1:j-i+1,1:j-i+1,:)  ! irreducible projective representations of litle co-group
    lit_cha_t(:,nirr_t) = cha_tem(:)
    nirr_t = nirr_t + 1
100 continue
    !
    if(j == lit_ord_t) exit
    i = j + 1
  end do
  nirr_t = nirr_t - 1
  !
  ! sort the irreps
  do i = 1, nirr_t-1
    do j = i+1, nirr_t
      if(abs(lit_cha_t(1,i))-abs(lit_cha_t(1,j)) > err_tol .or. &
        (abs(lit_cha_t(1,i)-lit_cha_t(1,j)) < err_tol .and. abs(sum(lit_cha_t(:,i)))-abs(sum(lit_cha_t(:,j))) < -err_tol)) then
        cha_tem = lit_cha_t(:,i); rep_tem = lit_irr_t(:,:,:,i)
        lit_cha_t(:,i) = lit_cha_t(:,j); lit_irr_t(:,:,:,i) = lit_irr_t(:,:,:,j)
        lit_cha_t(:,j) = cha_tem; lit_irr_t(:,:,:,j) = rep_tem
      end if
    end do
  end do
  !      
  return
end subroutine small_rep_t