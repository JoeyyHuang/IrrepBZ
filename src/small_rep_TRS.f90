subroutine small_rep_t(funit, k_crys, lit_ord, lit_ord_t, lit_spa_t, lit_ele_t, lit_symb_t, lit_axis_t, nirr_t, lit_cha_t, soc)
!=========================================================================
! get the regular projective representation of the given extended little group 
! and further construct the irreducible representation (small representation) by 
! numerically diagonalizing a random matrix commute with regular representation
!=========================================================================
  use constants, only : im, tpi, tol, tol3, tol22, tol33
  implicit none
!-------------------------------------------------------------------------
! arguments
  integer(kind=4),intent(in) :: funit, lit_ord, lit_ord_t
  integer(kind=4),intent(in) :: lit_spa_t(lit_ord_t), lit_axis_t(3,lit_ord_t)
  real(kind=8),intent(in) :: k_crys(3), lit_ele_t(3,4,lit_ord_t)
  character(len=2),intent(in) :: lit_symb_t(lit_ord_t)
  logical,intent(in) :: soc
  integer(kind=4),intent(out) :: nirr_t
  complex(kind=8),intent(out) :: lit_cha_t(lit_ord_t,lit_ord_t)   ! #element,#nirrep
! local variables
  integer(kind=4) :: i, j, k, l, m, n, p, ia, ib, ic, loop
  real(kind=8) :: eva_tem
  real(kind=8) :: k_crys2(3), fac_spin(lit_ord_t,lit_ord_t), g(3,lit_ord_t), eva(lit_ord_t), &
                  rmat(lit_ord_t,lit_ord_t), imat(lit_ord_t,lit_ord_t), toln(lit_ord_t)
  complex(kind=8) :: fac_sys
  complex(kind=8) :: rand_mat(lit_ord_t,lit_ord_t), rand_herm(lit_ord_t,lit_ord_t), comm_herm(lit_ord_t,lit_ord_t), &
                     cha_tem(lit_ord_t), rep_tem(lit_ord_t,lit_ord_t, lit_ord_t), &
                     lit_reg(lit_ord_t,lit_ord_t,lit_ord_t), evc(lit_ord_t,lit_ord_t), &
                     evc_tem(lit_ord_t), Q(lit_ord_t,lit_ord_t), R(lit_ord_t,lit_ord_t), &
                     ! notice that the general identity of \sum_{Dim}^2 = lit_ord_t no longer hold here
                     lit_irr_t(lit_ord_t,lit_ord_t,lit_ord_t,lit_ord_t)
  logical :: is_diag
!-------------------------------------------------------------------------
  write (funit,'(a)') 'Alternatively, take the anti-unitary time-reversal T as one group element,'
  write (funit,'(a)') 'so the combination operation TR, where reversal operation R satisfies Rk = -k,&
                    & is also one symmetry operation belonging to an extended little group.'
  write (funit,'(a)') 'Thus, the Irreps of original unitary little group, when time-reversal symmetry holds,&
                    & can be obtained by the small Reps of the extended little group.'
  if(lit_ord == lit_ord_t) then
    write (funit,'(a,i2,a$)') 'For current k-point, there is no reversal operation and&
                               & the extended little group is exactly the original unitary little group'
  else if(lit_ord_t-lit_ord == 1) then
    write (funit,'(a,i2,a$)') 'For current k-point, there is 1 reversal operation: '
    write (funit,'(a2,a,i2,a)') lit_symb_t(i), '[',lit_spa_t(i),']'
  else if(lit_ord_t-lit_ord > 1) then
    write (funit,'(a,i2,a$)') 'For current k-point, there are ', lit_ord_t-lit_ord, ' reversal operations: '
    do i = lit_ord+1, lit_ord_t-1
      write (funit,'(a2,a,i2,a,a,$)') lit_symb_t(i), '[',lit_spa_t(i),']',', '
    end do
    write (funit,'(a2,a,i2,a)') lit_symb_t(lit_ord_t), '[',lit_spa_t(lit_ord_t),']'
  end if
  write (funit,*)
  !
  if(soc) call spinor_t(lit_ord, lit_ord_t, lit_symb_t, lit_axis_t, lit_ele_t, fac_spin)
  !
  ! left projective regular representation
  lit_reg = (0.d0,0.d0)
  do i = 1, lit_ord_t
    k_crys2 = matmul(k_crys,lit_ele_t(:,1:3,i))
    do ia = -3, 3
      do ib = -3, 3
        do ic = -3, 3
          if(i <= lit_ord) then  ! k*alpha_i = k + G_i
            if(all(abs(k_crys2-k_crys-(/ia,ib,ic/)) < tol3)) g(:,i) = dble((/ia,ib,ic/))
          else  ! k*alpha_i = -k + G_i
            if(all(abs(k_crys2+k_crys-(/ia,ib,ic/)) < tol3)) g(:,i) = dble((/ia,ib,ic/))
          end if
        end do
      end do
    end do
    !
    do j = 1, lit_ord_t
      fac_sys = exp(-im*tpi*dot_product(g(:,i),lit_ele_t(:,4,j)))   ! factor system e^(-i*g_i*tau_j)
      if((i <= lit_ord .and. j <= lit_ord) .or. (i > lit_ord .and. j > lit_ord)) then
        do k = 1, lit_ord
          if(all(abs(matmul(lit_ele_t(:,1:3,i),lit_ele_t(:,1:3,j))-lit_ele_t(:,1:3,k)) < tol33)) then
            if(soc) lit_reg(k,j,i) = fac_sys*fac_spin(i,j)
            if(.NOT. soc) lit_reg(k,j,i) = fac_sys  ! here regular rep is already unitary
            exit
          end if
        end do !k
      else if((i <= lit_ord .and. j > lit_ord) .or. (i > lit_ord .and. j <= lit_ord)) then
        do k = lit_ord+1, lit_ord_t
          if(all(abs(matmul(lit_ele_t(:,1:3,i),lit_ele_t(:,1:3,j))-lit_ele_t(:,1:3,k)) < tol33)) then
            if(soc) lit_reg(k,j,i) = fac_sys*fac_spin(i,j)
            if(.NOT. soc) lit_reg(k,j,i) = fac_sys  ! here regular rep is already unitary
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
  loop = 0
  do while(loop < 20)
    loop = loop + 1
    call random_number(rmat)
    call random_number(imat)
    rand_mat = rmat + im*imat
    rand_herm = rand_mat + conjg(transpose(rand_mat))
    !
    ! construct Hermite matrix commute with reg_rep
    comm_herm = (0.d0,0.d0)
    do i = 1, lit_ord_t
      if(i <= lit_ord) comm_herm(:,:) = comm_herm(:,:)+matmul(lit_reg(:,:,i),matmul(rand_herm,conjg(transpose(lit_reg(:,:,i)))))
      if(i>lit_ord) comm_herm(:,:)=comm_herm(:,:)+matmul(lit_reg(:,:,i),matmul(conjg(rand_herm),conjg(transpose(lit_reg(:,:,i)))))
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
        if(i /= j .and. abs(comm_herm(i,j)) > tol) is_diag = .false.
      end do
    end do
    if(is_diag) then
      write (funit,'(a)') 'On-the-fly generate the small Reps of the extended little group by diagonalizing a random Hermite matrix'
      write (funit,'(a)') 'Diagonalization succeeded'
      write (funit,'(a)') 'Dim denotes the dimensionality of the Irrep'
      write (funit,*)
      exit
    else
      if(loop == 20) then
        write (funit,'(a)') 'Error: still failed after trying 20 times, stop the program'
        return
      end if
      write (funit,'(a)') 'On-the-fly generate the small Reps of the extended little group by diagonalizing a random Hermite matrix'
      write (funit,'(a)') 'Diagonalization failed'
      write (funit,'(a)') 'Now it is retrying with one other random Hermite matrix...'
    end if
  end do ! random matrix
  !
  ! eigenvalues
  do i = 1, lit_ord_t
    eva(i) = dreal(comm_herm(i,i))
  end do
  ! sort the eigenvalues in ascending order
  if(lit_ord_t > 1) then
    do i = 1, lit_ord_t-1
      do j = i+1, lit_ord_t
        if(eva(i) > eva(j)) then
          eva_tem = eva(i); eva(i) = eva(j); eva(j) = eva_tem
          evc_tem(:) = evc(:,i); evc(:,i) = evc(:,j); evc(:,j) = evc_tem(:)
        end if
      end do
    end do
  end if
  !
  ! get the irreps
  i = 1
  nirr_t = 1
  toln = tol
  lit_irr_t = (0.d0,0.d0)
  lit_cha_t = (0.d0,0.d0)
  do while(.true.)
    eva_tem = eva(i)
    do p = i, lit_ord_t
      if(abs(eva_tem-eva(p)) < tol) j = p
    end do
    !
    cha_tem = (0.d0,0.d0)
    do k = 1, lit_ord_t
      if(k <= lit_ord) rep_tem(1:j-i+1,1:j-i+1,k) = matmul(conjg(transpose(evc(:,i:j))),matmul(lit_reg(:,:,k),evc(:,i:j)))
      if(k > lit_ord) rep_tem(1:j-i+1,1:j-i+1,k) = matmul(conjg(transpose(evc(:,i:j))),matmul(lit_reg(:,:,k),evc(:,i:j)))
      do l = 1, j-i+1
        cha_tem(k) = cha_tem(k) + rep_tem(l,l,k)*exp(-im*tpi*dot_product(k_crys,lit_ele_t(:,4,k))) ! trace including translation
      end do
    end do
    !
    ! ignore the duplicated irreps
    do m = 1, nirr_t
      ! here we only compare the character of elements of the unitary little group
      if(all(abs(cha_tem(1:lit_ord)-lit_cha_t(1:lit_ord,m)) < toln(1:lit_ord))) goto 100  ! abs(complex) = real
     !if(all(abs(cha_tem(:)-lit_cha_t(:,m)) < toln(:))) goto 100  ! abs(complex) = real
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
  ! here only concern the character of elements of the unitary little group
  do i = 1, nirr_t-1
    do j = i+1, nirr_t
      if(abs(lit_cha_t(1,i))-abs(lit_cha_t(1,j)) > tol) then  ! different Dim
        cha_tem = lit_cha_t(:,i); rep_tem = lit_irr_t(:,:,:,i)
        lit_cha_t(:,i) = lit_cha_t(:,j); lit_irr_t(:,:,:,i) = lit_irr_t(:,:,:,j)
        lit_cha_t(:,j) = cha_tem; lit_irr_t(:,:,:,j) = rep_tem
      else if(abs(lit_cha_t(1,i)-lit_cha_t(1,j)) < tol) then  ! same Dim
        !
        do k = 1, lit_ord
        !do k = 1, lit_ord_t
          ! real Irrep precede complex Irrep
          if(abs(aimag(lit_cha_t(k,i))) < tol .and. abs(aimag(lit_cha_t(k,j))) > tol) then
            goto 200
          else if(abs(aimag(lit_cha_t(k,i))) > tol .and. abs(aimag(lit_cha_t(k,j))) < tol) then
            cha_tem = lit_cha_t(:,i); rep_tem = lit_irr_t(:,:,:,i)
            lit_cha_t(:,i) = lit_cha_t(:,j); lit_irr_t(:,:,:,i) = lit_irr_t(:,:,:,j)
            lit_cha_t(:,j) = cha_tem; lit_irr_t(:,:,:,j) = rep_tem
            goto 200
          ! positive complex Irrep precede negtive complex Irrep
          else if(abs(aimag(lit_cha_t(k,i))) > tol .and. abs(aimag(lit_cha_t(k,j))) > tol) then
            if(aimag(lit_cha_t(k,i))-aimag(lit_cha_t(k,j)) < -tol) then
              cha_tem = lit_cha_t(:,i); rep_tem = lit_irr_t(:,:,:,i)
              lit_cha_t(:,i) = lit_cha_t(:,j); lit_irr_t(:,:,:,i) = lit_irr_t(:,:,:,j)
              lit_cha_t(:,j) = cha_tem; lit_irr_t(:,:,:,j) = rep_tem
              goto 200
            end if
          end if
        end do
        do k = 1, lit_ord
        !do k = 1, lit_ord_t
          ! positive real Irrep precede negtive real Irrep
          if(dreal(lit_cha_t(k,i))-dreal(lit_cha_t(k,j)) > tol) then
            goto 200
          else if(dreal(lit_cha_t(k,i))-dreal(lit_cha_t(k,j)) < -tol) then
            cha_tem = lit_cha_t(:,i); rep_tem = lit_irr_t(:,:,:,i)
            lit_cha_t(:,i) = lit_cha_t(:,j); lit_irr_t(:,:,:,i) = lit_irr_t(:,:,:,j)
            lit_cha_t(:,j) = cha_tem; lit_irr_t(:,:,:,j) = rep_tem
            goto 200
          end if
        end do
        !
  200 end if
    end do
  end do
  !
  ! print small representation of little group
  ! here only print the character of elements of the unitary little group
  write (funit,'(a,i3)') 'Number of Irreps of the little group :', nirr_t
  do i = 1, lit_ord
    write (funit,'(a5,a2,a,i2,a,a5,$)') ' ', lit_symb_t(i), '[',lit_spa_t(i),']', ' '
  end do
  !write (funit,'(a,$)') ' | '
  !do i = lit_ord+1, lit_ord_t
  !  write (funit,'(a5,a2,a,i2,a,a5,$)') ' ', lit_symb_t(i), '[',lit_spa_t(i),']', ' '
  !end do
  write (funit,*)
  do i = 1, nirr_t
    write (funit,'(a,i2,a,5x,a,i2,5x)') '#',i,'*','Dim:', nint(dreal(lit_cha_t(1,i)))
    do j = 1, lit_ord
      write (funit,'("(",f6.3,","f6.3,") ",$)') lit_cha_t(j,i)+(1.d-8,1.d-8)
    end do
    !write (funit,'(a,$)') ' | '
    !do j = lit_ord+1, lit_ord_t
    !  write (funit,'("(",f6.3,","f6.3,") ",$)') lit_cha_t(j,i)+(1.d-8,1.d-8)
    !end do
    write (funit,*)
    write (funit,*)
  end do
  !
  return
end subroutine small_rep_t