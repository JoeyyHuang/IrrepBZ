subroutine Herring(funit, k_crys, lit_ord, lit_spa, lit_ele, lit_symb, nirr, lit_cha, soc)
!=========================================================================
! Herring rules to determine the reality of irreps and further 
! give the additional degeneracy arising from time-reversal symmetry
!=========================================================================
  use variables, only : spa_ord, spa_ele, spa_symb, spa_axis
  use mat_trans, only : mat_inv
  use constants, only : im, tpi, tol, tol3, tol22, tol33
  implicit none
!-------------------------------------------------------------------------
! arguments
  integer(kind=4),intent(in) :: funit, lit_ord, nirr
  real(kind=8),intent(in) :: k_crys(3)
  real(kind=8),intent(in) :: lit_ele(3,4,lit_ord)
  integer(kind=4),intent(in) :: lit_spa(lit_ord)
  complex(kind=8),intent(in) :: lit_cha(lit_ord,lit_ord)   ! #operation,#nirrep
  character(len=2),intent(in) :: lit_symb(lit_ord)
  logical,intent(in) :: soc
! local variables
  integer(kind=4) :: i, j, k, p, q, q1, q2, q3, ia, ib, ic, nrev
  integer(kind=4) :: rev(spa_ord)
  real(kind=8) :: fac
  real(kind=8) :: k_crys2(3), square(3,4), toln(lit_ord), fac_spin(spa_ord,spa_ord)
  complex(kind=8) :: reality(nirr), vec_rev(3), cha_rev(lit_ord)
  character(len=2) :: which1, which2  ! if complex rep, which irrep will be in degeneracy with the irrep of study
  character(len=9):: time_deg(nirr)
!-------------------------------------------------------------------------
  if(soc) call spinor(spa_ord, spa_symb, spa_axis, spa_ele, fac_spin)
  !
  ! reversed opeartions that meet k*alpha_i^-1 = -k + G_i
  nrev = 0
  do i = 1, spa_ord
    k_crys2 = matmul(k_crys,mat_inv(spa_ele(:,1:3,i)))
    do ia = -3, 3
      do ib = -3, 3
        do ic = -3, 3
          if(all(abs(k_crys2+k_crys-(/ia,ib,ic/)) < tol3)) then
            nrev = nrev + 1
            rev(nrev) = i
          end if
        end do
      end do
    end do
  end do ! i
  !
  toln = tol
  !
  ! Herring rules to determine the reality of irreps
  time_deg = '??????'
  if(nrev == 0) then
    time_deg = 'none'
  else
    reality(:) = (0.d0,0.d0)
    do i = 1, nrev
      square(:,1:3) = matmul(spa_ele(:,1:3,rev(i)),spa_ele(:,1:3,rev(i)))
      square(:,4) = matmul(spa_ele(:,1:3,rev(i)),spa_ele(:,4,rev(i)))+spa_ele(:,4,rev(i))
      do j = 1, lit_ord
        if(all(abs(square(:,1:3)-lit_ele(:,1:3,j)) < tol33)) then
          do k = 1, nirr
            if(soc) reality(k) = reality(k) + fac_spin(rev(i),rev(i))*lit_cha(j,k)*&
                                 exp(-im*tpi*dot_product(k_crys,square(:,4)-lit_ele(:,4,j)))
            if(.NOT. soc) reality(k) = reality(k) + lit_cha(j,k)*&
                                       exp(-im*tpi*dot_product(k_crys,square(:,4)-lit_ele(:,4,j)))
          end do ! k
          exit
        end if
      end do ! j
    end do ! i
    !
    do i = 1, nirr
      write (which1,'(i2)') i
      if(abs(imag(reality(i))) < tol .and. abs(real(reality(i))-lit_ord) < tol) then
        if(soc) time_deg(i) = 'double'
        if(.NOT. soc) time_deg(i) = 'none'
      else if(abs(imag(reality(i))) < tol .and. abs(real(reality(i))+lit_ord) < tol) then
        if(soc) time_deg(i) = 'none'
        if(.NOT. soc) time_deg(i) = 'double'
      else if(abs(imag(reality(i))) < tol .and. abs(real(reality(i))) < tol) then
        do j = 1, lit_ord
          vec_rev(:) = matmul(mat_inv(spa_ele(:,1:3,rev(1))),&
                       matmul(lit_ele(:,1:3,j),spa_ele(:,4,rev(1)))+lit_ele(:,4,j)-spa_ele(:,4,rev(1)))
          !
          if(soc) then
            do q = 1, spa_ord
              if(all(abs(matmul(spa_ele(:,1:3,rev(1)),spa_ele(:,1:3,q))-spa_ele(:,1:3,1)) < tol33)) then
                q1 = q
                fac = fac_spin(rev(1),q1)
                exit
              end if
            end do
            do q = 1, spa_ord
              if(all(abs(lit_ele(:,1:3,j)-spa_ele(:,1:3,q)) < tol33)) then
                q2 = q
                fac = fac*fac_spin(q1,q2)
                exit
              end if
            end do
            do q = 1, spa_ord
              if(all(abs(matmul(spa_ele(:,1:3,q1),lit_ele(:,1:3,j))-spa_ele(:,1:3,q)) < tol33)) then
                q3 = q
                fac = fac*fac_spin(q3,rev(1))
                exit
              end if
            end do
          end if ! spinor
          !
          do p = 1, lit_ord
            if(all(abs(matmul(lit_ele(:,1:3,j),spa_ele(:,1:3,rev(1))) - &
                       matmul(spa_ele(:,1:3,rev(1)),lit_ele(:,1:3,p))) < tol33)) then
              if(soc) cha_rev(j) = conjg(fac*lit_cha(p,i)*exp(-im*tpi*dot_product(k_crys,vec_rev(:)-lit_ele(:,4,p))))
              if(.NOT. soc) cha_rev(j) = conjg(lit_cha(p,i)*exp(-im*tpi*dot_product(k_crys,vec_rev(:)-lit_ele(:,4,p))))
              exit
            end if
          end do
        end do
        do k = 1, nirr
          if(all(abs(cha_rev-lit_cha(:,k)) < toln)) then
            write (which2,'(i2)') k
            time_deg(i) = '#'//trim(adjustl(which1))//' + '//'#'//trim(adjustl(which2))
            exit
          end if
        end do
      end if ! reality
      !
    end do ! i
  end if ! nrev
  !
  ! print small representation of little group
  write (funit,'(a,i3)') 'Number of Irreps of the little group:', nirr
  do i = 1, lit_ord
    write (funit,'(a5,a2,a,i2,a,a5,$)') ' ',lit_symb(i),'[',lit_spa(i),']',' '
  end do
  write (funit,*)
  do i = 1, nirr
    write (funit,'(a,i2,5x,a,i2,5x,a,a)') '#',i,'Dim:', nint(dreal(lit_cha(1,i))), 'TRS: ',time_deg(i)
    do j = 1, lit_ord
      write (funit,'("(",f6.3,","f6.3,") ",$)') lit_cha(j,i)+(1.d-8,1.d-8)
    end do
    write (funit,*)
    write (funit,*)
  end do
  !
  return
end subroutine Herring