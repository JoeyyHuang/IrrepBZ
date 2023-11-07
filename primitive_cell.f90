subroutine primitive_cell()
!=========================================================================
! find from the origical cell three translation vectors
! that can serve as the primitive basis vectors
!=========================================================================
  use variables, only : stdout, ntype, natom_conv_type, natom_prim_type, a_prim, &
                        a_conv, natom_conv, natom_prim, atom_conv_crys, atom_prim_crys, &
                        ntype, omega_conv, omega_prim, atom_conv_cart, atom_prim_cart
  use mat_trans, only : mat_inv
  use constants, only : err_tol, err_tol3
  implicit none
!-------------------------------------------------------------------------
! local varibales
  integer(kind=4) :: i, j, k, m, n, p, itype, ia, ib, ic, ntrans, start, omega_ratio
  integer(kind=4) :: loc(1)
  real(kind=8) :: vec(3), diff(3), basis(3), trans_list(3,natom_conv_type(1)+2), atom_tem(3,natom_conv)
  real(kind=8),allocatable :: a_prim_list(:,:,:), a_prim_len(:)
  logical :: is_trans(natom_conv)
  logical,allocatable :: is_basis(:)
!-------------------------------------------------------------------------
  !
  allocate(natom_prim_type(ntype))
  ! if there is only one atom in type1, original cell is already primitive
  if(natom_conv_type(1) == 1) then
    a_prim = a_conv
    omega_prim = omega_conv
    natom_prim = natom_conv
    natom_prim_type = natom_conv_type
    allocate(atom_prim_crys(3,natom_prim),atom_prim_cart(3,natom_prim))
    atom_prim_crys = atom_conv_crys
    atom_prim_cart = atom_conv_cart
    write (stdout,*) 'original cell is already primitive'
    write (stdout,*)
    write (stdout,'(a,i3)') 'atoms number in primitive cell is:', natom_prim
    write (stdout,*)
    write (stdout,*) 'primitive basis vectors:'
    write (stdout,'(3f15.10)') a_prim(:,1)
    write (stdout,'(3f15.10)') a_prim(:,2)
    write (stdout,'(3f15.10)') a_prim(:,3)
    write (stdout,*)
    do i = 1, natom_prim
      write (stdout,'(3f15.10)') atom_prim_crys(:,i)
    end do
    return
  end if
  !
  ntrans = 0
  trans_list = 0.d0   ! (coords,counter)
  !
  ! vector bwtween one atom and other rest atoms in the same type, here type 1
  do i = 2, natom_conv_type(1)
    vec(:) = atom_conv_cart(:,i) - atom_conv_cart(:,1)
    !
    is_trans(:) = .false.
    ! check if this vector is a symmetry-allowable translation vector
    do itype = 1, ntype
      if(itype == 1) then
        start = 1
      else
        start = sum(natom_conv_type(1:(itype-1)))+1
      end if
      do j = start, sum(natom_conv_type(1:itype))
        do k = start, sum(natom_conv_type(1:itype))
          do ia = -3, 3
            do ib = -3, 3
              do ic = -3, 3
                diff(:) = atom_conv_cart(:,j) + vec(:) - atom_conv_cart(:,k) - &
                          ia*a_conv(:,1) - ib*a_conv(:,2) - ic*a_conv(:,3)
                if(all(abs(diff) < err_tol3)) is_trans(j) = .true.
              end do
            end do
          end do
        end do
      end do
    end do
    !
    ! if all true, then vec is one symmetry-allowable translation vector
    if(all(is_trans)) then
      ntrans = ntrans + 1
      trans_list(:,ntrans) = vec(:)
    end if
    !
  end do ! i
  !
  ! if there is no translation vector, original cell is already primitive
  if(ntrans == 0) then
    a_prim = a_conv
    omega_prim = omega_conv
    natom_prim = natom_conv
    natom_prim_type = natom_conv_type
    allocate(atom_prim_crys(3,natom_prim),atom_prim_cart(3,natom_prim))
    atom_prim_crys = atom_conv_crys
    atom_prim_cart = atom_conv_cart
    write (stdout,*) 'original cell is already primitive'
    write (stdout,*)
    write (stdout,'(a,i3)') 'atoms number in primitive cell is:', natom_prim
    write (stdout,*)
    write (stdout,*) 'primitive basis vectors:'
    write (stdout,'(3f15.10)') a_prim(:,1)
    write (stdout,'(3f15.10)') a_prim(:,2)
    write (stdout,'(3f15.10)') a_prim(:,3)
    write (stdout,*)
    do i = 1, natom_prim
      write (stdout,'(3f15.10)') atom_prim_crys(:,i)
    end do
    return
  else
    ! From the translation list, find three as primitive basis vectors
    trans_list(:,ntrans+1) = a_conv(:,1)
    trans_list(:,ntrans+2) = a_conv(:,2)
    trans_list(:,ntrans+3) = a_conv(:,3)
    !
    allocate(is_basis(ntrans+3))
    allocate(a_prim_list(3,3,((ntrans+3)*(ntrans+2)*(ntrans+1)/6)))
    allocate(a_prim_len((ntrans+3)*(ntrans+2)*(ntrans+1)/6))
    m = 0
    !
    do i = 1, ntrans+1
      do j = i+1, ntrans+2
        do k = j+1, ntrans+3
          !
          is_basis(:) = .false.
          !
          do p = 1, ntrans+3
            do ia = -3, 3
              do ib = -3, 3
                do ic = -3, 3
                  basis(:) = trans_list(:,p)-ia*trans_list(:,i)-ib*trans_list(:,j)-ic*trans_list(:,k)
                  if(all(abs(basis) < err_tol3)) is_basis(p) = .true.
                end do
              end do
            end do
          end do
          if(all(is_basis)) then
            m = m + 1
            a_prim_list(:,1,m) = trans_list(:,i)
            a_prim_list(:,2,m) = trans_list(:,j)
            a_prim_list(:,3,m) = trans_list(:,k)
            a_prim_len(m) = sqrt(dot_product(a_prim_list(:,1,m),a_prim_list(:,1,m)))*&
                            sqrt(dot_product(a_prim_list(:,2,m),a_prim_list(:,2,m)))*&
                            sqrt(dot_product(a_prim_list(:,3,m),a_prim_list(:,3,m)))
          end if
          !
        end do ! k
      end do ! j
    end do ! i
    !
    loc = minloc(a_prim_len(1:m))
    a_prim(:,1) = a_prim_list(:,1,loc(1))
    a_prim(:,2) = a_prim_list(:,2,loc(1))
    a_prim(:,3) = a_prim_list(:,3,loc(1))
    omega_prim = (a_prim(2,1)*a_prim(3,2)-a_prim(2,2)*a_prim(3,1))*a_prim(1,3) - &
                 (a_prim(1,1)*a_prim(3,2)-a_prim(1,2)*a_prim(3,1))*a_prim(2,3) + &
                 (a_prim(1,1)*a_prim(2,2)-a_prim(1,2)*a_prim(2,1))*a_prim(3,3)
    if(omega_prim < 0.d0) then
      a_prim(:,1) = a_prim_list(:,2,loc(1))
      a_prim(:,2) = a_prim_list(:,1,loc(1))
      a_prim(:,3) = a_prim_list(:,3,loc(1))
      omega_prim = -omega_prim
    end if
    deallocate(is_basis,a_prim_list,a_prim_len)
    !
    omega_ratio = nint(omega_conv/omega_prim)
    natom_prim = natom_conv/omega_ratio
    natom_prim_type = natom_conv_type/omega_ratio
    allocate(atom_prim_crys(3,natom_prim),atom_prim_cart(3,natom_prim))
    !
    ! get the atomic positions in crystal coords in primitive cell
    n = 0
    do itype = 1, ntype
      if(itype == 1) then
        start = 1
      else
        start = sum(natom_conv_type(1:(itype-1)))+1
      end if
      do i = start, sum(natom_conv_type(1:itype))
        atom_tem(:,i) = matmul(mat_inv(a_prim),matmul(a_conv,atom_conv_crys(:,i))) ! crystal coords
      end do
      !
      do j = start, sum(natom_conv_type(1:itype))-1
        do k = j+1, sum(natom_conv_type(1:itype))
          do ia = -3, 3
            do ib = -3, 3
              do ic = -3, 3
                diff = atom_tem(:,j)-atom_tem(:,k)-dble((/ia,ib,ic/))
                if(all(abs(diff) < err_tol3)) goto 100
              end do
            end do
          end do
        end do ! k
        n = n + 1
        atom_prim_crys(:,n) = atom_tem(:,j)
        100  continue
      end do ! j
      n = n + 1
      atom_prim_crys(:,n) = atom_tem(:,j)
    end do ! itype
    !
    if(n /= natom_prim) then
      write (stdout,*) 'wrong number of atoms in primitie cell'
      stop
    end if
    !
    ! convert the atomic coords to 0~1
    do i = 1, n
      do j = 1, 3
        if(atom_prim_crys(j,i) < -err_tol .or. atom_prim_crys(j,i) > 1.d0-err_tol) then
          do ia = -3, 3
            if(.NOT. (atom_prim_crys(j,i)+dble(ia) < -err_tol .or. atom_prim_crys(j,i)+dble(ia) > 1.d0-err_tol)) &
            atom_prim_crys(j,i) = atom_prim_crys(j,i)+dble(ia)
          end do
        end if
      end do
    end do
    !
    atom_prim_cart = matmul(a_prim,atom_prim_crys)
    !
    write (stdout,'(a,i3)') 'atoms number in primitive cell is:', natom_prim
    write (stdout,*)
    write (stdout,*) 'primitive basis vectors:'
    write (stdout,'(3f15.10)') a_prim(:,1)
    write (stdout,'(3f15.10)') a_prim(:,2)
    write (stdout,'(3f15.10)') a_prim(:,3)
    write (stdout,*)
    do i = 1, natom_prim
      write (stdout,'(3f15.10)') atom_prim_crys(:,i)
    end do
    return
  end if ! ntrans
  !
end subroutine primitive_cell
  
   