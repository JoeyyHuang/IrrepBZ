subroutine qe_bands()
!=========================================================================
! group analysis of QE electronic bands
!=========================================================================
  use variables, only : a_prim
  use constants, only : tpi, im, tol, tol3, qe_band
  use mat_trans, only : mat_inv
  use read_binary, only : read_a_wfc
  implicit none
!-------------------------------------------------------------------------
! local variables
  integer(kind=4) :: i, j, k, l, p, ia, ib, ic, ik, nk, num, tot, ibnd, nbnd, lit_ord, nirr, nirr_t, nirr_t2, lit_ord_t
  integer(kind=4),allocatable :: lit_spa(:), lit_spa_t(:), lit_axis(:,:), lit_axis_t(:,:), npw(:), gvector(:,:)
  real(kind=8) :: angle
  real(kind=8) :: k_crys(3), k_crys2(3), gvector2(3), vec(3)
  real(kind=8),allocatable :: eva(:,:), kpoints(:,:), lit_ele(:,:,:), lit_ele_t(:,:,:), toln(:)
  complex(kind=8) :: spinor1(2), spinor2(2), SU2_mat(2,2)
  complex(kind=8),allocatable :: lit_cha(:,:), lit_cha_t(:,:), lit_cha_t2(:,:), Cg(:,:), band_cha(:,:), cha(:)
  character(len=132) :: line
  character(len=10) klabel
  character(len=12) wfc
  character(len=14) fname
  character(len=2),allocatable :: lit_symb(:), lit_symb_t(:)
  logical :: soc
!-------------------------------------------------------------------------
  ! read k-points and eigenvalues from data-file-schema.xml
  ik = 0
  open(unit=10, file='qe_band/data-file-schema.xml', status='old')
  do while(.true.)
    read(10,'(a132)') line
    p = index(line,'<nbnd>',back=.false.,kind=4)
    if(p /= 0) then
      backspace(10)
      backspace(10)
      read(10,'(a132)') line
      i = index(line,'>',back=.false.,kind=4)
      j = index(line,'</',back=.false.,kind=4)
      if(line(i+1:j-1) == 'false') soc = .false.
      if(line(i+1:j-1) == 'true') soc = .true.
      read(10,'(a132)') line
      i = index(line,'>',back=.false.,kind=4)
      j = index(line,'</',back=.false.,kind=4)
      read(line(i+1:j-1),*) nbnd
    end if
    p = index(line,'<starting_k_points>',back=.false.,kind=4)
    if(p /= 0) then
      read(10,'(a132)') line
      i = index(line,'>',back=.false.,kind=4)
      j = index(line,'</',back=.false.,kind=4)
      read(line(i+1:j-1),*) nk
      allocate(kpoints(3,nk), npw(nk), eva(nbnd,nk))
      do k = 1, nk
        read(10,'(a132)') line
        i = index(line,'>',back=.false.,kind=4)
        j = index(line,'</',back=.false.,kind=4)
        read(line(i+1:j-1),*) kpoints(:,k)
      end do
    end if
    p = index(line,'<npw>',back=.false.,kind=4)
    if(p /= 0) then
      ik = ik + 1
      i = index(line,'>',back=.false.,kind=4)
      j = index(line,'</',back=.false.,kind=4)
      read(line(i+1:j-1),*) npw(ik)
      read(10,'(a132)') line
      tot = 1
      do while(.true.)
        p = 1
        num = 0
        read(10,'(a132)') line
        do while(.true.)
          i = verify(line(p:), ' ')      !-- Find next non-blank
          if (i == 0) exit               !-- No word found
          num = num + 1                  !-- Found something
          p = p + i - 1                  !-- Move to start of the word
          i = scan(line(p:), ' ')        !-- Find next blank
          if (i == 0) exit               !-- No blank found
          p = p + i - 1                  !-- Move to the blank
        end do
        read(line,*) eva(tot:tot+num-1,ik)
        if(tot+num-1 == nbnd) exit
        tot = tot + num
      end do
      if(ik == nk) exit
    end if
    !
  end do
  close(10)
  !
  ! convert to cryst. coord.
  kpoints = matmul(transpose(a_prim),kpoints) / sqrt(a_prim(1,1)**2+a_prim(2,1)**2+a_prim(3,1)**2)
  !
  ! loops on k-points
  do ik = 1, nk
    if(ik < 10) then
      write(wfc,'(a3,i1.1,a4)') 'wfc', ik, '.dat'
      write(fname,'(a7,i1.1,a4)') 'band_ik', ik, '.txt'
    else if(ik < 100) then
      write(wfc,'(a3,i2.2,a4)') 'wfc', ik, '.dat'
      write(fname,'(a7,i2.2,a4)') 'band_ik', ik, '.txt'
    else if(ik < 1000) then
      write(wfc,'(a3,i3.3,a4)') 'wfc', ik, '.dat'
      write(fname,'(a7,i3.3,a4)') 'band_ik', ik, '.txt'
    end if
    open(unit=qe_band, file='qe_band/'//trim(fname), status='replace')
    if(.NOT. soc) write (qe_band,'(a)') 'The electronic bands are calculated without considering SOC'
    if(soc) write (qe_band,'(a)') 'The electronic bands are calculated with considering SOC'
    write (qe_band,'(a,i3)') 'Number of k points:', nk
    write (qe_band,'(a,i3/)') 'Number of electronic bands at each k point:', nbnd
    !
    ! read wfc#.dat file
    if(.NOT. soc) allocate(gvector(3,npw(ik)), Cg(npw(ik),nbnd))
    if(soc) allocate(gvector(3,npw(ik)), Cg(2*npw(ik),nbnd))
    call read_a_wfc('qe_band/'//wfc, gvector, Cg)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Irreps of little group
    k_crys = kpoints(:,ik)
    if(ik == 1) then
      write(klabel,'(a)') '1st'
    else if(ik == 2) then
      write(klabel,'(a)') '2nd'
    else if(ik == 3) then
      write(klabel,'(a)') '3rd'
    else if(ik < 10) then
      write(klabel,'(i1,a)') ik, 'th'
    else if(ik < 100) then
      write(klabel,'(i2,a)') ik, 'th'
    else if(ik < 1000) then
      write(klabel,'(i3,a)') ik, 'th'
    end if
    write (qe_band,'(a,a,a,3f10.5)') 'This is the ', trim(adjustl(klabel)), ' k-point (cryst. coord.):', k_crys
    !!!!!!!!!!!!!!!!!!!!! vector representation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ordinary little group with all unitary operations
    call little_order(k_crys, lit_ord)
    allocate(lit_spa(lit_ord),lit_ele(3,4,lit_ord), lit_symb(lit_ord), lit_axis(3,lit_ord), toln(lit_ord))
    call little_group(qe_band, k_crys, lit_ord, lit_spa, lit_ele, lit_symb, lit_axis)
    write (qe_band, '(a/)') REPEAT('=', 100)
    allocate(lit_cha(lit_ord,lit_ord))
    call small_rep(qe_band,k_crys,lit_ord,lit_ele,lit_symb,lit_axis,nirr,lit_cha,.false.)
    call Herring(qe_band, k_crys, lit_ord, lit_spa, lit_ele, lit_symb, nirr, lit_cha, .false.)
    ! extended little group when consider the anti-unitary operation of TRS as en element
    call little_order_t(k_crys, lit_ord, lit_ord_t)
    allocate(lit_spa_t(lit_ord_t),lit_ele_t(3,4,lit_ord_t),lit_symb_t(lit_ord_t),lit_axis_t(3,lit_ord_t))
    call little_group_t(qe_band, k_crys, lit_ord, lit_ord_t, lit_spa_t, lit_ele_t, lit_symb_t, lit_axis_t)
    allocate(lit_cha_t(lit_ord_t,lit_ord_t),lit_cha_t2(lit_ord_t,lit_ord_t))
    if(.NOT. soc) then
      call small_rep_t(qe_band,k_crys,lit_ord,lit_ord_t,lit_spa_t,lit_ele_t,lit_symb_t,lit_axis_t,nirr_t,lit_cha_t,.false.)
    else
      call small_rep_t(qe_band,k_crys,lit_ord,lit_ord_t,lit_spa_t,lit_ele_t,lit_symb_t,lit_axis_t,nirr_t2,lit_cha_t2,.false.)
    end if
    write (qe_band, '(a/)') REPEAT('=', 100)
    !!!!!!!!!!!!!!!!!!!!!!! spinor representation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ordinary little group with all unitary operations
    call small_rep(qe_band,k_crys,lit_ord,lit_ele,lit_symb,lit_axis,nirr,lit_cha,.true.)
    call Herring(qe_band, k_crys, lit_ord, lit_spa, lit_ele, lit_symb, nirr, lit_cha, .true.)
    ! extended little group when consider the anti-unitary operation of TRS as en element
    if(.NOT. soc) then
      call small_rep_t(qe_band,k_crys,lit_ord,lit_ord_t,lit_spa_t,lit_ele_t,lit_symb_t,lit_axis_t,nirr_t2,lit_cha_t2,.true.)
    else
      call small_rep_t(qe_band,k_crys,lit_ord,lit_ord_t,lit_spa_t,lit_ele_t,lit_symb_t,lit_axis_t,nirr_t,lit_cha_t,.true.)
    end if
    write (qe_band, '(a/)') REPEAT('=', 100)
    !
    ! calculate Irreps of QE electronic bands
    allocate(band_cha(lit_ord,nbnd), cha(lit_ord))
    band_cha = (0.d0,0.d0)
    do i = 1, lit_ord
      if(soc) then
        if(lit_symb(i) == 'E ' .or. lit_symb(i) == 'I ') then
          SU2_mat(:,:) = reshape((/(1.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0), (1.d0,0.d0)/),(/2, 2/))
        else
          if(lit_symb(i) == 'C2' .or. lit_symb(i) == 'm ') angle = tpi/2.d0
          if(lit_symb(i) == 'C3' .or. lit_symb(i) == 'S6') angle = tpi/3.d0
          if(lit_symb(i) == 'C4' .or. lit_symb(i) == 'S4') angle = tpi/4.d0
          if(lit_symb(i) == 'C6' .or. lit_symb(i) == 'S3') angle = tpi/6.d0
          vec = matmul(a_prim,dble(lit_axis(:,i)))
          vec = vec / sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
          SU2_mat(:,:) = reshape((/cos(angle/2.d0)-im*vec(3)*sin(angle/2.d0), (vec(2)-im*vec(1))*sin(angle/2.d0), &
                      -(vec(2)+im*vec(1))*sin(angle/2.d0), cos(angle/2.d0)+im*vec(3)*sin(angle/2.d0)/),(/2, 2/))
        end if
      end if
      !
      k_crys2 = matmul(k_crys,mat_inv(lit_ele(:,1:3,i)))
      do ia = -3, 3
        do ib = -3, 3
          do ic = -3, 3
            if(all(abs(k_crys2-k_crys-(/ia,ib,ic/)) < tol3)) goto 100
          end do
        end do
      end do
      100 continue
      do j = 1, npw(ik)
        gvector2 = matmul(gvector(:,j),mat_inv(lit_ele(:,1:3,i)))
        do k = 1, npw(ik)
          if(all(abs((/ia,ib,ic/)+gvector2-gvector(:,k)) < tol3)) then
            do l = 1, nbnd
              if(.NOT. soc) then
                band_cha(i,l) = band_cha(i,l) + &
                conjg(Cg(k,l))*Cg(j,l)*exp(im*tpi*dot_product(k_crys+(/ia,ib,ic/)+gvector2,-lit_ele(:,4,i)))
              else ! with soc
                spinor1(1)=Cg(k,l); spinor1(2)=Cg(npw(ik)+k,l)
                spinor2(1)=Cg(j,l); spinor2(2)=Cg(npw(ik)+j,l)
                band_cha(i,l) = band_cha(i,l) + dot_product(spinor1, matmul(SU2_mat,spinor2))*&
                exp(im*tpi*dot_product(k_crys+(/ia,ib,ic/)+gvector2,-lit_ele(:,4,i)))
              end if
            end do ! l
            exit
          end if
        end do ! k
      end do ! j
    end do ! i
    ! print Irreps of QE electronic bands
    write (qe_band,'(a)') '   <<  Irreps of Electronic Bands  >>'
    if(.NOT. soc) write (qe_band,'(a)') 'Determine the Irrep by comparing the band characters &
                                        &with the small Reps of the little group in "Scenario Without Considering SOC"'
    if(soc) write (qe_band,'(a)') 'Determine the Irrep by comparing the band characters &
                                  &with the small Reps of the little group in "Scenario With Considering SOC"'
    write (qe_band,'(a/)') 'For electronic Hamiltonian equation, the time-reversal symmetry is intrinsically holded'
    write (qe_band,'(a)') 'Characters of electronic bands associated with each symmetry opeartion'
    do i = 1, lit_ord
      write (qe_band,'(a5,a2,a,i2,a,a5,$)') ' ', lit_symb_t(i), '[',lit_spa_t(i),']', ' '
    end do
    write (qe_band,*)
    p = 1; k = 0; toln = tol
    do while(.true.)
      cha = (0.d0,0.d0)
      do i = p, nbnd
        if(abs(eva(i,ik)-eva(p,ik)) < tol) then
          cha = cha + band_cha(:,i)
          if(i == nbnd) k = nbnd
        else
          exit
        end if
      end do
      ! compare the bands characters with Irreps of little group
      do l = 1, nirr_t
        if(all(abs(lit_cha_t(1:lit_ord,l)-cha(1:lit_ord)) < toln(1:lit_ord))) exit
      end do
      !
      ! TODO: accidental degeneracy
      if(l > nirr_t) then
        if(i-p == 2) then
          do ibnd = p, p+1
            do l  = 1, nirr_t
              if(all(abs(lit_cha_t(1:lit_ord,l)-cha(1:lit_ord)) < toln(1:lit_ord))) exit
            end do
            write(qe_band,'(a,$)') 'bands:'
            write(qe_band,'(i4,$)') ibnd
            write (qe_band,'(a5,a,a,i2,a,a5,a,2i4)') ' ', 'Irrep: ' , '#', l, '*', ' ', 'accidental degeneracy:', p, p+1
            do j = 1, lit_ord
              write(qe_band,'("(",f6.3,","f6.3,") ",$)') band_cha(j,ibnd)+(1.d-8,1.d-8)
            end do
            write (qe_band, *)
          end do ! ibnd
        end if ! 2
        if(i-p == 3) then
          do l  = 1, nirr_t
            if(all(abs(lit_cha_t(1:lit_ord,l)-band_cha(1:lit_ord,p)) < toln(1:lit_ord))) exit
          end do
          if(l <= nirr_t) then
            write(qe_band,'(a,$)') 'bands:'
            write(qe_band,'(i4,$)') p
            write (qe_band,'(a5,a,a,i2,a,a5,a,3i4)') ' ', 'Irrep: ' , '#', l, '*', ' ', 'accidental degeneracy:', p, p+1, p+2
            do j = 1, lit_ord
              write(qe_band,'("(",f6.3,","f6.3,") ",$)') band_cha(j,p)+(1.d-8,1.d-8)
            end do
            write (qe_band, *)
            do l  = 1, nirr_t
              if(all(abs(lit_cha_t(1:lit_ord,l)-band_cha(1:lit_ord,p+1)) < toln(1:lit_ord))) exit
            end do
            if(l <= nirr_t) then
              write(qe_band,'(a,$)') 'bands:'
              write(qe_band,'(i4,$)') p+1
              write (qe_band,'(a5,a,a,i2,a,a5,a,3i4)') ' ', 'Irrep: ' , '#', l, '*', ' ', 'accidental degeneracy:', p, p+1, p+2
              do j = 1, lit_ord
                write(qe_band,'("(",f6.3,","f6.3,") ",$)') band_cha(j,p+1)+(1.d-8,1.d-8)
              end do
              write (qe_band, *)
              do l  = 1, nirr_t
                if(all(abs(lit_cha_t(1:lit_ord,l)-band_cha(1:lit_ord,p+2)) < toln(1:lit_ord))) exit
              end do
              write(qe_band,'(a,$)') 'bands:'
              write(qe_band,'(i4,$)') p+2
              write (qe_band,'(a5,a,a,i2,a,a5,a,3i4)') ' ', 'Irrep: ' , '#', l, '*', ' ', 'accidental degeneracy:', p, p+1, p+2
              do j = 1, lit_ord
                write(qe_band,'("(",f6.3,","f6.3,") ",$)') band_cha(j,p+2)+(1.d-8,1.d-8)
              end do
              write (qe_band, *)
            else
              do l  = 1, nirr_t
                if(all(abs(lit_cha_t(1:lit_ord,l)-band_cha(1:lit_ord,p+1)-band_cha(1:lit_ord,p+2)) < toln(1:lit_ord))) exit
              end do
              write(qe_band,'(a,$)') 'bands:'
              write(qe_band,'(2i4,$)') p+1, p+2
              write (qe_band,'(a5,a,a,i2,a,a5,a,3i4)') ' ', 'Irrep: ' , '#', l, '*', ' ', 'accidental degeneracy:', p, p+1, p+2
              do j = 1, lit_ord
                write(qe_band,'("(",f6.3,","f6.3,") ",$)') band_cha(j,p+1)+band_cha(j,p+2)+(1.d-8,1.d-8)
              end do
              write (qe_band, *)
            end if
          else
            do l  = 1, nirr_t
              if(all(abs(lit_cha_t(1:lit_ord,l)-band_cha(1:lit_ord,p)-band_cha(1:lit_ord,p+1)) < toln(1:lit_ord))) exit
            end do
            if(l <= nirr_t) then
              write(qe_band,'(a,$)') 'bands:'
              write(qe_band,'(2i4,$)') p, p+1
              write (qe_band,'(a5,a,a,i2,a,a5,a,3i4)') ' ', 'Irrep: ' , '#', l, '*', ' ', 'accidental degeneracy:', p, p+1, p+2
              do j = 1, lit_ord
                write(qe_band,'("(",f6.3,","f6.3,") ",$)') band_cha(j,p)+band_cha(j,p+1)+(1.d-8,1.d-8)
              end do
              write (qe_band, *)
              do l  = 1, nirr_t
                if(all(abs(lit_cha_t(1:lit_ord,l)-band_cha(1:lit_ord,p+2)) < toln(1:lit_ord))) exit
              end do
              write(qe_band,'(a,$)') 'bands:'
              write(qe_band,'(i4,$)') p+2
              write (qe_band,'(a5,a,a,i2,a,a5,a,3i4)') ' ', 'Irrep: ' , '#', l, '*', ' ', 'accidental degeneracy:', p, p+1, p+2
              do j = 1, lit_ord
                write(qe_band,'("(",f6.3,","f6.3,") ",$)') band_cha(j,p+2)+(1.d-8,1.d-8)
              end do
              write (qe_band, *)
            else
              do l  = 1, nirr_t
                if(all(abs(lit_cha_t(1:lit_ord,l)-band_cha(1:lit_ord,p)-band_cha(1:lit_ord,p+2)) < toln(1:lit_ord))) exit
              end do
              write(qe_band,'(a,$)') 'bands:'
              write(qe_band,'(2i4,$)') p, p+2
              write (qe_band,'(a5,a,a,i2,a,a5,a,3i4)') ' ', 'Irrep: ' , '#', l, '*', ' ', 'accidental degeneracy:', p, p+1, p+2
              do j = 1, lit_ord
                write(qe_band,'("(",f6.3,","f6.3,") ",$)') band_cha(j,p)+band_cha(j,p+2)+(1.d-8,1.d-8)
              end do
              write (qe_band, *)
              do l  = 1, nirr_t
                if(all(abs(lit_cha_t(1:lit_ord,l)-band_cha(1:lit_ord,p+1)) < toln(1:lit_ord))) exit
              end do
              write(qe_band,'(a,$)') 'bands:'
              write(qe_band,'(i4,$)') p+1
              write (qe_band,'(a5,a,a,i2,a,a5,a,3i4)') ' ', 'Irrep: ' , '#', l, '*', ' ', 'accidental degeneracy:', p, p+1, p+2
              do j = 1, lit_ord
                write(qe_band,'("(",f6.3,","f6.3,") ",$)') band_cha(j,p+1)+(1.d-8,1.d-8)
              end do
              write (qe_band, *)
            end if
          end if
        end if !3
        if(i-p > 3) then
          write(qe_band,'(a,$)') 'bands:'
          do ibnd = p, i-1
            write(qe_band,'(i4,$)') ibnd
          end do
          write (qe_band,'(a5,a,a,i2,a,a5,a)') ' ', 'Irrep: ' , '#', l, '*', ' ', 'accidental degeneracy:'
          do j = 1, lit_ord
            write(qe_band,'("(",f6.3,","f6.3,") ",$)') cha(j)+(1.d-8,1.d-8)
          end do
          write (qe_band, *)
        end if
      else
        write(qe_band,'(a,$)') 'bands:'
        do ibnd = p, i-1
          write(qe_band,'(i4,$)') ibnd
        end do
        write (qe_band,'(a5,a,a,i2,a)') ' ', 'Irrep: ' , '#', l, '*'
        do j = 1, lit_ord
          write(qe_band,'("(",f6.3,","f6.3,") ",$)') cha(j)+(1.d-8,1.d-8)
        end do
        write (qe_band, *)
      end if
      if(k == nbnd) exit
      p = i
    end do
    write (qe_band, *)
    write (qe_band, '(a/)') REPEAT('=', 100)
    !
    deallocate(lit_spa, lit_symb, lit_axis, lit_ele, lit_cha, &
               lit_spa_t, lit_symb_t, lit_axis_t, lit_ele_t, lit_cha_t, lit_cha_t2, gvector, Cg, band_cha, cha, toln)
    close(qe_band)
  end do ! ik
  !
  deallocate(kpoints,eva)
  !
  return
end subroutine qe_bands