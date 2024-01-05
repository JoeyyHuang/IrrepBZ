subroutine qe_phonons()
!=========================================================================
! group analysis of QE phonons
!=========================================================================
  use variables, only : a_prim, natom_prim, atom_prim_crys
  use constants, only : tpi, im, tol, tol3, qe_phonon
  use mat_trans, only : mat_inv
  use read_binary, only : read_a_wfc
  implicit none
!-------------------------------------------------------------------------
! local variables
  integer(kind=4) :: i, j, k, l, p, ia, ib, ic, iq, nq, iph, nfreq, lit_ord, lit_ord_t, nirr, nirr_t, nirr_t2
  integer(kind=4),allocatable :: lit_spa(:), lit_spa_t(:), lit_axis(:,:), lit_axis_t(:,:)
  real(kind=8) :: q_crys(3), atom_trans(3)
  real(kind=8),allocatable :: qpoints(:,:), freq(:,:), lit_ele(:,:,:), lit_ele_t(:,:,:), toln(:)
  complex(kind=8),allocatable :: evc(:,:), lit_cha(:,:), lit_cha_t(:,:), lit_cha_t2(:,:), T(:,:,:), phonon_cha(:,:), cha(:)
  character(len=2),allocatable :: lit_symb(:), lit_symb_t(:)
  character(len=10) qlabel
  character(len=20) fname
!-------------------------------------------------------------------------
  ! read q-points and frequencies from matdyn.freq
  open(unit=10, file='qe_phonon/matdyn.freq', status='old')
  read(10,'(12x,i4,6x,i4)') nfreq, nq
  allocate(qpoints(3,nq),freq(nfreq,nq),evc(nfreq,nfreq))
  do iq = 1, nq
    read(10,'(10x,3f10.6)') qpoints(:,iq)
    read(10,'(6f10.4)') (freq(i,iq), i=1,nfreq)
  end do
  close(10)
  ! convert to cryst. coord.
  qpoints = matmul(transpose(a_prim),qpoints) / sqrt(a_prim(1,1)**2+a_prim(2,1)**2+a_prim(3,1)**2)
  !
  ! read eigenvectors from matdyn.modes
  open(unit=10, file='qe_phonon/matdyn.eig', status='old')
  ! loops on q-points
  do iq = 1, nq
    if(iq < 10) then
      write(fname,'(a9,i1.1,a4)') 'phonon_iq', iq, '.txt'
    else if(iq < 100) then
      write(fname,'(a9,i2.2,a4)') 'phonon_iq', iq, '.txt'
    else if(iq < 1000) then
      write(fname,'(a9,i3.3,a4)') 'phonon_iq', iq, '.txt'
    end if
    open(unit=qe_phonon, file='qe_phonon/'//trim(fname), status='replace')
    write (qe_phonon,'(a,i3)') 'Number of q points:', nq
    write (qe_phonon,'(a,i3/)') 'Number of phonon modes at each q point:', nfreq
    !
    do i = 1, 4
      read(10,*)
    end do
    do i = 1, nfreq
      read(10,*)
      do j = 1, nfreq/3
        read(10,'(2x,3(f10.6,1x,f10.6,3x))') (evc((j-1)*3+k,i),k=1,3) ! cart. coord.
      end do
    end do
    read(10,*)
    !
    ! Irreps of little group
    q_crys = qpoints(:,iq)
    if(iq == 1) then
      write(qlabel,'(a)') '1st'
    else if(iq == 2) then
      write(qlabel,'(a)') '2nd'
    else if(iq == 3) then
      write(qlabel,'(a)') '3rd'
    else if(iq < 10) then
      write(qlabel,'(i1,a)') iq, 'th'
    else if(iq < 100) then
      write(qlabel,'(i2,a)') iq, 'th'
    else if(iq < 1000) then
      write(qlabel,'(i3,a)') iq, 'th'
    end if
    write (qe_phonon,'(a,a,a,3f10.5)') 'This is the ', trim(adjustl(qlabel)), ' q-point (cryst. coord.):', q_crys
    !!!!!!!!!!!!!!!!!!!!! vector representation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ordinary little group with all unitary operations
    call little_order(q_crys, lit_ord)
    allocate(lit_spa(lit_ord),lit_ele(3,4,lit_ord), lit_symb(lit_ord), lit_axis(3,lit_ord), toln(lit_ord))
    call little_group(qe_phonon, q_crys, lit_ord, lit_spa, lit_ele, lit_symb, lit_axis)
    write (qe_phonon, '(a/)') REPEAT('=', 100)
    allocate(lit_cha(lit_ord,lit_ord))
    call small_rep(qe_phonon,q_crys,lit_ord,lit_ele,lit_symb,lit_axis,nirr,lit_cha,.false.)
    call Herring(qe_phonon, q_crys, lit_ord, lit_spa, lit_ele, lit_symb, nirr, lit_cha, .false.)
    ! extended little group when consider the anti-unitary operation of TRS as en element
    call little_order_t(q_crys, lit_ord, lit_ord_t)
    allocate(lit_spa_t(lit_ord_t),lit_ele_t(3,4,lit_ord_t),lit_symb_t(lit_ord_t),lit_axis_t(3,lit_ord_t))
    call little_group_t(qe_phonon, q_crys, lit_ord, lit_ord_t, lit_spa_t, lit_ele_t, lit_symb_t, lit_axis_t)
    allocate(lit_cha_t(lit_ord_t,lit_ord_t),lit_cha_t2(lit_ord_t,lit_ord_t))
    call small_rep_t(qe_phonon,q_crys,lit_ord,lit_ord_t,lit_spa_t,lit_ele_t,lit_symb_t,lit_axis_t,nirr_t,lit_cha_t,.false.)
    write (qe_phonon, '(a/)') REPEAT('=', 100)
    !!!!!!!!!!!!!!!!!!!!!!! spinor representation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ordinary little group with all unitary operations
    call small_rep(qe_phonon,q_crys,lit_ord,lit_ele,lit_symb,lit_axis,nirr,lit_cha,.true.)
    call Herring(qe_phonon, q_crys, lit_ord, lit_spa, lit_ele, lit_symb, nirr, lit_cha, .true.)
    ! extended little group when consider the anti-unitary operation of TRS as en element
    call small_rep_t(qe_phonon,q_crys,lit_ord,lit_ord_t,lit_spa_t,lit_ele_t,lit_symb_t,lit_axis_t,nirr_t2,lit_cha_t2,.true.)
    write (qe_phonon, '(a/)') REPEAT('=', 100)
    !
    ! calculate Irreps of QE phonon modes
    allocate(T(nfreq,nfreq,lit_ord),phonon_cha(lit_ord,nfreq),cha(lit_ord))
    T = (0.d0,0.d0)
    do i = 1, lit_ord
      do j = 1, natom_prim
        atom_trans(:) = matmul(lit_ele(:,1:3,i),atom_prim_crys(:,j)) + lit_ele(:,4,i)
        do k = 1, natom_prim
          do ia = -3, 3
            do ib = -3, 3
              do ic = -3, 3
                if(all(abs(atom_trans(:) - atom_prim_crys(:,k) - dble((/ia,ib,ic/))) < tol3)) goto 100
              end do
            end do
          end do
        end do ! k
    100 T(3*(k-1)+1:3*k,3*(j-1)+1:3*j,i) = matmul(a_prim,matmul(lit_ele(:,1:3,i),mat_inv(a_prim)))*&
                                           exp(-im*tpi*dot_product(q_crys,dble((/ia,ib,ic/))))
      end do ! j
      do j = 1, nfreq
        phonon_cha(i,j) = dot_product(evc(:,j),matmul(T(:,:,i),evc(:,j)))*exp(-im*tpi*dot_product(q_crys,lit_ele(:,4,i)))
      end do
    end do ! i
    !
    ! print Irreps of QE phonon modes
    write (qe_phonon,'(a)') '   <<  Irreps of Phonon Modes  >>'
    write (qe_phonon,'(a)') 'Determine the Irrep by comparing the mode characters &
                            &with the small Reps of the little group in "Scenario Without Considering SOC"'
    write (qe_phonon,'(a/)') 'For lattice dynamical equation, the time-reversal symmetry is intrinsically holded'
    write (qe_phonon,'(a)') 'Characters of phonon modes associated with each symmetry opeartion'
    do i = 1, lit_ord
      write (qe_phonon,'(a5,a2,a,i2,a,a5,$)') ' ', lit_symb_t(i), '[',lit_spa_t(i),']', ' '
    end do
    write (qe_phonon,*)
    p = 1; k = 0; toln = tol
    do while(.true.)
      cha = (0.d0,0.d0)
      do i = p, nfreq
        if(abs(freq(i,iq)-freq(p,iq)) < tol) then
          cha = cha + phonon_cha(:,i)
          if(i == nfreq) k = nfreq
        else
          exit
        end if
      end do
      ! compare the modes characters with Irreps of little group
      do l = 1, nirr_t
        if(all(abs(lit_cha_t(1:lit_ord,l)-cha(1:lit_ord)) < toln(1:lit_ord))) exit
      end do
      !
      ! TODO: accidental degeneracy
      if(l > nirr_t) then
        if(i-p == 2) then
          do iph = p, p+1
            do l  = 1, nirr_t
              if(all(abs(lit_cha_t(1:lit_ord,l)-cha(1:lit_ord)) < toln(1:lit_ord))) exit
            end do
            write(qe_phonon,'(a,$)') 'modes:'
            write(qe_phonon,'(i4,$)') iph
            write (qe_phonon,'(a5,a,a,i2,a,a5,a,2i4)') ' ', 'Irrep: ', '#', l, '*', ' ', 'accidental degeneracy:', p, p+1
            do j = 1, lit_ord
              write(qe_phonon,'("(",f6.3,","f6.3,") ",$)') phonon_cha(j,iph)+(1.d-5,1.d-5)
            end do
            write (qe_phonon, *)
          end do ! iph
        end if ! 2
        if(i-p == 3) then
          do l  = 1, nirr_t
            if(all(abs(lit_cha_t(1:lit_ord,l)-phonon_cha(1:lit_ord,p)) < toln(1:lit_ord))) exit
          end do
          if(l <= nirr_t) then
            write(qe_phonon,'(a,$)') 'modes:'
            write(qe_phonon,'(i4,$)') p
            write (qe_phonon,'(a5,a,a,i2,a,a5,a,3i4)') ' ', 'Irrep: ', '#', l, '*', ' ', 'accidental degeneracy:', p, p+1, p+2
            do j = 1, lit_ord
              write(qe_phonon,'("(",f6.3,","f6.3,") ",$)') phonon_cha(j,p)+(1.d-5,1.d-5)
            end do
            write (qe_phonon, *)
            do l  = 1, nirr_t
              if(all(abs(lit_cha_t(1:lit_ord,l)-phonon_cha(1:lit_ord,p+1)) < toln(1:lit_ord))) exit
            end do
            if(l <= nirr_t) then
              write(qe_phonon,'(a,$)') 'modes:'
              write(qe_phonon,'(i4,$)') p+1
              write (qe_phonon,'(a5,a,a,i2,a,a5,a,3i4)') ' ', 'Irrep: ', '#', l, '*', ' ', 'accidental degeneracy:', p, p+1, p+2
              do j = 1, lit_ord
                write(qe_phonon,'("(",f6.3,","f6.3,") ",$)') phonon_cha(j,p+1)+(1.d-5,1.d-5)
              end do
              write (qe_phonon, *)
              do l  = 1, nirr_t
                if(all(abs(lit_cha_t(1:lit_ord,l)-phonon_cha(1:lit_ord,p+2)) < toln(1:lit_ord))) exit
              end do
              write(qe_phonon,'(a,$)') 'modes:'
              write(qe_phonon,'(i4,$)') p+2
              write (qe_phonon,'(a5,a,a,i2,a,a5,a,3i4)') ' ', 'Irrep: ', '#', l, '*', ' ', 'accidental degeneracy:', p, p+1, p+2
              do j = 1, lit_ord
                write(qe_phonon,'("(",f6.3,","f6.3,") ",$)') phonon_cha(j,p+2)+(1.d-5,1.d-5)
              end do
              write (qe_phonon, *)
            else
              do l  = 1, nirr_t
                if(all(abs(lit_cha_t(1:lit_ord,l)-phonon_cha(1:lit_ord,p+1)-phonon_cha(1:lit_ord,p+2)) < toln(1:lit_ord))) exit
              end do
              write(qe_phonon,'(a,$)') 'modes:'
              write(qe_phonon,'(2i4,$)') p+1, p+2
              write (qe_phonon,'(a5,a,a,i2,a,a5,a,3i4)') ' ', 'Irrep: ', '#', l, '*', ' ', 'accidental degeneracy:', p, p+1, p+2
              do j = 1, lit_ord
                write(qe_phonon,'("(",f6.3,","f6.3,") ",$)') phonon_cha(j,p+1)+phonon_cha(j,p+2)+(1.d-5,1.d-5)
              end do
              write (qe_phonon, *)
            end if
          else
            do l  = 1, nirr_t
              if(all(abs(lit_cha_t(1:lit_ord,l)-phonon_cha(1:lit_ord,p)-phonon_cha(1:lit_ord,p+1)) < toln(1:lit_ord))) exit
            end do
            if(l <= nirr_t) then
              write(qe_phonon,'(a,$)') 'modes:'
              write(qe_phonon,'(2i4,$)') p, p+1
              write (qe_phonon,'(a5,a,a,i2,a,a5,a,3i4)') ' ', 'Irrep: ', '#', l, '*', ' ', 'accidental degeneracy:', p, p+1, p+2
              do j = 1, lit_ord
                write(qe_phonon,'("(",f6.3,","f6.3,") ",$)') phonon_cha(j,p)+phonon_cha(j,p+1)+(1.d-5,1.d-5)
              end do
              write (qe_phonon, *)
              do l  = 1, nirr_t
                if(all(abs(lit_cha_t(1:lit_ord,l)-phonon_cha(1:lit_ord,p+2)) < toln(1:lit_ord))) exit
              end do
              write(qe_phonon,'(a,$)') 'modes:'
              write(qe_phonon,'(i4,$)') p+2
              write (qe_phonon,'(a5,a,a,i2,a,a5,a,3i4)') ' ', 'Irrep: ', '#', l, '*', ' ', 'accidental degeneracy:', p, p+1, p+2
              do j = 1, lit_ord
                write(qe_phonon,'("(",f6.3,","f6.3,") ",$)') phonon_cha(j,p+2)+(1.d-5,1.d-5)
              end do
              write (qe_phonon, *)
            else
              do l  = 1, nirr_t
                if(all(abs(lit_cha_t(1:lit_ord,l)-phonon_cha(1:lit_ord,p)-phonon_cha(1:lit_ord,p+2)) < toln(1:lit_ord))) exit
              end do
              write(qe_phonon,'(a,$)') 'modes:'
              write(qe_phonon,'(2i4,$)') p, p+2
              write (qe_phonon,'(a5,a,a,i2,a,a5,a,3i4)') ' ', 'Irrep: ', '#', l, '*', ' ', 'accidental degeneracy:', p, p+1, p+2
              do j = 1, lit_ord
                write(qe_phonon,'("(",f6.3,","f6.3,") ",$)') phonon_cha(j,p)+phonon_cha(j,p+2)+(1.d-5,1.d-5)
              end do
              write (qe_phonon, *)
              do l  = 1, nirr_t
                if(all(abs(lit_cha_t(1:lit_ord,l)-phonon_cha(1:lit_ord,p+1)) < toln(1:lit_ord))) exit
              end do
              write(qe_phonon,'(a,$)') 'modes:'
              write(qe_phonon,'(i4,$)') p+1
              write (qe_phonon,'(a5,a,a,i2,a,a5,a,3i4)') ' ', 'Irrep: ', '#', l, '*', ' ', 'accidental degeneracy:', p, p+1, p+2
              do j = 1, lit_ord
                write(qe_phonon,'("(",f6.3,","f6.3,") ",$)') phonon_cha(j,p+1)+(1.d-5,1.d-5)
              end do
              write (qe_phonon, *)
            end if
          end if
        end if !3
        if(i-p > 3) then
          write(qe_phonon,'(a,$)') 'modes:'
          do iph = p, i-1
            write(qe_phonon,'(i4,$)') iph
          end do
          write (qe_phonon,'(a5,a,a,i2,a,a5,a)') ' ', 'Irrep: ', '#', l, '*', ' ', 'accidental degeneracy:'
          do j = 1, lit_ord
            write(qe_phonon,'("(",f6.3,","f6.3,") ",$)') cha(j)+(1.d-5,1.d-5)
          end do
          write (qe_phonon, *)
        end if
      else
        write(qe_phonon,'(a,$)') 'modes:'
        do iph = p, i-1
          write(qe_phonon,'(i4,$)') iph
        end do
        write (qe_phonon,'(a5,a,a,i2,a)') ' ', 'Irrep: ', '#', l, '*'
        do j = 1, lit_ord
          write(qe_phonon,'("(",f6.3,","f6.3,") ",$)') cha(j)+(1.d-5,1.d-5)
        end do
        write (qe_phonon, *)
      end if
      if(k == nfreq) exit
      p = i
    end do
    write (qe_phonon, *)
    write (qe_phonon, '(a/)') REPEAT('=', 100)
    !
    deallocate(lit_spa,lit_spa_t,lit_symb,lit_symb_t,lit_axis,lit_axis_t,lit_ele,lit_ele_t,&
               lit_cha,lit_cha_t,lit_cha_t2,T,phonon_cha,cha,toln)
    !
    close(qe_phonon)
  end do ! iq
  close(10)
  !
  deallocate(qpoints,freq,evc)
  !
  return
end subroutine qe_phonons