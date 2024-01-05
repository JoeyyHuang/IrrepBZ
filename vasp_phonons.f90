subroutine vasp_phonons()
!=========================================================================
! group analysis of VASP phonons
!=========================================================================
  use variables, only : a_prim, natom_prim, atom_prim_crys
  use constants, only : tpi, im, tol, tol3, vasp_phonon
  use mat_trans, only : mat_inv
  implicit none
!-------------------------------------------------------------------------
! local variables
  integer(kind=4) :: i, j, k, p, l, ia, ib, ic, iq, nq, iph, natom, nfreq, lit_ord, lit_ord_t, nirr, nirr_t, nirr_t2
  integer(kind=4),allocatable :: lit_spa(:), lit_spa_t(:), lit_axis(:,:), lit_axis_t(:,:)
  real(kind=8) :: q_crys(3), atom_trans(3), qpoints(3)
  real(kind=8),allocatable :: freq(:), lit_ele(:,:,:), lit_ele_t(:,:,:), toln(:)
  complex(kind=8),allocatable :: evc(:,:), lit_cha(:,:), lit_cha_t(:,:), lit_cha_t2(:,:), T(:,:,:), phonon_cha(:,:), cha(:)
  character(len=132) :: line
  character(len=20) fname
  character(len=10) qlabel
  character(len=2),allocatable :: lit_symb(:), lit_symb_t(:)
  logical :: exst_band, exst_qpoints
!-------------------------------------------------------------------------
  exst_band = .false.; exst_qpoints = .false.
  inquire(file = 'vasp_phonon/band.yaml', exist = exst_band)
  inquire(file = 'vasp_phonon/qpoints.yaml', exist = exst_qpoints)
  if(exst_band .and. exst_qpoints) then
    open(unit=vasp_phonon, file='vasp_phonon/phonon_ik1.txt', status='replace')
    write (vasp_phonon,'(a)') 'Error: detect the presence of both band.yaml and qpoints_yaml,'
    write (vasp_phonon,'(a)') '       make sure only one of them exists!'
    close(vasp_phonon)
    return
  end if
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! read q-points, frequencies and eigenvectors from band.yaml or qpoints.yaml
  if(exst_band) open(unit=10, file='vasp_phonon/band.yaml', status='old')
  if(exst_qpoints) open(unit=10, file='vasp_phonon/qpoints.yaml', status='old')
  do while(.true.)
    read(10,'(a132)') line
    p = index(line(1:7),'nqpoint',back=.false.,kind=4)
    if(p /= 0) then
      read(line(9:),*) nq
    end if
    p = index(line,'natom',back=.false.,kind=4)
    if(p /= 0) then
      read(line(7:),*) natom
      nfreq = 3*natom
      allocate(freq(nfreq), evc(nfreq,nfreq))
    end if
    p = index(line,'phonon',back=.false.,kind=4)
    if(p /= 0) then
      do iq = 1, nq
        if(iq < 10) then
          write(fname,'(a9,i1.1,a4)') 'phonon_iq', iq, '.txt'
        else if(iq < 100) then
          write(fname,'(a9,i2.2,a4)') 'phonon_iq', iq, '.txt'
        else if(iq < 1000) then
          write(fname,'(a9,i3.3,a4)') 'phonon_iq', iq, '.txt'
        end if
        open(unit=vasp_phonon, file='vasp_phonon/'//trim(fname), status='replace')
        write (vasp_phonon,'(a,i3)') 'Number of q points:', nq
        write (vasp_phonon,'(a,i3/)') 'Number of phonon modes at each q point:', nfreq
        !
        read(10,'(15x,f13.7,1x,f13.7,1x,f13.7)') qpoints(:)
        if(exst_band) read(10,*)
        read(10,*)
        do i = 1, nfreq
          read(10,*)
          read(10,'(14x,f16.10)') freq(i)
          read(10,*)
          do j = 1, natom
            read(10,*)
            do k = 1, 3
              read(10,'(9x,f18.14,1x,f18.14)') evc((j-1)*3+k,i)
            end do !k
          end do !j
        end do ! i
        read(10,*)
        !
        ! Irreps of little group
        q_crys = qpoints
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
        write (vasp_phonon,'(a,a,a,3f10.5)') 'This is the ', trim(adjustl(qlabel)), ' q-point (cryst. coord.):', q_crys
        !!!!!!!!!!!!!!!!!!!!! vector representation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! ordinary little group with all unitary operations
        call little_order(q_crys, lit_ord)
        allocate(lit_spa(lit_ord),lit_ele(3,4,lit_ord), lit_symb(lit_ord), lit_axis(3,lit_ord), toln(lit_ord))
        call little_group(vasp_phonon, q_crys, lit_ord, lit_spa, lit_ele, lit_symb, lit_axis)
        write (vasp_phonon, '(a/)') REPEAT('=', 100)
        allocate(lit_cha(lit_ord,lit_ord))
        call small_rep(vasp_phonon,q_crys,lit_ord,lit_ele,lit_symb,lit_axis,nirr,lit_cha,.false.)
        call Herring(vasp_phonon, q_crys, lit_ord, lit_spa, lit_ele, lit_symb, nirr, lit_cha, .false.)
        ! extended little group when consider the anti-unitary operation of TRS as en element
        call little_order_t(q_crys, lit_ord, lit_ord_t)
        allocate(lit_spa_t(lit_ord_t),lit_ele_t(3,4,lit_ord_t),lit_symb_t(lit_ord_t),lit_axis_t(3,lit_ord_t))
        call little_group_t(vasp_phonon, q_crys, lit_ord, lit_ord_t, lit_spa_t, lit_ele_t, lit_symb_t, lit_axis_t)
        allocate(lit_cha_t(lit_ord_t,lit_ord_t),lit_cha_t2(lit_ord_t,lit_ord_t))
        call small_rep_t(vasp_phonon,q_crys,lit_ord,lit_ord_t,lit_spa_t,lit_ele_t,lit_symb_t,lit_axis_t,nirr_t,lit_cha_t,.false.)
        write (vasp_phonon, '(a/)') REPEAT('=', 100)
        !!!!!!!!!!!!!!!!!!!!!!! spinor representation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! ordinary little group with all unitary operations
        call small_rep(vasp_phonon,q_crys,lit_ord,lit_ele,lit_symb,lit_axis,nirr,lit_cha,.true.)
        call Herring(vasp_phonon, q_crys, lit_ord, lit_spa, lit_ele, lit_symb, nirr, lit_cha, .true.)
        ! extended little group when consider the anti-unitary operation of TRS as en element
        call small_rep_t(vasp_phonon,q_crys,lit_ord,lit_ord_t,lit_spa_t,lit_ele_t,lit_symb_t,lit_axis_t,nirr_t2,lit_cha_t2,.true.)
        write (vasp_phonon, '(a/)') REPEAT('=', 100)
        !
        ! calculate Irreps of VASP phonon modes
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
        ! print Irreps of VASP phonon modes
        write (vasp_phonon,'(a)') '   <<  Irreps of Phonon Modes  >>'
        write (vasp_phonon,'(a)') 'Determine the Irrep by comparing the mode characters &
                                  &with the small Reps of the little group in "Scenario Without Considering SOC"'
        write (vasp_phonon,'(a/)') 'For lattice dynamical equation, the time-reversal symmetry is intrinsically holded'
        write (vasp_phonon,'(a)') 'Characters of phonon modes associated with each symmetry opeartion'
        do i = 1, lit_ord
          write (vasp_phonon,'(a5,a2,a,i2,a,a5,$)') ' ', lit_symb_t(i), '[',lit_spa_t(i),']', ' '
        end do
        write (vasp_phonon,*)
        p = 1; k = 0; toln = tol
        do while(.true.)
          cha = (0.d0,0.d0)
          do i = p, nfreq
            if(abs(freq(i)-freq(p)) < tol) then
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
                write(vasp_phonon,'(a,$)') 'modes:'
                write(vasp_phonon,'(i4,$)') iph
                write (vasp_phonon,'(a5,a,a,i2,a,a5,a,2i4)') ' ','Irrep: ','#',l,'*',' ','accidental degeneracy:',p,p+1
                do j = 1, lit_ord
                  write(vasp_phonon,'("(",f6.3,","f6.3,") ",$)') phonon_cha(j,iph)+(1.d-5,1.d-5)
                end do
                write (vasp_phonon, *)
              end do ! iph
            end if ! 2
            if(i-p == 3) then
              do l  = 1, nirr_t
                if(all(abs(lit_cha_t(1:lit_ord,l)-phonon_cha(1:lit_ord,p)) < toln(1:lit_ord))) exit
              end do
              if(l <= nirr_t) then
                write(vasp_phonon,'(a,$)') 'modes:'
                write(vasp_phonon,'(i4,$)') p
                write (vasp_phonon,'(a5,a,a,i2,a,a5,a,3i4)') ' ','Irrep: ' ,'#',l,'*',' ','accidental degeneracy:',p,p+1,p+2
                do j = 1, lit_ord
                  write(vasp_phonon,'("(",f6.3,","f6.3,") ",$)') phonon_cha(j,p)+(1.d-5,1.d-5)
                end do
                write (vasp_phonon, *)
                do l  = 1, nirr_t
                  if(all(abs(lit_cha_t(1:lit_ord,l)-phonon_cha(1:lit_ord,p+1)) < toln(1:lit_ord))) exit
                end do
                if(l <= nirr_t) then
                  write(vasp_phonon,'(a,$)') 'modes:'
                  write(vasp_phonon,'(i4,$)') p+1
                  write (vasp_phonon,'(a5,a,a,i2,a,a5,a,3i4)') ' ','Irrep: ','#',l,'*',' ','accidental degeneracy:',p,p+1,p+2
                  do j = 1, lit_ord
                    write(vasp_phonon,'("(",f6.3,","f6.3,") ",$)') phonon_cha(j,p+1)+(1.d-5,1.d-5)
                  end do
                  write (vasp_phonon, *)
                  do l  = 1, nirr_t
                    if(all(abs(lit_cha_t(1:lit_ord,l)-phonon_cha(1:lit_ord,p+2)) < toln(1:lit_ord))) exit
                  end do
                  write(vasp_phonon,'(a,$)') 'modes:'
                  write(vasp_phonon,'(i4,$)') p+2
                  write (vasp_phonon,'(a5,a,a,i2,a,a5,a,3i4)') ' ','Irrep: ','#',l,'*',' ','accidental degeneracy:',p,p+1,p+2
                  do j = 1, lit_ord
                    write(vasp_phonon,'("(",f6.3,","f6.3,") ",$)') phonon_cha(j,p+2)+(1.d-5,1.d-5)
                  end do
                  write (vasp_phonon, *)
                else
                  do l  = 1, nirr_t
                    if(all(abs(lit_cha_t(1:lit_ord,l)-phonon_cha(1:lit_ord,p+1)-phonon_cha(1:lit_ord,p+2))<toln(1:lit_ord))) exit
                  end do
                  write(vasp_phonon,'(a,$)') 'modes:'
                  write(vasp_phonon,'(2i4,$)') p+1, p+2
                  write (vasp_phonon,'(a5,a,a,i2,a,a5,a,3i4)') ' ','Irrep: ','#',l,'*',' ','accidental degeneracy:',p,p+1,p+2
                  do j = 1, lit_ord
                    write(vasp_phonon,'("(",f6.3,","f6.3,") ",$)') phonon_cha(j,p+1)+phonon_cha(j,p+2)+(1.d-5,1.d-5)
                  end do
                  write (vasp_phonon, *)
                end if
              else
                do l  = 1, nirr_t
                  if(all(abs(lit_cha_t(1:lit_ord,l)-phonon_cha(1:lit_ord,p)-phonon_cha(1:lit_ord,p+1)) < toln(1:lit_ord))) exit
                end do
                if(l <= nirr_t) then
                  write(vasp_phonon,'(a,$)') 'modes:'
                  write(vasp_phonon,'(2i4,$)') p, p+1
                  write (vasp_phonon,'(a5,a,a,i2,a,a5,a,3i4)') ' ','Irrep: ','#',l,'*',' ','accidental degeneracy:',p,p+1,p+2
                  do j = 1, lit_ord
                    write(vasp_phonon,'("(",f6.3,","f6.3,") ",$)') phonon_cha(j,p)+phonon_cha(j,p+1)+(1.d-5,1.d-5)
                  end do
                  write (vasp_phonon, *)
                  do l  = 1, nirr_t
                    if(all(abs(lit_cha_t(1:lit_ord,l)-phonon_cha(1:lit_ord,p+2)) < toln(1:lit_ord))) exit
                  end do
                  write(vasp_phonon,'(a,$)') 'modes:'
                  write(vasp_phonon,'(i4,$)') p+2
                  write (vasp_phonon,'(a5,a,a,i2,a,a5,a,3i4)') ' ','Irrep: ','#',l,'*',' ','accidental degeneracy:',p,p+1,p+2
                  do j = 1, lit_ord
                    write(vasp_phonon,'("(",f6.3,","f6.3,") ",$)') phonon_cha(j,p+2)+(1.d-5,1.d-5)
                  end do
                  write (vasp_phonon, *)
                else
                  do l  = 1, nirr_t
                    if(all(abs(lit_cha_t(1:lit_ord,l)-phonon_cha(1:lit_ord,p)-phonon_cha(1:lit_ord,p+2)) < toln(1:lit_ord))) exit
                  end do
                  write(vasp_phonon,'(a,$)') 'modes:'
                  write(vasp_phonon,'(2i4,$)') p, p+2
                  write (vasp_phonon,'(a5,a,a,i2,a,a5,a,3i4)') ' ','Irrep: ','#',l,'*',' ','accidental degeneracy:',p,p+1,p+2
                  do j = 1, lit_ord
                    write(vasp_phonon,'("(",f6.3,","f6.3,") ",$)') phonon_cha(j,p)+phonon_cha(j,p+2)+(1.d-5,1.d-5)
                  end do
                  write (vasp_phonon, *)
                  do l  = 1, nirr_t
                    if(all(abs(lit_cha_t(1:lit_ord,l)-phonon_cha(1:lit_ord,p+1)) < toln(1:lit_ord))) exit
                  end do
                  write(vasp_phonon,'(a,$)') 'modes:'
                  write(vasp_phonon,'(i4,$)') p+1
                  write (vasp_phonon,'(a5,a,a,i2,a,a5,a,3i4)') ' ','Irrep: ','#',l,'*',' ','accidental degeneracy:',p,p+1,p+2
                  do j = 1, lit_ord
                    write(vasp_phonon,'("(",f6.3,","f6.3,") ",$)') phonon_cha(j,p+1)+(1.d-5,1.d-5)
                  end do
                  write (vasp_phonon, *)
                end if
              end if
            end if !3
            if(i-p > 3) then
              write(vasp_phonon,'(a,$)') 'modes:'
              do iph = p, i-1
                write(vasp_phonon,'(i4,$)') iph
              end do
              write (vasp_phonon,'(a5,a,a,i2,a,a5,a)') ' ', 'Irrep: ', '#', l, '*', ' ', 'accidental degeneracy:'
              do j = 1, lit_ord
                write(vasp_phonon,'("(",f6.3,","f6.3,") ",$)') cha(j)+(1.d-5,1.d-5)
              end do
              write (vasp_phonon, *)
            end if
          else  ! l <= nirr_t
            write(vasp_phonon,'(a,$)') 'modes:'
            do iph = p, i-1
              write(vasp_phonon,'(i4,$)') iph
            end do
            write (vasp_phonon,'(a5,a,a,i2,a)') ' ', 'Irrep: ', '#', l, '*'
            do j = 1, lit_ord
              write(vasp_phonon,'("(",f6.3,","f6.3,") ",$)') cha(j)+(1.d-5,1.d-5)
            end do
            write (vasp_phonon, *)
          end if ! l
          if(k == nfreq) exit
          p = i
        end do ! while
        write (vasp_phonon, *)
        write (vasp_phonon, '(a/)') REPEAT('=', 100)
        !
        deallocate(lit_spa,lit_spa_t,lit_symb,lit_symb_t,lit_axis,lit_axis_t,lit_ele,lit_ele_t,lit_cha,lit_cha_t,&
                   lit_cha_t2,T,phonon_cha,cha,toln)
        close(vasp_phonon)
      end do ! iq
      exit
    end if ! match 'phonon'
  end do ! while
  close(10)
  !
  deallocate(freq,evc)
  !
  return
end subroutine vasp_phonons