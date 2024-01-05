subroutine vasp_bands()
!=========================================================================
! group analysis of VASP electronic bands
!=========================================================================
  use variables, only : a_prim
  use constants, only : pi, c, tpi, im, tol, tol3, vasp_band
  use mat_trans, only : mat_inv
  implicit none
!-------------------------------------------------------------------------
! local variables
  integer(kind=4) :: i, j, k, l, p, ia, ib, ic, lit_ord, lit_ord_t, nirr, nirr_t, nirr_t2, npmax, iband, nrecl, &
                     iost, nspin, nprec, nplane, iplane, nwk, nband, nb1max, nb2max, nb3max, &
                     npmaxA, npmaxB, npmaxC, iwk, irec, ncnt, ig1, ig2, ig3, ig1p, ig2p, ig3p
  integer(kind=4),allocatable :: lit_spa(:), lit_spa_t(:), lit_axis(:,:), lit_axis_t(:,:), gvector(:,:), igall(:,:)
  real(kind=8) :: angle, xnrecl, xnspin, xnprec, xnwk, xnband, ecut, Vcell, a3mag, gtot, etot, &
                  b1mag, b2mag, b3mag, phi12, vmag, sinphi123, phi13, phi123, phi23, xnplane, &
                  nb1maxA,nb2maxA,nb3maxA,nb1maxB,nb2maxB,nb3maxB,nb1maxC,nb2maxC,nb3maxC
  real(kind=8) :: k_crys(3), k_crys2(3), gvector2(3), vec(3), &
                  a1(3), a2(3), a3(3), a2xa3(3), b1(3), b2(3), b3(3), vtmp(3), wk(3), sumkg(3)
  real(kind=8),allocatable :: lit_ele(:,:,:), lit_ele_t(:,:,:), occ(:), toln(:)
  complex(kind=4) :: norm
  complex(kind=8) :: spinor1(2), spinor2(2), SU2_mat(2,2)
  complex(kind=4),allocatable :: coeff(:,:)
  complex(kind=8),allocatable :: lit_cha(:,:), lit_cha_t(:,:), lit_cha_t2(:,:), band_cha(:,:), cha(:), cener(:)
  character(len=10) klabel
  character(len=14) fname
  character(len=2),allocatable :: lit_symb(:), lit_symb_t(:)
  logical :: soc
!-------------------------------------------------------------------------
  ! read WAVECAR to first check if soc is on
  nrecl = 24
  open(unit=10,file='vasp_band/WAVECAR',access='direct',recl=nrecl,iostat=iost,status='old')
  if (iost.ne.0) write(*,'(a,i1)') 'open error - iostat =',iost
  read(unit=10,rec=1) xnrecl,xnspin,xnprec
  close(unit=10)
  !
  nrecl=nint(xnrecl); nspin=nint(xnspin); nprec=nint(xnprec)
  if(nprec.eq.45210) then
    write(*,'(a)') '*** error - WAVECAR_double requires complex*16'
    stop
  end if
  !
  open(unit=10,file='vasp_band/WAVECAR',access='direct',recl=nrecl,iostat=iost,status='old')
  read(unit=10,rec=2) xnwk,xnband,ecut,(a1(j),j=1,3),(a2(j),j=1,3),(a3(j),j=1,3)
  nwk=nint(xnwk); nband=nint(xnband)
  allocate(occ(nband)); allocate(cener(nband))
  !
  call vcross(a2xa3,a2,a3)
  Vcell=dot_product(a1,a2xa3)
  a3mag=dsqrt(dot_product(a3,a3))
  call vcross(b1,a2,a3)
  call vcross(b2,a3,a1)
  call vcross(b3,a1,a2)
  b1=2.*pi*b1/Vcell
  b2=2.*pi*b2/Vcell
  b3=2.*pi*b3/Vcell
  b1mag=dsqrt(b1(1)**2+b1(2)**2+b1(3)**2)
  b2mag=dsqrt(b2(1)**2+b2(2)**2+b2(3)**2)
  b3mag=dsqrt(b3(1)**2+b3(2)**2+b3(3)**2)
  !
  phi12=acos((b1(1)*b2(1)+b1(2)*b2(2)+b1(3)*b2(3))/(b1mag*b2mag))
  call vcross(vtmp,b1,b2)
  vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
  sinphi123=(b3(1)*vtmp(1)+b3(2)*vtmp(2)+b3(3)*vtmp(3))/(vmag*b3mag)
  nb1maxA=(dsqrt(ecut*c)/(b1mag*abs(sin(phi12))))+1
  nb2maxA=(dsqrt(ecut*c)/(b2mag*abs(sin(phi12))))+1
  nb3maxA=(dsqrt(ecut*c)/(b3mag*abs(sinphi123)))+1
  npmaxA=nint(4.*pi*nb1maxA*nb2maxA*nb3maxA/3.)
  !   
  phi13=acos((b1(1)*b3(1)+b1(2)*b3(2)+b1(3)*b3(3))/(b1mag*b3mag))
  call vcross(vtmp,b1,b3)
  vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
  sinphi123=(b2(1)*vtmp(1)+b2(2)*vtmp(2)+b2(3)*vtmp(3))/(vmag*b2mag)
  phi123=abs(asin(sinphi123))
  nb1maxB=(dsqrt(ecut*c)/(b1mag*abs(sin(phi13))))+1
  nb2maxB=(dsqrt(ecut*c)/(b2mag*abs(sinphi123)))+1
  nb3maxB=(dsqrt(ecut*c)/(b3mag*abs(sin(phi13))))+1
  npmaxB=nint(4.*pi*nb1maxB*nb2maxB*nb3maxB/3.)
  !
  phi23=acos((b2(1)*b3(1)+b2(2)*b3(2)+b2(3)*b3(3))/(b2mag*b3mag))
  call vcross(vtmp,b2,b3)
  vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
  sinphi123=(b1(1)*vtmp(1)+b1(2)*vtmp(2)+b1(3)*vtmp(3))/(vmag*b1mag)
  phi123=abs(asin(sinphi123))
  nb1maxC=(dsqrt(ecut*c)/(b1mag*abs(sinphi123)))+1
  nb2maxC=(dsqrt(ecut*c)/(b2mag*abs(sin(phi23))))+1
  nb3maxC=(dsqrt(ecut*c)/(b3mag*abs(sin(phi23))))+1
  npmaxC=nint(4.*pi*nb1maxC*nb2maxC*nb3maxC/3.)
  !
  nb1max=nint(max(nb1maxA,nb1maxB,nb1maxC))
  nb2max=nint(max(nb2maxA,nb2maxB,nb2maxC))
  nb3max=nint(max(nb3maxA,nb3maxB,nb3maxC))
  !
  iwk = 1
  irec=3+(iwk-1)*(nband+1)
  read(unit=10,rec=irec) xnplane,(wk(j),j=1,3),(cener(iband),occ(iband),iband=1,nband)
  nplane=nint(xnplane)
  !
  ncnt = 0
  do ig3 = 0,2*nb3max
    ig3p = ig3
    if(ig3.gt.nb3max) ig3p=ig3-2*nb3max-1
    do ig2 = 0,2*nb2max
      ig2p = ig2
      if(ig2.gt.nb2max) ig2p=ig2-2*nb2max-1
      do ig1 = 0,2*nb1max
        ig1p = ig1
        if(ig1.gt.nb1max) ig1p=ig1-2*nb1max-1
        do j = 1,3
          sumkg(j)=(wk(1)+ig1p)*b1(j)+(wk(2)+ig2p)*b2(j)+(wk(3)+ig3p)*b3(j)
        end do
        gtot=sqrt(dot_product(sumkg,sumkg))
        etot=gtot**2/c
        if(etot.lt.ecut) ncnt=ncnt+1
      end do
    end do
  end do
  !
  if(2*ncnt == nplane) then 
    soc = .true.
    npmax=2*min0(npmaxA,npmaxB,npmaxC)
  else if(ncnt == nplane) then
    soc = .false.
    npmax=min0(npmaxA,npmaxB,npmaxC)
  end if
  allocate(igall(3,npmax),coeff(npmax,nband))
  !
  ! loops on k-points
  do iwk = 1, nwk
    if(iwk < 10) then
      write(fname,'(a7,i1.1,a4)') 'band_ik', iwk, '.txt'
    else if(iwk < 100) then
      write(fname,'(a7,i2.2,a4)') 'band_ik', iwk, '.txt'
    else if(iwk < 1000) then
      write(fname,'(a7,i3.3,a4)') 'band_ik', iwk, '.txt'
    end if
    open(unit=vasp_band, file='vasp_band/'//trim(fname), status='replace')
    if(.NOT. soc) write (vasp_band,'(a)') 'The electronic bands are calculated without considering SOC'
    if(soc) write (vasp_band,'(a)') 'The electronic bands are calculated with considering SOC'
    write (vasp_band,'(a,i3)') 'Number of k points:', nwk
    write (vasp_band,'(a,i3/)') 'Number of electronic bands at each k point:', nband
    !
    irec=3+(iwk-1)*(nband+1)
    read(unit=10,rec=irec) xnplane,(wk(j),j=1,3),(cener(iband),occ(iband),iband=1,nband)
    nplane=nint(xnplane)
    !
    ! gvectors
    ncnt = 0
    do ig3 = 0,2*nb3max
      ig3p = ig3
      if(ig3.gt.nb3max) ig3p=ig3-2*nb3max-1
      do ig2 = 0,2*nb2max
        ig2p = ig2
        if(ig2.gt.nb2max) ig2p=ig2-2*nb2max-1
        do ig1 = 0,2*nb1max
          ig1p = ig1
          if(ig1.gt.nb1max) ig1p=ig1-2*nb1max-1
          do j = 1,3
            sumkg(j)=(wk(1)+ig1p)*b1(j)+(wk(2)+ig2p)*b2(j)+(wk(3)+ig3p)*b3(j)
          end do
          gtot=sqrt(dot_product(sumkg,sumkg))
          etot=gtot**2/c
          if(etot.lt.ecut) then
            ncnt=ncnt+1
            igall(1,ncnt)=ig1p
            igall(2,ncnt)=ig2p
            igall(3,ncnt)=ig3p
          end if
        end do
      end do
    end do
    !
    allocate(gvector(3,ncnt))
    gvector = igall
    !
    do iband = 1, nband
      irec = irec + 1
      read(unit=10,rec=irec) (coeff(iplane,iband), iplane=1,nplane)
    end do
    ! normalization
    do iband = 1, nband
      norm = (0.d0,0.d0)
      do iplane = 1, nplane
        norm = norm + conjg(coeff(iplane,iband))*coeff(iplane,iband)
      end do
      coeff(:,iband) = coeff(:,iband) / sqrt(norm)
    end do
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Irreps of little group
    k_crys = wk
    if(iwk == 1) then
      write(klabel,'(a)') '1st'
    else if(iwk == 2) then
      write(klabel,'(a)') '2nd'
    else if(iwk == 3) then
      write(klabel,'(a)') '3rd'
    else if(iwk < 10) then
      write(klabel,'(i1,a)') iwk, 'th'
    else if(iwk < 100) then
      write(klabel,'(i2,a)') iwk, 'th'
    else if(iwk < 1000) then
      write(klabel,'(i3,a)') iwk, 'th'
    end if
    write (vasp_band,'(a,a,a,3f10.5)') 'This is the ', trim(adjustl(klabel)), ' k-point (cryst. coord.):', k_crys
    !!!!!!!!!!!!!!!!!!!!! vector representation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ordinary little group with all unitary operations
    call little_order(k_crys, lit_ord)
    allocate(lit_spa(lit_ord),lit_ele(3,4,lit_ord), lit_symb(lit_ord), lit_axis(3,lit_ord), toln(lit_ord))
    call little_group(vasp_band, k_crys, lit_ord, lit_spa, lit_ele, lit_symb, lit_axis)
    write (vasp_band, '(a/)') REPEAT('=', 100)
    allocate(lit_cha(lit_ord,lit_ord))
    call small_rep(vasp_band,k_crys,lit_ord,lit_ele,lit_symb,lit_axis,nirr,lit_cha,.false.)
    call Herring(vasp_band, k_crys, lit_ord, lit_spa, lit_ele, lit_symb, nirr, lit_cha, .false.)
    ! extended little group when consider the anti-unitary operation of TRS as en element
    call little_order_t(k_crys, lit_ord, lit_ord_t)
    allocate(lit_spa_t(lit_ord_t),lit_ele_t(3,4,lit_ord_t),lit_symb_t(lit_ord_t),lit_axis_t(3,lit_ord_t))
    call little_group_t(vasp_band, k_crys, lit_ord, lit_ord_t, lit_spa_t, lit_ele_t, lit_symb_t, lit_axis_t)
    allocate(lit_cha_t(lit_ord_t,lit_ord_t),lit_cha_t2(lit_ord_t,lit_ord_t))
    if(.NOT. soc) then
      call small_rep_t(vasp_band,k_crys,lit_ord,lit_ord_t,lit_spa_t,lit_ele_t,lit_symb_t,lit_axis_t,nirr_t,lit_cha_t,.false.)
    else
      call small_rep_t(vasp_band,k_crys,lit_ord,lit_ord_t,lit_spa_t,lit_ele_t,lit_symb_t,lit_axis_t,nirr_t2,lit_cha_t2,.false.)
    end if
    write (vasp_band, '(a/)') REPEAT('=', 100)
    !!!!!!!!!!!!!!!!!!!!!!! spinor representation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ordinary little group with all unitary operations
    call small_rep(vasp_band,k_crys,lit_ord,lit_ele,lit_symb,lit_axis,nirr,lit_cha,.true.)
    call Herring(vasp_band, k_crys, lit_ord, lit_spa, lit_ele, lit_symb, nirr, lit_cha, .true.)
    ! extended little group when consider the anti-unitary operation of TRS as en element
    if(.NOT. soc) then
      call small_rep_t(vasp_band,k_crys,lit_ord,lit_ord_t,lit_spa_t,lit_ele_t,lit_symb_t,lit_axis_t,nirr_t2,lit_cha_t2,.true.)
    else
      call small_rep_t(vasp_band,k_crys,lit_ord,lit_ord_t,lit_spa_t,lit_ele_t,lit_symb_t,lit_axis_t,nirr_t,lit_cha_t,.true.)
    end if
    write (vasp_band, '(a/)') REPEAT('=', 100)
    !
    ! calculate Irreps of VASP electronic bands
    allocate(band_cha(lit_ord,nband), cha(lit_ord))
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
      do j = 1, ncnt
        gvector2 = matmul(gvector(:,j),mat_inv(lit_ele(:,1:3,i)))
        do k = 1, ncnt
          if(all(abs((/ia,ib,ic/)+gvector2-gvector(:,k)) < tol3)) then
            do l = 1, nband
              if(.NOT. soc) then
                band_cha(i,l) = band_cha(i,l) + &
                conjg(coeff(k,l))*coeff(j,l)*exp(im*tpi*dot_product(k_crys+(/ia,ib,ic/)+gvector2,-lit_ele(:,4,i)))
              else ! with soc
                spinor1(1)=coeff(k,l); spinor1(2)=coeff(ncnt+k,l)
                spinor2(1)=coeff(j,l); spinor2(2)=coeff(ncnt+j,l)
                band_cha(i,l) = band_cha(i,l) + dot_product(spinor1, matmul(SU2_mat,spinor2))*&
                exp(im*tpi*dot_product(k_crys+(/ia,ib,ic/)+gvector2,-lit_ele(:,4,i)))
              end if
            end do ! l
            exit
          end if
        end do ! k
      end do ! j
    end do ! i
    ! print Irreps of VASP electronic bands
    write (vasp_band,'(a)') '   <<  Irreps of Electronic Bands  >>'
    if(.NOT. soc) write (vasp_band,'(a)') 'Determine the Irrep by comparing the band characters &
                                          &with the small Reps of the little group in "Scenario Without Considering SOC"'
    if(soc) write (vasp_band,'(a)') 'Determine the Irrep by comparing the band characters &
                                    &with the small Reps of the little group in "Scenario With Considering SOC"'
    write (vasp_band,'(a/)') 'For electronic Hamiltonian equation, the time-reversal symmetry is intrinsically holded'
    write (vasp_band,'(a)') 'Characters of electronic bands associated with each symmetry opeartion'
    do i = 1, lit_ord
      write (vasp_band,'(a5,a2,a,i2,a,a5,$)') ' ', lit_symb_t(i), '[',lit_spa_t(i),']', ' '
    end do
    write (vasp_band,*)
    p = 1; k = 0; toln = tol
    do while(.true.)
      cha = (0.d0,0.d0)
      do i = p, nband
        if(abs(cener(i)-cener(p)) < tol) then
          cha = cha + band_cha(:,i)
          if(i == nband) k = nband
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
          do iband = p, p+1
            do l  = 1, nirr_t
              if(all(abs(lit_cha_t(1:lit_ord,l)-band_cha(1:lit_ord,iband)) < toln(1:lit_ord))) exit
            end do
            write(vasp_band,'(a,$)') 'bands:'
            write(vasp_band,'(i4,$)') iband
            write (vasp_band,'(a5,a,a,i2,a,a5,a,2i4)') ' ', 'Irrep: ', '#', l, '*', ' ', 'accidental degeneracy:', p, p+1
            do j = 1, lit_ord
              write(vasp_band,'("(",f6.3,","f6.3,") ",$)') band_cha(j,iband)+(1.d-8,1.d-8)
            end do
            write (vasp_band, *)
          end do ! ibnd
        end if ! 2
        if(i-p == 3) then
          do l  = 1, nirr_t
            if(all(abs(lit_cha_t(1:lit_ord,l)-band_cha(1:lit_ord,p)) < toln(1:lit_ord))) exit
          end do
          if(l <= nirr_t) then
            write(vasp_band,'(a,$)') 'bands:'
            write(vasp_band,'(i4,$)') p
            write (vasp_band,'(a5,a,a,i2,a,a5,a,3i4)') ' ', 'Irrep: ', '#', l, '*', ' ', 'accidental degeneracy:', p, p+1, p+2
            do j = 1, lit_ord
              write(vasp_band,'("(",f6.3,","f6.3,") ",$)') band_cha(j,p)+(1.d-8,1.d-8)
            end do
            write (vasp_band, *)
            do l  = 1, nirr_t
              if(all(abs(lit_cha_t(1:lit_ord,l)-band_cha(1:lit_ord,p+1)) < toln(1:lit_ord))) exit
            end do
            if(l <= nirr_t) then
              write(vasp_band,'(a,$)') 'bands:'
              write(vasp_band,'(i4,$)') p+1
              write (vasp_band,'(a5,a,a,i2,a,a5,a,3i4)') ' ', 'Irrep: ', '#', l, '*', ' ', 'accidental degeneracy:', p, p+1, p+2
              do j = 1, lit_ord
                write(vasp_band,'("(",f6.3,","f6.3,") ",$)') band_cha(j,p+1)+(1.d-8,1.d-8)
              end do
              write (vasp_band, *)
              do l  = 1, nirr_t
                if(all(abs(lit_cha_t(1:lit_ord,l)-band_cha(1:lit_ord,p+2)) < toln(1:lit_ord))) exit
              end do
              write(vasp_band,'(a,$)') 'bands:'
              write(vasp_band,'(i4,$)') p+2
              write (vasp_band,'(a5,a,a,i2,a,a5,a,3i4)') ' ', 'Irrep: ', '#', l, '*', ' ', 'accidental degeneracy:', p, p+1, p+2
              do j = 1, lit_ord
                write(vasp_band,'("(",f6.3,","f6.3,") ",$)') band_cha(j,p+2)+(1.d-8,1.d-8)
              end do
              write (vasp_band, *)
            else
              do l  = 1, nirr_t
                if(all(abs(lit_cha_t(1:lit_ord,l)-band_cha(1:lit_ord,p+1)-band_cha(1:lit_ord,p+2)) < toln(1:lit_ord))) exit
              end do
              write(vasp_band,'(a,$)') 'bands:'
              write(vasp_band,'(2i4,$)') p+1, p+2
              write (vasp_band,'(a5,a,a,i2,a,a5,a,3i4)') ' ', 'Irrep: ', '#', l, '*', ' ', 'accidental degeneracy:', p, p+1, p+2
              do j = 1, lit_ord
                write(vasp_band,'("(",f6.3,","f6.3,") ",$)') band_cha(j,p+1)+band_cha(j,p+2)+(1.d-8,1.d-8)
              end do
              write (vasp_band, *)
            end if
          else
            do l  = 1, nirr_t
              if(all(abs(lit_cha_t(1:lit_ord,l)-band_cha(1:lit_ord,p)-band_cha(1:lit_ord,p+1)) < toln(1:lit_ord))) exit
            end do
            if(l <= nirr_t) then
              write(vasp_band,'(a,$)') 'bands:'
              write(vasp_band,'(2i4,$)') p, p+1
              write (vasp_band,'(a5,a,a,i2,a,a5,a,3i4)') ' ', 'Irrep: ', '#', l, '*', ' ', 'accidental degeneracy:', p, p+1, p+2
              do j = 1, lit_ord
                write(vasp_band,'("(",f6.3,","f6.3,") ",$)') band_cha(j,p)+band_cha(j,p+1)+(1.d-8,1.d-8)
              end do
              write (vasp_band, *)
              do l  = 1, nirr_t
                if(all(abs(lit_cha_t(1:lit_ord,l)-band_cha(1:lit_ord,p+2)) < toln(1:lit_ord))) exit
              end do
              write(vasp_band,'(a,$)') 'bands:'
              write(vasp_band,'(i4,$)') p+2
              write (vasp_band,'(a5,a,a,i2,a,a5,a,3i4)') ' ', 'Irrep: ', '#', l, '*', ' ', 'accidental degeneracy:', p, p+1, p+2
              do j = 1, lit_ord
                write(vasp_band,'("(",f6.3,","f6.3,") ",$)') band_cha(j,p+2)+(1.d-8,1.d-8)
              end do
              write (vasp_band, *)
            else
              do l  = 1, nirr_t
                if(all(abs(lit_cha_t(1:lit_ord,l)-band_cha(1:lit_ord,p)-band_cha(1:lit_ord,p+2)) < toln(1:lit_ord))) exit
              end do
              write(vasp_band,'(a,$)') 'bands:'
              write(vasp_band,'(2i4,$)') p, p+2
              write (vasp_band,'(a5,a,a,i2,a,a5,a,3i4)') ' ', 'Irrep: ', '#', l, '*', ' ', 'accidental degeneracy:', p, p+1, p+2
              do j = 1, lit_ord
                write(vasp_band,'("(",f6.3,","f6.3,") ",$)') band_cha(j,p)+band_cha(j,p+2)+(1.d-8,1.d-8)
              end do
              write (vasp_band, *)
              do l  = 1, nirr_t
                if(all(abs(lit_cha_t(1:lit_ord,l)-band_cha(1:lit_ord,p+1)) < toln(1:lit_ord))) exit
              end do
              write(vasp_band,'(a,$)') 'bands:'
              write(vasp_band,'(i4,$)') p+1
              write (vasp_band,'(a5,a,a,i2,a,a5,a,3i4)') ' ', 'Irrep: ', '#', l, '*', ' ', 'accidental degeneracy:', p, p+1, p+2
              do j = 1, lit_ord
                write(vasp_band,'("(",f6.3,","f6.3,") ",$)') band_cha(j,p+1)+(1.d-8,1.d-8)
              end do
              write (vasp_band, *)
            end if
          end if
        end if !3
        if(i-p > 3) then
          write(vasp_band,'(a,$)') 'bands:'
          do iband = p, i-1
            write(vasp_band,'(i4,$)') iband
          end do
          write (vasp_band,'(a5,a,a,i2,a,a5,a)') ' ', 'Irrep: ', '#', l, '*', ' ', 'accidental degeneracy:'
          do j = 1, lit_ord
            write(vasp_band,'("(",f6.3,","f6.3,") ",$)') cha(j)+(1.d-8,1.d-8)
          end do
          write (vasp_band, *)
        end if
      else
        write(vasp_band,'(a,$)') 'bands:'
        do iband = p, i-1
          write(vasp_band,'(i4,$)') iband
        end do
        write (vasp_band,'(a5,a,a,i2,a)') ' ', 'Irrep: ', '#', l, '*'
        do j = 1, lit_ord
          write(vasp_band,'("(",f6.3,","f6.3,") ",$)') cha(j)+(1.d-8,1.d-8)
        end do
        write (vasp_band, *)
      end if
      if(k == nband) exit
      p = i
    end do
    write (vasp_band, *)
    write (vasp_band, '(a/)') REPEAT('=', 100)
    !
    deallocate(lit_spa, lit_symb, lit_axis, lit_ele, lit_cha, &
               lit_spa_t, lit_symb_t, lit_axis_t, lit_ele_t, lit_cha_t, lit_cha_t2, gvector, band_cha, cha, toln)
    close(vasp_band)
  end do ! ik
  !
  deallocate(cener,occ,igall,coeff)
  !
  return
end subroutine vasp_bands
!
subroutine vcross(a,b,c)
  implicit none
  real(kind=8),intent(in) :: b(3),c(3)
  real(kind=8),intent(out) :: a(3)
  !
  a(1)=b(2)*c(3)-b(3)*c(2)
  a(2)=b(3)*c(1)-b(1)*c(3)
  a(3)=b(1)*c(2)-b(2)*c(1)
  return
end subroutine vcross