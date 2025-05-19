      program gmxkit 
      use xdr, only: xtcfile
      implicit real(a-h,o-z)
      type(xtcfile) :: xtc
      type(xtcfile) :: xtc1
      type(xtcfile) :: xtc2
      parameter (nskip=1000)
      parameter (delr=0.01e0)
      parameter (delg=0.5e0)
      parameter (zero=0.e0)
      parameter (chgmul=1.e6)
      parameter (ang2bohr=0.529e0)
      parameter (tiny=2.e-1)
      parameter (boltz =1.3806503e-23)  ! boltzmann constant in (J.K^-1)
      parameter (elechge=1.6022e-19)    ! unit electron charge in C
      parameter (temp=303.15)           ! K  ---needs to be corrected
      parameter (chgscale=0.7)          ! charge scaling factor, needs to be corrected
      logical :: valid,found_one,first,vehicular
      data first/.true./
      save first
      dimension cell(3,3),rcell(3,3),v1(3),v2(3),r(3),qori(3)
      dimension qcm1(3), layer(400),layer3(400)
      dimension umol(3,3)
      real*8 eleconne1,rand,urand
      real*8 eleconne2
      real*8, allocatable :: amsd(:,:),slope(:),dmsd(:,:)
      real*8, allocatable :: cmsd(:,:,:)
      integer*8, allocatable :: nh0(:),nht(:,:)
      logical, allocatable :: lgeomol(:)
      real, dimension(:), allocatable :: cnm,sigma,epsilon,chgsum,pm1,pm2,cdcm,cdss,chgmol,chg
      real, dimension(:), allocatable :: rkj,freesol,freean,area
      real, dimension(:,:), allocatable :: rdfmin,rdfminan,coord,coord_ion
      real, dimension(:,:), allocatable :: mass,charge,qcmol,gcm,gss,chgdisp,qdisp,nvhf1
      real, dimension(:,:), allocatable :: ctn,ct,tau
      real, dimension(:,:,:), allocatable :: qcm,qcmk,coordt,coord_iont,geomol
      real, dimension(:,:,:,:), allocatable :: qcmn,qqq
      integer, dimension(:), allocatable :: imole,nat,natan,ianion,nfreesol,nfreean,ncip,nsol1
      integer, dimension(:), allocatable :: nli1_sol,nli2_sol,nli3_sol,nli4_sol,nli1_an,nli2_an,nli3_an,nli4_an,nosol,noanion
      integer, dimension(:), allocatable :: res_num_i,n1,n2,inm,iat,itp_unit,npi,npj,itpi,itpj
      integer, dimension(:), allocatable :: natmn,nmoln0,nhtp,ngeom,natmd,igeomol
      integer, dimension(:,:), allocatable :: iatsol,iatan,nsol2,nrdf,indi,indj,nvhf
      integer, dimension(:,:), allocatable :: nda,aggdis
      integer, dimension(:,:,:), allocatable :: couple,couple_anion,ngrid
      integer, dimension(:,:,:,:), allocatable :: couple_aniono
      integer, dimension(:,:,:,:,:), allocatable :: couplert
      integer :: confout,atomtype,split_pos,cgnr,random_num
      character :: aA,residue,atom,bar
      character(len=100) :: line,molecule,ion,refrence
      character(len=100), dimension(:), allocatable :: res_num_c,res_name,res_type,itpfile
      character(len=100), dimension(:), allocatable :: atom_type,atom_type_itp,element
      character(len=10)  , dimension(:), allocatable :: mole,anion
      
      pi=4.e0*atan(1.e0)
      confout=5
      open(confout,file='confout.gro',status='old')
      read(confout,*)
      read(confout,*) natoms
      k=1
      l=1
      allocate (res_num_c(natoms),res_num_i(natoms),res_name(natoms))
      do i=1,natoms
        read(confout,*) line
        len_line=len_trim(line)
          do j=1,len_line
            if (line(j:j)>='A' .and. line(j:j) <= 'Z') then
              split_pos=j
              exit
            endif
          enddo
        res_num_c(k)=line(1:split_pos-1)
        res_name(l) =line(split_pos:len_line)
        read (res_num_c(k),*) res_num_i(k)
        k=k+1
        l=l+1
      enddo

      nty=1
      do i=1,natoms-1
        if (res_name(i) /= res_name(i+1)) then
          nty=nty+1
        endif
      enddo
      allocate (res_type(nty),n1(nty),n2(nty),inm(nty),iat(nty))
      res_type(1)=res_name(1)
      j=1
      do i=1,natoms-1
        if (res_name(i) /= res_name(i+1)) then
          n1(j)=i
          res_type(j+1)=res_name(i+1)
          j=j+1
        endif
      enddo
      n1(nty)=natoms
      n2(1)=n1(1)
      inm(1)=res_num_i(n1(1))
      iat(1)=n2(1)/inm(1)
      do i=2,nty
        n2(i)=n1(i)-n1(i-1)
        inm(i)=res_num_i(n1(i))-res_num_i(n1(i-1))
        iat(i)=n2(i)/inm(i)
      enddo
      write(*,'(a,i9)') 'Total number of atoms:', natoms
      write(*,'(a,i8)') 'Total type of residues:', nty
      do i=1,nty
        write(*,'(a8,i8,i8)') res_type(i),inm(i),iat(i)
      enddo

      nmols=0
      do i=1,nty
        nmols=nmols+inm(i)
      enddo
      
      nsize=inm(1)
      do i=2,nty
        if (inm(i) > nsize) then
          nsize=inm(i)
        endif
      enddo

      isize=iat(1)
      do i=2,nty
        if (iat(i) > isize) then
          isize=iat(i)
        endif
      enddo
      allocate (atom_type_itp(isize),charge(nty,isize),mass(nty,isize),chgsum(nty))
      allocate (itp_unit(nty),itpfile(nty))

      atomtype=7
      open(atomtype,file='atomtype.itp',status='old')

      nto=0
      do
        read(atomtype,'(A)',iostat=ios) line
        if (ios/=0) exit
        if (trim(line)=='') cycle
        if (line(1:1)==';') cycle
        if (line(1:1)=='[') cycle
        nto=nto+1
      enddo
      write(*,'(a,i9)') 'Total type of atoms',nto
      allocate (atom_type(nto),element(nto),sigma(nto),epsilon(nto))
      rewind(atomtype)
      i=1
      do
        read(atomtype,'(A)',iostat=ios) line
        if (ios/=0) exit
        if (trim(line)=='') cycle
        if (line(1:1)==';') cycle
        if (line(1:1)=='[') cycle
        backspace(atomtype)
        read(atomtype,*) atom_type(i),element(i),amass,acharge,aA,sigma(i),epsilon(i)
        write(*,'(a4,a4,2F10.5)') atom_type(i),element(i),sigma(i),epsilon(i)
        i=i+1
      enddo
      
      do i=1,nty
        itpfile(i)=trim(res_type(i)) // '.itp'
        write(*,'(a,a,a)') trim('-------------'),trim(res_type(i)),trim('-------------')
        itp_unit(i)=2*i
        open(itp_unit(i),file=itpfile(i),status='old')
        do
          read(itp_unit(i),'(A)',iostat=ios) line
          if (ios/=0) exit
          if (trim(line)=='') cycle
          if (line=='[ atoms ]') then
            exit
          endif
        enddo
        read(itp_unit(i),*)
        do j=1,iat(i)
          read(itp_unit(i),*) idx,atom_type_itp(j),resnr,residue,atom,cgnr,charge(i,j),mass(i,j)
          write(*,'(a8,2f10.4)') atom_type_itp(j),charge(i,j),mass(i,j)
        enddo
        chgsum(i)=sum(charge(i,1:iat(i)))
        write(*,'(a,f9.4)') 'SUMCHARGE',chgsum(i)
      enddo
      deallocate (res_num_c,res_num_i,res_name)
      deallocate (n1,n2)
      close(confout)
      close(atomtype)
      do i=1,nty
        close(itp_unit(i))
      enddo
      deallocate (itp_unit,itpfile)
      write(*,*)
      write(*,*) '========== Molecular Dynamics Analysis =========='
      write(*,*) ' 1) Cluster         2) RDF           3) MSD'
      write(*,*) ' 4) Residence Time  5) Self Van Hove 6) Distinct Van Hove'
      write(*,*) ' 7) SDF             8) XRD           9) Hopping'
      write(*,*) '10) Distance Angle 11) Viscosity    12) TI'
      write(*,*) ' 0) Quit'
      write(*,*) '------------>>'
      read(*,'(i)') inst
      if (inst==1) then
        do i=1,nty
          write(*,'(i3,a,a)')  i,') ',res_type(i)
        enddo
        write(*,*) '------------>>'
        open (101,file='cluster.inp',status='old')
        read (101,*) refrence
        do k=1,nty
          if (trim(refrence)==trim(res_type(k))) kty=k
        enddo
        write(*,*) refrence
        read (101,*) nkindsol
        write(*,*) nkindsol
        allocate (mole(nkindsol),nat(nkindsol),imole(nkindsol))
        read(101,*) (mole(i),i=1,nkindsol)
        do i=1,nkindsol
          do k=1,nty
            if (trim(mole(i))==trim(res_type(k))) then
              imole(i)=k
            endif
          enddo
        enddo
        write(*,*) (imole(i),i=1,nkindsol)
        read(101,*)  (nat(j),j=1,nkindsol)
        write(*,*)   (nat(j),j=1,nkindsol)
        natmax=nat(1)
        do j=2,nkindsol
          if (natmax < nat(j)) natmax=nat(j)
        enddo
        allocate (iatsol(nkindsol,natmax),rdfmin(nkindsol,natmax))
        read (101,*) ((iatsol(i,j),j=1,nat(i)),i=1,nkindsol)
        write (*,*)  ((iatsol(i,j),j=1,nat(i)),i=1,nkindsol)
        read (101,*) ((rdfmin(i,j),j=1,nat(i)),i=1,nkindsol)
        write (*,*)  ((rdfmin(i,j),j=1,nat(i)),i=1,nkindsol)
        
        read (101,*)
        read (101,*) nkindan
        write(*,*) nkindan
        allocate (anion(nkindan),natan(nkindan),ianion(nkindan))
        read(101,*) (anion(i),i=1,nkindan)
        do i=1,nkindan
          do k=1,nty
            if (trim(anion(i))==trim(res_type(k))) then
              ianion(i)=k
            endif
          enddo
        enddo
        write(*,*)  (ianion(i),i=1,nkindan)
        read(101,*)  (natan(j),j=1,nkindan)
        write(*,*)   (natan(j),j=1,nkindan)
        natanmax=natan(1)
        do j=2,nkindan
          if (natanmax < natan(j)) natanmax=natan(j)
        enddo
        allocate (iatan(nkindan,natanmax),rdfminan(nkindan,natanmax))
        read (101,*) ((iatan(i,j),j=1,natan(i)),i=1,nkindan)
        write (*,*)  ((iatan(i,j),j=1,natan(i)),i=1,nkindan)
        read (101,*) ((rdfminan(i,j),j=1,natan(i)),i=1,nkindan)
        write (*,*)  ((rdfminan(i,j),j=1,natan(i)),i=1,nkindan)
        open (102,file='freemol.dat')
        open (103,file='cluster.dat')
        open (104,file='ssip.dat')
        if (nkindsol==1 .and. nkindan==1) then
          open (105,file='agg1.dat')
          write(105,'(7a16)') 'Frame','Li Number','Nlayer','NLi','NAnion','NLi(Total)','NAnion(Total)'
          open (106,file='agg2.dat')
          write(106,'(5a16)') 'Frame','Agg Number','NLi(Total)','NAnion(Total)','N(Total)'
          open (107,file='agg3.dat')
          open (108,file='agg4.dat')
          open (111,file='Li_nOsol')
          open (112,file='Li_nOanion')
        endif
        open (109,file='sol.dat')
        open (110,file='anion.dat')
        
        allocate (coord(natoms,3),coord_ion(inm(kty),3))
        allocate (nfreesol(nkindsol),nfreean(nkindan),freesol(nkindsol),freean(nkindan))
        allocate (nli1_sol(nkindsol),nli2_sol(nkindsol),nli3_sol(nkindsol),nli4_sol(nkindsol))
        allocate (nli1_an(nkindan),nli2_an(nkindan),nli3_an(nkindan),nli4_an(nkindan))
        allocate (nosol(0:6),noanion(0:6))
        allocate (ncip(nkindan),nsol1(6),nsol2(0:6,0:6))
        allocate (aggdis(0:inm(kty),0:inm(ianion(1))))
        coord(:,:)=0.e0
        coord_ion(:,:)=0.e0
        nfreesol(:)=0
        nfreean(:)=0
        freesol(:)=0.e0
        freean(:)=0.e0
        nli1_sol(:)=0
        nli2_sol(:)=0
        nli3_sol(:)=0
        nli4_sol(:)=0
        nli1_an(:)=0
        nli2_an(:)=0
        nli3_an(:)=0
        nli4_an(:)=0
        nosol(:)=0
        noanion(:)=0
        nssip=0
        ncip(:)=0
        nsol1(:)=0
        nsol2(:,:)=0
        aggdis(:,:)=0
        kmax=0
        jmax=0
        imax=0
        nlen=0
        call xtc % init('traj_comp.xtc')
        do while (xtc % STAT==0)
          call xtc % read
          cell=xtc % box
          cell=cell*10.e0
          call recip(cell,rcell,vol)
          nlen=nlen+1
          do i=1,natoms
            coord(i,:)=xtc % pos(:,i)*10.e0
          enddo
          natmk0=0
          do k=1,kty-1
            natmk0=natmk0+inm(k)*iat(k)
          enddo
          do k=1,inm(kty)
            coord_ion(k,:)=coord(natmk0+1+(k-1)*iat(kty),:)
          enddo

          allocate (couple(nkindsol,inm(kty),nsize),natmn(20),rkj(natmax))
          allocate (couple_anion(nkindan,inm(kty),nsize),couple_aniono(nkindan,inm(kty),nsize,natan(1)))
          couple(:,:,:)=0
          couple_anion(:,:,:)=0
          couple_aniono(:,:,:,:)=0
          do n=1,nkindsol
            natmn0=0
            do j=1,imole(n)-1
              natmn0=natmn0+inm(j)*iat(j)
            enddo
            do k=1,inm(kty)
              do j=1,inm(imole(n))
                do m=1,nat(n)
                  natmn(m)=natmn0+iatsol(n,m)+(j-1)*iat(imole(n))
                  rkj(m)=dist(coord_ion(k,1:3),coord(natmn(m),1:3),cell,rcell)
                  if (rkj(m) <= rdfmin(n,m)) then
                    couple(n,k,j)=1
                    exit
                  endif
                enddo
              enddo
            enddo
          enddo
          
          do n=1,nkindan
            natmn0=0
            do j=1,ianion(n)-1
              natmn0=natmn0+inm(j)*iat(j)
            enddo
            do k=1,inm(kty)
              do j=1,inm(ianion(n))
                do m=1,natan(n)
                  natmn(m)=natmn0+iatan(n,m)+(j-1)*iat(ianion(n))
                  rkj(m)=dist(coord_ion(k,1:3),coord(natmn(m),1:3),cell,rcell)
                  if (rkj(m) <= rdfminan(n,m)) then
                    couple_anion(n,k,j)=1
                    exit
                  endif
                enddo
                do m=1,natan(n)
                  natmn(m)=natmn0+iatan(n,m)+(j-1)*iat(ianion(n))
                  rkj(m)=dist(coord_ion(k,1:3),coord(natmn(m),1:3),cell,rcell)
                  if (rkj(m) <= rdfminan(n,m)) then
                    couple_aniono(n,k,j,m)=1
                  endif
                enddo
              enddo
            enddo
          enddo
          ! Freemol
          do n=1,nkindsol
            do i=1,inm(imole(n))
              if (sum(couple(n,:,i))==0) then
                nfreesol(n)=nfreesol(n)+1
              endif
              if (sum(couple(n,:,i))==1) nli1_sol(n)=nli1_sol(n)+1
              if (sum(couple(n,:,i))==2) nli2_sol(n)=nli2_sol(n)+1
              if (sum(couple(n,:,i))==3) nli3_sol(n)=nli3_sol(n)+1
              if (sum(couple(n,:,i))==4) nli4_sol(n)=nli4_sol(n)+1
            enddo
            freesol(n)=real(nfreesol(n))/real(nlen)
          enddo
          do n=1,nkindan
            do i=1,inm(ianion(n))
              if (sum(couple_anion(n,:,i))==0) then
                nfreean(n)=nfreean(n)+1
              endif
              if (sum(couple_anion(n,:,i))==1) nli1_an(n)=nli1_an(n)+1
              if (sum(couple_anion(n,:,i))==2) nli2_an(n)=nli2_an(n)+1
              if (sum(couple_anion(n,:,i))==3) nli3_an(n)=nli3_an(n)+1
              if (sum(couple_anion(n,:,i))==4) nli4_an(n)=nli4_an(n)+1
            enddo
            freean(n)=real(nfreean(n))/real(nlen)
          enddo

          if (mod(nlen,nskip)==0) then
            write(102,'(i10,20f20.4)') nlen,(freesol(n),n=1,nkindsol),(freean(m),m=1,nkindan)
            write(109,'(i10,20f20.4)') nlen,((real(nli1_sol(n))/real(nlen),real(nli2_sol(n))/real(nlen),real(nli3_sol(n))/real(nlen),real(nli4_sol(n))/real(nlen)),n=1,nkindsol)
            write(110,'(i10,20f20.4)') nlen,((real(nli1_an(n))/real(nlen),real(nli2_an(n))/real(nlen),real(nli3_an(n))/real(nlen),real(nli4_an(n))/real(nlen)),n=1,nkindan)
          endif
          
          if (nkindsol==1 .and. nkindan==1) then
            do k=1,inm(kty)
              l=sum(couple(1,k,:))
              nosol(l)=nosol(l)+1
              m=sum(couple_aniono(1,k,:,:))
              noanion(m)=noanion(m)+1
            enddo
          endif

          if (mod(nlen,nskip)==0) then
            write(111,'(i10,20f20.4)') nlen,(real(nosol(n))/real(nlen),n=0,6)
            write(112,'(i10,20f20.4)') nlen,(real(noanion(n))/real(nlen),n=0,6)
          endif
             
          ! SSIP
          nnssip=0
          do k=1,inm(kty)
            nanion=sum(couple_anion(:,k,:))
            if (nanion==0) then
              nssip=nssip+1
              nnssip=nnssip+1
              if (nkindsol==1) then
                l=sum(couple(1,k,:))
                nsol1(l)=nsol1(l)+1
              endif
              if (nkindsol==2) then
                l=sum(couple(1,k,:))
                m=sum(couple(2,k,:))
                nsol2(l,m)=nsol2(l,m)+1
              endif
            endif
          enddo
          ! CIP
          do k=1,inm(kty)
            if (sum(couple_anion(:,k,:))==1) then
              valid = .true.
              do n=1,nkindan
                do j=1,inm(ianion(n))
                  if (couple_anion(n,k,j)==1 .and. sum(couple_anion(n,:,j))/=1) then
                    valid = .false.
                    exit
                  endif
                enddo
                if (valid) ncip(n)=ncip(n)+1
              enddo
            endif
          enddo 

          if (mod(nlen,nskip)==0) then
            write(103,'(i10,20f20.4)')  nlen,real(nssip)/real(nlen),(real(ncip(n))/real(nlen),n=1,nkindan)
            if (nkindsol==1) then
              write(104,'(i10,20f20.4)') nlen, (real(nsol1(n))/real(nlen),n=1,6)
            endif
            if (nkindsol==2) then
              write(104,'(i10)') nlen
              write(104,'(20f16.4)') (real(nsol2(0,n))/real(nlen),n=0,6)
              write(104,'(20f16.4)') (real(nsol2(1,n))/real(nlen),n=0,6)
              write(104,'(20f16.4)') (real(nsol2(2,n))/real(nlen),n=0,6)
              write(104,'(20f16.4)') (real(nsol2(3,n))/real(nlen),n=0,6)
              write(104,'(20f16.4)') (real(nsol2(4,n))/real(nlen),n=0,6)
              write(104,'(20f16.4)') (real(nsol2(5,n))/real(nlen),n=0,6)
              write(104,'(20f16.4)') (real(nsol2(6,n))/real(nlen),n=0,6)
              write(104,*)
            endif
          endif
          
          ! Cluster
          if (nkindan==1) then
            iagg=0
            kmaxlen=0
            jmaxlen=0
            do k=1,inm(kty)
              nlayer=1
              nli_agg=1
              nli_aggTL=1
              nanion_agg=0
              nanion_aggTL=0
              n_aggTL=0
              if (sum(couple_anion(1,k,:)) > 0) then
                if (mod(nlen,nskip)==0) then
                  write(105,'(9i16)') nlen,k,nlayer,nli_agg,nanion_agg,nli_aggTL,nanion_aggTL
                endif
                nli_agg=0
                layer(:)=0
                do j=1,inm(ianion(1))
                  if (couple_anion(1,k,j)==1) then
                    nanion_agg=nanion_agg+1
                    layer(nanion_agg)=j
                  endif
                enddo

                nanion_aggTL=nanion_aggTL+nanion_agg
                couple_anion(1,k,:)=0
                nlayer=nlayer+1
                if (mod(nlen,nskip)==0) then
                  write(105,'(9i16)') nlen,k,nlayer,nli_agg,nanion_agg,nli_aggTL,nanion_aggTL
                endif 

                do
                  layer3(:)=0
                  do l=1,nanion_agg
                    do m=1,inm(kty)
                      if (couple_anion(1,m,layer(l))==1) then
                        nli_agg=nli_agg+1
                        layer3(nli_agg)=m
                      endif
                    enddo
                    couple_anion(1,:,layer(l))=0
                  enddo

                  call remove_duplicates(layer3,nli_agg)
                  layer(:)=0
                  nli_aggTL=nli_aggTL+nli_agg
                  if (nli_agg==0) exit
                  if (nli_agg >0) then
                    nlayer=nlayer+1
                    nanion_agg=0
                    if (mod(nlen,nskip)==0) then
                      write(105,'(9i16)') nlen,k,nlayer,nli_agg,nanion_agg,nli_aggTL,nanion_aggTL
                    endif
                    do n=1,nli_agg
                      do j3=1,inm(ianion(1))
                        if (couple_anion(1,layer3(n),j3)==1) then
                          nanion_agg=nanion_agg+1
                          layer(nanion_agg)=j3
                        endif
                      enddo
                      couple_anion(1,layer3(n),:)=0
                    enddo
                    call remove_duplicates(layer,nanion_agg)
                    nlayer=nlayer+1
                    nli_agg=0
                    nanion_aggTL=nanion_aggTL+nanion_agg
                    if (mod(nlen,nskip)==0) then
                      write(105,'(9i16)') nlen,k,nlayer,nli_agg,nanion_agg,nli_aggTL,nanion_aggTL
                    endif
                  endif
                enddo ! do

                iagg=iagg+1
                n_aggTL=nli_aggTL+nanion_aggTL
                if (mod(nlen,nskip)==0) then
                  write(106,'(6i16)') nlen,iagg,nli_aggTL,nanion_aggTL,n_aggTL
                endif
                if(nli_aggTL > kmax)    kmax=nli_aggTL
                if(nanion_aggTL > jmax) jmax=nanion_aggTL
                if(nli_aggTL > kmaxlen)    kmaxlen=nli_aggTL
                if(nanion_aggTL > jmaxlen) jmaxlen=nanion_aggTL
                aggdis(nli_aggTL,nanion_aggTL)=aggdis(nli_aggTL,nanion_aggTL)+1
              endif
            enddo
            aggdis(0,1)=nfreean(1)
            aggdis(1,0)=nssip
            sum_kmaxlen=sum_kmaxlen+kmaxlen
            sum_jmaxlen=sum_jmaxlen+jmaxlen
            sum_cluster=sum_cluster+nnssip+iagg
            if (mod(nlen,nskip)==0) then
              write(108,'(i10,5f20.4)') nlen,sum_cluster/real(nlen),sum_kmaxlen/real(nlen),sum_jmaxlen/real(nlen)
            endif
            nsum=0
            do k=0,kmax
              do j=0,jmax
                nsum=nsum+aggdis(k,j)
              enddo
            enddo
            msum=0
            do k=0,kmax
              do j=0,jmax
                msum=msum+aggdis(k,j)*k+aggdis(k,j)*j
              enddo
            enddo
            if (mod(nlen,nskip*50)==0) then
              write(107,*) 'nlen=',nlen
              write(107,*) 'kmax=',kmax
              write(107,*) 'jmax=',jmax
              do k=0,kmax
                write(107,'(500f10.6)') (real(aggdis(k,j))/real(nsum), j=0,jmax)
              enddo
              write(107,*)
              write(107,*)
              do k=0,kmax
                write(107,'(500f10.6)') (real(aggdis(k,j)*k+aggdis(k,j)*j)/real(msum), j=0,jmax)
              enddo
              write(107,*)
              write(107,*)
            endif
          endif
          deallocate (couple,couple_anion,couple_aniono,natmn,rkj)
        enddo ! while
        call xtc % close
        deallocate (mole,nat,imole,iatsol,rdfmin)
        deallocate (anion,natan,ianion,iatan,rdfminan)
        deallocate (coord,coord_ion)
        deallocate (nfreesol,nfreean,freesol,freean,ncip)
        deallocate (nli1_sol,nli2_sol,nli3_sol,nli4_sol)
        deallocate (nli1_an,nli2_an,nli3_an,nli4_an)
        deallocate (nosol,noanion)
        deallocate (nsol1,nsol2)
        deallocate (aggdis)
        close (101)
        close (102)
        close (103)
        close (104)
        close (105)
        close (106)
        close (107)
        close (108)
        close (109)
        close (110)
      endif !inst=1
      
      if (inst==2) then
        ! cm-cm
        open(201,file='rdf.inp',status='old')
        read(201,*) ncomrdf
        if (ncomrdf==0) ncom=0
        if (ncomrdf>0)  ncom=1
        do icom=1,ncom
          call xtc % init("traj_comp.xtc")
          call xtc % read
          cell=xtc % box
          cell=cell*10.e0
          call recip(cell,rcell,vol)
          call xtc1 % init("cmcoor1.xtc", 'w')
          allocate (qcmol(3,nmols))
          nlens=0
          do while (xtc % STAT == 0)
            iatm=0
            imol=0
            do i=1,nty
              allocate (coord(iat(i),3))
              do j=1,inm(i)
                imol=imol+1
                do k=1,iat(i)
                  iatm=iatm+1
                  coord(k,:)=xtc % pos(:,iatm)*10.e0
                enddo
                call cenmas(iat(i),coord,mass(i,1:iat(i)),qcm1,cell,rcell)
                qcmol(1:3,imol)=qcm1(1:3)
              enddo
              deallocate (coord)
            enddo
            call xtc1 % write(nmols,nlens,xtc % time,xtc % box*10.e0,qcmol,xtc % prec)
            nlens=nlens+1
            call xtc % read
          enddo
          write(*,*)'nlen ',nlen
          call xtc % close
          call xtc1 % close
          deallocate (qcmol)
        
          allocate (qcm(nty,nsize,3))
          celllen1=sqrt(cell(1,1)**2+cell(2,1)**2+cell(3,1)**2)
          celllen2=sqrt(cell(1,2)**2+cell(2,2)**2+cell(3,2)**2)
          celllen3=sqrt(cell(1,3)**2+cell(2,3)**2+cell(3,3)**2)
          celllen=celllen1
          if(celllen > celllen2) celllen=celllen2
          if(celllen > celllen3) celllen=celllen3
          nlenrdf=nint(celllen/delr/2.e0)
          
          allocate (cnm(ncomrdf))
          allocate (nrdf(ncomrdf,0:nlenrdf))
          allocate (indi(ncomrdf,1),indj(ncomrdf,1))
          allocate (gcm(nlenrdf,ncomrdf),pm1(ncomrdf),cdcm(ncomrdf))
          do i=1,ncomrdf
            read(201,*) indi(i,1),indj(i,1)
          enddo
          do i=1,ncomrdf
            do j=0,nlenrdf
              nrdf(i,j)=0.e0
            enddo
          enddo
          nlen=0
          call xtc1 % init("cmcoor1.xtc")
          do while (xtc1 % STAT==0)
            nlen=nlen+1
            call xtc1 % read
            imol=0
            do itp=1,nty
              do jnm=1,inm(itp)
                imol=imol+1
                qcm(itp,jnm,1:3)=xtc1 % pos(1:3,imol)
              enddo
            enddo
            do i=1,ncomrdf
              imo=indi(i,1)
              jmo=indj(i,1)
              if (imo == jmo) then
                do j=1,inm(imo)-1
                  do k=j+1,inm(imo)
                    rjk=dist(qcm(imo,j,1:3),qcm(imo,k,1:3),cell,rcell)
                    nbin=nint(rjk/delr)
                    if(nbin <= nlenrdf) then
                      nrdf(i,nbin)=nrdf(i,nbin)+2
                    endif
                  enddo
                enddo
              else
                do j=1,inm(imo)
                  do k=1,inm(jmo)
                    rjk=dist(qcm(imo,j,1:3),qcm(jmo,k,1:3),cell,rcell)
                    nbin=nint(rjk/delr)
                    if(nbin <= nlenrdf) then
                      nrdf(i,nbin)=nrdf(i,nbin)+1
                    endif
                  enddo
                enddo
              endif
              rho=real(inm(jmo)/vol)
              cnm(i)=4.e0/3.e0*pi*rho*real(inm(imo))
            enddo
            if (mod(nlen,1000)==0) then
              write (*,'(a10,f6.2,a)')  'Processs:', real(nlen*100)/real(nlens),'%'
              open(202,file='rdfcm.dat')
              r0=0.5e0*delr
              r13=(0.5e0*delr)**3
              do j=1,nlenrdf
                r0=r0+delr
                r23=r0**3
                dr3=r23-r13
                r13=r23
                write(202,'(f10.4,12f16.4)') r0,(real(nrdf(i,j)/cnm(i)/dr3/real(nlen)),i=1,ncomrdf)
              enddo
              close(202)
              open(202,file='rdfcm.dat',status='old')
              do i=1,nlenrdf
                read(202,*) dis,(gcm(i,j),j=1,ncomrdf)
              enddo
              para=4.e0*pi/vol
              do j=1,ncomrdf
                cdcm(j)=0.e0
                pm1(j)=para*real(inm(indj(j,1)))
              enddo
              open(203,file='cdfcm.dat')
              do i=1,nlenrdf
                do j=1,ncomrdf
                  cdcm(j)=gcm(i,j)*(real(i)*delr)**2*delr+cdcm(j)
                enddo
                write(203,'(13f12.6)') real(i)*delr,(cdcm(j)*pm1(j),j=1,ncomrdf)
              enddo
              close(202)
              close(203)
            endif
          enddo
          deallocate (qcm,cnm)
          deallocate (nrdf,indi,indj)
          deallocate (gcm,pm1,cdcm)
          call xtc1 % close
          call system('/bin/rm -f cmcoor1.xtc')
        enddo ! icom

        !site-site
        read(201,*) nssrdf
        if (nssrdf==0) nss=0
        if (nssrdf> 0) nss=1
        do iss=1,nss
          call xtc % init('traj_comp.xtc')
          call xtc % read
          nlens=0
          do while (xtc % STAT==0)
            nlens=nlens+1
            call xtc %read
          enddo
          cell=xtc % box
          cell=cell*10.e0
          call recip(cell,rcell,vol)
          call xtc % close
          celllen1=sqrt(cell(1,1)**2+cell(2,1)**2+cell(3,1)**2)
          celllen2=sqrt(cell(1,2)**2+cell(2,2)**2+cell(3,2)**2)
          celllen3=sqrt(cell(1,3)**2+cell(2,3)**2+cell(3,3)**2)
          celllen=celllen1
          if(celllen > celllen2) celllen=celllen2
          if(celllen > celllen3) celllen=celllen3
          nlenrdf=nint(celllen/delr/2.e0)
          allocate (nrdf(nssrdf,0:nlenrdf))
          allocate (npi(nssrdf),npj(nssrdf))
          allocate (indi(nssrdf,isize),indj(nssrdf,isize))
          allocate (itpi(nssrdf),itpj(nssrdf))
          allocate (cnm(nssrdf))
          allocate (gss(nlenrdf,nssrdf),pm2(nssrdf),cdss(nssrdf))
          do issrdf=1,nssrdf
            read(201,*) itpi(issrdf),npi(issrdf)
            read(201,*) (indi(issrdf,k),k=1,npi(issrdf))
            read(201,*) itpj(issrdf),npj(issrdf)
            read(201,*) (indj(issrdf,k),k=1,npj(issrdf))
          enddo
          close (201)
          do issrdf=1,nssrdf
            do k=0,nlenrdf
              nrdf(issrdf,k)=0
            enddo
          enddo
          allocate (coord(natoms,3))
          nlen=0
          call xtc % init("traj_comp.xtc")
          do while (xtc % STAT==0)
            call xtc % read
            nlen=nlen+1
            do i=1,natoms
              coord(i,:)=xtc % pos(:,i)*10.e0
            enddo
            do issrdf=1,nssrdf
              itp=itpi(issrdf)
              jtp=itpj(issrdf)
              iatm0=0
              do i=1,itp-1
                iatm0=iatm0+inm(i)*iat(i)
              enddo
              jatm0=0
              do j=1,jtp-1
                jatm0=jatm0+inm(j)*iat(j)
              enddo
              do ipi=1,npi(issrdf)
                do irep=0,inm(itp)-1
                  iatm=iatm0+irep*iat(itp)+indi(issrdf,ipi)
                  do jpj=1,npj(issrdf)
                    do jrep=0,inm(jtp)-1
                      jatm=jatm0+jrep*iat(jtp)+indj(issrdf,jpj)
                      rjk=dist(coord(iatm,1:3),coord(jatm,1:3),cell,rcell)
                      nbin=nint(rjk/delr)
                      if(nbin <= nlenrdf) then
                        nrdf(issrdf,nbin)=nrdf(issrdf,nbin)+1
                      endif
                    enddo
                  enddo
                enddo
              enddo
              rho=real(inm(jtp)*npj(issrdf))/vol
              cnm(issrdf)=4.e0/3.e0*pi*rho*real(inm(itp)*npi(issrdf))
            enddo
            if (mod(nlen,1000)==0) then
              write (*,'(a10,f6.2,a)')  'Processs:', real(nlen*100)/real(nlens),'%'
              open(204,file='rdfss.dat')
              r0=0.5e0*delr
              r13=(0.5e0*delr)*3
              do j=1,nlenrdf
                r0=r0+delr
                r23=r0**3
                dr3=r23-r13
                r13=r23
                write(204,'(13f12.6)') r0,(real(nrdf(i,j))/cnm(i)/dr3/real(nlen),i=1,nssrdf)
              enddo
              close(204)
              open(204,file='rdfss.dat',status='old')
              do i=1,nlenrdf
                read(204,*) dis,(gss(i,j),j=1,nssrdf)
              enddo
              para=4.e0*pi/vol
              do j=1,nssrdf
                cdss(j)=0.e0
                pm2(j)=para*real(inm(itpj(j))*npj(j))
              enddo
              open(205,file='cdfss.dat')
              do i=1,nlenrdf
                do j=1,nssrdf
                  cdss(j)=gss(i,j)*(real(i)*delr)**2*delr+cdss(j)
                enddo
                write(205,'(13f12.6)') real(i)*delr,(cdss(j)*pm2(j),j=1,nssrdf)
              enddo
              close(204)
              close(205)
            endif
          enddo
          deallocate (coord,cnm,nrdf)
          deallocate (indi,indj,itpi,itpj,npi,npj)
          deallocate (gss,pm2,cdss)
          call xtc % close
        enddo ! iss

      endif ! inst=2

      if (inst==3) then
        write(*,*) 'enter the correlated time (ps)'
        read(*,*) tmsd
        allocate (chgmol(nty))
        lcond=.false.
        do i=1,nty
          chgmol(i)=0.
          charge(i,1:iat(i))=charge(i,1:iat(i))/chgscale
          chgmol(i)=real(nint(sum(charge(i,1:iat(i)))))
          if (abs(chgmol(i)) > tiny) lcond=.true.
        enddo
        allocate (qcmol(3,nmols))
        call xtc % init('traj_comp.xtc')
        call xtc % read
        cell=xtc % box
        cell=cell*10.0
        call recip(cell,rcell,vol)
        elecon=1.e20*(elechge**2)/3./boltz/temp/vol
        call xtc1 % init("cmcoor1.xtc", 'w')
        nlen=0
      !  open(666,file='current.dat')
        open(301,file='time_discontinued.dat')
        do while ( xtc % STAT == 0 )
          t1=xtc % time
          nlen=nlen+1
          if (nlen == 2) dt=t1
          iatm=0
          imol=0
          do i=1,nty
            allocate (coord(iat(i),3))
            do j=1,inm(i)
              imol=imol+1
              do k=1,iat(i)
                iatm=iatm+1
                coord(k,:)=xtc % pos(:,iatm)*10.0
              enddo
              call cenmas(iat(i),coord,mass(i,1:iat(i)),qcm1,cell,rcell)
              qcmol(1:3,imol)=qcm1(1:3)
            enddo
            deallocate (coord)
          enddo
          call xtc1 % write(nmols,nlen,xtc % time,xtc % box*10.,qcmol,xtc % prec)
          call xtc % read
        enddo
        write(301,'(f20.7,3h ps)') dt
        write(*,*)'nlen ',nlen
        call xtc % close
        call xtc1 % close
        deallocate (qcmol)
        close(301)

        call xtc1 % init('cmcoor1.xtc')
        call xtc2 % init('cmdisp.xtc','w')
        if (lcond) then
          ichgdisp=302
          open(ichgdisp,file='charge_disp')
        endif
        if(lcond) allocate (chg(nmols))
        imol=0
        do i=1,nty
          do j=1,inm(i)
            imol=imol+1
            if(lcond) chg(imol)=chgmol(i)
          enddo
        enddo
        allocate (qcm(2,nmols,3),qdisp(3,nmols))
        allocate (chgdisp(nty,3))
        qdisp(:,:)=0.
        call xtc1 % read
        do i=1,nmols
          qcm(1,i,1:3)=xtc1 % pos(1:3,i)
        enddo
        do ilen=1,nlen
          do i=1,nmols
            qcm(2,i,1:3)=xtc1 % pos(1:3,i)
          enddo
          call pbc2(nmols,qcm,qdisp,cell,rcell)
          if (lcond) chgdisp(:,:)=0.
          imol=0
          do i=1,nty
            do j=1,inm(i)
              imol=imol+1
              do k=1,3
                qcm(1,imol,k)=qcm(2,imol,k)
              enddo
              if (lcond) then
                if (abs(chg(imol)) > tiny) then
                  do k=1,3
                    chgdisp(i,k)=chgdisp(i,k)+chg(imol)*qdisp(k,imol)
                  enddo
                endif
              endif
            enddo
            if(lcond) write(ichgdisp,'(f12.4,3g20.10)') real(ilen-1)*dt,(chgdisp(i,k),k=1,3)
          enddo
          call xtc2 % write(nmols,ilen-1,xtc1 % time,xtc1 % box,qdisp,xtc1 % prec)
          call xtc1 % read
        enddo
        deallocate (qdisp,qcm,chgdisp)
        if(lcond) deallocate (chg)
        call xtc1 % close
        call xtc2 % close
        if (lcond) close (ichgdisp)
        call system('/bin/rm -f cmcoor1.xtc')
        
        imsdcond=48
        open(imsdcond,file='msd_cond.dat')
        nmsd=nint(tmsd/dt)
        if(nmsd*2.gt.nlen) then
          write(*,*)'tmsd too long or MD run too short!'
          stop
        endif
        nread=int(nlen/nmsd)
        write(imsdcond,*)'nlen =',nlen
        write(imsdcond,*)'tmsd =',tmsd
        write(imsdcond,*)'dt =',dt
        write(imsdcond,*)'nmsd =',nmsd
        write(imsdcond,*)

!       calculate conductivity from total charge displacement chgdisp
        if (lcond) then
          allocate (qcm(nmsd*2,nty+1,3))
          allocate (cmsd(nty,nty,0:nmsd))
          cmsd(:,:,:)=0.d0
          open(ichgdisp,file='charge_disp',status='old')
          rewind(ichgdisp)
          do i=1,nmsd
            do ity=1,nty
              read(ichgdisp,*)ttt,(qcm(i,ity,j),j=1,3)
            enddo
          enddo
          do iread=2,nread
            do i=nmsd+1,nmsd*2
              do ity=1,nty
                read(ichgdisp,*)ttt,(qcm(i,ity,j),j=1,3)
              enddo
            enddo
            do i=1,nmsd
              do j=i,i+nmsd
                do ity=1,nty
                  do itp=ity,nty
                    dum=(qcm(j,ity,1)-qcm(i,ity,1))*(qcm(j,itp,1)-qcm(i,itp,1)) &
                       +(qcm(j,ity,2)-qcm(i,ity,2))*(qcm(j,itp,2)-qcm(i,itp,2)) &
                       +(qcm(j,ity,3)-qcm(i,ity,3))*(qcm(j,itp,3)-qcm(i,itp,3))
                    cmsd(ity,itp,j-i)=cmsd(ity,itp,j-i)+dble(dum)
                  enddo
                enddo
              enddo
            enddo
            do i=1,nmsd
              do j=1,nty
                do k=1,3
                  qcm(i,j,k)=qcm(i+nmsd,j,k)
                enddo
              enddo
            enddo
            call runave(nty,iread,nmsd,cmsd,dt,elecon,tmsd,imsdcond,res_type)            
          enddo
          call flush(imsdcond)
          deallocate (qcm,cmsd)
          close (ichgdisp)
          call system('/bin/rm -f charge_disp')
        endif

        ! calculate msd from center of mass displacement cmdisp
        write(imsdcond,'(/,a)')'start self MSD:'
        allocate (qcm(nmsd*2,nmols,3))
        allocate (amsd(nty,0:nmsd))
        amsd(:,:)=0.d0
        call xtc2 % init("cmdisp.xtc")
        do i=1,nmsd
          call xtc2 % read
          do imol=1,nmols
            qcm(i,imol,1:3)=xtc2 % pos(1:3,imol)
          enddo
        enddo
        do iread=2,nread
          do i=nmsd+1,nmsd*2 
            call xtc2 % read
            do imol=1,nmols
              qcm(i,imol,1:3)=xtc2 % pos(1:3,imol)
            enddo
          enddo
          imol=0
          do ity=1,nty
            do iml=1,inm(ity)
              imol=imol+1
              do i=1,nmsd
                do j=i,i+nmsd
                  dum=(qcm(j,imol,1)-qcm(i,imol,1))**2 &
                     +(qcm(j,imol,2)-qcm(i,imol,2))**2 &
                     +(qcm(j,imol,3)-qcm(i,imol,3))**2
                  amsd(ity,j-i)=amsd(ity,j-i)+dble(dum)
                enddo
              enddo
            enddo
          enddo
          do i=1,nmsd
            do imol=1,nmols
              do k=1,3
                qcm(i,imol,k)=qcm(i+nmsd,imol,k)
              enddo
            enddo
          enddo
          call runave2(nty,iread,nmsd,amsd,dt,elecon,tmsd,imsdcond,lcond,inm,chgmol)
        enddo
        deallocate (qcm,amsd)
        close (imsdcond)
        call xtc2 % close
        call system('/bin/rm -f cmdisp.xtc')
        deallocate (chgmol)
      endif ! inst=3

      if (inst==4) then
        open (401,file='rt.inp',status='old')
        read (401,*) refrence
        do k=1,nty
          if (trim(refrence)==trim(res_type(k))) kty=k
        enddo
        write(*,*) refrence
        read (401,*) nkind
        write(*,*) nkind
        allocate (mole(nkind),nat(nkind),imole(nkind))
        read(401,*) (mole(i),i=1,nkind)
        write(*,*)  (mole(i),i=1,nkind)
        do i=1,nkind
          do k=1,nty
            if (trim(mole(i))==trim(res_type(k))) then
              imole(i)=k
            endif
          enddo
        enddo
        write(*,*) (imole(i),i=1,nkind)
        read(401,*) (nat(j),j=1,nkind)
        write(*,*)  (nat(j),j=1,nkind)
        natmax=nat(1)
        do j=2,nkind
          if (natmax < nat(j)) natmax=nat(j)
        enddo
        allocate (iatsol(nkind,natmax),rdfmin(nkind,natmax))
        read (401,*) ((iatsol(i,j),j=1,nat(i)),i=1,nkind)
        write (*,*)  ((iatsol(i,j),j=1,nat(i)),i=1,nkind)
        read (401,*) ((rdfmin(i,j),j=1,nat(i)),i=1,nkind)
        write (*,*)  ((rdfmin(i,j),j=1,nat(i)),i=1,nkind)
        close(401)
        write(*,*) 'enter the maximum correlated time (ps):'
        read(*,*) trt
        write(*,*) 'enter the rt mode (Continuous(1) or Interrupted(2)):'
        read(*,*) mode
        nlen=0
        call xtc % init('traj_comp.xtc')
        do while (xtc % STAT==0)
          call xtc % read
          t1=xtc % time
          nlen=nlen+1
          if (nlen == 2) dt=t1
        enddo
        call xtc % close
        
        nrt=nint(trt/dt)
        write(*,*) dt,nrt
        allocate (coord(natoms,3),coord_ion(inm(kty),3))
        allocate (couplert(nrt+1,nkind,inm(kty),nsize,natmax))
        allocate (natmn(natmax),rkj(natmax))
        allocate (nh0(nkind),nht(nkind,nrt),ctn(nkind,0:nrt),ct(nkind,0:nrt),tau(nkind,0:nrt))
        ctn(:,:)=0.e0
        ct(:,:)=0.e0
        tau(:,:)=0.e0
        call xtc % init('traj_comp.xtc')
        do ilen=1,nlen
          call xtc % read
          cell=xtc % box
          cell=cell*10.e0
          call recip(cell,rcell,vol)
          do iatoms=1,natoms
            coord(iatoms,1:3)=xtc % pos(1:3,iatoms)*10.e0
          enddo
          natmk0=0
          do k=1,kty-1
            natmk0=natmk0+inm(k)*iat(k)
          enddo
          do k=1,inm(kty)
            coord_ion(k,:)=coord(natmk0+1+(k-1)*iat(kty),:)
          enddo
          do n=1,nkind
            natmn0=0
            do j=1,imole(n)-1
              natmn0=natmn0+inm(j)*iat(j)
            enddo
            do k=1,inm(kty)
              do j=1,inm(imole(n))
                do m=1,nat(n)
                  natmn(m)=natmn0+iatsol(n,m)+(j-1)*iat(imole(n))
                  rkj(m)=dist(coord_ion(k,1:3),coord(natmn(m),1:3),cell,rcell)
                  if (rkj(m) <= rdfmin(n,m)) then
                    couplert(mod(ilen-1,nrt+1)+1,n,k,j,m)=1
                  else
                    couplert(mod(ilen-1,nrt+1)+1,n,k,j,m)=0
                  endif
                enddo
              enddo
            enddo
          enddo
          
          if (ilen > nrt) then
            if (mode==1) then
              nh0(:)=0
              nht(:,:)=0
              ih0=mod(ilen-nrt-1,nrt+1)+1
              do n=1,nkind
                do k=1,inm(kty)
                  do j=1,inm(imole(n))
                    do m=1,nat(n)
                      if (couplert(ih0,n,k,j,m)==1) then
                        nh0(n)=nh0(n)+1
                        do l=1,nrt
                          iht=mod(ih0+l-1,nrt+1)+1
                          if (couplert(iht,n,k,j,m)==1) then
                            nht(n,l)=nht(n,l)+1
                          endif
                        enddo
                      endif
                    enddo
                  enddo
                enddo
                do l=1,nrt
                  ctn(n,l)=ctn(n,l)+real(nht(n,l))/real(nh0(n))
                  ct(n,l)=ctn(n,l)/real(ilen-nrt)
                enddo
              enddo
            endif
            if (mode==2) then
              nh0(:)=0
              nht(:,:)=0
              ih0=mod(ilen-nrt-1,nrt+1)+1
              do n=1,nkind
                do k=1,inm(kty)
                  do j=1,inm(imole(n))
                    do m=1,nat(n)
                      if (couplert(ih0,n,k,j,m)==1) then
                        nh0(n)=nh0(n)+1
                        do l=1,nrt
                          iht=mod(ih0+l-1,nrt+1)+1
                          if (iht > ih0) then
                            if (all(couplert(ih0:iht,n,k,j,m)==1)) then
                              nht(n,l)=nht(n,l)+1
                            else
                              exit  
                            endif
                          endif
                          if (iht < ih0) then
                            if (all(couplert(ih0:nrt+1,n,k,j,m)==1) .and. all(couplert(1:iht,n,k,j,m)==1)) then
                              nht(n,l)=nht(n,l)+1
                            else
                              exit
                            endif
                          endif
                        enddo
                      endif
                    enddo
                  enddo
                enddo
                do l=1,nrt
                  ctn(n,l)=ctn(n,l)+real(nht(n,l))/real(nh0(n))
                  ct(n,l)=ctn(n,l)/real(ilen-nrt)
                enddo
              enddo
            endif

            if (mod(ilen,100)==1 .and. ilen > nrt) then
              write (*,'(a10,f6.2,a)')  'Processs:', real((ilen-1)*100)/real(nlen),'%'
              open(402,file='rt.dat')
              open(403,file='tau.dat')
              ct(:,0)=1.e0
              tau(:,0)=0.e0
              do l=1,nrt
                write(402,'(f10.2,10f10.6)') real(l*dt),(ct(n,l),n=1,nkind)
              enddo
              do l=1,nrt
                do n=1,nkind
                  tau(n,l)=tau(n,l-1)+(ct(n,l-1)+ct(n,l))*0.5*dt
                enddo
                write(403,'(f10.2,10f20.6)') real(l*dt),(tau(n,l),n=1,nkind)
              enddo
              close(402)
              close(403)
            endif
          endif
        enddo 
        deallocate (coord,coord_ion,couplert)
        deallocate (nh0,nht,ctn,ct,tau)
        deallocate (natmn,rkj)
        deallocate (mole,imole,nat,iatsol,rdfmin)
      endif ! inst=4

      if (inst==5) then
        write(*,*) 'enter the correlated time (ps):'
        read(*,*) tvh
        call xtc % init("traj_comp.xtc")
        call xtc % read
        cell=xtc % box
        cell=cell*10.e0
        call recip(cell,rcell,vol)
        call xtc1 % init("cmcoor1.xtc", 'w')
        allocate (qcmol(3,nmols))
        nlen=0
        do while (xtc % STAT == 0)
          t1=xtc % time
          nlen=nlen+1
          if (nlen == 2) dt=t1
          iatm=0
          imol=0
          do i=1,nty
            allocate (coord(iat(i),3))
            do j=1,inm(i)
              imol=imol+1
              do k=1,iat(i)
                iatm=iatm+1
                coord(k,:)=xtc % pos(:,iatm)*10.e0
              enddo
              call cenmas(iat(i),coord,mass(i,1:iat(i)),qcm1,cell,rcell)
              qcmol(1:3,imol)=qcm1(1:3)
            enddo
            deallocate (coord)
          enddo
          call xtc1 % write(nmols,nlen,xtc % time,xtc % box*10.e0,qcmol,xtc % prec)
          call xtc % read
        enddo
        write(*,*)'nlen ',nlen
        call xtc % close
        call xtc1 % close
        deallocate (qcmol)
        
        call xtc1 % init("cmcoor1.xtc")
        call xtc2 % init("cmdisp.xtc", 'w')
        allocate (qcm(2,nmols,3),qdisp(3,nmols))
        qdisp(:,:)=0.
        call xtc1 % read
        do i=1,nmols
          qcm(1,i,1:3)=xtc1 % pos(1:3,i)
        enddo
        do ilen=1,nlen
          do i=1,nmols
            qcm(2,i,1:3)=xtc1 % pos(1:3,i)
          enddo
          call pbc2(nmols,qcm,qdisp,cell,rcell)
          imol=0
          do i=1,nty
            do j=1,inm(i)
              imol=imol+1
              do k=1,3
                qcm(1,imol,k)=qcm(2,imol,k)
              enddo
            enddo
          enddo
          call xtc2 % write(nmols,ilen-1,xtc1 % time,xtc1 % box,qdisp,xtc1 % prec)
          call xtc1 % read
        enddo
        call xtc1 % close
        call xtc2 % close
        deallocate (qdisp,qcm)
        call system('/bin/rm -f cmcoor1.xtc')
        
        nvh=nint(tvh/dt)
        nread=int(nlen/nvh)
        allocate (qcm(nvh*2,nmols,3))
        rijmax=0.
        call xtc2 % init("cmdisp.xtc")
        do i=1,nvh
          call xtc2 % read
          do imol=1,nmols
            qcm(i,imol,1:3)=xtc2 % pos(1:3,imol)
          enddo
        enddo
        do iread=2,nread
          do i=nvh+1,nvh*2
            call xtc2 % read
            do imol=1,nmols
              qcm(i,imol,1:3)=xtc2 % pos(1:3,imol)
            enddo
          enddo
          imol=0
          do ity=1,nty
            do iml=1,inm(ity)
              imol=imol+1
              do i=1,nvh
                rij=sqrt((qcm(i+nvh,imol,1)-qcm(i,imol,1))**2 &
                        +(qcm(i+nvh,imol,2)-qcm(i,imol,2))**2 &
                        +(qcm(i+nvh,imol,3)-qcm(i,imol,3))**2)
                if (rij >= rijmax) rijmax=rij
              enddo
            enddo
          enddo
          do i=1,nvh
            do imol=1,nmols
              do k=1,3
                qcm(i,imol,k)=qcm(i+nvh,imol,k)
              enddo
            enddo
          enddo
        enddo
        write(*,*) 'rijmax', rijmax
        call xtc2 % close
        
        nbinmax=nint(rijmax/delr)        
        allocate (area(nty),nvhf(nty,nbinmax))
        allocate (nvhf1(nty,nbinmax))
        nvhf(:,:)=0
        nvhf1(:,:)=0
        call xtc2 % init("cmdisp.xtc")
        do i=1,nvh
          call xtc2 % read
          do imol=1,nmols
            qcm(i,imol,1:3)=xtc2 % pos(1:3,imol)
          enddo
        enddo
        do iread=2,nread
          do i=nvh+1,nvh*2
            call xtc2 % read
            do imol=1,nmols
              qcm(i,imol,1:3)=xtc2 % pos(1:3,imol)
            enddo
          enddo
          imol=0
          do ity=1,nty
            do iml=1,inm(ity)
              imol=imol+1
              do i=1,nvh
                rij=sqrt((qcm(i+nvh,imol,1)-qcm(i,imol,1))**2 &
                        +(qcm(i+nvh,imol,2)-qcm(i,imol,2))**2 &
                        +(qcm(i+nvh,imol,3)-qcm(i,imol,3))**2)
                nbin=nint(rij/delr)
                if(nbin <= nbinmax) then
                  nvhf(ity,nbin)=nvhf(ity,nbin)+1
                endif
              enddo
            enddo
          enddo
          do i=1,nvh
            do imol=1,nmols
              do k=1,3
                qcm(i,imol,k)=qcm(i+nvh,imol,k)
              enddo
            enddo
          enddo
          write (*,'(a10,f6.2,a)')  'Processs:', real(iread*100)/real(nread),'%'
          if (mod(iread,10)==0 .or. iread == nread) then
            open (501,file='self_vh.dat')
            write(501,'(i10,a10,i10)') iread,'/' ,nread
            r0=0.5*delr
            area(:)=0.e0
            nvhf1(:,0)=0.e0
            do ity=1,nty
              do ibin=1,nbinmax
                r0=r0+delr
                nvhf1(ity,ibin)=4*pi*r0*r0*real(nvhf(ity,ibin))/real(nvh)/real(inm(ity))/real(iread-1)
                area(ity)=area(ity)+(nvhf1(ity,ibin)+nvhf1(ity,ibin-1))*delr/2
              enddo
            enddo
            r0=0.5*delr
            do ibin=1,nbinmax
              r0=r0+delr
              write (501,'(f10.4,12f16.8)') r0, (nvhf1(ity,ibin)/area(ity),ity=1,nty)
            enddo
            close(501)
          endif  
        enddo
        call xtc2 % close
        call system('/bin/rm -f cmdisp.xtc')
        deallocate (qcm,area,nvhf1)
        !deallocate (nvhf)
      endif ! inst=5

      if (inst==6) then
        write(*,*) 'enter the correlated time (ps):'
        read(*,*) tvh
        do i=1,nty
          write(*,'(i3,a,a)')  i,') ',res_type(i)
        enddo
        write(*,*) '------------>>'
        write(*,*) 'enter the kinds number:'
        read(*,*) nkind
        write(*,*) 'enter the refrenced ion or molecule:'
        read(*,*) ion
        do k=1,nty
          if (trim(ion)==trim(res_type(k))) kty=k
        enddo
        nmol0=0
        do k=1,kty-1
          nmol0=nmol0+inm(k)
        enddo  
        write(*,*) 'enter the calculated ion or molecule:'
        allocate (mole(nkind),imole(nkind),nmoln0(nkind))
        read(*,*) (mole(n),n=1,nkind)
        do n=1,nkind
          do k=1,nty
            if (trim(mole(n))==trim(res_type(k))) then
              imole(n)=k
            endif
          enddo
        enddo
        write(*,*) (imole(n),inm(imole(n)),n=1,nkind)
        nmoln0(:)=0
        do n=1,nkind
          do k=1,imole(n)-1
            nmoln0(n)=nmoln0(n)+inm(k)
          enddo
        enddo
        write(*,*) (nmoln0(n),n=1,nkind)

        call xtc % init("traj_comp.xtc")
        call xtc % read
        cell=xtc % box
        cell=cell*10.e0
        call recip(cell,rcell,vol)
        call xtc1 % init("cmcoor1.xtc", 'w')
        allocate (qcmol(3,nmols))
        nlen=0
        do while (xtc % STAT == 0)
          t1=xtc % time
          nlen=nlen+1
          if (nlen == 2) dt=t1
          iatm=0
          imol=0
          do i=1,nty
            allocate (coord(iat(i),3))
            do j=1,inm(i)
              imol=imol+1
              do k=1,iat(i)
                iatm=iatm+1
                coord(k,:)=xtc % pos(:,iatm)*10.e0
              enddo
              call cenmas(iat(i),coord,mass(i,1:iat(i)),qcm1,cell,rcell)
              qcmol(1:3,imol)=qcm1(1:3)
            enddo
            deallocate (coord)
          enddo
          call xtc1 % write(nmols,nlen,xtc % time,xtc % box*10.e0,qcmol,xtc % prec)
          call xtc % read
          cell=xtc % box
          cell=cell*10.e0
          call recip(cell,rcell,vol)
        enddo
        write(*,*)'nlen ',nlen
        call xtc % close
        call xtc1 % close
        deallocate (qcmol)

        nvh=nint(tvh/dt)
        nread=int(nlen/nvh)
        celllen1=sqrt(cell(1,1)**2+cell(2,1)**2+cell(3,1)**2)
        celllen2=sqrt(cell(1,2)**2+cell(2,2)**2+cell(3,2)**2)
        celllen3=sqrt(cell(1,3)**2+cell(2,3)**2+cell(3,3)**2)
        celllen=celllen1
        if(celllen > celllen2) celllen=celllen2
        if(celllen > celllen3) celllen=celllen3
        nlenrdf=nint(celllen/delr/2.e0)
        allocate (qcm(nvh*2,nmols,3),qcmk(nvh*2,inm(kty),3),qcmn(nkind,nvh,nsize,3))
        allocate (cnm(nkind),nvhf(nkind,nlenrdf))
        nvhf(:,:)=0
        call xtc1 % init("cmcoor1.xtc")
        do i=1,nvh
          call xtc1 % read
          do imol=1,nmols
            qcm(i,imol,1:3)=xtc1 % pos(1:3,imol)
          enddo
          do k=1,inm(kty)
            qcmk(i,k,1:3)=qcm(i,nmol0+k,1:3)
          enddo
        enddo
        
        do iread=2,nread
          do i=nvh+1,nvh*2
            call xtc1 % read
            do imol=1,nmols
              qcm(i,imol,1:3)=xtc1 % pos(1:3,imol)
            enddo
            do k=1,inm(kty)
              qcmk(i,k,1:3)=qcm(i,nmol0+k,1:3)
              do n=1,nkind
                do j=1,inm(imole(n))
                  qcmn(n,i-nvh,j,1:3)=qcm(i,nmoln0(n)+j,1:3)
                  rrkj=dist(qcmk(i-nvh,k,1:3),qcmn(n,i-nvh,j,1:3),cell,rcell)
                  nbin=nint(rrkj/delr)
                  if(nbin <= nlenrdf) then
                    nvhf(n,nbin)=nvhf(n,nbin)+1
                  endif
                enddo
              enddo
            enddo
          enddo
          do i=1,nvh
            do k=1,inm(kty)
              qcmk(i,k,1:3)=qcmk(i+nvh,k,1:3)
            enddo
          enddo
          write (*,'(a10,f6.2,a)')  'Processs:', real(iread*100)/real(nread),'%'
          if (mod(iread,10)==0 .or. iread==nread) then
            open(601,file='distinct_vh.dat')
            write(601,'(i10,a10,i10)') iread,'/',nread
            do n=1,nkind
              rho=real(inm(imole(n)))/vol
              cnm(n)=4.e0/3.e0*pi*rho*real(inm(kty))*real(nvh)
            enddo
            r0=0.5e0*delr
            r13=(0.5e0*delr)**3
            do ibin=1,nlenrdf
              r0=r0+delr
              r23=r0**3
              dr3=r23-r13
              r13=r23
              write(601,'(10f12.6)') r0,(real(nvhf(n,ibin))/cnm(n)/dr3/real(iread-1),n=1,nkind)
            enddo
            close(601)
          endif
        enddo
        call xtc1 % close
        call system('/bin/rm -f cmcoor1.xtc')
        deallocate (mole,imole,cnm,nmoln0)
        deallocate (qcm,qcmk,qcmn)
        !deallocate (nvhf)
      endif ! inst=6
      
      if (inst==7) then
        call xtc % init("traj_comp.xtc")
        call xtc % read
        cell=xtc % box
        cell=cell*10.e0
        call recip(cell,rcell,vol)
        call xtc % close
        do i=1,3
          r(i)=sqrt(dotp(cell(i,1:3),cell(i,1:3)))
        enddo
        rmax=r(1)
        if(rmax > r(2)) rmax=r(2)
        if(rmax > r(3)) rmax=r(3)
        ngrd=int(rmax/delg/2.e0)
        if(ngrd > ngrdmax) ngrd=ngrdmax
        allocate (ngrid(-ngrd:ngrd,-ngrd:ngrd,-ngrd:ngrd))
        allocate (igeomol(nty))
        allocate (lgeomol(nty))
        call random_seed()
        call random_number(urand)
        print *,urand
        do i=1,nty
          igeomol(i)=int(real(inm(i))*urand)+1
          lgeomol(i)=.false.
        enddo
        allocate (natmd(isize))
        allocate (ngeom(isize))
        allocate (geomol(nty,isize,3))
        allocate (qqq(nty,nsize,isize,3))
        iprt=900
        open(701,file='grid.inp',status='old')
        open(702,file='grid.checkp',form='unformatted')
        open(703,file='grid.checkp.2',form='unformatted',status='new')
        
        do job=1,3
          read(701,*) ncalc
          if (job==1 .and. ncalc > 0 .and. first) then
            call xtc % init("traj_comp.xtc")
            call xtc1 % init("cmcoor.xtc",'w')
            call xtc % read
            nlen=0
            allocate(qcmol(3,nmols))
            do while (xtc % STAT == 0)
              imol=0
              iatm=0
              do i=1,nty
                allocate(coord(iat(i),3))
                do j=1,inm(i)
                  imol=imol+1
                  do k=1,iat(i)
                    iatm=iatm+1
                    coord(k,1:3)=xtc % pos(1:3,iatm)*10.e0
                  enddo
                  call cenmas(iat(i),coord,mass(i,1:iat(i)),qcm1,cell,rcell)
                  qcmol(:,imol)=qcm1(:)
                enddo
                deallocate(coord)
              enddo
              call xtc1 % write (nmols, nlen, xtc % time, xtc % box*10.e0, qcmol, xtc % prec)
              nlen=nlen+1
              call xtc % read
            enddo
            call xtc % close
            call xtc1 % close
            deallocate(qcmol)
            allocate (qcm(nty,nsize,3))
            first=.false.
          else if (job==2 .and. .not.first) then
            deallocate (qcm)
            call system('/bin/rm -f cmcoor.xtc')
          endif

          do icalc=1,ncalc
            read(701,*) itp,jtp,ngm
            read(701,*) m1,m2,m3,m4
            read(701,*) (ngeom(i),i=1,ngm)
            if (job==3) then
              read(701,*) natd
              read(701,*)(natmd(i),i=1,natd)
            endif

            ! read intermediate results from checkpoint file
            nfm=0
            read(702,end=109) nfm
   109      continue
            do i=-ngrd,ngrd
              do j=-ngrd,ngrd
                do k=-ngrd,ngrd
                  if (nfm==0) then
                    ngrid(i,j,k)=0
                  else
                    read(702) ngrid(i,j,k)
                  endif
                enddo
              enddo
            enddo

            call xtc % init("traj_comp.xtc")
            call xtc % read
            if (job==1 .and. .not.first) then
              call xtc1 % init("cmcoor.xtc")
              call xtc1 % read
            endif

            do while (xtc % STAT == 0)
              imol=0
              iatm=0
              t1=xtc % time
              do i=1,nty
                do j=1,inm(i)
                  imol=imol+1
                  do k=1,iat(i)
                    iatm=iatm+1
                    qqq(i,j,k,:)=xtc % pos(:,iatm)*10.e0
                  enddo
                  if (job==1 .and. .not.first) qcm(i,j,1:3)=xtc1 % pos(1:3,imol)
                enddo
              enddo
              ! pbc corrections for individual molecules
              do i=1,nty
                allocate (coord(iat(i),3))
                do j=1,inm(i)
                  do k=1,iat(i)
                    do l=1,3
                      coord(k,l)=qqq(i,j,k,l)
                    enddo
                  enddo
                  call pbc1(iat(i),coord,cell,rcell)
                  do k=1,iat(i)
                    do l=1,3
                      qqq(i,j,k,l)=coord(k,l)
                    enddo
                  enddo
                enddo
                deallocate (coord)
              enddo
              ! grid calculation
              nfm=nfm+1
              do imol=1,inm(itp)
                do idrt=1,3
                  qori(idrt)=0.e0
                  do j=1,ngm
                    qori(idrt)=qori(idrt)+qqq(itp,imol,ngeom(j),idrt)
                  enddo
                  qori(idrt)=qori(idrt)/real(ngm)
                enddo
                call distvec(qqq(itp,imol,m1,1:3),qqq(itp,imol,m2,1:3),v1)
                call distvec(qqq(itp,imol,m3,1:3),qqq(itp,imol,m4,1:3),v2)
                dnorm=sqrt(dotp(v1,v1))
                do idrt=1,3
                  umol(1,idrt)=v1(idrt)/dnorm ! normalize
                enddo
                call crossp(umol(1,1:3),v2,umol(3,1:3))
                dnorm=sqrt(dotp(umol(3,1:3),umol(3,1:3)))
                do idrt=1,3
                  umol(3,idrt)=umol(3,idrt)/dnorm
                enddo
                call crossp(umol(3,1:3),umol(1,1:3),umol(2,1:3))
                if(abs(umol(2,1)**2+umol(2,2)**2+umol(2,3)**2-1.e0) &
                  .gt.1.e-6) stop 'unit vectors are ill defined!'
                ! store molecular coor in mol frame for grid file
                if (imol==igeomol(itp) .and. .not.lgeomol(itp)) then
                  allocate (coord(iat(itp),3))
                  do k=1,iat(itp)
                    do l=1,3
                      coord(k,l)=qqq(itp,imol,k,l)-qori(l)
                    enddo
                  enddo
                  call pbc1(iat(itp),coord,cell,rcell)
                  do k=1,iat(itp)
                    call transform(coord(k,1:3),umol)
                    do l=1,3
                      geomol(itp,k,l)=coord(k,l)
                    enddo
                  enddo
                  lgeomol(itp)=.true.
                  deallocate (coord)
                endif
              ! grid
                n=0
                if (job==1 .or. job==3) then
                  do jmol=1,inm(jtp)
                    if(itp==jtp .and. imol==jmol) goto 94
                    if (job==1) then
                      call distvec(qcm(jtp,jmol,1:3),qori,r)
                      call image(r,cell,rcell)
                      call transform(r,umol)
                      call grid(ngrd,1,ngrid,r,delg)
                    else if (job==3) then
                      do jatm=1,natd
                        call distvec(qqq(jtp,jmol,natmd(jatm),1:3),qori,r)
                        call image(r,cell,rcell)
                        call transform(r,umol)
                        call grid(ngrd,1,ngrid,r,delg)
                      enddo
                    endif
   94               n=n+1
                  enddo
                else if (job==2) then
                  do ity=1,nty
                    do jmol=1,inm(ity)
                      if(itp==ity .and. imol==jmol) goto 95
                      do jatm=1,iat(ity)
                        if (jtp==1 .and. charge(ity,jatm) > 0.e0) then
                          nchg=nint(charge(ity,jatm)*chgmul)
                          call distvec(qqq(ity,jmol,jatm,1:3),qori,r)
                          call image(r,cell,rcell)
                          call transform(r,umol)
                          call grid(ngrd,nchg,ngrid,r,delg)
                        else if(jtp.eq.-1 .and. charge(ity,jatm).lt.0.e0)then
                          nchg=-nint(charge(ity,jatm)*chgmul)
                          call distvec(qqq(ity,jmol,jatm,1:3),qori,r)
                          call image(r,cell,rcell)
                          call transform(r,umol)
                          call grid(ngrd,nchg,ngrid,r,delg)
                        endif
                      enddo
   95                 n=n+1
                    enddo
                  enddo
                endif
              enddo ! imol
              call xtc % read
              call xtc1 % read
            enddo

            ! write check point file
   96       write(703) nfm
            do i=-ngrd,ngrd
              do j=-ngrd,ngrd
                do k=-ngrd,ngrd
                  write(703) ngrid(i,j,k)
                enddo
              enddo
            enddo

            iprt=iprt+1
            dum=real(nfm)*real(inm(itp))*delg**3
            if(job==2) dum=dum*chgmul
            if(job==3) dum=dum*real(natd)
            ngd=2*ngrd+1
            ori=real(ngd)*delg/2.e0/ang2bohr
            write(iprt,*)
            write(iprt,*)
            write(iprt,'(i5,3f12.6)') iat(itp),-ori,-ori,-ori
            write(iprt,'(i5,3f12.6)') ngd,delg/ang2bohr,zero,zero
            write(iprt,'(i5,3f12.6)') ngd,zero,delg/ang2bohr,zero
            write(iprt,'(i5,3f12.6)') ngd,zero,zero,delg/ang2bohr
            do iatm=1,iat(itp) 
              write(iprt,'(i5,4f12.6)') &
                namu(mass(itp,iatm)),zero,(geomol(itp,iatm,l)/ang2bohr,l=1,3)
            enddo
            do i=-ngrd,ngrd
              do j=-ngrd,ngrd
                do k=-ngrd,ngrd,6
                  if (k+5+ngrd+1 > 2*ngrd+1) then
                    write(iprt,'(6e15.6)')(real(ngrid(i,j,l))/dum,l=k,ngrd)
                  else
                    write(iprt,'(6e15.6)')(real(ngrid(i,j,l))/dum,l=k,k+5)
                  endif
                enddo
              enddo
            enddo

            call xtc % close
            if(job==1 .and. .not.first) call xtc1 % close
          enddo  ! icalc
        enddo    ! job
        close(702)
        close(703)
        call system('mv grid.checkp.2 grid.checkp')
        call system('rm cmcoor.xtc')
      endif ! inst=7

      if (inst==9) then
        open (901,file='rt.inp',status='old')
        read (901,*) refrence
        do k=1,nty
          if (trim(refrence)==trim(res_type(k))) kty=k
        enddo
        write(*,*) refrence
        read (901,*) nkind
        write(*,*) nkind
        allocate (mole(nkind),nat(nkind),imole(nkind))
        read(901,*) (mole(i),i=1,nkind)
        write(*,*)  (mole(i),i=1,nkind)
        do i=1,nkind
          do k=1,nty
            if (trim(mole(i))==trim(res_type(k))) then
              imole(i)=k
            endif
          enddo
        enddo
        write(*,*) (imole(i),i=1,nkind)
        read(901,*) (nat(j),j=1,nkind)
        write(*,*)  (nat(j),j=1,nkind)
        natmax=nat(1)
        do j=2,nkind
          if (natmax < nat(j)) natmax=nat(j)
        enddo
        allocate (iatsol(nkind,natmax),rdfmin(nkind,natmax))
        iatsol(:,:)=0
        read (901,*) ((iatsol(i,j),j=1,nat(i)),i=1,nkind)
        write (*,*)  ((iatsol(i,j),j=1,nat(i)),i=1,nkind)
        read (901,*) ((rdfmin(i,j),j=1,nat(i)),i=1,nkind)
        write (*,*)  ((rdfmin(i,j),j=1,nat(i)),i=1,nkind)
        close(901)
        write(*,*) 'enter the correlated time (ps):'
        read(*,*) thop
        call xtc % init('traj_comp.xtc')
        nlen=0
        do while (xtc % STAT==0)
          call xtc % read
          t1=xtc % time
          nlen=nlen+1
          if (nlen == 2) dt=t1
        enddo
        call xtc % close
        
        nhop=nint(thop/dt)
        nread=int(nlen/nhop)
        allocate (natmn(natmax),rkj(natmax),coordt(nhop*2,natoms,3),coord_iont(nhop*2,inm(kty),3))
        open(902,file='hopping.dat')
        allocate (couplert(nhop*2,nkind,inm(kty),nsize,natmax),nhtp(nkind))
        couplert(:,:,:,:,:)=0
        nhtp(:)=0
        call xtc % init('traj_comp.xtc')
        do i=1,nhop
          call xtc % read
          cell=xtc % box
          cell=cell*10.e0
          call recip(cell,rcell,vol)
          do iatoms=1,natoms
            coordt(i,iatoms,1:3)=xtc % pos(1:3,iatoms)*10.e0
          enddo
          natmk0=0
          do k=1,kty-1
            natmk0=natmk0+inm(k)*iat(k)
          enddo
          do k=1,inm(kty)
            coord_iont(i,k,:)=coordt(i,natmk0+1+(k-1)*iat(kty),:)
          enddo
          do n=1,nkind
            natmn0=0
            do j=1,imole(n)-1
              natmn0=natmn0+inm(j)*iat(j)
            enddo
            do k=1,inm(kty)
              do j=1,inm(imole(n))
                do m=1,nat(n)
                  natmn(m)=natmn0+iatsol(n,m)+(j-1)*iat(imole(n))
                  rkj(m)=dist(coord_iont(i,k,1:3),coordt(i,natmn(m),1:3),cell,rcell)
                  if (rkj(m) <= rdfmin(n,m)) then
                    couplert(i,n,k,j,m)=1
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo

        do iread=2,nread
          do i=nhop+1,nhop*2
            call xtc % read
            cell=xtc % box
            cell=cell*10.e0
            call recip(cell,rcell,vol)
            do iatoms=1,natoms
              coordt(i,iatoms,1:3)=xtc % pos(1:3,iatoms)*10.e0
            enddo
            natmk0=0
            do k=1,kty-1
              natmk0=natmk0+inm(k)*iat(k)
            enddo
            do k=1,inm(kty)
              coord_iont(i,k,:)=coordt(i,natmk0+1+(k-1)*iat(kty),:)
            enddo
            do n=1,nkind
              natmn0=0
              do j=1,imole(n)-1
                natmn0=natmn0+inm(j)*iat(j)
              enddo
              do k=1,inm(kty)
                do j=1,inm(imole(n))
                  do m=1,nat(n)
                    natmn(m)=natmn0+iatsol(n,m)+(j-1)*iat(imole(n))
                    rkj(m)=dist(coord_iont(i,k,1:3),coordt(i,natmn(m),1:3),cell,rcell)
                    if (rkj(m) <= rdfmin(n,m)) then
                      couplert(i,n,k,j,m)=1
                    endif
                  enddo
                enddo
              enddo
            enddo
          enddo
          call random_number(rand)
          random_num = ceiling(rand*inm(kty))
          do n=1,nkind
            do k=random_num,random_num
              do j=1,inm(imole(n))
                do m=1,nat(n)
                  do i=1,nhop
                    if (all(couplert(i:i+nhop,n,k,j,m)==1)) then
                      nhtp(n)=nhtp(n)+0
                    endif
                    if (couplert(i,n,k,j,m)==0 .and. couplert(i+nhop,n,k,j,m)==1) then
                      nhtp(n)=nhtp(n)+1
                    endif
                    if (couplert(i,n,k,j,m)==1 .and. couplert(i+nhop,n,k,j,m)==1) then
                      if (.not. all(couplert(i:i+nhop,n,k,j,m)==1)) then
                        nhtp(n)=nhtp(n)-1
                      endif
                    endif
                  enddo
                enddo
              enddo
            enddo
          enddo
          write (902,'(11f20.4)') iread*nhop*dt,(real(nhtp(n))/real(nhop),n=1,nkind)
          do n=1,nkind
            do k=1,inm(kty)
              do j=1,inm(imole(n))
                do m=1,nat(n)
                  do i=1,nhop
                    couplert(i,n,k,j,m)=couplert(i+nhop,n,k,j,m)
                  enddo
                  couplert(nhop+1:nhop*2,n,k,j,m)=0
                enddo
              enddo
            enddo
          enddo
        enddo
        call xtc % close
        close (902)
        deallocate (mole,imole,nat,iatsol,rdfmin)
        deallocate (natmn,rkj,nhtp)
      endif ! inst=9

      if (inst==10) then
        do i=1,nty
          write(*,'(i3,a,a)')  i,') ',res_type(i)
        enddo
        write(*,*) '------------>>'
        open(1001,file='da.inp',status='old')
        read(1001,*) ion
        write(*,*) ion
        do k=1,nty
          if (trim(ion)==trim(res_type(k))) kty=k
        enddo
        nmolk0=0
        do k=1,kty-1
          nmolk0=nmolk0+inm(k)
        enddo
        read(1001,*) molecule
        write(*,*)   molecule
        read(1001,*) matm
        write(*,*)   matm
        allocate (iatan(1,matm),rdfmin(1,matm))
        read(1001,*) (iatan(1,m),m=1,matm)
        write(*,*)   (iatan(1,m),m=1,matm)
        read(1001,*) (rdfmin(1,m),m=1,matm)
        write(*,*)   (rdfmin(1,m),m=1,matm)
        write(*,*) 'enter the 2 atomic number'
        read(*,*) na1,na2
        write(*,*) 'enter the deld and dela:'
        read(*,*) deld,dela

        do k=1,nty
          if (trim(ion)==trim(res_type(k))) kty=k
        enddo
        nmolk0=0
        do k=1,kty-1
          nmolk0=nmolk0+inm(k)
        enddo
        do j=1,nty
          if (trim(molecule)==trim(res_type(j))) jty=j
        enddo
        nmolj0=0
        do j=1,jty-1
          nmolj0=nmolj0+inm(j)
        enddo
        natoma0=0
        do j=1,jty-1
          natoma0=natoma0+inm(j)*iat(j)
        enddo

        call xtc % init("traj_comp.xtc")
        call xtc % read
        cell=xtc % box
        cell=cell*10.e0
        call recip(cell,rcell,vol)
        call xtc1 % init("cmcoor1.xtc", 'w')
        allocate (qcmol(3,nmols))
        nlen=0
        do while (xtc % STAT == 0)
        nlen=nlen+1
          iatm=0
          imol=0
          do i=1,nty
            allocate (coord(iat(i),3))
            do j=1,inm(i)
              imol=imol+1
              do k=1,iat(i)
                iatm=iatm+1
                coord(k,:)=xtc % pos(:,iatm)*10.e0
              enddo
              call cenmas(iat(i),coord,mass(i,1:iat(i)),qcm1,cell,rcell)
              qcmol(1:3,imol)=qcm1(1:3)
            enddo
            deallocate (coord)
          enddo
          call xtc1 % write(nmols,nlen,xtc % time,xtc % box*10.e0,qcmol,xtc % prec)
          call xtc % read
        enddo
        write(*,*)'nlen ',nlen
        call xtc % close
        call xtc1 % close
        deallocate (qcmol)
        
        call xtc % init("traj_comp.xtc")
        call xtc1 % init("cmcoor1.xtc")
        allocate (nda(0:nint(20/deld),0:nint(180/dela)))
        nda(:,:)=0
        do ilen=1,nlen
          call xtc % read
          call xtc1 % read
          allocate (coord(natoms,3),qcm(1,nmols,3),qcmk(1,inm(kty),3))
          allocate (natmn(matm),rkj(matm),couple_anion(1,inm(kty),inm(jty)))
          couple_anion(1,:,:)=0
          do iatm=1,natoms
            coord(iatm,1:3)=xtc % pos(1:3,iatm)*10.e0
          enddo
          do imol=1,nmols
            qcm(1,imol,1:3)=xtc1 % pos(1:3,imol)
          enddo
          do k=1,inm(kty)
            qcmk(1,k,1:3)=qcm(1,nmolk0+k,1:3)
            do j=1,inm(jty)
              do m=1,matm
                natmn(m)=natoma0+iatan(1,m)+(j-1)*iat(jty)
                rkj(m)=dist(qcmk(1,k,1:3),coord(natmn(m),1:3),cell,rcell)
                if (rkj(m) <= rdfmin(1,m)) then
                  couple_anion(1,k,j)=1
                  exit
                endif
              enddo
              if (couple_anion(1,k,j)==1) then
                dkj=dist(qcmk(1,k,1:3),qcm(1,nmolj0+j,1:3),cell,rcell)
                nbin=nint(dkj/deld)
                natoma1=natoma0+na1+(j-1)*iat(jty)
                natoma2=natoma0+na2+(j-1)*iat(jty)
                da1a2=dist(coord(natoma1,1:3),coord(natoma2,1:3),cell,rcell)
                v1(1)=qcmk(1,k,1)-qcm(1,nmolj0+j,1)
                v1(2)=qcmk(1,k,2)-qcm(1,nmolj0+j,2)
                v1(3)=qcmk(1,k,3)-qcm(1,nmolj0+j,3)
                v2(1)=coord(natoma1,1)-coord(natoma2,1)
                v2(2)=coord(natoma1,2)-coord(natoma2,2)
                v2(3)=coord(natoma1,3)-coord(natoma2,3)
                call pbc3(v1,cell,rcell)
                call pbc3(v2,cell,rcell)
                product=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
                cos_angle=product/(dkj*da1a2)
                cos_angle=max(-1.0,min(1.0, cos_angle))
             !   write(*,*) product,dkj*da1a2,cos_angle,dkj
                angle=acos(cos_angle)*180.0/pi
             !   if (angle > 90.0) angle=180.0-angle
                mbin=nint(angle/dela)
                nda(nbin,mbin)=nda(nbin,mbin)+1
              endif
            enddo
          enddo
          
          nsum=0
          do n=0,nint(6/deld)
            do m=0,nint(180/dela)
              nsum=nsum+nda(n,m)
            enddo
          enddo
          if (mod(ilen,nskip)==0) then
            write (*,'(a10,f6.2,a)')  'Processs:', real((ilen)*100)/real(nlen),'%'
            open(1002,file='da.dat')
            do n=0,nint(6/deld)
              write(1002,'(1000f10.6)') (real(nda(n,m))/nsum,m=0,nint(180/dela))
            enddo
            close(1002)
          endif
          deallocate (coord,qcm,qcmk,natmn,couple_anion,rkj)
        enddo
        deallocate (iatan,rdfmin)
        deallocate (nda)
        call xtc % close
        call xtc1 % close
      endif ! inst=10

      if (inst==12) then
        do i=1,nty
          write(*,'(i3,a,a)')  i,') ',res_type(i)
        enddo
        write(*,*) '------------>>'
        open (1201,file='ti.inp',status='old')
        read(1201,*) refrence
        do k=1,nty
          if (trim(refrence)==trim(res_type(k))) kty=k
        enddo
        write(*,*) refrence
        read (1201,*) nkind
        write(*,*) nkind
        allocate (mole(nkind),nat(nkind),imole(nkind))
        read(1201,*) (mole(i),i=1,nkind)
        do i=1,nkind
          do k=1,nty
            if (trim(mole(i))==trim(res_type(k))) then
              imole(i)=k
            endif
          enddo
        enddo
        write(*,*) (imole(i),i=1,nkind)
        read(1201,*) (nat(j),j=1,nkind)
        write(*,*)   (nat(j),j=1,nkind)
        natmax=nat(1)
        do j=2,nkind
          if (natmax < nat(j)) natmax=nat(j)
        enddo
        allocate (iatsol(nkind,natmax),rdfmin(nkind,natmax))
        read (1201,*) ((iatsol(i,j),j=1,nat(i)),i=1,nkind)
        write (*,*)   ((iatsol(i,j),j=1,nat(i)),i=1,nkind)
        read (1201,*) ((rdfmin(i,j),j=1,nat(i)),i=1,nkind)
        write (*,*)   ((rdfmin(i,j),j=1,nat(i)),i=1,nkind)
        write(*,*) 'enter the correlated time (ps):'
        read(*,*) tti
        nlen=0
        call xtc % init("traj_comp.xtc")
        call xtc % read
        cell=xtc % box
        cell=cell*10.e0
        call recip(cell,rcell,vol)
        call xtc1 % init("cmcoor1.xtc", 'w')

        do while (xtc % STAT == 0)
          t1=xtc % time
          nlen=nlen+1
          if (nlen == 2) dt=t1
          allocate (coord(3,natoms),coord_ion(3,inm(kty)))
          do iatoms=1,natoms
            coord(1:3,iatoms)=xtc % pos(1:3,iatoms)*10.e0
          enddo
          natmk0=0
          do k=1,kty-1
            natmk0=natmk0+inm(k)*iat(k)
          enddo
          do k=1,inm(kty)
            coord_ion(:,k)=coord(:,natmk0+1+(k-1)*iat(kty))
          enddo
          call xtc1 % write(inm(kty),nlen,xtc % time,xtc % box*10.e0,coord_ion,xtc % prec)
          call xtc % read
          deallocate (coord,coord_ion)
        enddo
        write(*,*)'nlen ',nlen
        write(*,*)'dt   ',dt
        call xtc % close
        call xtc1 % close

        call xtc1 % init("cmcoor1.xtc")
        call xtc2 % init("cmdisp.xtc", 'w')
        allocate (qcm(2,inm(kty),3),qdisp(3,inm(kty)))
        qdisp(:,:)=0.
        call xtc1 % read
        do k=1,inm(kty)
          qcm(1,k,1:3)=xtc1 % pos(1:3,k)
        enddo
        do ilen=1,nlen
          do k=1,inm(kty)
            qcm(2,k,1:3)=xtc1 % pos(1:3,k)
          enddo
          call pbc2(inm(kty),qcm,qdisp,cell,rcell)
          do k=1,inm(kty)
            do i=1,3
              qcm(1,k,i)=qcm(2,k,i)
            enddo
          enddo
          call xtc2 % write(inm(kty),ilen-1,xtc1 % time,xtc1 % box,qdisp,xtc1 % prec)
          call xtc1 % read
        enddo
        call xtc1 % close
        call xtc2 % close
        deallocate (qdisp,qcm)

        tin=0.e0
        open(1202,file='ti.dat')
        nti=nint(tti/dt)
        nread=int(nlen/nti)
        allocate (coord(natoms,3),coord_ion(inm(kty),3),qcm(nti*2,inm(kty),3))
        allocate (natmn(natmax),rkj(natmax))
        allocate (couplert(nti*2,nkind,inm(kty),nsize,natmax))
        call xtc % init("traj_comp.xtc")
        call xtc1 % init("cmcoor1.xtc")
        call xtc2 % init("cmdisp.xtc")
        do i=1,nti
          call xtc  % read
          call xtc1 % read
          call xtc2 % read
          do iatoms=1,natoms
            coord(iatoms,1:3)=xtc % pos(1:3,iatoms)*10.e0
          enddo
          do k=1,inm(kty)
            coord_ion(k,:)=xtc1 % pos(1:3,k)
            qcm(i,k,1:3)=xtc2 % pos(1:3,k)
          enddo
          do n=1,nkind
            natmn0=0
            do j=1,imole(n)-1
              natmn0=natmn0+inm(j)*iat(j)
            enddo
            do k=1,inm(kty)
              do j=1,inm(imole(n))
                do m=1,nat(n)
                  natmn(m)=natmn0+iatsol(n,m)+(j-1)*iat(imole(n))
                  rkj(m)=dist(coord_ion(k,1:3),coord(natmn(m),1:3),cell,rcell)
                  if (rkj(m) <= rdfmin(n,m)) then
                    couplert(i,n,k,j,m)=1
                  else
                    couplert(i,n,k,j,m)=0
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
        
        do iread=2,nread
          do i=nti+1,nti*2
            call xtc % read
            call xtc1 % read
            call xtc2 % read
            do iatoms=1,natoms
              coord(iatoms,1:3)=xtc % pos(1:3,iatoms)*10.e0
            enddo
            do k=1,inm(kty)
              coord_ion(k,:)=xtc1 % pos(1:3,k)
              qcm(i,k,1:3)=xtc2 % pos(1:3,k)
            enddo
            do n=1,nkind
              natmn0=0
              do j=1,imole(n)-1
                natmn0=natmn0+inm(j)*iat(j)
              enddo
              do k=1,inm(kty)
                do j=1,inm(imole(n))
                  do m=1,nat(n)
                    natmn(m)=natmn0+iatsol(n,m)+(j-1)*iat(imole(n))
                    rkj(m)=dist(coord_ion(k,1:3),coord(natmn(m),1:3),cell,rcell)
                    if (rkj(m) <= rdfmin(n,m)) then
                      couplert(i,n,k,j,m)=1
                    else
                      couplert(i,n,k,j,m)=0
                    endif
                  enddo
                enddo
              enddo
            enddo
          enddo
          
          vr2=0.e0
          sr2=0.e0
          do i=1,nti
            do k=1,inm(kty)
              rkk2=(qcm(i+nti,k,1)-qcm(i,k,1))**2 &
                  +(qcm(i+nti,k,2)-qcm(i,k,2))**2 &
                  +(qcm(i+nti,k,3)-qcm(i,k,3))**2
              vehicular=.true.
              do n=1,nkind
                do j=1,inm(imole(n))
                  do m=1,nat(n)
                    do l=1,nti
                      if (couplert(i+l,n,k,j,m) /= couplert(i,n,k,j,m)) then
                        vehicular=.false.
                        exit
                      endif
                    enddo
                    if (.not. vehicular) exit
                  enddo
                  if (.not. vehicular) exit
                enddo
                if (.not. vehicular) exit
              enddo
              if (vehicular) vr2=vr2+rkk2
              if (.not. vehicular) sr2=sr2+rkk2
            enddo
          enddo
          tin=tin+real(sr2)/real((vr2+sr2))
          ti=tin/real(iread-1)
          do i=1,nti
            couplert(i,:,:,:,:)=couplert(i+nti,:,:,:,:)
            qcm(i,:,:)=qcm(i+nti,:,:)
          enddo
          if (mod(iread*nti,nskip)==0) then
            write (*,'(a10,f6.2,a)')  'Processs:', real(iread*100)/real(nread),'%'
            write(1202,'(i10,5f20.10)') int(iread*nti*dt),ti
          endif
        enddo
        call xtc  % close
        call xtc1 % close
        call xtc2 % close
        deallocate (mole,imole,nat)
        deallocate (coord,coord_ion,qcm,couplert)
        deallocate (natmn,rkj,iatsol,rdfmin)
        close(1201)
        close(1202)
        call system('/bin/rm -f cmcoor1.xtc')
        call system('/bin/rm -f cmdisp.xtc')
      endif

      deallocate (res_type,inm,iat)
      deallocate (atom_type,element,sigma,epsilon,atom_type_itp,mass,charge,chgsum)
      print *,"finished"
      end program gmxkit

      subroutine recip(cell,rcell,vol)
      implicit real (a-h,o-z)
      dimension cell(3,3),rcell(3,3)
      rcell(1,1)=cell(2,2)*cell(3,3)-cell(3,2)*cell(2,3)
      rcell(2,1)=cell(3,2)*cell(1,3)-cell(1,2)*cell(3,3)
      rcell(3,1)=cell(1,2)*cell(2,3)-cell(2,2)*cell(1,3)
      rcell(1,2)=cell(2,3)*cell(3,1)-cell(3,3)*cell(2,1)
      rcell(2,2)=cell(3,3)*cell(1,1)-cell(1,3)*cell(3,1)
      rcell(3,2)=cell(1,3)*cell(2,1)-cell(2,3)*cell(1,1)
      rcell(1,3)=cell(2,1)*cell(3,2)-cell(2,2)*cell(3,1)
      rcell(2,3)=cell(3,1)*cell(1,2)-cell(1,1)*cell(3,2)
      rcell(3,3)=cell(1,1)*cell(2,2)-cell(2,1)*cell(1,2)
      vol=cell(1,1)*rcell(1,1)+cell(2,1)*rcell(2,1)+cell(3,1)*rcell(3,1)
      do i=1,3
        do j=1,3
          rcell(i,j)=rcell(i,j)/vol
        enddo
      enddo
      return
      end

      function dist(q1,q2,cell,rcell)
      implicit real (a-h,o-z)
      dimension q1(3),q2(3),cell(3,3),rcell(3,3)
      dimension dr(3),fr(3)
      do i=1,3
        dr(i)=q1(i)-q2(i)
      enddo
      do i=1,3
        fr(i)=rcell(1,i)*dr(1)+rcell(2,i)*dr(2)+rcell(3,i)*dr(3)
      enddo
      do i=1,3
        fr(i)=fr(i)-anint(fr(i))
      enddo
      do i=1,3
        dr(i)=cell(i,1)*fr(1)+cell(i,2)*fr(2)+cell(i,3)*fr(3)
      enddo
      dist=sqrt(dr(1)**2+dr(2)**2+dr(3)**2)
      return
      end
      
      subroutine pbc1(n,coord,cell,rcell)
      implicit real (a-h,o-z)
      dimension coord(n,3)
      dimension cell(3,3),rcell(3,3),fr(3)
      do iatm=2,n
        do i=1,3
          fr(i)=rcell(1,i)*(coord(iatm,1)-coord(1,1)) &
               +rcell(2,i)*(coord(iatm,2)-coord(1,2)) &
               +rcell(3,i)*(coord(iatm,3)-coord(1,3))
        enddo
        do i=1,3
          fr(i)=fr(i)-anint(fr(i))
        enddo
        do i=1,3
          coord(iatm,i)=coord(1,i)+cell(i,1)*fr(1)+cell(i,2)*fr(2)+cell(i,3)*fr(3)
        enddo
      enddo
      return
      end
      
      subroutine pbc2(nmols,qcm,qdisp,cell,rcell)
      implicit real (a-h,o-z)
      dimension qcm(2,nmols,3),qdisp(3,nmols)
      dimension cell(3,3),rcell(3,3),fr(3)
      do imol=1,nmols
        do i=1,3
          fr(i)=rcell(1,i)*(qcm(2,imol,1)-qcm(1,imol,1)) &
               +rcell(2,i)*(qcm(2,imol,2)-qcm(1,imol,2)) &
               +rcell(3,i)*(qcm(2,imol,3)-qcm(1,imol,3))
        enddo
        do i=1,3
          fr(i)=fr(i)-anint(fr(i))
        enddo
        do i=1,3
          qdisp(i,imol)=qdisp(i,imol)+cell(i,1)*fr(1)+cell(i,2)*fr(2)+cell(i,3)*fr(3)
        enddo
      enddo
      return
      end

      subroutine pbc3(v,cell,rcell)
      implicit real (a-h,o-z)
      dimension v(3)
      dimension cell(3,3),rcell(3,3),fr(3)
      do i=1,3
        fr(i)=rcell(1,i)*v(1)+rcell(2,i)*v(2)+rcell(3,i)*v(3)
      enddo
      do i=1,3
        fr(i)=fr(i)-anint(fr(i))
      enddo
      do i=1,3
        v(i)=cell(i,1)*fr(1)+cell(i,2)*fr(2)+cell(i,3)*fr(3)
      enddo
      return
      end

      subroutine cenmas(n,coord,awt,qcm1,cell,rcell)
      implicit real (a-h,o-z)
      dimension coord(n,3),awt(n)
      dimension qcm1(3),cell(3,3),rcell(3,3)
      if (n == 1) then
        qcm1(1)=coord(1,1)
        qcm1(2)=coord(1,2)
        qcm1(3)=coord(1,3)
        return
      endif
      if (n > 1) then
        call pbc1(n,coord,cell,rcell)
        qcm1(1)=0.
        qcm1(2)=0.
        qcm1(3)=0.
        wt=0.
        do i=1,n
          wt=wt+awt(i)
          do j=1,3
            qcm1(j)=qcm1(j)+awt(i)*coord(i,j)
          enddo
        enddo
        qcm1(1)=qcm1(1)/wt
        qcm1(2)=qcm1(2)/wt
        qcm1(3)=qcm1(3)/wt
      endif
      return
      end

      subroutine distvec(v1,v2,v3)
      implicit real (a-h,o-z)
      dimension v1(3),v2(3),v3(3)
      do i=1,3
        v3(i)=v1(i)-v2(i)
      enddo
      return
      end
      
      function dotp(v1,v2)
      implicit real (a-h,o-z)
      dimension v1(3),v2(3)
      dotp=0.e0
      do i=1,3
        dotp=dotp+v1(i)*v2(i)
      enddo
      return
      end

      subroutine crossp(v1,v2,v3)
      implicit real (a-h,o-z)
      dimension v1(3),v2(3),v3(3)
      v3(1)=epsil(2,3,1)*v1(2)*v2(3)+epsil(3,2,1)*v1(3)*v2(2)
      v3(2)=epsil(3,1,2)*v1(3)*v2(1)+epsil(1,3,2)*v1(1)*v2(3)
      v3(3)=epsil(1,2,3)*v1(1)*v2(2)+epsil(2,1,3)*v1(2)*v2(1)
      return
      end

      function epsil(i,j,k)
      implicit real (a-h,o-z)
      if(i.eq.1 .and. j.eq.2 .and. k.eq.3) epsil= 1.e0
      if(i.eq.2 .and. j.eq.1 .and. k.eq.3) epsil=-1.e0
      if(i.eq.2 .and. j.eq.3 .and. k.eq.1) epsil= 1.e0
      if(i.eq.3 .and. j.eq.2 .and. k.eq.1) epsil=-1.e0
      if(i.eq.3 .and. j.eq.1 .and. k.eq.2) epsil= 1.e0
      if(i.eq.1 .and. j.eq.3 .and. k.eq.2) epsil=-1.e0
      return
      end
      
      subroutine image(r,cell,rcell)
      implicit real (a-h,o-z)
      dimension cell(3,3),rcell(3,3)
      dimension r(3),fr(3)
      do i=1,3
        fr(i)=rcell(i,1)*r(1)+rcell(i,2)*r(2)+rcell(i,3)*r(3)
      enddo
      do i=1,3
        fr(i)=fr(i)-anint(fr(i))
      enddo
      do i=1,3
        r(i)=cell(1,i)*fr(1)+cell(2,i)*fr(2)+cell(3,i)*fr(3)
      enddo
      return
      end

      subroutine transform(r,umol)
      implicit real (a-h,o-z)
      dimension r(3),rmol(3),umol(3,3)
      do i=1,3
        rmol(i)=0.e0
        do j=1,3
          rmol(i)=rmol(i)+umol(i,j)*r(j)
        enddo
      enddo
      do i=1,3
        r(i)=rmol(i)
      enddo
      return
      end
      
      subroutine grid(ngrd,nct,ngrid,r,delg)
      implicit real (a-h,o-z)
      dimension r(3)
      dimension n(3)
      dimension ngrid(-ngrd:ngrd,-ngrd:ngrd,-ngrd:ngrd)
      do i=1,3
        n(i)=nint(r(i)/delg)
        if(abs(n(i)) > ngrd) return
      enddo
      ngrid(n(1),n(2),n(3))=ngrid(n(1),n(2),n(3))+nct
      return
      end

      function namu(amu)
      implicit real (a-h,o-z)
      if (nint(amu)==1) then
        namu=1
      else if (nint(amu).eq.12) then
        namu=6
      else if (nint(amu).eq.14) then
        namu=7
      else if (nint(amu).eq.16) then
        namu=8
      else if (nint(amu).eq.19) then
        namu=9
      else
        stop 'unknown amu!'
      endif
      return
      end

      subroutine remove_duplicates(array,n)
      implicit real (a-h,o-z)
      integer :: array(n)
      integer :: unique_array(n)
      logical :: found
      k = 0
      do i=1, n
        found = .false.
        do j=1,k
          if (array(i) == unique_array(j)) then
            found = .true.
            exit
          end if
        end do
        if (.not. found) then
          k=k+1
          unique_array(k) = array(i)
        end if
      enddo
      n=k
      array=unique_array
      return
      end
  
      subroutine remove_2_duplicates(array1,n,array2,m)
      implicit real (a-h,o-z)
      integer :: array1(n),array2(m) 
      integer :: unique_array1(n),unique_array2(m)
      logical :: found
      i=1
      do while (i <=  m)
        found = .false.
        do j=1,n
          if (array1(j) == array2(i)) then
            found = .true.
            exit
          endif
        enddo
        if (found) then
          do k=i,m-1
            array2(k)=array2(k+1)
          enddo
          m=m-1
        else
          i=i+1
        endif
      enddo
      return
      end
      
      subroutine runave(nty,iread,nmsd,cmsd,dt,elecon,tmsd,imsdcond,res_type)
      implicit real (a-h,o-z)
      real*8 cmsd(nty,nty,0:nmsd)
      real*8 fac,cond_time(nty,nty),cond_slope(nty,nty),cond_ind(nty+1)
      real*8, allocatable :: cmsd0(:,:)
      character res_type(nty)
      logical first
      data first/.true./
      save first
      
      if (first) then
        call system('/bin/rm -f running_transfernumber_time.dat')
        call system('/bin/rm -f running_conductivity_time.dat')
        call system('/bin/rm -f running_transfernumber_slope.dat')
        call system('/bin/rm -f running_conductivity_slope.dat')
        call system('/bin/rm -f running_ind_cond_time1.dat')
        call system('/bin/rm -f running_ind_cond_slope1.dat')
        call system('/bin/rm -f running_ind_cond_time.dat')
        call system('/bin/rm -f running_ind_cond_slope.dat')
        first=.false.
      endif
      
      open(45,file='chgmsd.dat')
      open(55,file='chgmsd_ind.dat')
      open(56,file='running_transfernumber_time.dat',position='append')
      open(57,file='running_conductivity_time.dat',position='append')
      open(58,file='running_transfernumber_slope.dat',position='append')
      open(59,file='running_conductivity_slope.dat',position='append')
      open(60,file='running_ind_cond_time.dat',position='append')
      open(61,file='running_ind_cond_slope.dat',position='append')      
      ! time
      allocate (cmsd0(nty+1,0:nmsd))
      cmsd0(:,:)=0.d0
      fac=1.d0/dble(nmsd)/dble(iread-1)
      do i=0,nmsd
        do ity=1,nty
          do jty=ity,nty
            cmsd0(ity,i)=cmsd0(ity,i)+cmsd(ity,jty,i)
            if(ity /= jty) then
              cmsd0(jty,i)=cmsd0(jty,i)+cmsd(ity,jty,i)
            endif
          enddo
          cmsd0(nty+1,i)=cmsd0(nty+1,i)+cmsd0(ity,i)
        enddo
        write(45,'(10f20.8)') real(i*dt),(cmsd0(ity,i)*fac,ity=1,nty+1)
        write(55,'(27f20.8)') real(i*dt),((cmsd(ity,jty,i)*fac,jty=ity,nty),ity=1,nty)
      enddo
      close(45)
      close(55)

      fac=fac*dble(elecon)/dble(tmsd)/2.d0
      do ity=1,nty
        do jty=ity,nty
          cond_time(ity,jty)=cmsd(ity,jty,nmsd)*fac
        enddo
        cond_ind(ity)=cmsd0(ity,nmsd)*fac
      enddo
      cond_ind(nty+1)=cmsd0(nty+1,nmsd)*fac
      deallocate (cmsd0)

      write(imsdcond,'(/,a,i6)')'charge MSD, iread = ',iread
      write(imsdcond,'(99hconductivity by total charge displacement divided &
            by time (Smith et al. JPCB 110(2006)22773 eq.8) =,f12.8,5h S/cm)') cond_ind(nty+1)
      write(imsdcond,*)'individual ionic conductivity (S/cm) and transference number by time:'
      do ity=1,nty
        write(imsdcond,'(i5,3x,2f12.8)')ity, cond_ind(ity),cond_ind(ity)/cond_ind(nty+1)
      enddo
      write(57,'(10g20.10)') real(iread-1)*tmsd,cond_ind(nty+1),(cond_ind(ity),ity=1,nty)
      write(56,'(10g20.10)') real(iread-1)*tmsd,(cond_ind(ity)/cond_ind(nty+1),ity=1,nty)
      write(imsdcond,*)'self/distinct and cross conductivity by charge displacement by time:'
      do ity=1,nty
        do jty=ity,nty
          write(imsdcond,'(2i5,3x,f12.8, 5h S/cm)')ity, jty, cond_time(ity,jty)
        enddo
      enddo
      write(60,'(27g20.10)') real(iread-1)*tmsd,((cond_time(ity,jty),jty=ity,nty),ity=1,nty)

      n1=nmsd/2
      n2=nmsd*3/4
      nt=n2-n1
      cond_slope(:,:)=0.d0
      cond_ind(:)=0.d0

      fac=1.d0/dble(nmsd)/dble(iread-1)/dble(nt+1)/dble(nt)/dble(dt)*dble(elecon)/2.d0
      do ity=1,nty
        do jty=ity,nty
          do k=n1,n2
            cond_slope(ity,jty)=cond_slope(ity,jty)+(cmsd(ity,jty,k+nt)-cmsd(ity,jty,k))
          enddo
          cond_slope(ity,jty)=cond_slope(ity,jty)*fac
          cond_ind(ity)=cond_ind(ity)+cond_slope(ity,jty)
          if(ity /= jty) then
            cond_ind(jty)=cond_ind(jty)+cond_slope(ity,jty)
          endif
        enddo
        cond_ind(nty+1)=cond_ind(nty+1)+cond_ind(ity)
      enddo

      write(imsdcond,'(/,23hconductivity by slope =,f12.8,5h S/cm)') cond_ind(nty+1)
      write(imsdcond,*)'individual ionic conductivity (S/cm) and transference number by slope:'
      do ity=1,nty
        write(imsdcond,'(i5,3x,2f12.8)')ity, cond_ind(ity),cond_ind(ity)/cond_ind(nty+1)
      enddo
      write(59,'(10g20.10)')real(iread-1)*tmsd,cond_ind(nty+1),(cond_ind(ity),ity=1,nty)
      write(58,'(10g20.10)')real(iread-1)*tmsd,(cond_ind(ity)/cond_ind(nty+1),ity=1,nty)
      write(imsdcond,*)'self/distinct and cross conductivity by charge displacement by slope:'
      do ity=1,nty
        do jty=ity,nty
          write(imsdcond,'(2i5,3x,f12.8, 5h S/cm)')ity,jty,cond_slope(ity,jty)
        enddo
      enddo
      write(61,'(27g20.10)')real(iread-1)*tmsd,((cond_slope(ity,jty),jty=ity,nty),ity=1,nty)
      write(imsdcond,*)'conductivity by slope should be adopted!'

      call flush(imsdcond)
      close(56)
      close(57)
      close(58)
      close(59)
      close(60)
      close(61)

      return
      end
      
      subroutine runave2(nty,iread,nmsd,amsd,dt,elecon,tmsd,imsdcond,lcond,inm,chgmol)
      implicit real(a-h,o-z)
      dimension inm(nty)
      dimension chgmol(nty)
      real*8 amsd(nty,0:nmsd),dmsd(2,nty)
      real*8 fac,cond_time(nty,nty),cond_slope(nty,nty)
      real*8, allocatable :: cmsd(:,:,:)
      logical lcond,first
      data first/.true./
      save first
      if (first) then
        call system('/bin/rm -f running_self_diffusion_coeff.dat')
        if(lcond) call system('/bin/rm -f running_conductivity_NE.dat')
        if(lcond) call system('/bin/rm -f running_ind_cond_time1.dat')
        if(lcond) call system('/bin/rm -f running_ind_cond_slope1.dat')
        write(*,'(a,x,i10,x,f12.6)')'in runave2, nmsd, dt:',nmsd,dt
        first=.false.
      endif
      open(60,file='running_ind_cond_time.dat')
      open(61,file='running_ind_cond_slope.dat')
      open(65,file='running_self_diffusion_coeff.dat',position='append')
      open(66,file='msd.dat')
      open(67,file='running_ind_cond_time1.dat')
      open(68,file='running_ind_cond_slope1.dat')
      open(77,file='running_conductivity_NE.dat',position='append')
      ! msd and self 
      fac=dble(nmsd)*dble(iread-1)
      do i=0,nmsd
        write(66,'(10g20.10)')real(i)*dt,(amsd(j,i)/fac/dble(inm(j)),j=1,nty)
      enddo
      close(66)
      write(imsdcond,'(/,a,i6,//)')'self MSD, iread = ',iread
      n1=nmsd/2
      n2=nmsd*3/4
      nt=n2-n1
      write(imsdcond,*)'        mol      D(A^2/ps)      D(cm^2/s)      D(m^2/s)'
      do itp=1,nty
        dmsd(1,itp)=amsd(itp,nmsd)/fac/dble(inm(itp))/dble(tmsd)/6.d0
        write(imsdcond,'(i3,3g20.8,24h  MSD by divided by time)')itp, &
              dmsd(1,itp),dmsd(1,itp)/1.d4,dmsd(1,itp)/1.d8
        dmsd(2,itp)=0.d0
        do i=n1,n2
          dmsd(2,itp)=dmsd(2,itp)+amsd(itp,i+nt)-amsd(itp,i)
        enddo
        dmsd(2,itp)=dmsd(2,itp)/fac/dble(nt+1)/dble(nt)/dble(dt)/dble(inm(itp))/6.d0
        write(imsdcond,'(i3,3g20.8,14h  MSD by slope)')itp, &
              dmsd(2,itp),dmsd(2,itp)/1.d4,dmsd(2,itp)/1.d8
      enddo
      write(65,'(10g20.10)')real(iread-1)*tmsd,((dmsd(i,ity),i=1,2),ity=1,nty)
      write(imsdcond,'(/)')
      call flush(imsdcond)
      call flush(65)

      if(.not.lcond) return

      do kkk=1,iread-1
        read(60,*)dum,((cond_time(ity,jty),jty=ity,nty),ity=1,nty)
        read(61,*)dum,((cond_slope(ity,jty),jty=ity,nty),ity=1,nty)
      enddo
      close(60)
      close(61)

      write(imsdcond,*)'individual conductivity (S/cm)'
      eleconne1=0.d0
      eleconne2=0.d0
      do ity=1,nty
        fac=dble(chgmol(ity))**2*dble(elecon)*3.d0*dble(inm(ity))
        dmsd(1,ity)=dmsd(1,ity)*fac
        dmsd(2,ity)=dmsd(2,ity)*fac
        eleconne1=eleconne1+dmsd(1,ity)
        eleconne2=eleconne2+dmsd(2,ity)
        cond_time(ity,ity)=cond_time(ity,ity)-dmsd(1,ity)
        cond_slope(ity,ity)=cond_slope(ity,ity)-dmsd(2,ity)
        do jty=ity,nty
          if(jty.eq.ity) then
            write(imsdcond,'(2i5,3x,2f12.8,14h S/cm by time )')ity,ity,dmsd(1,ity),cond_time(ity,ity)
            write(imsdcond,'(2i5,3x,2f12.8,14h S/cm by slope)')ity,ity,dmsd(2,ity),cond_slope(ity,ity)
          else
            write(imsdcond,'(2i5,3x,f12.8, 26h             S/cm by time )')ity,jty,cond_time(ity,jty)
            write(imsdcond,'(2i5,3x,f12.8, 26h             S/cm by slope)')ity,jty,cond_slope(ity,jty)
          endif
        enddo
      enddo
      write(imsdcond,'(/)')
      write(imsdcond,'(92hconductivity by Nernst-Einstein relation (see &
              Smith et al. JPCB 110(2006)22773 eq.9, by t) =,f12.8,5h S/cm)') eleconne1
      write(imsdcond,'(93hconductivity by Nernst-Einstein relation (see &
              Smith et al. JPCB 110(2006)22773 eq.9, slope) =,f12.8,5h S/cm)') eleconne2
      call flush(imsdcond)
      write(77,'(15g20.10)')real(iread-1)*tmsd,((dmsd(i,ity),i=1,2),ity=1,nty),eleconne1,eleconne2
      call flush(77)
      write(67,'(25g20.10)')real(iread-1)*tmsd,(dmsd(1,ity),(cond_time(ity,jty),jty=ity,nty),ity=1,nty)
      call flush(67)
      write(68,'(25g20.10)')real(iread-1)*tmsd,(dmsd(2,ity),(cond_slope(ity,jty),jty=ity,nty),ity=1,nty)
      call flush(68)

      allocate (cmsd(nty,0:nty,0:nmsd))
      open(69,file='msd.dat')
      do i=0,nmsd
        read(69,*)dum,(cmsd(ity,0,i),ity=1,nty)
        do ity=1,nty
          cmsd(ity,0,i)=cmsd(ity,0,i)*dble(chgmol(ity))**2*dble(inm(ity))
        enddo
      enddo
      close(69)
      open(55,file='chgmsd_ind.dat')
      do i=0,nmsd
        read(55,*)dum,((cmsd(ity,jty,i),jty=ity,nty),ity=1,nty)
        do ity=1,nty
          do jty=ity+1,nty
            cmsd(jty,ity,i)=cmsd(ity,jty,i)
          enddo
          cmsd(ity,ity,i)=cmsd(ity,ity,i)-cmsd(ity,0,i)
        enddo
      enddo
      close(55)
      open(70,file='chgmsd1.dat')
      do i=0,nmsd
        write(70,'(47g20.10)')real(i)*dt,((cmsd(ity,jty,i),jty=0,nty),ity=1,nty)
      enddo
      close(70)
      deallocate (cmsd)
      return
      end

