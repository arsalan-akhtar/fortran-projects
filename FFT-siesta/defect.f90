!=======================================================================
! 21 August Adding FFT
!=======================================================================


Program Defect
!-----------------------------------------------------------------------
!For FFT
!Modules
!-----------------------------------------------------------------------
 use precision,   only : dp, grid_p
 use parallel,    only : Node, Nodes, ProcessorY
 use sys,         only : die
 use alloc,       only : re_alloc, de_alloc
 use m_fft,       only : fft     ! 3-D fast Fourier transform
 use cellsubs,    only : reclat  ! Finds reciprocal lattice vectors
 use cellsubs,    only : volcel  ! Finds unit cell volume  
 use mesh,         only : nsm
 use units,       only:pi
!-----------------------------------------------------------------------
!FOR iorho
implicit none
character(20)                       ::  fname_sub_neutral,task_sub_neutral
character(20)                       ::  fname_sub_charge ,task_sub_charge
integer                             ::  maxp_sub_neutral, mesh_sub_neutral(3)&
&, nspin_sub_neutral, nsm_sub_neutral
integer                             ::  maxp_sub_charge , mesh_sub_charge(3)&
&, nspin_sub_charge,nsm_sub_charge
real,allocatable,dimension(:,:)     ::  f_sub_neutral 
real,allocatable,dimension(:,:)     ::  f_sub_charge 
integer                             ::  i_x,i_s,i2,i3,ip,is,ind,iv 
double precision                    ::  cell_sub_neutral(3,3)!,dcell_sub(3,3)
double precision                    ::  cell_sub_charge(3,3)!,dcell_sub(3,3)
logical                             ::  found_sub_neutral,found_sub_charge,check
real(kind=8),parameter              ::  bohr_to_ang=0.529177
!-----------------------------------------------------------------------
real(dp)                :: VOLUME,rcell_sub_neutral_md(3),rcell_sub_charge_md(3),epsi
real(kind=8)            :: E_lat,g2,ecut,E_lat_z,deltag,deltax,gn2 
double precision        :: rcell_sub_charge(3,3),rcell_sub_neutral(3,3)
integer                 :: q_x,q_y,q_z,numpoint,h,ign
real(kind=8)          :: gamma2,beta2,x
!-----------------------------------------------------------------------
!FOR FFT
!-----------------------------------------------------------------------
integer                          :: n1,n2,n3,ic
real(kind=8) ,allocatable,dimension(:,:)  ::  f_sub_kind_neutral !(kind=8)
real(kind=8) ,allocatable,dimension(:,:)  ::  f_sub_kind_charge
real(grid_p), pointer :: CG_neutral(:,:)
real(grid_p), pointer :: CG_charge(:,:)
!integer               ::      nsm
!-----------------------------------------------------------------------
real(kind=8) ::q_model,ev_to_k
!-----------------------------------------------------------------------
!write(*,*) "Please Enter the File Name: "   
!read(*,*)  fname_sub  
!fname_sub_neutral='ZrO2-G.VT'
fname_sub_neutral='ZrO2.VT'
fname_sub_charge='ZrO2-2p.VT' ! +2VO charge
task_sub_neutral='read'
task_sub_charge='read'
nsm_sub_neutral=1
nsm_sub_charge=1
ecut=ev_to_k(2400.0_dp)
epsi=20.0
x=0.0_dp
gamma2=1.0_dp
beta2=1.0_dp
!=======================================================================
!READING Neutral VT
call iorho(task_sub_neutral,fname_sub_neutral,cell_sub_neutral,&
& mesh_sub_neutral,nsm_sub_neutral,maxp_sub_neutral,nspin_sub_neutral,f_sub_neutral,found_sub_neutral)
!call iorho(task_sub,fname_sub_neutral,cell_sub,mesh_sub,nsm_sub,maxp_sub,nspin_sub,f_sub,found_sub)!
!READING Charged VT
call iorho(task_sub_charge,fname_sub_charge,cell_sub_charge,&
& mesh_sub_charge,nsm_sub_charge,maxp_sub_charge,nspin_sub_charge,f_sub_charge,found_sub_charge)
! CHECKING FILE EXIST OR NOT
if ( found_sub_neutral ) then
    write (*,*) "Neutral Potential Founded",found_sub_neutral
    else
    write (*,*) "ERROR Couldn't read Neutral Potential File'"
    STOP
endif
if ( found_sub_charge ) then
    write (*,*) "Charged Potential Founded",found_sub_charge
else
    write (*,*) "ERROR Couldn't read Charge Potential File'"
    STOP
endif

write (*,*) "mode of (Neutral) :" , task_sub_neutral
write (*,*) "mode of (Charge) :" , task_sub_neutral

100 format (3(F10.5,5X))
! Neutral READING INFO
write(*,*) "============================================================"
write(*,*) "                      NEUTRAL INFO"
write(*,*) "============================================================"
write (*,*)  " The Neutral Lattice Vector is in Ang"
write (*,100) cell_sub_neutral(1,1)*bohr_to_ang, cell_sub_neutral(1,2)*bohr_to_ang, cell_sub_neutral(1,3)*bohr_to_ang
write (*,100) cell_sub_neutral(2,1)*bohr_to_ang, cell_sub_neutral(2,2)*bohr_to_ang, cell_sub_neutral(2,3)*bohr_to_ang
write (*,100) cell_sub_neutral(3,1)*bohr_to_ang, cell_sub_neutral(3,2)*bohr_to_ang, cell_sub_neutral(3,3)*bohr_to_ang
write (*,*) "the number of spins (nspin_sub) :" , nspin_sub_neutral
write (*,*) "the nubmer of mesh point (mesh_sub):" , mesh_sub_neutral
write (*,*) "the number of sub mesh point (nsm) :" , nsm_sub_neutral
write (*,*) "the number of maximum mesh point (maxp) :" , maxp_sub_neutral
! Charge READING INFO
write(*,*) "============================================================"
write(*,*) "                      CHARGE INFO"
write(*,*) "============================================================"
write (*,*)  " The Charge Lattice Vector is in Ang"
write (*,100) cell_sub_charge(1,1)*bohr_to_ang, cell_sub_charge(1,2)*bohr_to_ang, cell_sub_charge(1,3)*bohr_to_ang
write (*,100) cell_sub_charge(2,1)*bohr_to_ang, cell_sub_charge(2,2)*bohr_to_ang, cell_sub_charge(2,3)*bohr_to_ang
write (*,100) cell_sub_charge(3,1)*bohr_to_ang, cell_sub_charge(3,2)*bohr_to_ang, cell_sub_charge(3,3)*bohr_to_ang
write (*,*) "the number of spins (nspin_sub) :" , nspin_sub_charge
write (*,*) "the nubmer of mesh point (mesh_sub):" , mesh_sub_charge
write (*,*) "the number of sub mesh point (nsm) :" , nsm_sub_charge
write (*,*) "the number of maximum mesh point (maxp) :" , maxp_sub_charge
!-----------------------------------------------------------------------

allocate(f_sub_neutral(maxp_sub_neutral,nspin_sub_neutral))
allocate(f_sub_charge(maxp_sub_charge,nspin_sub_charge))
call iorho(task_sub_neutral,fname_sub_neutral,cell_sub_neutral,mesh_sub_neutral&
& ,nsm_sub_neutral,maxp_sub_neutral,nspin_sub_neutral,f_sub_neutral,found_sub_neutral)!
call iorho(task_sub_charge,fname_sub_charge,cell_sub_charge,&
& mesh_sub_charge,nsm_sub_charge,maxp_sub_charge,nspin_sub_charge,f_sub_charge,found_sub_charge)
!-----------------------------------------------------------------------
!Printing the readed VT
!-----------------------------------------------------------------------
!101 format ("For spin=",2 (5X,I10),E25.10)
do i_s=1,nspin_sub_neutral 
    do i_x=1,maxp_sub_neutral 
!        write(*,101) i_s,i_x, f_sub_neutral(i_x,i_s)
    end do
end do
do i_s=1,nspin_sub_charge 
    do i_x=1,maxp_sub_charge
!        write(*,101) i_s,i_x, f_sub_charge(i_x,i_s)
    end do
end do
!-----------------------------------------------------------------------
!Taking FFT
!-----------------------------------------------------------------------
!allocate(f_sub_kind(maxp_sub,nspin_sub))
!call iorho(task_sub,fname_sub,cell_sub,mesh_sub,nsm_sub,maxp_sub,nspin_sub,f_sub_kind,found_sub)!
!VOLUME = VOLCEL( cell_sub )
!write(*,*)"The Volume is  :",VOLUME ," A^3"
write(*,*) "============================================================"
write(*,*) "                      NEUTRAL FFT                           "
write(*,*) "============================================================"
allocate(f_sub_kind_neutral(maxp_sub_neutral,nspin_sub_neutral))
f_sub_kind_neutral(:,:)=f_sub_neutral(:,:)
nullify( CG_neutral )
call re_alloc( CG_neutral, 1, 2, 1, mesh_sub_neutral(1)*mesh_sub_neutral(2)*mesh_sub_neutral(3))
!call re_alloc( CG_neutral, 1, 2, 1, 80*80*80)
!C     Copy density to complex array
!!$OMP parallel do default(shared), private(I)
write(*,*) "Neutral: REAL MESH of "
do ic=1,maxp_sub_neutral
   CG_neutral(1,ic) = f_sub_kind_neutral(ic,1)
   CG_neutral(1,ic) = f_sub_neutral(ic,1)
   write(*,*)ic,f_sub_neutral(ic,1)    
   CG_neutral(2,ic) = 0.0_grid_p
end do
nsm=nsm_sub_neutral
!!$OMP end parallel do
call fft( CG_neutral, mesh_sub_neutral, -1 )

write(*,*) "FT"
do ic=1,maxp_sub_neutral
    write(*,*) ic,CG_neutral(1,ic)
end do
!write(*,*) "BFT"
!call fft( CG, mesh_sub, +1 )
!do ic=1,maxp_sub
!    write(*,*) ic,CG(1,ic)
!end do
!call fft(f_sub_kind,mesh_sub, -1 )
!101 format ("For spin=",2 (5X,I10),E25.10)
!write(*,*) "4 kind"
!do i_s=1,nspin_sub 
!    do i_x=1,maxp_sub 
!        write(*,101) i_s,i_x, f_sub_kind(i_x,i_s)
!    end do
!end do
!call de_alloc(f_sub_kind )
!call de_alloc(CG)
write(*,*) "============================================================"
write(*,*) "                      Charged FFT                           "
write(*,*) "============================================================"
allocate(f_sub_kind_charge(maxp_sub_charge,nspin_sub_charge))
f_sub_kind_charge(:,:)=f_sub_charge(:,:)
nullify( CG_charge )
call re_alloc( CG_charge, 1, 2, 1, mesh_sub_charge(1)*mesh_sub_charge(2)*mesh_sub_charge(3))
!call re_alloc( CG_neutral, 1, 2, 1, 80*80*80)
!C     Copy density to complex array
!!$OMP parallel do default(shared), private(I)
write(*,*) "Charge: REAL MESH of "
do ic=1,maxp_sub_charge
   CG_charge(1,ic) = f_sub_kind_charge(ic,1)
   CG_charge(1,ic) = f_sub_charge(ic,1)
   write(*,*)ic,f_sub_charge(ic,1)    
   CG_charge(2,ic) = 0.0_grid_p
end do
nsm=nsm_sub_charge
!!$OMP end parallel do
call fft( CG_neutral, mesh_sub_neutral, -1 )

write(*,*) "FT"
do ic=1,maxp_sub_neutral
    write(*,*) ic,CG_neutral(1,ic)
end do
!-----------------------------------------------------------------------
!Calculating V^{lr}
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!Calculating E_lat
!-----------------------------------------------------------------------
!write (*,*)"Charge", q_model(0.0_dp,1.0_dp,1.0_dp,1.0_dp)
call reclat( cell_sub_neutral,rcell_sub_neutral,1 )


write(*,*) ,"The G VECTOR is  :",rcell_sub_neutral(1,1),rcell_sub_neutral(1,2),rcell_sub_neutral(1,3)        
write(*,*) ,"The G VECTOR is  :",rcell_sub_neutral(2,1),rcell_sub_neutral(2,2),rcell_sub_neutral(2,3)         
write(*,*) ,"The G VECTOR is  :",rcell_sub_neutral(3,1),rcell_sub_neutral(3,2),rcell_sub_neutral(3,3)          
!write(*,*) "============================================"
VOLUME = VOLCEL( cell_sub_neutral )
write(*,*)"The Volume is  :",VOLUME ," A^3"
!Finding norm (modules) of Reciprocal Cell-vector 
 do iv = 1,3
   rcell_sub_neutral_md(iv) = dot_product(rcell_sub_neutral(:,iv),rcell_sub_neutral(:,iv))
   rcell_sub_neutral_md(iv) = sqrt(rcell_sub_neutral_md(iv))
 enddo
write(*,*),"The parameter G Vector: ",rcell_sub_neutral_md(1),rcell_sub_neutral_md(2),rcell_sub_neutral_md(3)
!write(*,*),"The Norm of G Vector: ", sqrt(rcell_sub_md(1)**2+rcell_sub_md(2)**2+rcell_sub_md(3)**2)    
write(*,*),"The Norm of G Vector: ",rcell_sub_neutral_md(1),rcell_sub_neutral_md(2),rcell_sub_neutral_md(3)
write(*,*),"ECut k :", ecut!ev_to_k(500.0)
!200 format (4(A3,10X))
!write(*,*),"             G2             G_x                  G_y                  G_z "
!Finding norm (modules) of G-vectors within cutoff 
do q_x=1,int(ecut)
     do q_y=1,int(ecut)
          do q_z=1,int(ecut)
               if ((sqrt((rcell_sub_neutral_md(1)*q_x)**2+(rcell_sub_neutral_md(1)*q_y)**2&
&               +(rcell_sub_neutral_md(1)*q_z )**2)) < ecut ) then
                    g2=(sqrt((rcell_sub_neutral_md(1)*q_x)**2+(rcell_sub_neutral_md(1)*q_y)**2&
&                    +(rcell_sub_neutral_md(1)*q_z )**2))
!                    write(*,*)  g2,rcell_sub_md(1)*q_x , rcell_sub_md(1)*q_y,rcell_sub_md(1)*q_z 
               end if
          enddo
     enddo
end do
!Calculating E_lat
E_lat=0.0
E_lat_z=0.0
deltag=0.0
do q_x=1,int(ecut)
     do q_y=1,int(ecut)
          do q_z=1,int(ecut)
               if ( (sqrt((rcell_sub_neutral_md(1)*q_x)**2+(rcell_sub_neutral_md(1)*q_y)**2 &
&               +(rcell_sub_neutral_md(1)*q_z )**2)) < ecut ) then
                    g2=(sqrt((rcell_sub_neutral_md(1)*q_x)**2+(rcell_sub_neutral_md(1)*q_y)**2 &
&                    +(rcell_sub_neutral_md(1)*q_z )**2))
                    E_lat= E_lat+((q_model(x,gamma2,beta2,g2)**2)/(g2**2))
!                    V_lr(q_)                    
!                    write(*,*) q_model(x,gamma2,beta2,g2),E_lat*((2*pi)/(epsi*VOLUME ))
                    
               end if
          enddo
     enddo
end do
numpoint=1000000
deltag=ecut/numpoint
do ign=0,numpoint
     E_lat_z= E_lat_z + deltag*(q_model(x,gamma2,beta2,deltag)**2)              
     deltag=ign+deltag
!     write(*,*) ign,(E_lat_z*(1/pi*(epsi)))
enddo
write(*,*),"TOTAL E_Lattice", ((E_lat*((2*pi)/(epsi*VOLUME )))-(E_lat_z*(1/pi*(epsi))))



!                    E_lat_z= E_lat_z + deltag*(q_model(x,gamma2,beta2,g2)**2)              
!                    deltag=g2-deltag
!,(E_lat_z*(1/pi*(epsi)))

!write(*,100)
!f=q_model(0.0_dp,1.0_dp,1.0_dp,g2)







deallocate(f_sub_neutral)
deallocate(f_sub_charge)
deallocate(f_sub_kind_neutral)
deallocate(f_sub_kind_charge)
!end subroutine test
End Program Defect



! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
subroutine iorho( task, fname, cell, mesh, nsm, maxp, nspin, f, found )!

! *********************************************************************
! Saves/recovers the electron density at the mesh points.
! This simplified version reads only files written in serial mode.
! Writen by J.Soler July 1997.
! *************************** INPUT **********************************
! character*(*) task      : 'read'/'READ' or 'write'/'WRITE'
! character*(*) fname     : File name for input or output
! integer nsm             : Number of sub-mesh points per mesh point
!                           (not used in this version)
! integer maxp            : First dimension of array rho
! integer nspin           : Second dimension of array rho
! ************************** OUTPUT **********************************
! integer maxp            : Required first dimension of array rho,
!                           equal to mesh(1)*mesh(2)*mesh(3)
!                           Set only when task='read' and required
!                           value is larger than input value
! integer nspin           : Number of spin polarizations (1 or 2)
! logical found           : Were data found? (only when task='read')
! ******************** INPUT or OUTPUT (depending on task) ***********
! real*8  cell(3,3)       : Lattice vectors
! integer mesh(3)         : Number of mesh divisions of each
!                           lattice vector
! real    f(maxp,nspin)   : Electron density
!                           Notice single precision in this version
! *************************** UNITS ***********************************
! Units should be consistent between task='read' and 'write'
! ******************** BEHAVIOUR **************************************
! If task='read', and the values of maxp or nspin on input are less than
! those required to copy the array f from the file, then the required
! values of maxp and nspin are returned on output, but f is not read.
! *********************************************************************
      implicit          none
! Arguments
      character*(*)     fname, task
      integer           maxp, mesh(3), nspin, nsm
      real             f(maxp,nspin)
      double precision  cell(3,3)
      logical           found
! Internal variables and arrays
      character fform*11
      integer   i2, i3, ind, ip, is, np, ns
! Fix whether formatted or unformatted files will be used
      fform = 'unformatted'
! Look for data file
      inquire( file=fname, exist=found )
      if (.not.found) return
      !write (*,*)"Not Found"
      
! Read unit cell vectors, number of mesh points and spin components
      open( unit=1, file=fname, status='old', form=fform )
      if (fform .eq. 'formatted') then
        read(1,*) cell
        read(1,*) mesh, ns
      else
        read(1) cell
        read(1) mesh, ns
      endif
!      write (*,*) "Founded"
! Read density (only if array f is large enough)
      np = mesh(1) * mesh(2) * mesh(3)
      if (ns.gt.nspin .or. np.gt.maxp) then
        maxp = np
      else
        if (fform .eq. 'formatted') then
!          ind = 0  !THIS IS THE BUG FOR READING
          do is = 1,ns
            ind = 0
            do i3 = 1,mesh(3)
              do i2 = 1,mesh(2)
                read(1,*) (f(ind+ip,is),ip=1,mesh(1))
                ind = ind + mesh(1)
              enddo
            enddo
          enddo
        else
!          ind = 0  !THIS IS THE BUG FOR READING
          do is = 1,ns
            ind = 0
            do i3 = 1,mesh(3)
              do i2 = 1,mesh(2)
                read(1) (f(ind+ip,is),ip=1,mesh(1))
                ind = ind + mesh(1)
              enddo
            enddo
          enddo
        endif
      endif
      close(1)
      nspin = ns
      maxp = np
      end subroutine iorho
!  DO LOOP FOR READING
!  Spin(ns)|  z  |  y  | x 
!-------------------------------
!        1 |  1  |  1  |   1 ...  80
!        1 |  1  |  2  |  81 ... 160    
!        1 |  1  |  3  | 161 ... 240
!        1 | ... | ... | 241 ... 320
!        1 | ... | ... | 321 ... 400
!        1 | ... | ... | 401 ... 480
!        1 |  1  |  80 | ..  ...  ..



!-----------------------------------------------------------------------
!Q MODEL Functions
!-----------------------------------------------------------------------
real(kind=8) Function q_model(x,gamma2,beta2,g2)
!     use precision,   only : dp, grid_p
!    "Gaussian Model in PyCDT"
!    "q(r)=q[x exp(-r/gamma) +(1-x) exp(-r^2/beta^2)]"
!    "q(g2)=x/sqrt(1+gamma^2*g2)+*(1-x)*exp^(-beta2*g2) "
!     INPUT          x: Weight
!               gamma2: Exponential Decay Constant
!                beta2: Gaussian Decay    
!                   g2: square of reciprocal Vector  
!    OUTPUT    q_model: Charge densitu at reciprocal magnitude 
     implicit none
     real(kind=8)          :: gamma2,beta2,x
     real(kind=8)       :: g2 
     q_model=x/sqrt(1+beta2*g2)+(1-x)*exp(-0.25*beta2*g2)
end Function q_model 
!-----------------------------------------------------------------------
!eV to k  Functions
!-----------------------------------------------------------------------
real(kind=8) Function ev_to_k(eV)
!    Convert energy to reciprocal vector magnitude k via hbar*k^2/2m
!    Args:
!        a: Energy in eV.
!
!    Returns:   Reciprocal vector magnitude (units of 1/Bohr).
     implicit none
     real(kind=8), parameter :: ang_to_bohr=1.889726878
     real(kind=8) :: eV
     ev_to_k=sqrt(eV/3.80986)*ang_to_bohr   
end Function


!=======================================================================
!allocate(f_normal_sub(mesh_sub(1),mesh_sub(2),mesh_sub(3),nspin_sub))


!ind = 0  I THink bug
!outer_counter=0
!do is = 1,nspin_sub
!    ind = 0  !I THink its correct
!    inner_counter=0  
!    do i3 = 1,mesh_sub(3)
!        do i2 = 1,mesh_sub(2)
!            do ip=1,mesh_sub(1)        
!                f_normal_sub(ip,i2,i3,is)=f_sub(ind+ip,is)
!                inner_counter=inner_counter+1
!                outer_counter=outer_counter+1
!                write(*,*)ip,i2,i3,is,inner_counter,outer_counter,f_sub(inner_counter,is)
!            end do
!        enddo
!    enddo
!enddo


!do i_x=1,mesh_sub(1)
!    write(*,101) i_x, f_normal_sub(i_x,i_x,1)
!end do
!write(*,*) f_normal_sub(1,2,1,1)
