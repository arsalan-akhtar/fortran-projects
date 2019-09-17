!=======================================================================
! 21 August Adding FFT
!=======================================================================
!VO Position  = (0.25,0.25,0.25) 

Program Defect
!-----------------------------------------------------------------------
! Loading Modules
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
implicit none
!-----------------------------------------------------------------------
! Declaration of iorho Neutral Reading
!-----------------------------------------------------------------------
character(20)  ::  fname_sub_neutral,task_sub_neutral
integer        ::  maxp_sub_neutral,mesh_sub_neutral(3)&
&,nspin_sub_neutral, nsm_sub_neutral
real,allocatable,dimension(:,:)   ::  f_sub_neutral,f_sub_charge 
double precision                  ::  cell_sub_neutral(3,3)
logical                           ::  found_sub_neutral
real(dp)       :: volume_neutral,rcell_sub_neutral_md(3)
!-----------------------------------------------------------------------
!Declaration of iorho Charge Reading
!-----------------------------------------------------------------------
character(20)  ::  fname_sub_charge ,task_sub_charge
integer        ::  maxp_sub_charge , mesh_sub_charge(3)&
&, nspin_sub_charge,nsm_sub_charge
double precision  ::  cell_sub_charge(3,3)
logical           ::  found_sub_charge,check
real(dp)          :: volume_charge,rcell_sub_charge_md(3)
!-----------------------------------------------------------------------
!Declaration For Iterations
!-----------------------------------------------------------------------
integer  ::  iv,ns,n1,n2,n3,ic

!-----------------------------------------------------------------------
!Declaration Of Parameters
!-----------------------------------------------------------------------
real(kind=8),parameter  :: bohr_to_ang=0.529177
real(kind=8), parameter :: ang_to_bohr=1.889726878
real(kind=8),parameter  :: ry_to_ev=13.6056980659
!-----------------------------------------------------------------------
real(kind=8)            :: E_q_lat,E_q_lat_z,ecut,deltax,gn2 
double precision        :: rcell_sub_charge(3,3),rcell_sub_neutral(3,3)!,r(3)
integer                 :: q_x,q_y,q_z,numpoint,h,ign

!-----------------------------------------------------------------------
!Declaration For E_lattice
!-----------------------------------------------------------------------
real(kind=4)            :: epsi
integer                 :: charge
real(kind=8)            :: E_r_lat,delta
!-----------------------------------------------------------------------
!Declaration For V_lr
!-----------------------------------------------------------------------
real(kind=4),allocatable,dimension(:,:)     ::  V_lr 
!real(grid_p), pointer                       :: V_lr(:,:)
!-----------------------------------------------------------------------
!FOR FFT
!-----------------------------------------------------------------------
real(kind=8) ,allocatable,dimension(:,:)  ::  f_sub_kind_neutral
real(kind=8) ,allocatable,dimension(:,:)  ::  f_sub_kind_charge
real(grid_p), pointer :: CG_neutral(:,:)
real(grid_p), pointer :: CG_charge(:,:)
real(grid_p), pointer :: V_lr_p(:,:),V_sr(:,:),deltaV(:,:)
!integer               ::      nsm
!-----------------------------------------------------------------------
!Declaration For Functions
!-----------------------------------------------------------------------
real(kind=8)   :: ev_to_k,q_model,d_position(3)
real(kind=4)   :: q_model_r,r_norm,gamma_d,beta_d,x_d,n
real(kind=8)   :: gamma2,beta2,x,deltag,g2
integer(kind=8) :: Vl_d,i_x,i_s,i2,i3,ip,is,ind,counter,r_0(3),r_x,r_y,r_z
!=======================================================================
! INPUTS & User Defined Values
!=======================================================================
!write(*,*) "Please Enter the File Name: "   
!read(*,*)  fname_sub  
!fname_sub_neutral='ZrO2-G.VT'
fname_sub_neutral='ZrO2.VT'
fname_sub_charge='ZrO2-2p.VT' ! +2VO charge
task_sub_neutral='read'
task_sub_charge='read'
nsm_sub_neutral=1
nsm_sub_charge=1
ecut=ev_to_k(700000.0_dp)
epsi=20.0
x_d=0.0
x=0.0
gamma_d=0.5!1.0
gamma2=0.5!1.0
beta_d=1.0 !*(10**(-5))!*ang_to_bohr
beta2=1.0!*(10**(-5))
d_position(1)=0.25
d_position(2)=0.25
d_position(3)=0.25
r_0(1)=1  ! going to use internaly
r_0(2)=1
r_0(3)=1
ns=1
charge=-2

!=======================================================================
!-----------------------------------------------------------------------
!READING Neutral VT
!-----------------------------------------------------------------------
call iorho(task_sub_neutral,fname_sub_neutral,cell_sub_neutral,&
& mesh_sub_neutral,nsm_sub_neutral,maxp_sub_neutral,nspin_sub_neutral,f_sub_neutral,found_sub_neutral)
!-----------------------------------------------------------------------
!READING Charged VT
call iorho(task_sub_charge,fname_sub_charge,cell_sub_charge,&
& mesh_sub_charge,nsm_sub_charge,maxp_sub_charge,nspin_sub_charge,f_sub_charge,found_sub_charge)
!-----------------------------------------------------------------------
! CHECKING FILE EXIST OR NOT
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
! Neutral READING INFO
100 format (3(F10.5,5X))
!-----------------------------------------------------------------------
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
volume_neutral = VOLCEL( cell_sub_neutral )
write(*,*)"The Volume is  :",volume_neutral ," A^3"

!-----------------------------------------------------------------------
! Charge READING INFO
!-----------------------------------------------------------------------
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
volume_charge = VOLCEL( cell_sub_charge )
write(*,*)"The Volume is  :",volume_charge ," A^3"
!-----------------------------------------------------------------------
! Finding Location of Defect
!-----------------------------------------------------------------------
write(*,*) "============================================================"
write(*,*) "                      Position of  Defect                   "
write(*,*) "============================================================"
write(*,*)int(mesh_sub_charge(1)*d_position(1))
write(*,*)int(mesh_sub_charge(2)*d_position(2))
write(*,*)int(mesh_sub_charge(3)*d_position(3))
!-----------------------------------------------------------------------
! Allocating For V_neutral,V_charge
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
!   write(*,*)ic,f_sub_neutral(ic,1)    
   CG_neutral(2,ic) = 0.0_grid_p
end do
nsm=nsm_sub_neutral
!!$OMP end parallel do
call fft( CG_neutral, mesh_sub_neutral, -1 )

write(*,*) "FT"
do ic=1,maxp_sub_neutral
!    write(*,*) ic,CG_neutral(1,ic)
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
!   write(*,*)ic,f_sub_charge(ic,1)    
   CG_charge(2,ic) = 0.0_grid_p
end do
nsm=nsm_sub_charge
!!$OMP end parallel do
call fft( CG_neutral, mesh_sub_neutral, -1 )
write(*,*) "FT"
do ic=1,maxp_sub_neutral
!    write(*,*) ic,CG_neutral(1,ic)
end do
write(*,*) "============================================================"
!deallocate(f_sub_kind_charge)
!-----------------------------------------------------------------------
!Calculating V^{lr}
!-----------------------------------------------------------------------
write(*,*) "============================================================"
write(*,*) "Calculating V^{lr}"
write(*,*) "============================================================"
!call de_alloc(V_lr)
!Vl_d=(1600)!maxp_sub_charge!(160*160*160*160)
!deallocate(V_lr)
write(*,*) maxp_sub_neutral,nspin_sub_charge ! For testing *************
allocate(V_lr(maxp_sub_neutral,nspin_sub_charge))
!nullify( V_lr )
!call re_alloc( V_lr, 1, 2, 1, mesh_sub_charge(1)*mesh_sub_charge(2)*mesh_sub_charge(3))
!call re_alloc( V_lr, 1, 2, 1, 160*160*160*160)
!deallocate(V_lr)
!V_lr_d=160*160*160
!allocate(V_lr(4096000,2)) Worked!!!!!!!!
counter=1
do is = 1,nspin_sub_charge
     ind = 0
     do i3 = 1,mesh_sub_charge(3)
          do i2 = 1,mesh_sub_charge(2)
               do ip=1,mesh_sub_charge(1)
                    r_x=abs(ip-r_0(1))
                    r_y=abs(i2-r_0(2))
                    r_z=abs(i3-r_0(3))
                    r_norm = (r_x**2+r_y**2+r_z**2)**(1.0/2.0)
!                    r_norm = (r_x+r_y+r_z)
                    !V_lr(ind+ip,is)=0.0  ! For testing *************
                    if ((i3 .EQ. r_0(3)) .AND. (i2.EQ.r_0(2)) .AND. (ip.EQ.r_0(1)))then
                         V_lr(ind+ip,is)=0.0
                    else 
					V_lr(ind+ip,is)= (1/epsi)*q_model_r(charge,x_d,gamma_d,beta_d,ip,i2,i3,r_0)/r_norm
                    endif                     
!                    write (*,*) ind+ip,ip,i2,i3,is, V_lr(ind+ip,is)
                    counter=counter+1
               end do
               ind = ind + mesh_sub_charge(1)
          end do
     end do    
end do
!-----------------------------------------------------------------------
!Calculating V^{lr}_periodic  
!-----------------------------------------------------------------------
nullify( V_lr_p )
call re_alloc( V_lr_p, 1, 2, 1, mesh_sub_charge(1)*mesh_sub_charge(2)*mesh_sub_charge(3))
!write(*,*) "Test"
do ic=1,maxp_sub_charge
   V_lr_p(1,ic) = V_lr(ic,1)
   V_lr_p(1,ic) = V_lr(ic,1)
!   write(*,*)ic,V_lr(ic,1)    !For testing
   V_lr_p(2,ic) = 0.0_grid_p
end do
write(*,*) "============================================================"
write(*,*) "Calculating FFT of V^{lr}                                   "
write(*,*) "============================================================"
call fft( V_lr_p, mesh_sub_charge, -1 )
do ic=1,maxp_sub_charge
!    write(*,*) ic,V_lr_p(1,ic)    
end do
call reclat( cell_sub_charge,rcell_sub_charge,1 )
!-----------------------------------------------------------------------
!Finding norm (modules) of Reciprocal Cell-vector 
!-----------------------------------------------------------------------
do iv = 1,3
   rcell_sub_charge_md(iv) = dot_product(rcell_sub_charge(:,iv),rcell_sub_charge(:,iv))
   rcell_sub_charge_md(iv) = sqrt(rcell_sub_charge_md(iv))
enddo
counter=0
write(*,*) "============================================================"
write(*,*) "Calculating V^{lr}_{periodic} Divide by |G|^2"
write(*,*) "============================================================"
do q_x=1,mesh_sub_charge(1)
     do q_y=1,mesh_sub_charge(2)
          do q_z=1,mesh_sub_charge(3)
               g2=(sqrt((rcell_sub_charge_md(1)*q_x)**2+(rcell_sub_charge_md(2)*q_y)**2&
&                    +(rcell_sub_charge_md(1)*q_z )**2))              
               counter=counter+1
               V_lr_p(1,counter)=(4*pi*V_lr_p(1,counter))/((epsi)*g2**2)
!               write(*,*) counter,g2,V_lr_p(1,counter)
          enddo
     enddo
end do
!-----------------------------------------------------------------------
!Back FFT of  V_lr_p
!-----------------------------------------------------------------------
write(*,*) "============================================================"
write(*,*) "Back FFT of V^{lr}_{periodic}"
write(*,*) "============================================================"
call fft( V_lr_p, mesh_sub_charge, +1 )
do ic=1,maxp_sub_charge
!    write(*,*) ic,V_lr_p(1,ic)    
end do
!-----------------------------------------------------------------------
!Fixing the dimension of V_lr_p
!-----------------------------------------------------------------------
counter=0
do is = 1,nspin_sub_charge
     ind = 0
     counter=0
     do q_x=1,mesh_sub_charge(3)
          do q_y=1,mesh_sub_charge(2)
               do q_z=1,mesh_sub_charge(1)
                    counter=counter+1
                    V_lr_p(ind+q_z,is)=V_lr_p(1,counter) 
!               write(*,*) counter,V_lr_p(1,counter)
!               write(*,*) counter,V_lr_p(ind+q_z,is)               
               enddo
               ind = ind + mesh_sub_charge(1)
          enddo
     enddo
end do
!call de_alloc(V_lr_p)
!-----------------------------------------------------------------------
!Calculating G Vector
!-----------------------------------------------------------------------
!write (*,*)"Charge", q_model(0.0_dp,1.0_dp,1.0_dp,1.0_dp)
call reclat( cell_sub_neutral,rcell_sub_neutral,1 )
write(*,*) "============================================================"
write(*,*) "The G VECTOR is  :"
write(*,100) rcell_sub_neutral(1,1),rcell_sub_neutral(1,2),rcell_sub_neutral(1,3)        
write(*,100) rcell_sub_neutral(2,1),rcell_sub_neutral(2,2),rcell_sub_neutral(2,3)         
write(*,100) rcell_sub_neutral(3,1),rcell_sub_neutral(3,2),rcell_sub_neutral(3,3)          
volume_neutral = VOLCEL( cell_sub_neutral )
write(*,*)"The Volume is  :",volume_neutral ," A^3"
write(*,*) "============================================================"
!-----------------------------------------------------------------------
!Finding norm (modules) of Reciprocal Cell-vector 
!-----------------------------------------------------------------------
 do iv = 1,3
   rcell_sub_neutral_md(iv) = dot_product(rcell_sub_neutral(:,iv),rcell_sub_neutral(:,iv))
   rcell_sub_neutral_md(iv) = sqrt(rcell_sub_neutral_md(iv))
 enddo
write(*,*) "The parameter G Vector: ",rcell_sub_neutral_md(1),rcell_sub_neutral_md(2),rcell_sub_neutral_md(3)
!write(*,*),"The Norm of G Vector: ", sqrt(rcell_sub_md(1)**2+rcell_sub_md(2)**2+rcell_sub_md(3)**2)    
write(*,*) "ECut k :", ecut
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
!                    write(*,*)  g2,rcell_sub_neutral_md(1)*q_x ,&
!&                    rcell_sub_neutral_md(1)*q_y,rcell_sub_neutral_md(1)*q_z 
               end if
          enddo
     enddo
end do
!-----------------------------------------------------------------------
!Calculating E_lat(q) 
!-----------------------------------------------------------------------
write(*,*) "============================================================"
write(*,*) "Calculating  E_lattice(q)"
write(*,*) "============================================================"
E_q_lat=0.0
E_q_lat_z=0.0
deltag=0.0
do q_x=1,int(ecut)
     do q_y=1,int(ecut)
          do q_z=1,int(ecut)
               if ( (sqrt((rcell_sub_neutral_md(1)*q_x)**2+(rcell_sub_neutral_md(1)*q_y)**2 &
&               +(rcell_sub_neutral_md(1)*q_z )**2)) < ecut ) then
                    g2=(sqrt((rcell_sub_neutral_md(1)*q_x)**2+(rcell_sub_neutral_md(1)*q_y)**2 &
&                    +(rcell_sub_neutral_md(1)*q_z )**2))
                    E_q_lat= E_q_lat+((q_model(x,gamma2,beta2,g2)**2)/(g2**2))
!                    write(*,*) q_model(x,gamma2,beta2,g2)!,E_q_lat*((2*pi)/(epsi*volume_neutral ))
!                    write(*,*) E_q_lat !*((2*pi)/(epsi*volume_neutral ))
!                    write(*,*) q_model(x,gamma2,beta2,g2)**2                    
               end if
          enddo
     enddo
end do
!write(*,*) E_q_lat*((2*pi)/(epsi*volume_neutral))
numpoint=1000
deltag=ecut/numpoint
write(*,*) "E_cut is ",ecut, "Delta g is",deltag
do ign=0,numpoint
     E_q_lat_z= E_q_lat_z + (q_model(x,gamma2,beta2,deltag)**2)              
!     write(*,*) ign,deltag,q_model(x,gamma2,beta2,deltag)**2,(E_q_lat_z)
     deltag=ign+deltag
!     write(*,*) ign,deltag,(E_q_lat_z*(1/pi*(epsi)))
enddo
!write(*,*) "Total E_Lattice(g) =", ((E_q_lat*((2*pi)/(epsi*volume_neutral ))))
write(*,*) "Total E_Lattice(g) with g=0 =", E_q_lat*((2*pi)/(epsi*volume_neutral))!*ang_to_bohr**3
write(*,*) "Total Removing g=0 with it self ",E_q_lat_z*(1.0/pi*(epsi))! *deltag
write(*,*) "Total E_Lattice(g)-without g=0 =", (E_q_lat*((2*pi)/(epsi*volume_neutral))- E_q_lat_z*(1/pi*(epsi)))*ry_to_ev!*deltag
!write(*,*) "Total E_Lattice(g) =", ((E_q_lat*((2*pi)/(epsi*volume_neutral )))-(E_q_lat_z*(1/pi*(epsi))))
!                    E_lat_z= E_lat_z + deltag*(q_model(x,gamma2,beta2,g2)**2)              
!                    deltag=g2-deltag
!,(E_lat_z*(1/pi*(epsi)))
!-----------------------------------------------------------------------
!Calculating E_lat(r) 
!-----------------------------------------------------------------------
write(*,*) "============================================================"
write(*,*) "Calculating  E_lattice(r)"
write(*,*) "============================================================"
ind = 0
E_r_lat=0.0
deltax=(volume_neutral/(maxp_sub_charge))
n=(-1)*(charge/volume_neutral)
!write(*,*) volume_charge
do i3 = 1,mesh_sub_charge(3)
     do i2 = 1,mesh_sub_charge(2)
          do ip=1,mesh_sub_charge(1)
               E_r_lat=E_r_lat+(0.5*(q_model_r(charge,x_d,gamma_d,beta_d,ip,i2,i3,r_0) &
& + (n*(V_lr_p(ip+ind,1)-V_lr(ip+ind,1)))) &
& + (n*V_lr(ip+ind,1)))               
!         write(*,*) E_r_lat
          end do
          ind = ind + mesh_sub_charge(1)
     end do
end do    
write(*,*) "Total E_lattice (r) =", E_r_lat*deltax*ry_to_ev !*(ang_to_bohr**3)
!-----------------------------------------------------------------------
!Calculating Delta 
!-----------------------------------------------------------------------
write(*,*) "============================================================"
write(*,*) "Calculating  Delta V"
write(*,*) "============================================================"
allocate(V_sr(maxp_sub_neutral,nspin_sub_charge))
allocate(deltaV(maxp_sub_neutral,nspin_sub_charge))
ind=0
do i3 = 1,mesh_sub_charge(3)
     do i2 = 1,mesh_sub_charge(2)
          do ip=1,mesh_sub_charge(1)              
               deltaV(ind+ip,1)=-f_sub_kind_charge(ind+ip,is)+f_sub_kind_charge(ind+ip,is)+V_lr_p(ind+ip,is)
               V_sr(ind+ip,1)=f_sub_kind_charge(ind+ip,is)-f_sub_kind_charge(ind+ip,is)-V_lr_p(ind+ip,is)
!               write(*,*) "Counter = ",ind+ip,"V_sr_not_alligned= ",V_sr(ind+ip,1),"DeltaV = ",deltaV(ind+ip,is)
               V_sr(ind+ip,1)=f_sub_kind_charge(ind+ip,is)-f_sub_kind_charge(ind+ip,is)-V_lr_p(ind+ip,is)+deltaV(ind+ip,is)
!               write(*,*) "Counter = ",ind+ip,"V_sr_alligned= ",V_sr(ind+ip,1)
          end do
          ind = ind + mesh_sub_charge(1)
     end do
end do   
!-----------------------------------------------------------------------
!Calculating Delta 
!-----------------------------------------------------------------------
write(*,*) "============================================================"
write(*,*) "Calculating  Delta"
write(*,*) "============================================================"
delta=0.0
ind=0
do i3 = 1,mesh_sub_charge(3)
     do i2 = 1,mesh_sub_charge(2)
          do ip=1,mesh_sub_charge(1)              
               delta=delta+ V_sr(ind+ip,1)
!               write(*,*) "Counter = ",ind+ip,"V_sr_alligned= ",V_sr(ind+ip,1)
          end do
          ind = ind + mesh_sub_charge(1)
     end do
end do   
write(*,*) "Delta = ",delta*(1/volume_neutral)

deallocate(f_sub_neutral)
deallocate(f_sub_charge)
call de_alloc(V_lr_p)
deallocate(f_sub_kind_neutral)
deallocate(f_sub_kind_charge)
deallocate(V_lr)
deallocate(V_sr)
deallocate(deltaV)
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
      real              f(maxp,nspin)
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
!Q MODEL(g2) Functions
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
!Q MODEL(r-r_0) Functions
!-----------------------------------------------------------------------
real(kind=4) Function q_model_r(charge_d,x,gamma_d,beta_d,r1,r2,r3,r_0)
     use units,   only : pi
!    "Gaussian Model in PyCDT"
!    "q(r)=q[x exp(-r/gamma) +(1-x) exp(-r^2/beta^2)]"
!    "q(g2)=x/sqrt(1+gamma^2*g2)+*(1-x)*exp^(-beta2*g2) "
!     INPUT          x: Weight
!               gamma_d: Exponential Decay Constant
!                beta_d: Gaussian Decay    
!                   r: real Vector  
!                   r_0: center of defect real Vector  
!    OUTPUT    q_model: Charge density at real magnitude 
     implicit none
     real(kind=4)          :: gamma_d,beta_d,x,N_gamma,N_beta
     real(kind=4)          :: a,b,c
     integer(kind=8)       :: r_0(3),r1,r2,r3
     integer               :: charge_d 
     N_gamma=8*pi*gamma_d**3
     N_beta=pi**(3.0/2.0)*beta_d**3
     a=abs(r1-r_0(1))
     b=abs(r2-r_0(2))
     c=abs(r3-r_0(3))
     q_model_r=charge_d*x*N_gamma**(-1)*exp((-1)*sqrt(a**2+b**2+c**2)/beta_d)+&
&    charge_d*(1-x)*N_beta**(-1)*exp((-1)*(sqrt(a**2+b**2+c**2))**2/beta_d**2)
     !q_model_r=(exp((-1)*a)**2)/beta_d**2
end Function q_model_r



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
!integer Function defect_position(d_p,mesh_points)
!	implicit none
!	integer      :: mesh_points(3),d_p(3),defect_position(3)
	
!	defect_position(1)=int(mesh_points(1)/d_p(1))
!	defect_position(2)=int(mesh_points(2)/d_p(2))
!	defect_position(3)=int(mesh_points(3)/d_p(3))

!end Function defect_position
