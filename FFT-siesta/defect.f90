!=======================================================================
! 21 August Adding FFT
!=======================================================================


Program Defect
!-----------------------------------------------------------------------
!For FFT
!     Modules
 !subroutine test()
 use precision,   only : dp, grid_p
 use parallel,    only : Node, Nodes, ProcessorY
 use sys,         only : die
 !use alloc,       only : re_alloc, de_alloc
 !use m_fft,       only : fft     ! 3-D fast Fourier transform
 !use cellsubs,    only : reclat  ! Finds reciprocal lattice vectors
 use cellsubs,    only : volcel  ! Finds unit cell volume
!-----------------------------------------------------------------------
      



!-----------------------------------------------------------------------
!use iorho
implicit none
character(20)                    ::  fname_sub, task_sub
integer                          ::  maxp_sub, mesh_sub(3), nspin_sub, nsm_sub,i_x,i_s
real ,allocatable,dimension(:,:)  ::  f_sub !(kind=8)
real ,allocatable,dimension(:,:,:,:)  ::  f_normal_sub
integer                          ::  i2,i3,ip,is,ind,inner_counter,outer_counter
double precision                 ::  cell_sub(3,3)!,dcell_sub(3,3)
logical                          ::  found_sub,check
real,parameter                   ::  bohr_to_ang=0.529177
!-----------------------------------------------------------------------
real(dp)        :: VOLUME   
!-----------------------------------------------------------------------

!write(*,*) "Please Enter the File Name: "   
!read(*,*)  fname_sub  
fname_sub='ZrO2-G.VT'
task_sub='read'
nsm_sub=1
!=.TRUE.
call iorho(task_sub,fname_sub,cell_sub,mesh_sub,nsm_sub,maxp_sub,nspin_sub,f_sub,found_sub)!

if ( found_sub ) then
    write (*,*) found_sub
endif

write (*,*) "mode (task_sub) :" , task_sub


100 format (3(F10.5,5X))
write (*,*)  " The Lattice Vector is in Ang"
write (*,100) cell_sub(1,1)*bohr_to_ang, cell_sub(1,2)*bohr_to_ang, cell_sub(1,3)*bohr_to_ang
write (*,100) cell_sub(2,1)*bohr_to_ang, cell_sub(2,2)*bohr_to_ang, cell_sub(2,3)*bohr_to_ang
write (*,100) cell_sub(3,1)*bohr_to_ang, cell_sub(3,2)*bohr_to_ang, cell_sub(3,3)*bohr_to_ang
write (*,*)  " The Lattice Vector is in Bohr"
write (*,100) cell_sub(1,1), cell_sub(1,2), cell_sub(1,3)
write (*,100) cell_sub(2,1), cell_sub(2,2), cell_sub(2,3)
write (*,100) cell_sub(3,1), cell_sub(3,2), cell_sub(3,3)

write (*,*) "the number of spins (nspin_sub) :" , nspin_sub
write (*,*) "the nubmer of mesh point (mesh_sub):" , mesh_sub
!write (*,*) "the nubmer of mesh point (mesh_sub):" , mesh_sub(10)
write (*,*) "the number of sub mesh point (nsm) :" , nsm_sub
write (*,*) "the number of maximum mesh point (maxp) :" , maxp_sub

allocate(f_sub(maxp_sub,nspin_sub))
call iorho(task_sub,fname_sub,cell_sub,mesh_sub,nsm_sub,maxp_sub,nspin_sub,f_sub,found_sub)!

101 format ("For spin=",2 (5X,I10),E25.10)
do i_s=1,nspin_sub 
    do i_x=1,maxp_sub 
        !write(*,101) i_s,i_x, f_sub(i_x,i_s)
    end do
end do
!-----------------------------------------------------------------------
!Taking FFT
VOLUME = VOLCEL( cell_sub )
write(*,*)"The Volume is :",VOLUME

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



deallocate(f_sub)
!deallocate(f_normal_sub)
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


