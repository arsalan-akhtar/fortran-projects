subroutine iorho( task, fname, cell, mesh, nsm, maxp, nspin, f, found , fnormal)!

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
      real             f(maxp,nspin),fnormal(mesh(1),mesh(1),mesh(1),nspin) ! new format
!      real       f_normal(mesh(1),mesh(1),mesh(1),nspin) ! new format
!      real             f_normal(80,80,80,2) ! new format
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
      !write (*,*) "Founded"
! Read density (only if array f is large enough)
      np = mesh(1) * mesh(2) * mesh(3)
      if (ns.gt.nspin .or. np.gt.maxp) then
        maxp = np
      else
        if (fform .eq. 'formatted') then
          ind = 0
          do is = 1,ns
            do i3 = 1,mesh(3)
              do i2 = 1,mesh(2)
!                read(1,*) (f(ind+ip,is),ip=1,mesh(1))
!                read(1,*) (f_normal(ip,i2,i3,is),ip=1,mesh(1))
                ind = ind + mesh(1)
              enddo
            enddo
          enddo
        else
          ind = 0
          do is = 1,ns
            do i3 = 1,mesh(3)
              do i2 = 1,mesh(2)
                read(1) (f(ind+ip,is),ip=1,mesh(1))
                read(1) (fnormal(ip,i2,i3,is),ip=1,mesh(1))
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



