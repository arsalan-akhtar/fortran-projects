! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      MODULE m_siesta_init
      private
      public :: siesta_init

      CONTAINS

      subroutine siesta_init()
      use Kpoint_grid,        only: setup_Kpoint_grid, gamma_scf, nkpnt
      USE Kpoint_pdos,        only: gamma_pdos
      use Band,               only: gamma_bands, setup_bands
      use m_ksvinit,          only: gamma_polarization,
     &                              estimate_pol_kpoints
      use m_struct_init,      only: struct_init
      use siesta_options
      use timer_options,      only: use_tree_timer, use_parallel_timer
      use sparse_matrices
      use class_Fstack_Pair_Geometry_dSpData2D, only: new, max_size
      use siesta_geom
      use atomlist,           only: no_u, rmaxkb, amass, lasto, qtot,
     &                              iza, rmaxo, zvaltot, superc,
     &                              initatomlists, no_l, rmaxldau
      use fdf
      use sys,                only: die, bye
#ifdef BSC_CELLXC
      use bsc_xcmod,          only: setXC
#endif /* BSC_CELLXC */
      use molecularmechanics, only: inittwobody
      use metaforce,          only: initmeta
      use m_mpi_utils,        only: broadcast
      use alloc,              only: re_alloc, alloc_report
      use parallelsubs,       only: getnodeorbs
      use m_iostruct,         only: write_struct, read_struct
      use zmatrix,            only: lUseZmatrix
      use zmatrix,            only: write_canonical_ucell_and_Zmatrix
      use m_supercell,        only: exact_sc_ag
      use files,              only: slabel
      use siesta_cmlsubs,     only: siesta_cml_init
      use m_timestamp,        only: timestamp
      use m_wallclock,        only: wallclock
      use parallel,           only: Node, Nodes, PEXSINodes, IOnode
      use parallel,           only: SIESTA_worker
      use parallel,           only: SIESTA_Group, SIESTA_Comm
      use densematrix,        only: Haux, Saux, Psi
      use m_energies
      use m_steps
      use m_spin,             only: init_spin, spin
      use m_spin,             only: qs, init_spiral
      use m_rmaxh
      use m_forces 
      use m_eo
      use m_gamma
      use m_fixed,            only: init_fixed, print_fixed
      use m_ioxv,             only: xv_file_read
      use m_projected_DOS,    only: init_projected_DOS
      use writewave,          only: gamma_wavefunctions,
     &                              setup_wf_kpoints
      use m_new_dm,           only: get_allowed_history_depth
      use m_object_debug,     only: set_object_debug_level
#ifdef DEBUG_XC
      USE siestaXC, only: setDebugOutputUnit
#endif /* DEBUG_XC */
      USE m_timer, only: timer_report ! Sets options for report of CPU times
      use m_check_walltime
#ifdef MPI
      use mpi_siesta
#endif
      use m_cite, only: add_citation, init_citation

#ifdef TRANSIESTA
      use m_ts_init, only : ts_init
#endif

#ifdef SIESTA__FLOOK
      use siesta_dicts, only : dict_populate
      use flook_siesta, only : slua_init, slua_call, LUA_INITIALIZE
#endif

#ifdef NCDF_4
      use netcdf_ncdf
#endif
#ifdef _OPENMP
      use omp_lib, only : omp_get_num_threads, omp_get_schedule
      use omp_lib, only : OMP_SCHED_STATIC, OMP_SCHED_DYNAMIC
      use omp_lib, only : OMP_SCHED_GUIDED, OMP_SCHED_AUTO
#else
!$    use omp_lib, only : omp_get_num_threads, omp_get_schedule
!$    use omp_lib, only : OMP_SCHED_STATIC, OMP_SCHED_DYNAMIC
!$    use omp_lib, only : OMP_SCHED_GUIDED, OMP_SCHED_AUTO
#endif

      implicit none

#ifdef MPI
      logical :: initialized
      integer :: MPIerror  ! Return error code in MPI routines
#endif
      real(dp):: veclen       ! Length of a unit-cell vector
      integer :: i, is, ispin, n_dm_items
      integer :: neigmin      ! Min. number of eigenstates (per k point)
      integer :: ns           ! Number of species
      logical :: user_wants_supercell
      real(dp) :: walltime_m, walltime_s

      integer :: World_Group
      integer :: npSIESTA
      external :: read_options
      external :: reset_messages_file

!---------------------------------------------------------------------- BEGIN
! Initialise MPI unless siesta is running as a subroutine 
! of a driver program that may have initialized MPI already

#ifdef MPI
C     Initialise MPI and set processor number

      call MPI_Initialized( initialized, MPIerror )
      if (.not.initialized) then
#ifdef _OPENMP
         call MPI_Init_Thread(MPI_Thread_Funneled, i, MPIerror)
         if ( MPI_Thread_Funneled /= i ) then
            ! the requested threading level cannot be asserted
            ! Notify the user
           write(*,'(a)') '!!! Could not assert funneled threads'
         end if
#else
         call MPI_Init( MPIerror )
#endif
#ifdef _TRACE_
        call MPItrace_shutdown( )
#endif
      endif ! (.not.initialized)
      
      call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
      call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )
      
      ! Keeper for PEXSI-runs..
      PEXSINodes = Nodes

      call reset_messages_file()
      

      if (.not. fdf_parallel()) then
        call die('siesta_init: ERROR: FDF module ' //
     &           'doesn''t have parallel support')
      endif
#else
      call reset_messages_file()
#endif
      if (Node==0) then
         call system("/bin/rm -f 0_NORMAL_EXIT")
      endif

#ifdef DEBUG_XC
! Set output file unit for debug info
      call setDebugOutputUnit()
#endif /* DEBUG_XC */

!     Initialise input/output...........................................
      call reinit( sname )

#ifdef MPI
      ! Optionally, restrict the number of processors that SIESTA uses, out
      ! of the total pool
      npSIESTA  = fdf_get("MPI.Nprocs.SIESTA",Nodes)
      call MPI_Comm_Group(MPI_Comm_World, World_Group, MPIerror)
      call MPI_Group_incl(World_Group, npSIESTA,
     $                    (/ (i,i=0,npSIESTA-1) /),
     $                    SIESTA_Group, MPIerror)
      call MPI_Comm_create(MPI_Comm_World, SIESTA_Group,
     $                     SIESTA_Comm, MPIerror)

      SIESTA_worker = (Node < npSIESTA)

      ! Swap communicator
      ! This is needed since MPI_COMM_WORLD is used implicitly by
      ! Siesta for all operations. 
      MPI_COMM_WORLD = SIESTA_Comm
      if (SIESTA_worker) then
         call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
         call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )
      endif
#else
      SIESTA_worker = .true.
      Node = 0
      Nodes = 1
      PEXSINodes = 1
#endif
      IOnode = SIESTA_worker .and. (Node .eq. 0)

      ! Be sure to initialize the spin-configuration data
      ! for all processors.

      call init_spin( )

      if (.not. SIESTA_worker) RETURN

#ifdef DEBUG
!     Generates a debug file for every process to track the execution.
!     The file is called debug.$(PID) Where PID is the process number.
!     It also works in secuencial mode.
      call debugMpiOn( )
#endif
#ifdef NCDF_4
      call ncdf_IONode( IONode )
#endif

!     Print version information
      if (IOnode) then
        call prversion()
#ifdef MPI
        if ( PEXSINodes /= Nodes ) then
           write(6,'(/,a,2(i0,a))')
     &          '* Running on ', Nodes, ' SIESTA-nodes and ',
     &          PEXSINodes, ' PEXSI-nodes in parallel'
        else if ( Nodes > 1 ) then           
           write(6,'(/,a,i0,a)')
     &          '* Running on ', Nodes, ' nodes in parallel'
        else
           write(6,'(/,a)') '* Running in serial mode with MPI'
        endif
#else
        write(6,'(/,a)') '* Running in serial mode'
#endif
!$OMP parallel default(shared)
!$OMP master
!$    i = omp_get_num_threads()
!$    write(*,'(a,i0,a)') '* Running ',i,' OpenMP threads.'
!$    write(*,'(a,i0,a)') '* Running ',Nodes*i,' processes.'
#ifdef _OPENMP
!$    write(*,'(a,i0)') '* OpenMP version ', _OPENMP
#endif
!$    call omp_get_schedule(i,is)
!$    select case ( i )
!$    case ( OMP_SCHED_STATIC ) 
!$    write(*,'(a,i0)') '* OpenMP runtime schedule STATIC, chunks ',is
!$    case ( OMP_SCHED_DYNAMIC ) 
!$    write(*,'(a,i0)') '* OpenMP runtime schedule DYNAMIC, chunks ',is
!$    if ( is == 1 ) then
!$     ! this is the default scheduling, probably the user
!$     ! have not set the value, predefine it to 32
!$     is = 32
!$     write(*,'(a,i0)')'** OpenMP runtime schedule DYNAMIC, chunks ',is
!$    end if
!$    case ( OMP_SCHED_GUIDED ) 
!$    write(*,'(a,i0)') '* OpenMP runtime schedule GUIDED, chunks ',is
!$    case ( OMP_SCHED_AUTO ) 
!$    write(*,'(a,i0)') '* OpenMP runtime schedule AUTO, chunks ',is
!$    case default
!$    write(*,'(a,i0)') '* OpenMP runtime schedule UNKNOWN, chunks ',is
!$    end select
!$OMP end master
!$OMP end parallel
!$    call omp_set_schedule(i,is)
         call timestamp('Start of run')
         call wallclock('Start of run')
      endif

!
!     Nullify arrays
      nullify(Haux,Saux,psi)
      nullify(Dscf,Dold,Escf)
      nullify(listh)
      nullify(H,S)
      nullify(xijo)

!     Initialize some variables
      call init_Energies()
      
!
!     Set object debugging level
!
      call set_object_debug_level(0)
      i = fdf_get('DebugObjects.Node',0)
      if (Node==i) then
         if (fdf_get('DebugObjects',.false.)) then
            call set_object_debug_level(1)
         endif
      endif
         
!     Initialize CML (relies on reinit)
      call siesta_cml_init( )

! Set timer report file and threshold .................................
      threshold = fdf_get('timer_report_threshold', 0._dp)
      call timer_report( file=trim(slabel)//'.times', 
     .                   threshold=threshold )
      ! Note that the parallel timer might be inefficient for reports
      ! when large numbers of processors are used
      if ( PEXSINodes /= Nodes ) then
         use_parallel_timer = fdf_get('UseParallelTimer', .false.)
         use_tree_timer = fdf_get('UseTreeTimer', .true.)
      else
         use_parallel_timer = fdf_get('UseParallelTimer', .true.)
         use_tree_timer = fdf_get('UseTreeTimer', .false.)
      end if

!     Start time counter
!     Note new placement of this first use, so that
!     initialization is done after the setup of relevant variables

      call timer( 'siesta', 0 )
      call timer( 'siesta', 1 )
      call timer( 'Setup', 1 )

! Set allocation report level .........................................
! variables level and threshold imported from module siesta_options
      level = fdf_get('alloc_report_level', 0)
      threshold = fdf_get('alloc_report_threshold', 0._dp)
      call alloc_report( level=level, file=trim(slabel)//'.alloc',
     .                   threshold=threshold,
     &                   printNow=.false. )
      
      if ( IONode ) then
         ! Add citations
         call init_citation(trim(slabel))
         call add_citation("10.1088/0953-8984/14/11/302")
      end if


!     Initialise exchange-correlation functional information
#ifndef BSC_CELLXC
      call read_xc_info()
#else /* BSC_CELLXC */
      call setXC( )
#endif /* BSC_CELLXC */

!     Initialise force field component
      call inittwobody( )

!     Initialize pseudopotentials and atomic orbitals
      if (IOnode) call initatom( ns )
      call broadcast( ns )

      call broadcast_basis( )

      call register_rfs()

      atmonly = fdf_get( 'Atom-Setup-Only', .false. )
      if (atmonly) call bye( 'End of atom setup' )

!     Read geometry
      call struct_init( )      ! Sets na_u, isa, ucell

      ! Initialize the spin-spiral settings
      call init_spiral( ucell )
      
!     Initialize atom lists
      call initatomlists( )    ! Sets iza

!     early exit if only checking the structure
      struct_only = fdf_get( 'Output-Structure-Only', .false. )
      if (IONode) then
         call write_struct( ucell, na_u, isa, iza, xa )
         if (fdf_boolean('WriteCoorXmol',.false.)) then
            call coxmol(iza, xa, na_u)
         endif
         if (lUseZmatrix) then
            call write_canonical_ucell_and_Zmatrix(
     &                        filename="OUT.UCELL.ZMATRIX.INITIAL")
         endif
      endif
      if (struct_only) then
         call bye('End of structure processing')
      endif

!     Walltime control (numbers in seconds)
      if ( fdf_isphysical("MaxWalltime") ) then
         walltime_m = fdf_get("MaxWalltime", huge(1._dp), 's')
      else
         walltime_m = fdf_get("MaxWalltime", huge(1._dp))
      end if
      ! Period for clean-up operations
      if ( fdf_isphysical("MaxWalltime.Slack") ) then
         walltime_s = fdf_get("MaxWalltime.Slack", 5._dp, 's')
      else
         walltime_s = fdf_get("MaxWalltime.Slack", 5._dp)
      end if
      ! Note that the slack detracts from net available time
      walltime_warning = walltime_m - walltime_s

!-------------- Now we have the initial geometry

!     End of Initial Structure Processing
      if (Node.eq.0) then
        write(6,'(/,a,20("*"),a,28("*"))')
     &    'siesta: ', ' Simulation parameters '
        write(6,'(a)')  'siesta:'
        write(6,'(a)')  'siesta: The following are some of the '//
     &                           'parameters of the simulation.'
        write(6,'(a)')  'siesta: A complete list of the parameters '//
     &                           'used, including default values,'
        write(6,'(a,a)')'siesta: can be found in file out.fdf'
        write(6,'(a)')  'siesta:'
      endif

!     Allocate other arrays based on read sizes
      nullify(fa,cfa)
      call re_alloc( fa, 1, 3, 1, na_u, 'fa', 'siesta_init' )
      call re_alloc( cfa, 1, 3, 1, na_u, 'cfa', 'siesta_init' )

!     Read simulation data
      call read_options( na_u, ns, spin%Grid )

      call get_allowed_history_depth(n_dm_items)
      call new(DM_history,n_dm_items,"(DM history stack)")
      if (ionode) print "(a,i0)", "Size of DM history Fstack: ",
     $                            max_size(DM_history)

      qtot = qtot - charnet     ! qtot set in initatomlists
                                ! charnet set in redata
      if (IOnode) then
        write(6,fmt="(a,f12.6)") 'Total number of electrons: ', qtot
        write(6,fmt="(a,f12.6)") 'Total ionic charge: ', zvaltot
      endif
!
!     Warn the user: if not doing a direct optimization, the Zmatrix
!     coordinates are no longer updated. Only coordinates are treated.
!
      if (lUseZmatrix) then
         if (idyn .ne. 0) then
            write(6,"(a)")
     &      " *** WARNING: Zmatrix form will be used only for input !!"
            write(0,"(a)")
     &      " *** WARNING: Zmatrix form will be used only for input !!"
         endif
      endif

! Calculate spin populations for fixed spin case...
      if (fixspin) then
        if (.not.spin%Col)
     &    call die( 'siesta: ERROR: ' //
     &              'You can only fix the spin of the system' //
     &              ' for collinear spin polarized calculations.' )
        qs(1) = (qtot + ts) / 2.0_dp
        qs(2) = (qtot - ts) / 2.0_dp
      else
        qs(1:spin%spinor) = 0.0_dp
        do ispin = 1,spin%spinor
           qs(ispin) = qtot/spin%spinor
        end do
      endif

!     Find maximum interaction range
      if ( negl ) then
        rmaxh = 2.0_dp*(rmaxo + rmaxldau)
      else
        rmaxh = 2.0_dp*(rmaxo + max(rmaxkb,rmaxldau))
      endif

!     Madelung correction for charged systems
      if (charnet .ne. 0.0_dp) then
        call madelung(ucell, shape, charnet, Emad)
      endif

!     Parallel initialisation
      call initparallel( no_u, na_u, lasto, xa, ucell, rmaxh )
      if (IOnode) call show_distribution( )

!     Find number of locally stored orbitals and allocated related arrays
      call GetNodeOrbs(no_u,Node,Nodes,no_l)

!     Initialise arrays
      nullify( listhptr, numh )
      call re_alloc( listhptr, 1, max(1,no_l), 'listhptr', 'sparseMat' )
      call re_alloc( numh,     1, max(1,no_l), 'numh',     'sparseMat' )
      listhptr    = 0
      numh        = 0

!     Get number of eigenstates that need to be calculated
      neigwanted = fdf_get('NumberOfEigenStates',no_u)

!     Check number of eigenstates - cannot be larger than number of
!     basis functions or smaller than number of occupied states + 1
!     so that the Fermi level can be estimated
      do is = 1 , spin%spinor
         neigmin = nint(qs(is)/real(3 - spin%spinor, dp)) + 1
         neigwanted = max(neigwanted,neigmin)
      enddo
      neigwanted = min(neigwanted,no_u)

!     Find k-grid for Brillouin zone integration 
!     NOTE: We need to know whether gamma is .true. or
!     not early, in order to decide whether to use an 
!     auxiliary supercell for the calculation of matrix elements.
      call setup_Kpoint_grid( ucell )
      gamma = gamma_scf

!     Call initialisation of PDOS here since we need to check if 
!     the auxiliary supercell is needed for a non-gamma calculation
      call init_projected_DOS( )
      if (do_pdos) then
        gamma = gamma .and. gamma_pdos
      endif

      nullify(eo,qo)
      call re_alloc(eo, 1, no_u, 1, spin%spinor, 1, nkpnt, 'eo', 
     &              'siesta_init')
      call re_alloc(qo, 1, no_u, 1, spin%spinor, 1, nkpnt, 'qo', 
     &              'siesta_init')

      call setup_bands( )
      gamma = gamma .and. gamma_bands

      call setup_wf_kpoints( )
      gamma = gamma .and. gamma_wavefunctions

      call estimate_pol_kpoints( ucell )
      gamma = gamma .and. gamma_polarization
!      print *, "gamma, gamma_pol: ", gamma, gamma_polarization
!
!     User can request that the calculation is done with an explicit
!     auxiliary supercell and hermitian version, even if using only
!     the gamma point
!
      user_wants_supercell = fdf_get( 'ForceAuxCell', .false. )
      if (user_wants_supercell) gamma = .false.

!     Find required supercell
!     2*rmaxh is used to guarantee that two given orbitals in the
!     supercell can only overlap once
      if (gamma) then
        nsc(1:3) = 1
      else
        ! Discussion with Jose-Soler we get this size ! NPA
        do i = 1 , 3
          veclen = sqrt(ucell(1,i)**2+ucell(2,i)**2+ucell(3,i)**2)
          nsc(i) = 1 + 2 * ceiling( rmaxh / veclen )
        end do
        if ( .not. naiveauxcell ) 
     &       call exact_sc_ag(negl,ucell,na_u,isa,xa,nsc)
      end if

      mscell = 0.0_dp
      do i = 1, 3
        mscell(i,i) = nsc(i)
        nscold(i)   = nsc(i)
      enddo

!     Find auxiliary supercell (required only for k sampling)
      call superc( ucell, scell, nsc )

!     Initialise metadynamic forces if required
      call initmeta( )

      if (idyn .eq. 0) then
        inicoor = 0
        fincoor = nmove
      else if (idyn .ge. 1 .and. idyn .le. 5) then
        inicoor = istart
        fincoor = ifinal
      else if (idyn .eq. 6) then
        inicoor = 0
        fincoor = (ia2-ia1+1)*3*2
      else if (idyn .eq. 7) then
        call die( "'PHONON' support is deprecated" )
      else if (idyn .eq. 8) then
        inicoor = 0
        fincoor = huge(1)
#ifdef NCDF_4
      else if (idyn == 9) then
         if ( IONode ) then
            write(*,'(/,a)')'expcoord: '//repeat('*',69)
            write(*,'(a,i0)')
     &           'expcoord: Currently reached coordinate step = ',
     &           inicoor
            write(*,'(a,i0)')
     &           'expcoord: Final coordinate step = ',
     &           fincoor
            write(*,'(a,/)')'expcoord: '//repeat('*',69)
         end if
#endif
#ifdef SIESTA__FLOOK
      else if (idyn == 10) then
         inicoor = 0
         ! Controlled in external LUA machine
         fincoor = huge(1)
#endif
      else
        call die( 'siesta: wrong idyn' )
      endif

      ! initialize the fixed positions (we need to check whether the 
      ! fixed positions make sense)
      call init_fixed( ucell, na_u , isa, iza )
      call print_fixed( )

      ! Build initial velocities according to Maxwell-Bolzmann distribution....
      if (idyn .ne. 0 .and. idyn .ne. 6 .and. (.not. xv_file_read)) 
     &    call vmb( na_u, tempinit, amass, xa, isa, va )

      istp = 0
      
#ifdef TRANSIESTA
      call ts_init(spin%Grid,ucell,na_u,xa,lasto,no_u,inicoor,fincoor)
#endif

#ifdef SIESTA__FLOOK

      ! Populate dictionaries
      call dict_populate()

      ! Initialize LUA handle
      call slua_init(LUA)

      ! Call a preprocess step right after initialization
      call slua_call(LUA,LUA_INITIALIZE)

#endif

!     Output memory use before main loop
!!      call printmemory( 6, 0 )

!     Initialization now complete. Flush stdout.
      if (ionode) call pxfflush( 6 )

      call timer( 'Setup', 2 )

!--------------------------------------------------------------------------- END

      END subroutine siesta_init

      END MODULE m_siesta_init
