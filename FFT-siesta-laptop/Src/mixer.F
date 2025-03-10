! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
      MODULE m_mixer
      private
      public :: mixer

      CONTAINS
 
      subroutine mixer( iscf )

      use precision,        only: dp
      use siesta_options,   only: mixH, mix_scf_first
      use siesta_options,   only: mullipop, muldeb
      ! Spin mixing options
      use m_mixing_scf, only: MIX_SPIN_ALL, MIX_SPIN_SPINOR
      use m_mixing_scf, only: MIX_SPIN_SUM, MIX_SPIN_SUM_DIFF
      use m_mixing_scf, only: mix_spin

      use sparse_matrices,  only: Dold, Dscf, Hold, H, S
      use sparse_matrices,  only: maxnh, numh, listhptr, listh
      use siesta_geom,      only: na_u, isa
      use atomlist,         only: iaorb, iphorb, lasto, no_u, no_l
      use atomlist,         only: indxuo

      use m_mixing_scf, only: scf_mix
      use m_mixing, only: mixing
      
      use m_spin,           only: h_spin_dim, SpOrb
      use parallel,         only: IONode

      implicit none

      real(dp), pointer :: Xin(:,:), Xout(:,:)
      real(dp), allocatable :: F(:,:)

      integer,  intent(in)  :: iscf

      real(dp)              :: dtmp
      integer               :: iiscf, i, is

      external :: mulliken

!-------------------------------------------------------------------- BEGIN

      call timer( 'MIXER', 1 )

      if ( mixH ) then
         ! Mix Hamiltonian
         Xin => Hold
         Xout => H
      else
         ! Mix density matrix
         Xin => Dold
         Xout => Dscf
      end if


      ! Create residual function to minimize
      allocate(F(maxnh,h_spin_dim))

!$OMP parallel default(shared), private(dtmp,i)

      ! prepare input...
      select case ( mix_spin )
      case ( MIX_SPIN_ALL , MIX_SPIN_SPINOR )
         ! Xin and Xout are unchanged...
      case ( MIX_SPIN_SUM , MIX_SPIN_SUM_DIFF )

         ! transfer spin density to
         ! spin-sum and spin-difference
!$OMP do
         do i = 1 , size(Xin,1)
            dtmp = (Xin(i,1) + Xin(i,2)) * 0.5_dp
            Xin(i,2) = (Xin(i,1) - Xin(i,2)) * 0.5_dp
            Xin(i,1) = dtmp
         end do
!$OMP end do nowait

!$OMP do
         do i = 1 , size(Xout,1)
            dtmp = (Xout(i,1) + Xout(i,2)) * 0.5_dp
            Xout(i,2) = (Xout(i,1) - Xout(i,2)) * 0.5_dp
            Xout(i,1) = dtmp
         end do
!$OMP end do

      end select

      ! Calculate the residual
!$OMP workshare
      F = Xout - Xin
!$OMP end workshare

!$OMP end parallel

      ! Call mixing routine
      ! Xin contains the input element, F contains
      !  F = Xout - Xin
      ! Upon exit Xout contains the mixed quantity
      select case ( mix_spin )
      case ( MIX_SPIN_ALL )
         call mixing( scf_mix, maxnh, h_spin_dim, Xin, F, Xout)
      case ( MIX_SPIN_SPINOR )
         call mixing( scf_mix, maxnh, h_spin_dim, Xin, F, Xout,
     &        nsub=2)
      case ( MIX_SPIN_SUM )
         call mixing( scf_mix, maxnh, h_spin_dim, Xin, F, Xout,
     &        nsub=1)
      case ( MIX_SPIN_SUM_DIFF )
         call mixing( scf_mix, maxnh, h_spin_dim, Xin, F, Xout,
     &        nsub=2)
      end select

!$OMP parallel default(shared), private(dtmp,i,is)

      ! correct output...
      select case ( mix_spin )
      case ( MIX_SPIN_ALL , MIX_SPIN_SPINOR )
         ! Xin and Xout are unchanged...
      case ( MIX_SPIN_SUM , MIX_SPIN_SUM_DIFF )

         ! transfer spin-sum and difference to
         ! spin-up and spin-down

!$OMP do
         do i = 1 , size(Xin,1)
            dtmp = Xin(i,1) + Xin(i,2)
            Xin(i,2) = Xin(i,1) - Xin(i,2)
            Xin(i,1) = dtmp
         end do
!$OMP end do nowait

!$OMP do
         do i = 1 , size(Xout,1)
            dtmp = Xout(i,1) + Xout(i,2)
            Xout(i,2) = Xout(i,1) - Xout(i,2)
            Xout(i,1) = dtmp
         end do
!$OMP end do

      end select

      
      ! Correctly handle mixed quantity
      ! move over mixed quantity to "input"
      if ( mix_scf_first ) then
         
         ! always allow mixing
!$OMP do
         do is = 1 , size(Xin, 2)
            do i = 1 , size(Xin, 1)
               Xin(i,is) = Xout(i,is)
            end do
         end do
!$OMP end do nowait
         
      else if ( iscf == 1 ) then
         
         ! We are not allowed to mix the first SCF
         ! Remember that Xout contains the MIXED
         ! quantity and F is the current residual
         ! hence:
         !   Xout = Xin + F
!$OMP do
         do is = 1 , size(Xin, 2)
            do i = 1 , size(Xin, 1)
               Xin(i,is) = Xin(i,is) + F(i,is)
               ! we require that Xin == Xout:
               Xout(i,is) = Xin(i,is)
            end do
         end do
!$OMP end do nowait

      end if

!$OMP end parallel
      
      deallocate(F)

      ! Print populations at each SCF step, if requested
      ! Note that this is after mixing, which is not
      ! entirely correct. It should be moved to the top,
      ! or done somewhere else.

      if (muldeb) then 
         if (IONode)
     &        write (6,"(/a)")
     &        'siesta: Mulliken populations after mixing'
         if ( SpOrb ) then
            call moments( mullipop, na_u, no_u, maxnh, numh, listhptr,
     .           listh, S, Dscf, isa, lasto, iaorb, iphorb,
     .           indxuo )
         else
            call mulliken( mullipop, h_spin_dim, na_u, no_u, maxnh,
     &           numh, listhptr, listh, S, Dscf, isa,
     &           lasto, iaorb, iphorb )
         endif
      endif

      call timer( 'MIXER', 2 )

!-------------------------------------------------------- END
      END subroutine mixer

      End MODULE m_mixer


