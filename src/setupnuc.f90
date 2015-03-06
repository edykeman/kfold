! ==============================================================================
! Subroutine: SETUPNUC (N)
! 
! Purpose: Creates a table of nucleation probabilites between pairs of
!          nucleotides based on a worm like chain model.
!          See:
!          (1) Toan et al. J. Phys. Chem. B 112, 6094-6106 (2008).
!          (2) S. Kuznetsov and A. Ansari "A kinetic zipper model with
!          interchain interactions applied to nucleic acid hairpin
!          folding kinetics", Biophys J 102, 1001-111 (2012).
!
! Method: Uses the formula from:
!
! Arguments:
!
!             N - Number of nucleotides in the sequence.
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            01/01/2015   Original Code
!
! Dependancies:
!
! Modules -
! Functions -
! Subroutines -
!
! Author(s): Eric Dykeman
!            Copyright (c) 2015 (Please Refer to LICENCE)
!
! ==============================================================================

      SUBROUTINE SETUPNUC (N)

        USE RNAVar, ONLY : pnuc,beta

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        INTEGER, INTENT(IN) :: n

        !=== VARIABLES ===!

        INTEGER :: i

        DOUBLE PRECISION :: x,xi,e

        DOUBLE PRECISION, PARAMETER :: c = 0.1785714290d0
!       DOUBLE PRECISION, PARAMETER :: c2= 3.9274668195d8  !rate in (1/s)
        DOUBLE PRECISION, PARAMETER :: c2= 3.9274668195d2  !rate in (1/uS)
        DOUBLE PRECISION, PARAMETER :: xp= 4.0000000000d0


        IF ( ALLOCATED (pnuc) ) DEALLOCATE (pnuc)

        ALLOCATE (pnuc(n))

        pnuc(:) = 0.0d0

        DO i=5,n

          x  = c * DBLE(i-1)
          xi = 1.0d0 / x

          IF ( x <= xp ) THEN

            e = -7.0270d0 * xi + 0.4920d0 * x
            x = 84.90d0 * ( xi ** 5.50d0 )
            x = c2 * x / beta

            pnuc(i) = x * DEXP(e)

          ELSE

            x = xi * xi

            e = 1.0d0 - 0.6250d0 * xi - 0.12343750d0 * x
            x = c2 * x / beta

            pnuc(i) = x * e

          ENDIF

        ENDDO

        RETURN

      END SUBROUTINE SETUPNUC
