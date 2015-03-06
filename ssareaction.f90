! ==============================================================================
! Subroutine: SSAREACTION (RNA,ISEED,TIME,TOUT)
! 
! Purpose: Calculates an RNA folding reaction to fire based on the
!          (s)tochastic (s)imulation (a)lgorithm of Gillespie.
!
! Method: See Gillespie, Daniel T. (1977). "Exact Stochastic Simulation
!         of Coupled Chemical Reactions". J. Phys. Chem. 81, 2340.
!
! Arguments:
!
!           RNA - Class structure containing information on the
!                 RNA secondary structure and possible reactions.
!         ISEED - Integer seed for the random number generator.
!          TIME - Current Time
!          TOUT - Time to write next trajectory output.
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            01/01/2015   Original Code
!
! Dependancies:
!
! Modules - CLASS_RNAFOLD RNAVar
! Functions - RANDOM
! Subroutines - ESTRUC V2CT
!
! Author(s): Eric Dykeman
!            Copyright (c) 2015 (Please Refer to LICENCE)
!
! ==============================================================================

      SUBROUTINE SSAREACTION (RNA,ISEED,TIME,TOUT)

        USE RNAVar, ONLY : mxnt

        USE Class_RNAFold

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE (RNA_STRUC), INTENT(INOUT) :: rna

        INTEGER, INTENT(INOUT) :: iseed

        DOUBLE PRECISION, INTENT(INOUT) :: time
        DOUBLE PRECISION, INTENT(IN) :: tout

        !=== VARIABLES ===!

        REAL :: e

        INTEGER :: i,j,k,n
        INTEGER :: n1,n2,indx

        CHARACTER :: fld(mxnt)

        DOUBLE PRECISION :: r,tau,random
        DOUBLE PRECISION :: atot,amax


        !=== Total Transition Rate ===!

        n = rna% nsum / 2

        atot = rna% psum(n)


        !=== Compute Time Increment ===!

        r = RANDOM(iseed)

        tau = DLOG(1.0d0/r)
        tau = tau / atot

        time = time + tau

        !=== Output Current Structure? ===!

        IF ( time > tout ) THEN

          CALL V2CT (rna%ibsp,fld,'V',rna%n)

          CALL ESTRUC (rna%iseq,rna%ibsp,rna%n,e)

          WRITE(2,'(2E16.8)')tout,e
          WRITE(2,'(10000A1)')(rna%seq(i),i=1,rna%n)
          WRITE(2,'(10000A1)')(    fld(i),i=1,rna%n)

        ENDIF


        !=== Fire Reaction ===!

        r = RANDOM(iseed)
        amax = r * atot

        !=== Find Reaction to Fire ===!

        i = n

        DO WHILE ( MOD(n,2) == 0 )

          n = n / 2
          j = i - n

          r = rna% psum(j)

          IF ( r >= amax ) THEN

            i = j

          ELSE

            i = i + n
            amax = amax - r

          ENDIF

        ENDDO

        !=== Choose Between i and i+1 ===!

        r = rna% ptot(i)

        IF ( r >= amax ) THEN

          indx = i

        ELSE

          indx = i + 1
          amax = amax - r

        ENDIF

        !=== Fire Reaction ===!

        CALL LOOP_FIRE (rna,indx,amax)

        RETURN

      END SUBROUTINE SSAREACTION
