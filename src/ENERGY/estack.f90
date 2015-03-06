! ==============================================================================
! Subroutine: ESTACK (ISEQ,I,J,IP,JP,N,ES)
! 
! Purpose: Computes the energy of helix stacking between two bp using
!          the empirical MFOLD 3.0 energy function.
!
! Method: Uses the MFOLD 3.0 energy function for RNA @ T=37 given by:
!
!         ES = E_stack
!
!               5' (I) X W (IP) 3'
!               3' (J) Y Z (JP) 5'
!
!               NOTE: I < IP and JP < J
!
! Arguments:
!
!          ISEQ - Array of length N containing the sequence
!                 in numerical code (A=1,C=2,G=3,U=4)
!             I - Nucleotide position of the first basepair 5'.
!             J - Nucleotide position of the first basepair 3'.
!            IP - Nucleotide position of the second basepair 5'.
!            JP - Nucleotide position of the second basepair 3'.
!             N - Number of nucleotides in the sequence.
!            ES - (OUTPUT) MFOLD 3.0 energy of the stacking.
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            01/01/2015   Original Code
!
! Dependencies:
!
! Modules -
! Functions -
! Subroutines - TSTACK
!
! Author(s): Eric Dykeman
!            Copyright (c) 2015 (Please Refer to LICENCE)
!
! ==============================================================================

      SUBROUTINE ESTACK (ISEQ,I,J,IP,JP,N,ES)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        INTEGER, INTENT(IN) :: i,j,ip,jp,n
        INTEGER, INTENT(IN) :: iseq(n)

        REAL, INTENT(OUT) :: es

        !=== VARIABLES ===!

        INTEGER :: list(4)


        es = 0.0e0

        !5' (i) A X (ip) 3'
        !3' (j) U Y (jp) 5'

        list(1) = iseq(i)
        list(2) = iseq(j)
        list(3) = iseq(ip)
        list(4) = iseq(jp)

        CALL TSTACK (list,es)

        RETURN

      END SUBROUTINE ESTACK
