! ==============================================================================
! Subroutine: DELTAG_HI (R,I,J,DG)
! 
! Purpose: Computes the difference in free energy of an RNA Helix due
!          to a deletion of the internal i-j base pair using the emperical
!          MFOLD 3.0 energy function.
!
! Method: Uses the MFOLD 3.0 energy function for RNA @ T=37.
!
! Arguments:
!
!             R - Class containing information about the RNA fold and
!                 the RNA sequence.
!             I - Nucleotide position of the 5' nucleotide.
!             J - Nucleotide position of the 3' nucleotide.
!            DG - (OUTPUT) The energy difference from deletion of the
!                 base-pair i-j
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            01/01/2015   Original Code
!
! Dependencies:
!
! Modules - Class_RNAFold, RNAVar
! Functions -
! Subroutines - ESTACK, EBULGE
!
! Author(s): Eric Dykeman
!            Copyright (c) 2015 (Please Refer to LICENCE)
!
! ==============================================================================

      SUBROUTINE DELTAG_HI (R,I,J,DG)

        USE Class_RNAFold

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE(RNA_STRUC), INTENT(IN) :: r
        INTEGER, INTENT(IN) :: i,j

        DOUBLE PRECISION, INTENT(OUT) :: dg

        !=== VARIABLES ===!

        INTEGER :: n

        REAL :: ei,ef,es


        n = r% n

        !=== INITIAL ENERGY ===!

        CALL ESTACK (r%iseq,i-1,j+1,i,j,n,ei)
        CALL ESTACK (r%iseq,i,j,i+1,j-1,n,es)

        ei = ei + es


        !=== FINAL ENERGY ===!

        CALL EBULGE (r%iseq,i-1,j+1,i+1,j-1,n,ef)

        dg = DBLE(ef) - DBLE(ei)

        RETURN

      END SUBROUTINE DELTAG_HI
