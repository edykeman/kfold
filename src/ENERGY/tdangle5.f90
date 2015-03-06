! ==============================================================================
! Subroutine: TDANGLE5 (LIST,ES)
! 
! Purpose: Performs a table lookup for the stacking interaction
!          between a dangling base and a basepair in a helix.
!
! Method: Uses the MFOLD 3.0 energy function table for RNA @ T=37.
!
! Arguments:
!
!          LIST - Array of length 3 containing the nucleotides in
!                 numerical code (A=1,C=2,G=3,U=4) for the
!                 following locations:
!
!                 5' (1) A       3'
!                 3' (2) U X (3) 5'
!
!                 where LIST(1) = letter code for position 1 etc.
!
!            ES - (OUTPUT) MFOLD 3.0 stacking energy of the sequence
!                 provided in LIST.
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
! Subroutines -
!
! Author(s): Eric Dykeman
!            Copyright (c) 2015 (Please Refer to LICENCE)
!
! ==============================================================================

      SUBROUTINE TDANGLE5 (LIST,ES)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        INTEGER, INTENT(IN) :: list(3)
        REAL, INTENT(INOUT) :: es

        !=== VARIBLES ===!

        INTEGER :: i,i1,i2,i3
        REAL :: au(4),cg(4),gc(4)
        REAL :: ua(4),gu(4),ug(4)

        DATA (au(i),i=1,4) / -0.30e0,-0.10e0,-0.20e0,-0.20e0 /

        DATA (cg(i),i=1,4) / -0.20e0,-0.30e0,-0.00e0,-0.00e0 /

        DATA (gc(i),i=1,4) / -0.50e0,-0.30e0,-0.20e0,-0.10e0 /

        DATA (ua(i),i=1,4) / -0.30e0,-0.30e0,-0.40e0,-0.20e0 /

        DATA (gu(i),i=1,4) / -0.30e0,-0.10e0,-0.20e0,-0.20e0 /

        DATA (ug(i),i=1,4) / -0.30e0,-0.30e0,-0.40e0,-0.20e0 /


        ! 5' (1) A       3' !
        ! 3' (2) U X (3) 5' !

        i1 = list(1)
        i2 = list(2)
        i3 = list(3)

        IF ( i1 == 1 .and. i2 == 4 ) es = es + au(i3)
        IF ( i1 == 2 .and. i2 == 3 ) es = es + cg(i3)
        IF ( i1 == 3 .and. i2 == 2 ) es = es + gc(i3)
        IF ( i1 == 4 .and. i2 == 1 ) es = es + ua(i3)
        IF ( i1 == 3 .and. i2 == 4 ) es = es + gu(i3)
        IF ( i1 == 4 .and. i2 == 3 ) es = es + ug(i3)

        RETURN

      END SUBROUTINE TDANGLE5
