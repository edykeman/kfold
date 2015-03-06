! ==============================================================================
! Subroutine: TSTACKI (LIST,ES)
! 
! Purpose: Performs a table lookup for the stacking interaction of
!          the two nucleotides positioned over a closing basepair in
!          an internal loop.
!
! Method: Uses the MFOLD 3.0 energy function table for RNA @ T=37.
!
! Arguments:
!
!          LIST - Array of length 4 containing the nucleotides in
!                 numerical code (A=1,C=2,G=3,U=4) for the
!                 following locations:
!
!                 5' (1) A X (3) 3' INTERNAL LOOP!
!                 3' (2) U Y (4) 5' INTERNAL LOOP!
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

      SUBROUTINE TSTACKI (LIST,ES)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        INTEGER, INTENT(IN) :: list(4)
        REAL, INTENT(INOUT) :: es

        !=== VARIABLES ===!

        INTEGER :: i,i1,i2,i3,i4
        REAL :: au(4,4),cg(4,4),gc(4,4)
        REAL :: ua(4,4),gu(4,4),ug(4,4)

        DATA (au(1,i),i=1,4) /  0.70e0, 0.70e0,-0.40e0, 0.70e0 /
        DATA (au(2,i),i=1,4) /  0.70e0, 0.70e0, 0.70e0, 0.70e0 /
        DATA (au(3,i),i=1,4) / -0.40e0, 0.70e0, 0.70e0, 0.70e0 /
        DATA (au(4,i),i=1,4) /  0.70e0, 0.70e0, 0.70e0, 0.00e0 /

        DATA (cg(1,i),i=1,4) /  0.00e0, 0.00e0,-1.10e0, 0.00e0 /
        DATA (cg(2,i),i=1,4) /  0.00e0, 0.00e0, 0.00e0, 0.00e0 /
        DATA (cg(3,i),i=1,4) / -1.10e0, 0.00e0, 0.00e0, 0.00e0 /
        DATA (cg(4,i),i=1,4) /  0.00e0, 0.00e0, 0.00e0,-0.70e0 /

        DATA (gc(1,i),i=1,4) /  0.00e0, 0.00e0,-1.10e0, 0.00e0 /
        DATA (gc(2,i),i=1,4) /  0.00e0, 0.00e0, 0.00e0, 0.00e0 /
        DATA (gc(3,i),i=1,4) / -1.10e0, 0.00e0, 0.00e0, 0.00e0 /
        DATA (gc(4,i),i=1,4) /  0.00e0, 0.00e0, 0.00e0,-0.70e0 /

        DATA (ua(1,i),i=1,4) /  0.70e0, 0.70e0,-0.40e0, 0.70e0 /
        DATA (ua(2,i),i=1,4) /  0.70e0, 0.70e0, 0.70e0, 0.70e0 /
        DATA (ua(3,i),i=1,4) / -0.40e0, 0.70e0, 0.70e0, 0.70e0 /
        DATA (ua(4,i),i=1,4) /  0.70e0, 0.70e0, 0.70e0, 0.00e0 /

        DATA (gu(1,i),i=1,4) /  0.70e0, 0.70e0,-0.40e0, 0.70e0 /
        DATA (gu(2,i),i=1,4) /  0.70e0, 0.70e0, 0.70e0, 0.70e0 /
        DATA (gu(3,i),i=1,4) / -0.40e0, 0.70e0, 0.70e0, 0.70e0 /
        DATA (gu(4,i),i=1,4) /  0.70e0, 0.70e0, 0.70e0, 0.00e0 /

        DATA (ug(1,i),i=1,4) /  0.70e0, 0.70e0,-0.40e0, 0.70e0 /
        DATA (ug(2,i),i=1,4) /  0.70e0, 0.70e0, 0.70e0, 0.70e0 /
        DATA (ug(3,i),i=1,4) / -0.40e0, 0.70e0, 0.70e0, 0.70e0 /
        DATA (ug(4,i),i=1,4) /  0.70e0, 0.70e0, 0.70e0, 0.00e0 /


        ! 5' (1) A X (3) 3' INTERNAL LOOP!
        ! 3' (2) U Y (4) 5' INTERNAL LOOP!

        i1 = list(1)
        i2 = list(2)
        i3 = list(3)
        i4 = list(4)

        IF ( i1 == 1 .and. i2 == 4 ) es = es + au(i3,i4)
        IF ( i1 == 2 .and. i2 == 3 ) es = es + cg(i3,i4)
        IF ( i1 == 3 .and. i2 == 2 ) es = es + gc(i3,i4)
        IF ( i1 == 4 .and. i2 == 1 ) es = es + ua(i3,i4)
        IF ( i1 == 3 .and. i2 == 4 ) es = es + gu(i3,i4)
        IF ( i1 == 4 .and. i2 == 3 ) es = es + ug(i3,i4)

        RETURN

      END SUBROUTINE TSTACKI
