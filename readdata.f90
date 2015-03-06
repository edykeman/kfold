! ==============================================================================
! Subroutine: READDATA
! 
! Purpose: Reads in the RNA energy data files and sets up some tables.
!
! Method:
!
! Arguments:
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            01/01/2015   Original Code
!
! Dependancies:
!
! Modules - RNAVar
! Functions -
! Subroutines -
!
! Author(s): Eric Dykeman
!            Copyright (c) 2015 (Please Refer to LICENCE)
!
! ==============================================================================

      SUBROUTINE READDATA

        USE RNAVar, ONLY : eau,eaup,iwc

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        !=== VARIABLES ===!

        INTEGER :: i,j,k


        !=== Read Data Files ===!

        !=== Setup Tables ===!

        !=== A=1,C=2,G=3,U=4 ===!

        iwc(:,:) = 0
        iwc(1,4) = 1
        iwc(4,1) = 1
        iwc(2,3) = 1
        iwc(3,2) = 1
        iwc(3,4) = 1
        iwc(4,3) = 1

        eaup(:,:) = 0.0e0
        eaup(1,4) = eau
        eaup(4,1) = eau
        eaup(3,4) = eau
        eaup(4,3) = eau

        RETURN

      END SUBROUTINE READDATA
