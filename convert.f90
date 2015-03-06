! ==============================================================================
! Subroutine: CONVERT (GSEQ,ISEQ,NN)
! 
! Purpose: Converts a character sequence of A C G U into numerical code.
!
! Method: Converts a nucleic acid sequence to the following numerical
!         code:
!
!              A = 1
!              C = 2
!              G = 3
!            T/U = 4
!
! Arguments:
!
!          GSEQ - (INPUT) Array of length NN containing the sequence
!                 of characters to be converted into numbers.
!          ISEQ - (OUTPUT) Array of length NN containing the sequence
!                 in the numerical code.
!            NN - Number of characters in the sequence.
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

      SUBROUTINE CONVERT (GSEQ,ISEQ,NN)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        INTEGER, INTENT(IN) :: nn
        INTEGER, INTENT(OUT) :: iseq(nn)

        CHARACTER, INTENT(IN) :: gseq(nn)

        !=== VARIABLES ===!

        INTEGER :: i


        DO i=1,nn

          iseq(i) = 0

          IF ( gseq(i) == 'A' ) iseq(i) = 1
          IF ( gseq(i) == 'a' ) iseq(i) = 1

          IF ( gseq(i) == 'C' ) iseq(i) = 2
          IF ( gseq(i) == 'c' ) iseq(i) = 2

          IF ( gseq(i) == 'G' ) iseq(i) = 3
          IF ( gseq(i) == 'g' ) iseq(i) = 3

          IF ( gseq(i) == 'U' ) iseq(i) = 4
          IF ( gseq(i) == 'u' ) iseq(i) = 4

          IF ( gseq(i) == 'T' ) iseq(i) = 4
          IF ( gseq(i) == 't' ) iseq(i) = 4

        ENDDO

        RETURN

      END SUBROUTINE CONVERT
