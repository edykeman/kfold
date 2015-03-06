! ==============================================================================
! Subroutine: TLOOP (LIST,EL)
! 
! Purpose: Performs a table lookup for the special bonus energies
!          for various tetra-loops.
!
! Method: Uses the MFOLD 3.0 energy function table for RNA @ T=37.
!
! Arguments:
!
!          LIST - Array of length 6 containing the nucleotides in
!                 numerical code (A=1,C=2,G=3,U=4) for the
!                 following locations:
!
!                     1 2 3 4 5 6
!                 5'  A W X Y Z U   3'
!                       L O O P
!
!                 where LIST(1) = letter code for position 1 etc.
!
!            EL - (OUTPUT) MFOLD 3.0 tetra-loop bonus energy of the
!                 sequence provided in LIST.
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

      SUBROUTINE TLOOP (LIST,EL)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        INTEGER, INTENT(IN) :: list(6)
        REAL, INTENT(INOUT) :: el

        !=== VARIABLES ===!

        INTEGER :: i
        CHARACTER :: cl(6)

        REAL :: ener(30)
        CHARACTER (LEN=6) :: cwrk,cseq(30)


        DATA (cseq(i),i=1,30) / 'GGGGAC','GGUGAC','CGAAAG','GGAGAC',&
            & 'CGCAAG','GGAAAC','CGGAAG','CUUCGG','CGUGAG','CGAAGG',&
            & 'CUACGG','GGCAAC','CGCGAG','UGAGAG','CGAGAG','AGAAAU',&
            & 'CGUAAG','CUAACG','UGAAAG','GGAAGC','GGGAAC','UGAAAA',&
            & 'AGCAAU','AGUAAU','CGGGAG','AGUGAU','GGCGAC','GGGAGC',&
            & 'GUGAAC','UGGAAA' /

        DATA (ener(i),i=1,30) / -3.00e0,-3.00e0,-3.00e0,-3.00e0,&
              & -3.00e0,-3.00e0,-3.00e0,-3.00e0,-3.00e0,-2.50e0,&
              & -2.50e0,-2.50e0,-2.50e0,-2.50e0,-2.00e0,-2.00e0,&
              & -2.00e0,-2.00e0,-2.00e0,-1.50e0,-1.50e0,-1.50e0,&
              & -1.50e0,-1.50e0,-1.50e0,-1.50e0,-1.50e0,-1.50e0,&
              & -1.50e0,-1.50e0 /


        !      1 2 3 4 5 6      !
        !  5'  A W X Y Z U   3' !
        !        L O O P        !

        DO i=1,6

          IF ( list(i) == 1 ) cl(i) = 'A'
          IF ( list(i) == 2 ) cl(i) = 'C'
          IF ( list(i) == 3 ) cl(i) = 'G'
          IF ( list(i) == 4 ) cl(i) = 'U'

        ENDDO

        WRITE(cwrk,'(6A1)')(cl(i),i=1,6)


        DO i=1,30
        IF ( cseq(i) == cwrk ) el = el + ener(i)
        ENDDO

        RETURN

      END SUBROUTINE TLOOP
