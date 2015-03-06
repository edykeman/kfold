! ==============================================================================
! Subroutine: EDANGLE (ISEQ,I,J,K,N,ED)
! 
! Purpose: Computes the energy of a dangling nucleotide over a basepair
!          using the empirical MFOLD 3.0 energy function.
!
! Method: Uses the MFOLD 3.0 energy function for RNA @ T=37 given by:
!
!         ED = E_dangle
!
!               5' (I) X       3'
!               3' (J) Y Z (K) 5'
!
!
! Arguments:
!
!          ISEQ - Array of length N containing the sequence
!                 in numerical code (A=1,C=2,G=3,U=4)
!             I - Nucleotide position of the basepair 5'.
!             J - Nucleotide position of the basepair 3'.
!             K - Nucleotide position of the dangling nucleotide.
!             N - Number of nucleotides in the sequence.
!            ED - (OUTPUT) MFOLD 3.0 energy of the dangling nucleotide.
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
! Subroutines - TDANGLE5, TDANGLE3
!
! Author(s): Eric Dykeman
!            Copyright (c) 2015 (Please Refer to LICENCE)
!
! ==============================================================================

      SUBROUTINE EDANGLE (ISEQ,I,J,K,N,ED)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        INTEGER, INTENT(IN) :: i,j,k,n
        INTEGER, INTENT(IN) :: iseq(n)

        REAL, INTENT(OUT) :: ed

        !=== VARIABLES ===!

        INTEGER :: list(3)


        ed = 0.0e0

        IF ( k < 1 ) RETURN
        IF ( k > n ) RETURN


        IF ( k == i+1 ) THEN

          !5' (i) A X (k) 3'
          !3' (j) U       5'

          list(1) = iseq(i)
          list(2) = iseq(j)
          list(3) = iseq(k)

          CALL TDANGLE3 (list,ed)

        ENDIF

        IF ( k == j-1 ) THEN

          !5' (i) A       3'
          !3' (j) U X (k) 5'

          list(1) = iseq(i)
          list(2) = iseq(j)
          list(3) = iseq(k)

          CALL TDANGLE5 (list,ed)

        ENDIF

        IF ( k == j+1 ) THEN

          !5'       A (i) 3'
          !3' (k) X U (j) 5'

          list(1) = iseq(j)
          list(2) = iseq(i)
          list(3) = iseq(k)

          CALL TDANGLE3 (list,ed)

        ENDIF

        IF ( k == i-1 ) THEN

          !5' (k) X A (i) 3'
          !3'       U (j) 5'

          list(1) = iseq(j)
          list(2) = iseq(i)
          list(3) = iseq(k)

          CALL TDANGLE5 (list,ed)

        ENDIF

        RETURN

      END SUBROUTINE EDANGLE
