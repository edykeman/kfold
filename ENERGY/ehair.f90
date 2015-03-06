! ==============================================================================
! Subroutine: EHAIR (ISEQ,I,J,N,EH)
! 
! Purpose: Computes the energy of an RNA hairpin turn using the
!          empirical MFOLD 3.0 energy function.
!
! Method: Uses the MFOLD 3.0 energy function for RNA @ T=37 given by:
!
!         EH = E_entropic + E_stack + E_bonus + E_penalty
!
!               5' (I) X ... loop 3'
!               3' (J) Y ... loop 5'
!
!               NOTE: I < J
!
! Arguments:
!
!          ISEQ - Array of length N containing the sequence
!                 in numerical code (A=1,C=2,G=3,U=4)
!             I - Nucleotide position of the loop basepair 5'.
!             J - Nucleotide position of the loop basepair 3'.
!             N - Number of nucleotides in the sequence.
!            EH - (OUTPUT) MFOLD 3.0 energy of the hairpin sequence.
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            01/01/2015   Original Code
!
! Dependencies:
!
! Modules - RNAvar
! Functions -
! Subroutines - TSTACKH, TLOOP
!
! Author(s): Eric Dykeman
!            Copyright (c) 2015 (Please Refer to LICENCE)
!
! ==============================================================================

      SUBROUTINE EHAIR (ISEQ,I,J,N,EH)

        USE RNAvar, ONLY : eau,beta

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        INTEGER, INTENT(IN) :: i,j,n
        INTEGER, INTENT(IN) :: iseq(n)

        REAL, INTENT(OUT) :: eh

        !=== VARIABLES ===!

        INTEGER :: list(6)
        INTEGER :: k,ic,nl

        REAL :: x,c,elh(30)

        DATA (elh(k),k=1,30) / 0.00e0,0.00e0,5.70e0,5.60e0,5.60e0,&
                             & 5.40e0,5.90e0,5.60e0,6.40e0,6.50e0,&
                             & 6.60e0,6.70e0,6.80e0,6.90e0,6.90e0,&
                             & 7.00e0,7.10e0,7.10e0,7.20e0,7.20e0,&
                             & 7.30e0,7.30e0,7.40e0,7.40e0,7.50e0,&
                             & 7.50e0,7.50e0,7.60e0,7.60e0,7.70e0 /


        eh = 0.0e0

        nl = j - i - 1

        c = 1.750e0 / REAL(beta)


        !=== TERM 1 --> Entropic Term ===!

        IF ( nl > 30 ) THEN

          x = REAL(nl) / 30.0e0
          x = c * LOG(x)

          eh = elh(30) + x

        ELSE

          eh = elh(nl)

        ENDIF


        !=== TERM 2 --> Stacking Energy ===!

        IF ( nl > 3 ) THEN

          ! 5' (i) A X (i+1) LOOP
          ! 3' (j) U Y (j-1) LOOP

          list(1) = iseq(i)
          list(2) = iseq(j)
          list(3) = iseq(i+1)
          list(4) = iseq(j-1)

          CALL TSTACKH (list,eh)

        ENDIF


        !=== TERM 3 --> Bonuses ===!

        !=== Tetra-loop Bonus ===!

        IF ( nl == 4 ) THEN

          list(1) = iseq(i)
          list(2) = iseq(i+1)
          list(3) = iseq(i+2)
          list(4) = iseq(i+3)
          list(5) = iseq(i+4)
          list(6) = iseq(i+5)

          CALL TLOOP (list,eh)

        ENDIF

        !=== GGG Hairpin Bonus ===!

        IF ( iseq(i) == 3 .and. iseq(j) == 4 ) THEN

          ic = 0

          DO k=MAX(1,i-2),i
          IF ( iseq(k) == 3 ) ic = ic + 1
          ENDDO

          IF ( ic == 3 ) eh = eh - 2.20e0

        ENDIF


        !=== TERM 4 --> Penalties ===!

        !=== Poly C Penalty ===!

        ic = 0

        DO k=i+1,j-1
        IF ( iseq(k) == 2 ) ic = ic + 1
        ENDDO

        IF ( ic == nl ) THEN
        IF ( nl == 3 ) THEN
          eh = eh + 1.40e0
        ELSE
          eh = eh + 1.60e0
          eh = eh + 0.30e0 * REAL(ic)
        ENDIF
        ENDIF

        !=== A-U / G-U closing a Tri-loop ===!

        IF ( nl == 3 ) THEN

          IF ( iseq(i) == 4 ) eh = eh + eau
          IF ( iseq(j) == 4 ) eh = eh + eau

        ENDIF

        RETURN

      END SUBROUTINE EHAIR
