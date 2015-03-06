! ==============================================================================
! Subroutine: ESTRUC (ISEQ,IBSP,N,E)
! 
! Purpose: Computes the energy of an RNA SS using the emperical
!          MFOLD 3.0 energy function.
!
! Method: Uses the MFOLD 3.0 energy function for RNA @ T=37 given by:
!            
!            E = SUM E_loop + E_stack
!
! Arguments:
!
!          ISEQ - Array of length N containing the sequence
!                 in numerical code (A=1,C=2,G=3,U=4)
!          IBSP - Array of dimension (N) containing the information
!                 on base pairs in the RNA fold.
!                 IBSP(i) = j [i base pairs with j]
!                 IBSP(i) = 0 [i is single stranded]
!             N - Number of nucleotides in the sequence.
!             E - (OUTPUT) MFOLD 3.0 energy of the RNA fold.
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
! Subroutines - ELOOP, ESTACK
!
! Author(s): Eric Dykeman
!            Copyright (c) 2015 (Please Refer to LICENCE)
!
! ==============================================================================

      SUBROUTINE ESTRUC (ISEQ,IBSP,N,E)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(IN) :: iseq(n),ibsp(n)

        REAL, INTENT(OUT) :: e

        !=== VARIABLES ===!

        INTEGER :: i,j,ip,jp
        INTEGER :: il,nl,loop(n)

        REAL :: el,es


        e  = 0.0e0
        el = 0.0e0
        es = 0.0e0

        !=== Find Loops ===!

        nl = 1
        loop(1) = n

        DO i=1,n

          j = ibsp(i)

          IF ( j > i ) THEN

            ip = i + 1
            jp = j - 1

            IF ( ibsp(ip) /= jp ) THEN
              nl = nl + 1
              loop(nl) = i
            ENDIF

          ENDIF

        ENDDO
           write(99,*)'nl = ',nl
        !=== Compute Energy ===!

        DO il=1,nl

          i = loop(il)
          j = ibsp(i)

          IF ( i == n ) j = 1

          !=== Loop Energy ===!

          CALL ELOOP (iseq,ibsp,i,j,n,el)
          write(99,*)'eloop #',il,el
          e = e + el

          !=== Stacking Energy ===!

          ip = i - 1
          jp = j + 1

          IF ( ip < 1 ) CYCLE
          IF ( jp > n ) CYCLE
          IF ( i  > j ) CYCLE

          DO WHILE ( ibsp(ip) == jp )

            CALL ESTACK (iseq,ip,jp,i,j,n,es)

            e = e + es

            i = ip
            j = jp

            ip = ip - 1
            jp = jp + 1

            IF ( ip < 1 ) EXIT
            IF ( jp > n ) EXIT

          ENDDO

        ENDDO

        RETURN

      END SUBROUTINE ESTRUC
