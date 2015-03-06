! ==============================================================================
! Subroutine: LOOP_RESUM (R,INDX)
! 
! Purpose: Resums the partial sum table of transition rates and total
!          flux when the transition rate for a single loop element
!          (#indx) has changed.
!
! Method: Recomputes the total flux in LOG_2(N) time using a partial
!         sum table.
!
! Arguments:
!
!             R - Class structure containing information on the
!                 RNA secondary structure and possible reactions.
!          INDX - The indx number of the loop element.
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

      SUBROUTINE LOOP_RESUM (R,INDX)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE(RNA_STRUC), INTENT(INOUT) :: r

        INTEGER, INTENT(IN) :: indx

        !=== VARIABLES ===!

        INTEGER :: i,j,k,nsum
        INTEGER :: n,n1,n2


        nsum = r% nsum

        !=== Resum Partial Sum Table ===!

        n = 1
        n1= 2
        n2= 4

        IF ( MOD(indx,2) == 1 ) i = indx
        IF ( MOD(indx,2) == 0 ) i = indx - 1

        r% psum(i) = r% ptot(i) + r% ptot(i+1)

        DO WHILE ( n1 < nsum )

          i = INT(i/n2) * n2 + n1

          j = i - n
          k = i + n

          r% psum(i) = r% psum(j) + r% psum(k)

          n  = n1
          n1 = n2
          n2 = 2 * n2

        ENDDO

        RETURN

      END SUBROUTINE LOOP_RESUM
