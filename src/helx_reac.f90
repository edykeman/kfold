! ==============================================================================
! Subroutine: HELX_REAC (R,I)
! 
! Purpose: Recomputes the reaction rates for opening one of the internal
!          base-pairs within a continuous section of helix.
!
! Method:
!
!               A-U           A-U
!               x-x   --->   x   x
!               U-A           U-A
!
! Arguments:
!
!             R - Class structure containing information on the
!                 RNA secondary structure and possible reactions.
!             I - The indx number of one of the nucelotides which
!                 belongs to the helix.
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

      SUBROUTINE HELX_REAC (R,I)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE(RNA_STRUC), INTENT(INOUT) :: r

        INTEGER, INTENT(IN) :: i

        !=== VARIABLES ===!

        INTEGER :: ke,ip,jp,indx

        DOUBLE PRECISION :: x,dg,atot


        jp = i

        ip = r% ibsp(jp)
        jp = MIN(ip,jp)

        DO WHILE ( r% link(jp) == 0 )

          jp = jp + 1

        ENDDO

        indx = r% link(jp)


        !=== Open BP Inside Helix ===!

        ke = r% ibsp(jp)
        ip = r% ibsp(jp)

        atot = r% wrk1(jp)

        IF ( r% link(ke) == 0 ) THEN

          ip = ip + 1
          jp = jp - 1

          DO WHILE ( r% link(ip) == 0 )

            CALL DELTAG_HI (r,ip,jp,dg)

            dg = dg / 2.0d0

            x = beta * dg
            x = DEXP(-x) * rateh

            r% wrk1(ip) = x
            r% wrk1(jp) = 0.0d0

            r% wrk2(ip) = 0.0d0
            r% wrk2(jp) = 0.0d0

            atot = atot + x

            ip = ip + 1
            jp = jp - 1

          ENDDO

        ENDIF

        r% ptot(indx) = atot

        CALL LOOP_RESUM (r,indx)

        RETURN

      END SUBROUTINE HELX_REAC
