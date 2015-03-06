! ==============================================================================
! Subroutine: LOOP_INIT (R)
! 
! Purpose: Initializes the data structures (loop elements) required
!          for RNA kinetics.
!
! Method:
!
! Arguments:
!
!             R - Class structure containing information on the
!                 RNA secondary structure and possible reactions.
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

      SUBROUTINE LOOP_INIT (R)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE(RNA_STRUC), INTENT(INOUT) :: r

        !=== VARIABLES ===!

        INTEGER :: i,j,n,nl,ns,nh
        INTEGER :: ip,jp,kp,ks,ke
        INTEGER :: nsum


        !=== Initialize RNA Data ===!

        n = r% n

        r% link(:) = 0
        r% loop(:) = 0
        r% nhlx(:) = 0
        r% nsgl(:) = 0

        r% wrk1(:) = 0.0d0
        r% wrk2(:) = 0.0d0

        r% psum(:) = 0.0d0
        r% ptot(:) = 0.0d0


        !=== Find Loops ===!

        nl = 1
        r% loop(1) = n

        DO i=1,n

          j = r% ibsp(i)

          IF ( j > i ) THEN

            ip = i + 1
            jp = j - 1

            IF ( r% ibsp(ip) /= jp ) THEN
              nl = nl + 1
              r% loop(nl) = i
            ENDIF

          ENDIF

        ENDDO

        r% nl = nl


        !=== Make Links ===!

        DO i=1,nl

          ip = r% loop(i)
          jp = r% ibsp(ip)

          IF ( ip == n ) jp = 1

          r% link(ip) = i

          IF ( ip < jp ) THEN
            ks = ip + 1
            ke = jp
            nh = 1
            ns = 0
          ELSE
            ks = jp
            ke = ip
            nh = 0
            ns = 0
          ENDIF

          kp = ks

          DO WHILE ( kp <= ke )

            IF ( r% ibsp(kp) > kp ) nh = nh + 1
            IF ( r% ibsp(kp) == 0 ) ns = ns + 1

            IF ( r% ibsp(kp) > kp ) THEN
              kp = r% ibsp(kp)
              r% link(kp) = i
            ELSE
              kp = kp + 1
            ENDIF

          ENDDO

          r% nhlx(i) = nh
          r% nsgl(i) = ns

        ENDDO


        !=== Compute Size of Partial Sum Table ===!

        nsum = 2
        DO WHILE ( nsum < nl )
        nsum = 2 * nsum
        ENDDO

        r% nsum = nsum


        !=== Compute Reactions for Loops ===!

        DO i=1,nl

          CALL LOOP_REAC (r,i)

        ENDDO

        RETURN

      END SUBROUTINE LOOP_INIT
