! ==============================================================================
! Subroutine: V2CT (IBSP,FLD,JOB,N)
! 
! Purpose: Converts a secondary structure from either Vienna format
!          to CT format, or CT to Vienna.
!
! Method:
!
! Arguments:
!
!          IBSP - Base-pair information for the RNA.
!           FLD - Vienna fold for the RNA.
!           JOB - Conversion job to perform.
!                 JOB = 'V' Convert to Vienna Format
!                 JOB = 'C' Convert to CT Format.
!             N - Number of nucleotides in the RNA.
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

      SUBROUTINE V2CT (IBSP,FLD,JOB,N)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(INOUT) :: ibsp(n)

        CHARACTER, INTENT(IN) :: job
        CHARACTER, INTENT(INOUT) :: fld(n)


        !=== VARIABLES ===!

        INTEGER :: i,j,is
        INTEGER :: istack(n)


        !=== Convert to Vienna ===!

        IF ( job == 'V' ) THEN

          DO i=1,n
          IF ( ibsp(i) == 0 ) THEN

            fld(i) = '.'

          ELSE

            j = ibsp(i)

            IF ( j > i ) THEN
              fld(i) = '('
              fld(j) = ')'
            ELSE
              fld(i) = ')'
              fld(j) = '('
            ENDIF

          ENDIF
          ENDDO

        ENDIF

        !=== Convert to CT ===!

        IF ( job == 'C' ) THEN

          is = 0

          ibsp(:) = 0
          istack(:) = 0

          DO i=1,n
          IF ( fld(i) == '.' ) THEN

            ibsp(i) = 0

          ELSEIF ( fld(i) == '(' ) THEN

            is = is + 1
            istack(is) = i

          ELSEIF ( fld(i) == ')' ) THEN

            j = istack(is)

            ibsp(j) = i
            ibsp(i) = j

            is = is - 1

          ENDIF
          ENDDO

        ENDIF

        RETURN

      END SUBROUTINE V2CT
