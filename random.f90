! ==============================================================================
! Function: RANDOM (ISEED)
! 
! Purpose: Generates a random number on the interval [0,1] given an
!          integer seed.
!
! Method: Efficent Fortran Programming, John Wiley and Sons, New York (1990)
!         pp. 17-18 Library of Congress code QA76.73.F25 K78
!
! Arguments:
!
!           ISEED - Four byte integer seed.
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

      DOUBLE PRECISION FUNCTION RANDOM (ISEED)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        INTEGER, INTENT(INOUT) :: iseed

        !=== VARIABLES ===!

        INTEGER :: hi,lo,test

        INTEGER, PARAMETER :: a = 16807
        INTEGER, PARAMETER :: m = 2147483647
        INTEGER, PARAMETER :: q = 127773
        INTEGER, PARAMETER :: r = 2836


        hi = INT(iseed/q)
        lo = MODULO(iseed,q)

        test = a * lo - r * hi

        IF ( test > 0 ) THEN

          iseed = test

        ELSE

          iseed = test + m

        ENDIF

        RANDOM = DBLE(iseed) / DBLE(m)

        RETURN

      END FUNCTION RANDOM
