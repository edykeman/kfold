! ==============================================================================
! Module: CLASS_RNAFOLD
! 
! Purpose: A FORTRAN class structure containing subroutines and data
!          elements required for computing transition probabilites between
!          different RNA secondary structures.
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            01/01/2015   Original Code
!
! Dependancies:
!
! Modules - RNAVar
! Functions -
! Subroutines - 
!
! Author(s): Eric Dykeman
!            Copyright (c) 2015 (Please Refer to LICENCE)
!
! ==============================================================================

      MODULE CLASS_RNAFOLD

        USE RNAVar, ONLY : rateh,ratem,rated,beta,pnuc,&
                         & iwc,eaup,mxnt

        IMPLICIT NONE

        PRIVATE

        PUBLIC :: LOOP_INIT, LOOP_FIRE

        TYPE, PUBLIC :: RNA_STRUC

          CHARACTER :: seq(mxnt)

          INTEGER :: iseq(mxnt)
          INTEGER :: ibsp(mxnt)
          INTEGER :: link(mxnt)

          INTEGER :: loop(mxnt)
          INTEGER :: nhlx(mxnt)
          INTEGER :: nsgl(mxnt)

          INTEGER :: n
          INTEGER :: nl
          INTEGER :: nsum

          DOUBLE PRECISION :: psum(mxnt)
          DOUBLE PRECISION :: ptot(mxnt)

          DOUBLE PRECISION :: wrk1(mxnt)
          DOUBLE PRECISION :: wrk2(mxnt)

        END TYPE RNA_STRUC

        CONTAINS

        INCLUDE 'loop_init.f90'
        INCLUDE 'loop_resum.f90'

        INCLUDE 'helx_reac.f90'
        INCLUDE 'loop_reac.f90'
        INCLUDE 'loop_fire.f90'

      END MODULE CLASS_RNAFOLD
