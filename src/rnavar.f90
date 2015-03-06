! ==============================================================================
! Module: RNAVAR
! 
! Purpose: Contains the variables, work arrays, and parameters needed for 
!          computing the kinetics of folding for an RNA molecule.
!
!          For a list of parameters see:
!
!          A. M. Zuker, "Algorithms and thermodynamics for RNA secondary
!          structure predicitions: A practical guide".
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            01/01/2015   Original Code
!
! Contains:
!
! Modules -
! Functions -
! Subroutines -
!
! Author(s): Eric Dykeman
!            Copyright (c) 2015 (Please Refer to LICENCE)
!
! ==============================================================================

      MODULE RNAVAR

        IMPLICIT NONE

        !=== ALLOCATABLE ARRAYS ===!

        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: pnuc

        !=== VARIABLES ===!

        DOUBLE PRECISION :: beta = 0.16225023135094183147d1

        INTEGER :: iwc(4,4)

        REAL :: eaup(4,4)

        !=== PARAMETERS ===!

        DOUBLE PRECISION, PARAMETER :: gcons = 1.987206d-3
        DOUBLE PRECISION, PARAMETER :: temp = 310.15d0

        DOUBLE PRECISION, PARAMETER :: rateh = 1.0d2 !Helix Extension (1/uS)!
        DOUBLE PRECISION, PARAMETER :: ratem = 5.0d0 !Helix Moprhing  (1/uS)!
        DOUBLE PRECISION, PARAMETER :: rated = 5.0d0 !Helix Diffusion (1/uS)!

        INTEGER, PARAMETER :: mxnt = 10000

        REAL, PARAMETER :: em  = 10.10e0
        REAL, PARAMETER :: eh  = -0.30e0
        REAL, PARAMETER :: es  = -0.30e0
        REAL, PARAMETER :: eau = +0.50e0

      END MODULE
