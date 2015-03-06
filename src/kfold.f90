! ==============================================================================
! Program: KFOLD
! 
! Description: A program for computing the folding kinetics of an RNA
!              sequence using Turner energies.
!
!              Please See/Cite:
!              Dykeman,E.C.
!
! Notes:
!        EXECUTION PATH
!
!        Section 1 - Read In Files and Prepare
!        Section 2 - Perform the simulation
!        Section 3 - Clean Up and Finish Program
!        END
!
!        FILE TREE
!
!        Unit = 1  --- Sequence File
!        Unit = 2  --- Trajectory Output File
!        Unit = 3  --- Log File
!        Unit = 5  --- Standard Out (NOT USED)
!        Unit = 6  --- Standard In (NOT USED)
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            01/01/2015   Original Code
!
! Dependancies:
!
! Modules - Class_RNAFold, RNAVar
! Functions -
! Subroutines - READDATA SSAREACTION V2CT
!
! Author(s): Eric Dykeman
!            Copyright (c) 2015 (Please Refer to LICENCE)
!
! ==============================================================================

      PROGRAM KFOLD

        USE RNAVar, ONLY : mxnt

        USE Class_RNAFold

        IMPLICIT NONE

        !=== VARIABLES ===!

        TYPE(RNA_STRUC) :: rna

        CHARACTER :: seq(mxnt)
        CHARACTER :: fld(mxnt)

        INTEGER :: iseq(mxnt)
        INTEGER :: ibpi(mxnt)
        INTEGER :: ibpf(mxnt)

        CHARACTER (LEN=mxnt+1) :: fasta

        CHARACTER (LEN=70) :: seqfile,outfile,logfile
        CHARACTER (LEN=70) :: xdata,arg

        INTEGER :: i,j,k,n,nn,is,io
        INTEGER :: isim,nsim,narg,iseed

        DOUBLE PRECISION :: random,tstart,time
        DOUBLE PRECISION :: tout,tmax,dt

        LOGICAL :: istart,istop


        !=== DEFAULT SETTINGS ===!

        !=== Max Time (microseconds) ===!

        tstart = 0.0d0
        tmax = 1.0d0

        nsim = 1
        iseed = 61928712

        !=== Inital / Final Structure ===!

        ibpi(:) = 0
        ibpf(:) = 0

        istart= .FALSE.
        istop = .FALSE.

        seqfile = 'seq.fasta'
        outfile = 'seq.traj'
        logfile = 'seq.log'


        !=== SECTION 0 - Input From Command Line ===!

        narg = IARGC ()

        DO i=1,narg,2

          CALL GETARG (i,arg)

          SELECT CASE (arg)

            CASE ('-i')
              CALL GETARG (i+1,seqfile)
            CASE ('-o')
              CALL GETARG (i+1,outfile)
            CASE ('-l')
              CALL GETARG (i+1,logfile)
            CASE ('-n')
              CALL GETARG (i+1,xdata)
              READ(xdata,*)nsim
            CASE ('-s')
              CALL GETARG (i+1,xdata)
              READ(xdata,*)iseed
            CASE ('-t')
              CALL GETARG (i+1,xdata)
              READ(xdata,*)tmax
            CASE DEFAULT

              WRITE(*,*)arg,'Invalid Line Argument'
              STOP

          END SELECT

        ENDDO


        !=== SECTION 1 - Read in data ===!

        CALL READDATA

        OPEN (UNIT=1, FILE=seqfile, STATUS='Unknown')
        OPEN (UNIT=2, FILE=outfile, STATUS='Unknown')
        OPEN (UNIT=3, FILE=logfile, STATUS='Unknown')

        READ(1,*)fasta

        nn = LEN_TRIM(fasta)

        IF ( nn > mxnt ) THEN
          WRITE(*,*)'ERROR: Maximum number of nt = ',mxnt
          WRITE(*,*)'Increase mxnt in rnavar.f90'
          STOP
        ENDIF

        READ(fasta,'(10000A1)')(seq(k),k=1,nn)

        READ(1,*,IOSTAT=io)fasta

        IF ( io == 0 ) THEN

          istart = .TRUE.

          !=== Start Structure Specified ===!

          READ(fasta,'(10000A1)')(fld(k),k=1,nn)

          !=== Convert Vienna to CT ===!

          CALL V2CT (ibpi,fld,'C',nn)

        ENDIF

        READ(1,*,IOSTAT=io)fasta

        IF ( io == 0 ) THEN

          istop = .TRUE.

          !=== Stop Structure Specified ===!

          READ(fasta,'(10000A1)')(fld(k),k=1,nn)

          !=== Convert Vienna to BP ===!

          CALL V2CT (ibpf,fld,'C',nn)

        ENDIF

        !=== Setup RNA ===!

        CALL CONVERT (seq,iseq,nn)
        CALL SETUPNUC (nn)

        rna% seq(:) = seq(:)
        rna% iseq(:) = iseq(:)
        rna% n = nn


        !=== SECTION 2 - Perform RNA Kinetics ===! 

        DO isim=1,nsim

          io = 1
          dt = 1.0d-2
          time = tstart

          IF ( istart ) THEN
            rna% ibsp(:) = ibpi(:)
          ELSE
            rna% ibsp(:) = 0
          ENDIF
 
          CALL LOOP_INIT (rna)

          !=== STOCHASTIC SIMULATION ===!

          DO WHILE ( time < tmax )

            CALL SSAREACTION (rna,iseed,time,tout)

            !=== Increment tout ===!

            IF ( time > tout ) THEN

              tout = tout + dt

              io = io + 1

              IF ( io > 9 ) THEN
                io = 1
                dt = dt * 10.0d0
              ENDIF

            ENDIF

            !=== Check for stop structure ===!

            IF ( istop ) THEN

              j = 0

              DO i=1,nn
              IF ( ibpf(i) == rna%ibsp(i) ) j = j + 1
              ENDDO

              IF ( j == nn ) THEN
                WRITE(3,'(E16.8)')time
                EXIT
              ENDIF

            ENDIF

          ENDDO

        ENDDO

        CLOSE (UNIT=1)
        CLOSE (UNIT=2)
        CLOSE (UNIT=3)

      END PROGRAM KFOLD
