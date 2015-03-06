! ==============================================================================
! Subroutine: EBULGE (ISEQ,I,J,IP,JP,N,EB)
! 
! Purpose: Computes the energy of an RNA buldge with two helices using the
!          empirical MFOLD 3.0 energy function.
!
! Method: Uses the MFOLD 3.0 energy function for RNA @ T=37 given by:
!
!         EB = E_entropic + E_stack + E_stack + E_asymmetry
!
!               5' (I) X ... W (IP) 3'
!               3' (J) Y ... Z (JP) 5'
!
!               NOTE: I < IP and JP < J
!
! Arguments:
!
!          ISEQ - Array of length N containing the sequence
!                 in numerical code (A=1,C=2,G=3,U=4)
!             I - Nucleotide position of the starting basepair 5'.
!             J - Nucleotide position of the starting basepair 3'.
!            IP - Nucleotide position of the ending basepair 3'.
!            JP - Nucleotide position of the ending basepair 5'.
!             N - Number of nucleotides in the sequence.
!            EB - (OUTPUT) MFOLD 3.0 bulge energy of the sequence.
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
! Subroutines - TSTACK, TSTACKI, TINT11, TINT12, TINT22
!
! Author(s): Eric Dykeman
!            Copyright (c) 2015 (Please Refer to LICENCE)
!
! ==============================================================================

      SUBROUTINE EBULGE (ISEQ,I,J,IP,JP,N,EB)

        USE RNAvar, ONLY : eau,beta

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        INTEGER, INTENT(IN) :: i,j,ip,jp,n
        INTEGER, INTENT(IN) :: iseq(n)

        REAL, INTENT(OUT) :: eb

        !=== VARIABLES ===!

        INTEGER :: k,ibul,imin,imax
        INTEGER :: n1,n2,nt,na,list(8)

        REAL :: x,c,f(4)
        REAL :: elb(30),eli(30)

        DATA (  f(k),k=1,4 ) / 0.50e0,0.50e0,0.50e0,0.50e0 /

        DATA (elb(k),k=1,30) / 3.80e0,2.80e0,3.20e0,3.60e0,4.00e0,&
                             & 4.40e0,4.60e0,4.70e0,4.80e0,4.90e0,&
                             & 5.00e0,5.10e0,5.20e0,5.30e0,5.40e0,&
                             & 5.40e0,5.50e0,5.50e0,5.60e0,5.70e0,&
                             & 5.70e0,5.80e0,5.80e0,5.80e0,5.90e0,&
                             & 5.90e0,6.00e0,6.00e0,6.00e0,6.10e0 /

        DATA (eli(k),k=1,30) / 0.00e0,0.00e0,0.00e0,1.70e0,1.80e0,&
                             & 2.00e0,2.20e0,2.30e0,2.40e0,2.50e0,&
                             & 2.60e0,2.70e0,2.80e0,2.90e0,3.00e0,&
                             & 3.00e0,3.10e0,3.10e0,3.20e0,3.30e0,&
                             & 3.30e0,3.40e0,3.40e0,3.40e0,3.50e0,&
                             & 3.50e0,3.60e0,3.60e0,3.60e0,3.70e0 /


        eb = 0.0e0

        n1 = ip - i - 1
        n2 = j - jp - 1

        nt = n1 + n2
        na = IABS(n1-n2)

        imin = MIN(n1,n2)
        imax = MAX(n1,n2)

        c = 1.750e0 / REAL(beta)


        !=== Get Buldge Type ===!

        ibul = 6
        IF ( n1 == 0 .and. n2 == 1 ) ibul = 0
        IF ( n1 == 1 .and. n2 == 0 ) ibul = 0
        IF ( n1 == 0 .and. n2 >= 2 ) ibul = 1
        IF ( n1 >= 2 .and. n2 == 0 ) ibul = 1
        IF ( n1 == 1 .and. n2 == 1 ) ibul = 2
        IF ( n1 == 1 .and. n2 == 2 ) ibul = 3
        IF ( n1 == 2 .and. n2 == 1 ) ibul = 4
        IF ( n1 == 2 .and. n2 == 2 ) ibul = 5


        !=== TERM 1 --> Entropic Term ===!

        IF ( nt > 30 ) THEN

          x = REAL(nt) / 30.0d0
          x = c * LOG(x)

          IF ( ibul <= 1 ) eb = elb(30) + x
          IF ( ibul == 6 ) eb = eli(30) + x

        ELSE

          IF ( ibul <= 1 ) eb = elb(nt)
          IF ( ibul == 6 ) eb = eli(nt)

        ENDIF


        !=== TERM 2 & 3 --> Stacking Energy ===!

        SELECT CASE (ibul)

          CASE (0)

            ! 5' (i) A . X (ip) 3'
            ! 3' (j) U   Y (jp) 5'

            list(1) = iseq(i)
            list(2) = iseq(j)
            list(3) = iseq(ip)
            list(4) = iseq(jp)

            CALL TSTACK (list,eb)

          CASE (1)

            ! 5' (i) A .. X (ip) 3'
            ! 3' (j) U    Y (jp) 5'

            !=== Closing A-U / G-U Penalty ===!

            IF ( iseq(i) ==  4 ) eb = eb + eau
            IF ( iseq(j) ==  4 ) eb = eb + eau
            IF ( iseq(ip) == 4 ) eb = eb + eau
            IF ( iseq(jp) == 4 ) eb = eb + eau

          CASE (2)

            ! 5' (i) A . X (ip) 3'
            ! 3' (j) U . Y (jp) 5'

            list(1) = iseq(i)
            list(2) = iseq(j)
            list(3) = iseq(i+1)
            list(4) = iseq(j-1)
            list(5) = iseq(ip)
            list(6) = iseq(jp)

            CALL TINT11 (list,eb)

          CASE (3)

            ! 5' (i) A .   X (ip) 3'
            ! 3' (j) U . . Y (jp) 5'

            list(1) = iseq(i)
            list(2) = iseq(j)
            list(3) = iseq(i+1)
            list(4) = iseq(j-1)
            list(5) = iseq(j-2)
            list(6) = iseq(ip)
            list(7) = iseq(jp)

            CALL TINT12 (list,eb)

          CASE (4)

            ! 5' (i) A . . X (ip) 3'
            ! 3' (j) U   . Y (jp) 5'

            list(1) = iseq(jp)
            list(2) = iseq(ip)
            list(3) = iseq(jp+1)
            list(4) = iseq(ip-1)
            list(5) = iseq(ip-2)
            list(6) = iseq(j)
            list(7) = iseq(i)

            CALL TINT12 (list,eb)

          CASE (5)

            ! 5' (i) A . . X (ip) 3'
            ! 3' (j) U . . Y (jp) 5'

            list(1) = iseq(i)
            list(2) = iseq(j)
            list(3) = iseq(i+1)
            list(4) = iseq(j-1)
            list(5) = iseq(ip-1)
            list(6) = iseq(jp+1)
            list(7) = iseq(ip)
            list(8) = iseq(jp)

            CALL TINT22 (list,eb)

          CASE (6)

            ! 5' (i) A X .. G (ip) 3'
            ! 3' (j) U Y .. C (jp) 5'

            list(1) = iseq(i)
            list(2) = iseq(j)
            list(3) = iseq(i+1)
            list(4) = iseq(j-1)

            !=== GAIL Rule ===!

            IF ( imin == 1 .and. imax > 2 ) THEN

              list(3) = 1
              list(4) = 1

            ENDIF

            CALL TSTACKI (list,eb)

            ! 5' (i) A .. X G (ip) 3'
            ! 3' (j) U .. Y C (jp) 5'

            list(1) = iseq(jp)
            list(2) = iseq(ip)
            list(3) = iseq(jp+1)
            list(4) = iseq(ip-1)

            !=== GAIL Rule ===!

            IF ( imin == 1 .and. imax > 2 ) THEN

              list(3) = 1
              list(4) = 1
            
            ENDIF

            CALL TSTACKI (list,eb)

          CASE DEFAULT

            WRITE(*,*)'ERROR: ibul out of range [0-6]'
            STOP

        END SELECT


        !=== PART 4 ---> Asymmetry Penalty ===!

        IF ( ibul == 6 ) THEN

          k = MIN(4,n1,n2)

          x = REAL(na) * f(k)

          x = MIN(x,3.0d0)

          eb = eb + x

        ENDIF

        RETURN

      END SUBROUTINE EBULGE
