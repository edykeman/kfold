! ==============================================================================
! Subroutine: ELOOP (ISEQ,IBSP,I,J,N,EL)
! 
! Purpose: Computes the energy of an RNA loop using the emperical
!          MFOLD 3.0 energy function.
!
! Method: Uses the MFOLD 3.0 energy function for RNA @ T=37 given by:
!            
!            EL = EHAIR                --> FOR k = 1
!            EL = EBULGE               --> FOR k = 2
!                                      --> FOR k > 2
!
!            EL = a + b*NS + c*NH + GS               IF NS <= 6
!               = a + b *6 + c*NH + GS + d*LOG(NS/6) IF NS  > 6
!
!            where k = number of helices in loop
!                  a = Multi-loop Penality
!                  b = Single Stranded Penalty
!                  c = Helix Penality
!                  d = 1.75 KT
!                  GS= Stacking Energy of Dangling Bases
!
! Arguments:
!
!          ISEQ - Array of length N containing the sequence
!                 in numerical code (A=1,C=2,G=3,U=4)
!          IBSP - Array of dimension (N) containing the information
!                 on base pairs in the RNA fold.
!                 IBSP(i) = j [i base pairs with j]
!                 IBSP(i) = 0 [i is single stranded]
!             I - Nucleotide position of the 5' most nucleotide.
!             J - Nucleotide position of the 3' most nucleotide.
!             N - Number of nucleotides in the sequence.
!            EL - (OUTPUT) MFOLD 3.0 energy of the RNA loop.
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
! Subroutines - EDANGLE
!
! Author(s): Eric Dykeman
!            Copyright (c) 2015 (Please Refer to LICENCE)
!
! ==============================================================================

      SUBROUTINE ELOOP (ISEQ,IBSP,I,J,N,EL)

        USE RNAvar, ONLY : em,eh,es,eaup,beta

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        INTEGER, INTENT(IN) :: i,j,n
        INTEGER, INTENT(IN) :: iseq(n),ibsp(n)

        REAL, INTENT(OUT) :: el

        !=== VARIABLES ===!

        INTEGER :: k,ip,jp,kp,ks,ke
        INTEGER :: nh,ns,iloop,ilast

        REAL :: x,c,ed,e3,e5


        el = 0.0e0

        nh = 0
        ns = 0

        c = 1.750e0 / REAL(beta)

        !=== Internal Loop iloop = 1 ===!
        !=== External Loop iloop = 0 ===!

        IF ( i < j ) iloop = 1
        IF ( i > j ) iloop = 0

        ks = MIN(i,j)
        ke = MAX(i,j)


        !=== Count number of Helices in Loop ===!

        k = ks

        DO WHILE ( k <= ke )

          IF ( ibsp(k) == 0 ) ns = ns + 1
          IF ( ibsp(k) >  k ) nh = nh + 1

          IF ( ibsp(k) > k ) THEN
          IF ( iloop == 0 .or. k /= ks ) THEN
            k = ibsp(k)
          ENDIF
          ENDIF

          k = k + 1

        ENDDO


        !=== Compute Loop Energy ===!

        IF ( nh == 1 .and. iloop == 1 ) THEN

          CALL EHAIR (iseq,i,j,n,el)

        ELSEIF ( nh == 2 .and. iloop == 1 ) THEN

          ip = i + 1
          DO WHILE ( ibsp(ip) == 0 )
          ip = ip + 1
          ENDDO

          jp = ibsp(ip)

          CALL EBULGE (iseq,i,j,ip,jp,n,el)

        ELSE

          e3 = 0.0e0

          IF ( iloop == 0 ) THEN

            k = ks
            ilast = 0

          ELSEIF ( iloop == 1 ) THEN

            ip = i
            jp = ibsp(i)

            k = ip + 1
            ilast = k

            IF ( k >= 1 .and. k <= n ) THEN
            IF ( ibsp(k) == 0 ) THEN
              CALL EDANGLE (iseq,ip,jp,k,n,e3)
            ENDIF
            ENDIF

            IF ( ns <= 6 ) THEN
              el = em + es * REAL(ns) + eh * REAL(nh)
            ELSE
              x = REAL(ns) / 6.0e0
              el = em + es * 6.0e0 + eh * REAL(nh)
              el = el + c * LOG(x)
            ENDIF

          ENDIF

          DO WHILE ( k <= ke )

            IF ( ibsp(k) /= 0 ) THEN

              ip = k
              jp = ibsp(k)

              kp = ip - 1

              e5 = 0.0e0

              IF ( kp >= 1 .and. kp <= n ) THEN
              IF ( ibsp(kp) == 0 ) THEN
                CALL EDANGLE (iseq,ip,jp,kp,n,e5)
              ENDIF
              ENDIF

              IF ( ilast == kp ) THEN
                ed = MIN(e3,e5)
              ELSE
                ed = e3 + e5
              ENDIF

              el = el + ed + eaup(iseq(ip),iseq(jp))

              kp = jp + 1
              ilast = kp

              e3 = 0.0e0

              IF ( kp >= 1 .and. kp <= n ) THEN
              IF ( ibsp(kp) == 0 ) THEN
                CALL EDANGLE (iseq,ip,jp,kp,n,e3)
              ENDIF
              ENDIF

              IF ( k /= ke ) THEN
                k = ibsp(k)
              ENDIF

            ENDIF

            k = k + 1

          ENDDO

          IF ( iloop == 0 ) THEN
            el = el + e3
          ENDIF

        ENDIF

        RETURN

      END SUBROUTINE ELOOP
