! ==============================================================================
! Subroutine: DELTAG_HD (R,II,JJ,KK,DG)
! 
! Purpose: Computes the difference in free energy of an RNA loop due
!          to a nucleotide diffusion along the helix i.e. the ii-jj
!          base pair will shift to either ii-kk or kk-jj. The energy
!          change is calculated using the emperical MFOLD 3.0 energy
!          function.
!
! Method: Uses the MFOLD 3.0 energy function for RNA @ T=37.
!
! Arguments:
!
!             R - Class containing information about the RNA fold and
!                 the RNA sequence.
!            II - Nucleotide position of the 5' most nucleotide.
!            JJ - Nucleotide position of the 3' most nucleotide.
!            KK - Nucleotide position of the single-stranded nucleotiode
!                 that either ii/jj in the base-pair ii-jj will swap with.
!            DG - (OUTPUT) The energy difference from diffusion of kk.
!
! History:
!
! Version    Date         Comment
! --------   ----------   -----------------------
!            01/01/2015   Original Code
!
! Dependencies:
!
! Modules - Class_RNAFold, RNAVar
! Functions -
! Subroutines - EHAIR, EBULGE, ESTACK, EDANGLE
!
! Author(s): Eric Dykeman
!            Copyright (c) 2015 (Please Refer to LICENCE)
!
! ==============================================================================

      SUBROUTINE DELTAG_HD (R,II,JJ,KK,DG)

        USE RNAvar, ONLY : em,eh,es,eaup,beta

        USE Class_RNAFold

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE(RNA_STRUC), INTENT(IN) :: r
        INTEGER, INTENT(IN) :: ii,jj,kk

        DOUBLE PRECISION, INTENT(OUT) :: dg

        !=== VARIABLES ===!

        INTEGER :: i,j,k,l,n,ip,jp
        INTEGER :: nh,ns,mh,ms,is,js
        INTEGER :: indx,iloop,kloop

        REAL :: e1,e2,e3,e4,e5,ed
        REAL :: x,c,ei,ef


        indx = r% link(jj)

        i = r% loop(indx)
        j = r% ibsp(i)
        n = r% n

        nh = r% nhlx(indx)
        ns = r% nsgl(indx)

        c = 1.750e0 / REAL(beta)

        !=== Internal Loop iloop = 1 ===!
        !=== External Loop iloop = 0 ===!

        iloop = 1
        kloop = 1

        IF ( i > j ) THEN
          iloop = 0
          i = 1
          j = n
        ENDIF

        IF ( r% link(ii) /= 0 ) THEN

          k = r% link(ii)

          mh = r% nhlx(k)
          ms = r% nsgl(k)

          k = r% loop(k)
          l = r% ibsp(k)

          IF ( k > l ) THEN
            kloop = 0
            k = 1
            l = n
          ENDIF

        ELSE

          mh = 2
          ms = 0

          k = ii
          l = jj

        ENDIF


        e4 = 0.0e0

        !=== INITIAL ENERGY ===!

        IF ( nh == 1 .and. iloop == 1 ) THEN

          CALL EHAIR (r%iseq,i,j,n,e1)

        ELSEIF ( nh == 2 .and. iloop == 1 ) THEN

          ip = i + 1
          DO WHILE ( r% ibsp(ip) == 0 )
          ip = ip + 1
          ENDDO

          jp = r% ibsp(ip)

          CALL EBULGE (r%iseq,i,j,ip,jp,n,e1)

        ELSE

          e1 = 0.0e0
          e5 = 0.0e0
          e3 = 0.0e0

          is = r% iseq(ii)
          js = r% iseq(jj)

          IF ( iloop == 1 ) THEN

            IF ( ns <= 6 ) THEN
              e1 = em + es * REAL(ns) + eh * REAL(nh)
            ELSE
              x = REAL(ns) / 6.0e0
              e1 = em + es * 6.0e0 + eh * REAL(nh)
              e1 = e1 + c * LOG(x)
            ENDIF

          ENDIF

          IF ( ii > 1 ) THEN
          IF ( r% ibsp(ii-1) == 0 ) THEN

            CALL EDANGLE (r%iseq,ii,jj,ii-1,n,e5)

            IF ( ii > 2 ) THEN

              ip = ii - 2
              jp = r% ibsp(ip)

              IF ( jp /= 0 ) THEN
                CALL EDANGLE (r%iseq,ip,jp,ii-1,n,ed)
                IF ( kk == ii+1 ) e4 = e4 + ed
                e5 = MIN(e5,ed) 
              ENDIF

            ENDIF

            IF ( ii > 3 .and. kk == ii-1 ) THEN
            IF ( r% ibsp(ii-2) == 0 ) THEN

              ip = ii - 3
              jp = r% ibsp(ip)

              IF ( jp /= 0 ) THEN
                CALL EDANGLE (r%iseq,ip,jp,ii-2,n,ed)
                e1 = e1 + ed
              ENDIF

            ENDIF
            ENDIF

          ENDIF
          ENDIF

          IF ( jj < n ) THEN
          IF ( r% ibsp(jj+1) == 0 ) THEN

            CALL EDANGLE (r%iseq,ii,jj,jj+1,n,e3)

            IF ( jj < n-1 ) THEN

              ip = jj + 2
              jp = r% ibsp(ip)

              IF ( jp /= 0 ) THEN
                CALL EDANGLE (r%iseq,ip,jp,jj+1,n,ed)
                IF ( kk == jj-1 ) e4 = e4 + ed
                e3 = MIN(e3,ed)
              ENDIF

            ENDIF

            IF ( jj < n-2 .and. kk == jj+1 ) THEN
            IF ( r% ibsp(jj+2) == 0 ) THEN

              ip = jj + 3
              jp = r% ibsp(ip)

              IF ( jp /= 0 ) THEN
                CALL EDANGLE (r%iseq,ip,jp,jj+2,n,ed)
                e1 = e1 + ed
              ENDIF

            ENDIF
            ENDIF

          ENDIF
          ENDIF

          e1 = e1 + e3 + e5 + eaup(is,js)

        ENDIF


        IF ( mh == 1 .and. kloop == 1 ) THEN

          CALL EHAIR (r%iseq,ii,jj,n,e2)

        ELSEIF ( mh == 2 .and. kloop == 1 ) THEN

          IF ( ms == 0 ) THEN

            CALL ESTACK (r%iseq,ii,jj,ii+1,jj-1,n,e2)

          ELSE

            ip = ii + 1
            DO WHILE ( r% ibsp(ip) == 0 )
            ip = ip + 1
            ENDDO

            jp = r% ibsp(ip)

            CALL EBULGE (r%iseq,ii,jj,ip,jp,n,e2)

          ENDIF

        ELSE

          e2 = 0.0e0
          e5 = 0.0e0
          e3 = 0.0e0

          is = r% iseq(ii)
          js = r% iseq(jj)

          IF ( kloop == 1 ) THEN

            IF ( ms <= 6 ) THEN
              e2 = em + es * REAL(ms) + eh * REAL(mh)
            ELSE
              x = REAL(ms) / 6.0e0
              e2 = em + es * 6.0e0 + eh * REAL(mh)
              e2 = e2 + c * LOG(x)
            ENDIF

          ENDIF

          IF ( jj > 1 ) THEN
          IF ( r% ibsp(jj-1) == 0 ) THEN

            CALL EDANGLE (r%iseq,ii,jj,jj-1,n,e5)

            IF ( jj > 2 ) THEN

              ip = jj - 2
              jp = r% ibsp(ip)

              IF ( jp /= 0 ) THEN
                CALL EDANGLE (r%iseq,ip,jp,jj-1,n,ed)
                IF ( kk == jj+1 ) e4 = e4 + ed
                e5 = MIN(e5,ed)
              ENDIF

            ENDIF

            IF ( jj > 3 .and. kk == jj-1 ) THEN
            IF ( r% ibsp(jj-2) == 0 ) THEN

              ip = jj - 3
              jp = r% ibsp(ip)

              IF ( jp /= 0 ) THEN
                CALL EDANGLE (r%iseq,ip,jp,jj-2,n,ed)
                e2 = e2 + ed
              ENDIF

            ENDIF
            ENDIF

          ENDIF
          ENDIF

          IF ( ii < n ) THEN
          IF ( r% ibsp(ii+1) == 0 ) THEN

            CALL EDANGLE (r%iseq,ii,jj,ii+1,n,e3)

            IF ( ii < n-1 ) THEN

              ip = ii + 2
              jp = r% ibsp(ip)

              IF ( jp /= 0 ) THEN
                CALL EDANGLE (r%iseq,ip,jp,ii+1,n,ed)
                IF ( kk == ii-1 ) e4 = e4 + ed
                e3 = MIN(e3,ed)
              ENDIF

            ENDIF

            IF ( ii < n-2 .and. kk == ii+1 ) THEN
            IF ( r% ibsp(ii+2) == 0 ) THEN

              ip = ii + 3
              jp = r% ibsp(ip)

              IF ( jp /= 0 ) THEN
                CALL EDANGLE (r%iseq,ip,jp,ii+2,n,ed)
                e2 = e2 + ed
              ENDIF

            ENDIF
            ENDIF

          ENDIF
          ENDIF

          e2 = e2 + e3 + e5 + eaup(is,js)

        ENDIF

        ei = e1 + e2


        !=== FINAL ENERGY ===!

        IF ( kk == ii-1 .or. kk == jj+1 ) THEN
          ns = ns - 1
          ms = ms + 1
        ELSE
          ns = ns + 1
          ms = ms - 1
        ENDIF

        IF ( nh == 1 .and. iloop == 1 ) THEN

          IF ( kk == i-1 .or. kk == i+1 ) THEN
            CALL EHAIR (r%iseq,kk,j,n,e1)
          ENDIF

          IF ( kk == j-1 .or. kk == j+1 ) THEN
            CALL EHAIR (r%iseq,i,kk,n,e1)
          ENDIF

        ELSEIF ( nh == 2 .and. iloop == 1 ) THEN

          IF ( ns == 0 ) THEN

            IF ( kk == ii-1 ) THEN
              CALL ESTACK (r%iseq,kk-1,jj+1,kk,jj,n,e1)
            ELSE
              CALL ESTACK (r%iseq,ii-1,kk+1,ii,kk,n,e1)
            ENDIF

          ELSE

            ip = i + 1
            DO WHILE ( r% ibsp(ip) == 0 )
            ip = ip + 1
            ENDDO

            jp = r% ibsp(ip)

            IF ( ii == j ) THEN
            IF ( kk == ii+1 .or. kk == ii-1 ) THEN
              CALL EBULGE (r%iseq,jj,kk,ip,jp,n,e1)
            ELSE
              CALL EBULGE (r%iseq,kk,ii,ip,jp,n,e1)
            ENDIF
            ELSE
            IF ( kk == ii+1 .or. kk == ii-1 ) THEN
              CALL EBULGE (r%iseq,i,j,kk,jj,n,e1)
            ELSE
              CALL EBULGE (r%iseq,i,j,ii,kk,n,e1)
            ENDIF
            ENDIF

          ENDIF

        ELSE

          e1 = 0.0e0
          e5 = 0.0e0
          e3 = 0.0e0

          IF ( iloop == 1 ) THEN

            IF ( ns <= 6 ) THEN
              e1 = em + es * REAL(ns) + eh * REAL(nh)
            ELSE
              x = REAL(ns) / 6.0e0
              e1 = em + es * 6.0e0 + eh * REAL(nh)
              e1 = e1 + c * LOG(x)
            ENDIF

          ENDIF

          IF ( kk == jj+1 .or. kk == jj-1 ) THEN

            is = r% iseq(ii)
            js = r% iseq(kk)

            IF ( ii > 1 ) THEN
            IF ( r% ibsp(ii-1) == 0 ) THEN

              CALL EDANGLE (r%iseq,ii,kk,ii-1,n,e5)

              IF ( ii > 2 ) THEN

                ip = ii - 2
                jp = r% ibsp(ip)

                IF ( jp /= 0 ) THEN
                  CALL EDANGLE (r%iseq,ip,jp,ii-1,n,ed)
                  e5 = MIN(e5,ed)
                ENDIF

              ENDIF

            ENDIF
            ENDIF

            IF ( kk < n ) THEN
            IF ( r% ibsp(kk+1) == 0 .or. kk == jj-1 ) THEN

              CALL EDANGLE (r%iseq,ii,kk,kk+1,n,e3)

              IF ( kk < n-1 ) THEN

                ip = kk + 2
                jp = r% ibsp(ip)

                IF ( jp /= 0 ) THEN
                  CALL EDANGLE (r%iseq,ip,jp,kk+1,n,ed)
                  e3 = MIN(e3,ed)
                ENDIF

              ENDIF

            ENDIF
            ENDIF

          ENDIF

          IF ( kk == ii-1 .or. kk == ii+1 ) THEN

            is = r% iseq(kk)
            js = r% iseq(jj)

            IF ( kk > 1 ) THEN
            IF ( r% ibsp(kk-1) == 0 .or. kk == ii+1 ) THEN

              CALL EDANGLE (r%iseq,kk,jj,kk-1,n,e5)

              IF ( kk > 2 ) THEN

                ip = kk - 2
                jp = r% ibsp(ip)

                IF ( jp /= 0 ) THEN
                  CALL EDANGLE (r%iseq,ip,jp,kk-1,n,ed)
                  e5 = MIN(e5,ed)
                ENDIF

              ENDIF

            ENDIF
            ENDIF

            IF ( jj < n ) THEN
            IF ( r% ibsp(jj+1) == 0 ) THEN

              CALL EDANGLE (r%iseq,kk,jj,jj+1,n,e3)

              IF ( jj < n-1 ) THEN

                ip = jj + 2
                jp = r% ibsp(ip)

                IF ( jp /= 0 ) THEN
                  CALL EDANGLE (r%iseq,ip,jp,jj+1,n,ed)
                  e3 = MIN(e3,ed)
                ENDIF

              ENDIF

            ENDIF
            ENDIF

          ENDIF

          e1 = e1 + e3 + e5 + eaup(is,js)

        ENDIF


        IF ( mh == 1 .and. kloop == 1 ) THEN

          IF ( kk == ii-1 .or. kk == ii+1 ) THEN
            CALL EHAIR (r%iseq,kk,jj,n,e2)
          ENDIF

          IF ( kk == jj-1 .or. kk == jj+1 ) THEN
            CALL EHAIR (r%iseq,ii,kk,n,e2)
          ENDIF

        ELSEIF ( mh == 2 .and. kloop == 1 ) THEN

          IF ( ms == 0 ) THEN

            IF ( kk == ii+1 ) THEN
              CALL ESTACK (r%iseq,kk,jj,kk+1,jj-1,n,e2)
            ELSE
              CALL ESTACK (r%iseq,ii,kk,ii+1,kk-1,n,e2)
            ENDIF

          ELSE

            ip = ii + 1
            DO WHILE ( r% ibsp(ip) == 0 )
            ip = ip + 1
            ENDDO

            jp = r% ibsp(ip)

            IF ( kk == ii+1 .or. kk == ii-1 ) THEN
              CALL EBULGE (r%iseq,kk,jj,ip,jp,n,e2)
            ELSE
              CALL EBULGE (r%iseq,ii,kk,ip,jp,n,e2)
            ENDIF

          ENDIF

        ELSE

          e2 = 0.0e0
          e5 = 0.0e0
          e3 = 0.0e0

          IF ( kloop == 1 ) THEN

            IF ( ms <= 6 ) THEN
              e2 = em + es * REAL(ms) + eh * REAL(mh)
            ELSE
              x = REAL(ms) / 6.0e0
              e2 = em + es * 6.0e0 + eh * REAL(mh)
              e2 = e2 + c * LOG(x)
            ENDIF

          ENDIF

          IF ( kk == jj+1 .or. kk == jj-1 ) THEN

            is = r% iseq(ii)
            js = r% iseq(kk)

            IF ( kk > 1 ) THEN
            IF ( r% ibsp(kk-1) == 0 .or. kk == jj+1 ) THEN

              CALL EDANGLE (r%iseq,ii,kk,kk-1,n,e5)

              IF ( kk > 2 ) THEN

                ip = kk - 2
                jp = r% ibsp(ip)

                IF ( jp /= 0 ) THEN
                  CALL EDANGLE (r%iseq,ip,jp,kk-1,n,ed)
                  e5 = MIN(e5,ed)
                ENDIF

              ENDIF

            ENDIF
            ENDIF

            IF ( ii < n ) THEN
            IF ( r% ibsp(ii+1) == 0 ) THEN

              CALL EDANGLE (r%iseq,ii,kk,ii+1,n,e3)

              IF ( ii < n-1 ) THEN

                ip = ii + 2
                jp = r% ibsp(ip)

                IF ( jp /= 0 ) THEN
                  CALL EDANGLE (r%iseq,ip,jp,ii+1,n,ed)
                  e3 = MIN(e3,ed)
                ENDIF

              ENDIF

            ENDIF
            ENDIF

          ENDIF

          IF ( kk == ii-1 .or. kk == ii+1 ) THEN

            is = r% iseq(kk)
            js = r% iseq(jj)

            IF ( jj > 1 ) THEN
            IF ( r% ibsp(jj-1) == 0 ) THEN

              CALL EDANGLE (r%iseq,kk,jj,jj-1,n,e5)

              IF ( jj > 2 ) THEN

                ip = jj - 2
                jp = r% ibsp(ip)

                IF ( jp /= 0 ) THEN
                  CALL EDANGLE (r%iseq,ip,jp,jj-1,n,ed)
                  e5 = MIN(e5,ed)
                ENDIF

              ENDIF

            ENDIF
            ENDIF

            IF ( kk < n ) THEN
            IF ( r% ibsp(kk+1) == 0 .or. kk == ii-1 ) THEN

              CALL EDANGLE (r%iseq,kk,jj,kk+1,n,e3)

              IF ( kk < n-1 ) THEN

                ip = kk + 2
                jp = r% ibsp(ip)

                IF ( jp /= 0 ) THEN
                  CALL EDANGLE (r%iseq,ip,jp,kk+1,n,ed)
                  e3 = MIN(e3,ed)
                ENDIF

              ENDIF

            ENDIF
            ENDIF

          ENDIF

          e2 = e2 + e3 + e5 + eaup(is,js)

        ENDIF

        ef = e1 + e2 + e4

        dg = DBLE(ef) - DBLE(ei)

        RETURN

      END SUBROUTINE DELTAG_HD
