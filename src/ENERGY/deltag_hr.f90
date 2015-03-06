! ==============================================================================
! Subroutine: DELTAG_HR (R,II,JJ,DG)
! 
! Purpose: Computes the difference in free energy of an RNA loop due
!          to a deletion of the ii-jj base pair using the emperical
!          MFOLD 3.0 energy function.
!
! Method: Uses the MFOLD 3.0 energy function for RNA @ T=37.
!
! Arguments:
!
!             R - Class containing information about the RNA fold and
!                 the RNA sequence.
!            II - Nucleotide position of the 5' most nucleotide.
!            JJ - Nucleotide position of the 3' most nucleotide.
!            DG - (OUTPUT) The energy difference from deletion of the
!                 base-pair ii-jj
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
! Subroutines - ESTACK, EHAIR, EBULGE, EDANGLE
!
! Author(s): Eric Dykeman
!            Copyright (c) 2015 (Please Refer to LICENCE)
!
! ==============================================================================

      SUBROUTINE DELTAG_HR (R,II,JJ,DG)

        USE RNAvar, ONLY : em,eh,es,eaup,beta

        USE Class_RNAFold

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE(RNA_STRUC), INTENT(IN) :: r
        INTEGER, INTENT(IN) :: ii,jj

        DOUBLE PRECISION, INTENT(OUT) :: dg

        !=== VARIABLES ===!

        INTEGER :: i,j,k,n,ip,jp,is,js
        INTEGER :: nh,ns,mh,ms,lh,ls
        INTEGER :: iloop,indx

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

        IF ( i > j ) THEN
          iloop = 0
          i = 1
          j = n
        ENDIF

        IF ( r% ibsp(ii+1) /= jj-1 ) THEN

          k  = r% link(ii)
          mh = r% nhlx(k)
          ms = r% nsgl(k)

        ELSE

          mh = 2
          ms = 0

        ENDIF

        !=== Final Loop Size ===!

        lh = nh + mh - 2
        ls = ns + ms + 2


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

          IF ( mh > 2 ) THEN

            is = r% iseq(i)
            js = r% iseq(j)

            IF ( r% ibsp(i+1) == 0 ) THEN
              CALL EDANGLE (r%iseq,i,j,i+1,n,ed)
              e4 = e4 + ed
            ENDIF

            IF ( r% ibsp(j-1) == 0 ) THEN
              CALL EDANGLE (r%iseq,i,j,j-1,n,ed)
              e4 = e4 + ed
            ENDIF

            e4 = e4 + eaup(is,js)

          ENDIF

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
                e5 = MIN(e5,ed)
                e4 = e4 + ed
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
                e3 = MIN(e3,ed)
                e4 = e4 + ed
              ENDIF

            ENDIF

          ENDIF
          ENDIF

          e1 = e1 + e3 + e5 + eaup(is,js)

        ENDIF


        IF ( mh == 1 ) THEN

          CALL EHAIR (r%iseq,ii,jj,n,e2)

        ELSEIF ( mh == 2 ) THEN

          IF ( ms == 0 ) THEN

            CALL ESTACK (r%iseq,ii,jj,ii+1,jj-1,n,e2)

          ELSE

            ip = ii + 1
            DO WHILE ( r% ibsp(ip) == 0 )
            ip = ip + 1
            ENDDO

            jp = r% ibsp(ip)

            CALL EBULGE (r%iseq,ii,jj,ip,jp,n,e2)

            IF ( nh > 2 .or. iloop == 0 ) THEN

              is = r% iseq(ip)
              js = r% iseq(jp)

              IF ( r% ibsp(ip-1) == 0 ) THEN
                CALL EDANGLE (r%iseq,ip,jp,ip-1,n,ed)
                e4 = e4 + ed
              ENDIF

              IF ( r% ibsp(jp+1) == 0 ) THEN
                CALL EDANGLE (r%iseq,ip,jp,jp+1,n,ed)
                e4 = e4 + ed
              ENDIF

              e4 = e4 + eaup(is,js)

            ENDIF

          ENDIF

        ELSE

          e2 = 0.0e0
          e5 = 0.0e0
          e3 = 0.0e0

          is = r% iseq(ii)
          js = r% iseq(jj)

          IF ( ms <= 6 ) THEN
            e2 = em + es * REAL(ms) + eh * REAL(mh)
          ELSE
            x = REAL(ms) / 6.0e0
            e2 = em + es * 6.0e0 + eh * REAL(mh)
            e2 = e2 + c * LOG(x)
          ENDIF

          IF ( r% ibsp(jj-1) == 0 ) THEN

            CALL EDANGLE (r%iseq,ii,jj,jj-1,n,e5)

            ip = jj - 2
            jp = r% ibsp(ip)

            IF ( jp /= 0 ) THEN
              CALL EDANGLE (r%iseq,ip,jp,jj-1,n,ed)
              e5 = MIN(e5,ed)
              e4 = e4 + ed
            ENDIF

          ENDIF

          IF ( r% ibsp(ii+1) == 0 ) THEN

            CALL EDANGLE (r%iseq,ii,jj,ii+1,n,e3)

            ip = ii + 2
            jp = r% ibsp(ip)

            IF ( jp /= 0 ) THEN
              CALL EDANGLE (r%iseq,ip,jp,ii+1,n,ed)
              e3 = MIN(e3,ed)
              e4 = e4 + ed
            ENDIF

          ENDIF

          e2 = e2 + e3 + e5 + eaup(is,js)

        ENDIF

        ei = e1 + e2


        !=== FINAL ENERGY ===!

        IF ( lh == 1 .and. iloop == 1 ) THEN

          IF ( j == ii ) THEN
            CALL EHAIR (r%iseq,i-1,j+1,n,ef)
          ELSE
            CALL EHAIR (r%iseq,i,j,n,ef)
          ENDIF

        ELSEIF ( lh == 2 .and. iloop == 1 ) THEN

          ip = ii + 1
          IF ( mh == 1 ) ip = i + 1
          IF ( j == ii ) ip = i + 1

          DO WHILE ( r% ibsp(ip) == 0 )
          ip = ip + 1
          ENDDO

          jp = r% ibsp(ip)

          IF ( nh > 2 ) THEN

            e3 = 0.0e0
            e5 = 0.0e0

            IF ( ip == ii ) THEN

              ip = jj + 1
              DO WHILE ( r% ibsp(ip) == 0 )
              ip = ip + 1
              ENDDO

              jp = r% ibsp(ip)

            ENDIF

            is = r% iseq(i)
            js = r% iseq(j)

            IF ( i+2 /= ii ) THEN
            IF ( r% ibsp(i+1) == 0 ) THEN

              CALL EDANGLE (r%iseq,i,j,i+1,n,e3)

            ENDIF
            ENDIF

            IF ( j-2 /= jj ) THEN
            IF ( r% ibsp(j-1) == 0 ) THEN

              CALL EDANGLE (r%iseq,i,j,j-1,n,e5)

            ENDIF
            ENDIF

            e1 = e1 + eaup(is,js)

            is = r% iseq(ip)
            js = r% iseq(jp)

            IF ( ip-2 /= jj ) THEN
            IF ( r% ibsp(ip-1) == 0 ) THEN

              CALL EDANGLE (r%iseq,ip,jp,ip-1,n,ed)

              IF ( ip-2 == i ) THEN
                e3 = MIN(e3,ed)
              ELSE
                e3 = e3 + ed
              ENDIF

            ENDIF
            ENDIF

            IF ( jp+2 /= ii ) THEN
            IF ( r% ibsp(jp+1) == 0 ) THEN

              CALL EDANGLE (r%iseq,ip,jp,jp+1,n,ed)

              IF ( jp+2 == j ) THEN
                e5 = MIN(e5,ed)
              ELSE
                e5 = e5 + ed
              ENDIF

            ENDIF
            ENDIF

            e1 = e1 + eaup(is,js)
            e1 = e1 + e3 + e5

            ei = e1 + e2

          ENDIF

          IF ( j == ii ) THEN
            CALL EBULGE (r%iseq,i-1,j+1,ip,jp,n,ef)
          ELSE
            CALL EBULGE (r%iseq,i,j,ip,jp,n,ef)
          ENDIF

        ELSE

          ef = 0.0e0
          e5 = 0.0e0
          e3 = 0.0e0

          IF ( iloop == 1 ) THEN

            IF ( ls <= 6 ) THEN
              ef = em + es * REAL(ls) + eh * REAL(lh)
            ELSE
              x = REAL(ls) / 6.0e0
              ef = em + es * 6.0e0 + eh * REAL(lh)
              ef = ef + c * LOG(x)
            ENDIF

          ENDIF

          ip = ii + 1
          jp = r% ibsp(ip)

          IF ( jp /= 0 ) THEN
            CALL EDANGLE (r%iseq,ip,jp,ii,n,e5)
          ENDIF

          IF ( ii > 1 ) THEN

            ip = ii - 1
            jp = r% ibsp(ip)

            IF ( jp /= 0 ) THEN
              CALL EDANGLE (r%iseq,ip,jp,ii,n,ed)
              e5 = MIN(e5,ed)
            ENDIF

          ENDIF

          ip = jj - 1
          jp = r% ibsp(ip)

          IF ( jp /= 0 ) THEN
            CALL EDANGLE (r%iseq,ip,jp,jj,n,e3)
          ENDIF

          IF ( jj < n ) THEN

            ip = jj + 1
            jp = r% ibsp(ip)

            IF ( jp /= 0 ) THEN
              CALL EDANGLE (r%iseq,ip,jp,jj,n,ed)
              e3 = MIN(e3,ed)
            ENDIF

          ENDIF

          IF ( mh == 2 .and. ms == 0 ) THEN

            is = r% iseq(ii+1)
            js = r% iseq(jj-1)

            ef = ef + eaup(is,js)

          ENDIF

          ef = ef + e3 + e4 + e5

        ENDIF

        dg = DBLE(ef) - DBLE(ei)

        RETURN

      END SUBROUTINE DELTAG_HR
