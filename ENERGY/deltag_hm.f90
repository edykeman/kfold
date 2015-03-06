! ==============================================================================
! Subroutine: DELTAG_HM (R,II,JJ,DG)
! 
! Purpose: Computes the difference in free energy of an RNA loop due
!          to the swapping of base-pairs between two helices using the
!          emperical MFOLD 3.0 energy function. The swapping of base-pairs
!          occurs as a result of the extension of the helix below the ii-jj
!          base-pair.
!
! Method: Uses the MFOLD 3.0 energy function for RNA @ T=37.
!
! Arguments:
!
!             R - Class containing information about the RNA fold and
!                 the RNA sequence.
!            II - Nucleotide position of the 5' most nucleotide.
!            JJ - Nucleotide position of the 3' most nucleotide.
!            DG - (OUTPUT) The energy difference from extension of the
!                 helix below base-pair ii-jj.
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
! Subroutines - ESTACK, EDANGLE
!
! Author(s): Eric Dykeman
!            Copyright (c) 2015 (Please Refer to LICENCE)
!
! ==============================================================================

      SUBROUTINE DELTAG_HM (R,II,JJ,DG)

        USE RNAvar, ONLY : em,eh,es,eaup,beta

        USE Class_RNAFold

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE(RNA_STRUC), INTENT(IN) :: r
        INTEGER, INTENT(IN) :: ii,jj

        DOUBLE PRECISION, INTENT(OUT) :: dg

        !=== VARIABLES ===!

        INTEGER :: i,j,k,n,ip,jp,kp,lp
        INTEGER :: is,js,nh,ns,mh,ms
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

        !=== Final Loop Size ===!

        mh = nh
        ms = ns

        ip = r% ibsp(ii-1)
        jp = r% ibsp(jj+1)

        IF ( ip /= 0 .and. jp /= 0 ) ms = ms + 2


        e4 = 0.0e0

        !=== INITIAL ENERGY ===!

        ei = 0.0e0
        e1 = 0.0e0
        e2 = 0.0e0
        e3 = 0.0e0
        e5 = 0.0e0

        is = r% iseq(ii)
        js = r% iseq(jj)

        IF ( iloop == 1 ) THEN

          IF ( ns <= 6 ) THEN
            ei = em + es * REAL(ns) + eh * REAL(nh)
          ELSE
            x = REAL(ns) / 6.0e0
            ei = em + es * 6.0e0 + eh * REAL(nh)
            ei = ei + c * LOG(x)
          ENDIF

        ENDIF

        ei = ei + eaup(is,js)


        !=== 5' Side ===!

        IF ( r% ibsp(ii-1) /= 0 ) THEN

          jp = ii - 1
          ip = r% ibsp(jp)

          is = r% iseq(ip)
          js = r% iseq(jp)

          CALL ESTACK (r%iseq,ip,jp,ip+1,jp-1,n,e1)

          IF ( ip > 1 ) THEN
          IF ( r% ibsp(ip-1) == 0 ) THEN

            CALL EDANGLE (r%iseq,ip,jp,ip-1,n,e5)

            IF ( ip > 2 ) THEN

              kp = ip - 2
              lp = r% ibsp(kp)

              IF ( lp /= 0 ) THEN
                CALL EDANGLE (r%iseq,kp,lp,ip-1,n,ed)
                IF ( lp /= jj+1 ) e4 = e4 + ed
                e5 = MIN(e5,ed)
              ENDIF

            ENDIF

          ENDIF
          ENDIF

          e1 = e1 + e5 + eaup(is,js)

          is = r% iseq(ip+1)
          js = r% iseq(jp-1)

          CALL EDANGLE (r%iseq,ip+1,jp-1,ip,n,e5)

          IF ( ip > 1 ) THEN
          IF ( r% ibsp(ip-1) /= 0 ) THEN

            kp = ip - 1
            lp = r% ibsp(kp)

            IF ( lp /= jj+1 ) THEN
              CALL EDANGLE (r%iseq,kp,lp,ip,n,ed)
              e5 = MIN(e5,ed)
            ENDIF

          ENDIF
          ENDIF

          e4 = e4 + e5 + eaup(is,js)

        ELSE

          CALL EDANGLE (r%iseq,ii,jj,ii-1,n,e5)

          IF ( ii > 2 ) THEN

            ip = ii - 2
            jp = r% ibsp(ip)

            IF ( jp /= 0 ) THEN

              CALL EDANGLE (r%iseq,ip,jp,ii-1,n,ed)
              e5 = MIN(e5,ed)

            ELSEIF ( ii > 3 ) THEN

              kp = ii - 3
              lp = r% ibsp(kp)

              IF ( lp /= 0 ) THEN
                CALL EDANGLE (r%iseq,kp,lp,ii-2,n,ed)
                e5 = e5 + ed
              ENDIF

            ENDIF

          ENDIF

          e1 = e5

        ENDIF


        !=== 3' Side ===!

        IF ( r% ibsp(jj+1) /= 0 ) THEN

          ip = jj + 1
          jp = r% ibsp(ip)

          is = r% iseq(ip)
          js = r% iseq(jp)

          CALL ESTACK (r%iseq,ip,jp,ip+1,jp-1,n,e2)

          IF ( jp < n ) THEN
          IF ( r% ibsp(jp+1) == 0 ) THEN

            CALL EDANGLE (r%iseq,ip,jp,jp+1,n,e3)

            IF ( jp < n-1 ) THEN

              kp = jp + 2
              lp = r% ibsp(kp)

              IF ( lp /= 0 ) THEN

                CALL EDANGLE (r%iseq,kp,lp,jp+1,n,ed)

                IF ( lp /= ii-1 ) THEN
                  e4 = e4 + ed
                  e3 = MIN(e3,ed)
                ELSE
                  e3 = 0.0e0
                ENDIF

              ENDIF

            ENDIF

          ENDIF
          ENDIF

          e2 = e2 + e3 + eaup(is,js)

          is = r% iseq(ip+1)
          js = r% iseq(jp-1)

          CALL EDANGLE (r%iseq,ip+1,jp-1,jp,n,e3)

          IF ( jp < n ) THEN
          IF ( r% ibsp(jp+1) /= 0 ) THEN

            kp = jp + 1
            lp = r% ibsp(kp)

            IF ( lp /= ii-1 ) THEN
              CALL EDANGLE (r%iseq,kp,lp,jp,n,ed)
              e3 = MIN(e3,ed)
            ENDIF

          ENDIF
          ENDIF

          e4 = e4 + e3 + eaup(is,js)

        ELSE

          CALL EDANGLE (r%iseq,ii,jj,jj+1,n,e3)

          IF ( jj < n-1 ) THEN

            ip = jj + 2
            jp = r% ibsp(ip)

            IF ( jp /= 0 ) THEN

              CALL EDANGLE (r%iseq,ip,jp,jj+1,n,ed)
              e3 = MIN(e3,ed)

            ELSEIF ( jj < n-2 ) THEN

              kp = jj + 3
              lp = r% ibsp(kp)

              IF ( lp /= 0 ) THEN
                CALL EDANGLE (r%iseq,kp,lp,jj+2,n,ed)
                e3 = e3 + ed
              ENDIF

            ENDIF

          ENDIF

          e2 = e3

        ENDIF

        ei = ei + e1 + e2


        !=== FINAL ENERGY ===!

        ef = 0.0e0
        e5 = 0.0e0
        e3 = 0.0e0

        ip = ii - 1
        jp = jj + 1

        is = r% iseq(ip)
        js = r% iseq(jp)

        IF ( iloop == 1 ) THEN

          IF ( ms <= 6 ) THEN
            ef = em + es * REAL(ms) + eh * REAL(mh)
          ELSE
            x = REAL(ms) / 6.0e0
            ef = em + es * 6.0e0 + eh * REAL(mh)
            ef = ef + c * LOG(x)
          ENDIF

        ENDIF

        CALL ESTACK (r%iseq,ip,jp,ii,jj,n,ed)

        ef = ef + ed + eaup(is,js)

        IF ( ip > 1 ) THEN
        IF ( r% ibsp(ip-1) == 0 ) THEN

          CALL EDANGLE (r%iseq,ip,jp,ip-1,n,e5)

          IF ( ip > 2 ) THEN

            kp = ip - 2
            lp = r% ibsp(kp)

            IF ( lp /= 0 ) THEN
              CALL EDANGLE (r%iseq,kp,lp,ip-1,n,ed)
              e5 = MIN(e5,ed)
            ENDIF

          ENDIF

        ENDIF
        ENDIF

        IF ( jp < n ) THEN
        IF ( r% ibsp(jp+1) == 0 ) THEN

          CALL EDANGLE (r%iseq,ip,jp,jp+1,n,e3)

          IF ( jp < n-1 ) THEN

            kp = jp + 2
            lp = r% ibsp(kp)

            IF ( lp /= 0 ) THEN
              CALL EDANGLE (r%iseq,kp,lp,jp+1,n,ed)
              e3 = MIN(e3,ed)
            ENDIF

          ENDIF

        ENDIF
        ENDIF

        ef = ef + e5 + e3 + e4

        dg = DBLE(ef) - DBLE(ei)

        RETURN

      END SUBROUTINE DELTAG_HM
