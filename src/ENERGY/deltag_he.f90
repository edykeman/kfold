! ==============================================================================
! Subroutine: DELTAG_HE (R,II,JJ,DG)
! 
! Purpose: Computes the difference in free energy of an RNA loop due
!          to extension of the helix below the ii-jj base pair using the
!          emperical MFOLD 3.0 energy function.
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
!                 helix below the ii-jj base-pair (ii-1 base-pairs with jj+1)
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

      SUBROUTINE DELTAG_HE (R,II,JJ,DG)

        USE RNAvar, ONLY : em,eh,es,eaup,beta

        USE Class_RNAFold

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE(RNA_STRUC), INTENT(IN) :: r
        INTEGER, INTENT(IN) :: ii,jj

        DOUBLE PRECISION, INTENT(OUT) :: dg

        !=== VARIABLES ===!

        INTEGER :: i,j,k,n,ip,jp,is,js
        INTEGER :: nh,ns,mh,ms,iloop,indx

        REAL :: e1,e3,e5,ei,ef
        REAL :: x,c,ed


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
        ms = ns - 2


        !=== INITIAL ENERGY ===!

        IF ( nh == 1 .and. iloop == 1 ) THEN

          CALL EHAIR (r%iseq,i,j,n,ei)

        ELSEIF ( nh == 2 .and. iloop == 1 ) THEN

          ip = i + 1
          DO WHILE ( r% ibsp(ip) == 0 )
          ip = ip + 1
          ENDDO

          jp = r% ibsp(ip)

          CALL EBULGE (r%iseq,i,j,ip,jp,n,ei)

        ELSE

          ei = 0.0e0
          e5 = 0.0e0
          e3 = 0.0e0

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

          IF ( ii > 1 ) THEN
          IF ( r% ibsp(ii-1) == 0 ) THEN

            CALL EDANGLE (r%iseq,ii,jj,ii-1,n,e5)

            IF ( ii > 2 ) THEN

              ip = ii - 2
              jp = r% ibsp(ip)

              IF ( jp /= 0 ) THEN
                CALL EDANGLE (r%iseq,ip,jp,ii-1,n,ed)
                e5 = MIN(e5,ed)
              ENDIF

            ENDIF

            IF ( ii > 3 ) THEN
            IF ( r% ibsp(ii-2) == 0 ) THEN

              ip = ii - 3
              jp = r% ibsp(ip)

              IF ( jp /= 0 ) THEN
                CALL EDANGLE (r%iseq,ip,jp,ii-2,n,ed)
                ei = ei + ed
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
                e3 = MIN(e3,ed)
              ENDIF

            ENDIF

            IF ( jj < n-2 ) THEN
            IF ( r% ibsp(jj+2) == 0 ) THEN

              ip = jj + 3
              jp = r% ibsp(ip)

              IF ( jp /= 0 ) THEN
                CALL EDANGLE (r%iseq,ip,jp,jj+2,n,ed)
                ei = ei + ed
              ENDIF

            ENDIF
            ENDIF

          ENDIF
          ENDIF

          ei = ei + e3 + e5 + eaup(is,js)

        ENDIF


        !=== FINAL ENERGY ===!

        CALL ESTACK (r%iseq,ii-1,jj+1,ii,jj,n,ef)

        IF ( mh == 1 .and. iloop == 1 ) THEN

          CALL EHAIR (r%iseq,i+1,j-1,n,e1)

        ELSEIF ( mh == 2 .and. iloop == 1 ) THEN

          IF ( ms == 0 ) THEN

            CALL ESTACK (r%iseq,ii-2,jj+2,ii-1,jj+1,n,e1)

          ELSE

            ip = i + 1
            DO WHILE ( r% ibsp(ip) == 0 )
            ip = ip + 1
            ENDDO

            jp = r% ibsp(ip)

            IF ( j == ii ) THEN
              CALL EBULGE (r%iseq,i+1,j-1,ip,jp,n,e1)
            ELSE
              CALL EBULGE (r%iseq,i,j,ip-1,jp+1,n,e1)
            ENDIF

          ENDIF

        ELSE

          e1 = 0.0e0
          e5 = 0.0e0
          e3 = 0.0e0

          is = r% iseq(ii-1)
          js = r% iseq(jj+1)

          IF ( iloop == 1 ) THEN

            IF ( ms <= 6 ) THEN
              e1 = em + es * REAL(ms) + eh * REAL(mh)
            ELSE
              x = REAL(ms) / 6.0e0
              e1 = em + es * 6.0e0 + eh * REAL(mh)
              e1 = e1 + c * LOG(x)
            ENDIF

          ENDIF

          IF ( ii > 2 ) THEN
          IF ( r% ibsp(ii-2) == 0 ) THEN

            CALL EDANGLE (r%iseq,ii-1,jj+1,ii-2,n,e5)

            IF ( ii > 3 ) THEN

              ip = ii - 3
              jp = r% ibsp(ip)

              IF ( jp /= 0 ) THEN
                CALL EDANGLE (r%iseq,ip,jp,ii-2,n,ed)
                e5 = MIN(e5,ed)
              ENDIF

            ENDIF

          ENDIF
          ENDIF

          IF ( jj < n-1 ) THEN
          IF ( r% ibsp(jj+2) == 0 ) THEN

            CALL EDANGLE (r%iseq,ii-1,jj+1,jj+2,n,e3)

            IF ( jj < n-2 ) THEN

              ip = jj + 3
              jp = r% ibsp(ip)

              IF ( jp /= 0 ) THEN
                CALL EDANGLE (r%iseq,ip,jp,jj+2,n,ed)
                e3 = MIN(e3,ed)
              ENDIF

            ENDIF

          ENDIF
          ENDIF

          e1 = e1 + e3 + e5 + eaup(is,js)

        ENDIF

        ef = ef + e1

        dg = DBLE(ef) - DBLE(ei)

        RETURN

      END SUBROUTINE DELTAG_HE
