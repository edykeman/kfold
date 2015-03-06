! ==============================================================================
! Subroutine: LOOP_REAC (R,INDX)
! 
! Purpose: Recomputes the possible reactions within a loop element and
!          their corrisponding rates. 
!
! Method:  There are 6 possible reactions around the loop:
!
!          (1) Nucleation                   A
!                                          A A
!              A A A A A A U A A  -->  A A A-U A A
!
!          (2) Helix Extension
!                                        x-x
!              A A A x-x U A A  -->  A A A-U A A
!
!          (3) Helix Retraction
!                    x-x
!                A A A-U A A   -->  A A A x-x U A A
!
!          (4) Helix Morphing               x-x
!                  x-x  x-x                 U-A
!                A A-U  U-A A  -->  A A x-x U-A
!
!          (5) Defect Diffusion         x-x
!                    x-x                  U
!                A A A-U U A A -->  A A A-U A A
!
!          (6) Helix Open
!                         U-A       U-A
!                         x-x      x   x
!                         A-U  -->  A-U
!
! Arguments:
!
!             R - Class structure containing information on the
!                 RNA secondary structure and possible reactions.
!          INDX - The indx number of the loop element that reactions
!                 will be calculated for.
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

      SUBROUTINE LOOP_REAC (R,INDX)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE(RNA_STRUC), INTENT(INOUT) :: r

        INTEGER, INTENT(IN) :: indx

        !=== VARIABLES ===!

        INTEGER :: i,j,k,n,ip,jp,kp,is,js
        INTEGER :: ks,ke,l,lmx,icnt,iloop
        INTEGER :: nt,nh,ns,mt,mh,ms,icase

        DOUBLE PRECISION :: x,dg,atot,rate


        i = r% loop(indx)
        j = r% ibsp(i)
        n = r% n

        IF ( i == n ) j = 1

        nh = r% nhlx(indx)
        ns = r% nsgl(indx)

        nt = ns + 2 * nh

        !=== Internal Loop iloop = 1 ===!
        !=== External Loop iloop = 0 ===!

        IF ( i < j ) iloop = 1
        IF ( i > j ) iloop = 0

        IF ( iloop == 1 ) THEN
          ks = i + 1
          ke = j
        ELSE
          ks = j
          ke = i
        ENDIF


        !=== COMPUTE REACTIONS ===!

        atot = 0.0d0

        k = ks
        icnt = 0

        DO WHILE ( k <= ke )

          !=== Nucleation Events ===!

          IF ( r% ibsp(k) == 0 ) THEN

            r% wrk1(k) = 0.0d0
            r% wrk2(k) = 0.0d0

            x = 0.0d0

            l = 2
            lmx = nt / 2 + 1

            IF ( MOD(nt,2) == 0 ) THEN
            IF ( icnt+1 > lmx-1 ) lmx = lmx - 1
            ENDIF

            IF ( iloop == 0 ) lmx = nt - icnt

            kp = k + 1
            is = r% iseq(k)

            DO WHILE ( l <= lmx )

              IF ( r% ibsp(kp) == 0 ) THEN

                js = r% iseq(kp)

                IF ( l > 4 .and. iwc(is,js) == 1 ) THEN
                  x = x + pnuc(l)
                ENDIF

              ELSE

                l = l + 1
                kp = r% ibsp(kp)

              ENDIF

              l = l + 1
              kp = kp + 1

            ENDDO
 
            r% wrk1(k) = x
            atot = atot + x

          ENDIF


          !=== Helix Events ===!

          IF ( r% ibsp(k) > 0 ) THEN

            ip = k
            jp = r% ibsp(k)

            r% wrk1(jp) = 0.0d0
            r% wrk2(ip) = 0.0d0

            IF ( r% link(ip) == 0 ) THEN
            r% wrk1(ip) = 0.0d0
            r% wrk2(jp) = 0.0d0
            ENDIF

            icase = 0

            !=== Helix Extension ===!

            IF ( ip > 1 .and. jp < n ) THEN
            IF ( nh > 1  .or. ns > 4 ) THEN
            IF ( r%ibsp(ip-1) == 0 .and. r%ibsp(jp+1) == 0 ) THEN

              is = r% iseq(ip-1)
              js = r% iseq(jp+1)

              IF ( iwc(is,js) == 1 ) THEN

                icase = 1

                IF ( iloop == 1 ) THEN
                IF ( nh == 2 .and. ns == 2 ) THEN
                  icase = 2
                  IF ( k == ke ) icase = 0
                ELSEIF ( k == ke ) THEN
                  icase = 3
                ENDIF
                ENDIF

              ENDIF

            ENDIF
            ENDIF
            ENDIF

            IF ( icase > 0 ) THEN

              CALL DELTAG_HE (r,ip,jp,dg)

              dg = dg / 2.0d0

              x = beta * dg
              x = DEXP(-x) * rateh

              r% wrk2(ip) = x
              atot = atot + x

            ENDIF

            icase = 0

            !=== Helix Retraction ===!

            IF ( ip /= n .and. jp /= 1 ) THEN

              IF ( r% ibsp(ip+1) == jp-1 ) THEN

                icase = 1

                IF ( iloop == 1 ) THEN
                IF ( k == ke ) icase = 2
                ENDIF

                rate = rateh

              ELSEIF ( iloop == 0 .or. k /= ke ) THEN

                icase = 3

                l  = r% link(ip)
                mh = r% nhlx(l)
                ms = r% nsgl(l)

                mt = ms + 2 * mh

                IF ( iloop == 1 ) THEN
                  l = MIN(nt,mt)
                ELSE
                  l = mt
                ENDIF

                rate = pnuc(l)

              ENDIF

            ENDIF

            IF ( icase > 0 ) THEN

              CALL DELTAG_HR (r,ip,jp,dg)

              IF ( icase /= 3 ) THEN
                dg = dg / 2.0d0
              ENDIF

              x = beta * dg
              x = DEXP(-x) * rate

              IF ( icase == 2 ) THEN
                r% wrk1(ip) = x
              ELSE
                r% wrk1(jp) = x
              ENDIF

              atot = atot + x

            ENDIF

            icase = 0

            !=== Helix Morphing ===!

            IF ( iloop == 0 .or. nh > 2 ) THEN
            IF ( ip > 1 .and. jp < n ) THEN

              is = r% iseq(ip-1)
              js = r% iseq(jp+1)

              IF ( iwc(is,js) == 1 ) icase = 1

              is = r% ibsp(ip-1)
              js = r% ibsp(jp+1)

              IF ( is /= 0 ) THEN
              IF ( r% link(is) /= 0 ) icase = 0
              ENDIF

              IF ( js /= 0 ) THEN
              IF ( r% link(jp+1) /= 0 ) icase = 0
              ENDIF

              IF ( is == 0 .and. js == 0 ) icase = 0

            ENDIF
            ENDIF

            IF ( icase > 0 ) THEN

              CALL DELTAG_HM (r,ip,jp,dg)

              dg = dg / 2.0d0

              x = beta * dg
              x = DEXP(-x) * ratem

              atot = atot + x

            ENDIF

            !=== Defect Diffusion ===!

            !=== PUSH ===!

            icase = 0

            IF ( r% link(ip) == 0 ) THEN

              icase = 2

              IF ( iloop == 1 ) THEN
              IF ( nh == 2 .and. ns == 1 ) icase = 3
              IF ( nh == 1 .and. ns == 3 ) icase = 0
              ENDIF

            ELSEIF ( iloop == 0 .or. k /= ke ) THEN

              icase = 1

              IF ( iloop == 1 ) THEN
              IF ( nh == 2 .and. ns == 1 ) icase = 4
              ENDIF

            ENDIF

            IF ( icase > 0 ) THEN

              !=== Push 5' End ===!

              kp = ip - 1

              IF ( kp >= 1 ) THEN
              IF ( r% ibsp(kp) == 0 ) THEN

                is = r% iseq(kp)
                js = r% iseq(jp)

                IF ( iwc(is,js) == 1 ) THEN

                  CALL DELTAG_HD (r,ip,jp,kp,dg)

                  dg = dg / 2.0d0

                  x = beta * dg
                  x = DEXP(-x) * rated

                  atot = atot + x

                ENDIF

              ENDIF
              ENDIF

              !=== Push 3' End ===!

              kp = jp + 1

              IF ( kp <= n ) THEN
              IF ( r% ibsp(kp) == 0 ) THEN

                is = r% iseq(ip)
                js = r% iseq(kp)

                IF ( iwc(is,js) == 1 ) THEN

                  CALL DELTAG_HD (r,ip,jp,kp,dg)

                  dg = dg / 2.0d0

                  x = beta * dg
                  x = DEXP(-x) * rated

                  atot = atot + x

                ENDIF

              ENDIF
              ENDIF

            ENDIF

            !=== PULL === !

            icase = 0

            IF ( r% link(ip) /= 0 ) THEN
            IF ( iloop == 0 .or. k /= ke ) THEN

               l = r% link(ip)
              mh = r% nhlx(l)
              ms = r% nsgl(l)

              icase = 1

              IF ( mh == 1 .and. ms == 3 ) icase = 0
              IF ( mh == 2 .and. ms == 1 ) icase = 4

            ENDIF
            ENDIF

            IF ( icase > 0 ) THEN

              !=== Pull 5' End ===!

              kp = ip + 1

              IF ( r% ibsp(kp) == 0 ) THEN

                is = r% iseq(kp)
                js = r% iseq(jp)

                IF ( iwc(is,js) == 1 ) THEN

                  CALL DELTAG_HD (r,ip,jp,kp,dg)

                  dg = dg / 2.0d0

                  x = beta * dg
                  x = DEXP(-x) * rated

                  atot = atot + x

                ENDIF

              ENDIF

              !=== Pull 3' End ===!

              kp = jp - 1

              IF ( r% ibsp(kp) == 0 ) THEN

                is = r% iseq(ip)
                js = r% iseq(kp)

                IF ( iwc(is,js) == 1 ) THEN

                  CALL DELTAG_HD (r,ip,jp,kp,dg)

                  dg = dg / 2.0d0

                  x = beta * dg
                  x = DEXP(-x) * rated

                  atot = atot + x

                ENDIF

              ENDIF

            ENDIF

            !=== Save Loop Reaction Rate ===!

            IF ( iloop == 1 ) r% wrk1(i) = atot

            !=== Open Internal Helix BP ===!

            IF ( r% link(ip) == 0 ) THEN
            IF ( iloop == 1 .and. k == ke ) THEN

              is = ip + 1
              js = jp - 1

              DO WHILE ( r% link(is) == 0 )

                CALL DELTAG_HI (r,is,js,dg)

                dg = dg / 2.0d0

                x = beta * dg
                x = DEXP(-x) * rateh

                r% wrk1(is) = x
                r% wrk1(js) = 0.0d0

                r% wrk2(is) = 0.0d0
                r% wrk2(js) = 0.0d0

                atot = atot + x

                is = is + 1
                js = js - 1

              ENDDO

            ENDIF
            ENDIF


            IF ( k /= ke ) THEN
              k = r% ibsp(k)
              icnt = icnt + 1
            ENDIF

          ENDIF

          k = k + 1
          icnt = icnt + 1

        ENDDO

        r% ptot(indx) = atot

        CALL LOOP_RESUM (r,indx)

        RETURN

      END SUBROUTINE LOOP_REAC
