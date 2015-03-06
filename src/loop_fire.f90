! ==============================================================================
! Subroutine: LOOP_FIRE (R,INDX,AMAX)
! 
! Purpose: Finds the reaction J in the loop such that S >= AMAX where
!          S is the partial sum of the reaction rates for reactions
!          {1,J} in the loop. AMAX is determined from the SSA protocol.
!          Once reaction J is found, this reaction is "fired" and the
!          reactions for neighboring loop elements updated.
!
! Method:
!
! Arguments:
!
!             R - Class structure containing information on the
!                 RNA secondary structure and possible reactions.
!          INDX - The indx number of the loop element that a reaction
!                 will be choosen.
!          AMAX - The reaction "threshold" determined from the SSA protocol.
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

      SUBROUTINE LOOP_FIRE (R,INDX,AMAX)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        TYPE(RNA_STRUC), INTENT(INOUT) :: r

        DOUBLE PRECISION, INTENT(IN) :: amax
        INTEGER, INTENT(IN) :: indx

        !=== VARIABLES ===!

        INTEGER :: i,j,k,n,nl,ip,jp,kp,is,js
        INTEGER :: ks,ke,l,lmx,icnt,iloop
        INTEGER :: nt,nh,ns,mt,mh,ms,icase
        INTEGER :: jndx,kndx,nsum

        DOUBLE PRECISION :: x,dg,atot


        i = r% loop(indx)
        j = r% ibsp(i)

        n = r% n
        nl= r% nl

        nsum = r% nsum

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


        !=== FIND REACTION TO FIRE ===!

        atot = 0.0d0
        jndx = 0
        kndx = 0

        k = ks
        icnt = 0

        DO WHILE ( k <= ke )

          !=== Nucleation Events ===!

          IF ( r% ibsp(k) == 0 ) THEN
          IF ( atot + r% wrk1(k) >= amax ) THEN

            l = 2
            lmx = nt / 2 + 1

            IF ( MOD(nt,2) == 0 ) THEN
            IF ( icnt+1 > lmx-1 ) lmx = lmx - 1
            ENDIF

            IF ( iloop ==  0 ) lmx = nt - icnt

            kp = k + 1
            is = r% iseq(k)

            DO WHILE ( l <= lmx )

              IF ( r% ibsp(kp) == 0 ) THEN

                js = r% iseq(kp)

                IF ( l > 4 .and. iwc(is,js) == 1 ) THEN

                  atot = atot + pnuc(l)

                  IF ( atot >= amax ) THEN

                    nl = nl + 1

                    r% ibsp(k) = kp
                    r% ibsp(kp)= k
                    r% nl = nl

                    IF ( nl > nsum ) THEN
                      r% nsum = 2 * nsum
                    ENDIF

                    IF ( k < kp ) THEN
                      r% loop(nl) = k
                      r% link(k)  = nl
                      r% link(kp) = indx
                    ELSE
                      r% loop(nl) = kp
                      r% link(k)  = indx
                      r% link(kp) = nl
                    ENDIF

                    !=== Fix Links in New Loop ===!

                    mh = 1
                    ms = 0

                    ip = MIN(k,kp)
                    jp = ip + 1

                    jndx = r% link(ip)

                    DO WHILE ( jp < r% ibsp(ip) )

                      IF ( r% link(jp) == indx ) THEN
                        r% link(jp) = jndx
                      ENDIF

                      IF ( r% ibsp(jp) > jp ) mh = mh + 1
                      IF ( r% ibsp(jp) == 0 ) ms = ms + 1

                      IF ( r% ibsp(jp) > jp ) THEN
                        jp = r% ibsp(jp)
                      ELSE
                        jp = jp + 1
                      ENDIF

                    ENDDO

                    r% nhlx(indx) = nh - mh + 2
                    r% nsgl(indx) = ns - ms - 2

                    r% nhlx(jndx) = mh
                    r% nsgl(jndx) = ms

                    CALL LOOP_REAC (r,indx)
                    CALL LOOP_REAC (r,jndx)

                    !=== Recalc Lower Loop? ===!

                    IF ( iloop == 0 ) kndx = 0
                    IF ( iloop == 1 ) kndx = r% link(j)

                    IF ( kndx /= 0 ) CALL LOOP_REAC (r,kndx)

                    RETURN

                  ENDIF

                ENDIF

              ELSE

                l = l + 1
                kp = r% ibsp(kp)

              ENDIF

              l = l + 1
              kp = kp + 1

            ENDDO

          ENDIF
          atot = atot + r% wrk1(k)
          ENDIF


          !=== Helix Events ===!

          IF ( r% ibsp(k) > 0 ) THEN

            ip = k
            jp = r% ibsp(k)

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

              atot = atot + r% wrk2(ip)

              IF ( atot >= amax ) THEN

                IF ( icase == 1 ) THEN

                  jndx = r% link(ip)
                  r% nsgl(indx) = ns - 2

                ELSEIF ( icase == 2 ) THEN

                  jndx = r% link(i+2)

                  IF ( jndx == nl ) THEN
                    jndx = indx
                  ENDIF

                  !=== Delete Loop indx ===!

                  !=== Copy Loop nl to indx ===!

                  IF ( indx /= nl ) THEN

                    r% loop(indx) = r% loop(nl)
                    r% nhlx(indx) = r% nhlx(nl)
                    r% nsgl(indx) = r% nsgl(nl)
                    r% ptot(indx) = r% ptot(nl)

                    CALL LOOP_RESUM (r,indx)

                    kp = r% loop(nl)

                    r% link(kp) = indx

                    l = kp + 1

                    DO WHILE ( l < r% ibsp(kp) )

                      IF ( r% link(l) == nl ) THEN
                        r% link(l) = indx
                      ENDIF

                      IF ( r% ibsp(l) > l ) THEN
                        l = r% ibsp(l)
                      ELSE
                        l = l + 1
                      ENDIF

                    ENDDO

                  ENDIF

                  r% link(i)  = 0
                  r% link(j-2)= 0

                  r% loop(nl) = 0
                  r% nhlx(nl) = 0
                  r% nsgl(nl) = 0
                  r% ptot(nl) = 0.0d0

                  CALL LOOP_RESUM (r,nl)

                  r% nl = nl - 1

                  IF ( nsum > 2 ) THEN
                  IF ( r% nl <= nsum / 2 ) THEN
                    nsum = nsum / 2
                    r% nsum = nsum
                    r% psum(nsum) = 0.0d0
                  ENDIF
                  ENDIF

                ELSEIF ( icase == 3 ) THEN

                  r% loop(indx) = jp + 1
                  r% nsgl(indx) = ns - 2

                ENDIF

                !=== Adjust Base-Pairs ===!

                r% ibsp(ip-1) = jp+1
                r% ibsp(jp+1) = ip-1

                !=== Fix Links ===!

                r% link(jp+1) = r% link(jp)
                r% link(jp) = 0

                !=== Recalc Main Loop? ===!

                IF ( icase /= 2 ) THEN
                  CALL LOOP_REAC (r,indx)
                ENDIF

                !=== Recalc Upper Loop? ===!

                IF ( icase /= 3 ) THEN
                IF ( jndx == 0 ) CALL HELX_REAC (r,ip)
                IF ( jndx /= 0 ) CALL LOOP_REAC (r,jndx)
                ENDIF

                !=== Recalc Lower Loop? ===!

                IF ( iloop == 0 ) kndx = 0
                IF ( iloop == 1 ) kndx = r% link(j)

                IF ( kndx /= 0 ) CALL LOOP_REAC (r,kndx)

                RETURN

              ENDIF

            ENDIF

            !=== Helix Retraction ===!

            icase = 0

            IF ( ip /= n .and. jp /= 1 ) THEN

              IF ( r% ibsp(ip+1) == jp-1 ) THEN

                icase = 1

                IF ( iloop == 1 ) THEN
                IF ( k == ke ) icase = 2
                ENDIF

              ELSEIF ( iloop == 0 .or. k /= ke ) THEN

                l  = r% link(ip)
                mh = r% nhlx(l)
                ms = r% nsgl(l)

                icase = 3

              ENDIF

            ENDIF

            IF ( icase > 0 ) THEN

              IF ( icase == 2 ) THEN
                atot = atot + r% wrk1(ip)
              ELSE
                atot = atot + r% wrk1(jp)
              ENDIF

              IF ( atot >= amax ) THEN

                r% ibsp(ip) = 0
                r% ibsp(jp) = 0

                IF ( icase == 1 ) THEN

                  r% nsgl(indx) = ns + 2

                  r% link(jp-1) = r% link(jp)
                  IF ( jp /= n ) r% link(jp) = 0

                  !=== Recalc Main Loop ===!

                  CALL LOOP_REAC (r,indx)

                  !=== Recalc Upper/Lower Loops? ===!

                  jndx = r% link(ip+1)

                  IF ( iloop == 0 ) kndx = 0
                  IF ( iloop == 1 ) kndx = r% link(j)

                  IF ( jndx == 0 ) CALL HELX_REAC (r,ip+1)
                  IF ( jndx /= 0 ) CALL LOOP_REAC (r,jndx)
                  IF ( kndx /= 0 ) CALL LOOP_REAC (r,kndx)

                ELSEIF ( icase == 2 ) THEN

                  r% nsgl(indx) = ns + 2
                  r% loop(indx) = jp - 1

                  r% link(jp-1) = r% link(jp)
                  IF ( jp /= n ) r% link(jp) = 0

                  !=== Recalc Main Loop ===!

                  CALL LOOP_REAC (r,indx)

                  !=== Recalc Lower Loop? ===!

                  kndx = r% link(ip+1)

                  IF ( kndx /= 0 ) CALL LOOP_REAC (r,kndx)

                ELSEIF ( icase == 3 ) THEN

                  r% nhlx(indx) = nh + mh - 2
                  r% nsgl(indx) = ns + ms + 2

                  jndx = r% link(ip)

                  !=== Fix Links in Loop Being Deleted ===!

                  l = ip + 1

                  DO WHILE ( l < jp )

                    IF ( r% link(l) == jndx ) THEN
                      r% link(l) = indx
                    ENDIF

                    IF ( r% ibsp(l) > l ) THEN
                      l = r% ibsp(l)
                    ELSE
                      l = l + 1
                    ENDIF

                  ENDDO

                  !=== Copy Loop nl to jndx ===!

                  IF ( jndx /= nl ) THEN

                    r% loop(jndx) = r% loop(nl)
                    r% nhlx(jndx) = r% nhlx(nl)
                    r% nsgl(jndx) = r% nsgl(nl)
                    r% ptot(jndx) = r% ptot(nl)

                    CALL LOOP_RESUM (r,jndx)

                    kp = r% loop(nl)

                    r% link(kp) = jndx

                    l = kp + 1

                    DO WHILE ( l < r% ibsp(kp) )

                      IF ( r% link(l) == nl ) THEN
                        r% link(l) = jndx
                      ENDIF

                      IF ( r% ibsp(l) > l ) THEN
                        l = r% ibsp(l)
                      ELSE
                        l = l + 1
                      ENDIF

                    ENDDO

                  ENDIF

                  IF ( indx /= nl ) jndx = indx

                  r% link(ip) = 0
                  IF ( jp /= n ) r% link(jp) = 0

                  r% loop(nl) = 0
                  r% nhlx(nl) = 0
                  r% nsgl(nl) = 0
                  r% ptot(nl) = 0.0d0

                  CALL LOOP_RESUM (r,nl)

                  r% nl = nl - 1

                  IF ( nsum > 2 ) THEN
                  IF ( r% nl <= nsum / 2 ) THEN
                    nsum = nsum / 2
                    r% nsum = nsum
                    r% psum(nsum) = 0.0d0
                  ENDIF
                  ENDIF

                  !=== Recalc Main Loop ===!

                  CALL LOOP_REAC (r,jndx)

                  !=== Recalc Lower Loop? ===!

                  IF ( iloop == 0 ) kndx = 0
                  IF ( iloop == 1 ) kndx = r% link(j)

                  IF ( kndx /= 0 ) CALL LOOP_REAC (r,kndx)

                ENDIF

                RETURN

              ENDIF

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

              IF ( atot >= amax ) THEN

                is = r% ibsp(ip-1)
                js = r% ibsp(jp+1)

                IF ( is /= 0 ) THEN

                  r% ibsp(is)  = 0
                  r% link(ip-1)= 0
                  r% link(ip-2)= indx

                  IF ( iloop == 1 .and. is == ke ) THEN
                    r% loop(indx) = ip - 2
                  ENDIF

                ENDIF

                IF ( js /= 0 ) THEN

                  r% ibsp(js)  = 0
                  r% link(js)  = 0
                  r% link(js-1)= indx

                  IF ( iloop == 1 .and. js ==  i ) THEN
                    r% loop(indx) = js - 1
                  ENDIF

                ENDIF

                IF ( is /= 0 .and. js /= 0 ) THEN
                  ns = ns + 2
                  r% nsgl(indx) = ns
                ENDIF

                !=== Adjust Base-Pairs ===!

                r% ibsp(ip-1) = jp+1
                r% ibsp(jp+1) = ip-1

                !=== Fix Links ===!

                r% link(jp+1) = r% link(jp)
                r% link(jp) = 0

                IF ( iloop == 1 .and. k == ke ) THEN
                  r% loop(indx) = jp + 1
                ENDIF

                !=== Recalc Main Loop ===!

                CALL LOOP_REAC (r,indx)

                jndx = r% link(ip)
                IF ( jndx == 0 ) CALL HELX_REAC (r,ip)
                IF ( jndx /= 0 ) CALL LOOP_REAC (r,jndx)

                !=== Recalc 5' 3' Loops? ===!

                IF ( is /= 0 ) THEN
                  jndx = r% link(is+1)
                  IF ( jndx == 0 ) CALL HELX_REAC (r,is+1)
                  IF ( jndx /= 0 ) CALL LOOP_REAC (r,jndx)
                ENDIF

                IF ( js /= 0 ) THEN
                  jndx = r% link(jp+2)
                  IF ( jndx == 0 ) CALL HELX_REAC (r,jp+2)
                  IF ( jndx /= 0 ) CALL LOOP_REAC (r,jndx)
                ENDIF

                !=== Recalc Lower Loop? ===!

                IF ( iloop == 0 ) kndx = 0
                IF ( iloop == 1 ) kndx = r% link(j)

                IF ( kndx /= 0 ) CALL LOOP_REAC (r,kndx)

                RETURN

              ENDIF

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

                  IF ( atot >= amax ) GOTO 1

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

                  IF ( atot >= amax ) GOTO 1

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

                  IF ( atot >= amax ) GOTO 1

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

                  IF ( atot >= amax ) GOTO 1

                ENDIF

              ENDIF

            ENDIF

 1          IF ( atot >= amax ) THEN

              !=== Adjust Base-Pairs ===!

              IF ( kp == ip+1 .or. kp == ip-1 ) THEN
                r% ibsp(ip) = 0
                r% ibsp(jp) = kp
                r% ibsp(kp) = jp
              ELSE
                r% ibsp(ip) = kp
                r% ibsp(jp) = 0
                r% ibsp(kp) = ip
              ENDIF

              IF ( icase == 1 ) THEN

                jndx = r% link(ip)

                IF ( kp == ip+1 .or. kp == ip-1 ) THEN
                  r% loop(jndx) = kp
                  r% link(kp) = r% link(ip)
                  r% link(ip) = 0
                ELSE
                  r% link(kp) = r%link(jp)
                  r% link(jp) = 0
                ENDIF

                IF ( kp == ip+1 .or. kp == jp-1 ) THEN
                  r% nsgl(indx) = r% nsgl(indx) + 1
                  r% nsgl(jndx) = r% nsgl(jndx) - 1
                ELSE
                  r% nsgl(indx) = r% nsgl(indx) - 1
                  r% nsgl(jndx) = r% nsgl(jndx) + 1
                ENDIF

                CALL LOOP_REAC (r,indx)
                CALL LOOP_REAC (r,jndx)

                !=== Recalc Lower Loop? ==!

                IF ( iloop == 0 ) kndx = 0
                IF ( iloop == 1 ) kndx = r% link(j)

                IF ( kndx /= 0 ) CALL LOOP_REAC (r,kndx)

              ELSEIF ( icase == 2 ) THEN

                !=== Add Loop ===!

                nl = nl + 1

                r% nl = nl

                r% link(jp-1) = nl

                IF ( nl > nsum ) THEN
                  r% nsum = 2 * nsum
                ENDIF

                IF ( kp == jp+1 ) THEN
                  r% link(ip) = nl
                  r% link(jp+1) = r% link(jp)
                  r% link(jp) = 0
                ELSE
                  r% link(kp) = nl
                ENDIF

                IF ( iloop == 1 .and. k == ke ) THEN
                  r% loop(nl) = i - 1
                  IF ( kp == jp+1 ) THEN
                  r% loop(indx) = i + 1
                  ENDIF
                ELSEIF ( kp == ip-1 ) THEN
                  r% loop(nl) = kp
                ELSE
                  r% loop(nl) = ip
                ENDIF

                r% nsgl(indx) = ns - 1
                r% nhlx(nl) = 2
                r% nsgl(nl) = 1

                CALL LOOP_REAC (r,indx)
                CALL LOOP_REAC (r,nl)

                !=== Recalc Upper Loop? ===!

                jndx = r% link(ip+1)
                IF ( jndx /= 0 ) CALL LOOP_REAC (r,jndx)

                IF ( iloop == 0 .or. k /= ke ) THEN
                IF ( jndx == 0 ) CALL HELX_REAC (r,ip+1)
                ENDIF

                !=== Recalc Lower Loop? ===!

                IF ( iloop == 0 ) kndx = 0
                IF ( iloop == 1 ) kndx = r% link(j)

                IF ( kndx /= 0 ) CALL LOOP_REAC (r,kndx)

              ELSEIF ( icase == 3 ) THEN

                !=== Move Loop ===!

                r% link(jp-1) = r% link(jp)
                r% link(jp) = 0

                IF ( kp == jp+1 ) THEN
                  r% link(ip) = r% link(ip-1)
                  r% link(ip-1) = 0
                ELSE
                  r% link(ip-1) = r% link(ip-2)
                  r% link(ip-2) = 0
                ENDIF

                IF ( k == ke ) THEN

                  r% loop(indx) = i-1

                  kndx = r% link(ip+1)

                  IF ( kp == ip-1 ) jndx = r% link(i+1)
                  IF ( kp == jp+1 ) jndx = r% link(i+2)

                  IF ( jndx == 0 ) CALL HELX_REAC (r,j-1)

                ELSE

                  r% loop(indx) = i+1

                  jndx = r% link(ip+1)
                  kndx = r% link(j)

                  IF ( jndx == 0 ) CALL HELX_REAC (r,ip+1)

                ENDIF

                CALL LOOP_REAC (r,indx)

                IF ( jndx /= 0 ) CALL LOOP_REAC (r,jndx)
                IF ( kndx /= 0 ) CALL LOOP_REAC (r,kndx)

              ELSEIF ( icase == 4 ) THEN

                !=== Delete Loop ===!

                IF ( kp == ip+1 ) THEN

                  jndx = r% link(ip)

                  r% nsgl(indx) = r% nsgl(indx) + 1

                  r% link(ip)  = 0
                  r% link(jp-1)= 0

                  CALL LOOP_REAC (r,indx)

                  !=== Recalc Upper Loop? ===!

                  kndx = r% link(ip+2)

                  IF ( kndx /= 0 ) CALL LOOP_REAC (r,kndx)
                  IF ( kndx == 0 ) CALL HELX_REAC (r,ip+2)

                  !=== Recalc Lower Loop? ===!

                  IF ( iloop == 0 ) kndx = 0
                  IF ( iloop == 1 ) kndx = r% link(j)

                  IF ( kndx /= 0 ) CALL LOOP_REAC (r,kndx)

                ELSEIF ( kp == jp-1 ) THEN

                  jndx = r% link(ip)

                  r% nsgl(indx) = r% nsgl(indx) + 1

                  r% link(jp-1)= r% link(jp)
                  r% link(ip)  = 0
                  r% link(jp)  = 0
                  r% link(jp-2)= 0

                  CALL LOOP_REAC (r,indx)

                  !=== Recalc Upper Loop? ===!

                  kndx = r% link(ip+1)

                  IF ( kndx /= 0 ) CALL LOOP_REAC (r,kndx)
                  IF ( kndx == 0 ) CALL HELX_REAC (r,ip+1)

                  !=== Recalc Lower Loop? ===!

                  IF ( iloop == 0 ) kndx = 0
                  IF ( iloop == 1 ) kndx = r% link(j)

                  IF ( kndx /= 0 ) CALL LOOP_REAC (r,kndx)

                ELSE

                  jndx = r% link(ip)
                  kndx = r% link(j)

                  r% nsgl(jndx) = r% nsgl(jndx) + 1

                  r% link(i) = 0
                  r% link(jp)= 0

                  IF ( kp == ip-1 ) THEN
                    r% loop(jndx) = ip - 1
                    r% link(ip-1) = r% link(ip)
                    r% link(ip) = 0
                  ENDIF

                  !=== Recalc Loops ===!

                  CALL LOOP_REAC (r,jndx)

                  IF ( kndx /= 0 ) CALL LOOP_REAC (r,kndx)

                  jndx = indx

                ENDIF

                !=== Delete Loop jndx ===!

                !=== Copy Loop nl to indx ===!

                IF ( jndx /= nl ) THEN

                  r% loop(jndx) = r% loop(nl)
                  r% nhlx(jndx) = r% nhlx(nl)
                  r% nsgl(jndx) = r% nsgl(nl)
                  r% ptot(jndx) = r% ptot(nl)

                  CALL LOOP_RESUM (r,jndx)

                  is = r% loop(nl)

                  r% link(is) = jndx

                  l = is + 1

                  DO WHILE ( l < r% ibsp(is) )

                    IF ( r% link(l) == nl ) THEN
                      r% link(l) = jndx
                    ENDIF

                    IF ( r% ibsp(l) > l ) THEN
                      l = r% ibsp(l)
                    ELSE
                      l = l + 1
                    ENDIF

                  ENDDO

                ENDIF

                r% loop(nl) = 0
                r% nhlx(nl) = 0
                r% nsgl(nl) = 0
                r% ptot(nl) = 0.0d0

                CALL LOOP_RESUM (r,nl)

                r% nl = nl - 1

                IF ( nsum > 2 ) THEN
                IF ( r% nl <= nsum / 2 ) THEN
                  nsum = nsum / 2
                  r% nsum = nsum
                  r% psum(nsum) = 0.0d0
                ENDIF
                ENDIF

                !=== Recalc Main Loop? ===!

                IF ( jndx /= indx ) CALL LOOP_REAC (r,indx)

              ENDIF

              RETURN

            ENDIF

            !=== Open BP Inside Helix ===!

            IF ( r% link(ip) == 0 ) THEN
            IF ( iloop == 1 .and. k == ke ) THEN

              is = ip + 1
              js = jp - 1

              DO WHILE ( r% link(is) == 0 )

                atot = atot + r% wrk1(is)

                is = is + 1
                js = js - 1

                IF ( atot >= amax ) THEN

                  nl = nl + 1

                  r% ibsp(is-1) = 0
                  r% ibsp(js+1) = 0
                  r% nl = nl

                  r% loop(nl) = js
                  r% link(js) = nl
                  r% link(is-2) = nl

                  r% nhlx(nl) = 2
                  r% nsgl(nl) = 2

                  IF ( nl > nsum ) THEN
                    r% nsum = 2 * nsum
                  ENDIF

                  !=== Recalc Upper Loop? ===!

                  jndx = r% link(js+2)

                  IF ( jndx == 0 ) CALL HELX_REAC (r,js+2)
                  IF ( jndx /= 0 ) CALL LOOP_REAC (r,jndx)

                  !=== Calculate New Loop ===!

                  jndx = r% link(js)

                  CALL LOOP_REAC (r,jndx)

                  !=== Recalc Lower Loop? ===!

                  kndx = r% link(is)

                  IF ( kndx /= 0 ) CALL LOOP_REAC (r,kndx)

                  RETURN

                ENDIF

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

        RETURN

      END SUBROUTINE LOOP_FIRE
