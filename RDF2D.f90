!     A simple 2D example for computing the Radial Distribution Function, g(r),
!     of a random (almost uniform) distribution of points on a 2D square plate.
!     Developed by Ahmad Ababaei (ahmad.ababaei@imgw.pl)
!     Institute of Meteorology and Water Management — National Research Institute
!     Podleśna 61, 01-673 Warsaw, Poland

      PROGRAM RadialDistributionFunction
      REAL, ALLOCATABLE, DIMENSION(:) :: x,y,rdf

      plt = 1.0  ! square plate size
      nsh = 10   ! no. of shells (rings)
      dsc = 0.1  ! scan radius: 0 < r < dsc
      dsh = dsc / real(nsh) ! shell thickness
       pi = 4.0 * atan(1.0)
      ALLOCATE (rdf(nsh))
      rdf = 0.0

!---- READ -------------------------------
!     OPEN (1, FILE = 'input.dat')
!     np = 0
!     DO
!       READ (1,*, END=1)
!     np = np + 1
!     ENDDO
!1    CLOSE (1)
!     dnp = real(np)*real(np-1)/2.0/plt**2
!     WRITE(*,*) 'Number of points:', np
!     WRITE(*,*) 'Pair density:',  dnp
!     nsi = int ( 3.0 * 4.0 * real(np) / plt**2 * dsc * plt )
!     nco = int ( 3.0 * 4.0 * real(np) / plt**2 * dsc * dsc )
!     WRITE(*,*) "Allocation size:", np+nsi+nco
!     ALLOCATE ( x(np+nsi+nco), y(np+nsi+nco) )
!     ne = np
!     OPEN (1, FILE = 'input.dat')
!     DO i = 1, np
!        READ (1,*) x(i), y(i)
!     ENDDO
!     CLOSE (1)
!     GOTO 4
!-----------------------------------------

!---- RANDOM DISTRIBUTION ----------------
       np = 10000
      dnp = real(np)*real(np-1)/2.0/plt**2
      WRITE(*,*) 'Number of points:', np
      WRITE(*,*) 'Pair density:',  dnp
      nsi = int ( 3.0 * 4.0 * real(np) / plt**2 * dsc * plt )
      nco = int ( 3.0 * 4.0 * real(np) / plt**2 * dsc * dsc )
      WRITE(*,*) "Allocation size:", np+nsi+nco
      ALLOCATE ( x(np+nsi+nco), y(np+nsi+nco) )
      ne = np

      OPEN (3, FILE = 'distribution.dat')
      DO n = 1, np
         CALL RANDOM_NUMBER(r)
!        r = r**2 ! non-uniform (disable for a uniform) distribution
         x(n) = r * plt
         CALL RANDOM_NUMBER(r)
!        r = r**2 ! non-uniform (disable for a uniform) distribution
         y(n) = r * plt
         WRITE(3,*) x(n), y(n)

!------- EXTENSION OF DISTRIBUTION -------
! For the points whose distances to the edges of the plate are less than the scanning radius, there is a point-free region (outside the square boundaries) where a neighboring point never exists. Computing RDF would result in g(r) being a decreasing function of r. If we are to judge the distribution in the domain with no influence from the boundaries (for example, when we have a sample distribution from a larger domain), the remedy is to periodically repeat the distribution in the domain. In such a case, an ideal uniform distribution would result in g(r) = 1. Since the scanning radius is limited, in order to periodically extend the domain, we only need to repeat the parts within a distance of scanning radius next to the boundaries. This is as if we have infinitely many domains next to the current one from all eight cardinal directions.

         IF ( x(n) .LT. dsc ) THEN ! W
               ne = ne + 1
            x(ne) = x(n) + plt
            y(ne) = y(n)
            IF ( y(n) .LT. dsc ) THEN ! SW
                  ne = ne + 1
               x(ne) = x(n) + plt
               y(ne) = y(n) + plt
            ELSEIF ( y(n) .GT. (plt-dsc) ) THEN ! NW
                  ne = ne + 1
               x(ne) = x(n) + plt
               y(ne) = y(n) - plt
            ENDIF
         ELSEIF ( x(n) .GT. (plt-dsc) ) THEN ! E
               ne = ne + 1
            x(ne) = x(n) - plt
            y(ne) = y(n)
            IF ( y(n) .LT. dsc ) THEN ! SE
                  ne = ne + 1
               x(ne) = x(n) - plt
               y(ne) = y(n) + plt
            ELSEIF ( y(n) .GT. (plt-dsc) ) THEN ! NE
                  ne = ne + 1
               x(ne) = x(n) - plt
               y(ne) = y(n) - plt
            ENDIF
         ENDIF

         IF ( y(n) .LT. dsc ) THEN ! S
               ne = ne + 1
            x(ne) = x(n)
            y(ne) = y(n) + plt
         ELSEIF ( y(n) .GT. (plt-dsc) ) THEN ! N
               ne = ne + 1
            x(ne) = x(n)
            y(ne) = y(n) - plt
         ENDIF
      ENDDO

      WRITE(*,*) "np extended:", ne
      WRITE(*,*)
      OPEN (4, FILE = 'extension.dat')
      DO n = np+1, ne
         WRITE(4,*) x(n), y(n)
      ENDDO
      CLOSE (4)
!-----------------------------------------

!---- RADIAL DISTANCE FUNCTION -----------
4     DO i = 1, np
         DO j = i+1, np ! Pairs (i,j) both inside the square
            dst = sqrt ( ( x(j)-x(i) )**2 + ( y(j)-y(i) )**2 )
            IF ( dst .GE. dsc ) GOTO 2
            ish = int( dst / dsc * real (nsh) ) + 1
            rdf(ish) = rdf(ish) + 1.0
2        ENDDO
         DO j = np+1, ne ! Pairs (i,j): i ‒> inside & j ‒> outside (extended)
            dst = sqrt ( ( x(j)-x(i) )**2 + ( y(j)-y(i) )**2 )
            IF ( dst .GE. dsc ) GOTO 3
            ish = int( dst / dsc * real (nsh) ) + 1
            rdf(ish) = rdf(ish) + 0.5 ! half-counted because i-j is counted but j-i is not!
3        ENDDO
      ENDDO
      DEALLOCATE (x,y)
!-----------------------------------------

!------ N O R M A L I Z A T I O N --------
!---- RADIAL DISTRIBUTION FUNCTION -------
      rdf_av = 0.0
      OPEN (2, FILE = 'r_g.dat')
      ri = 0.0
      WRITE(*,*) " r:               g(r):"
      DO ish = 1, nsh
          ro = ri + dsh
          sr = pi * ( ro**2 - ri**2 )    ! surface area of rings (shells)
          rdf(ish) = rdf(ish) / sr / dnp ! normalization
          rdf_av = rdf_av + rdf(ish) * dsh
          WRITE(2,*) (ri+ro)/2.0, rdf(ish)
          WRITE(*,*) (ri+ro)/2.0, rdf(ish)
          ri = ro
      ENDDO
      CLOSE (2)
      DEALLOCATE (rdf)
      WRITE(*,*)
      WRITE(*,*) 'Average RDF =', rdf_av / dsc
!-----------------------------------------

      END
