       program inzone

c Program to compute innermost zone effect on gravity anomaly or geoid   
c due to geoid gradient.
c
c
c -N: north-south geoid gradient (in micro-rad)
c -E: west-east geoid gradient (in micro-rad) 
! -G: output .grd3 file
c -O: for geoidal height [default: for gravity anomalies]
      implicit none
      integer*4 idim2,type,i,j,ix,iy,nx0,ny0,ldf,kx,ky,
     1 iright,ileft,jdown,jup,ii,jj,k,span
      parameter (idim2=34*60+1,span=6,ldf=2*span+1)
      real*4 north(idim2,idim2),east(idim2,idim2),out(idim2,idim2),
     1 xa(ldf),ya(ldf),data1(ldf,ldf),data2(ldf,ldf),xi_dy,eta_dx,gama,
     2 qd2dr
      real*8 tmp,phi,rl,pi,sqrt_pi,re,dtor,gnorm,factor,dx,dy,s,s2,
     1 lata(idim2),dlat(idim2),time0,time1,
     2 xmin,xmax,ymin,ymax,gridx,gridy,west1,east1,south1,north1,gdx,gdy

c variables for input arguments
      integer*4 nargs,iargc,par
      logical io(10),check
      character*80 cha,file_n, file_e,fileanom
      data par/3/
ccccccccccccccccccccccccccc

c First executable
      call second(time0)
C get arguments
      nargs=iargc()
      do i=1,par
      io(i)=.true.
      end do
      type=0
      i=0
!      do ii=1,nargs
!      call getarg(ii,cha)
!      if(cha(1:1).eq.'-') then
!
!		if(cha(2:2).eq.'N' .or. cha(2:2).eq.'n') then
!		file_n=cha(3:)
!		i=i+1
!                io(i)=.false.
!
!		elseif(cha(2:2).eq.'E' .or. cha(2:2).eq.'e') then
!		file_e=cha(3:)
!		i=i+1
!                io(i)=.false.
!
!		elseif(cha(2:2).eq.'G' .or. cha(2:2).eq.'g') then
!		fileanom=cha(3:)
!		i=i+1
!                io(i)=.false.
!
!
!		elseif(cha(2:2).eq.'O' .or. cha(2:2).eq.'o') then
!		type=1
!		else
!	        call write_error
!                stop
!		end if
!      else
!      call write_error
!      stop
!      end if
!      end do

! check necessary parameters
!      do i=1,par
!      if(io(i)) then
!      call write_error
!      stop
!      end if
!      end do
      file_n="re_gradients_B6_L12_n.grd3"
      file_e="re_gradients_B6_L12_e.grd3"
      fileanom="dg_inzone_B6_L12.grd3"
C Define some constants
      re=6371790.0d0 ! mean earth radius in meter
      pi=datan(1.d0)*4.d0
      sqrt_pi=dsqrt(pi)
      dtor=pi/180.0
      open(10,file=file_n,form='unformatted',status='old')
      read(10)ix,iy,xmin,xmax,ymin,ymax,gridx,gridy
      nx0=ix
      ny0=iy
      if(nx0.gt.idim2) stop'column out of bound'
      if(ny0.gt.idim2) stop'row out of bound'
      open(11,file=file_e,form='unformatted',status='old')
! Output file
      open(60,file=fileanom,form='unformatted')
C read north component
      read(10)((north(i,j),j=1,iy),i=1,ix)
c read east component
      read(11)ix,iy,west1,east1,south1,north1,gdx,gdy
      read(11)((east(i,j),j=1,iy),i=1,ix)
      if (ix.ne.nx0)stop'x dim not mathced'
      if (iy.ne.ny0)stop'y dim not mathced'
      if(west1.ne.xmin) stop'west not matched'
      if(east1.ne.xmax) stop'east not matched'
      if(south1.ne.ymin) stop'south not matched'
      if(north1.ne.ymax) stop'north not matched'
      if(gridx.ne.gdx)stop'grid size x not mathched'
      if(gridy.ne.gdy)stop'grid size y not mathched'

      check=.false.
      do j=1,ny0
      lata(j)=(ymin+dfloat(j-1)*gridy)*dtor
      dlat(j)=gridy*dtor
      end do

      do 1 jj=1,iy
      phi=lata(jj)
      if(type.eq.0)gnorm=gama(phi)
      factor=dcos(phi)*re*dtor
c dx and dy are in meters
      dy=dlat(jj)*re
      dx=gridx*factor
      s=dsqrt(dx*dy)/sqrt_pi
      s2=s*s
      jdown=jj-span
      if(jdown.lt.1)jdown=1
      jup=jj+span
      if(jup.gt.iy)jup=iy
      ky=jup-jdown+1
      k=0
      do j=jdown, jup
      k=k+1
      tmp=lata(j)
      ya(k)=(tmp-phi)*re
      end do

      do 1 ii=1,ix
      rl=xmin+(ii-1)*gridx

      ileft=ii-span
      if(ileft.lt.1)ileft=1
      iright=ii+span
      if(iright.gt.ix)iright=ix
      kx= iright-ileft+1
      k=0
      do i=ileft, iright
      k=k+1
      tmp=xmin+(i-1)*gridx
      xa(k)=(tmp-rl)*factor
      end do

      do j=jdown, jup
      do i=ileft, iright
      data1(i-ileft+1, j-jdown+1)=north(i,j)
      data2(i-ileft+1, j-jdown+1)=east(i,j)
      end do
      end do

      xi_dy= qd2dr(0, 1, 0.0, 0.0, kx, xa, ky,
     & ya, data1, ldf, check)
      eta_dx= qd2dr(1, 0, 0.0, 0.0, kx, xa, ky,
     & ya, data2, ldf, check)
      if(type.eq.0) then
      tmp=-0.5*gnorm*1.e-6*s*(xi_dy+eta_dx)
      else
      tmp=-0.25*s2*1.e-6*(xi_dy+eta_dx)
      end if
      out(ii,jj)=tmp
1     continue  
c output innermost zone effect on gravity anomalies or geoidal heights
      write(60)ix,iy,xmin,xmax,ymin,ymax,gridx,gridy
      write(60)((out(i,j),j=1,iy),i=1,ix)
      call second(time1)
      write(6,103)time1-time0
 103  format('total execution time (seconds):',F20.3)
      END

      subroutine second(time)
      real*8 time
      real*4 T(2)
      data TOT/0.D0/
      TIME=DTIME(T)+TOT
      TOT=TIME
      RETURN
      END

C-----------------------------------------------------------------------
C
C  Purpose:    Evaluate the derivative of a function defined on a
C              rectangular grid using quadratic interpolation.
C
C  Usage:      QD2DR(IXDER, IYDER, X, Y, NXDATA, XDATA, NYDATA, YDATA,
C                    FDATA, LDF, CHECK)
C
C  Arguments:
C     IXDER  - Order of the X-derivative.  (Input)
C     IYDER  - Order of the Y-derivative.  (Input)
C     X      - X-coordinate of the point at which the function is to
C              be evaluated.  (Input)
C     Y      - Y-coordinate of the point at which the function is to
C              be evaluated.  (Input)
C     NXDATA - Number of data points in the X-direction.  (Input)
C              NXDATA must be at least 3.
C     XDATA  - Array of length NXDATA containing the location of the
C              data points in the X-direction.  (Input)
C              XDATA must be increasing.
C     NYDATA - Number of data points in the Y-direction.  (Input)
C              NYDATA must be at least 3.
C     YDATA  - Array of length NYDATA containing the location of the
C              data points in the Y-direction.  (Input)
C              YDATA must be increasing.
C     FDATA  - Array of size NXDATA by NYDATA containing function
C              values.  (Input)
C              FDATA(I,J) is the value of the function at
C              (XDATA(I),YDATA(J)).
C     LDF    - Leading dimension of FDATA exactly as specified in the
C              dimension statement of the calling program.  (Input)
C              LDF must be at least as large as NXDATA.
C     CHECK  - Logical variable that is .TRUE. if checking of XDATA and
C              YDATA is required or .FALSE. if checking is not required.
C              (Input)
C     QD2DR  - Value of the (IXDER,IYDER) derivative of the function
C              at (X,Y).  (Output)
C
C  Remark:
C     Because quadratic interpolation is used, if the order of any
C     derivative is greater than 2 then the returned value is zero.
C
C
C-----------------------------------------------------------------------
C
      REAL FUNCTION QD2DR (IXDER, IYDER, X, Y, NXDATA, XDATA, NYDATA,
     &      YDATA, FDATA, LDF, CHECK)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    IXDER, IYDER, NXDATA, NYDATA, LDF
      REAL       X, Y, XDATA(*), YDATA(*), FDATA(LDF,*)
      LOGICAL    CHECK
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, IC, J, JC
      REAL       C1, C2, C3, C4, C5, DX0, DXB, DY0, DYB, F0, F1, F2,
     &           F3, F4, FC, HXA, HXB, HXC, HXCB, HYA, HYB, HYC, HYCB,
     &           T1, T3, VALUE, X0, XA, XB, XC, Y0, YA, YB, YC
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  INT,SIGN
      INTRINSIC  INT, SIGN
      INTEGER    INT
      REAL       SIGN
C                                  SPECIFICATIONS FOR SUBROUTINES
c      EXTERNAL   E1MES, E1POP, E1PSH, E1STI, E1STR
C                                  SPECIFICATIONS FOR FUNCTIONS
c      EXTERNAL   N1RTY, Q3DER
      INTEGER    Q3DER
C
c      CALL E1PSH ('QD2DR ')
C                                  Set default VALUE
      VALUE = 0.0
C                                  Check NXDATA
      IF (NXDATA .LT. 3) stop'NXDATA must be greater than 3'
C                                  Check NYDATA
      IF (NYDATA .LT. 3) stop'NYDATA must be greater than 3'
C                                  Check IXDER
      IF (IXDER .LT. 0) stop'ixder must be greater than 0'
      IF (IYDER .LT. 0) stop'iyder must be greater than 0'
C                                  Check LDF
      IF (LDF .LT. NXDATA) THEN
      write(6,*) 'The leading dimension of FDATA must be ',
     &               'at least as large as NXDATA '
      stop
      END IF
C                                  Check for errors
c      IF (N1RTY(0) .NE. 0) GO TO 9000
      IF (CHECK) THEN
         DO 10  I=2, NXDATA
            IF (XDATA(I-1) .GE. XDATA(I)) 
     &stop'Points in the X direction not strictly increasing'
   10    CONTINUE
         DO 20  I=2, NYDATA
            IF (YDATA(I-1) .GE. YDATA(I))
     & stop'Points in the Y direction not strictly increasing'
   20    CONTINUE
      END IF
C                                   If IX+IY > 2 then VALUE=0.0
      IF (IXDER+IYDER .GT. 2) GO TO 9000
C
C                    F3    FC       F0, F1, F2, F3 are the values
C                     +             at the grid points.  X is the
C                     + X           point at which the function
C              F2++++F0++++F1       is to be estimated. (It need
C                     +             not be in the First quadrant).
C                     +             FC - the outer corner point
C                    F4             nearest X.
C
C                                   F0 is the value of the FDATA at
C                                   FDATA(I,J), it is the interior mesh
C                                   point nearest  X.
C                                   The coordinates of F0 are (X0,Y0),
C                                   The coordinates of F1 are (XB,Y0),
C                                   The coordinates of F2 are (XA,Y0),
C                                   The coordinates of F3 are (X0,YB),
C                                   The coordinates of F4 are (X0,YA),
C                                   The coordinates of FC are (XC,YC),
C
      I = Q3DER(X,NXDATA,XDATA)
      J = Q3DER(Y,NYDATA,YDATA)
C
      X0 = XDATA(I)
      Y0 = YDATA(J)
      XA = XDATA(I-1)
      XB = XDATA(I+1)
      YA = YDATA(J-1)
      YB = YDATA(J+1)
C
      DX0 = X - X0
      DY0 = Y - Y0
C
      DXB = X - XB
      DYB = Y - YB
C
      F0 = FDATA(I,J)
      F1 = FDATA(I+1,J)
      F2 = FDATA(I-1,J)
      F3 = FDATA(I,J+1)
      F4 = FDATA(I,J-1)
C
      IC = I + INT(SIGN(1.0,DX0))
      JC = J + INT(SIGN(1.0,DY0))
C
      XC = XDATA(IC)
      YC = YDATA(JC)
      FC = FDATA(IC,JC)
C                       O               HXA, HXB are the mesh spacings
C                       +               in the X-direction to the left
C                      HYB              and right of the center point.
C                       +
C               O++HXA++O++HXB++O       HYB, HYA are the mesh spacings
C                       +               in the Y-direction.
C                      HYA
C                       +               HXC equals either  HXB  or  HXA
C                       O               depending on where the corner
C                                       point is located.
C
      HXA = X0 - XA
      HXB = XB - X0
      HYA = Y0 - YA
      HYB = YB - Y0
C
      HXC = XC - X0
      HYC = YC - Y0
C
      HXCB = XC - XB
      HYCB = YC - YB
C                                       Construct the interpolant
C                                       F = F0 + C1*(X-X0) +
C                                           C2*(X-X0)*(X-X1) +
C                                           C3*(Y-Y0) + C4*(Y-Y0)*(Y-Y1)
C                                           + C5*(X-X0)*(Y-Y0)
      C1 = (F1-F0)/HXB
      T1 = (F0-F2)/HXA
      C2 = (C1-T1)/(HXA+HXB)
      C3 = (F3-F0)/HYB
      T3 = (F0-F4)/HYA
      C4 = (C3-T3)/(HYA+HYB)
      C5 = (FC-F0-HXC*C1-HXC*HXCB*C2-HYC*C3-HYC*HYCB*C4)/(HXC*HYC)
C
      IF (IXDER.EQ.0 .AND. IYDER.EQ.0) THEN
         VALUE = F0 + DX0*(C1+DXB*C2+DY0*C5) + DY0*(C3+DYB*C4)
      ELSE IF (IXDER.EQ.1 .AND. IYDER.EQ.0) THEN
         VALUE = C1 + (DXB+DX0)*C2 + DY0*C5
      ELSE IF (IXDER.EQ.0 .AND. IYDER.EQ.1) THEN
         VALUE = C3 + (DYB+DY0)*C4 + DX0*C5
      ELSE IF (IXDER.EQ.1 .AND. IYDER.EQ.1) THEN
         VALUE = C5
      ELSE IF (IXDER.EQ.2 .AND. IYDER.EQ.0) THEN
         VALUE = 2.0*C2
      ELSE IF (IXDER.EQ.0 .AND. IYDER.EQ.2) THEN
         VALUE = 2.0*C4
      ELSE
         VALUE = 0.0
      END IF
C
 9000 CONTINUE
      QD2DR = VALUE
      RETURN
      END
C-----------------------------------------------------------------------
C
C  Purpose:    Find the index of the point in DATA(2:NDATA-1) closest
C              to X.
C
C  Usage:      Q3DER(X, NDATA, DATA)
C
C  Arguments:
C     X      - A given point.  (Input)
C     NDATA  - Number of points in DATA.  (Input)
C              It must be at least 3.
C     DATA   - Array of data points.  (Input)
C              DATA must be strictly increasing.
C     Q3DER  - The index of the point in DATA(2:NDATA-1) closest to X.
C              (Output)
C
C-----------------------------------------------------------------------
C
      INTEGER FUNCTION Q3DER (X, NDATA, DATA)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    NDATA
      REAL       X, DATA(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    IVAL
C                                  SPECIFICATIONS FOR INTRINSICS
C     INTRINSIC  MAX0,MIN0
      INTRINSIC  MAX0, MIN0
      INTEGER    MAX0, MIN0
C                                  SPECIFICATIONS FOR SUBROUTINES
c      EXTERNAL   P3DER
C
      CALL P3DER (1, NDATA-1, DATA, X, IVAL)
C
      IF (IVAL .EQ. 1) IVAL = 2
      IF ((X-DATA(IVAL-1)) .LT. (DATA(IVAL)-X)) IVAL = IVAL - 1
C
      Q3DER = MAX0(MIN0(IVAL,NDATA-1),2)
C
      RETURN
      END
C-----------------------------------------------------------------------
C
C  Purpose:    Compute MAX (I, 1 .LE. NINTV. .AND. BREAK(I) .LE. X)
C
C  Usage:      CALL P3DER (KORD, NINTV, BREAK, X, LEFT)
C
C  Arguments:
C     KORD   - Order of the polynomial.  (Input)
C     NINTV  - Number of polynomial pieces.  (Input)
C     BREAK  - Vector of length NINTV+1 containing the breakpoints
C              of the piecewise polynomial representation.  (Input)
C     X      - The point whose location in BREAK is to be found.
C     LEFT   - Integer whose value is
C                LEFT
C                  1      IF                       X .LT.  BREAK(1)
C                  I      IF         BREAK(I).LE. X .LT. BREAK(I+1)
C                 NINTV   IF                    BREAK(NINTV) .LE. X
C              The asymmetric treatment of the interval is due to the
C              decision to make all PP functions continuous from the
C              right.  (Output)
C
C
C-----------------------------------------------------------------------
C
      SUBROUTINE P3DER (KORD, NINTV, BREAK, X, LEFT)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER    KORD, NINTV, LEFT
      REAL       X, BREAK(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    IHI, ISTEP, MIDDLE
C                                  SPECIFICATIONS FOR SAVE VARIABLES
      INTEGER    ILO
      SAVE       ILO
C
      DATA ILO/1/
C
      IHI = ILO + 1
      IF (IHI .GE. NINTV) THEN
         IF (X .GE. BREAK(NINTV)) THEN
            LEFT = NINTV
            GO TO 9000
         ELSE IF (NINTV .LE. 1) THEN
            LEFT = 1
            GO TO 9000
         END IF
         ILO = NINTV - 1
         IHI = NINTV
      END IF
C
      IF (X .LT. BREAK(IHI)) THEN
         IF (X .GE. BREAK(ILO)) THEN
            LEFT = ILO
            GO TO 9000
         END IF
C                                  Now X .LT. BREAK(ILO) - decrease ILO
C                                  to capture X
         ISTEP = 1
   10    CONTINUE
         IHI = ILO
         ILO = IHI - ISTEP
         IF (ILO .GT. 1) THEN
            IF (X .GE. BREAK(ILO)) GO TO 30
            ISTEP = ISTEP*2
            GO TO 10
         END IF
         ILO = 1
         IF (X .LT. BREAK(1)) THEN
            LEFT = 1
            GO TO 9000
         END IF
         GO TO 30
      END IF
C                                  Now X .GE. BREAK(IHI) - increase IHI
C                                  to capture X
      ISTEP = 1
   20 CONTINUE
      ILO = IHI
      IHI = ILO + ISTEP
      IF (IHI .LT. NINTV) THEN
         IF (X .LT. BREAK(IHI)) GO TO 30
         ISTEP = ISTEP*2
         GO TO 20
      END IF
      IF (X .GE. BREAK(NINTV)) THEN
         LEFT = NINTV
         GO TO 9000
      END IF
      IHI = NINTV
C                                  Now BREAK(ILO) .LE. X .LT.
C                                  BREAK(IHI) - narrow the inteval
   30 CONTINUE
      MIDDLE = (ILO+IHI)/2
      IF (MIDDLE .EQ. ILO) THEN
         LEFT = ILO
         GO TO 9000
      END IF
C                                  It is assumed that MIDDLE = ILO in
C                                  case IHI = ILO+1
      IF (X .LT. BREAK(MIDDLE)) THEN
         IHI = MIDDLE
      ELSE
         ILO = MIDDLE
      END IF
      GO TO 30
 9000 CONTINUE
      RETURN
      END
      
       subroutine mergrid(south,north,dlon,lata,dlata,n)

      implicit none
      integer n
      real*8 south,north,dlat,dlon,lat,lat0,dtor,lata(*),dlata(*)
      data dtor/1.7453292519943D-02/
      n=1
      dlat=dlon*dcos(south*dtor)
      dlata(1)=dlat*dtor
      lat0=south
      lata(1)=south*dtor
      lat=south
      do while (lat.le.north) 
      lat=lat0+dlat
      n=n+1
      lata(n)=lat*dtor
c      write(60,'(2f10.6)')lat,dlat
      lat0=lat
      dlat=dlon*dcos(lat0*dtor)
      dlata(n)=dlat*dtor
      end do
      north=lata(n-1)/dtor
      n=n-1
      return
      end
c
      real function gama(glat)
 
*** compute normal gravity  (returned in milligals (not m/(s*s) !)
***                         (input in radians -- single pre.)
 
      implicit double precision(a-h,o-z)
      real*8 glat
 
 
      parameter(e2 =0.00669438002290d0)
      parameter(cay=0.001931851353d0  )
      parameter(ge =978032.67715d0    )
 
*** somigilana's formula
 
      slat=dsin(glat)
      slat2=slat*slat
      gamma =ge*((1.d0+cay*slat2)/dsqrt(1.d0-e2*slat2))
      gama=sngl(gamma)
 
      return
      end

      subroutine write_error
      write(0,*)' '
      write(0,*)'Usage: inzone -Nnorth.grd3 -Eeast.grd3 -Gout.grd3 [-O]'
      write(0,*)' '
      write(0,*)'-N: input file of north component of gradient'
      write(0,*)'-E: input file of east component of gradient'
      write(0,*)'-G: utput file'
      write(0,*)'-O: compute effect on geoid [default: on gravity 
     &anomaly]'
      return
      end
