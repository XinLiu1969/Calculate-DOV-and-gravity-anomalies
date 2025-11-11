       program gift
c
c Gravity-geoid Inversion by FFT (gift)
c
c Program to calculate gravity anomalies or geoidal heights from geoid gradients 
c (or dov) using an 1-D FFT method. 
c
c Usage: gift -Nnorth.grd3 -Eeast.grd3 -Gout_file.grd3 
c        [-D -M -O -V -Cradius -L -Sstdx/stdy -Fwidth]
c
c See subroutine write_error for more info.
c
c *** PROGRAMMING NOTES ***
c
c 1. Use single precison for FFT computation
c 2. Use complex algebra for the FFT to get a concise code
c     

      implicit none
      integer*4 idim1,idim2,type,i,j,jy,nx2,ix,iy,nx0,ny0,lent,irept
     &	,iper,jmin,jmax,nx,ny

c     parameter (idim1=1962,idim2=idim1/2)
      parameter (idim1=(60*34+1)*2,idim2=idim1/2)

      integer*4 iwk(6*idim1+150)

      real*8 lata(idim2),dlat(idim2),dlon,time0,time1,
     1 xmin,xmax,ymin,ymax,gridx,gridy,west1,east1,south1,north1,gdx,gdy
     2 ,dtor,pi,re,pi4,gnorm,du,dv,rmax, dr,lon,lat

      real*4 north(idim2,idim2),east(idim2,idim2),
     1 half_lon(0:idim1),sin_lon(0:idim1),
     2 cos_lat(idim2),sin_lat(idim2),wk(6*idim1+150),
     3 out(idim2,idim2),gama,coslat,radius,smax,rad_max,
     4 factor,s_north,ss_north,s_east,ss_east,width

      complex*8 cwk(idim1),c_hcos(idim1),c_hsin(idim1),csum(idim1)
     &	,cn(idim1), ce(idim1),cnorth(idim1,idim2), ceast(idim1,idim2)
     &	,cwork(idim1,idim1) 
!      logical ir_isnan
c variables for input arguments
      integer*4 ii,nargs,iargc,i1,i2
      logical lg1, lg2, lg3, lgrad,lwt,lmer,lverb,lodd,lfilter,lnoise
      character*80 cha,file_n, file_e,fileanom,tbuf,file
ccccccccccccccccccccccccccc
      equivalence (north(1,1),out(1,1)),(cwork(1,1),cnorth(1,1))
	write(0,*)'New version of GIFT'

c First executable
C get arguments
      nargs=iargc()
      lg1=.true.
      lg2=.true.
      lg3=.true.
      lgrad=.true.
      lmer=.false.
      lverb=.true.
      lfilter=.false.
      lnoise=.false.
      type=0
      radius=110.0 
      width=999999.0

      file_n="re_gradients_B6_L12_n.grd3"
      file_e="re_gradients_B6_L12_e.grd3" 
      fileanom="dg_nofilter_B6_L12.grd3"
      file="dg_B3_L1.txt"
      call second(time0)
C Define some constants
      re=6371790.0d0 ! mean earth radius in meter
      pi=datan(1.d0)*4.d0
      pi4=4.d0*pi
      dtor=pi/180.0
      open(3,file=file_n,form='unformatted')
      read(3)ix,iy,xmin,xmax,ymin,ymax,gridx,gridy
      nx0=ix
      ny0=iy
      if(nx0.gt.idim2) stop'column out of bound'
      if(ny0.gt.idim2) stop'row out of bound'
      close(3)
      lent=56+nx0*ny0*4+4 
      open(10,file=file_n,form='unformatted')
      open(11,file=file_e,form='unformatted')
      
      lent=56+nx0*ny0*4
      open(20,file=fileanom,form='unformatted')
C read north component
      read(10)ix,iy,xmin,xmax,ymin,ymax,gridx,gridy 
      read(10)((north(i,j),j=1,iy),i=1,ix)
c      read(10,rec=1),
c     & ((north(i,j),j=1,iy),i=1,ix),s_north
      s_north=5.22
      if(lnoise) then
       s_north=ss_north**2
       else
       s_north=s_north**2
       end if
c read east component
      s_east=5.07
      read(11)ix,iy,west1,east1,south1,north1,gdx,gdy 
      read(11)((east(i,j),j=1,iy),i=1,ix)
c      read(11,rec=1)ix,iy,west1,east1,south1,north1,gdx,gdy,
c     & ((east(i,j),j=1,iy),i=1,ix),s_east
      if(lnoise) then
      s_east=ss_east**2 
      else
      s_east=s_east**2
      end if
      if (ix.ne.nx0)stop'x dim not mathced'
      if (iy.ne.ny0)stop'y dim not mathced'
      if(west1.ne.xmin) stop'west not matched'
      if(east1.ne.xmax) stop'east not matched'
      if(south1.ne.ymin) stop'south not matched'
      if(north1.ne.ymax) stop'north not matched'
      if(gridx.ne.gdx)stop'grid size x not mathched'
      if(gridy.ne.gdy)stop'grid size y not mathched'

      lodd=.false.
      if(mod(ix,2).ne.0) then
      lodd=.true.
      nx0=ix-1 
      end if

      if(lmer) then
      call mergrid(ymin,ymax,gridx,lata,dlat,i)
      if(i.ne.ny0)stop'mercator grid is in trouble'
      else
      do j=1,ny0
      lata(j)=(ymin+dfloat(j-1)*gridy)*dtor
      dlat(j)=gridy*dtor
      end do
      end if
      dlon=gridx*dtor


      if(lfilter) then
      factor=dcos((ymin+ymax)*0.5*dtor)
      nx=2*nx0
      ny=2*ny0

      s_north=s_north*dfloat(nx*ny)
      s_east=s_east*dfloat(nx*ny)
C fundamental frequencies
      du=1.0/(nx*gridx*factor)
      dv=1.0/(ny*gridy)
      dr=min(du,dv)
      rmax=dsqrt( 1/(2*gridx*factor)**2+ 1/(2*gridy)**2)
      call filter(north,cwork,idim2,idim1, nx,ny,nx0,ny0, du,dv,
     & rmax, dr, s_north,width)
      call filter(east,cwork,idim2,idim1,nx,ny,nx0,ny0, du,dv,
     & rmax, dr, s_east,width)
      end if


C Get tables. Use 100% zero padding
      nx2 =2*nx0
      call table(half_lon,sin_lon,cos_lat,sin_lat,lata,gridx,nx2,ny0)
C Compute FT of data for  each parallel
      do j=1,ny0
c Multiply each parallel by cos(lat)*dlat*dlon (all in radian)
      coslat=dcos(lata(j))*dlat(j)*dlon 
      do i=1,nx0
      cwk(i)=cmplx(north(i,j),east(i,j))*cmplx(coslat,0.0)
      cwk(i+nx0)=0.0
      end do
      call fftcc2(cwk,cn,ce,nx2,wk,iwk)
      do i=1,nx2
      cnorth(i,j)=cn(i) !(F(yita))
      ceast(i,j)=ce(i)  ! F(kesi)
      end do
      end do
      if(type.eq.1) then
      factor=-re*1.e-6/pi4/nx2
      if(.not.lgrad)factor=-factor
      end if
C Compute geoid or gravity for all grids along a parallel, from
c south to north
      radius=radius/111.1949266*dtor
      rad_max=20000.0/111.1949266*dtor
      smax=sin(radius/2.0)**2
c accumulate the contributions from the involved parallels  
c If -V then report progress every 10% 
      irept=ny0/10.0
      iper=0
      do 1 j=1,ny0
      do i=1,nx2
      csum(i)=0.0
      end do

      if(radius.lt.rad_max) then
      call limit(j,jmin,jmax,lata,ny0,radius)
      else
      jmin=1
      jmax=ny0
      end if

      do 2 jy=jmin,jmax
      call kernel(c_hcos,c_hsin,half_lon,sin_lon,cos_lat,sin_lat,lata,
     & j,jy,nx2,type,wk,iwk,cwk,smax)
      do i=1,nx2
      csum(i)=csum(i)+c_hcos(i)*cnorth(i,jy)+c_hsin(i)*ceast(i,jy)
      end do
2     continue
c get ready for inverse fft by fftcc
      do i=1,nx2
      csum(i)=conjg(csum(i))
      end do
      call fftcc(csum,nx2,iwk,wk)
C factor includes sign for dov or gradient, normal g, and nx2
      if(type.eq.0) then
      gnorm= gama(lata(j))
      factor = -gnorm*1.e-6/pi4/nx2
      if(.not.lgrad)factor=-factor
      end if
      do i=1,nx0
      out(i,j)=real(csum(i))*factor
!!	if(ir_isnan(out(i,j)))stop'NaN found in array, no output!!'
      end do
C If ix is odd then linearly extrapolate to get the last column
      if(lodd)out(ix,j)=2.0*out(nx0,j)-out(nx0-1,j)
      if(lverb) then
      if(mod(j,irept).eq.0) then
      iper=iper+10
      write(0,"('Percentage completed:', i3)")iper
      end if
      end if
1      continue
c output gravity anomalies or geoidal heights
 
      write(20)ix,iy,xmin,xmax,ymin,ymax,gridx,gridy
      write(20)((out(i,j),j=1,iy),i=1,ix)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*)"write .txt file....." 
      open(62,file=trim(file))
      do i=1,ix
        lon=xmin+DFLOAT(i-1)*gridx
         do j=1,iy
            lat=ymin+DFLOAT(j-1)*gridy
             write(62,100) lon,lat,out(i,j)
         enddo
      enddo
100       format(3f12.4)
      close(62)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call second(time1)
      WRITE(6,103)time1-time0
 103  FORMAT('total execution time (seconds):',F20.3)
      pause
      END
      subroutine second(time)
      real*8 TIME
      real*4 T(2)
      data TOT/0.D0/
      TIME=DTIME(T)+TOT
      TOT=TIME
      RETURN
      END

      subroutine limit(j,jmin,jmax,lata,ny0,radius)
c Program to find the upper and lower parallels for integration 
      real*8 lata(*)
      real*4 lat,lo,hi,radius
      integer j,jmin,jmax,ny0
      lat=lata(j)
      lo=lat-radius
      if(lo.lt.lata(1))lo=lata(1)
      hi=lat+radius
      if(hi.gt.lata(ny0))hi=lata(ny0)
      do n=j,1,-1
      if(lata(n).lt.lo) go to 1
      end do
1     jmin=n
      do n=j,ny0
      if(lata(n).gt.hi) go to 2
      end do
2     jmax=n
      if(jmin.lt.1)jmin=1
      if(jmax.gt.ny0)jmax=ny0
      return
      end

      subroutine kernel(c_hcos,c_hsin,half_lon,sin_lon,cos_lat
     &	,sin_lat,lata,p,q,n,type,wk,iwk,cwk,smax)
c Program to compute the kernel functions for geoid-gradient to gravity
c anomaly or geoid-gradient to geoid conversion
c
c c_hcos ... FT of dH/dpsi*cos(alpha) or dC/dpsi*cos(alpha) 
c c_hsin ... FT of dH/dpsi*sin(alpha)  or dC/dpsi*sin(alpha)
c half_lon ... 1-d array containing sin(dlon*i/2)**2, for i=0, ...
c cos_lat ... cos(lat(i)), for i=1, ...
c sin_lat ... sin(lat(i)), for i=1, ...
c lata ... 1-d array of latitude (in radian) from south to north
c p, q ... kernel function values are between latitude bands p and q.
c n ... number of points along a parallel, n must be even
c type ... 0 for converting geoid gradient to gravity anomaly, 1 for 
c          converting geoid gradient to geoid.
      implicit none
      integer type,p,q,n,iwk(*),nlo,nlo1,nlo2,iiu,iu2,iu
      real*4 half_lon(0:n),sin_lon(0:n),cos_lat(*),sin_lat(*),
     1 cder,hder,s,s2,c,sc,sin_ph2,sin_ph,a1,a2,alpha,cosa,sina,
     2 tmp,wk(*),smax
      real*8 lata(*)
      complex*8 cwk(n),c_hcos(n),c_hsin(n)
C psi-derivative of H-function (geoid gradient to gravity anomaly)
      hder(s,s2)=0.5*sqrt(1.0-s2)/s*(2.0*s2+2.0*s-1.0)/(s+s2)
c psi-derivative of C-function (geoid gradient to geoid)
c c=cos(psi), sc=sin(psi)
      cder(c,sc)=sc*(1.0/(c-1.0)+1.0+1.5*c)

C First executable
      nlo=n+1
      nlo1=n
      nlo2=n/2
*** note: wavenumber progression as i goes from 1 -> n
***  i    1, 2, 3, ...,  n/2   , (n/2)+1,  (n/2)+2, ..., n-2, n-1,  n
*** iiu   0, 1, 2, ..., (n/2)-1,  -n/2  , -(n/2)+1, ..., -3 , -2 , -1
***       ^-------------------^  ^----------------------------------^
***              positive                      negative
      sin_ph2=(dsin( (lata(p)-lata(q))/2.0)  )**2
      sin_ph=-dsin((lata(q)-lata(p)))
      do 1 iu=1,(nlo2+1)
        iiu=iu-1
        iu2=(nlo+1)-iu
        if(iu2.gt.nlo1) iu2=1


          s2= sin_ph2+ half_lon(iiu)*cos_lat(p)*cos_lat(q)
			        if(s2.gt.0.0 .and. s2.lt.smax) then
          s =sqrt(s2)
c Compute azimuth (alpha) from q to p (not from p to q !!)
         a1=cos_lat(p)*sin_lon(iiu)
         a2=sin_ph+2.0*sin_lat(q)*cos_lat(p)*half_lon(iiu)
         alpha=atan2(a1,a2)
         cosa=cos(alpha)
         sina=sin(alpha)
         if(type.eq.0) then
         tmp=hder(s,s2)
         else
          c=1.0-2.0*s2 
          sc=sqrt(1.0-c**2)
         tmp=cder(c,sc)
         end if
         cwk(iu)=cmplx(tmp*cosa,tmp*sina)
c The azimuth in the second half is negative to the azimuth
c in the first half, so we use conjugate 
         cwk(iu2)=conjg(cwk(iu))

			        else
            cwk(iu)=cmplx(0.0,0.0)
            cwk(iu2)=cmplx(0.0,0.0)
			        end if
   1  continue
c Forward Fourier transform
      call fftcc2(cwk,c_hcos,c_hsin,n,wk,iwk)
      return
      end
! H kernel
	function hder(s,s2)
	implicit none
	real*4 s,s2,hder
C psi-derivative of H-function (geoid gradient to gravity anomaly)
      if(s.eq.0) then
	hder=0
	return
        end if
      hder=0.5*sqrt(1.0-s2)/s*(2.0*s2+2.0*s-1.0)/(s+s2)
      return
	end
! 
	function cder(c,sc)
	implicit none
	real*4 c,sc,cder
c psi-derivative of C-function (geoid gradient to geoid)
c c=cos(psi), sc=sin(psi)
      if(c.eq.1.0) then
	cder=0
	return
	end if
      cder=sc*(1.0/(c-1.0)+1.0+1.5*c)
      return
	end


      subroutine table(half_lon,sin_lon,cos_lat,sin_lat,lata,dlon,nx,ny)
C Compute table of sin(dlon*i/2)**2, and sin (dlon*i), etc.
      implicit none
      integer nx,ny,i
      real*4 half_lon(0:nx),cos_lat(ny),sin_lat(ny),sin_lon(0:nx)
     1 ,lat,s,del
      real*8 dlon,lata(*),dtor
      data dtor/1.7453292519943D-02/

      del=0.5d0*dlon*dtor
      do 12 i=1,nx
      s=sin(i*del)
      half_lon(i)=s*s
      sin_lon(i)=sin(i*2.0*del)
   12 continue
      half_lon(0)=0.0
      sin_lon(0)=0.0

c Compute tables of cosine and sine latitude

      do 14 i=1,ny
        lat=lata(i)
        cos_lat(i)=cos(lat)
        sin_lat(i)=sin(lat)
   14 continue
      return
      end
      
       subroutine mergrid(south,north,dlon,lata,dlata,n)
c *** NOTE ***
c . north may be changed, south is unchanged.
c . the latitude in array lata(.) always increases as the index increases.

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
      real function gama(glat)
 
*** compute normal gravity  (returned in milligals (not m/(s*s) !)
***                         (input in radians -- single pre.)
 
      implicit double precision(a-h,o-z)
      real*8 glat
 
 
      parameter(e2 =0.00669438002290d0)
      parameter(cay=0.001931851353d0  )
      parameter(ge =978032.67715d0    )
 
*** somigilana's formula
 
      slat=sin(glat)
      slat2=slat*slat
      gamma =ge*((1.d0+cay*slat2)/dsqrt(1.d0-e2*slat2))
      gama=sngl(gamma)
 
      return
      end

      SUBROUTINE FFTCC2(ARR,A1,A2,NLON,WK,IWK)
C Compute ft of two real-valued functions at a time
c If   arr=cmplx(f1,f2), then
c A1 = FT of f1
c A2 = FT of f2
      REAL*4 WK(6*NLON+150),RN,RN1,AIN,AIN1
      INTEGER IWK(6*NLON+150)
      COMPLEX*8 ARR(NLON),A1(NLON),A2(NLON)
      CALL FFTCC(ARR,NLON,IWK,WK)
      A1(1)=REAL(ARR(1))
      A2(1)=AIMAG(ARR(1))
      N2=NLON+2 
      DO 2 N=2,NLON
      N1=N2-N
      RN=REAL(ARR(N))
      RN1=REAL(ARR(N1))
      AIN=AIMAG(ARR(N))
      AIN1=AIMAG(ARR(N1))
      A1(N)=CMPLX(0.5*(RN+RN1),0.5*(AIN-AIN1))
      A2(N)=CMPLX(0.5*(AIN+AIN1),0.5*(RN1-RN))
 2    CONTINUE
      RETURN
      END 
C   IMSL ROUTINE NAME   - FFTCC
C
C-----------------------------------------------------------------------

C
C   PURPOSE             - COMPUTE THE FAST FOURIER TRANSFORM OF A
C                           COMPLEX VALUED SEQUENCE
C
C   USAGE               - CALL FFTCC (A,N,IWK,WK)
C
C   ARGUMENTS    A      - COMPLEX VECTOR OF LENGTH N. ON INPUT A
C                           CONTAINS THE COMPLEX VALUED SEQUENCE TO BE
C                           TRANSFORMED. ON OUTPUT A IS REPLACED BY THE
C                           FOURIER TRANSFORM.
C                N      - INPUT NUMBER OF DATA POINTS TO BE
C                           TRANSFORMED. N MAY BE ANY POSITIVE
C                           INTEGER.
C                IWK    - INTEGER WORK VECTOR OF LENGTH 6*N+150.
C                           (SEE PROGRAMMING NOTES FOR FURTHER DETAILS)
C                WK     - REAL WORK VECTOR OF LENGTH 6*N+150.
C                           (SEE PROGRAMMING NOTES FOR FURTHER DETAILS)
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  FFTCC COMPUTES THE FOURIER TRANSFORM, X, ACCORDING
C                TO THE FOLLOWING FORMULA;
C
C                  X(K+1) = SUM FROM J = 0 TO N-1 OF
C                           A(J+1)*CEXP((0.0,(2.0*PI*J*K)/N))
C                  FOR K=0,1,...,N-1 AND PI=3.1415...
C
C                NOTE THAT X OVERWRITES A ON OUTPUT.
C            2.  FFTCC CAN BE USED TO COMPUTE
C
C                  X(K+1) = (1/N)*SUM FROM J = 0 TO N-1 OF
C                           A(J+1)*CEXP((0.0,(-2.0*PI*J*K)/N))
C                  FOR K=0,1,...,N-1 AND PI=3.1415...
C
C                BY PERFORMING THE FOLLOWING STEPS;
C
C                     DO 10 I=1,N
C                        A(I) = CONJG(A(I))
C                  10 CONTINUE
C                     CALL FFTCC (A,N,IWK,WK)
C                     DO 20 I=1,N
C                        A(I) = CONJG(A(I))/N
C                  20 CONTINUE

C
C-----------------------------------------------------------------------
C
      SUBROUTINE FFTCC (A,N,IWK,WK)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IWK(1)
      REAL*4		   WK(1)
      COMPLEX*8         A(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,IAM,IAP,IBM,IBP,IC,ICC,ICF,ICK,ID,IDM1,II,
     1                   IJA,IKB,IKT,ILL,IM,IRD,ISF,ISK,ISP,ISS,ITA,ITB,
     2                   J,JA,JF,JJ,JK,K,K0,K1,K2,K3,KA,KB,KD2,KF,KH,KN,
     3                   KT,KTP,L,L1,M,MM,MM1,MP
      real*4             CM,SM,C1,C2,C3,S1,S2,S3,C30,RAD,A0,A1,A4,B4,
     1                   A2,A3,B0,B1,B2,B3,ZERO,HALF,ONE,TWO,Z0(2),
     2                   Z1(2),Z2(2),Z3(2),Z4(2)
      COMPLEX*8          ZA0,ZA1,ZA2,ZA3,ZA4,AK2
      EQUIVALENCE        (ZA0,Z0(1)),(ZA1,Z1(1)),(ZA2,Z2(1)),
     1                   (ZA3,Z3(1)),(A0,Z0(1)),(B0,Z0(2)),(A1,Z1(1)),
     2                   (B1,Z1(2)),(A2,Z2(1)),(B2,Z2(2)),(A3,Z3(1)),
     3                   (B3,Z3(2)),(ZA4,Z4(1)),(Z4(1),A4),(Z4(2),B4)
      DATA               RAD/6.283185307179586D0/,
     1                   C30/.8660254037844386D0/
      DATA               ZERO,HALF,ONE,TWO/0.0D0,0.5D0,1.0D0,2.0D0/
C                                  FIRST EXECUTABLE STATEMENT
      IF (N .EQ. 1) GO TO 9005
      K = N
      M = 0
      J = 2
      JJ = 4
      JF = 0
C                                  DETERMINE THE SQUARE FACTORS OF N
      IWK(1) = 1
    5 I = K/JJ
      IF (I*JJ .NE. K) GO TO 10
      M = M+1
      IWK(M+1) = J
      K = I
      GO TO 5
   10 J = J + 2
      IF (J .EQ. 4) J = 3
      JJ = J * J
      IF (JJ .LE. K) GO TO 5
      KT = M
C                                  DETERMINE THE REMAINING FACTORS OF N
      J = 2
   15 I = K / J
      IF (I*J .NE. K) GO TO 20
      M = M + 1
      IWK(M+1) = J
      K = I
      GO TO 15
   20 J = J + 1
      IF (J .EQ. 3) GO TO 15
      J = J + 1
      IF(J.LE.K) GO TO 15
      K = IWK(M+1)
      IF (IWK(KT+1) .GT. IWK(M+1)) K = IWK(KT+1)
      IF(KT.LE.0) GO TO 30
      KTP = KT + 2
      DO 25  I = 1,KT
         J = KTP - I
         M = M+1
         IWK(M+1) = IWK(J)
   25 CONTINUE
   30 MP = M+1
      IC = MP+1
      ID = IC+MP
      ILL = ID+MP
      IRD = ILL+MP+1
      ICC = IRD+MP
      ISS = ICC+MP
      ICK = ISS+MP
      ISK = ICK+K
      ICF = ISK+K
      ISF = ICF+K
      IAP = ISF+K
      KD2 = (K-1) / 2 + 1
      IBP = IAP + KD2
      IAM = IBP + KD2
      IBM = IAM + KD2
      MM1 = M-1
      I = 1
   35 L = MP - I
      J = IC - I
      IWK(ILL+L) = 0
      IF ((IWK(J-1) + IWK(J)) .EQ. 4) IWK(ILL+L) = 1
      IF (IWK(ILL+L) .EQ. 0) GO TO 40
      I = I + 1
      L = L - 1
      IWK(ILL+L) = 0
   40 I = I + 1
      IF(I.LE.MM1) GO TO 35
      IWK(ILL+1) = 0
      IWK(ILL+MP) = 0
      IWK(IC) = 1
      IWK(ID) = N
      DO 45  J = 1,M
         K = IWK(J+1)
         IWK(IC+J) = IWK(IC+J-1) * K
         IWK(ID+J) = IWK(ID+J-1) / K
         WK(IRD+J) = RAD/IWK(IC+J)
         C1 = RAD/K
         IF (K .LE. 2) GO TO 45
         WK(ICC+J) = COS(C1)
         WK(ISS+J) = SIN(C1)
   45 CONTINUE
      MM = M
      IF (IWK(ILL+M) .EQ. 1) MM = M - 1
      IF (MM .LE. 1) GO TO 50
      SM = IWK(IC+MM-2) * WK(IRD+M)
      CM = COS(SM)
      SM = SIN(SM)
   50 KB = 0
      KN = N
      JJ = 0
      I = 1
      C1 = ONE
      S1 = ZERO
      L1 = 1
   55 IF (IWK(ILL+I+1) .EQ. 1) GO TO 60
      KF = IWK(I+1)
      GO TO 65
   60 KF = 4
      I = I+1

   65 ISP = IWK(ID+I)
      IF (L1 .EQ. 1) GO TO 70
      S1 = JJ * WK(IRD+I)
      C1 = COS(S1)
      S1 = SIN(S1)
C                                  FACTORS OF 2, 3, AND 4 ARE
C                                  HANDLED SEPARATELY.
   70 IF (KF .GT. 4) GO TO 140
      GO TO (75,75,90,115), KF
   75 K0 = KB + ISP
      K2 = K0 + ISP
      IF (L1 .EQ. 1) GO TO 85
   80 K0 = K0 - 1
      IF (K0 .LT. KB) GO TO 190
      K2 = K2 - 1
      ZA4 = A(K2+1)
      A0 = A4*C1-B4*S1
      B0 = A4*S1+B4*C1
      A(K2+1) = A(K0+1)-ZA0
      A(K0+1) = A(K0+1)+ZA0
      GO TO 80
   85 K0 = K0 - 1
      IF (K0 .LT. KB) GO TO 190
      K2 = K2 - 1
      AK2 = A(K2+1)
      A(K2+1) = A(K0+1)-AK2
      A(K0+1) = A(K0+1)+AK2
      GO TO 85
   90 IF (L1 .EQ. 1) GO TO 95
      C2 = C1 * C1 - S1 * S1
      S2 = TWO * C1 * S1
   95 JA = KB + ISP - 1
      KA = JA + KB
      IKB = KB+1
      IJA = JA+1
      DO 110 II = IKB,IJA
         K0 = KA - II + 1
         K1 = K0 + ISP
         K2 = K1 + ISP
         ZA0 = A(K0+1)
         IF (L1 .EQ. 1) GO TO 100
         ZA4 = A(K1+1)
         A1 = A4*C1-B4*S1
         B1 = A4*S1+B4*C1
         ZA4 = A(K2+1)
         A2 = A4*C2-B4*S2
         B2 = A4*S2+B4*C2
         GO TO 105
  100    ZA1 = A(K1+1)
         ZA2 = A(K2+1)
  105    A(K0+1) = CMPLX(A0+A1+A2,B0+B1+B2)
         A0 = -HALF * (A1+A2) + A0
         A1 = (A1-A2) * C30
         B0 = -HALF * (B1+B2) + B0
         B1 = (B1-B2) * C30
         A(K1+1) = CMPLX(A0-B1,B0+A1)
         A(K2+1) = CMPLX(A0+B1,B0-A1)
  110 CONTINUE
      GO TO 190
  115 IF (L1 .EQ. 1) GO TO 120
      C2 = C1 * C1 - S1 * S1
      S2 = TWO * C1 * S1
      C3 = C1 * C2 - S1 * S2
      S3 = S1 * C2 + C1 * S2
  120 JA = KB + ISP - 1
      KA = JA + KB
      IKB = KB+1
      IJA = JA+1
      DO 135 II = IKB,IJA
         K0 = KA - II + 1
         K1 = K0 + ISP
         K2 = K1 + ISP
         K3 = K2 + ISP
         ZA0 = A(K0+1)
         IF (L1 .EQ. 1) GO TO 125
         ZA4 = A(K1+1)
         A1 = A4*C1-B4*S1
         B1 = A4*S1+B4*C1
         ZA4 = A(K2+1)
         A2 = A4*C2-B4*S2
         B2 = A4*S2+B4*C2
         ZA4 = A(K3+1)
         A3 = A4*C3-B4*S3
         B3 = A4*S3+B4*C3
         GO TO 130
  125    ZA1 = A(K1+1)
         ZA2 = A(K2+1)
         ZA3 = A(K3+1)
  130    A(K0+1) = CMPLX(A0+A2+A1+A3,B0+B2+B1+B3)
         A(K1+1) = CMPLX(A0+A2-A1-A3,B0+B2-B1-B3)
         A(K2+1) = CMPLX(A0-A2-B1+B3,B0-B2+A1-A3)
         A(K3+1) = CMPLX(A0-A2+B1-B3,B0-B2-A1+A3)
  135 CONTINUE
      GO TO 190
  140 JK = KF - 1
      KH = JK/2
      K3 = IWK(ID+I-1)
      K0 = KB + ISP
      IF (L1 .EQ. 1) GO TO 150
      K = JK - 1
      WK(ICF+1) = C1
      WK(ISF+1) = S1
      DO 145 J = 1,K
         WK(ICF+J+1) = WK(ICF+J) * C1 - WK(ISF+J) * S1
         WK(ISF+J+1) = WK(ICF+J) * S1 + WK(ISF+J) * C1
  145 CONTINUE
  150 IF (KF .EQ. JF) GO TO 160
      C2 = WK(ICC+I)
      WK(ICK+1) = C2
      WK(ICK+JK) = C2
      S2 = WK(ISS+I)
      WK(ISK+1) = S2
      WK(ISK+JK) = -S2
      DO 155 J = 1,KH
         K = JK - J
         WK(ICK+K) = WK(ICK+J) * C2 - WK(ISK+J) * S2
         WK(ICK+J+1) = WK(ICK+K)
         WK(ISK+J+1) = WK(ICK+J) * S2 + WK(ISK+J) * C2
         WK(ISK+K) = -WK(ISK+J+1)
  155 CONTINUE
  160 K0 = K0 - 1
      K1 = K0
      K2 = K0 + K3
      ZA0 = A(K0+1)
      A3 = A0
      B3 = B0
      DO 175 J = 1,KH
         K1 = K1 + ISP
         K2 = K2 - ISP
         IF (L1 .EQ. 1) GO TO 165
         K = KF - J
         ZA4 = A(K1+1)
         A1 = A4*WK(ICF+J)-B4*WK(ISF+J)
         B1 = A4*WK(ISF+J)+B4*WK(ICF+J)
         ZA4 = A(K2+1)
         A2 = A4*WK(ICF+K)-B4*WK(ISF+K)
         B2 = A4*WK(ISF+K)+B4*WK(ICF+K)
         GO TO 170
  165    ZA1 = A(K1+1)
         ZA2 = A(K2+1)
  170    WK(IAP+J) = A1 + A2
         WK(IAM+J) = A1 - A2
         WK(IBP+J) = B1 + B2
         WK(IBM+J) = B1 - B2
         A3 = A1 + A2 + A3
         B3 = B1 + B2 + B3
  175 CONTINUE
      A(K0+1) = CMPLX(A3,B3)
      K1 = K0
      K2 = K0 + K3
      DO 185 J = 1,KH
         K1 = K1 + ISP
         K2 = K2 - ISP
         JK = J
         A1 = A0
         B1 = B0
         A2 = ZERO
         B2 = ZERO
         DO 180  K = 1,KH
            A1 = A1 + WK(IAP+K) * WK(ICK+JK)
            A2 = A2 + WK(IAM+K) * WK(ISK+JK)
            B1 = B1 + WK(IBP+K) * WK(ICK+JK)
            B2 = B2 + WK(IBM+K) * WK(ISK+JK)
            JK = JK + J
            IF (JK .GE. KF) JK = JK - KF
  180    CONTINUE
         A(K1+1) = CMPLX(A1-B2,B1+A2)
         A(K2+1) = CMPLX(A1+B2,B1-A2)
  185 CONTINUE
      IF (K0 .GT. KB) GO TO 160
      JF = KF
  190 IF ( I .GE. MM ) GO TO 195
      I = I + 1
      GO TO 55
  195 I = MM
      L1 = 0
      KB = IWK(ID+I-1) + KB
      IF (KB .GE. KN) GO TO 215
  200 JJ = IWK(IC+I-2) + JJ
      IF (JJ .LT. IWK(IC+I-1)) GO TO 205
      I = I - 1
      JJ = JJ - IWK(IC+I)
      GO TO 200
  205 IF (I .NE. MM) GO TO 210
      C2 = C1
      C1 = CM * C1 - SM * S1
      S1 = SM * C2 + CM * S1
      GO TO 70
  210 IF (IWK(ILL+I) .EQ. 1) I = I + 1
      GO TO 55
  215 I = 1
      JA = KT - 1
      KA = JA + 1
      IF(JA.LT.1) GO TO 225
      DO 220  II = 1,JA
         J = KA - II
         IWK(J+1) = IWK(J+1) - 1
         I = IWK(J+1) + I
  220 CONTINUE
C                                  THE RESULT IS NOW PERMUTED TO
C                                  NORMAL ORDER.
  225 IF (KT .LE. 0) GO TO 270
      J = 1
      I = 0
      KB = 0
  230 K2 = IWK(ID+J) + KB
      K3 = K2
      JJ = IWK(IC+J-1)
      JK = JJ
      K0 = KB + JJ
      ISP = IWK(IC+J) - JJ
  235 K = K0 + JJ
  240 ZA4 = A(K0+1)
      A(K0+1) = A(K2+1)
      A(K2+1) = ZA4
      K0 = K0 + 1
      K2 = K2 + 1
      IF (K0 .LT. K) GO TO 240
      K0 = K0 + ISP
      K2 = K2 + ISP
      IF (K0 .LT. K3) GO TO 235
      IF (K0 .GE. K3 + ISP) GO TO 245
      K0 = K0 - IWK(ID+J) + JJ
      GO TO 235
  245 K3 = IWK(ID+J) + K3
      IF (K3 - KB .GE. IWK(ID+J-1)) GO TO 250
      K2 = K3 + JK
      JK = JK + JJ
      K0 = K3 - IWK(ID+J) + JK
      GO TO 235
  250 IF (J .GE. KT) GO TO 260
      K = IWK(J+1) + I
      J = J + 1
  255 I = I + 1
      IWK(ILL+I) = J
      IF (I .LT. K) GO TO 255
      GO TO 230
  260 KB = K3
      IF (I .LE. 0) GO TO 265
      J = IWK(ILL+I)
      I = I - 1
      GO TO 230
  265 IF (KB .GE. N) GO TO 270
      J = 1
      GO TO 230
  270 JK = IWK(IC+KT)
      ISP = IWK(ID+KT)
      M = M - KT
      KB = ISP/JK-2
      IF (KT .GE. M-1 ) GO TO 9005
      ITA = ILL+KB+1
      ITB = ITA+JK
      IDM1 = ID-1
      IKT = KT+1
      IM = M+1
      DO 275 J = IKT,IM
         IWK(IDM1+J) = IWK(IDM1+J)/JK
  275 CONTINUE
      JJ = 0
      DO 290 J = 1,KB
         K = KT
  280    JJ = IWK(ID+K+1) + JJ
         IF (JJ .LT. IWK(ID+K)) GO TO 285
         JJ = JJ - IWK(ID+K)
         K = K + 1
         GO TO 280
  285    IWK(ILL+J) = JJ
         IF (JJ .EQ. J) IWK(ILL+J) = -J
  290 CONTINUE
C                                  DETERMINE THE PERMUTATION CYCLES
C                                  OF LENGTH GREATER THAN OR EQUAL
C                                  TO TWO.
      DO 300  J = 1,KB
         IF (IWK(ILL+J) .LE. 0) GO TO 300
         K2 = J
  295    K2 = IABS(IWK(ILL+K2))
         IF (K2 .EQ. J) GO TO 300
         IWK(ILL+K2) = -IWK(ILL+K2)
         GO TO 295
  300 CONTINUE
C                                  REORDER A FOLLOWING THE
C                                  PERMUTATION CYCLES
      I = 0
      J = 0
      KB = 0
      KN = N
  305 J = J + 1
      IF (IWK(ILL+J) .LT. 0) GO TO 305
      K = IWK(ILL+J)
      K0 = JK * K + KB
  310 ZA4 = A(K0+I+1)
      WK(ITA+I) = A4
      WK(ITB+I) = B4
      I = I + 1
      IF (I .LT. JK) GO TO 310
      I = 0
  315 K = -IWK(ILL+K)
      JJ = K0
      K0 = JK * K + KB
  320 A(JJ+I+1) = A(K0+I+1)
      I = I + 1
      IF (I .LT. JK) GO TO 320
      I = 0
      IF (K .NE. J) GO TO 315
  325 A(K0+I+1) = CMPLX(WK(ITA+I),WK(ITB+I))
      I = I + 1
      IF (I .LT. JK) GO TO 325
      I = 0
      IF (J .LT. K2) GO TO 305
      J = 0
      KB = KB + ISP
      IF (KB .LT. KN) GO TO 305
 9005 RETURN
      END

      SUBROUTINE FFT2D  (A,IA1,N1,N2,IJOB,IWK,RWK,CWK)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IA1,N1,N2,IJOB,IWK(1)
      real*4   RWK(1)
      COMPLEX*8         A(IA1,N2),CWK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,J,L,M
      real*4   R12
      COMPLEX*8         C12
C                                  FIRST EXECUTABLE STATEMENT
      IF (IJOB.GT.0) GO TO 10
C                                  INVERSE TRANSFORM
      DO 5 I=1,N1
      DO 5 J=1,N2
         A(I,J) = CONJG(A(I,J))
    5 CONTINUE

C                                  TRANSFORM SECOND SUBSCRIPT
10     DO 40 L=1,N1
         DO 30 M=1,N2
            CWK(M) = A(L,M)
   30    CONTINUE
         CALL FFTCC (CWK,N2,IWK,RWK)
         DO 35 J=1,N2
            A(L,J) = CWK(J)
   35    CONTINUE
   40 CONTINUE
C                                  TRANSFORM FIRST SUBSCRIPT
      DO 45 J=1,N2
         CALL FFTCC (A(1,J),N1,IWK,RWK)
   45 CONTINUE
      IF (IJOB.GT.0) GO TO 55
      R12 = N1*N2
      C12 = CMPLX(R12,0.0D0)
      DO 50 I=1,N1
      DO 50 J=1,N2
         A(I,J) = CONJG(A(I,J))/C12
   50 CONTINUE
   55 RETURN
      END

      subroutine filter(arr,carra,idim,idim1, nx,ny,nyqx,nyqy, du,dv,
     & rmax, dr, noise,width)
C Subroutine to estimate the Wiener filter given a white noise spectrum.
C The filter is fitted by the Gaussian: exp(-a* f**2)
c 
c
      implicit none
      integer*4 idimx,idim1
      parameter(idimx=4000)
      integer*4 is(3000),iwk(6*idimx+150),
     & nx,ny,nyqx,nyqy, idim,iter,
     & i, j, i1, j1,k,icount,n
      real*8  pw(3000), saa, tmp, tem1, a,du,dv,rmax,dr,r,u,v,da,
     & sum, sum1,f2
      real*4 noise,work(3000),aa,arr(idim,nyqy),rwk(6*idimx+150),q2
     &	,width

      complex*8 carra(idim1,*),cwk(idimx)

      if(idim1.gt.idimx) stop
      do j=1,ny
      do i=1,nx
      carra(i,j)=0
      end do
      end do
      do j=1,nyqy
      do i=1,nyqx
      carra(i,j)=arr(i,j)
      end do
      end do
      call fft2d (carra,idim1,nx,ny,1,iwk,rwk,cwk)

      if(width.eq.999999.0) then

C  
      k= int(rmax/dr)+1
      if(k.gt.3000)stop'filter: k greater than 3000' ! k取值依据
c      write(6,*)rmax, dr, dr*k
      do i=1,k
      pw(i)=0.d0
      is(i)=0
      end do

      icount=0
      do 31 j = 0, ny-1
      j1=j+1

      if(j.le.nyqy) then
      v=dfloat(j)*dv
      else
      v=dfloat(j-ny)*dv
      end if

      do 31 i = 0, nx-1
      i1=i+1
      if(i.le.nyqx)then
      u=dfloat(i)*du
      else
      u=dfloat(i-nx)*du
      end if

      r=dsqrt(u**2 + v**2 )
      n=nint(r/dr)+1
      if(n.gt.k) go to 31

      saa=real( conjg(carra(i1,j1)) * carra(i1,j1) )
      if(saa.gt.0.0) then
      tem1=saa-noise
    	  if(tem1.gt.0) then
	      tem1=tem1/saa
              pw(n)=pw(n)+tem1
              is(n)=is(n)+1 
 	     end if
      end if


   31 continue

C Average the filter
      do i=1,k
      if(is(i).ne.0) then
      pw(i) = pw(i)/ is(i)
c      write(iu,'(2f7.4)')(i-1)*dr, pw(i)
      else
      pw(i)=9999.0
      end if
      end do

C Fit the filter by a Gaussian having the form exp(- a f**2)
      do i=1,k
      work(i)=((i-1)*dr)**2
      end do
      a=0

10    do iter=1,10
      sum=0.d0
      sum1=0.0
       do i=1,k

      if(pw(i).ne.9999.0) then
      f2=work(i)
      tmp=dexp(-a*f2)
      sum=sum+ f2*f2 * tmp*tmp
      sum1=sum1-f2*tmp*(pw(i)-tmp)
      end if
      end do

      da=sum1/sum
c      write(6,*)'da,i1=', da,i1
      a=da+a
      if(dabs(da).lt.1.d-6) go to 11 
      end do
11    continue
      write(6,100)dsqrt(a/3.141592654)*111.1949
100   format('Filter width of the Gaussian filter:', f7.2, ' km')
c      do i=1,k
c      write(iu1,'(2f7.4)')(i-1)*dr, dexp(-a*work(i))
c      end do
      aa=a
      else
      aa= (width/111.1949)**2*3.141592654
      end if

      do 32 j=0,ny-1
 
      j1=j+1
      if(j.le.nyqy) then
      v=dfloat(j)*dv
      else
      v=dfloat(j-ny)*dv
      end if

      do 32 i=0,nx-1
      i1=i+1
      if(i.le.nyqx)then
      u=dfloat(i)*du
      else
      u=dfloat(i-nx)*du
      end if
      q2=u*u+v*v
C Filter gradient only when the frequency is greater than 1 cycle/degree
      if(q2.gt.1)carra(i1,j1)=carra(i1,j1) * exp(-aa*q2)
32    continue
      call fft2d (carra,idim1,nx,ny,-1,iwk,rwk,cwk)
      do i=1,nyqx
      do j=1,nyqy
      arr(i,j)=carra(i,j)
      end do
      end do
      return
      end 

      subroutine write_error
      write(0,*)
      write(0,*)'Usage: gift -Nnorth -Eeast -Gout_file 
     &	[-D -W -M -C -O -V -Sstdx/stdy -L]'
      write(0,*)' '
      write(0,*)'-N	input file of north component of gradient (or dov)'
      write(0,*)'-E	input file of east component of gradient (or dov)'
      write(0,*)'-G	output file'
      write(0,*)'-D	input data type is deflection of the vertical
     &	[default: geoid gradient]'
      write(0,*)'-M	use Mercator grid [default: .grd1 grid]'
      write(0,*)'-C	influence radius in km [default: 20000 km]'
      write(0,*)'-O	compute geoid [default: gravity anomaly]'
      write(0,*)'-V	report progress every 10% [default: silent]'
      write(0,*)'-S	std of gradient in east and north directions 
     &	[default: from the gradient files]'
      write(0,*)'-L	use adaptive filters [default: no filter]'
      return
      end

