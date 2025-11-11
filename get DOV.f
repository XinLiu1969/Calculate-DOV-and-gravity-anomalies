! compute north and east gradients and errors on a grid 
!         using collocation                                                              !           
!
! DESCRIPTION                                                                                      
!  f.xyg: file containing geoid gradients (in arc-sec) in        
!                    ascii form  
! -C: covaraince table                               
! -I: grid interval in long and lat directions in       
!                    minutes                                           
! -N:ouput geoid gradient in north-south direction      
!                    in micro-radian                                      
! -E:ouput geoid gradient in west-east direction in     
!                    micro-radian
! -W: size of window (in minutes) for data selection
!
! OPTIONS         
! -S: Scale global covariance function [default: no]
! -O: Compute error estimates of gradients [default: no]
!
! NOTES
! (1) The input file (f.xyg) is created by rmgrad.f  or bingrad.f, and
!      each record contains
!
!     longitude  : degree					        
!     latitude   : degree                            
!     SSG        : arc-second					        
!     azimuth    : radian					        
!     STD.       : arc-second 
!                                         
! (2) The output files are                                                   
!     north.grd3 : north direction gradient of ssg and  std error 
!                 (unit : microrad)    
!     east.grd3   : east  direction gradient os ssg and std error 
!                  (unit : microrad)
!    north_err.grd3: error estimate of north gradient  
!    east_err.grd3 : error estimate of east gradient  
!                                                                     
!========1=========2=========3=========4=========5=========6=========7==
      implicit none
      integer*4  I, J,n, N2,IWIN, IWINX,npp,win0,sum_npt,
     &  NX,NY,NARGS,IARGC,ik1,ik2,ik3,ik4,iik,par,k,
     &  idim,nda,num,icov,ncov,npred,npt,i1,j1,i2,j2,i3,
     &  il,ir,jd,ju,nsel,ii1,m,nskip
      integer nsing
      parameter(idim=60*34+1,NDA=20,NUM=6000) 
      
      integer*2 count(idim,idim)
      integer*4 indx(num*(num+1)/2)

      real*8 west,east,south,north,dx,dy,fx,fy,window,time0,sum_fact,
     & PSI,PSIIV,VAR,dist,tmp1,tmp2,tmp3,zero,pi,redtor,gridx,gridy,
     & c0,c2,tmp,factor1,factor2,halfpi,cosfi,dtor,covee,xi,yi,gridxx,
     & distance,timex,azm,zpred1,zpred2,toldst

      real*8 C(NUM*(NUM+1)/2),XTMP(NUM),YTMP(NUM),ZTMP(NUM),
     &  AZTMP(NUM),STDTMP(NUM),CLL(1001),CMM(1001),
     & ARR1(NUM),ARR2(NUM),xARR1(NUM),xARR2(NUM),work(num),
     & dov_n,dov_e
    
      real*4 LAT,LON,DOV,AZ,STD,xi_err,eta_err

      real*4 DA(IDIM,IDIM,5,nda),ARRX(IDIM,IDIM),ARRY(IDIM,IDIM),
     & ARRX_err(IDIM,IDIM),ARRY_err(IDIM,IDIM),tmpa(nda),record(20,20),
     & result(4)

 
      character*80 filename,cha,ofile2,ofile1,file_i,tbuf,ifile1,
     & ofile3,ofile4

      logical io(10),lscale,lerr
      COMMON /COVTAB/CLL,CMM,PSI,PSIIV,VAR
      COMMON/GRDINFO/NX,NY,WEST,EAST,SOUTH,NORTH,FX,FY
      data icov/1001/
	data toldst/0.d0/ 
      lscale=.true.
      lerr=.true.
      filename="cov_xgm2019e(2159).txt"                    
      file_i="re_geiod_gradients_B6_L8.dat"             
      ofile1="re_gradients_B6_L8_n.grd3"
      ofile2="re_gradients_B6_L8_e.grd3"
      ofile3="re_gradients_B6_L8_n.txt"
      ofile4="re_gradients_B6_L8_e.txt"
      DX=1.0
      DY=1.0
      east=-119
      west=-151
      north=88
      south=59
      FX=DX/60.d0
      FY=DY/60.d0
      NX=(east-west)/fx+1.01
      NY=(north-south)/fy+1.01
      window=1 
	write(*,*)Nx,NY,idim
      if(NX.gt.idim .or. NY.gt.idim) stop'increase idim'      

      OPEN(60,FILE=ofile1,form='unformatted') ! north component
      OPEN(61,FILE=ofile2,form='unformatted') ! east component
      OPEN(62,FILE=ofile3) ! north component
      OPEN(63,FILE=ofile4) ! east component

      WRITE(*,*)'Reading covariance table ...'
      OPEN(100,FILE=filename,STATUS='OLD')
      READ(100,*)NCOV,PSI
      if(ncov.gt.icov) stop'increase icov'
      DO  I=1,NCOV
          READ(100,253)DIST,TMP1,TMP2,TMP3, CLL(I),CMM(I)
253       FORMAT(F6.2,5D20.13)
      END DO
      close(100)
! check table
      IF(abs(CLL(1)-CMM(1)) .gt. 1.d-10)  
     &STOP'COVARAINCE TABLE NOT CORRECT'
      VAR=CLL(1)
      PSIIV=1.D0/PSI
      
 !     WRITE(0,*)'Reading landmask ... '

 !     CALL READMASK
   
      ZERO=0.D0
      DTOR=0.01745329252D0

      IWIN= nint(window/dx/2.0)

      win0=window/dx/2.0/dcos(dmax1(abs(north),abs(south))*dtor)

      if((2*iwin+1)*(2*win0+1)*5.gt.num) stop'increase NUM'    

      PI=DATAN(1.D0)*4.D0
      HALFPI=PI/2.D0
      factor2=1.d0/0.2062648062d0
      REDTOR=111.1949266D0
	toldst=(toldst/6371.d0/dtor)**2
      WRITE(6,*)'Output GRD3 info:',NX,NY,west,east,south,north,fx,fy     
      OPEN(11,FILE=file_i)
   
      DO  I=1,idim
          DO  J=1,idim
	        count(i,j)=0
          ENDDO
      ENDDO 

      write(0,*)'Reading gradient data ...'      
      npt=0
      
102   read(11,*,end=103)lon,lat,dov, az, std
      
      IF(LAT.LT.south.OR.LAT.GT.north) GOTO 102
      if(lon.lt.west .or. lon.gt.east) go to 102
      I= (LON-WEST)/FX+1.01                       
      J= (LAT-south)/FY+ 1.01
      npt=npt+1
! count(i,j) records the number of gradients at grid (i,j)
      count(i,j)=count(i,j)+1
      if(count(i,j).gt.nda) then
          write(6,*)'Cell with clustered data:',lon,lat,i,j,count(i,j)
          count(i,j)=nda
          go to 102
!       stop'increase NDA'
      endif
      da(i,j,1,count(i,j))=lat
      da(i,j,2,count(i,j))=lon
      da(i,j,3,count(i,j))=dov
      da(i,j,4,count(i,j))=az
      da(i,j,5,count(i,j))=std**2
      go to 102
       
103   continue
      
      WRITE(6,*)'Number of selected gradients:',npt

! Select data if count(i,j) > 5
       do j=1,ny
           do i=1,nx
               if(count(i,j).gt.5) then
                   do k=1, count(i,j)
                      tmpa(k)=da(i,j,5,k)
                      indx(k)=k
                   end do
                   n=count(i,j)
                   call indexx(n,tmpa,indx)       
                   do k=1,5
                      record(k,1)=da(i,j,1,indx(k))
                      record(k,2)=da(i,j,2,indx(k))
                      record(k,3)=da(i,j,3,indx(k))
                      record(k,4)=da(i,j,4,indx(k))
                      record(k,5)=da(i,j,5,indx(k))
                   end do
                   count(i,j)=5
                   do k=1,count(i,j)
                      da(i,j,1,k)=record(k,1)
                      da(i,j,2,k)=record(k,2)
                      da(i,j,3,k)=record(k,3)
                      da(i,j,4,k)=record(k,4)
                      da(i,j,5,k)=record(k,5)
                   end do
               end if
           end do
       end do


       CALL SECOND(TIME0)

       write(0,*)'Start gridding ...'      
       npred=0
       npp=0
       sum_npt=0
       sum_fact=0
        DO J1=1,ny  
           GRIDY=south+DFLOAT(J1-1)*FY
           COSFI=DCOS(GRIDY*DTOR)
           IWINX=NINT(IWIN/COSFI)
          if(dmod(gridy,1.d0).eq.0.d0)
     &    write(*,*)'Gridding is now at latitude=',nint(gridy)
           JD=max0(J1-IWIN,1)  !lat
           JU=min0(J1+IWIN,ny)
        DO I1=1,NX !lon
           GRIDX=WEST+DFLOAT(I1-1)*FX
           GRIDXX=GRIDX*COSFI

           IL=max0(I1-IWINX,1) !lon
           IR=min0(I1+IWINX,NX)         
           N2=0           
           DO J2=JD,JU
               DO I2=IL,IR
                   DO I3=1,count(I2,J2)
!                      Search data in a (2*iwinx+1) by (2*iwin+1)  window
!                      dist=sqrt(((i1-i2)*cosfi)**2+(j1-j2)**2)
!                      if(dist.le.iwin) then
 			         N2=N2+1 
			         if(n2.gt.num) then
	                     n2=n2-1
			             go to 107
	                 end if
                 
                       YTMP(N2)= DA(I2,J2,1,I3)
                       XTMP(N2)= DA(I2,J2,2,I3)*COSFI
                       ZTMP(N2)=DA(I2,J2,3,I3)
                       AZTMP(N2)=DA(I2,J2,4,I3)
                       STDTMP(N2)=DA(I2,J2,5,I3)
!               end if
                   END DO 

105           END DO
106       END DO
107    continue

C Predictions are made only no. of points greater than 4

      if(n2.le.4) then
          ARRY(I1,J1)=999.0
          ARRX(I1,J1)=999.0
	    arrx_err(i1,J1)=999.0
	    arry_err(i1,J1)=999.0
          go to 110
      end if    
! Remove data points that are too close to avoid singulars in
! lsc 
      call  SLTPTV(Xtmp,Ytmp,N2,TOLDST,NSKIP,INDX)
      IF(NSKIP.GT.0) THEN
          m=0
          DO 6 I=1,N2
              DO 7 j=1,NSKIP
                  IF(I.EQ.INDX(j)) GO TO 6
7                 CONTINUE
                  m=m+1
                  xtmp(m)=xtmp(I)
                  ytmp(m)=ytmp(I)
                  ztmp(m)=ztmp(I)
	            AZTMP(m)=AZTMP(i)
                  stdtmp(m)=stdtmp(I)
6                 CONTINUE
                  N2=m
      END IF
C Predicting and filtering by collocation       
C Compute covariance matrices

      DO I=1,N2
C NORTH COMPONENT
          arr1(I)=COVEE(GRIDXX,GRIDY,XTMP(I),YTMP(I),ZERO,aztmp(i),
     &                  distance)
          xarr1(I)=arr1(i)
C EAST COMPONENT
C====&===1=========2=========3=========4=========5=========6=========7==
         ARR2(I)=COVEE(GRIDXX,GRIDY,XTMP(I),YTMP(I),HALFPI,aztmp(i)
     &           ,distance)
          xarr2(I)=ARR2(i)
      END DO 
      nsel=n2

C Scaling factor will affect lsc filtering.
      if (lscale) then
          C0=0.D0
          C2=0.D0
          DO I=1,NSEL
              TMP=ztmp(i)
              C0=C0+TMP
              C2=C2+TMP*TMP
          END DO
          C0=C0/NSEL
          C2=(C2/NSEL-C0*C0)*nsel/(nsel-1)
          factor1=c2/var
2         IF(factor1.LT.0.1) factor1=0.1
          IF(factor1.GT.6.0) factor1=6.0
      else
          factor1=1.0
      end if

C Compute covariance matrix of observations
      DO  I = 1,NSEL
          II1=I*(I-1)/2
          XI = XTMP(I)
          YI = YTMP(I)
          AZM= aztmp(i)
          DO  J = 1,I-1
              tmp=aztmp(j)
              C(II1+J)=COVEE(XI,YI,XTMP(J),YTMP(J),AZM,tmp,distance)
          END DO
C DIAGONAL ELEMENT
          C(II1+I)=VAR+stdtmp(i)/factor1
c      C(II1+I)=VAR+stdtmp(i)
      END DO
      ZPRED1 = 0
      ZPRED2 = 0

      if(lerr) then
 
          CALL SOLUTION(C,ARR1,ARR2,NSEL,NSING)

          if(nsing.eq.1) then
              npp=npp+1
              write(100,*)real(gridx),real(gridy)
              ARRY(I1,J1)=999.0
              ARRX(I1,J1)=999.0
              ARRY_err(I1,J1)=999.0
              ARRX_err(I1,J1)=999.0
              go to 110
          end if

          xi_err=0
          eta_err=0
          DO I = 1,NSEL
              TMP=ztmp(i)
! Compute gradients     
              ZPRED1 = ZPRED1 + ARR1(I)*TMP
              ZPRED2 = ZPRED2 + ARR2(I)*TMP
! Compute error estimates of gradients
              xi_err=xi_err+arr1(i)*xarr1(i)
              eta_err=eta_err+arr2(i)*xarr2(i)
     
          END DO
          xi_err=sqrt(abs(var-xi_err)*factor1) 
          eta_err=sqrt(abs(var-eta_err)*factor1)
 
      else   
    
          call  solve(c,ztmp,nsel,0,c,work,nsing)
          if(nsing.eq.1) then
              npp=npp+1
              write(100,*)real(gridx),real(gridy)
              ARRY(I1,J1)=999.0
              ARRX(I1,J1)=999.0
              ARRY_err(I1,J1)=999.0
              ARRX_err(I1,J1)=999.0
              go to 110
          end if
        


          do i=1,nsel
              TMP=ztmp(i)
              ZPRED1 = ZPRED1 + ARR1(I)*TMP
              ZPRED2 = ZPRED2 + ARR2(I)*TMP
          end do

      end if

C Accept only reasonable values

      if(dabs(ZPRED1).gt.100.0 .or.dabs(ZPRED2).gt.100.0 ) then
            ARRY(I1,J1)=999.0
            ARRX(I1,J1)=999.0
            ARRY_err(I1,J1)=999.0
            ARRX_err(I1,J1)=999.0
            go to 110
      end if

! Convert arc-sec to micro-radian
      ARRY(I1,J1) = ZPRED1*factor2
      ARRX(I1,J1) = ZPRED2*factor2
      if(lerr) then
          ARRY_err(I1,J1) = xi_err*factor2
          ARRX_err(I1,J1) = eta_err*factor2
      end if
      sum_npt=sum_npt+nsel
      sum_fact=sum_fact+factor1
      npred=npred+1
110   ENDDO
109   ENDDO
 
301   CONTINUE  ! This is the main loop

! END OF PREDICTION 
      write(*,*) trim(file_i)
      N = npred
      WRITE(6,250) N,nx*ny-n 
      write(6,*)'Averaged number of points used for computation:',
     & sum_npt/npred
      write(6,*)'Averaged scale factor:', sum_fact/npred  
      WRITE(6,*)'Number of singular cases:', npp
250   FORMAT('predicted/empty grids: ',2I10)

       write(0,*)'Filling data gaps over the oceans ...'
      CALL FILL(ARRX,IDIM,IDIM)
      CALL FILL(ARRY,IDIM,IDIM)

      if(lerr) then
      call fill(arrx_err,IDIM,IDIM)
      call fill(arry_err,IDIM,IDIM)
      end if       
      
! Compute statistics
      call stat(arrx,idim,nx,ny,result)
      write(6,*)'mean,std dev, max, min (micro-rad) of north component:'
      write(6,'(4f10.4)')(result(i),i=1,4)
      call stat(arry,idim,nx,ny,result)
      write(6,*)'mean, std dev, max, min (micro-rad) of east component:'
      write(6,'(4f10.4)')(result(i),i=1,4)

      if(lerr) then
      call stat(arrx_err,idim,nx,ny,result)
      write(6,*)'mean, std dev, max, min (micro-rad) of north component
     &(error):'
      write(6,'(4f10.4)')(result(i),i=1,4)
      call stat(arry_err,idim,nx,ny,result)
      write(6,*)'mean, std dev, max, min (micro-rad) of east component 
     &(error):'
      write(6,'(4f10.4)')(result(i),i=1,4)
      end if

      WRITE(0,*)'Writing .grd3 files ... '  

! North component
      WRITE(60)NX,NY,WEST,EAST,SOUTH,NORTH,FX,FY
      WRITE(60)((ARRX(I,J),J=1,NY),I=1,NX)
      if(lerr)WRITE(60)((ARRX_err(I,J),J=1,NY),I=1,NX)

! east component
      WRITE(61)NX,NY,WEST,EAST,SOUTH,NORTH,FX,FY
      WRITE(61)((ARRY(I,J),J=1,NY),I=1,NX)
      if(lerr)WRITE(61)((ARRY_err(I,J),J=1,NY),I=1,NX)
      
      WRITE(0,*)'Writing .txt files ... '  
!NX=(east-west)/fx+1.01
!NY=(north-south)/fy+1.01
      do i=1,nx
          lon=west+DFLOAT(i-1)*FX
          do j=1,ny
              lat=south+DFLOAT(j-1)*FY
              dov_n=ARRX(I,J)
              dov_e=ARRY(I,J)
              write(62,1000) lon,lat,dov_n
              write(63,1000) lon,lat,dov_e
          enddo
      enddo
1000  format(3f12.4)

 
      CALL SECOND(TIMEX)
      WRITE(6,1004)TIMEX-TIME0
 1004 FORMAT('TOTAL EXECUTION TIME (SECONDS):',F12.3)
      pause
      STOP
      END

      subroutine stat(a,dim,nx,ny,result)
! Compute max, min, std dev, average of data in array a
      implicit none
      integer nx,ny,i,j,n,dim
      real*4 a(dim,*),maxg,ming,result(*),flag
      real*8 sum,sum2
!
      do i=1,4
      result(i)=999.0
      end do
      maxg=-99999.0
      ming=99999.0
      n=0
      sum=0
      sum2=0
      do j=1,ny
      do i=1,nx
       n=n+1
      sum=sum+a(i,j)
      sum2=sum2+a(i,j)**2
      maxg = amax1(maxg,a(i,j))
      ming = amin1(ming,a(i,j))
       end do
      end do
      if(n.ne.0) then
      result(1)=sum/n ! mean
      result(2)=sqrt((sum2-result(1)**2/n)/(n-1)) ! std dev
      result(3)=maxg ! maximum
      result(4)=ming ! minimum
      end if
      return
      end 
      
C
C
      SUBROUTINE SOLUTION(A,X,Y,N,nsing)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C PROGRAM TO SOLVE FOR A LINEAR SYSTEM USING CHOLESKY DECOMPOSTION.
C TWO SOLUTIONS ARE DONE AT A TIME.
C
C       A ... COVARIANCE MATRIX BETWEEN OBSERVATIONS.
C       X ... COVARIANCE VECTOR BETWEEN ETA (NORTH COMPONENT OF
C             GEOID GRADIENTS) AND THE OBS.
C       Y ... COVARIANCE VECTOR BETWEEN XI (EAST COMPONENT OF 
C             GEOID GRADIENTS ) AND THE OBS.
C       N ... NUMBER OF OBS.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8  A(*),SUM, SUMX,SUMY,X(*),Y(*)
      integer nsing
C CHOLESKY DECOMPOSITION
      nsing=0
      DO 1 K=1,N
      kk=k*(k-1)/2
      DO 2 I=1,K-1
      II=I*(I-1)/2
      SUM=0.D0
      DO 3 J=1,I-1
 3    SUM=SUM+A(II+J)*A(KK+J)
      II1=KK+I
 2    A(II1)=(A(II1)-SUM)/A(II+I)
      SUM=0.D0
      DO 4 J=1,K-1
      IK=KK+J
 4    SUM=SUM+A(IK)*A(IK)
      IK=KK+K
      A(IK)=A(IK)-SUM
      IF(A(IK).LE.0.0) THEN
c      STOP'COVARIANCE MATRIX IS NOT POSITIVE DEFINITE'
c      write(6,*)'a(ik)=',k,a(ik)
!      A(IK)=1.D10
      nsing=1
      return
      END IF
1     A(IK)=DSQRT(A(IK))
c FORWARD SUBSTITUTION
      DO I=1,N
      II=I*(I-1)/2
      SUMX=0.D0
      SUMY=0.D0
      DO J=1,I-1
      SUMX=SUMX+A(II+J)*X(J)
      SUMY=SUMY+A(II+J)*Y(J)
      END DO
      X(I)=(X(I)-SUMX)/A(II+I)
      Y(I)=(Y(I)-SUMY)/A(II+I)
      END DO
C BACKWARD SUBSTITUTION
      DO J=N,1,-1
      SUMX=0.D0
      SUMY=0.D0
      DO I=J+1,N
      II=I*(I-1)/2
      SUMX=SUMX+A(II+J)*X(I)
      SUMY=SUMY+A(II+J)*Y(I)
      END DO
      JJ=J*(J+1)/2
      X(J)=(X(J)-SUMX)/A(JJ)
      Y(J)=(Y(J)-SUMY)/A(JJ)
      END DO
      RETURN
      END
C
      FUNCTION COVEE( X1,Y1,X2,Y2,A1,A2,dist)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C FUNCTION TO CALCULATE THE COVARIANCE BETWEEN
C GEOID GRADIENT WITH AZIMUTH A1 AT X1, Y1,  AND GEOID GRADIENT WITH 
c AZIMUTH A2 AT X2, Y2 .
C
C IF THE TWO POINTS ARE IDENTICAL THEN THE COVARIANCE IS 
C CLL(0) * COS (A1-A2).
C
C X1,Y1,X2,Y2 ARE IN DEGREES WITH X1 AND X2
C BEING MULTIPLIED BY A FACTOR COS (LAT),
C AND A1 AND A2 IN RADIANS
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,O-Z)
c
C COVARIANCES CLL AND CMM ARE STORED IN ARRAYS  CLL(.) AND CMM(.) WITH
C A STEPSIZE PSI (IN DEGREES)
      COMMON /COVTAB/CLL(1001),CMM(1001),PSI,PSIIV,VAR
      REDTOR=111.1949266D0
      DTOR=0.01745329252D0
C      COSFI=DCOS((Y2+Y1)*0.5D0*DTOR)      
c      DX=(X2-X1)*REDTOR
c      DY=(Y2-Y1)*REDTOR
      DX=X2-X1
      DY=Y2-Y1
      DIST=DSQRT(DX*DX+DY*DY)
C IF THE DISTANCE BETWEEN THE TWO POINTS ARE LESS THAN 
C 1 METER= 9.E-6 DEGREE,
C THEY ARE REGARDED AS THE SAME POINT
      IF(DIST.LT.(9.D-6)) THEN
c      write(6,*)'dist=',dist
      COVEE=VAR*DCOS(A1-A2)
      RETURN
      END IF
C
      AL=DATAN2(DX,DY)
      ID=DIST*PSIIV+1.01
c      write(6,*)dist,id
      IF(id.gt.501)then
      covee=0
      return
      end if
c      write(6,*)id,dist
      TEMP=PSIIV*(DIST-(ID-1)*PSI)
      COVLL=CLL(ID)+(CLL(ID+1)-CLL(ID))*TEMP
      COVMM=CMM(ID)+(CMM(ID+1)-CMM(ID))*TEMP
C PLANAR APPROXIMATION: ALPHA (P,Q) = ALPHA (Q,P) + PI
      DALPHA1=A1-AL
      DALPHA2=A2-AL
      COVEE=COVLL*DCOS(DALPHA1)*DCOS(DALPHA2)+COVMM*
     & DSIN(DALPHA1)*DSIN(DALPHA2)
      RETURN
      END
c
      SUBROUTINE SECOND(TIME)
      REAL*8 TIME
      REAL*4 T(2)
      DATA TOT/0.D0/
      TIME=DTIME(T)+TOT
      TOT=TIME
      RETURN
      END
C================================================================
C
      function landmask(lat,lon)
      real*8 lat,lon,lon1
      integer*2 i,j
      logical landmask
      logical*1 mask(4320,2160)
      common /landocean/mask
      lon1=lon
      if(lon1.gt.360.0)lon1=lon1-360.0
      if(lon1.lt.0.0)lon1=lon1+360.0
      i=idnint(lon1*12.d0)+1
      if(i.gt.4320) i=4320
      j=idnint((lat+90.d0)*12.d0)+1
      if(j.gt.2160)j=2160
C land=1, water=0
      landmask=mask(i,j)
      return
      end

C=================================================================
C
 !     subroutine readmask
 !     logical*1 mask(4320,2160)
 !     common /landocean/mask
!      lent=4320*2160
! glon_mask_01 is coded in sun binary
!      open(88,file='glob_mask_01'
!     &,RECL=LENT,ACCESS='DIRECT')
!      read(88,rec=1)((mask(i,j),i=1,4320), j=1,2160)
!      return
!      end

C=================================================================
   
      subroutine fill(arr,idim,idim1)
C Subroutine to fill the data gaps over the oceans
      implicit real*8(a-h,o-z)
      integer*4 idim, ix, iy, k, l, i, j,window,idim1,
     & iend,istart,jend,jstart
      real*4 arr(idim,idim1)
!      LOGICAL LANDMASK
!      logical*1 mask(4320,2160)
      real*8 sum,sum1,xmin,xmax,ymin,ymax,gridx,gridy,lat,lon
      common /grdinfo/ix,iy,xmin,xmax,ymin,ymax,gridx,gridy
      common /landocean/mask

! initial window size
      window=3
      REDTOR=111.1949266D0
      DTOR=0.01745329252D0
      d=100 !km
      do 12 i=1,ix
      lon=xmin+dfloat(i-1)*gridx
!      if(lon.gt.360.d0) lon=lon-360.d0
      do 12 j=1,iy
      lat=ymin+dfloat(j-1)*gridy
      if(arr(i,j).ne.999.0) go to 12
  
14    sum=0.d0
      sum1=0.d0
      istart=i-window
      if(istart.le.0)istart=1
      iend=i+window-1
      if(iend.gt.ix)iend=ix
      jstart=j-window
      if(jstart.le.0)jstart=1
      jend=j+window-1
      if(jend.gt.iy)jend=iy
      x1=lon
      y1=lat
      n=0
      do 13 l=jstart, jend
      y2=ymin+dfloat(l-1)*gridy
      cosfi=dcos((y1+y2)/2.d0*dtor)
      do 13 k=istart, iend  
      if(k.eq.i .and. l.eq.j) goto 13 
      if(abs(arr(k,l)).eq.999.0) goto 13
        x3=xmin+dfloat(k-1)*gridx
        x2=x3
        if(x2.gt.360.d0) x2=x2-360.d0
        if(x2.lt.0.d0) x2=x2+360.d0
 !       if(landmask(y2,x2)) goto 13
        xdis=(x1-x3)*redtor*cosfi
        ydis=(y1-y2)*redtor
        dist2=xdis**2+ydis**2
        
        wt=1.d0/dist2
        sum=sum+ arr(k,l)*wt
        sum1=sum1+wt
        n=n+1
13    continue
      if(n.lt.4 .and. window.lt.6) then
        window=window+1
        goto 14
      endif
      
       if(sum1.ne. 0.) then
       arr(i,j)=sngl(sum/sum1)
       else
       arr(i,j)=0.
      end if

12    continue
      return
      end

      SUBROUTINE solve(A,X,N,mode,C,work,nsing)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C PROGRAM TO SOLVE FOR A LINEAR SYSTEM USING CHOLESKY DECOMPOSTION.
C
C    A ... A psd matrix, or the normal matrix. On output, it
c          becomes a lower triangular matrix from the Cholesky
c          decomposition or its inverse, C
C    X ... on input it is the observation vector; on output
c          it is the solution vector
c    C ... inverse of A
c    work .. a work array with dim=N
C    N ... order of A
c    mode .. 0, solve for x only
c            1, invert A (the inverse is C) and solve for x
c            2, invert A only. X is not used in this mode
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      REAL*8  A(*),SUM,X(*),work(*),C(*)
      integer order,nsing
C CHOLESKY DECOMPOSITION of A
      nsing=0
      if(mode.lt.0 .or. mode.gt.3) stop'incorrect mode'
      DO 1 K=1,N
      kk=k*(k-1)/2
      DO 2 I=1,K-1
      II=I*(I-1)/2
      SUM=0.D0
      DO 3 J=1,I-1
 3    SUM=SUM+A(II+J)*A(KK+J)
      II1=KK+I
 2    A(II1)=(A(II1)-SUM)/A(II+I)
      SUM=0.D0
      DO 4 J=1,K-1
      IK=KK+J
 4    SUM=SUM+A(IK)*A(IK)
      IK=KK+K
      A(IK)=A(IK)-SUM
      IF(A(IK).LE.0.0) THEN
      write(0,*)'negative diagonal element a(ik)=',k,a(ik)
c      A(IK)=1.D10
      nsing=1
!      stop'matrix is not positive definite'
      return
      END IF
1     A(IK)=DSQRT(A(IK))
      if(mode.eq.0 .or.(mode.eq.1)) then
      call LTRISOL(A,X,N)
      call UTRISOL(A,X,N)
      end if
c Invert A
      if((mode.eq.1) .or. (mode.eq.2)) then
      do order=1,n
      do i=1,n
      work(i)=0
      end do
      work(order)=1
      call LTRISOL(A,work,N)
      Call UTRISOL(A,work,N)
      do i=order,n
      ii=i*(i-1)/2
      c(ii+order)=work(i)
      end do
      end do
      do i=1,n*(n+1)/2
      a(i)=c(i)
      end do
      end if
c      write(0,*)'Nsing=',nsing
      RETURN
      END

      SUBROUTINE LTRISOL(A,X,N)
C SUBROUTINE TO SOLVE A LINEAR SYSTEM. A IS A LOWER TRIANGULAR MATRIX.
C INPUT :
C        A ----- A NONSINGULAR LOWER TRIANGULAR MATRIX
C        N ----- ORDER OF T
C        X  ---- OBSERVATION VECTOR
C OUTPUT :
C        X ----- SOLUTION VECTOR
C
       REAL*8 A(*),X(*)
       REAL*8 SUM
       DO 1 I=1,N
       II=I*(I-1)/2
       SUM=0.D0
       DO 2 J=1,I-1
 2     SUM=SUM+A(II+J)*X(J)
1      X(I)=(X(I)-SUM)/A(II+I)
       RETURN
       END
      SUBROUTINE UTRISOL(A,X,N)
C SUBROUTINE TO SOLVE A LINEAR SYSTEM. A IS AN UPPER
C TRIANGULAR MATRIX FROM THE CHOLESKY DECOMPOSITION.
C INPUT :
C        A ----- A NONSINGULAR LOWER TRIANGULAR MATRIX
C        N ----- ORDER OF A
C        X  ---- OBSERVATION VECTOR
C OUTPUT :
C        X ----- SOLUTION VECTOR
C
       REAL*8 A(*),X(*)
       REAL*8 SUM
C BACKWARD SUBSTITUTION
      DO J=N,1,-1
      SUM=0.D0
      DO I=J+1,N
      II=I*(I-1)/2
      SUM=SUM+A(II+J)*X(I)
      END DO
      JJ=J*(J+1)/2
      X(J)=(X(J)-SUM)/A(JJ)
      END DO
      RETURN
      END

      SUBROUTINE INDEXX(N,ARRIN,INDX)
      DIMENSION ARRIN(N),INDX(N)
      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      IF(N.EQ.1)RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
      END

      SUBROUTINE SLTPTV(X,Y,N,TOLDST,NSKIP,INDXX)
C A VERSION OF SUBROUTINE SLTPT, BUT IT EXECUTES FASTER THAN SLTPT
C DOES. MORE POINTS MIGHT BE DELETED USING THIS ALGORITHM.
C ** TOLDST IS THE SQUARE OF THE TOLERABLE DISTANCE.
      IMPLICIT REAL*8(A-H,O-Z)
	real*8 x(*),y(*)
      DIMENSION INDXX(N*(N-1)/2)
      NSKIP=0
      DO   I=1,N-1
      DO   J=I+1,N
	dist2=(X(I)-X(J))**2+(Y(I)-Y(J))**2
      IF(dist2.le.TOLDST)then
      NSKIP=NSKIP+1
      INDXX(NSKIP)=J
	end if
      end do
      end do
      RETURN
      END
c=========================================================================
      subroutine write_error
      write(0,*)' '
      write(0,*)'Usage: gridg f.xyg  -Rwest/east/south/north -Idx/dy 
     & -Wwindow -Nnorth.grd1 -Eeast.grd1 [-S -O]> gridg.sta' 
      write(0,*)'f.xyg file containing geoid gradients (in sec) 
     & in ascii form (see rmgrad.f)'
      write(0,*)'-R for west, east, south, north'
      write(0,*)'-I grid interval in long and lat directions in minutes'
      write(0,*)'-N ouput gridded geoid gradient in north-south 
     &	direction'
      write(0,*)'-E ouput gridded geoid gradient in west-east direction'
      write(0,*)'-S scale global covariance function[default: No scale]'
      write(0,*)'-O compute error estimates of gradients[default: no]'
	write(0,*)'-W: size of window (in minutes) for data selection'
      return
      end

