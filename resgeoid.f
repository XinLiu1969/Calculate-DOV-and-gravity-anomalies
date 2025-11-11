      program resgeoid
	use DFLIB


! remove reference geoid and compute along-track residual gradients 去除参考大地水准面并计算沿轨剩余大地水准面梯度
!
!
! DESCRITPION
! -F: contain input files
! -G: file of residual geoid gradients 
! -M: reference geoid model in .grd3   
! -I: sst model in .grd3               
! -O: file of outliers                 
! -T: threshold                        
! -N: Satellite ID                    
!
! The output ascii file (-G) contains:  
!
!     latitude (degree)
!     longitude (degree)
!     geoid gradient (arc-second)
!     azimuth of geoid gradient (radian)
!     std dev of geoid gradient (arc-second)
!     satellite id
!
!
      implicit none 
	integer*4 narg,i,npass,sat_id,np,iargc,nmax,dim,dim1,nx,ny,j,deg
     &,numg,n,Nm,ii,ib,ie,npt,nx1,ny1,maxindex,k,total,ncout,fcout,numf,
     &numr
      parameter(dim=40*361+1,nmax=800000,dim1=30*42+1) !dim=40*361+1
	real*4 geoidsst(dim,dim),geoidsst1(dim,dim)
      character*200 buf,filelist,ofile1,ofile2,INF,outf,filename,ifile1
     &,outf1,ifile2
      real*8 plo,pla,ssh,sstd,dd,mg,std,thresh,summ,summ1,mean1
     &,lon(nmax),lat(nmax),h(nmax),stdh(nmax),refgeoi(nmax)
     &,xmin,xmax,ymin,ymax,dx,dy,val,az(nmax),dista(nmax),smin
     &,xmin1,xmax1,ymin1,ymax1,dx1,dy1,val1,gh
     &,grad_std(nmax),glon(nmax),glat(nmax),grad_ht(nmax)
     &,grad_res(nmax),p(nmax),delh_res(nmax)
     &,geoid(nmax),grad_az(nmax),grad_geoid(nmax),diff
     &,grad_refgeoi(nmax),maxdiff,delh(nmax),delh_geoid(nmax)
     &, rmss,stdd,stddsum,grad_res_total(30000000),factor2
c	REAL*8, ALLOCATABLE ::  correct(:)
	logical EXISTS,P1(nmax),ret
	data deg/6/                                  
	data smin/8000.d0/ ! minimum acceptable distance
	                                                 

! get command-line arguments 
 !     narg=iargc()                                     
 !     if(narg.lt.1) then                             
 !     call write_error
 !     stop
 !     end if
 !以上为命令行下运行Fortran命令行语句，可改为windows

!      read (*,*) buf
!	do i=1,narg
!		!call getarg(i,buf)
!		select case (buf(1:2))
!		case ('-F','-f')
!			filelist=buf(3:)
!	    case ('-N','-n')
!			read(buf(3:),*)N
!		case ('-G','-g')
!			ofile1=buf(3:)
!		CASE ('-O','-o')  
!			ofile2=buf(3:)
!	    case ('-M','-m')
!			ifile1=buf(3:)
!         case ('-I','-i')
!			ifile2=buf(3:)
!		case ('-T','-t')
!			read(buf(3:),*)thresh
!		case default
!	write(*,*)'resgeoid -Ffilelist -N -M -Goutfile -Tthresh -O'
!			stop
!	    end select
!     end do
      ncout=0
	fcout=0
      filelist="cyclepasspathlist.dat"              
	sat_id=11
	ofile1="re_geiod_gradients_L12.dat"
	ofile2="outliers_IS2.dat"
!	ifile1="EGM2008gxinm1.grd3"
	thresh=10                                             
      factor2=1.d0/0.2062648062d0
! Output file of residual geoid gradients
      open(60,file=ofile1)              
! Output file of outliers
      open(61,file=ofile2)
! read geoid file in .grd3 format
 !-------------------------------------------------------------------     
 !     open(10,file=ifile1,form='unformatted')
 !     read(10)nx,ny,xmin,xmax,ymin,ymax,dx,dy
 !     write(0,*)nx,ny,xmin,xmax,ymin,ymax,dx,dy
 !     if(nx.gt.dim)stop'increase dim'
 !     if(ny.gt.dim)stop'increase dim'
 !     read(10)((geoidsst(i,j),j=1,ny),i=1,nx)
 !     close(10)
	!write(*,*)filelist,ofile1,ofile2,smin,nx
      open(50,file=filelist,status='old')
      stddsum=0.0 
      total=0
21    continue
	read(50,*,end=166)filename               
      inf=TRIM(filename)
      !if(filename(68:69).NE.'3r') GOTO 21   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!I:\altimeter_data\4-removeMDT\GM\IS-2\cycle_0015\C0015_P0161_3L.dat
      write(*,*)filename
      npass=npass+1
      INQUIRE(FILE = INF, EXIST = EXISTS)
      IF(.NOT. EXISTS) GOTO 3000
      open(12,file=INF,status='old')
      np=0
1     continue	
      read(12,*,end=66)plo,pla,ssh,gh,sstd         !,sat_id
      sstd=0.089    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      call interp2d(geoidsst,nx,ny,dim,dx,dy,xmin,ymin,deg,pla,plo,val) 
c         call interp2d(geoidsst1,nx1,ny1,dim,dx1,dy1,xmin1,ymin1,deg          
c     &     ,plo,pla,val1)
      if(gh.ne.999999) then      
          np=np+1 
	    lon(np)=plo
	    lat(np)=pla
          h(np)=ssh
          geoid(np)=gh
          stdh(np)=sstd  
      endif
      goto 1
66    continue
      close(12)
	if(np.lt.1) goto 3000
	call along_dist(lon,lat,np,dista,az)!
	call filter1d(dista,h,np,8.d0,1,p,numf) 
	call gradient(az,dista,lon,lat,h,stdh,np,smin,
     &     grad_ht,delh,glon,glat,grad_az,grad_std,numg) 
      call gradient(az,dista,lon,lat,geoid,stdh,np,smin,
     &     grad_geoid,delh_geoid,glon,glat,grad_az,grad_std,numg)

	! compute gradient of reference geoid

      grad_res=grad_ht -grad_geoid 
      delh_res=delh-delh_geoid 
      ncout=ncout+numg
      do  i=1,numg
	    if(abs(grad_res(i)).lt.thresh)then
              write(60,'(2f12.6,3f12.6,I6)') 
     &glon(i),glat(i),grad_res(i),grad_az(i),grad_std(i),
     % sat_id
	        fcout=fcout+1
              grad_res_total(fcout)=grad_res(i)!*factor2
	    else
	        write(61,'(2f12.6,2f12.3)')glon(i),glat(i),grad_res(i),
     &	    delh_res(i) 
	    endif
      enddo 
      total= total+1
3000  continue 
      goto 21
166   continue
      close(50)
      close(61)
      close(60)
      write(*,*) ncout,fcout
      call tongji(grad_res_total,fcout)
	pause
      end

      subroutine along_dist(lon,lat,n,dista,az)
! Program to compute along-track distance (in km) using spherical formula
! Input:
! lon,lat: longitude and latitude
! n: number of data points
! Output:
! dista: along-track distance in meter
! az: along-track azimuth in radian (up to index=n-1)
! 
      implicit none
      integer n,i
      real*8 lon(n),lat(n),dista(n),az(n)
      real*8 latm,dlon,dlat,d2r,re,dd,a,e2
      data d2r,re/0.017453292d0,6371000.d0/
! TOPEX ellipsoidal parameters
   !   data a,e2/6378136.3d0,0.00669438002290d0/
! WGS84 ellipsoidal parameters
      data a,e2/6378137.0d0,0.00669437999013d0/
      dista(1)=0
      do i=2,n
      latm=(lat(i)+lat(i-1))*0.5*d2r
! radius = gauss mean radius at average latitude
      re=a*dsqrt(1.d0-e2)/(1.d0-e2*(dsin(latm))**2)
      dlat=(lat(i)-lat(i-1))*d2r
      dlon=(lon(i)-lon(i-1))*dcos(latm)*d2r
!      dd=dsqrt(dlat**2+dlon**2)*re
! spherical distance       
      dd=re*dacos( dsin(lat(i)*d2r)*dsin(lat(i-1)*d2r)+
     &  dcos(lat(i)*d2r)*dcos(lat(i-1)*d2r)*dcos((lon(i)-lon(i-1))*d2r))
      dista(i)=dista(i-1)+dd
      az(i-1)=datan2(dlon,dlat) !max i is n-1 
      end do
      return
      end

      subroutine gradient(az,dista,lon,lat,h,err,n,smin,
     & g,delg,glon,glat,gaz,gstd,k)
! Direct computation of slope as gradient
! Input:
! az: along-track azimuth in radian
! dista: along-track distance in meter
! lon,lat: longitude and latitude in degree
! h: along-track ssh in meter
! mh: mean ssh in meter
! err: std dev of h  in meter
! n: number of ssh
! smin: minimum acceptable distance  (in meter) for computing g
!
! Output:
! g: gradient in arc-second
! glon,glat: lon and lat of g in degree
! gaz: azimuth of g in radian
! gstd: std of g in arc-second
! k: no. of accepted gradients

      implicit none
      integer i,k,n
      real*8 dista(n),lon(n),lat(n),h(n),g(n),glon(n),glat(n),
     &  az(n),gaz(n),gstd(n),err(n),mh(n),mg(n),delg(n)
      real*8 d2r,re,dlon,dlat,diff,smin,rho
      data re,d2r/6371000.d0,0.017453292d0/
      data rho/206264.8062D0/
      k=0
      
	do 11 i=1,n-1
          diff=dista(i+1)-dista(i)
          if(diff.gt.smin) go to 11
          k=k+1
          glat(k)=(lat(i)+lat(i+1))/2.d0
          glon(k)=(lon(i)+lon(i+1))/2.d0
          g(k)=(h(i+1)-h(i))/diff*rho
          delg(k)=(h(i+1)-h(i))
c mg is reference gradient

c      mg(k)=(mh(i+1)-mh(i))/diff*rho
          dlon=(glon(k)-lon(i))*cos(glat(k)*d2r)*d2r
          dlat=(glat(k)-lat(i))*d2r
          gaz(k)=az(i)
          gstd(k)=dsqrt(err(i+1)**2+err(i)**2)/diff*rho
11    end do
      return
      end

      subroutine interp2d(a,nx,ny,dimx,dx,dy,x1,y1,deg,x,y,vint)

c Porgram to do two dimensional interpolation using 
c divided difference with a regular grid
c
c Input:
c a ... data array in .grd3 format. SP
c nx,ny ... used dimensions of a along x and y. INT
c dimx ... dimension of a exactly as declared in the calling program. INT
c dx,dy ... grid intervals in a. DP
c x1,y1 ... minimum values of x and y coordinates (or west and south). DP
c x,y ... x and y coordinates where interpolation is wanted. DP
c degree ... number of data points used; the higher 
c            the more accurate. DP
c 
c Output: 
c vint ... interpolated value at x, y. DP
c
c Notes: 
c SP ... single precision
c DP ... double precision
c INT ... integer*4
c SINT ... integer*2

      implicit none 
      integer nx,ny,dimx,degree,deg,k,k1,k2,kx,ky,k1x,k1y,i,j
      real*4 a(dimx,*)
      real*8 dx,dy,x1,y1,x2,y2,x,y,vint,tmp,x0,y0,xa(50),ya(50)
      x2=x1+(nx-1)*dx
      y2=y1+(ny-1)*dy
      if(x.lt.x1 .or. x.gt.x2) then
      write(0,*)'station is not inside the grid, int. value =999999.0'
      vint=999999.d0
      return
      end if
      if(y.lt.y1 .or. y.gt.y2) then
      write(0,*)'station is not inside the grid, int. value =999999.0'
      vint=999999.d0
      return
      end if
      degree=(deg/2)*2
      k=(x-x1)/dx+1.01
      k1=k-degree/2+1
      if(k1.lt.1)k1=1
      k2=k+degree/2
      if(k2.gt.nx)k2=nx
      kx=k2-k1+1
      x0=x1+(k1-1)*dx
      k1x=k1-1
      k=(y-y1)/dy+1.01
      k1=k-degree/2+1
      if(k1.lt.1)k1=1
      k2=k+degree/2
      if(k2.gt.ny)k2=ny
      ky=k2-k1+1
      y0=y1+(k1-1)*dy
      k1y=k1-1
c Interpolate to y
      do j=1,ky
      do i=1,kx
      xa(i)=a(k1x+i,k1y+j)
      end do
      call divint(x0,xa,kx,dx,x,tmp)
      ya(j)=tmp
      end do
c one final interpolation
      call divint(y0,ya,ky,dy,y,tmp)
      vint=tmp
      return
      end

      SUBROUTINE DIVINT(X0,YA,N,H,X,Y)
C==============================================================
C   Polynomial interpolation using divided Difference
C==============================================================
C   Veriables :
C       X0    ==> coordinate of the first point ( Input )
c                 ie, its value is ya(1)
C       YA(N) ==> Data             ( Input )
C       N     ==> Number of points ( Input )
C       H     ==> Stepsize         ( Input )
C       X     ==> coordinate where int is needed ( Input )
C       Y     ==> Result           ( Output )
C       DF    ==> Work array
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION YA(N),DF(50,50)
      if(n.gt.50) stop'increase dim of DF'
      S=(X-X0)/H
c If x happens to be on the knot, then s is an integer
c
      if(abs(int(s)-s).lt.1.0d-7) then
      y=ya(int(s)+1)
      return
      end if

      S1=1.D0
      S2=1.D0
      DO J=1,N
         DF(1,J)=YA(J)
      ENDDO
      Y=YA(1)
      DO  I=2,N
      DO  J=1,N-I+1
         DF(I,J)=DF(I-1,J+1)-DF(I-1,J)
      ENDDO
      di=dfloat(i)
         S1=S1*(S-di+2.d0)
         S2=S2*(di-1.d0)
         DY=S1/S2*DF(I,1)
         Y=Y+DY
      ENDDO
      RETURN
      END


      subroutine filter1d(x,y,np,window,mode,p,num)
 
c
c A space (time) domain filter and outlier cleaner for a one-d time series.
c
c Input:
c x, y: time series of x, y pair, x must be in the increasing order. y=y(x)
c npt: number of data points, may be changed if mode=2
c window: size of window (same unit as x's)
! mode: 1= filter+reject outliers, 2= reject outliers only
 
c filter_type: 1 = gauss, 2 = average
c
c Output:
c y: y is replaced by a filtered time series 
c p: index of y for outliers
c Notes:
c
c 1. The Gaussian and average filter formulae used are:
c 
c y_filter(k) = [SUM (wt(i)*y(i))] / [SUM (wt(i)]
c
c where 
c wt(i) = exp(-(ds/sigma)**2) or 1, ds = abs(x(k)-x(i)), 
c i = i_min,...,k-1,k,k+1,..., i_max. 
c i_min and i_max are indices fulfilling ds<window/2, sigma=window/6


      implicit none
      integer*4 npt,i,k,mode,max_i,j,number,np,num
	parameter(npt=80000)
       real*8 x(npt),y(npt),p(npt),sum2,diff,weight,max_diff,
     & wt,sum,ds,dels,std,std3,sigma2,tmp,window
      logical loutlier
      real*8 tdat(npt)
c     weight(filter_type,ds,sigma)=exp(-0.5*(ds/sigma)**2)

      if(npt.gt.100000) stop'increase the dim of p'
      number=np
      dels=0.5*window*1000
      sigma2=(window*1000/6.d0)**2
c	write(*,*) window,np,number
c       mostnum=int(0.1*npt)
      do i=1,number
      p(i)=1.0
      end do     

200   continue
      do k=1,number
		 if(p(k).gt.2.0) goto 7
		j=0
		wt=0.0
		sum=0.0
		do i=k,1,-1
			ds=abs(x(k)-x(i))
			if(ds.gt.dels) go to 5
			if(p(i).lt.2.0) then
c	 tmp=0.5+0.5*dcos(3.1415926d0*ds/dels)
				tmp=dexp(-0.5*ds*ds/sigma2)!?*0.5
				wt=wt+tmp
				sum=sum+y(i)*tmp
				j=j+1
			end if
		end do
5         do i=k+1,number
			ds=abs(x(i)-x(k))
			if(ds.gt.dels) go to 6
			if(p(i).lt.2.0) then
				tmp=dexp(-0.5*ds*ds/sigma2)!
c      tmp=0.5+0.5*dcos(3.1415926d0*ds/dels)
				wt=wt+tmp
				sum=sum+y(i)*tmp
				j=j+1
			end if
		end do
6         continue
		if(j>1)then
		 tdat(k)=sum/wt
		else
		 p(k)=9999.0
		endif
7         continue
      enddo
  
      loutlier=.false.
      sum2=0.0
	max_diff=-9999.0
	k=0
      
	if(mode.eq.1) then
! Filter + reject outliers
		num=0
		do i=1,number
			if(p(i).le.2.0) then
			 num=num+1
			 y(i)=tdat(i)
			endif
		end do
      end if

      return
	end
      

      subroutine write_error
      write(0,*) 'error'
	! SYNOPSIS
	write(*,*)'h2g -Mgeoid_mdoel.grd3 -Gout.xyg -Ooutlier -Tthresh'
!
! DESCRITPION
! file_list: list of .xyh files. A .xyh file contains:
!           longitude, latitude, ssh, std dev of ssh, satellite id
!
	write(*,*)' -G: file of residual geoid gradients'
	write(*,*)' -M: reference geoid model in .grd3'
	write(*,*)' -O: file of outliers'
! -P: file that altimetry geoid gradients minus reference geoid gradients
	write(*,*)' -T: threshold'
      return
      end
      
      
      subroutine tongji(dssh,n)
       integer count,n,num
       real*8 dssh(n),kssh(30000000)
       real*8 lat,lon,rmax,rmin,rmean,std,sum1,rms
       real*8 rmax1,rmin1,sumssh,dssh1,dssh0
       rmax=0
       rmin=0
       rmean=0
       sum1=0
       count=n
      do i=1,n
            rmean=dssh(i)+rmean
            sum1=sum1+dssh(i)*dssh(i)
            if(dssh(i)>rmax) rmax=dssh(i)
            if(dssh(i)<rmin) rmin=dssh(i)
      end do
      rmean=rmean/n
      rms=dsqrt(sum1/n)
      std=dsqrt(rms**2-rmean**2)
      write(*,*) '---------the total crossover discrepency-------------'
      write(6,500) n
500   format(/,7x,'参与比较点数:',I10,' (Input from File)')
      write(6,501) rmax,rmin,rmean,rms,std
501    format(/,7x,'最大差值：',f10.4,//
     &,7x,'最小差值：',f10.4,//,7x,'平均值：',f10.4,//
     &,7x,'rms：',f10.4,//,7x,'std：',f10.4)
      
      write(*,*) '-----------------------------------------------'  
      num=1
      dssh1=0
     
      sumssh=0
      rmax1=0
      rmin1=0
      do i1=1,count
         if (dabs(dssh(i1)-rmean)<3*STD) then
             kssh(num)=dssh(i1)
             dssh1=dssh1+kssh(num)
             sumssh=sumssh+kssh(num)*kssh(num)
             if(dssh(i1)>rmax1) rmax1=dssh(i1)
             if(dssh(i1)<rmin1) rmin1=dssh(i1)
             num=num+1
         endif
      enddo

      rmean=dssh1/num
      rms=sqrt(sumssh/num)
      std=sqrt(sumssh/num-rmean*rmean)
       !write(*,*) n,num
      write(*,*)'---the crossover discrepency after removes 3STD-------'
      write(6,600) num
600    format(/,7x,'参与比较点数:',I10,' (Input from File)')
      write(6,601) rmax1,rmin1,rmean,rms,std
601   format(/,7x,'最大差值：',f10.4,//,7x,'最小差值：'
     &,f10.4,//,7x,'平均值：',f10.4,//,7x,'rms：',f10.4
     &,//,7x,'std：',f10.4)
      write(*,*)'---------------------------------------'
      end 
      