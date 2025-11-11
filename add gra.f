       program xgridadd
!
! DESCRIPTION
! 
! 
! -A: input .grd3 file
! -B: input .grd3 file
! -G: output .grd3 file
! 
      implicit none
      integer*4 i,j,nx,ny,dim,type,nargs,iargc,par,ia
      integer*4 nn1(2),nn2(2)
      parameter(dim=34*60+1)
      real*4 indata1(dim,dim),indata2(dim,dim),outdata(dim,dim)
      real*8 x1,x2,y1,y2,dx,dy,lon,lat
      real*8 head1(6),head2(6)
      character*80 ifile1,ifile2,ofile1,ofile2,tbuf,cha
      logical io(10)
! Define constants
      data par/3/
! Get command-line arguments
!      nargs=iargc()
!      if(nargs.eq.0) then
!      call write_error
!      stop
!      end if
!      do i=1,par
!      io(i)=.true.
!      end do
      i=0
!      do ia=1,nargs
!      call getarg(ia,cha)
!	        if(cha(1:1).eq.'-') then
!! Necessary parametes
!                if(cha(2:2).eq.'A' .or. cha(2:2).eq.'a') then
!                ifile1=cha(3:)
! 		        i=i+1
!                io(i)=.false.
!                elseif(cha(2:2).eq.'B' .or. cha(2:2).eq.'b') then
!                ifile2=cha(3:)
! 		        i=i+1
!                io(i)=.false.
!                elseif(cha(2:2).eq.'G' .or. cha(2:2).eq.'g') then
!		        ofile1=cha(3:)
! 		        i=i+1
!                io(i)=.false.
!                else
!                call write_error
!                stop
!                end if
!      else
!      call write_error
!      stop
!      end if

!      end do
!! check necessary parameters
!      do i=1,par
!      if(io(i)) then
!      call write_error
!      stop
!      end if
!      end do
      ifile1="dg_nofilter_B3_L12.grd3"
      ifile2="dg_inzone_B3_L12.grd3"
      ofile1="ddg_final_B3_L12.grd3"
      ofile2="ddg_final_B3_L12.dat"
! open input files
      write(0,*) ifile1,ifile2
      open(10,file=ifile1,form='unformatted')
      open(11,file=ifile2,form='unformatted')
! Open output file
      open(60,file=ofile1,form='unformatted') ! signal
      open(61,file=ofile2)
! Read header
      read(10)(nn1(i),i=1,2),(head1(i),i=1,6)
      read(11)(nn2(i),i=1,2),(head2(i),i=1,6)
! check numbers of grids
      do i=1,2
      if(nn1(i).ne.nn2(i)) then
      write(0,*)'Inconsistent no. of grids. nn1, nn2=', nn1(i),nn2(i)
      stop
      end if
      end do
! check borders
      do i=1,6
      if(head1(i).ne.head2(i)) then
      write(0,*)'Inconsistent data borders. head1, head2=', 
     &         head1(i),head2(i)
      stop
      end if
      end do
      nx=nn1(1)
      ny=nn1(2)
      x1=head1(1)
      x2=head1(2)
      y1=head1(3)
      y2=head1(4)
      dx=head1(5)
      dy=head1(6)
! chose the data window for the given point
      	if(nx.gt.dim)stop'increase dim'
	if(ny.gt.dim)stop'increase dim'
      write(0,*)'GRD3 info:',nx,ny,x1,x2,y1,y2,dx,dy
      write(60)nx,ny,x1,x2,y1,y2,dx,dy
!Read data
      read(10)((indata1(i,j),J=1,NY),I=1,NX)
      read(11)((indata2(i,j),J=1,NY),I=1,NX)
      
      do i=1,nx
         lon=x1+(i-1)*dx
         do j=1,ny           
           lat=y1+(j-1)*dy
           outdata(i,j)=indata1(i,j)+indata2(i,j)
           write(61,*) lon,lat,outdata(i,j)
         end do
      
      end do
      write(60)((outdata(i,j),j=1,ny),i=1,nx)
      close(60)
      close(61)
      end       
      
      
           
      subroutine write_error
      write(0,*)'xgridadd -Afile1.grd3 -Bfile2.grd3 -Gout.grd3'
      return
      end   
  
