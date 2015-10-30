!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Cubic spline interpolation
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012) -- Ch.1
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Aug. 19, 2015
!-----------------------------------------------------------------------------!

program spline_test
implicit none
integer::j,nd,np,io
real*8 ::x,f,fexact,xp,fp,pi
real*8,allocatable::xd(:),fd(:)


!---------------------------------------------!
!Given data set on equally spaced grid points
!---------------------------------------------!
nd = 10
allocate(xd(0:nd),fd(0:nd))
do j=0,nd
xd(j) =-1.0d0 + dfloat(j)*2.0d0/(dfloat(nd))
fd(j) = 1.0d0/(1.0d0+25.0d0*xd(j)**2)
end do

!Clustered towards to ends
!Chebyshev-Lobatto points (cosine distribution)
io=0   !set io=1 for clustered points
if (io.eq.1) then
	pi = 4.0d0*datan(1.0d0)
	do j=0,nd
	xd(j) =-dcos((j)*pi/dfloat(nd))
	fd(j) = 1.0d0/(1.0d0+25.0d0*xd(j)**2)
	end do
end if


! Writing data to text file
open(11, file="dataset.plt")
write(11,*)'variables ="x","f"'
do j=0,nd
write(11,*) xd(j),fd(j)
end do
close(11)


! Writing exact solution using np points
np = 2000 !use 2000 points to plot curve
open(12, file="exact_curve.plt")
write(12,*)'variables ="x","f"'
	do j=0,np
		xp =-1.0d0 + dfloat(j)*2.0d0/(dfloat(np))
		fp = 1.0d0/(1.0d0+25.0d0*xp**2)
		write(12,*) xp,fp
	end do
close(12)

!compare solutions at a desired point
x = 0.7d0
fexact = 1.0d0/(1.0d0+25.0d0*x**2)
call spline(nd,xd,fd,x,f)
write(*,100)"numerical:",f
write(*,100)"exact    :",fexact

100 format (A20,F12.4)

! Writing numerical solution using np points 
open(13, file="numerical_curve.plt")
write(13,*)'variables ="x","f"'
	do j=0,np
		xp =-1.0d0 + dfloat(j)*2.0d0/(dfloat(np))
		call spline(nd,xd,fd,xp,f)
		write(13,*) xp,f
	end do
close(13)

       
end


!-----------------------------------------------------------------------------!
!Construct cubic spline interpolation from a given data set fd(xd) 
!fd(i): given data values 
!xd(i): given points
!     where i=0,1,2,...nd 
!     nd: maximum number of data point (total number of points is nd+1)
!
!
!f: interpolated value
!x: interpolation location 
!
!note: use natural spline boundary conditions
!-----------------------------------------------------------------------------!
subroutine spline(nd,xd,fd,x,f)
implicit none
integer::nd,i
real*8 ::xd(0:nd),fd(0:nd),g(0:nd)
real*8 ::x,f,d
real*8, dimension(:),allocatable::a,b,c,r,z

!compute second derivatives at data points
g(0) = 0.0d0 !natural spline b.c.
g(nd)= 0.0d0 !natural spline b.c.

!construct tridiagonal system:
allocate(a(nd-1),b(nd-1),c(nd-1),r(nd-1),z(nd-1))
!---------------------------------------------!
!Given system
!---------------------------------------------!
do i=1,nd-1
a(i) = (xd(i)-xd(i-1))/6.0d0
b(i) = (xd(i+1)-xd(i-1))/3.0d0
c(i) = (xd(i+1)-xd(i))/6.0d0 
r(i) = (fd(i+1)-fd(i))/(xd(i+1)-xd(i)) &
     - (fd(i)-fd(i-1))/(xd(i)-xd(i-1))
end do
!following two lines can be deactivated for natural spline b.c.
!since g(0)=0 and g(nd)=0
r(1)    = r(1) - a(1)*g(0)        !b.c.
r(nd-1) = r(nd-1) - c(nd-1)*g(nd) !b.c.

!Thomas Algorithm
call tdma(a,b,c,r,z,1,nd-1)

!second derivatives
do i=1,nd-1
g(i) = z(i)
end do  
 
!construct the cubic spline polynomials inside the domain
do i=0,nd-1
  
	if (x.ge.xd(i).and.x.le.xd(i+1)) then
    d = xd(i+1)-xd(i) 
    
	f = g(i)/6.0d0*(((xd(i+1)-x)**3)/d - d*(xd(i+1)-x)) &
      + g(i+1)/6.0d0*(((x-xd(i))**3)/d - d*(x-xd(i))) &
      + fd(i)*(xd(i+1)-x)/d + fd(i+1)*(x-xd(i))/d 
    
	end if
end do

! outside of the domain
if (x.lt.xd(0)) then  !lower off-bound
i=0
    
	d = xd(i+1)-xd(i)

	f = g(i)/6.0d0*(((xd(i+1)-x)**3)/d - d*(xd(i+1)-x)) &
      + g(i+1)/6.0d0*(((x-xd(i))**3)/d - d*(x-xd(i))) &
      + fd(i)*(xd(i+1)-x)/d + fd(i+1)*(x-xd(i))/d 

else if (x.gt.xd(nd)) then !upper off-bound
i=nd-1

	d = xd(i+1)-xd(i)

	f = g(i)/6.0d0*(((xd(i+1)-x)**3)/d - d*(xd(i+1)-x)) &
      + g(i+1)/6.0d0*(((x-xd(i))**3)/d - d*(x-xd(i))) &
      + fd(i)*(xd(i+1)-x)/d + fd(i+1)*(x-xd(i))/d
   
end if

    

end



!------------------------------------------------------------------!
!Tridiagonal matrix algorithm (TDMA)
!Thomas algorithm
!solution tridiagonal systems
!a: lower diagonal
!b: main diagonal
!c: upper diagonal
!r: source vector
!x: solution vector
!   for indices s(start) to e(end)
!   i: s,s+1,s+2, ....,i,....,e 
!
!Note: a(s) and c(e) are dummy coefficients, not used.
!------------------------------------------------------------------!

subroutine tdma(a,b,c,r,x,s,e)
implicit none
integer s,e,i
real*8, dimension(s:e) ::a,b,c,r,x    

! forward elimination phase
do i=s+1,e
b(i) = b(i) - a(i)/b(i-1)*c(i-1)
r(i) = r(i) - a(i)/b(i-1)*r(i-1)
end do
! backward substitution phase 
x(e) = r(e)/b(e)
do i=e-1,s,-1
x(i) = (r(i)-c(i)*x(i+1))/b(i)
end do

return
end
