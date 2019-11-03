program parametric_curves
use types
implicit none

real(dp) :: wa=3._dp, xa, ya
real(dp) :: wb=2._dp**.5_dp, xb, yb
real(dp) :: wc=2._dp/3._dp, xc, yc
real(dp) :: wd=pi, xd, yd
real(dp) :: we=5._dp/3._dp, xe, ye
real(dp) :: wf=.5_dp*(1._dp+5._dp**.5_dp), xf, yf
real(dp) :: t, dt
integer :: fua=11,fub=22,fuc=33,fud=44,fue=55,fuf=66
integer :: i, N

rewind(fua)
rewind(fub)
rewind(fuc)
rewind(fud)
rewind(fue)
rewind(fuf)

open(unit=fua,file='parametric_curves_a.dat')
open(unit=fub,file='parametric_curves_b.dat')
open(unit=fuc,file='parametric_curves_c.dat')
open(unit=fud,file='parametric_curves_d.dat')
open(unit=fue,file='parametric_curves_e.dat')
open(unit=fuf,file='parametric_curves_f.dat') 

print*, "Number of points, N, for parametric_curves plots:"
read(*,*) N

dt = (40._dp*pi)/N


do i=1,N
t = -20._dp*pi + dt*(i-1)
call func_para(t,wa,xa,ya)
write(fua,*) xa, ya
call func_para(t,wb,xb,yb)
write(fub,*) xb, yb
call func_para(t,wc,xc,yc)	
write(fuc,*) xc, yc
call func_para(t,wd,xd,yd)	
write(fud,*) xd, yd
call func_para(t,we,xe,ye)
write(fue,*) xe, ye
call func_para(t,wf,xf,yf)	
write(fuf,*) xf, yf		
enddo
close(fua)
close(fub)
close(fuc)
close(fud)
close(fue)
close(fuf)

contains

subroutine func_para(t,w,x,y)
    implicit none
    real(dp),intent(in) :: t,w
    real(dp),intent(out) :: x,y
    y = sin(t*w)
    x = sin(t)
end subroutine func_para
end program parametric_curves