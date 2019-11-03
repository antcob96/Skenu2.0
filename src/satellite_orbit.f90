program satellite_orbit
use types
implicit none

real(dp) :: dx, a_o, v_o, a_n, v_n, b, e_n
real(dp) :: yu, yb, x, delta_a, xb, ap

real(dp) :: pe, e_o, delta ,G, mass_e, mu
real(dp) :: time, v_0, v_f

integer :: i, N, file_unit = 10,file_unit2 = 11, k, l
integer :: file_unit3 = 12

open(unit=file_unit,file='Orbit.dat')
open(unit=file_unit2,file='Ratio_semi_maj_min.dat')
open(unit=file_unit3,file='Velocity.dat')


print*, "Number of points for satellite orbit graph, N:"
read*, N 



call parameters_sat(pe, e_o, delta, G, mass_e)

mu = G*mass_e
! Initialize 
a_o = pe/(1._dp-e_o)

ap = (1._dp+e_o)/(1._dp-e_o)*pe

a_n = pe/(1._dp-e_o)



v_o = sqrt(((G*mass_e)/pe)*(1._dp+e_o))

k = 1

dx = 2._dp*a_o/float(N)

delta_a = a_n - a_o
!print*, delta_a
 v_0 = sqrt(((G*mass_e)/ap)*(1._dp-e_o))
	! Write velocity at r_max in a data file
    write(file_unit3,*) k, v_0


time = 2._dp*pi*sqrt((a_o**3)/mu)

b = sqrt(pe*ap)


write(file_unit2,*) k, b/a_o

	! Plot orbit
	do i=1,N+1
	x = a_o - dx*float(i-1)
	yu = (b/a_o)*sqrt(a_o**2-(x-delta_a)**2)

	write(file_unit,*) x, yu
	enddo

	do i=1,N+1
	x = -a_o + dx*float(i-1)
	yb = -(b/a_o)*sqrt(a_o**2-(x-delta_a)**2)

	write(file_unit,*) x, yb
	enddo






do while (e_o>1._dp/7._dp)

v_n = (v_o-v_o*delta)




	! Find new eccentricity
e_n = sqrt(1._dp+((pe*v_n/mu)**2)*(v_n**2-(2._dp*mu/pe)))



	! Find semi-major axis value
	a_o = pe/(1._dp - e_n)


	
	! Find the nubmer of orbits it takes to be a circular orbit
	k = k+1



dx = abs((a_n-a_o))/float(N)

time = time + 2._dp*pi*sqrt((a_o**3)/mu)


!print*, delta_a

ap = pe*(1._dp+e_n)/(1._dp-e_n)
b = sqrt(ap*pe)

delta_a = (a_n-a_o)

print*, delta_a, a_n, a_o

write(file_unit2,*) k, b/a_o

	! Plot orbit
	do i=1,N+1
	x = a_o - dx*float(i-1)
	yu = (b/a_o)*sqrt(a_o**2-(x-delta_a)**2)
	
	write(file_unit,*) x, yu
	enddo


	! Plot orbit
	do i=1,N+1
	x = -a_o + dx*float(i-1)
	yb = -(b/a_o)*sqrt(a_o**2-(x-delta_a)**2)

	write(file_unit,*) x, yb
	enddo

	v_o = v_n

	v_f = sqrt(((G*mass_e)/ap)*(1._dp-e_n))

	! Write velocity in a data file

		write(file_unit3,*) k, v_f


	e_o = e_n
enddo

v_f = v_f

time = ((time/60._dp)/60._dp)/24._dp

print*, "Number of orbits it takes to circularize:"
print*, k

print*, "Time for the satellite to reach a circular orbit in [days]:"
print*, time 

print*, "Velocity of 1st orbit at r_max"
print*, v_0

print*, "Velocity of last orbit at r_max"
print*, v_f


print*, "Difference between the 1st orbit and last orbit velocity:"
print*, v_f-v_0

contains

subroutine parameters_sat(perihelion, eccentricity_init, delta_size, Gravity, mass_of_planet)
    implicit none
    real(dp), intent(out) :: perihelion, eccentricity_init, delta_size, Gravity, mass_of_planet
    perihelion = 6.6_dp*10._dp**6 ! Perihelion distance Mm
    eccentricity_init = 0.9656_dp ! Inital eccentricity
    delta_size = 0.01_dp ! Fraction reduction of velocity
    mass_of_planet = 5.98_dp*10._dp**(24)! Kg
    Gravity = 6.67_dp*10._dp**(-11) ! Gravitational cosntant
end subroutine parameters_sat



end program satellite_orbit