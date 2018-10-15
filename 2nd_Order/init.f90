program init
implicit none

integer :: i,n
double precision :: m(1:9), Theta(1:9), vc(1:9), ax(1:9)
double precision :: r(1:9,1:3), v(1:9,1:3) 
double precision :: CoM(1:3), CoV(1:3)
double precision pi, G, AU, GM, Mtot

n=9

G=6.67259e-11
AU=149597870700.
GM=1.32712440018e20

Theta=0
pi=4.*atan(1.)

!Sun
!mass = 1·9884e30
!Planet 	a (AU)      	e           	i (deg) 	Omega(deg)	~omega(deg)	L (deg) 	Mass (kg)
! Mercury	0.38709893  	0.20563069  	7.00487 	48.33167	77.45645	252.25084	3.30114e23
! Venus 	0.72333199  	0.00677323  	3.39471 	76.68069	131.53298	181.97973	4.86747e24
! Earth 	1.00000011  	0.01671022  	0.00005 	-11.26064	102.94719	100.46435	5.97237e24
! Mars  	1.52366231  	0.09341233  	1.85061 	49.57854	336.04084	355.45332	6.41712e23
! Jupiter	5.20336301  	0.04839266  	1.30530 	100.55615	14.75385	34.40438	1.898187e27
! Saturn	9.53707032  	0.05415060  	2.48446 	113.71504	92.43194	49.94432	5.68336e26
! Uranus	19.19126393 	0.04716771  	0.76986 	74.22988	170.96424	313.23218	8.68127e25
! Neptune	30.06896348 	0.00858587  	1.76917 	131.72169	44.97135	304.88003	1.024126e26

r=0; v=0

ax(1)=0
ax(2)=0.38709893*AU
ax(3)=0.72333199*AU
ax(4)=1.00000011*AU
ax(5)=1.52366231*AU
ax(6)=5.20336301*AU
ax(7)=9.53707032*AU
ax(8)=19.19126393*AU
ax(9)=30.06896348*AU

m(1)=1.9884e30
m(2)=3.30114e23
m(3)=4.86747e24
m(4)=5.97237e24
m(5)=6.41712e23
m(6)=1.898187e27
m(7)=5.68336e26
m(8)=8.68127e25
m(9)=1.024126e26

do i=2,n
    Theta(i)=(i-2)*(2*pi)/(n-1)
end do
do i=2,n
    vc(i)=sqrt(GM/ax(i))
end do

do i=2,n
    r(i,1)=ax(i)*cos(theta(i))
    r(i,2)=ax(i)*sin(theta(i))
    v(i,1)=-vc(i)*sin(theta(i))
    v(i,2)=vc(i)*cos(theta(i))
end do

CoM=0
CoV=0
Mtot=0

do i=1,n
    CoM(1:3) = CoM(1:3) + m(i)*r(i,1:3)
    CoV(1:3) = CoV(1:3) + m(i)*v(i,1:3)
    Mtot=Mtot+m(i)
end do

CoM=CoM/Mtot
CoV=CoV/Mtot

do i=1,n
    r(i,:)=r(i,:)-CoM(:)
    v(i,:)=v(i,:)-CoV(:)
end do

open(10,file='init.txt',status='unknown')
write(10,*) m(1:n)
do i=1,3
    write(10,*) r(1:n,i)
end do
do i=1,3
    write(10,*) v(1:n,i)
end do
close(10)

end program

























