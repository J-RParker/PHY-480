program main
implicit none

integer :: i,j,n
double precision, allocatable, dimension(:,:) :: init !index, mass/x/y/z/vx/vy/vz
double precision, allocatable, dimension(:) :: m
double precision, allocatable, dimension(:,:) :: r, a0 !index, x/y/z)
double precision, allocatable, dimension(:,:,:) :: v, a !index, x/y/z , times
double precision, dimension(1:3) :: dr !xyz
double precision :: t_init !Time stuff
double precision :: kE, pE, E0, En, dE !Energy stuff
double precision :: r2
double precision :: G !Constants

G=6.67e-11

open(10,file='init.txt',status='old')
read(10,*) n
write(6,*) n

allocate(init(1:n,1:7))
allocate(m(1:n))
allocate(r(1:n,1:3), a0(1:n,1:3))
allocate(v(1:n,1:3,1:8), a(1:n,1:3,1:8)) !Change later?

read(10,*) init

v=0; a=0

m=init(:,1)
r=init(:,2:4)
v(:,:,1)=init(:,5:7)

a=0; a0=0; kE=0; pE=0; E0=0; En=0; dE=0; dr=0

deallocate(init)

!Calculate initial energy

do i =1,n
    kE=kE+0.5*m(i)*(v(i,1,1)**2 + v(i,2,1)**2 + v(i,3,1)**2)
end do

do i=1,n-1
    do j=i+1,n
        dr(1:3) = r(j,1:3)-r(i,1:3)
        r2=sqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
        pE=pE-G*m(i)*m(j)/r2
    end do
end do

E0 = kE+pE
call cpu_time(t_init)
write(6,*) E0
write(6,'(A,F6.3,A)') 't_int =', t_init,'s'

end program