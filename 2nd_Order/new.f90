program new
implicit none

integer :: i,j,n
double precision :: r(1:3,1:3), v(1:3,1:3), a(1:3,1:3), a0(1:3,1:3), init(1:3,1:7)
double precision :: dr(1:3), m(1:3)
double precision :: dt, t, logt, G, AU, simt, r2, t_init, t_run, t_out

n=3
dt=1000
simt=1e10
t_out=1e6

G=6.67e-11
AU=149597870700.
t=0
logt=0

open(10,file='init.txt',status='old')
read(10,*) init
close(10)

m=init(:,1)
r=init(:,2:4)
v=init(:,5:7)
a0=0
a=0

call cpu_time(t_init)
write(6,'(A,F6.3,A)') 't_int =', t_init,'s'

open(11,file='sun.txt',status='unknown')
open(12,file='earth.txt',status='unknown')
open(13,file='jupiter.txt',status='unknown')
do i=1,n
    write(i+10,*) r(i,1:3)/AU
end do
do i=1,n
    do j=1,n
        if (i==j) cycle
        dr(1:3) = r(j,1:3)-r(i,1:3)
        r2=sqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
        a0(i,1:3) = a0(i,1:3) + G*m(j)*dr(1:3)/(r2**3)
    end do
end do
do
    r = r + v*dt + 0.5*a0*dt**2
    do i=1,n
        do j=1,n
            if (i==j) cycle
            dr(1:3) = r(j,1:3)-r(i,1:3)
            r2=sqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
            a(i,1:3) = a(i,1:3) + G*m(j)*dr(1:3)/(r2**3)
        end do
    end do
    v=v+0.5*(a+a0)*dt
    a0=a
    a=0
    t=t+dt
    logt=logt+dt
    if (logt>t_out) then
        do i=1,n
            write(i+10,*) r(i,1:3)/AU
        end do
        logt=0
    end if
    
    if (t>simt) exit
end do
do i=1,n
    close(i+10)
end do
call cpu_time(t_run)
write(6,'(A,F6.3,A)') 't_run =', t_run,'s'
end program