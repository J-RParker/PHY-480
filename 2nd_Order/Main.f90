program Main
implicit none

integer :: i,j,n
double precision, Allocatable :: init(:,:) !(body index, mass/xyz/vx vy vz data)
double precision, Allocatable :: r(:,:), v(:,:), a(:,:), a0(:,:) !(body index, xyz data)
double precision, Allocatable :: m(:) !(body index)
double precision :: dr(1:3) !(xyz data)
double precision :: dt, t, t_count_pos, t_count_en, simt, t_init, t_run, t_log_pos, t_log_en
double precision :: G, AU, r2, kE, pE, E0, Et, dE

n=9

dt=1000
simt=1e10
t_log_pos=1e5
t_log_en=1e7

G=6.67e-11
AU=149597870700.
t=0
t_count_pos=0
t_count_en=0

Allocate(r(1:n,1:3), v(1:n,1:3), a(1:n,1:3), a0(1:n,1:3))
Allocate(m(1:n))
Allocate(init(1:n,1:7))

open(10,file='init.txt',status='old')
read(10,*) init
close(10)

m=init(:,1)
r=init(:,2:4)
v=init(:,5:7)

deallocate(init)

a0=0; a=0; kE=0; pE=0; E0=0; En=0; dE=0

do i =1,n
    kE=kE+0.5*m(i)*(v(i,1)**2 + v(i,2)**2 + v(i,3)**2)
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
write(6,'(A,F6.3,A)') 't_int =', t_init,'s'

open(10,file='dE.txt',status='unknown')

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
    t_count_pos=t_count_pos+dt
    t_count_en=t_count_en+dt
    if (t_count_pos>t_log_pos) then
        do i=1,n
            write(i+10,*) r(i,1:3)/AU
        end do
        t_count_pos=0
    end if
    if (t_count_en>t_log_en) then
        kE=0
        pE=0
        do i =1,n
            kE=kE+0.5*m(i)*(v(i,1)**2 + v(i,2)**2 + v(i,3)**2)
        end do
        do i=1,n-1
            do j=i+1,n
                dr(1:3) = r(j,1:3)-r(i,1:3)
                r2=sqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
                pE=pE-G*m(i)*m(j)/r2
            end do
        end do
        Et=kE+pE
        dE=(E0-Et)/E0
        write(10,*) dE
        t_count_en=0
    end if
    
    if (t>simt) exit
end do
close(10)
deallocate(r,v,a,a0,m)
call cpu_time(t_run)
write(6,'(A,F6.3,A)') 't_run =', t_run,'s'
end program