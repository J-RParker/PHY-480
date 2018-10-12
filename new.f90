program new
implicit none

integer :: i,j,n
double precision :: init(1:3,1:7)
double precision :: r(1:3,1:3), v(1:3,1:3), a(1:3,1:3), a0(1:3,1:3)
double precision :: dr(1:3), m(1:3) CoM(1:3), CoV(1:3) 
double precision :: dt, t, t_count_pos,t_count_en, simt, t_init, t_run, t_log_pos, t_log_en
double precision :: G, AU, r2, Mtot, kE, pE, E0, En, dE
n=3
dt=1000
simt=1e10
t_log_pos=1e6
t_log_en=1e7

G=6.67e-11
AU=149597870700.
t=0
t_count_pos=0
t_count_en=0

open(10,file='init.txt',status='old')
read(10,*) init
close(10)

m=init(:,1)
r=init(:,2:4)
v=init(:,5:7)
a0=0
a=0

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

kE=0
pE=0
E0=0
En=0
dE=0

do i =1,n
	kE=kE+0.5*m(i)*(v(i,1)**2 + v(i,2)**2 + v(i,3)**2)
end do

do i=1,n-1
	do j=i,n
		dr(1:3) = r(j,1:3)-r(i,1:3)
        r2=sqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
		pE=pE-G*m(i)*m(j)/r2
	end do
end do

E0 = kE=pE


call cpu_time(t_init)
write(6,'(A,F6.3,A)') 't_int =', t_init,'s'


open(10,file='dE.txt',status='unknown')
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
		En=kE+pE
		dE=(E0-En)/E0
		write(10,*) dE
		t_count_en=0
	end if
    
    if (t>simt) exit
end do
do i=0,n
    close(i+10)
end do
call cpu_time(t_run)
write(6,'(A,F6.3,A)') 't_run =', t_run,'s'
end program