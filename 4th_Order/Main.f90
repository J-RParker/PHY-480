program main
implicit none

integer :: i,j,n,bs,ssc
character(len=120) :: filename
double precision, allocatable, dimension(:,:) :: init !index, mass/x/y/z/vx/vy/vz
double precision, allocatable, dimension(:) :: m
double precision, allocatable, dimension(:,:) :: r, rp, vp ,err !index, x/y/z)
double precision, allocatable, dimension(:,:,:) :: v, a, v_interp, a_interp !index, x/y/z , times
double precision, dimension(1:3) :: dr !xyz
double precision :: t_init, t_run, dt, dt_max, t, simt, t_count_pos, t_log_pos !Time stuff
double precision :: kE, pE, E0, En, dE !Energy stuff
double precision :: small, relErr !Error stuff
double precision :: r2
double precision :: G, AU !Constants

G=6.67e-11
AU=149597870700.
t=0
dt=100
dt_max=dt*2**8
simt=1e10

t_count_pos=0
t_log_pos=1e5

ssc=0
small=1e-5
relErr=5e-6

open(10,file='init.txt',status='old')
read(10,*) n
write(6,*) n

allocate(init(1:n,1:7))
allocate(m(1:n))
allocate(r(1:n,1:3), rp(1:n,1:3), vp(1:n,1:3))
allocate(err(1:n,1:6))
allocate(v(1:n,1:3,-6:1), a(1:n,1:3,-6:1))
allocate(v_interp(1:n,1:3,-3:0), a_interp(1:n,1:3,-3:0))

read(10,*) init

r=0; v=0; a=0

m=init(:,1)
r=init(:,2:4)
v(:,:,-6)=init(:,5:7)

kE=0; pE=0; E0=0; En=0; dE=0; dr=0

deallocate(init)
close(10)
!allocate(bodies(1:n))
do i=1,n
    write(filename,'(A,I0,A)') 'Data/0',i,'_pos.txt'
    open(i+10,file=filename,status='unknown')
end do
open(9,file='Data/dt.txt',status='unknown')

!Calculate initial energy

do i =1,n
    kE=kE+0.5*m(i)*(v(i,1,-6)**2 + v(i,2,-6)**2 + v(i,3,-6)**2)
end do

do i=1,n-1
    do j=i+1,n
        dr(1:3) = r(j,1:3)-r(i,1:3)
        r2=dsqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
        pE=pE-G*m(i)*m(j)/r2
    end do
end do

E0 = kE+pE
call cpu_time(t_init)
write(6,*) E0
write(6,'(A,F6.3,A)') 't_int =', t_init,'s'





do i=1,n
    do j=1,n
        if (i==j) cycle
        dr(:)=r(j,:)-r(i,:)
        r2=dsqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
        a(i,:,-6) = a(i,:,-6) + G*m(j)*dr(:)/(r2**3)
    end do
end do

do bs= -6,-1
    r(:,:) = r(:,:) + v(:,:,bs)*dt + 0.5*a(:,:,bs)*dt**2
    do i=1,n
        do j=1,n
            if (i==j) cycle
            dr(:)=r(j,:)-r(i,:)
            r2=dsqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
            a(i,:,bs+1) = a(i,:,bs+1) + G*m(j)*dr(:)/(r2**3)
        end do
    end do
    v(:,:,bs+1) = v(:,:,bs) + 0.5*(a(:,:,bs+1) + a(:,:,bs))*dt  
end do
!---------------------------------------------------------------
do
    v(:,:,1)=0; a(:,:,1)=0
    rp=0; vp=0  
    
    !Predict r and v
    rp(:,:)=r(:,:) + (dt/24.)*((-9*v(:,:,-3))+(37*v(:,:,-2))+(-59*v(:,:,-1))+(55*v(:,:,0)))
    vp(:,:)=v(:,:,0) + (dt/24.)*((-9*a(:,:,-3))+(37*a(:,:,-2))+(-59*a(:,:,-1))+(55*a(:,:,0)))
    
    ! Calc new acc    
    do i=1,n
        do j=1,n
            if (i==j) cycle
            dr(:)=r(j,:)-r(i,:)
            r2=dsqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
            a(i,:,1) = a(i,:,1) + G*m(j)*dr(:)/(r2**3)
        end do
    end do

    !Correct r and v
    r(:,:)=r(:,:) + (dt/24.)*( (v(:,:,-2))+(-5*v(:,:,-1))+(19*v(:,:,0))+(9*vp(:,:)) )
    v(:,:,1)=v(:,:,0) + (dt/24.)*( (a(:,:,-2))+(-5*a(:,:,-1))+(19*a(:,:,0))+(9*a(:,:,1)) )
    
    t=t+dt
    t_count_pos=t_count_pos+dt
    
    if (t_count_pos>t_log_pos) then
        do i=1,n
            write(i+10,*) r(i,:)/AU
        end do
        write(9,*) dt
        t_count_pos=0
    end if
    
    if (t>simt) exit
    
    !Shuffle v and a arrays    
    do i=-5,1
        v(:,:,i-1)=v(:,:,i)
        a(:,:,i-1)=a(:,:,i)
    end do
    
    !Increment steps since change   
    ssc=ssc+1
    
    !Calc errors between p and c
    err(:,1:3)= (19./270.) * abs(r(:,:)-rp(:,:)) / (abs(r(:,:)) + small)
    err(:,4:6)= (19./270.) * abs(v(:,:,0)-vp(:,:)) / (abs(v(:,:,0)) + small)
    
    if (any(err > relErr)) then
        !Halve dt
        dt=0.5*dt
        v_interp(:,:,0)=v(:,:,0)
        v_interp(:,:,-1)=( (-5*v(:,:,-4)) + (28*v(:,:,-3)) + (-70*v(:,:,-2)) + (140*v(:,:,-1)) + (35*v(:,:,0)) )/128
        v_interp(:,:,-2)=v(:,:,-1)
        v_interp(:,:,-3)=( (3*v(:,:,-4)) + (-20*v(:,:,-3)) + (90*v(:,:,-2)) + (60*v(:,:,-1)) + (-5*v(:,:,0)) )/128
        a_interp(:,:,0)=a(:,:,0)
        a_interp(:,:,-1)=( (-5*a(:,:,-4)) + (28*a(:,:,-3)) + (-70*a(:,:,-2)) + (140*a(:,:,-1)) + (35*a(:,:,0)) )/128
        a_interp(:,:,-2)=a(:,:,-1)
        a_interp(:,:,-3)=( (3*a(:,:,-4)) + (-20*a(:,:,-3)) + (90*a(:,:,-2)) + (60*a(:,:,-1)) + (-5*a(:,:,0)) )/128
        
        v(:,:,-3:0)=v_interp(:,:,:)
        a(:,:,-3:0)=a_interp(:,:,:)
        ssc=0
    end if
    if (all(err < relErr/100.) .and. ssc>=3 .and. dt /= dt_max) then
        !Double dt (set ssd=0)
        dt=2*dt
        do i=1,3
            v(:,:,-i)=v(:,:,-2*i)
            a(:,:,-i)=a(:,:,-2*i)
        end do
        ssc=0
    end if
    
    
end do

call cpu_time(t_run)
write(6,'(A,F6.3,A)') 't_run =', t_run,'s'


end program