PROGRAM test
IMPLICIT NONE
INTEGER :: i
DOUBLE PRECISION :: n(1000000),nn(1000000),init_t

open (10,file='nums.txt', status='old')
open (11,file='f_out.txt',status='unknown')
read(10,*) n

nn=n**n

write(11,*) nn
close(10)
close(11)
call cpu_time(init_t)
write(6,*) init_t
end program
