program lab11
 real*8 :: start,finish
 integer :: time1,time2,n_min,n,error
 real :: init 
 real*8,allocatable :: a(:,:),b(:,:),c0(:,:),c1(:,:),c2(:,:),c3(:,:)
 n=2
 do while (n<=1024)
  print *,'n=',n
  call init_random_seed()
  allocate(a(n,n),b(n,n),c0(n,n),c1(n,n),c2(n,n),c3(n,n),stat=error)
  do i=1,n
   do j=1,n
    call random_number(init)
    a(i,j)=init
    call random_number(init)
    b(i,j)=init
    c0(i,j)=0
    c1(i,j)=0
    c2(i,j)=0
    c3(i,j)=0
   end do
  end do
  
  CALL cpu_time(start)
  c0 = matmul(a, b)
  CALL cpu_time(finish)
  print *,'стандартный-алгоритм:'
  print *,'Время работы = ',(finish-start),' с'
  
  CALL cpu_time(start)
  call jik(a,b,c1,n)
  CALL cpu_time(finish)
  print *,'jik-алгоритм:'
  print *,'Время работы = ',(finish-start),' с'
  
  CALL cpu_time(start)
  call jki(a,b,c2,n)
  CALL cpu_time(finish)
  print *,'jki-алгоритм:'
  print *,'Время работы = ',(finish-start),' с'
  
  CALL cpu_time(start)
  call strass(a,b,c3,n)
  CALL cpu_time(finish)
  print *,'Штрассен-алгоритм:'
  print *,'Время работы = ',(finish-start),' с'
  n=n*2
  deallocate(a,b,c0,c1,c2,c3)
 end do
end

subroutine init_random_seed()
    integer :: i, n, clock
    integer, DIMENSION(:), allocatable :: seed
    CALL RANDOM_SEED(size = n)
    allocate(seed(n))
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)
    deallocate(seed)
end

subroutine jik(a,b,c,n)
 real*8 a(n,n),b(n,n),c(n,n)
 do j=1,n
  do i=1,n
   do k=1,n
    c(i,j)=c(i,j)+a(i,k)*b(k,j)
   end do
  end do
 end do
end 

subroutine jki(a,b,c,n)
 real*8 a(n,n),b(n,n),c(n,n)
 do j=1,n
  do k=1,n
   do i=1,n
    c(i,j)=c(i,j)+a(i,k)*b(k,j)
   end do
  end do
 end do
end 

RECURSIVE subroutine strass(a,b,c,n)
 real*8 a(n,n),b(n,n),c(n,n)
 integer n_min,error
 real*8,allocatable :: P1(:,:),P2(:,:),P3(:,:),P4(:,:),P5(:,:),P6(:,:),P7(:,:)
 if (n.le.64) then 
  c = matmul(a, b)
 else
  m=n/2
  allocate(P1(m,m),P2(m,m),P3(m,m),P4(m,m),P5(m,m),P6(m,m),P7(m,m),stat=error)
  call strass(a(1:m,1:m)+a(m+1:n,m+1:n),b(1:m,1:m)+b(m+1:n,m+1:n),P1,m)
  call strass(a(m+1:n,1:m)+a(m+1:n,m+1:n),b(1:m,1:m),P2,m)
  call strass(a(1:m,1:m),b(1:m,m+1:n)-b(m+1:n,m+1:n),P3,m)
  call strass(a(m+1:n,m+1:n),b(m+1:n,1:m)-b(1:m,1:m),P4,m)
  call strass(a(1:m,1:m)+a(1:m,m+1:n),b(m+1:n,m+1:n),P5,m)
  call strass(a(m+1:n,1:m)-a(1:m,1:m),b(1:m,1:m)+b(1:m,m+1:n),P6,m)
  call strass(a(1:m,m+1:n)-a(m+1:n,m+1:n),b(m+1:n,1:m)+b(m+1:n,m+1:n),P7,m)
  c(1:m,1:m)=P1+P4-P5+P7
  c(1:m,m+1:n)=P3+P5
  c(m+1:n,1:m)=P2+P4
  c(m+1:n,m+1:n)=P1+P3-P2+P6
  deallocate(P1,P2,P3,P4,P5,P6,P7)
 end if
end 
