program lab12
 parameter (n=1024)
 real*8 :: start,finish
 integer :: time1,time2,block_size,n_min
 real init 
 real*8 a(n,n),b(n,n),c4(n,n)
 call init_random_seed()
 do i=1,n
  do j=1,n
   call random_number(init)
   a(i,j)=init
   call random_number(init)
   b(i,j)=init
   c4(i,j)=0
  end do
 end do
 
 print *,'блочный-алгоритм:'
 print *,n
 block_size=2
 do while (block_size<=1024)
  print *,'размер блоков=',block_size
  CALL cpu_time(start)
  call block(a,b,c4,n,block_size)
  CALL cpu_time(finish)
  print *,'Время работы = ',(finish-start),' с'
  block_size=block_size*2
 end do
end

subroutine timems(t,time)
 parameter(n=8)
 integer t(n),time
 time=((t(5)*60+t(6))*60+t(7))*10**3+t(8)
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

subroutine block(a,b,c,n,block_size)
 real*8 a(n,n),b(n,n),c(n,n)
 integer block_kol,block_size
 block_kol=n/block_size
 do jb=1,block_kol
  do j=(jb-1)*block_size+1,jb*block_size
   do kb=1,block_kol
    do k=(kb-1)*block_size+1,kb*block_size
     do ib=1,block_kol
      do i=(ib-1)*block_size+1,ib*block_size
       c(i,j)=c(i,j)+a(i,k)*b(k,j)
	  end do
     end do
    end do
   end do
  end do
 end do
end 
