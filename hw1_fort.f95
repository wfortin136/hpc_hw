program main
implicit none

!Structure to store complex numbers. Decided to use this
!instead of using the COMPLEX library
type cplex
  sequence
  real*8 :: real
  real*8 :: imag
end type cplex

  real*8, dimension(:), allocatable :: data_l, data_tmp
  integer c_count, nx, ny, i, j
  type(cplex) = c

  allocate(data_l(nx*ny))

  do i=0, 6
    c%real = i/(nx*1.) * 4. - 2. 
    do j=0, 6
      c%imag = j/(ny*1.)*4. - 2.
      c_count = calc_pixel(c)
      *data_l+=1
      *data_l=c_count
  
      print *, i, j, c_count
    enddo
  enddo


  deallocate(data_l)
!end program main

contains

pure FUNCTION calc_pixel(type(complex) c)
  type(cplex), intent(in) :: c
  integer :: counter

  type(cplex) = z
  integer max_it
  real*8 temp, len_squrd=0

  max_iter =500
  counter=0
  c%real = 0
  c%imag=0

  do while((len_sqrd <4) && (counter<max_it))
    !real and imaginary parts of Zn+1 where Zn+1=Zn^2+c
    temp = c%real * c%real - c%imag * c%imag + c%real
    z%imag = 2*c%reak*c%imag +c%imag
    z%real=temp
    len_sqrd = (z%real*z%real) +(z%imag * z%imag)
    counter++
  end do
RETURN
END





!integer :: x, y, i, j, k, N=2048, N2=1024,blocklen=4, nn
!double precision, DIMENSION(:,:), ALLOCATABLE :: A, B, C
!real*8 mat_size, time, flops
!real rand_1, rand_2
!integer start, finish, rate, flop
!
!allocate(A(N,N))
!allocate(B(N,N))
!allocate(C(N,N))
!
!mat_size=sizeof(A)*3/1024.
!
!write(6, 5) "Matricies: 3 ", N, "x", N, ", Total Size: ", mat_size, "KB"
!5 format(A13,I4,A1,I4,1X, A13, F10.3, 1x, A2) 
!
!do x=1, N
!  do y=1, N
!    call RANDOM_NUMBER(rand_1)
!    call RANDOM_NUMBER(rand_2)
!    A(x,y)=rand_1*10
!    B(x,y)=rand_2*10
!    !print *, x, y, A(x,y), B(x,y)
!  enddo
!enddo
!
!print *, "***Naive***"
!do nn=1, 11
!  call naive_mmx(A,B,C, 2**nn)
!  !call printC(C,N,"Naive")
!enddo
!
!print *, "***Tiling***"
!do nn=2, 11
!  call tiling_mmx(A,B,C, 2**nn, blocklen)
!  !call printC(C,N,"Tiling")
!enddo
!
!C(:,:)=0.
!print *, "***Recursion***"
!do nn=2, 11
!  call recursion_mmx(A,B,C, 2**nn, blocklen)
!  !call printC(C,N,"Recursion")
!enddo
!
!print *, "***DGEMM***"
!do nn=2, 11
!  call blas_mmx(A,B,C, 2**nn)
!enddo
!
!print *, "***Tiling***"
!do nn=2, 10
!  call tiling_mmx(A,B,C, N2, 2**nn)
!  !call printC(C,N,"Tiling")
!enddo
!
!print *, "***Recursion***"
!do nn=2, 10
!  call recursion_mmx(A,B,C, N2, 2**nn)
!  !call printC(C,N,"Tiling")
!enddo
!
!deallocate (A)
!deallocate (B)
!deallocate (C)
!!end program main
!
!contains
!
!  SUBROUTINE calc_pixel(TYPE(complex) 
!
!SUBROUTINE calc_speed(start, finish, rate, n, flop, time, flops)
!  integer, intent(in)::n
!  integer*8, intent(in) :: start, finish, rate
!  integer*8, intent(out) ::flop
!  real*8, intent(out) :: time, flops 
!
!  time = (finish-start)/(1.* rate)
!  flop=2*n**3
!  flops=(flop*1.)/time
!
!END SUBROUTINE calc_speed
!
!SUBROUTINE printstats(flop, time, flops, n)
!  real*8, intent(in) :: time, flops
!  integer*8, intent(in) ::  flop
!  integer, intent(in) :: n
!  !write(6, 15) "FLOPs: ", flop
!  !15 format(A7,I12)
!  write(6, 25) "Time: ", time, "sec", "N", n
!  25 format(A6,E10.1E2, 1x,A3,1x, A1, I5)
!  !write(6, 35) "FLOPs/sec: ", flops 
!  !35 format(A11,E10.1E2)
!
!END SUBROUTINE printstats
!
!SUBROUTINE naive_mmx(A,B,C, n)
!  implicit none
!  double precision, intent(in) :: A(:,:), B(:,:)
!  double precision, intent(out) :: C(:,:)
!  integer, intent(in) :: n
!
!  integer :: i,j,k
!  integer*8 :: start, finish, rate, flop
!  double precision, dimension(:), ALLOCATABLE :: rowA, colB
!  double precision :: cval
!  real*8 :: time, flops
!  allocate(rowA(n))
!  allocate(colB(n))
!
!  call system_clock(COUNT_RATE=rate)
!  call system_clock(COUNT=start)
!
!  do i=1, n
!    rowA=A(i,1:n)
!    do j=1, n
!      cval=C(i,j)
!      colB=B(1:n,j)
!      do k=1, n
!        cval= cval+(rowA(k)*colB(k))
!        !print *, C(i,j)
!      enddo
!      C(i,j)=cval
!    enddo
!  enddo
!
!  call system_clock(COUNT=finish)
!
!  call calc_speed(start,finish, rate, n, flop, time, flops)
!
!  call printstats(flop, time, flops, n)
!
!  deallocate(rowA)
!  deallocate(colB)
!END SUBROUTINE naive_mmx
!
!SUBROUTINE tiling_mmx(A,B,C, n,bs)
!  implicit none
!  double precision, intent(in) :: A(:,:), B(:,:)
!  double precision, intent(out) :: C(:,:)
!  integer, intent(in) :: n, bs
!
!  integer*8 :: start, finish, rate, flop
!  real*8 :: time, flops, mat_size
!  double precision, dimension(:,:), ALLOCATABLE :: subA, subB, subC
!  integer :: i,j,k,blocks, x1,y1, x2, y2, a3, b3
!
!  C(:,:)=0.
!  if(MOD(n,bs) .NE. 0) then
!    print *, "ERROR: Block size does not evenly divide matricies"
!  END IF
!  blocks = n/bs
!
!  allocate(subA(bs,bs))
!  allocate(subB(bs,bs))
!  allocate(subC(bs,bs))
!
!  mat_size=sizeof(subA)*3/1024.
!
!!  write(6, 5) "Block Matrix: 3 ", bs, "x", bs, ", Total Size: ", mat_size, "KB"
!!    5 format(A16,I4,A1,I4,1X, A13, F10.3, 1x, A2)
!
!  call system_clock(COUNT_RATE=rate)
!  call system_clock(COUNT=start)
!
!  do i=1, blocks
!    do j=1, blocks
!      x1=bs*(i-1)+1
!      x2=bs*i
!      y1=bs*(j-1)+1
!      y2=bs*j
!      subC=C(x1:x2,y1:y2)
!      !print *,x1,x2,y1,y2
!      do k=1, blocks
!	a3=bs*(k-1)+1
!	b3=bs*k
!	subA=A(x1:x2,a3:b3)
!	subB=B(a3:b3,y1:y2)
!	!print *,a3,b3
!	subC= subC+MATMUL(subA,subB)
!	!print *, C(i,j)
!      enddo
!      C(x1:x2,y1:y2)=subC
!    enddo
!  enddo
!
!  call system_clock(COUNT=finish)
!
!  call calc_speed(start,finish, rate, n, flop, time, flops)
!
!  call printstats(flop, time, flops, n)
!
!  deallocate(subA)
!  deallocate(subB)
!  deallocate(subC)
!
!END SUBROUTINE tiling_mmx
!
!RECURSIVE SUBROUTINE RMM(A,B,C,n,bs)
!  double precision, intent(in) :: A(:,:), B(:,:)
!  double precision, intent(out) :: C(:,:)
!  integer, intent(in) :: n, bs
!
!  integer x1,y1,x2,y2
! !print *, n/2
!  if( n==bs) THEN 
!    C=C+MATMUL(A,B)
!!    print *, n
!  else
!    x1=1
!    y1=bs/2
!    x2=bs/2+1
!    y2=bs
!  
!    call RMM(A(x1:y1,x1:y1),B(x1:y1,x1:y1),C(x1:y1,x1:y1),&
!      n/2,bs)
!    call RMM(A(x1:y1,x2:y2),B(x2:y2,x1:y1),C(x1:y1,x1:y1),&
!      n/2,bs)
!
!    call RMM(A(x1:y1,x1:y1),B(x1:y1,x2:y2),C(x1:y1,x2:y2),&
!      n/2,bs)
!    call RMM(A(x1:y1,x2:y2),B(x2:y2,x2:y2),C(x1:y1,x2:y2),&
!      n/2,bs)
!
!    call RMM(A(x2:y2,x1:y1),B(x1:y1,x1:y1),C(x2:y2,x1:y1),&
!      n/2,bs)
!    call RMM(A(x2:y2,x2:y2),B(x2:y2,x1:y1),C(x2:y2,x1:y1),&
!      n/2,bs)
!
!    call RMM(A(x2:y2,x1:y1),B(x1:y1,x2:y2),C(x2:y2,x2:y2),&
!      n/2,bs)
!    call RMM(A(x2:y2,x2:y2),B(x2:y2,x2:y2),C(x2:y2,x2:y2),&
!      n/2,bs)
!  end if
!
!END SUBROuTINE RMM
!
!SUBROUTINE recursion_mmx(A,B,C,n,bs)
!  double precision, intent(in) :: A(:,:), B(:,:)
!  double precision, intent(out) :: C(:,:)
!  integer, intent(in) :: n, bs
!
!  integer*8 :: start, finish, rate, flop
!  real*8 :: time, flops
!
!  call system_clock(COUNT_RATE=rate)
!  call system_clock(COUNT=start)
!
!  call RMM(A,B,C, n, bs)
!  
!  call system_clock(COUNT=finish)
!  call calc_speed(start,finish, rate, n, flop, time, flops)
!  call printstats(flop, time, flops,n)
!
!END SUBROUTINE recursion_mmx
!
!subroutine printC(C,n,mattype)
!  integer, intent(in)::n
!  double precision, intent(in) :: C(:,:)
!  character(LEN=*) :: mattype
!
!  integer x,y
!
!  print *, "***********"
!  print *, mattype
!
!  do x=1, n
!    do y=1, n
!      print *, C(x,y)
!    enddo
!  enddo 
!end subroutine printC
!
!SUBROUTINE blas_mmx(A,B,C, N)
!  double precision, intent(in) :: A(:,:), B(:,:)
!  double precision, intent(out) :: C(:,:)
!  integer, intent(in) :: N
!  
!  integer*8 :: start, finish, rate, flop
!  real*8 :: time, flops
!
!  call system_clock(COUNT_RATE=rate)
!  call system_clock(COUNT=start)
!
!  call DGEMM('n','n', N, N, N, 1, A, N, B, N, 1, C, N)
!  
!  call system_clock(COUNT=finish)
!  call calc_speed(start,finish, rate, n, flop, time, flops)
!  call printstats(flop, time, flops,n)
!end subroutine
!
!end program main
