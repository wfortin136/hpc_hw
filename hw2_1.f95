program main
  use mpi

  implicit none

  !include 'mpif.h'
  !Structure to store complex numbers. Decided to use this
  !instead of using the COMPLEX library
  type cplex
    sequence
    real*8 :: real
    real*8 :: imag
    integer*8 :: mand
  end type cplex

  character(10) :: prog_name, arg1_str, arg2_str, arg3_str
  character(50) :: filename
  integer arg1, arg2, i, j, x
  integer nx, ny, nrows_l, lstart, lend, data_size
  type(cplex), dimension(:), allocatable :: static_l_s, dyn_l_s
  integer, dimension(:), allocatable :: status_array
  integer ierr, rank, numrank, last_rank
  logical flag
  integer stat(MPI_STATUS_SIZE)

  double precision :: time1, time2, dur,gdur

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numrank, ierr)
  
  !Get command line arguments. This is kerazy compared to c!!
  call getarg(0, prog_name)
  call getarg(1, arg1_str)
  call getarg(2, arg2_str)
  call getarg(3, arg3_str)

  filename = prog_name(:len_trim(prog_name)-2)//"_"//trim(arg1_str)//"x"//trim(arg2_str)//".out"
  open(99, file=trim(filename))
  
  read (arg1_str, *) arg1
  read (arg2_str, *) arg2
  nx=arg1
  ny=arg2

  if(arg3_str == "static") then
    nrows_l = (nx/numrank) !

    allocate(static_l_s((nrows_l)*(ny)))

    lstart = rank*nrows_l
    lend = lstart + nrows_l
    time1= MPI_Wtime()
    call calc_mandel(lstart, lend, nx, ny, static_l_s)
    call static_run(static_l_s, size(static_l_s,1), numrank, rank)
    time2= MPI_Wtime()
    deallocate(static_l_s)
  end if
  
  if(arg3_str == "dynamic") then
    nrows_l = 1
    allocate(dyn_l_s((nrows_l)*(ny)))
    time1= MPI_Wtime()
    call dynamic_run(dyn_l_s, size(dyn_l_s,1), numrank, rank)
    time2= MPI_Wtime()
    deallocate(dyn_l_s)
  end if
  
  dur=time2-time1
  print *, dur
  call MPI_Reduce(dur, gdur, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  
  if(rank==0) then
    print *, arg3_str, numrank, dur
  end if
  call MPI_FINALIZE(ierr)

contains

  integer FUNCTION calc_pixel(c)
    implicit none
    type(cplex), intent(in) :: c
    integer :: counter

    type(cplex) :: z
    integer max_it
    real*8 temp, len_sqrd

    len_sqrd=0
    max_it =256
    counter=0
    z%real = 0
    z%imag=0

    do while((len_sqrd <4) .AND. (counter<max_it))
      !real and imaginary parts of Zn+1 where Zn+1=Zn^2+c
      temp = c%real * c%real - c%imag * c%imag + c%real
      z%imag = 2*c%real*c%imag +c%imag
      z%real=temp
      len_sqrd = (z%real*z%real) +(z%imag * z%imag)
      counter=counter+1
    end do
    if(counter==max_it) then
      calc_pixel=1
    else
      calc_pixel=0
    end if
  end function calc_pixel


  subroutine calc_mandel(start, fin, xcount,ycount, mand_struc)
    integer, intent(in) :: start, fin, xcount, ycount
    type(cplex), intent(inout) :: mand_struc(*)
    
    type(cplex) c, mand_ar
    integer x, y, indx
    !count across values of x.
    !for each value of x, iterate through values of y
    !use the equation x/nx*4 - 2 to normalize
    !array between bounds of -2 to 2. Be definition
    !we know anything >|2| will approach infinity
    
    indx=0
    do x=start, fin-1
      
      c%real = x/(xcount*1.) * 4. - 2.
      do y=0, (ny-1)
        c%imag = y/(ycount*1.)*4. - 2.
        indx=indx+1
        mand_struc(indx)%imag=c%imag
        mand_struc(indx)%real=c%real
        mand_struc(indx)%mand=calc_pixel(c)
      enddo
    enddo
  end subroutine calc_mandel
  
  subroutine static_run(data_set, data_size, nums, r)
    type(cplex), intent(inout) :: data_set(*)
    integer, intent(in) :: nums, r, data_size

    integer i,j
    integer stat(MPI_STATUS_SIZE)
  
    if(r==0) then
      !print *, size(data_l_s,1) 
      do j=1, data_size
        !print *, data_set(j)
        write(99,*) data_set(j)
      end do

      do i=1, nums-1
        !instead of creating a new MPI_data_type struct, I cheated and used an 8byte integer
        !and two 8 byte reals for the array size
        call MPI_RECV(data_set, 3*data_size, MPI_REAL8, i, 0, MPI_COMM_WORLD, stat ,ierr)
        do j=1, data_size
          !print *, data_set(j)
          write(99,*) data_set(j)
        end do
      end do
    
    else
      call MPI_SEND(data_set, 3*data_size, MPI_REAL8, 0, 0, MPI_COMM_WORLD, ierr)
    end if

  end subroutine static_run

  subroutine dynamic_run(data_set, data_size, nums, r)
    type(cplex), intent(inout) :: data_set(*)
    integer, intent(in) :: nums, r, data_size

    integer i,j, start, fin
    integer stat(MPI_STATUS_SIZE)
  
    if(r==0) then
    
    x=0
    do i=1, nums-1
      call MPI_SEND(x, 1, MPI_INTEGER, i, 0, MPI_COMM_WORLD, ierr)
      x=x+1
    end do

    do while(x<nx)
          call MPI_RECV(data_set, 3*data_size, MPI_REAL8, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, stat, ierr)
          call MPI_SEND(x, 1, MPI_INTEGER, stat(MPI_SOURCE), 0, MPI_COMM_WORLD, ierr)
          x=x+1
          do j=1, data_size
            !print *, data_set(j)
            write(99,*) data_set(j)
          end do
    end do

    do i=1, nums-1
      call MPI_RECV(data_set, 3*data_size, MPI_REAL8, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, stat, ierr)
      do j=1, data_size
        !print *, data_set(j)
        write(99,*) data_set(j)
      end do
    end do 

    do i=1, nums-1
      call MPI_SEND(0, 1, MPI_INTEGER, i, 1, MPI_COMM_WORLD, ierr)
    end do
  else
    do while(.true.)
        call MPI_RECV(start, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, stat ,ierr)
        
        if(stat(MPI_TAG) == 1) then
          exit
        end if
        fin = start+nrows_l
        call calc_mandel(start, fin, nx, ny, data_set)
        call MPI_SEND(data_set, 3*data_size, MPI_REAL8, 0, 0, MPI_COMM_WORLD, ierr)
    end do
  end if
  end subroutine dynamic_run
end program main
