program main
  use mpi

  implicit none

  integer, dimension(1000000) :: temp1, temp2
  integer ierr, rank, numrank, i
  integer stat(MPI_STATUS_SIZE)

  double precision :: time1, time2, dur,gdur

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numrank, ierr)
  
  temp1=4
  temp2=5
  do i=0, 1000000, 100000 
    if(rank==0) then
      time1=MPI_Wtime()
      call MPI_SEND(temp1(:i), i, MPI_INTEGER, 1, 0, MPI_COMM_WORLD, ierr)
      call MPI_RECV(temp2(:i), i, MPI_INTEGER, 1, 0, MPI_COMM_WORLD, stat, ierr) 
      time2=MPI_Wtime()
      print *, sizeof(temp1(:i)), time2-time1
    else if(rank==1) then
      !time1=MPI_Wtime()
      call MPI_RECV(temp2(:i), i, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, stat ,ierr)
      time2=MPI_Wtime()
      call MPI_SEND(temp1(:i), i, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
    end if
  end do

  !if(rank==0) then
  !  print *, rank, time2-time1
  !end if

  call MPI_FINALIZE(ierr)
end program main
