module MpiWrapper

  use ConstantsModule, only: LENMEMPATH, LENVARNAME
  use KindModule, only: DP, I4B

  implicit none

#ifdef MPI_PARALLEL
  include 'mpif.h'
#endif

  private

  integer, parameter :: mpiwrplblk_size = 64
#ifdef MPI_PARALLEL
  character(len=MPI_MAX_PROCESSOR_NAME) :: mpiwrpmyname
  integer, dimension(MPI_STATUS_SIZE,mpiwrplblk_size) :: mpiwrpstats
#endif
  integer :: mpiwrpcomm_world, mpiwrpcomm_null
  integer :: mpiwrpnrproc, mpiwrpmyrank
  integer, dimension(1000) :: rreq, sreq
  integer :: lenbuf

  ! -- public functions
  public :: mpiwrpcomm_size
  public :: mpiwrpcomm_rank
  public :: mpiwrpinit
  public :: mpiwrpbarrier
  public :: mpiwrpfinalize
  public :: mpiwrpcommworld
  public :: mpiwrpallgather
  public :: mpiwrpallgatherv
  public :: mpiwrpreduce
  public :: mpiwrpallreduce
  public :: mpiwrpallreduceloc
  public :: mpiwrpcommgroup
  public :: mpiwrpgroupincl
  public :: mpiwrpcommcreate
  public :: mpiwrpgrouprank
  public :: mpiwrpstats
  public :: mpiwrpwaitall
  public :: mpiwrpprobe
  public :: mpiwrptypefree
  public :: mpiwrpgetcount
  ! -- public variables
  public :: mpiwrpnrproc
  public :: mpiwrpmyrank

  save

  interface mpiwrpreduce
    module procedure mpiwrpreduced
  end interface

  interface mpiwrpallreduceloc
    module procedure mpiwrpallreducedloc
  end interface

  interface mpiwrpallreduce
    module procedure mpiwrpallreducei, mpiwrpallreducer, mpiwrpallreduced
  end interface

  interface mpiwrpallgather
    module procedure mpiwrpallgatheri
  end interface

  interface mpiwrpallgatherv
    module procedure mpiwrpallgathervi, mpiwrpallgathervd
  end interface

  contains

  integer function mpiwrpcomm_size(comm)
! ******************************************************************************
! Wrapper subroutine around MPI_COMM_SIZE te determine number of
! processes.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in) :: comm ! (I) communicator
    ! -- local
    integer :: ierr, size
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call MPI_Comm_size(comm, size, ierr)
    mpiwrpcomm_size = size
#else
    mpiwrpcomm_size = 1
#endif
    ! -- return
    return
  end function mpiwrpcomm_size

  integer function mpiwrpcomm_rank(comm)
! ******************************************************************************
! Wrapper function around MPI_COMM_RANK to determine own
! process rank ID.
! @return rank ID of my own process.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in) :: comm ! (I) communicator
    ! -- local
    integer :: ierr, rank
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
   if (comm /= MPI_COMM_NULL) then
     call MPI_Comm_rank(comm, rank, ierr)
     mpiwrpcomm_rank = rank
   else
     mpiwrpcomm_rank = -1
   endif
#else
   mpiwrpcomm_rank = 0
#endif
    ! -- return
    return
  end function mpiwrpcomm_rank

   subroutine mpiwrpinit()
! ******************************************************************************
! Wrapper subroutine around MPI_INIT for initialization.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    ! -- local
    integer :: ierr, i, required, provided
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    required = MPI_THREAD_FUNNELED
    call MPI_Init_thread(required, provided, ierr)
    if (required /= provided) then
      if (mpiwrpmyrank == 0) then
        write(*,*) 'Warning, could not guarantee thread safety.'
      end if
    end if
#endif
    ! -- return
    return
  end subroutine mpiwrpinit

  subroutine mpiwrpfinalize()
! ******************************************************************************
! Wrapper subroutine around MPI_FINALIZE termination.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    ! -- local
    integer :: ierr
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call MPI_Finalize(ierr)
#endif
    ! -- return
    return
  end subroutine mpiwrpfinalize

  subroutine mpiwrpbarrier(comm)
! ******************************************************************************
! Wrapper subroutine around MPI_BARRIER barrier.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in) :: comm
    ! -- local
    integer :: ierr
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call MPI_Barrier(comm, ierr)
#endif
    ! -- return
    return
  end subroutine mpiwrpbarrier

  integer function mpiwrpcommworld()
! ******************************************************************************
! Wrapper subroutine around MPI_BARRIER barrier.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    ! -- local
    integer :: ierr
! ------------------------------------------------------------------------------

#ifdef MPI_PARALLEL
    mpiwrpcommworld = MPI_COMM_WORLD
#else
    mpiwrpcommworld = 0
#endif
    ! -- return
    return
  end function mpiwrpcommworld

  subroutine mpiwrpwaitall(count, reqs, status)
! ******************************************************************************
! Waits for a given set of communication requests to complete.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    integer, intent(in)                                  :: count   ! (I)   number of request handles
    integer, dimension(count), intent(inout)             :: reqs   ! (I/O) array with request handles
    integer, dimension(MPI_STATUS_SIZE,*), intent(inout) :: status ! (I) status
    ! -- local
    integer :: ierr, i, n
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    if ( count <= mpiwrplblk_size ) then
       call MPI_Waitall(count, reqs, status, ierr)
       return
    end if

    do i = 1, count, mpiwrplblk_size
       n = min(mpiwrplblk_size, count - i + 1)
       call MPI_Waitall(n, reqs(i), status, ierr)
    end do
#endif
    ! -- return
    return
  end subroutine mpiwrpwaitall

  subroutine mpiwrpprobe(source, tag, comm, status)
! ******************************************************************************
! Waits for a given set of communication requests to complete.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    integer, intent(in)                                  :: source ! (I) source rank
    integer, intent(in)                                  :: tag    ! (I) tag
    integer, intent(in)                                  :: comm   ! (I) communicator
    integer, dimension(MPI_STATUS_SIZE,*), intent(inout) :: status ! (I/O) status
    ! -- local
    integer :: ierr
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
     call MPI_Probe(source, tag, comm, status, ierr)
#else
     call mpiwrperror(comm, 'mpiwrpprobe', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpprobe

  subroutine mpiwrpgetcount(status, datatype, count)
! ******************************************************************************
! Waits for a given set of communication requests to complete.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    integer, dimension(MPI_STATUS_SIZE), intent(in) :: status   ! (I) status
    integer, intent(in)                             :: datatype ! (I) datatype
    integer, intent(out)                            :: count    ! (O) count
    ! -- local
    integer :: ierr
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
     call MPI_Get_count(status, datatype, count, ierr)
#else
     call mpiwrperror(comm, 'mpiwrpgetcount', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpgetcount

 subroutine mpiwrpreduced(sendbuf, count, op, root, comm)
! ******************************************************************************
! Combines values (doubles) from all processes and distributes the result back
! to root processes.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in)                               :: count   ! Number of elements in send buffer
    double precision, dimension(count), intent(inout) :: sendbuf ! Send buffer
    character(len=*)                                  :: op      ! Operation
    integer, intent(in)                               :: root    ! Rank of root process
    integer, intent(in)                               :: comm    ! Communicator
    ! -- local
    integer  :: ierr, i
    double precision, dimension(count) :: recvbuf ! receive buffer.
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    if (comm == MPI_COMM_NULL) then
      return
    endif
    select case(op)
      case('mpi_min')
        call MPI_Reduce(sendbuf, recvbuf, count, mpi_double_precision, mpi_min, &
                           root, comm, ierr)
      case('mpi_max')
        call MPI_Reduce(sendbuf, recvbuf, count, mpi_double_precision, mpi_max, &
                           root, comm, ierr)
      case('mpi_sum')
        call MPI_Reduce(sendbuf, recvbuf, count, mpi_double_precision, mpi_sum, &
                           root, comm, ierr)
      case default
        call mpiwrperror(comm, 'mpiwrpreduced', 'invalid operation')
    end select
    !
    do i = 1, count
      sendbuf(i) = recvbuf(i)
    end do
    !
#else
    call mpiwrperror(comm, 'mpiwrpreduced', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpreduced

  subroutine mpiwrpallreducei(sendbuf, count, op, comm)
! ******************************************************************************
! Combines values (integers) from all processes and distributes the result back
! to all processes.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in)                   :: count   ! Number of elements in send buffer
    integer, dimension(count), intent(inout) :: sendbuf ! Send buffer
    character(len=*)                      :: op      ! Operation
    integer, intent(in)                   :: comm    ! Communicator
    ! -- local
    integer  :: ierr, i
    integer, dimension(count) :: recvbuf ! receive buffer.
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    if (comm == MPI_COMM_NULL) then
      return
    endif
    select case(op)
      case('mpi_min')
        call MPI_Allreduce(sendbuf, recvbuf, count, mpi_integer, mpi_min,       &
                           comm, ierr)
      case('mpi_max')
        call MPI_Allreduce(sendbuf, recvbuf, count, mpi_integer, mpi_max,       &
                           comm, ierr)
      case('mpi_sum')
        call MPI_Allreduce(sendbuf, recvbuf, count, mpi_integer, mpi_sum,       &
                           comm, ierr)
      case default
        call mpiwrperror(comm, 'mpiwrpallreducei', 'invalid operation')
    end select
    !
    do i = 1, count
      sendbuf(i) = recvbuf(i)
    end do
    !
#else
    call mpiwrperror(comm, 'mpiwrpallreducei', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpallreducei

  subroutine mpiwrpallreducer(sendbuf, count, op, comm)
! ******************************************************************************
! Combines values (reals) from all processes and distributes the result back
! to all processes.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in)                   :: count   ! Number of elements in send buffer
    real, dimension(count), intent(inout) :: sendbuf ! Send buffer
    character(len=*)                      :: op      ! Operation
    integer, intent(in)                   :: comm    ! Communicator
    ! -- local
    integer  :: ierr, i
    real, dimension(count) :: recvbuf ! receive buffer.
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    if (comm == MPI_COMM_NULL) then
      return
    endif
    select case(op)
      case('mpi_min')
        call MPI_Allreduce(sendbuf, recvbuf, count, mpi_real, mpi_min,          &
                           comm, ierr)
      case('mpi_max')
        call MPI_Allreduce(sendbuf, recvbuf, count, mpi_real, mpi_max,          &
                           comm, ierr)
      case('mpi_sum')
        call MPI_Allreduce(sendbuf, recvbuf, count, mpi_real, mpi_sum,          &
                           comm, ierr)
      case default
        call mpiwrperror(comm, 'mpiwrpallreducer', 'invalid operation')
    end select
    !
    do i = 1, count
      sendbuf(i) = recvbuf(i)
    end do
    !
#else
    call mpiwrperror(comm, 'mpiwrpallreducer', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpallreducer

  subroutine mpiwrpallreduced(sendbuf, count, op, comm)
! ******************************************************************************
! Combines values (doubles) from all processes and distributes the result back
! to all processes.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in)                               :: count   ! Number of elements in send buffer
    double precision, dimension(count), intent(inout) :: sendbuf ! Send buffer
    character(len=*)                                  :: op      ! Operation
    integer, intent(in)                               :: comm    ! Communicator
    ! -- local
    integer  :: ierr, i
    double precision, dimension(count) :: recvbuf ! receive buffer.
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    if (comm == MPI_COMM_NULL) then
      return
    endif
    select case(op)
      case('mpi_min')
        call MPI_Allreduce(sendbuf, recvbuf, count, mpi_double_precision,       &
                           mpi_min, comm, ierr)
      case('mpi_max')
        call MPI_Allreduce(sendbuf, recvbuf, count, mpi_double_precision,       &
                           mpi_max, comm, ierr)
      case('mpi_sum')
        call MPI_Allreduce(sendbuf, recvbuf, count, mpi_double_precision,       &
                           mpi_sum, comm, ierr)
      case default
        call mpiwrperror(comm, 'mpiwrpallreduced', 'invalid operation')
    end select
    !
    do i = 1, count
      sendbuf(i) = recvbuf(i)
    end do
    !
#else
    call mpiwrperror(comm, 'mpiwrpallreduced', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpallreduced

  subroutine mpiwrpallreducedloc(sendbuf, count, op, comm)
! ******************************************************************************
! Combines values (doubles) from all processes and distributes the result back
! to all processes.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in)                                 :: count   ! Number of elements in send buffer
    double precision, dimension(2,count), intent(inout) :: sendbuf ! Send buffer
    character(len=*)                                    :: op      ! Operation
    integer, intent(in)                                 :: comm    ! Communicator
    ! -- local
    integer  :: ierr, i
    double precision, dimension(2,count) :: recvbuf ! receive buffer.
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    if (comm == MPI_COMM_NULL) then
      return
    endif
    select case(op)
      case('mpi_min')
        call MPI_Allreduce(sendbuf, recvbuf, count, mpi_2double_precision,      &
                           mpi_minloc, comm, ierr)
      case('mpi_max')
        call MPI_Allreduce(sendbuf, recvbuf, count, mpi_2double_precision,      &
                           mpi_maxloc, comm, ierr)
      case default
        call mpiwrperror(comm, 'mpiwrpallreducedloc', 'invalid operation')
    end select
    !
    do i = 1, count
      sendbuf(1,i) = recvbuf(1,i)
      sendbuf(2,i) = recvbuf(2,i)
    end do
    !
#else
    call mpiwrperror(comm, 'mpiwrpallreduced', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpallreducedloc

  subroutine mpiwrperror(comm, subname, message)
! ******************************************************************************
! Prints an error message and aborts the current program.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in)           :: comm    ! (I) communicator
    character(len=*), intent(in)  :: subname ! (I) the name of the subroutine in which the error occured
    character(len=*), intent(in)  :: message ! (I) the error message to be printed
    ! -- local
    integer :: ierr, i
    integer, dimension(1) :: idummy(1)
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    write (6,*)
    write (6,*) '*** fatal runtime error in subroutine '//subname
    write (6,*)
    write (6,*) '    '//message
    write (6,*)
    write (6,*) '*** trying to terminate all processes...'
    write (6,*)
    !
    call mpi_abort(comm, 1, ierr)
#else
    write (6,*)
    write (6,*) '*** fatal runtime error in subroutine '//subname
    write (6,*)
    write (6,*) '    '//message
    write (6,*)
    write (6,*) '*** aborting program...'
    write (6,*)
    !
    ! --   try to force a core dump
    i         = -666666666
    idummy(i) =  666666666
    !
    stop
    !
#endif
    ! -- return
    return
  end subroutine mpiwrperror

  subroutine mpiwrpallgatheri(comm, gsbuf, gscnt, grbuf, grcnt)
! ******************************************************************************
! Gathers integer values into specified memory
! locations from a group of processes and broadcasts the
! gathered data to all processes.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer :: comm, gscnt
    integer :: grcnt
    integer, dimension(*) :: gsbuf, grbuf
    ! -- local
    integer :: ierr, i, k
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call MPI_Allgather(gsbuf, gscnt, mpi_integer,                               &
                       grbuf, grcnt, mpi_integer,                               &
                       comm, ierr)
#endif
    ! -- return
    return
  end subroutine mpiwrpallgatheri

  subroutine mpiwrpallgathervi(comm, gsbuf, gscnt, grbuf, grcnt, offsets)
! ******************************************************************************
! Gathers integer values into specified memory
! locations from a group of processes and broadcasts the
! gathered data to all processes.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer :: comm, gscnt
    integer, dimension(*) :: grcnt, offsets
    integer, dimension(*) :: gsbuf, grbuf
    ! -- local
    integer :: ierr, i, k
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call MPI_Allgatherv(gsbuf, gscnt, mpi_integer,                             &
                        grbuf, grcnt, offsets, mpi_integer,                    &
                        comm, ierr)
#else
    if ( grcnt(1) < gscnt ) then
      call mpiwrperror(comm, 'mpiwrpallgathervi',                              &
                        'receive buffer too small')
    endif
    !
    k = offsets(1)
    !
    do i = 1, gscnt
      grbuf(k+i) = gsbuf(i)
    enddo
#endif
    ! -- return
    return
  end subroutine mpiwrpallgathervi

  subroutine mpiwrpallgathervd(comm, gsbuf, gscnt, grbuf, grcnt, offsets)
! ******************************************************************************
! Gathers double values into specified memory
! locations from a group of processes and broadcasts the
! gathered data to all processes.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer :: comm, gscnt
    integer, dimension(*) :: grcnt, offsets
    double precision, dimension(*) :: gsbuf, grbuf
    ! -- local
    integer :: ierr, i, k
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call MPI_Allgatherv(gsbuf, gscnt, mpi_double_precision,                    &
                        grbuf, grcnt, offsets, mpi_double_precision,           &
                        comm, ierr)
#else
    if ( grcnt(1) < gscnt ) then
      call mpiwrperror(comm, 'mpiwrpallgathervd',                              &
                       'receive buffer too small')
    endif
    !
    k = offsets(1)
    !
    do i = 1, gscnt
      grbuf(k+i) = gsbuf(i)
    enddo
#endif
    ! -- return
    return
  end subroutine mpiwrpallgathervd

  subroutine mpiwrpcommgroup(comm, group)
! ******************************************************************************
! Accesses the group associated with given communicator
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in) :: comm
    integer, intent(out) :: group
    ! -- local
    integer :: ierr
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call MPI_Comm_group(comm, group, ierr)
#else
    call mpiwrperror(comm, 'mpiwrpcommgroup', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpcommgroup

  subroutine mpiwrpgroupincl(old_group, n, rnks, new_group)
! ******************************************************************************
! Produces a group by reordering an existing group and taking only listed
! members.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in) :: old_group
    integer, intent(in) :: n
    integer, dimension(n), intent(in) :: rnks
    integer, intent(out) :: new_group
    ! -- local
    integer :: ierr
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call MPI_Group_incl(old_group, n, rnks, new_group, ierr)
#else
    call mpiwrperror(comm, 'mpiwrpgroupincl', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpgroupincl

  subroutine mpiwrpcommcreate(old_comm, group, new_comm)
! ******************************************************************************
! Creates a new communicator.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in) :: old_comm
    integer, intent(in) :: group
    integer, intent(out) :: new_comm
    ! -- local
    integer :: ierr
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call MPI_Comm_create(old_comm, group, new_comm, ierr)
#else
    call mpiwrperror(comm, 'mpiwrpcommcreate', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpcommcreate

  subroutine mpiwrpgrouprank(group, rank)
! ******************************************************************************
! Returns the rank of this process in the given group.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in) :: group
    integer, intent(out) :: rank
    ! -- local
    integer :: ierr
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call mpi_group_rank(group, rank, ierr)
#else
    call mpiwrperror(comm, 'mpiwrpgrouprank', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpgrouprank

  subroutine mpiwrptypefree(newtype)
! ******************************************************************************
! Create a new MPI data type.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in) :: newtype
    ! -- local
    integer :: ierr
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL

   call MPI_Type_free(newtype, ierr)
#endif
    ! -- return
    return
  end subroutine mpiwrptypefree

end module MpiWrapper