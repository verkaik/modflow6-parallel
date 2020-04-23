module MpiWrapper
  
  use ConstantsModule, only: LENORIGIN, LENVARNAME
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
  
  type :: MetaMemoryType
    character(len=LENVARNAME) :: name
    character(len=LENORIGIN)  :: origin
    integer(I4B)              :: memitype
    integer(I4B)              :: isize
    integer(I4B)              :: ncol = 0
    integer(I4B)              :: nrow = 0
    integer(I4B)              :: nlay = 0
  end type MetaMemoryType  
  
  type ColMemoryType
    character(len=LENVARNAME)  :: name        !name of the array
    character(len=LENORIGIN)   :: origin      !name of origin
    integer(I4B)               :: memitype    !integer type
    logical                    :: logicalsclr !logical
    integer(I4B)               :: intsclr     !integer
    real(DP)                   :: dblsclr     !double
    integer(I4B), dimension(3) :: aint1d      !1d integer array
  end type ColMemoryType
  
  ! -- public functions
  public :: mpiwrpcomm_size
  public :: mpiwrpcomm_rank
  public :: mpiwrpinit
  public :: mpiwrpbarrier
  public :: mpiwrpfinalize
  public :: mpiwrpcommworld
  public :: mpiwrpisend
  public :: mpiwrpirecv
  public :: mpiwrpallgather
  public :: mpiwrpallgatherv
  public :: mpiwrpreduce
  public :: mpiwrpallreduce
  public :: mpiwrpcommgroup
  public :: mpiwrpgroupincl
  public :: mpiwrpcommcreate
  public :: mpiwrpcolstruct
  public :: mpiwrpmmtstruct
  public :: mpiwrpmtstruct
  public :: mpiwrpgrouprank
  public :: mpiwrpstats
  public :: mpiwrpwaitall
  public :: mpiwrpprobe
  public :: mpiwrptypefree
  public :: mpiwrpgetcount
  ! -- public variables
  public :: mpiwrpnrproc
  public :: mpiwrpmyrank
  ! -- public types
  public :: ColMemoryType
  public :: MetaMemoryType
  
  save
  
  interface mpiwrpisend
    module procedure mpiwrpisendmmt, mpiwrpisendmt
  end interface

  interface mpiwrpirecv
    module procedure mpiwrpirecvmmt, mpiwrpirecvmt
  end interface
  
  interface mpiwrpreduce
    module procedure mpiwrpreduced
  end interface
  
  interface mpiwrpallreduce
    module procedure mpiwrpallreducei, mpiwrpallreducer, mpiwrpallreduced
  end interface
  
  interface mpiwrpallgather
    module procedure mpiwrpallgatheri
  end interface 
  
  interface mpiwrpallgatherv
    module procedure mpiwrpallgathervi, mpiwrpallgathervd,                      &
                     mpiwrpallgathervcol 
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

  subroutine mpiwrpisendmmt(buf, count, datatype, dest, tag, comm, req)
! ******************************************************************************
! Begins a non-blocking send of a meta memory type structure to the
! process with rankID dest.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in)                                :: count    ! number of integers to be sent
    type(MetaMemoryType), dimension(count), intent(in) :: buf      ! buffer to be sent
    integer, intent(in)                                :: datatype ! data type
    integer, intent(in)                                :: dest     ! rank id of destination process
    integer, intent(in)                                :: tag      ! message tag
    integer, intent(in)                                :: comm     ! communicator
    integer, intent(out)                               :: req      ! request handle
    ! -- local
    integer :: ierr      
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call MPI_Isend(buf, count, datatype, dest, tag, comm, req, ierr)
#else
    req = -1
    call mpiwrperror(comm, 'mpiwrpisendmmt', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpisendmmt

  subroutine mpiwrpisendmt(buf, count, datatype, dest, tag, comm, req)
! ******************************************************************************
! Begins a non-blocking send of a memory type structure to the
! process with rankID dest.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryTypeModule, only: MemoryType
    ! -- dummy
    integer, intent(in)                            :: count    ! number of objects to be sent
    type(MemoryType), dimension(count), intent(in) :: buf      ! buffer to be sent
    integer, intent(in)                            :: datatype ! data type
    integer, intent(in)                            :: dest     ! rank id of destination process
    integer, intent(in)                            :: tag      ! message tag
    integer, intent(in)                            :: comm     ! communicator
    integer, intent(out)                           :: req      ! request handle
    ! -- local
    integer :: ierr      
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call MPI_Isend(buf, count, datatype, dest, tag, comm, req, ierr)
#else
    req = -1
    call mpiwrperror(comm, 'mpiwrpisendmt', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpisendmt

  subroutine mpiwrpirecvmmt(buf, count, datatype, dest, tag, comm, req)
! ******************************************************************************
! Begins a non-blocking receive of a meta memory type structure from the
! process with rankID dest.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in)                                   :: count    ! number of objects to be received
    type(MetaMemoryType), dimension(count), intent(inout) :: buf      ! buffer to be received
    integer, intent(in)                                   :: datatype ! buffer data type
    integer, intent(in)                                   :: dest     ! rank id of destination process
    integer, intent(in)                                   :: tag      ! message tag
    integer, intent(in)                                   :: comm     ! communicator
    integer, intent(out)                                  :: req      ! request handle
    ! -- local
    integer :: ierr      
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call MPI_Irecv(buf, count, datatype, dest, tag, comm, req, ierr)
#else
    req = -1
    call mpiwrperror(comm, 'mpiwrpirecvmmt', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpirecvmmt
  
  subroutine mpiwrpirecvmt(buf, count, datatype, dest, tag, comm, req)
! ******************************************************************************
! Begins a non-blocking receive of a  memory type structure from the
! process with rankID dest.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryTypeModule, only: MemoryType
    ! -- dummy
    integer, intent(in)                               :: count    ! number of objects to be received
    type(MemoryType), dimension(count), intent(inout) :: buf      ! buffer to be received
    integer, intent(in)                               :: datatype ! buffer data type
    integer, intent(in)                               :: dest     ! rank id of destination process
    integer, intent(in)                               :: tag      ! message tag
    integer, intent(in)                               :: comm     ! communicator
    integer, intent(out)                              :: req      ! request handle
    ! -- local
    integer :: ierr      
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call MPI_Irecv(buf, count, datatype, dest, tag, comm, req, ierr)
#else
    req = -1
    call mpiwrperror(comm, 'mpiwrpirecvmt', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpirecvmt
  
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
        call MPI_Reduce(sendbuf, recvbuf, count, mpi_double, mpi_min,            &
                           root, comm, ierr)
      case('mpi_max')
        call MPI_Reduce(sendbuf, recvbuf, count, mpi_double, mpi_max,            &
                           root, comm, ierr)
      case('mpi_sum')
        call MPI_Reduce(sendbuf, recvbuf, count, mpi_double, mpi_sum,            &
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
        call MPI_Allreduce(sendbuf, recvbuf, count, mpi_double, mpi_min,          &
                           comm, ierr)
      case('mpi_max')
        call MPI_Allreduce(sendbuf, recvbuf, count, mpi_double, mpi_max,          &
                           comm, ierr)
      case('mpi_sum')
        call MPI_Allreduce(sendbuf, recvbuf, count, mpi_double, mpi_sum,          &
                           comm, ierr)
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

  subroutine mpiwrpallgathervcol(comm, gsbuf, gscnt, gstype, grbuf, grcnt,      &
                                 grtype, offsets)
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
    integer :: comm, gscnt, gstype, grtype
    integer, dimension(*) :: grcnt
    type(ColMemoryType), dimension(*) :: gsbuf
    type(ColMemoryType), dimension(*) :: grbuf
    integer, dimension(*) :: offsets

    ! -- local
    integer :: ierr, i, k
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call MPI_Allgatherv(gsbuf, gscnt, gstype, grbuf, grcnt, offsets, grtype,    &
                        comm, ierr)
    !
#endif
    ! -- return
    return
  end subroutine mpiwrpallgathervcol
  
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

  subroutine mpiwrpcolstruct(newtype)
! ******************************************************************************
! Create a new MPI data type.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ConstantsModule, only: LENORIGIN, LENVARNAME
    use KindModule, only: DP, I4B    
    ! -- dummy
    integer, intent(out) :: newtype
    ! -- local
    type(ColMemoryType) :: cmtdum
    
    integer, dimension(7) :: types, blocklengths
    integer(KIND=MPI_ADDRESS_KIND), dimension(7) :: displacements
    integer(KIND=MPI_ADDRESS_KIND) :: base
    
    integer :: i, ierr 
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    ! -- set-up derived MPI data type 
    types(1) = MPI_CHARACTER
    types(2) = MPI_CHARACTER
    types(3) = MPI_INTEGER
    types(4) = MPI_INTEGER
    types(5) = MPI_REAL
    types(6) = MPI_DOUBLE
    types(7) = MPI_INTEGER
    blocklengths(1)  = LENVARNAME
    blocklengths(2)  = LENORIGIN
    blocklengths(3)  = 1
    blocklengths(4)  = 1
    blocklengths(5)  = 1
    blocklengths(6)  = 1
    blocklengths(7)  = 3
    call MPI_GET_ADDRESS(cmtdum%name,       displacements(1), ierr)
    call MPI_GET_ADDRESS(cmtdum%origin,     displacements(2), ierr)
    call MPI_GET_ADDRESS(cmtdum%memitype,   displacements(3), ierr)
    call MPI_GET_ADDRESS(cmtdum%logicalsclr,displacements(4), ierr)
    call MPI_GET_ADDRESS(cmtdum%intsclr,    displacements(5), ierr)
    call MPI_GET_ADDRESS(cmtdum%dblsclr,    displacements(6), ierr)
    call MPI_GET_ADDRESS(cmtdum%aint1d,     displacements(7), ierr)
    base = displacements(1)
    do i = 1, 7
      displacements(i) = displacements(i) - base
    enddo
    ! -- create and commit datatype
    call MPI_TYPE_CREATE_STRUCT(7, blocklengths, displacements, types, newtype, ierr)
    call MPI_TYPE_COMMIT(newtype, ierr)
#else
    call mpiwrperror(comm, 'mpiwrpcolstruct', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpcolstruct

  subroutine mpiwrpmmtstruct(newtype)
! ******************************************************************************
! Create a new MPI data type.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ConstantsModule, only: LENORIGIN, LENVARNAME
    use KindModule, only: DP, I4B    
    ! -- dummy
    integer, intent(out) :: newtype
    ! -- local
    type(MetaMemoryType) :: mmtdum
    
    integer, dimension(7) :: types, blocklengths
    integer(KIND=MPI_ADDRESS_KIND), dimension(7) :: displacements
    integer(KIND=MPI_ADDRESS_KIND) :: base
    
    integer :: i, ierr 
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    ! -- set-up derived MPI data type 
    types(1) = MPI_CHARACTER
    types(2) = MPI_CHARACTER
    types(3) = MPI_INTEGER
    types(4) = MPI_INTEGER
    types(5) = MPI_INTEGER
    types(6) = MPI_INTEGER
    types(7) = MPI_INTEGER
    blocklengths(1)  = LENVARNAME
    blocklengths(2)  = LENORIGIN
    blocklengths(3)  = 1
    blocklengths(4)  = 1
    blocklengths(5)  = 1
    blocklengths(6)  = 1
    blocklengths(7)  = 1
    call MPI_GET_ADDRESS(mmtdum%name,     displacements(1), ierr)
    call MPI_GET_ADDRESS(mmtdum%origin,   displacements(2), ierr)
    call MPI_GET_ADDRESS(mmtdum%memitype, displacements(3), ierr)
    call MPI_GET_ADDRESS(mmtdum%isize,    displacements(4), ierr)
    call MPI_GET_ADDRESS(mmtdum%ncol,     displacements(5), ierr)
    call MPI_GET_ADDRESS(mmtdum%nrow,     displacements(6), ierr)
    call MPI_GET_ADDRESS(mmtdum%nlay,     displacements(7), ierr)
    base = displacements(1)
    do i = 1, 7
      displacements(i) = displacements(i) - base
    enddo
    ! -- create and commit datatype
    call MPI_TYPE_CREATE_STRUCT(7, blocklengths, displacements, types, newtype, ierr)
    call MPI_TYPE_COMMIT(newtype, ierr)
#else
    call mpiwrperror(comm, 'mpiwrpmmtstruct', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpmmtstruct
  
  subroutine mpiwrpmtstruct(mt, mmt, nmt, newtype)
  ! ******************************************************************************
  ! ******************************************************************************
  !
  !    SPECIFICATIONS:
  ! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryTypeModule, only: MemoryType,                                       &
                                ilogicalsclr, iintsclr, idblsclr,                 & 
                                iaint1d, iaint2d, iaint3d,                        & 
                                iadbl1d, iadbl2d, iadbl3d
    ! -- dummy
    integer, intent(in) :: nmt
    type(MemoryType), dimension(nmt), intent(inout) :: mt
    type(MetaMemoryType), dimension(nmt), intent(in) :: mmt
    integer, intent(out) :: newtype
    ! -- local
    integer :: i, j, isize, ncol, nrow, nlay, ierr, ntypes
    integer, dimension(:), allocatable :: types, blocklengths
    integer(KIND=MPI_ADDRESS_KIND), dimension(:), allocatable :: displacements
    integer(KIND=MPI_ADDRESS_KIND) :: base
  ! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    ntypes = 13
    allocate(types(nmt*ntypes))
    allocate(blocklengths(nmt*ntypes))
    allocate(displacements(nmt*ntypes))
    !
    j = 0
    do i = 1, nmt
      ! -- name
      j = j + 1
      types(j)        = MPI_CHARACTER
      blocklengths(j) = LENVARNAME
      call MPI_GET_ADDRESS(mt(i)%name, displacements(j), ierr)
      ! -- origin
      j = j + 1
      types(j)        = MPI_CHARACTER
      blocklengths(j) = LENORIGIN
      call MPI_GET_ADDRESS(mt(i)%origin, displacements(j), ierr)
      ! -- memitype
      j = j + 1
      types(j)        = MPI_INTEGER
      blocklengths(j) = 1
      call MPI_GET_ADDRESS(mt(i)%memitype, displacements(j), ierr)
      ! -- isize
      j = j + 1
      types(j)        = MPI_INTEGER
      blocklengths(j) = 1
      call MPI_GET_ADDRESS(mt(i)%isize, displacements(j), ierr)
      ! -- logicalsclr
      j = j + 1
      types(j) = MPI_LOGICAL
      if (mmt(i)%memitype == ilogicalsclr) then
        if (.not.associated(mt(i)%logicalsclr)) then
          allocate(mt(i)%logicalsclr)
        endif
        blocklengths(j) = 1
        call MPI_GET_ADDRESS(mt(i)%logicalsclr, displacements(j), ierr)
      else
        blocklengths(j) = 0
        displacements(j) = displacements(j-1)
      endif
      ! -- intsclr
      j = j + 1
      types(j) = MPI_INTEGER
      if (mmt(i)%memitype == iintsclr) then
        if (.not.associated(mt(i)%intsclr)) then
          allocate(mt(i)%intsclr)
          mt(i)%intsclr = 0
        endif
        blocklengths(j) = 1
        call MPI_GET_ADDRESS(mt(i)%intsclr, displacements(j), ierr)
      else
        blocklengths(j) = 0
        displacements(j) = displacements(j-1)
      endif
      ! -- idblsclr
      j = j + 1
      types(j) = MPI_DOUBLE_PRECISION
      if (mmt(i)%memitype == idblsclr) then
        if (.not.associated(mt(i)%dblsclr)) then
          allocate(mt(i)%dblsclr)
        endif
        blocklengths(j) = 1
        call MPI_GET_ADDRESS(mt(i)%dblsclr, displacements(j), ierr)
      else
        blocklengths(j) = 0
        displacements(j) = displacements(j-1)
      endif
      ! -- aint1d
      j = j + 1
      types(j) = MPI_INTEGER
      if (mmt(i)%memitype == iaint1d) then
        isize = mmt(i)%isize
        if (.not.associated(mt(i)%aint1d)) then
          allocate(mt(i)%aint1d(isize))
        endif  
        blocklengths(j) = isize
        call MPI_GET_ADDRESS(mt(i)%aint1d, displacements(j), ierr)
      else
        blocklengths(j) = 0
        displacements(j) = displacements(j-1)
      endif
      ! -- aint2d
      j = j + 1
      types(j) = MPI_INTEGER
      if (mmt(i)%memitype == iaint2d) then
        isize = mmt(i)%isize
        ncol  = mmt(i)%ncol
        nrow  = mmt(i)%nrow
        if (.not.associated(mt(i)%aint2d)) then
          allocate(mt(i)%aint2d(ncol,nrow))
        endif
        blocklengths(j) = isize
        call MPI_GET_ADDRESS(mt(i)%aint2d, displacements(j), ierr)
      else
        blocklengths(j) = 0
        displacements(j) = displacements(j-1)
      endif
      ! -- aint3d
      j = j + 1
      types(j) = MPI_INTEGER
      if (mmt(i)%memitype == iaint3d) then
        isize = mmt(i)%isize
        ncol  = mmt(i)%ncol
        nrow  = mmt(i)%nrow
        nlay  = mmt(i)%nlay
        if (.not.associated(mt(i)%aint3d)) then
          allocate(mt(i)%aint3d(ncol,nrow,nlay))
        endif
        blocklengths(j) = isize
        call MPI_GET_ADDRESS(mt(i)%aint3d, displacements(j), ierr)
      else
        blocklengths(j) = 0
        displacements(j) = displacements(j-1)
      endif
      ! -- adbl1d
      j = j + 1
      types(j) = MPI_DOUBLE_PRECISION
      if (mmt(i)%memitype == iadbl1d) then
        isize = mmt(i)%isize
        if (.not.associated(mt(i)%adbl1d)) then
          allocate(mt(i)%adbl1d(isize))
        endif
        blocklengths(j) = isize
        call MPI_GET_ADDRESS(mt(i)%adbl1d, displacements(j), ierr)
      else
        blocklengths(j) = 0
        displacements(j) = displacements(j-1)
      endif
      ! -- adbl2d
      j = j + 1
      types(j) = MPI_DOUBLE_PRECISION
      if (mmt(i)%memitype == iadbl2d) then
        isize = mmt(i)%isize
        ncol  = mmt(i)%ncol
        nrow  = mmt(i)%nrow
        if (.not.associated(mt(i)%adbl2d)) then
          allocate(mt(i)%adbl2d(ncol,nrow))
        endif
        blocklengths(j) = isize
        call MPI_GET_ADDRESS(mt(i)%adbl2d, displacements(j), ierr)
      else
        blocklengths(j) = 0
        displacements(j) = displacements(j-1)
      endif
      ! -- adbl3d
      j = j + 1
      types(j) = MPI_DOUBLE_PRECISION
      if (mmt(i)%memitype == iadbl3d) then
        isize = mmt(i)%isize
        ncol  = mmt(i)%ncol
        nrow  = mmt(i)%nrow
        nlay  = mmt(i)%nlay
        if (.not.associated(mt(i)%adbl3d)) then
          allocate(mt(i)%adbl3d(ncol,nrow,nlay))
        endif
        blocklengths(j) = isize
        call MPI_GET_ADDRESS(mt(i)%adbl3d, displacements(j), ierr)
      else
        blocklengths(j) = 0
        displacements(j) = displacements(j-1)
      endif
    enddo
    !
    base = displacements(1)
    do i = 1, nmt*ntypes
      displacements(i) = displacements(i) - base
    enddo
    !
    call MPI_TYPE_CREATE_STRUCT(nmt*ntypes, blocklengths, displacements, types, newtype, ierr)
    call MPI_TYPE_COMMIT(newtype, ierr)
    !
    deallocate(types, blocklengths, displacements)
#else
    call mpiwrperror(comm, 'mpiwrpmtstruct', 'invalid operation')
#endif    
    !
    ! -- return
    return
  end subroutine mpiwrpmtstruct
  
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