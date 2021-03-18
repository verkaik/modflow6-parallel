module MpiWrapperMemory
  
  use ConstantsModule, only: LENMEMPATH, LENVARNAME
  use KindModule, only: DP, I4B    
  
  implicit none
 
#ifdef MPI_PARALLEL
  include 'mpif.h'
#endif
  
  private
  
  type :: MetaMemoryType
    character(len=LENVARNAME) :: name
    character(len=LENMEMPATH) :: path
    integer(I4B)              :: memitype
    integer(I4B)              :: isize
    integer(I4B)              :: ncol = 0
    integer(I4B)              :: nrow = 0
    integer(I4B)              :: nlay = 0
  end type MetaMemoryType  
  
  type ColMemoryType
    character(len=LENVARNAME)  :: name        !name of the array
    character(len=LENMEMPATH)  :: path        !memory path
    integer(I4B)               :: memitype    !integer type
    logical                    :: logicalsclr !logical
    integer(I4B)               :: intsclr     !integer
    real(DP)                   :: dblsclr     !double
    integer(I4B), dimension(3) :: aint1d      !1d integer array
  end type ColMemoryType
  
  ! -- public functions
  public :: mpiwrpisend
  public :: mpiwrpirecv
  public :: mpiwrpcolstruct
  public :: mpiwrpmmtstruct
  public :: mpiwrpmtstruct
  public :: mpiwrpallgathervcol
  ! -- public variables
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
  
  contains
  
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

  subroutine mpiwrpcolstruct(newtype)
! ******************************************************************************
! Create a new MPI data type.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ConstantsModule, only: LENMEMPATH, LENVARNAME
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
    types(6) = MPI_DOUBLE_PRECISION
    types(7) = MPI_INTEGER
    blocklengths(1)  = LENVARNAME
    blocklengths(2)  = LENMEMPATH
    blocklengths(3)  = 1
    blocklengths(4)  = 1
    blocklengths(5)  = 1
    blocklengths(6)  = 1
    blocklengths(7)  = 3
    call MPI_GET_ADDRESS(cmtdum%name,       displacements(1), ierr)
    call MPI_GET_ADDRESS(cmtdum%path,       displacements(2), ierr)
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
    use ConstantsModule, only: LENMEMPATH, LENVARNAME
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
    blocklengths(2)  = LENMEMPATH
    blocklengths(3)  = 1
    blocklengths(4)  = 1
    blocklengths(5)  = 1
    blocklengths(6)  = 1
    blocklengths(7)  = 1
    call MPI_GET_ADDRESS(mmtdum%name,     displacements(1), ierr)
    call MPI_GET_ADDRESS(mmtdum%path,     displacements(2), ierr)
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
      blocklengths(j) = LENMEMPATH
      call MPI_GET_ADDRESS(mt(i)%path, displacements(j), ierr)
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
  
end module MpiWrapperMemory