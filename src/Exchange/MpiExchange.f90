  ! -- groups functionality for MPI exchanges
  module MpiExchangeModule
  use TimerModule, only: code_timer
  use KindModule, only: DP, I4B  
  use SimModule, only: ustop, store_error
  use ConstantsModule, only: LENPACKAGENAME, LENMODELNAME, LENORIGIN,          &
                             LENVARNAME, LINELENGTH
  use ArrayHandlersModule, only: ExpandArray, ifind
  use ListModule, only: ListType
  use MemoryTypeModule, only: MemoryType
  use MpiExchangeGenModule, only: serialrun, parallelrun, writestd, partstr
  use MpiWrapper, only: mpiwrpinit, mpiwrpfinalize, mpiwrpnrproc, mpiwrpmyrank,&
                        mpiwrpbarrier, mpiwrpcommworld, mpiwrpcomm_size,       &
                        mpiwrpcomm_rank, mpiwrpallgather, mpiwrpallgatherv,    &
                        mpiwrpcolstruct, mpiwrptypefree, ColMemoryType,        &
                        mpiwrpstats, mpiwrpisend, mpiwrpirecv,                 &
                        mpiwrpmmtstruct, mpiwrpmtstruct, mpiwrpwaitall,        &
                        mpiwrpprobe, mpiwrpgetcount, mpiwrpallreduce,          &
                        mpiwrpcommgroup, mpiwrpgroupincl, mpiwrpcommcreate
  use MpiWrapper, only: MetaMemoryType

  implicit none
  
  private
  
  ! -- Public types
  public :: MpiExchangeType
  public :: VarGroupType
  ! -- Public variables
  public :: MpiWorld
  public :: isrcmvr, ipckmvr, itgtmvr, iunpmvr
  ! -- Public functions
  public :: mpi_initialize_world
  public :: mpi_world_da
  !public :: mpi_add_halo_model
  public :: mpi_to_colmem
  
  integer, parameter :: MAXNVG  = 10  ! maximum number of variable groups 
  integer, parameter :: MAXNVAR = 100
  integer, parameter :: MAXNEX  = 5
  
  integer, parameter :: isrcsol = 1
  integer, parameter :: isrcgwf = 2
  integer, parameter :: isrcmvr = 3
  integer, parameter :: isrchal = 4
  integer, parameter :: isrchll = 5
  !
  integer, parameter :: ipckmm1  = 1
  integer, parameter :: ipckhm1  = 2
  integer, parameter :: ipckall  = 3
  integer, parameter :: ipckmvr  = 4
  !
  integer, parameter :: itgtsol = 1
  integer, parameter :: itgtgwf = 2
  integer, parameter :: itgtmvr = 3
  integer, parameter :: itgthal = 4
  integer, parameter :: itgthll = 5
  !
  integer, parameter :: iunpmem = 1 !memory manager
  integer, parameter :: iunphal = 2
  integer, parameter :: iunpmvr = 3
  !
  type MpiGwfBuf
    integer(I4B)                                :: nsnd = 0         ! total number of variables to send
    integer(I4B)                                :: nrcv = 0         ! total number of variables to receive
    type(MemoryType), dimension(:), pointer     :: sndmt => null()  ! send data
    type(MemoryType), dimension(:), pointer     :: rcvmt => null()  ! receive data
    type(MetaMemoryType), dimension(:), pointer :: sndmmt => null() ! send meta data
    type(MetaMemoryType), dimension(:), pointer :: rcvmmt => null() ! receive meta data
  end type MpiGwfBuf

  type VarType
    logical                   :: lsnd = .true.
    logical                   :: lrcv = .true.
    character(len=LENORIGIN)  :: origin
    character(len=LENVARNAME) :: name
    character(len=LENVARNAME) :: nameext
    integer(I4B)              :: srctype
    integer(I4B)              :: pcktype
    integer(I4B)              :: tgttype
    integer(I4B)              :: unptype
    !logical                   :: lhalonode = .true.
    integer(I4B), dimension(1) :: id = 0
  end type Vartype

  type VarGroupType
    integer(I4B)                      :: nvar = 0
    type(VarType), dimension(MAXNVAR) :: var
  end type VarGroupType
  
  type ExchangeType
    type(VarGroupType), dimension(:), pointer :: vgvar => null()
    character(len=LENVARNAME)                 :: name
    character(len=LENMODELNAME)               :: halo_name
    character(len=LENMODELNAME)               :: m1_name
    character(len=LENMODELNAME)               :: m2_name
    character(len=4)                          :: m1_dis
    character(len=4)                          :: m2_dis
    integer(I4B)                              :: halo_offset
    integer(I4B)                              :: m1_id !CGC
    integer(I4B)                              :: m2_id !CGC
  end type ExchangeType
  
  ! -- Local types
  type :: MpiGwfCommInt
    integer(I4B), pointer                     :: xprnk     => null() ! neighboring rank
    integer(I4B), pointer                     :: nexchange => null() ! number of exchanges 
    type(ExchangeType), dimension(:), pointer :: exchange  => null() ! exchanges
    type(MpiGwfBuf), dimension(:), pointer    :: vgbuf     => null() ! variable group buffers over all exchanges
  end type MpiGwfCommInt
  
  type :: MpiExchangeType
    logical                                                :: linit = .false.
    character(len=LENPACKAGENAME)                          :: name          ! name (origin)
    character(len=LENPACKAGENAME)                          :: solname       ! solution name (origin)
    integer(I4B)                                           :: gnmodel = 0   ! number of global models
    character(len=LENMODELNAME), dimension(:), allocatable :: gmodelnames   ! global model names
    integer(I4B)                                           :: gnsub = 0     ! number of global subdomains
    integer(I4B), dimension(:), allocatable                :: gsubs         ! global subdomains
    integer(I4B)                                           :: lnmodel = 0   ! number of local models
    character(len=LENMODELNAME), dimension(:), allocatable :: lmodelnames   ! local model names
    integer(I4B)                                           :: lnsub = 0     ! number of local subdomains
    integer(I4B), dimension(:), allocatable                :: lsubs         ! model local subdomains
    integer(I4B)                                           :: hnmodel = 0   ! number of halo models
    character(len=LENMODELNAME), dimension(:), allocatable :: hmodelnames   ! halo model names
    character(len=LENMODELNAME), dimension(:), allocatable :: hmodelm1names   ! halo model m1 names
    character(len=LENMODELNAME), dimension(:), allocatable :: hmodelm2names   ! halo model m2 names
    integer(I4B), pointer                                :: comm   => null() ! MPI communicator
    integer(I4B), pointer                                :: nrproc => null() ! number of MPI process for this communicator
    integer(I4B), pointer                                :: myrank => null() ! MPI rank in for this communicator
    integer(I4B), pointer                                :: myproc => null() ! MPI proc in for this communicator
    integer(I4B), dimension(:), pointer                  :: procmap => null() ! Mapping
    integer(I4B), pointer                                :: nrxp => null() ! Number of exchange partners
    type(MpiGwfCommInt), dimension(:), pointer           :: lxch => null() ! Point to-point communication structure
    integer(I4B)                                         :: nvg  ! number of variable groups
    character(len=LINELENGTH), dimension(MAXNVG)         :: vg   ! variable groups
    logical,  dimension(MAXNVG)                          :: lxchmeta = .true. ! exchange meta data 
    character(len=50), pointer                           :: nrprocstr => null() ! Number of processes string
    integer(I4B), pointer                                :: npdigits  => null() ! Number of digits for nrproc
    character(len=50), pointer                           :: partstr   => null() ! Partition string
    logical                                              :: lp2p = .true. ! flag indicating if point-to-point communication is necessary
    real(DP)                                             :: ttgsum = 0.d0
    real(DP)                                             :: ttgmax = 0.d0
    real(DP)                                             :: ttgmin = 0.d0
    real(DP)                                             :: ttbarr = 0.d0
    real(DP)                                             :: ttpack = 0.d0
    real(DP)                                             :: ttupck = 0.d0
  contains
    procedure :: mpi_barrier
    procedure :: mpi_create_output_str
    procedure :: mpi_is_iproc
    procedure :: mpi_addmodel
    procedure :: mpi_addhmodel
    procedure :: mpi_getmodel
    procedure :: mpi_addsub
    procedure :: mpi_local_exchange_init
    procedure :: mpi_local_exchange
    procedure :: mpi_add_vg
    procedure :: mpi_add_vmt
    procedure :: mpi_init_vg
    procedure :: mpi_mv_halo
    procedure :: mpi_global_exchange_sum1
    procedure :: mpi_global_exchange_sum2
    generic   :: mpi_global_exchange_sum => mpi_global_exchange_sum1,          &
                                            mpi_global_exchange_sum2
    procedure :: mpi_global_exchange_absmin1
    procedure :: mpi_global_exchange_absmin2
    generic   :: mpi_global_exchange_absmin => mpi_global_exchange_absmin1,    &
                                               mpi_global_exchange_absmin2
    procedure :: mpi_global_exchange_absmax1
    procedure :: mpi_global_exchange_absmax2
    generic   :: mpi_global_exchange_absmax => mpi_global_exchange_absmax1,    &
                                               mpi_global_exchange_absmax2
    procedure :: mpi_global_exchange_max => mpi_global_exchange_max_int
    procedure :: mpi_global_exchange_cgc !CGC
    procedure :: mpi_debug
    procedure :: mpi_da
    procedure :: mpi_not_supported
    procedure :: mpi_get_halo_rcvmt
    procedure :: mpi_copy_int_to_halo
    generic   :: mpicopyinttohalo => mpi_copy_int_to_halo
    procedure :: mpi_copy_dbl_to_halo
    generic   :: mpicopydbltohalo => mpi_copy_dbl_to_halo
  end type MpiExchangeType
  
  ! -- World communicator
  type(MpiExchangeType), pointer :: MpiWorld => null()

  character(len=LINELENGTH) :: errmsg
  
  real(DP) :: t
  
  save
  
  contains
  
  subroutine mpi_initialize_world()
! ******************************************************************************
! MPI initialization for the world communicator.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryManagerModule, only: mem_allocate !, mem_get_info
    use MpiExchangeGenModule,   only: mpi_initialize
    ! -- dummy
    ! -- local
    character(len=LENORIGIN) :: origin
! ------------------------------------------------------------------------------
    !
    ! -- Initialize MPI
    call mpi_initialize()
    !
    ! -- Allocate MpiWorld object
    allocate(MpiWorld)
    !
    ! -- Set name
    MpiWorld%name = 'MPI_WORLD'
    !
    ! -- Allocate scalars
    origin = MpiWorld%name
    call mem_allocate(MpiWorld%comm,   'COMM',   origin)
    call mem_allocate(MpiWorld%nrproc, 'NRPROC', origin)
    call mem_allocate(MpiWorld%myrank, 'MYRANK', origin)
    call mem_allocate(MpiWorld%myproc, 'MYPROC', origin)
    call mem_allocate(MpiWorld%nrxp,   'NRXP',   origin)
    !
    MpiWorld%comm   = mpiwrpcommworld()
    MpiWorld%nrproc = mpiwrpcomm_size(MpiWorld%comm)
    MpiWorld%myrank = mpiwrpcomm_rank(MpiWorld%comm)
    MpiWorld%myproc = MpiWorld%myrank + 1
    MpiWorld%nrxp   = 0
    !
    if (MpiWorld%nrproc == 1) then
      serialrun = .true.
    else
      serialrun = .false.
    end if
    parallelrun = .not.serialrun
    ! @@@@@@ DEBUG
    !serialrun = .false.
    ! @@@@@@ DEBUG
    !
    if (MpiWorld%myrank == 0) then
      writestd = .true.
    else
      writestd = .false.
    end if
    !
    ! -- output strings
    call MpiWorld%mpi_create_output_str()
    partstr = MpiWorld%partstr
    !
    ! -- return
    return
  end subroutine mpi_initialize_world
  
  subroutine mpi_world_da()
! ******************************************************************************
! Deallocate for the world communicator.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    ! -- local
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    end if
    !
    call MpiWorld%mpi_barrier()
    call MpiWorld%mpi_da()
    deallocate(MpiWorld)
    !
    ! -- return
    return
  end subroutine mpi_world_da
  
  subroutine mpi_barrier(this)
! ******************************************************************************
! MPI barrier.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiExchangeType) :: this
    ! -- local
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    end if
    !
#ifdef MPI_TIMER
    call code_timer(0, t, this%ttbarr)
#endif
    call mpiwrpbarrier(this%comm)
#ifdef MPI_TIMER
    call code_timer(1, t, this%ttbarr)
#endif
    !
    ! -- return
    return
  end subroutine mpi_barrier
  
  subroutine mpi_create_output_str(this)
! ******************************************************************************
! Create several strings for output.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(MpiExchangeType) :: this
    ! -- local
    character(len=LENORIGIN) :: origin
    character(len=20) :: fmt  
! ------------------------------------------------------------------------------
    ! -- Allocate
    origin = this%name
    allocate(this%nrprocstr)
    call mem_allocate(this%npdigits, 'NPDIGITS', origin)
    allocate(this%partstr)
    !
    write(this%nrprocstr,*) this%nrproc
    this%nrprocstr = adjustl(this%nrprocstr)
    this%npdigits = len_trim(this%nrprocstr)
    write(fmt,*) this%npdigits
    fmt = adjustl(fmt)
    fmt = '(a,i'//trim(fmt)//'.'//trim(fmt)//')'
    write(this%partstr,fmt) 'p',this%myrank
    !    
    ! -- return
    return
  end subroutine mpi_create_output_str
  
  subroutine mpi_debug(this)
! ******************************************************************************
! Debug subroutine.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiExchangeType) :: this
    ! -- local
! ------------------------------------------------------------------------------
    if (this%myrank == 0) then
      write(*,*) 'Press a key...' 
      pause
    end if
    call mpi_barrier(this)
    ! -- return
    return
  end subroutine mpi_debug
    
  subroutine mpi_local_exchange_init(this)
! ******************************************************************************
! This subroutine initializes the local exchange data structure for the
! solution communicator.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(MpiExchangeType) :: this
    ! -- local
    character(len=LENORIGIN) :: origin
    integer(I4B) :: im, i, j, isub
    integer(I4B), dimension(:), allocatable :: ranks
    integer(I4B) :: world_group, new_group
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    end if
    !
    ! -- Allocate scalars
    origin = this%name
    call mem_allocate(this%comm,   'COMM',   origin)
    call mem_allocate(this%nrproc, 'NRPROC', origin)
    call mem_allocate(this%myrank, 'MYRANK', origin)
    call mem_allocate(this%myproc, 'MYPROC', origin)
    call mem_allocate(this%nrxp,   'NRXP',   origin)
    !
    ! -- loop over the models associated with this solution
    allocate(this%procmap(MpiWorld%nrproc))
    this%procmap = 0
    do im = 1, this%gnmodel
      ! -- find the global subdomain number for this model
      i = ifind(MpiWorld%gmodelnames, this%gmodelnames(im))
      ! -- add subdomain to local MPI subdomain list
      if (i > 0) then
        isub = MpiWorld%gsubs(i)
        this%procmap(isub) = 1
      else
        call store_error('Program error: mpi_local_exchange.')
        call ustop()
      end if
    end do
    this%nrproc = sum(this%procmap)
    allocate(ranks(max(this%nrproc,1)))
    j = 0
    do i = 1, MpiWorld%nrproc
      if (this%procmap(i) == 1) then
        j = j + 1
        ranks(j) = i-1
        this%procmap(i) = j
      end if
    end do
    !
    call mpiwrpcommgroup(MpiWorld%comm, world_group)
    call mpiwrpgroupincl(world_group, this%nrproc, ranks, new_group)
    call mpiwrpcommcreate(MpiWorld%comm, new_group, this%comm)
    !
    this%myrank = mpiwrpcomm_rank(this%comm)
    this%myproc = this%myrank + 1
    !
    ! clean up
    deallocate(ranks)
    !
    !this%comm   = MpiWorld%comm
    !this%nrproc = MpiWorld%nrproc
    !this%myrank = MpiWorld%myrank
    !this%myproc = MpiWorld%myproc
    !
    ! -- return
    return
  end subroutine mpi_local_exchange_init

  subroutine mpi_add_vg(this, vgname)
! ******************************************************************************
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiExchangeType) :: this
    character(len=*), intent(in) :: vgname
    ! -- local
    integer(I4B) :: ivg, n, ixp, iex
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    end if
    !
    if (.not.this%linit) then
      call store_error('Program error: mpi_add_vmt')
      call ustop()
    end if
    !
    ivg = ifind(this%vg, trim(vgname))
    if (ivg < 0) then
      this%nvg = this%nvg + 1
      n = this%nvg
      this%vg(n) = trim(vgname)
    else
      write(errmsg,'(a)') 'Program error mpi_add_vg'
      call store_error(errmsg)
      call ustop()
    end if
    !
    ! -- Allocate arrays
    do ixp = 1, this%nrxp
      if (.not.associated(this%lxch(ixp)%vgbuf)) then
        allocate(this%lxch(ixp)%vgbuf(MAXNVG))
      end if
      do iex = 1, this%lxch(ixp)%nexchange
        if (.not.associated(this%lxch(ixp)%exchange(iex)%vgvar)) then
          allocate(this%lxch(ixp)%exchange(iex)%vgvar(MAXNVG))
        end if
      end do
    end do
    !
    ! -- return
    return
  end subroutine mpi_add_vg

  subroutine mpi_add_vmt(this, vgname, name, nameext, srctype, pcktype,         &
     tgttype, unptype)
! ******************************************************************************
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiExchangeType) :: this
    character(len=*), intent(in) :: vgname
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: nameext
    character(len=*), intent(in) :: srctype
    character(len=*), intent(in) :: pcktype
    character(len=*), intent(in) :: tgttype
    character(len=*), intent(in) :: unptype
    !logical, intent(in), optional :: lhalonode
    ! -- local
    type(VarGroupType), pointer :: vgvar
    integer(I4B) :: ivg, ixp, iex, iv
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    end if
    !
    if (.not.this%linit) then
      call store_error('Program error 1: mpi_add_vmt')
      call ustop()
    end if
    !
    ivg = ifind(this%vg, trim(vgname))
    if (ivg < 0) then
      write(errmsg,'(a)') 'Program error 2: mpi_add_vmt'
      call store_error(errmsg)
      call ustop()
    end if
    !
    do ixp = 1, this%nrxp
      do iex = 1, this%lxch(ixp)%nexchange
        vgvar => this%lxch(ixp)%exchange(iex)%vgvar(ivg)
        iv = vgvar%nvar
        iv = iv + 1 
        if (iv > MAXNVAR) then
          call store_error('Program error 3: mpi_add_vmt')
          call ustop()
        end if
        vgvar%nvar = iv
        vgvar%var(iv)%name    = trim(name)
        vgvar%var(iv)%nameext = trim(nameext)
        !if (present(lhalonode)) then
        !  vgvar%var(iv)%lhalonode = lhalonode
        !end if
        select case(srctype)
          case('SOL')
            vgvar%var(iv)%srctype = isrcsol
          case('GWF')
            vgvar%var(iv)%srctype = isrcgwf
          case('MVR')        
            vgvar%var(iv)%srctype = isrcmvr
          case('HAL')        
            vgvar%var(iv)%srctype = isrchal
          case('HLL')        
            vgvar%var(iv)%srctype = isrchll
          case default
            call store_error('Program error 4: mpi_add_vmt')
            call ustop()
        end select
        select case(pcktype)
          case('MM1')
            vgvar%var(iv)%pcktype = ipckmm1
          case('HM1')
            vgvar%var(iv)%pcktype = ipckhm1
          case('ALL')
            vgvar%var(iv)%pcktype = ipckall
          case default
            call store_error('Program error 5: mpi_add_vmt')
            call ustop()
        end select
        select case(tgttype)
          case('SOL')
            vgvar%var(iv)%tgttype = itgtsol
          case('GWF')
            vgvar%var(iv)%tgttype = itgtgwf
          case('MVR')        
            vgvar%var(iv)%tgttype = itgtmvr
          case('HAL')        
            vgvar%var(iv)%tgttype = itgthal
          case('HLL')        
            vgvar%var(iv)%tgttype = itgthll
          case default
            call store_error('Program error 6: mpi_add_vmt')
            call ustop()
        end select
        select case(unptype)
        case('MEM')
            vgvar%var(iv)%unptype = iunpmem
          case('HAL')
            vgvar%var(iv)%unptype = iunphal
          case('MVR')
            vgvar%var(iv)%unptype = iunpmvr
         case default 
            call store_error('Program error 7: mpi_add_vmt')
            call ustop()
        end select
      end do
    end do
    !
    ! -- return
    return
  end subroutine mpi_add_vmt
  
 subroutine mpi_init_vg(this, vgname)
! ******************************************************************************
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiExchangeType) :: this
    character(len=*), intent(in) :: vgname
    ! -- local
    type(VarGroupType), pointer :: vgvar
    type(MpiGwfBuf), pointer :: buf
    integer(I4B) :: ivg, ixp, iex, iv, nsnd
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    end if
    !
    ivg = ifind(this%vg, trim(vgname))
    if (ivg < 0) then
      write(errmsg,'(a)') 'Program error mpi_init_vg'
      call store_error(errmsg)
      call ustop()
    end if
    !
    do ixp = 1, this%nrxp
      nsnd = 0
      do iex = 1, this%lxch(ixp)%nexchange
        vgvar => this%lxch(ixp)%exchange(iex)%vgvar(ivg)
        do iv = 1, vgvar%nvar
          if (vgvar%var(iv)%lsnd) then
            nsnd = nsnd + 1
          end if
        end do  
      end do
      buf => this%lxch(ixp)%vgbuf(ivg)
      buf%nsnd = nsnd 
      allocate(buf%sndmt(max(nsnd,1)))
      allocate(buf%sndmmt(max(nsnd,1)))
    end do
    !
    ! -- return
    return
  end subroutine mpi_init_vg  
  
  subroutine mpi_local_exchange(this, solname, vgname, lunpack)
! ******************************************************************************
! Local point-to-point exchange (wrapper).
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryTypeModule, only: iintsclr, idblsclr, iaint1d, iadbl1d !@@@@DEBUG
    use MemoryManagerModule, only: mem_get_ptr, mem_setptr, mem_setval,        &
                                   mem_setval_id
    use MpiWrapper, only: mpiwrpstats
    ! -- dummy
    class(MpiExchangeType) :: this
    character(len=*), intent(in) :: solname
    character(len=*), intent(in) :: vgname
    logical, intent(in) :: lunpack 
    ! -- local
    !type(VarMemoryType), pointer :: vmt
    !
    type(MemoryType), pointer :: mt
    type(ExchangeType), pointer :: ex
    type(VarGroupType), pointer :: vgvar
    type(VarType), pointer :: var
    type(MetaMemoryType), pointer :: rcvmmt
    type(MpiGwfBuf), pointer :: vgbuf

    character(len=LENORIGIN) :: mod_origin, sol_origin, src_origin, tgt_origin
    integer(I4B) :: ixp, iex, nrcv, iv, is, i, j, istat
    character(len=LENMODELNAME) :: mname, id, m1_name, m2_name, halo_name
    !
    integer(I4B) :: newtype
    integer(I4B), dimension(:), allocatable :: snd_newtype, rcv_newtype
    integer(I4B), dimension(:), allocatable :: sreq, rreq
    !
    integer(I4B) :: irank, memitype !@@@DEBUG
    !
    integer(I4B) :: ivg
    integer(I4B) :: nexg
    integer(I4B), dimension(:), pointer :: nodem1
    integer(I4B) ::  moffset
    integer(I4B), pointer :: tmp
    
    logical :: lpack, lpackmeta
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    end if
    !
    if (.not.this%lp2p) then
      return
    end if
    !
    lpack = .true.
    lpackmeta = .true.
    !
    ! -- find the variable group name
    ivg = ifind(this%vg, trim(vgname))
    if (ivg < 0) then
      call store_error('Program error: mpi_local_exchange for '//trim(vgname))
      call ustop()
    end if
    !
    ! ---
#ifdef MPI_TIMER
    call code_timer(0, t, this%ttpack)
#endif
    if (lpack) then
      ! -- Loop over the exchange partners
      do ixp = 1, this%nrxp
        vgbuf => this%lxch(ixp)%vgbuf(ivg)
        is = 0
        ! -- Loop over the exchanges
        do iex = 1, this%lxch(ixp)%nexchange
          ex => this%lxch(ixp)%exchange(iex)
          vgvar => ex%vgvar(ivg)
          m1_name = ex%m1_name
          halo_name = ex%halo_name
          ! -- Loop over the exchange send variables
          do iv = 1, vgvar%nvar
            var => vgvar%var(iv)
            if (var%lsnd) then
              is = is + 1
              ! -- check
              if (is > vgbuf%nsnd) then
                !write(*,*) '@@@ packing '//trim(vgname)//' '//trim(var%name),is
                call store_error('Program error 1: mpi_local_exchange')
                call ustop()
              end if
              !
              ! source and offset
              select case(var%srctype)
                case(isrcgwf)
                  if (trim(var%nameext) == 'DIS') then
                    write(mod_origin,'(a,1x,a)') trim(m1_name), trim(ex%m1_dis)
                  else
                    write(mod_origin,'(a,1x,a)') trim(m1_name), trim(var%nameext)
                  endif
                  src_origin = mod_origin
                  moffset = 0
                case(isrcsol)
                  read(solname,*) sol_origin
                  write(sol_origin,'(a,1x,a)') trim(sol_origin), trim(var%nameext)
                  src_origin = sol_origin
                  call mem_setptr(tmp, 'MOFFSET', trim(m1_name))
                  moffset = tmp
                case(isrchal)
                  src_origin = trim(halo_name)
                  moffset = 0
                case(isrchll)
                  src_origin = trim(halo_name)//'_M1 '//trim(var%nameext)
                  moffset = 0
                case(isrcmvr)
                  src_origin = var%origin
                  moffset = 0
              end select
              !
              !write(*,*) '@@@@ Getting "'//trim(var%name)//'" "'//trim(src_origin)//'" for rank',this%myrank
              call mem_get_ptr(var%name, src_origin, mt)
              !
              select case(var%pcktype)
                case(ipckmm1)
                  call mem_setptr(tmp, 'NEXG', trim(ex%name))
                  call mem_setptr(nodem1, 'NODEM1', trim(ex%name))
                  if (.not.associated(tmp).or..not.associated(nodem1)) then
                    call store_error('Program error 2: mpi_local_exchange')
                    call ustop()
                  end if
                  nexg = tmp
                  call mpi_pack_mt(src_origin,  mt, vgbuf%sndmt(is), nodem1, nexg, moffset)
                case(ipckhm1)
                  call mem_setptr(nodem1, 'IMAPMTOHALO', trim(halo_name)//'_M1')
                  nexg = size(nodem1)
                  call mpi_pack_mt(src_origin, mt, vgbuf%sndmt(is), nodem1, nexg, moffset)
                case(ipckmvr)
                  nexg = 1
                  nodem1 => var%id
                  call mpi_pack_mt(src_origin, mt, vgbuf%sndmt(is), nodem1, nexg, moffset)
                case(ipckall)
                  call mpi_pack_mt(tgt_origin, mt, vgbuf%sndmt(is))
                  !write(*,*) '@@@@ Packed "'//trim(var%name)//'" "'//trim(src_origin)//'" for rank',this%myrank
              end select
              !
            end if
          end do
        end do
      end do
    end if
    !
#ifdef MPI_TIMER
    call code_timer(1, t, this%ttpack)
#endif
    !-- Pack the meta data
    if (lpackmeta) then
      ! -- Loop over the exchange partners
      do ixp = 1, this%nrxp
        vgbuf => this%lxch(ixp)%vgbuf(ivg)
        ! -- Loop over all send variables
        do is = 1, vgbuf%nsnd
          call mpi_pack_mmt(vgbuf%sndmt(is), vgbuf%sndmmt(is))
        end do
      end do
    end if
    !
    ! -- First point-to-point communication to exchange memory type meta-data
    !
    ! -- Allocate request arrays
    allocate(sreq(this%nrxp), rreq(this%nrxp))
    if (this%lxchmeta(ivg)) then
      ! -- Create MPI derived datatype
      call mpiwrpmmtstruct(newtype)
      ! -- Non-blocking send
      do ixp = 1, this%nrxp
        call mpiwrpisend(this%lxch(ixp)%vgbuf(ivg)%sndmmt,                      &
                         this%lxch(ixp)%vgbuf(ivg)%nsnd,                        &
                         newtype, this%lxch(ixp)%xprnk, 0, this%comm, sreq(ixp))
      end do
      ! -- Probe for the data sizes and allocate
      do ixp = 1, this%nrxp
        vgbuf => this%lxch(ixp)%vgbuf(ivg)
        call mpiwrpprobe(this%lxch(ixp)%xprnk, 0, this%comm, mpiwrpstats(1,ixp))
        call mpiwrpgetcount(mpiwrpstats(1,ixp), newtype, nrcv)
        vgbuf%nrcv = nrcv
        ! -- Allocate receive buffers
        if (.not.associated(vgbuf%rcvmmt)) then
          allocate(vgbuf%rcvmmt(max(1,nrcv)))
        end if
        if (.not.associated(vgbuf%rcvmt)) then
          allocate(vgbuf%rcvmt(max(1,nrcv)))
        end if
      end do
      ! -- Non-blocking receive
      do ixp = 1, this%nrxp
        vgbuf => this%lxch(ixp)%vgbuf(ivg)
        call mpiwrpirecv(vgbuf%rcvmmt, vgbuf%nrcv,                              &
                         newtype, this%lxch(ixp)%xprnk, 0, this%comm, rreq(ixp))
      end do
      call mpiwrpwaitall(this%nrxp, sreq, mpiwrpstats)
      call mpiwrpwaitall(this%nrxp, rreq, mpiwrpstats)
      ! -- Clean up MPI datatype
      call mpiwrptypefree(newtype)
      ! -- Debug
      if (.false.) then
        do irank = 0, this%nrproc-1
          if (irank == this%myrank) then
            write(*,*) '====myrank',this%myrank
            do ixp = 1, this%nrxp
              write(*,*) '=======received from',this%lxch(ixp)%xprnk
              do i = 1, this%lxch(ixp)%vgbuf(ivg)%nrcv
                rcvmmt => this%lxch(ixp)%vgbuf(ivg)%rcvmmt(i)
                write(*,*) 'name:     ', trim(rcvmmt%name)
                write(*,*) 'origin:   ', trim(rcvmmt%origin)
                write(*,*) 'memitype: ', rcvmmt%memitype
                write(*,*) 'isize:    ', rcvmmt%isize
                write(*,*) 'ncol:     ', rcvmmt%ncol
                write(*,*) 'nrow:     ', rcvmmt%nrow
              end do
            end do
          end if
          call mpiwrpbarrier(this%comm)
        end do
      end if
      !
      this%lxchmeta(ivg) = .false.
    end if
    !
    ! -- Second point-to-point communication to exchange actual data
    
    ! -- Create MPI send and receive datatypes
    allocate(snd_newtype(this%nrxp), rcv_newtype(this%nrxp))
    do ixp = 1, this%nrxp
      vgbuf => this%lxch(ixp)%vgbuf(ivg)
      call mpiwrpmtstruct(vgbuf%sndmt, vgbuf%sndmmt,                            &
                          vgbuf%nsnd, snd_newtype(ixp))
      call mpiwrpmtstruct(vgbuf%rcvmt, vgbuf%rcvmmt,                            &
                          vgbuf%nrcv, rcv_newtype(ixp))
    end do
    !
    ! -- Non-blocking send
    do ixp = 1, this%nrxp
      vgbuf => this%lxch(ixp)%vgbuf(ivg)
      call mpiwrpisend(vgbuf%sndmt, 1, snd_newtype(ixp),                        &
                       this%lxch(ixp)%xprnk, 0, this%comm, sreq(ixp))
    end do
    ! -- Non-blocking receive
    do ixp = 1, this%nrxp
      vgbuf => this%lxch(ixp)%vgbuf(ivg)
      call mpiwrpirecv(vgbuf%rcvmt, 1, rcv_newtype(ixp),                        &
                       this%lxch(ixp)%xprnk, 0, this%comm, rreq(ixp))
    end do
    call mpiwrpwaitall(this%nrxp, sreq, mpiwrpstats)
    call mpiwrpwaitall(this%nrxp, rreq, mpiwrpstats)
    ! -- Clean up MPI datatypes
    do ixp = 1, this%nrxp
      call mpiwrptypefree(snd_newtype(ixp))
      call mpiwrptypefree(rcv_newtype(ixp))
    end do
    deallocate(snd_newtype, rcv_newtype)
    deallocate(sreq, rreq)
    !
    ! -- Debug
    if (.false. .and. trim(vgname) == 'MOVER') then
      write(*,*) '@@@@@ VARGROUP = '//trim(vgname)
      do irank = 0, this%nrproc-1
        if (irank == this%myrank) then
          write(*,*) '=================myrank',this%myrank
          do ixp = 1, this%nrxp
            vgbuf => this%lxch(ixp)%vgbuf(ivg)
            write(*,*) '---sent to',this%lxch(ixp)%xprnk, vgbuf%nsnd
            do i = 1, vgbuf%nsnd
              !if (trim(vgbuf%sndmt(i)%name) /= 'X') cycle
              write(*,*) 'name:      ', trim(vgbuf%sndmt(i)%name)
              write(*,*) 'origin:    ', trim(vgbuf%sndmt(i)%origin)
              write(*,*) 'isize:     ', vgbuf%sndmt(i)%isize
              memitype = vgbuf%sndmt(i)%memitype
              !write(*,*) 'memitype:  ', memitype
              select case(memitype)
                case(iadbl1d)
                  do j = 1, vgbuf%sndmt(i)%isize
                    write(*,*) 'dval:      ', vgbuf%sndmt(i)%adbl1d(j)
                  end do
                case(iaint1d)
                  do j = 1, vgbuf%sndmt(i)%isize
                    write(*,*) 'ival:      ', vgbuf%sndmt(i)%aint1d(j)
                  end do
              end select
            end do
            write(*,*) '---received from',this%lxch(ixp)%xprnk, vgbuf%nrcv
            do i = 1,vgbuf%nrcv
              !if (trim(vgbuf%rcvmt(i)%name) /= 'X') cycle
              write(*,*) 'name:      ', trim(vgbuf%rcvmt(i)%name)
              write(*,*) 'origin:    ', trim(vgbuf%rcvmt(i)%origin)
              write(*,*) 'isize:     ', vgbuf%rcvmt(i)%isize
              !write(*,*) 'memitype:  ', vgbuf%rcvmt(i)%memitype
              memitype = vgbuf%rcvmt(i)%memitype
              !write(*,*) 'memitype:  ', memitype
              select case(memitype)
                case(iadbl1d)
                  do j = 1, vgbuf%rcvmt(i)%isize
                    write(*,*) 'dval:      ', vgbuf%rcvmt(i)%adbl1d(j)
                  end do
                case(iaint1d)
                  do j = 1, vgbuf%rcvmt(i)%isize
                    write(*,*) 'ival:      ', vgbuf%rcvmt(i)%aint1d(j)
                  end do
               end select
            end do
          end do
        end if
        call mpiwrpbarrier(this%comm)
      end do
      !call ustop('@@@@@debug')
    end if
    !
    ! -- Set the received data for the halo (m2) models
#ifdef MPI_TIMER
    call code_timer(0, t, this%ttupck)
#endif
    if (lunpack) then
      do ixp = 1, this%nrxp
        vgbuf => this%lxch(ixp)%vgbuf(ivg)
        is = 0
        do iex = 1, this%lxch(ixp)%nexchange
          ex => this%lxch(ixp)%exchange(iex)
          ! -- halo model name
          m2_name = ex%m2_name
          vgvar => ex%vgvar(ivg)
          halo_name = ex%halo_name
          do iv = 1, vgvar%nvar
            var => vgvar%var(iv)
            if (var%lrcv) then
              is = is + 1
              ! check
              if (vgbuf%rcvmt(is)%name /= var%name) then
                call store_error('Program error 5: mpi_local_exchange')
                call ustop()
              end if
              !
              ! target origin
              select case (var%tgttype)
                case(itgtgwf,itgtsol)
                  read(vgbuf%rcvmt(is)%origin,*,iostat=istat) mname, id
                  if (istat /= 0) then
                    read(vgbuf%rcvmt(is)%origin,*,iostat=istat) mname
                    id = ''
                  end if
                  if (index(id,'HALO') > 0) then
                    id = ''
                  end if
                  if (trim(id) == 'DIS') then
                    id = ex%m2_dis
                  end if
                  tgt_origin = trim(m2_name)//' '//trim(id)
                case(itgthal)
                  read(vgbuf%rcvmt(is)%origin,*,iostat=istat) mname, id
                  if (istat /= 0) then
                    read(vgbuf%rcvmt(is)%origin,*,iostat=istat) mname
                    id = ''
                  end if
                  if (trim(id) == 'DIS') then
                    id = 'DISU'
                  end if
                  tgt_origin = trim(halo_name)//' '//trim(id)
                case(itgthll)  
                  tgt_origin = trim(halo_name)//'_M2 '//trim(var%nameext)
                case(itgtmvr)
                  tgt_origin = var%origin
              end select
              !
              vgbuf%rcvmt(is)%origin = tgt_origin
              !
              select case(var%unptype)
                case(iunpmem)
                  !if (trim(vgname)=='HALO_INIT_CON_A' .and. trim(var%name)=='CL1') then
                  !  write(*,*) '-->Unpacking for HALO_A: '//trim(var%name)//' '//trim(vgbuf%rcvmt(is)%origin), var%unptype
                  !end if
                  call mem_setval(vgbuf%rcvmt(is))
                case(iunpmvr)
                  call mem_setval_id(vgbuf%rcvmt(is), var%id, 1)
                case(iunphal)
                  call mem_setptr(nodem1, 'IMAPMTOHALO', trim(halo_name)//'_M1')
                  if (.not.associated(nodem1)) then
                    call store_error('Program error 6: mpi_local_exchange')
                    call ustop()
                  endif
                  call mem_setval(vgbuf%rcvmt(is), size(nodem1))
              end select
            end if
          end do
        end do
      end do
    end if
#ifdef MPI_TIMER
    call code_timer(1, t, this%ttupck)
#endif
    !
    ! -- return
    return
  end subroutine mpi_local_exchange
  
  subroutine mpi_mv_halo(this, solname, vgname, x, m2_id) !CGC
! ******************************************************************************
! Add terms for halo matrix-vector product.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ConstantsModule, only: LENORIGIN, DZERO
    use MemoryManagerModule, only: mem_setptr
    use MemoryTypeModule, only: iadbl1d
    ! -- dummy
    class(MpiExchangeType) :: this
    character(len=*), intent(in) :: solname
    character(len=*), intent(in) :: vgname
    real(DP), dimension(*), intent(inout) :: x
    integer(I4B), intent(in), optional :: m2_id !CGC
    ! -- local
    type(MpiGwfBuf), pointer :: vgbuf
    character(len=LENORIGIN) :: sol_origin
    integer(I4B) :: ivg, ixp, ix, iv, iexg, n
    integer(I4B), pointer :: nexg
    integer(I4B), dimension(:), pointer :: nodem1
    integer(I4B), dimension(:), pointer :: active
    integer(I4B), pointer :: moffset
    real(DP), dimension(:), pointer :: cond
    real(DP), dimension(:), pointer :: newtonterm
    real(DP) :: v
    logical :: lm2_id !CGC
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    end if
    !
    if (present(m2_id)) then !CGC
      lm2_id = .true. !CGC
    else !CGC
      lm2_id = .false. !CG
    end if !CGC
    !
    ! -- find the variable group name
    ivg = ifind(this%vg, trim(vgname))
    if (ivg < 0) then
      call store_error('Program error: mpi_mv_halo for '//trim(vgname))
      call ustop()
    end if
    !
    read(solname,*) sol_origin
    call mem_setptr(active, 'IACTIVE', trim(sol_origin))
    !
    do ixp = 1, this%nrxp
      vgbuf => this%lxch(ixp)%vgbuf(ivg)
      do ix = 1, this%lxch(ixp)%nexchange
        ! -- cycle in case global model 2 ID does not match the optional argument
        if (lm2_id) then !CGC
          if (this%lxch(ixp)%exchange(ix)%m2_id /= m2_id) then !CGC
            cycle !CGC
          end if !CGC
        end if !CGC
        ! 
        ! -- get exchange nodes and conductance
        call mem_setptr(nexg, 'NEXG', trim(this%lxch(ixp)%exchange(ix)%name))
        call mem_setptr(nodem1, 'NODEM1',                                      &
          trim(this%lxch(ixp)%exchange(ix)%name))
        call mem_setptr(cond, 'COND', trim(this%lxch(ixp)%exchange(ix)%name))
        call mem_setptr(newtonterm, 'NEWTONTERM',                              &
          trim(this%lxch(ixp)%exchange(ix)%name))
        call mem_setptr(moffset, 'MOFFSET',                                    &
          trim(this%lxch(ixp)%exchange(ix)%m1_name))
        !
        iv = ix
        if (vgbuf%rcvmt(iv)%memitype /= iadbl1d) then
          call store_error('Program error 1: mpi_mv_halo')
          call ustop()
        end if
        if (vgbuf%rcvmt(iv)%isize /= nexg) then
          call store_error('Program error 2: mpi_mv_halo')
          call ustop()
        end if
        do iexg = 1, nexg
          n = nodem1(iexg) + moffset
          v = vgbuf%rcvmt(iv)%adbl1d(iexg)
          if (active(n) > 0) then
            x(n) = x(n) + (cond(iexg) + newtonterm(iexg))*v
          end if
        end do
      end do 
    end do
    !
    ! -- return
    return
  end subroutine mpi_mv_halo
  
  function mpi_get_halo_rcvmt(this, vgname, varname, varext, hmname) result(mt)
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiExchangeType) :: this
    character(len=*), intent(in) :: vgname
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: varext
    character(len=*), intent(in) :: hmname
    type(MemoryType), pointer :: mt
    ! -- local
    type(MpiGwfBuf), pointer :: vgbuf
    type(ExchangeType), pointer :: ex
    type(VarGroupType), pointer :: vgvar
    type(VarType), pointer :: var
    logical :: lfound
    integer(I4B) :: ivg, ixp, iex, iv, is
! ------------------------------------------------------------------------------
    !
    mt => null()
    !
    if (serialrun) then
      return
    end if
    !
    ! -- find the variable group name
    ivg = ifind(this%vg, trim(vgname))
    if (ivg < 0) then
      call store_error('Program error: mpicopydbltohalo for '//trim(vgname))
      call ustop()
    end if
    !
    lfound = .false.
    do ixp = 1, this%nrxp
      vgbuf => this%lxch(ixp)%vgbuf(ivg)
      is = 0
      do iex = 1, this%lxch(ixp)%nexchange
        ex => this%lxch(ixp)%exchange(iex)
        vgvar => ex%vgvar(ivg)
        do iv = 1, vgvar%nvar
          var => vgvar%var(iv)
          if (var%lrcv) then
            is = is + 1
          end if
          if ((trim(hmname) == trim(ex%halo_name)) .and.                       &
              (trim(var%name) == trim(varname)) .and.                          &
              (trim(var%nameext) == trim(varext))) then
            mt => vgbuf%rcvmt(is)
            lfound = .true.
          end if
          if(lfound) exit
        end do
        if (lfound) exit
      end do
      if (lfound) exit
    end do
    !
    ! -- return
    return
  end function mpi_get_halo_rcvmt
  
  subroutine mpi_copy_int_to_halo(this, vgname, varname, varext, hmname,        &
    gwfhaloarray, offset)
! ******************************************************************************
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryTypeModule, only: iaint1d
    ! -- dummy
    class(MpiExchangeType) :: this
    character(len=*), intent(in) :: vgname
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: varext
    character(len=*), intent(in) :: hmname
    integer(I4B), dimension(:), intent(inout) :: gwfhaloarray
    integer(I4B), intent(in) :: offset
    ! -- local
    type(MemoryType), pointer :: mt
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    end if
    !
    mt => this%mpi_get_halo_rcvmt(vgname, varname, varext, hmname)
    !
    if (.not.associated(mt)) then
      call store_error('Program error 1: mpi_copy_int_to_halo')
      call ustop()
    end if
    !
    ! check
    if (mt%memitype /= iaint1d) then
      call store_error('Program error 2: mpi_copy_int_to_halo')
      call ustop()
    endif
    !
    do i = 1, mt%isize
      gwfhaloarray(i+offset) = mt%aint1d(i)
    end do
    !
    ! -- return
    return
  end subroutine mpi_copy_int_to_halo

  subroutine mpi_copy_dbl_to_halo(this, vgname, varname, varext, hmname,        &
    gwfhaloarray, offset)
! ******************************************************************************
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryTypeModule, only: iadbl1d
    ! -- dummy
    class(MpiExchangeType) :: this
    character(len=*), intent(in) :: vgname
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: varext
    character(len=*), intent(in) :: hmname
    real(DP), dimension(:), intent(inout) :: gwfhaloarray
    integer(I4B), intent(in) :: offset
    ! -- local
    type(MemoryType), pointer :: mt
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    end if
    !
    mt => this%mpi_get_halo_rcvmt(vgname, varname, varext, hmname)
    !
    if (.not.associated(mt)) then
      call store_error('Program error 1: mpi_copy_dbl_to_halo')
      call ustop()
    end if
    !
    ! check
    if (mt%memitype /= iadbl1d) then
      call store_error('Program error 2: mpi_copy_dbl_to_halo')
      call ustop()
    endif
    !
    do i = 1, mt%isize
      gwfhaloarray(i+offset) = mt%adbl1d(i)
    end do
    !
    ! -- return
    return
  end subroutine mpi_copy_dbl_to_halo

  subroutine mpi_global_exchange_sum1(this, dval)
! ******************************************************************************
! Collective sum over all processes for one double precision value.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiExchangeType) :: this
    real(DP), intent(inout) :: dval
    ! -- local
    real(DP), dimension(1) :: dwrk
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    end if
    !
    dwrk(1) = dval
#ifdef MPI_TIMER
    call code_timer(0, t, this%ttgsum)
#endif
    call mpiwrpallreduce(dwrk, 1, 'mpi_sum', this%comm)
#ifdef MPI_TIMER
    call code_timer(1, t, this%ttgsum)
#endif
    dval = dwrk(1)
    !
    ! -- return
    return
  end subroutine mpi_global_exchange_sum1

  subroutine mpi_global_exchange_sum2(this, dval1, dval2)
! ******************************************************************************
! Collective sum over all processes for two double precision values.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiExchangeType) :: this
    real(DP), intent(inout) :: dval1, dval2
    ! -- local
    real(DP), dimension(2) :: dwrk
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    end if
    !
    dwrk(1) = dval1
    dwrk(2) = dval2
#ifdef MPI_TIMER
    call code_timer(0, t, this%ttgsum)
#endif
    call mpiwrpallreduce(dwrk, 2, 'mpi_sum', this%comm)
#ifdef MPI_TIMER
    call code_timer(1, t, this%ttgsum)
#endif
    dval1 = dwrk(1)
    dval2 = dwrk(2)
    !
    ! -- return
    return
  end subroutine mpi_global_exchange_sum2
  
  subroutine mpi_global_exchange_absmax1(this, dval)
! ******************************************************************************
! Collective maximum over all processes for one double precision value.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiExchangeType) :: this
    real(DP), intent(inout) :: dval
    ! -- local
    real(DP), dimension(1) :: dwrk
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    end if
    !
    dwrk(1) = abs(dval)
#ifdef MPI_TIMER
    call code_timer(0, t, this%ttgmax)
#endif
    call mpiwrpallreduce(dwrk, 1, 'mpi_max', this%comm)
#ifdef MPI_TIMER
    call code_timer(0, t, this%ttgmax)
#endif
    dval = dwrk(1)
    !
    ! -- return
    return
  end subroutine mpi_global_exchange_absmax1
  
  subroutine mpi_global_exchange_absmax2(this, dval1, dval2)
! ******************************************************************************
! Collective maximum over all processes for two double precision values.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiExchangeType) :: this
    real(DP), intent(inout) :: dval1, dval2
    ! -- local
    real(DP), dimension(2) :: dwrk
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    end if
    !
    dwrk(1) = abs(dval1)
    dwrk(2) = abs(dval2)
#ifdef MPI_TIMER
    call code_timer(0, t, this%ttgmax)
#endif
    call mpiwrpallreduce(dwrk, 2, 'mpi_max', this%comm)
#ifdef MPI_TIMER
    call code_timer(1, t, this%ttgmax)
#endif
    dval1 = dwrk(1)
    dval2 = dwrk(2)
    !
    ! -- return
    return
  end subroutine mpi_global_exchange_absmax2
  
  subroutine mpi_global_exchange_max_int(this, ival)
! ******************************************************************************
! Collective maximum over all processes for one integer value.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiExchangeType) :: this
    integer(I4B), intent(inout) :: ival
    ! -- local
    integer(I4B), dimension(1) :: iwrk
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    end if
    !
    iwrk(1) = ival
#ifdef MPI_TIMER
    call code_timer(0, t, this%ttgmax)
#endif
    call mpiwrpallreduce(iwrk, 1, 'mpi_max', this%comm)
#ifdef MPI_TIMER
    call code_timer(1, t, this%ttgmax)
#endif
    ival = iwrk(1)
    !
    ! -- return
    return
  end subroutine mpi_global_exchange_max_int
  
  subroutine mpi_global_exchange_absmin1(this, dval)
! ******************************************************************************
! Collective minimum over all processes for one double precision value.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiExchangeType) :: this
    real(DP), intent(inout) :: dval
    ! -- local
    real(DP), dimension(1) :: dwrk
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    end if
    !
    dwrk(1) = abs(dval)
#ifdef MPI_TIMER
    call code_timer(0, t, this%ttgmin)
#endif
    call mpiwrpallreduce(dwrk, 1, 'mpi_min', this%comm)
#ifdef MPI_TIMER
    call code_timer(1, t, this%ttgmin)
#endif
    dval = dwrk(1)
    !
    ! -- return
    return
  end subroutine mpi_global_exchange_absmin1
  
  subroutine mpi_global_exchange_absmin2(this, dval1, dval2)
! ******************************************************************************
! Collective minimum over all processes for two double precision values.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiExchangeType) :: this
    real(DP), intent(inout) :: dval1, dval2
    ! -- local
    real(DP), dimension(2) :: dwrk
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    end if
    !
    dwrk(1) = abs(dval1)
    dwrk(2) = abs(dval2)
#ifdef MPI_TIMER
    call code_timer(0, t, this%ttgmin)
#endif
    call mpiwrpallreduce(dwrk, 2, 'mpi_min', this%comm)
#ifdef MPI_TIMER
    call code_timer(1, t, this%ttgmin)
#endif
    dval1 = dwrk(1)
    dval2 = dwrk(2)
    !
    ! -- return
    return
  end subroutine mpi_global_exchange_absmin2

  subroutine mpi_global_exchange_cgc(this, nla1d, nga1d, la1d, ga1d,            &
    gnia, gnja, gia, gja, iopt) !CGC
! ******************************************************************************
! Collective gather over all processes for coarse grid matrix
! coefficients.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiExchangeType) :: this
    integer(I4B), intent(in) :: nla1d
    integer(I4B), intent(in) :: nga1d
    real(DP), intent(in), dimension(nla1d)  :: la1d
    real(DP), intent(out), dimension(nga1d) :: ga1d
    integer(I4B), intent(in) :: gnia
    integer(I4B), intent(in) :: gnja
    integer(I4B), dimension(gnia), intent(in) :: gia
    integer(I4B), dimension(gnja), intent(in) :: gja
    integer(I4B), intent(in) :: iopt
    ! -- local
    integer(I4B) :: i, j, k, ip, im, ic0, ic1, m
    integer(I4B), dimension(:), allocatable :: recvcounts, displs
    real(I4B), dimension(:), allocatable :: iwrk
    real(DP), dimension(:), allocatable :: rbuf
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      do i = 1, nga1d
        ga1d(i) = la1d(i)
      end do
      return
    end if
    !
    ! -- allocate work arrays
    allocate(recvcounts(this%nrproc), displs(this%nrproc), iwrk(this%nrproc), rbuf(nga1d))
    !
    ! -- loop over global topology and determine number of coeffients to send
    do ip = 1, this%nrproc
      recvcounts(ip) = 0
    end do
    do i = 1, gnia-1 ! loop over all models
      ic0 = gia(i)
      if (iopt == 1) then
        ic1 = gia(i+1) - 1
      else
        ic1 = ic0
      end if  
      do j = ic0, ic1
        im = gja(j) ! global model
        ip = MpiWorld%gsubs(im) ! processor ID
        recvcounts(ip) = recvcounts(ip) + 1
      end do
    end do
    !
    ! -- determine displacements
    displs(1) = 0
    do i = 2, this%nrproc
      displs(i) = displs(i-1) + recvcounts(i-1)
    end do
    !
#ifdef MPI_TIMER
    call code_timer(0, t, this%ttgmax)
#endif
    ! - gather all matrix coefficients
    call mpiwrpallgatherv(this%comm, la1d, nla1d, rbuf, recvcounts, displs)
#ifdef MPI_TIMER
    call code_timer(1, t, this%ttgmax)
#endif
    !
    iwrk = 0
    m = 0
    do i = 1, gnia-1 ! loop over all models
      ic0 = gia(i)
      if (iopt == 1) then
        ic1 = gia(i+1) - 1
      else
        ic1 = ic0
      end if
      im = gja(ic0) ! global model
      ip = MpiWorld%gsubs(im) ! processor ID
      do j = ic0, ic1
        iwrk(ip) = iwrk(ip) + 1
        k = displs(ip) + iwrk(ip) ! receive buffer index
        m = m + 1
        ga1d(m) = rbuf(k)
      end do
    end do
    !
    ! -- cleanup
    deallocate(recvcounts, displs, rbuf, iwrk)
    !
    ! -- return
    return
  end subroutine mpi_global_exchange_cgc
    
  subroutine mpi_pack_mt(origin, mti, mto, node, nexg, moffset)
! ******************************************************************************
! Pack memory type for point-to-point communication.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryTypeModule, only: iintsclr, idblsclr, iaint1d, iadbl1d
    ! -- dummy
    character(len=*), intent(in) :: origin
    type(MemoryType), intent(in) :: mti
    type(MemoryType), intent(out) :: mto
    integer(I4B), intent(in), optional :: nexg
    integer(I4B), intent(in), optional :: moffset
    integer(I4B), dimension(:), intent(in), optional :: node
    ! -- local
    integer(I4B) :: i, n
! ------------------------------------------------------------------------------
    !
    write(errmsg,'(a)') 'Program error in mpi_pack_mt.'
    !
    mto%name     = mti%name
    mto%origin   = trim(origin)
    mto%memitype = mti%memitype
    mto%isize    = 0 
    !
    ! TODO: checks for K22, K33, ANGLE1, ANGLE2 and ANGLE3
    if (present(node)) then
      if ((mti%memitype == iaint1d) .or. (mti%memitype == iadbl1d)) then
        if (mti%isize < nexg) then
          return
        end if
      end if
    end if
    !
    select case(mti%memitype)
      case(iintsclr)
        allocate(mto%intsclr)
        mto%intsclr = mti%intsclr
      case(idblsclr)
        allocate(mto%dblsclr)
        mto%dblsclr = mti%dblsclr
      case(iaint1d)
        if (present(node)) then
          mto%isize = nexg
          allocate(mto%aint1d(mto%isize))
          do i = 1, mto%isize
            n = node(i) + moffset
            if ((n > 0) .and. (n <= size(mti%aint1d))) then
              mto%aint1d(i) = mti%aint1d(n)
            end if
          end do
        else
          mto%isize = mti%isize
          allocate(mto%aint1d(mto%isize))
          do i = 1, mto%isize
            mto%aint1d(i) = mti%aint1d(i)
          end do
        end if
        !write(*,*) '# int n, isize',trim(mti%name)//' '//trim(mti%origin), n,size(mti%aint1d)
      case(iadbl1d)
        if (present(node)) then
          mto%isize = nexg
          allocate(mto%adbl1d(mto%isize))
          do i = 1, mto%isize
            n = node(i) + moffset
            if ((n > 0) .and. (n <= size(mti%adbl1d))) then
              mto%adbl1d(i) = mti%adbl1d(n)
            end if
          end do
        else
          mto%isize = mti%isize
          allocate(mto%adbl1d(mto%isize))
          do i = 1, mto%isize
            mto%adbl1d(i) = mti%adbl1d(i)
          end do
        end if
        !write(*,*) '# dbl n, isize',trim(mti%name)//' '//trim(mti%origin), n,size(mti%adbl1d)
      case default
        call store_error(errmsg)
        call ustop()
    end select
    !
    ! -- return
    return
  end subroutine mpi_pack_mt
  
  subroutine mpi_pack_mmt(mt, mmt)
! ******************************************************************************
! Pack memory type for point-to-point communication.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryTypeModule, only: iaint2d, iadbl2d
    ! -- dummy
    type(MemoryType), intent(in) :: mt
    type(MetaMemoryType), intent(out) :: mmt
    ! -- local
! ------------------------------------------------------------------------------
    !
    mmt%name     = mt%name
    mmt%origin   = mt%origin
    mmt%memitype = mt%memitype
    mmt%isize    = mt%isize
    select case(mt%memitype)
    case(iaint2d)
        mmt%ncol = size(mt%aint2d, dim=1)
        mmt%nrow = size(mt%aint2d, dim=2)
      case(iadbl2d)
        mmt%ncol = size(mt%adbl2d, dim=1)
        mmt%nrow = size(mt%adbl2d, dim=2)
    end select
    !
    ! -- return
    return
  end subroutine mpi_pack_mmt

  subroutine mpi_is_iproc(this, isub, add)
! ******************************************************************************
! Add a model to my subdomain.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiExchangeType) :: this
    integer(I4B), intent(in)  :: isub
    logical, intent(out) :: add
    ! -- local
! ------------------------------------------------------------------------------
    if (serialrun) then
      !add = .true. !DEBUG 
      !return
    end if
    !
    if (isub == this%myproc) then
      add = .true.
    else
      add = .false.
    end if  
    !
    ! -- return
    return
  end subroutine mpi_is_iproc

  subroutine mpi_addmodel(this, iopt, mname)
! ******************************************************************************
! Add a model to my subdomain.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ArrayHandlersModule, only: ExpandArray
    use SimModule, only: store_error, ustop
    ! -- dummy
    class(MpiExchangeType) :: this
    integer(I4B), intent(in) :: iopt
    character(len=*), intent(in) :: mname
    ! -- local
    integer(I4B) :: n
! ------------------------------------------------------------------------------
    if (serialrun) then
       return
    end if
    !
    select case(iopt)
      ! store global
      case(1)
        call ExpandArray(this%gmodelnames)
        n = this%gnmodel
        n = n + 1
        this%gmodelnames(n) = mname
        this%gnmodel = n
      ! store local
      case(2)
        call ExpandArray(this%lmodelnames)
        n = this%lnmodel
        n = n + 1
        this%lmodelnames(n) = mname
        this%lnmodel = n
      ! store local
      case default
        write(errmsg,'(a)') 'Program error in mpi_addmodel.'
        call store_error(errmsg)
        call ustop()
    end select
    !
    ! -- return
    return
  end subroutine mpi_addmodel
  
  subroutine mpi_addhmodel(this, hmname, m1name, m2name)
! ******************************************************************************
! Add a model to my subdomain.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ArrayHandlersModule, only: ExpandArray
    use SimModule, only: store_error, ustop
    ! -- dummy
    class(MpiExchangeType) :: this
    character(len=*), intent(in) :: hmname, m1name, m2name
    ! -- local
    integer(I4B) :: n
! ------------------------------------------------------------------------------
    if (serialrun) then
      return
    end if
    !
    call ExpandArray(this%hmodelnames)
    call ExpandArray(this%hmodelm1names)
    call ExpandArray(this%hmodelm2names)
    n = this%hnmodel
    n = n + 1
    this%hmodelnames(n)   = hmname
    this%hmodelm1names(n) = m1name
    this%hmodelm2names(n) = m2name
    this%hnmodel = n
    !
    ! -- return
    return
  end subroutine mpi_addhmodel
  
  subroutine mpi_getmodel(this, mname, lok)
! ******************************************************************************
! Get the model name. In case this model turn out to be halo, return the first
! interface.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
    ! -- modules
    use ArrayHandlersModule, only: ExpandArray
    use SimModule, only: store_error, ustop
    ! -- dummy
    class(MpiExchangeType) :: this
    character(len=*), intent(inout) :: mname
    logical, intent(out) :: lok
    ! -- local
    type(ExchangeType), pointer :: ex
    integer(I4B) :: m, ixp, iex
! ------------------------------------------------------------------------------
    !
    lok = .false.
    if (serialrun) then
      lok = .true.
      return
    end if
    !
    m = ifind(this%gmodelnames, mname)
    if (m <= 0) then
      write(errmsg,'(a)') 'Error, mover model not found.'
      call store_error(errmsg)
      call ustop()
    end if
    !
    m = ifind(this%lmodelnames, mname)
      ! find the halo model
    if (m <= 0) then
      do ixp = 1, this%nrxp
        do iex = 1, this%lxch(ixp)%nexchange
          ex => this%lxch(ixp)%exchange(iex)
          if (trim(ex%m2_name) == trim(mname)) then
            mname = ex%halo_name
            lok = .true.
          end if
        end do
      end do
    else
      lok = .true.
    end if
    !
    ! -- return
    return
  end subroutine mpi_getmodel
  
  subroutine mpi_addsub(this, iopt, isub)
! ******************************************************************************
! Add a model to my subdomain.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ArrayHandlersModule, only: ExpandArray
    use SimModule, only: store_error, ustop
    ! -- dummy
    class(MpiExchangeType) :: this
    integer(I4B), intent(in) :: iopt
    integer(I4B), intent(in) :: isub
    ! -- local
    integer(I4B) :: n
! ------------------------------------------------------------------------------
    if (serialrun) then
      return
    end if
    !
    select case(iopt)
      ! store global
      case(1)
        call ExpandArray(this%gsubs)
        n = this%gnsub
        n = n + 1
        this%gsubs(n) = isub
        this%gnsub = n
      ! store local
      case(2)
        call ExpandArray(this%lsubs)
        n = this%lnsub
        n = n + 1
        this%lsubs(n) = isub
        this%lnsub = n
      case default
        write(errmsg,'(a)') 'Program error in mpi_addsub.'
        call store_error(errmsg)
        call ustop()
    end select
    !
    ! -- return
    return
  end subroutine mpi_addsub
      
!  subroutine mpi_add_halo_model(im, modelname)
!! ******************************************************************************
!! This subroutine sets the list of halo (m2) models
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    use MpiExchangeGenModule, only: nhalo, modelname_halo,                      &
!                                    mpi_create_modelname_halo
!    ! -- dummy
!    integer, intent(in) :: im
!    character(len=*), intent(inout) :: modelname
!    ! -- local
!    integer(I4B) :: m
!! ------------------------------------------------------------------------------
!    call mpi_create_modelname_halo(im, modelname)
!    m = ifind(modelname_halo, modelname)
!    if (m < 0) then
!      nhalo = nhalo + 1
!      call ExpandArray(modelname_halo)
!      modelname_halo(nhalo) = modelname
!    end if
!    !
!    ! -- return
!    return
!  end subroutine mpi_add_halo_model
  
  subroutine mpi_to_colmem(mt, is, cmt, iopt)
! ******************************************************************************
! Convert to collective MemoryType
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ConstantsModule, only: DZERO  
    use MemoryTypeModule, only: MemoryType
    ! -- dummy
    type(MemoryType), intent(in) :: mt
    integer(I4B), intent(inout) :: is
    type(ColMemoryType), dimension(*), intent(inout) :: cmt
    integer(I4B), intent(in) :: iopt
    ! -- local
    integer(I4B) :: i
! ------------------------------------------------------------------------------
    is = is + 1
    if (iopt == 1) then
      return
    end if
    
    cmt(is)%name     = mt%name
    cmt(is)%origin   = mt%origin
    cmt(is)%memitype = mt%memitype
    if (associated(mt%logicalsclr)) then
      cmt(is)%logicalsclr = mt%logicalsclr
    else
      cmt(is)%logicalsclr = .false.
    end if  
    if (associated(mt%intsclr)) then
      cmt(is)%intsclr = mt%intsclr
    else
      cmt(is)%intsclr = 0
    end if
    if (associated(mt%dblsclr)) then
      cmt(is)%dblsclr = mt%dblsclr
    else
      cmt(is)%dblsclr = DZERO
    end if
    if (associated(mt%aint1d)) then
      if (mt%isize > size(cmt(is)%aint1d)) then
        call store_error('Program error: mpi_to_colmem')
        call ustop()
      end if
      do i = 1, mt%isize
        cmt(is)%aint1d(i) = mt%aint1d(i)
      end do
    else
      cmt(is)%aint1d = 0
    end if
    ! -- return
    return
  end subroutine mpi_to_colmem
  
  subroutine mpi_not_supported(this, msg)
! ******************************************************************************
! Stop in case the functionality is not supported in parallel.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use SimModule, only: store_error, ustop
    ! -- dummy
    class(MpiExchangeType) :: this
    character(len=*), intent(in) :: msg
    ! -- local
! ------------------------------------------------------------------------------
    if (serialrun) then
      return
    end if
    !
    errmsg = 'Not (yet) supported in parallel: '//trim(msg)//'. Stopping.'
    call ustop(errmsg,0,writestd)
    !
    ! -- return
    return
  end subroutine mpi_not_supported
  
  subroutine mpi_da(this)
! ******************************************************************************
! mpi_da -- deallocate
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(MpiExchangeType) :: this
    ! -- local
    type(MpiGwfBuf), pointer :: vgbuf
    type(ExchangeType), pointer :: ex
    integer(I4B) :: ivg, ixp, iex
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    end if
    
    ! -- communicator
    if (associated(this%comm)) then
      call mem_deallocate(this%comm)
    end if
    if (associated(this%nrproc)) then
      call mem_deallocate(this%nrproc)
    end if  
    if (associated(this%myrank)) then
      call mem_deallocate(this%myrank)
    end if
    if (associated(this%myproc)) then
      call mem_deallocate(this%myproc)
    end if
    if(associated(this%procmap)) then
      deallocate(this%procmap)
    end if
    !
    ! -- model lists
    if(allocated(this%gmodelnames)) then
      deallocate(this%gmodelnames)
    end if
    if(allocated(this%lmodelnames)) then
      deallocate(this%lmodelnames)
    end if
    if(allocated(this%hmodelnames)) then
      deallocate(this%hmodelnames)
    end if
    if(allocated(this%gsubs)) then
      deallocate(this%gsubs)
    end if
    if(allocated(this%lsubs)) then
      deallocate(this%lsubs)
    end if
    !
    ! -- point-to-point
    if (associated(this%nrxp)) then
      do ixp = 1, this%nrxp
        do ivg = 1, this%nvg
          vgbuf => this%lxch(ixp)%vgbuf(ivg)
          if(associated(vgbuf%sndmt)) then
            deallocate(vgbuf%sndmt)
          end if
          if(associated(vgbuf%rcvmt)) then
            deallocate(vgbuf%rcvmt)
          end if
          if(associated(vgbuf%sndmmt)) then
            deallocate(vgbuf%sndmmt)
          end if
          if(associated(vgbuf%rcvmmt)) then
            deallocate(vgbuf%rcvmmt)
          end if
        end do
        do iex = 1, this%lxch(ixp)%nexchange
          ex => this%lxch(ixp)%exchange(iex)
          if (associated(ex%vgvar)) then
            deallocate(ex%vgvar)
          end if
        end do
        deallocate(this%lxch(ixp)%xprnk)
        deallocate(this%lxch(ixp)%exchange)
        deallocate(this%lxch(ixp)%nexchange)
        deallocate(this%lxch(ixp)%vgbuf)
      end do
    end if
    if (associated(this%lxch)) then
      deallocate(this%lxch)
    end if
    if (associated(this%nrxp)) then
      call mem_deallocate(this%nrxp)
    end if
    if(associated(this%nrprocstr)) then
      deallocate(this%nrprocstr)
    end if
    if(associated(this%npdigits)) then
      call mem_deallocate(this%npdigits)
    end if
    if(associated(this%partstr)) then
      deallocate(this%partstr)
    end if
    !
    ! -- Initialize constants
    this%name = ''
    this%linit = .false.
    this%lnmodel = 0
    this%gnmodel = 0
    this%gnsub = 0
    this%lnsub = 0
    this%lp2p = .true.
    this%nvg = 0
    this%vg = ''
    this%lxchmeta = .true.
    !
    ! -- return
    return
  end subroutine mpi_da
end module MpiExchangeModule