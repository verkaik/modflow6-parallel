module MpiExchangeModule
  
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
  public :: ivmtmvr
  ! -- Public functions
  public :: mpi_initialize_world
  public :: mpi_world_da
  public :: mpi_add_halo_model
  public :: mpi_to_colmem
  
  integer, parameter :: MAXNVG  = 10  ! maximum number of variable groups 
  integer, parameter :: MAXNVAR = 100
  integer, parameter :: MAXNEX  = 5
  integer, parameter :: ivmtsol = 1
  integer, parameter :: ivmtgwf = 2
  integer, parameter :: ivmtmvr = 3
  
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
    integer(I4B)              :: vmttype
    integer(I4B), dimension(1) :: id = 0
  end type Vartype

  type VarGroupType
    integer(I4B)                      :: nvar = 0
    type(VarType), dimension(MAXNVAR) :: var
  end type VarGroupType
  
  type ExchangeType
    type(VarGroupType), dimension(:), pointer :: vgvar => null()
    character(len=LENVARNAME)                 :: name
    character(len=LENMODELNAME)               :: m1_name
    character(len=LENMODELNAME)               :: m2_name
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
    integer(I4B)                                           :: gnmodel = 0   ! number of global models
    character(len=LENMODELNAME), dimension(:), allocatable :: gmodelnames   ! global model names
    integer(I4B)                                           :: gnsub = 0     ! number of global subdomains
    integer(I4B), dimension(:), allocatable                :: gsubs         ! global subdomains
    integer(I4B)                                           :: lnmodel = 0   ! number of local models
    character(len=LENMODELNAME), dimension(:), allocatable :: lmodelnames   ! local model names
    integer(I4B)                                           :: lnsub = 0     ! number of local subdomains
    integer(I4B), dimension(:), allocatable                :: lsubs         ! model local subdomains
    integer(I4B), pointer                                  :: comm      => null() ! MPI communicator
    integer(I4B), pointer                                  :: nrproc    => null() ! number of MPI process for this communicator
    integer(I4B), pointer                                  :: myrank    => null() ! MPI rank in for this communicator
    integer(I4B), pointer                                  :: myproc    => null() ! MPI proc in for this communicator
    integer(I4B), dimension(:), pointer                    :: procmap   => null() ! Mapping
    integer(I4B), pointer                                  :: nrxp => null() ! Number of exchange partners
    type(MpiGwfCommInt), dimension(:), pointer             :: lxch => null() ! Point to-point communication structure
    integer(I4B)                                           :: nvg  ! number of variable groups
    character(len=LINELENGTH), dimension(MAXNVG)           :: vg   ! variable groups
    logical,  dimension(MAXNVG)                            :: lxchmeta = .true. ! exchange meta data 
    character(len=50), pointer                             :: nrprocstr => null() ! Number of processes string
    integer(I4B), pointer                                  :: npdigits  => null() ! Number of digits for nrproc
    character(len=50), pointer                             :: partstr   => null() ! Partition string
    logical                                                :: lp2p = .true.
  contains
    procedure :: mpi_barrier
    procedure :: mpi_create_output_str
    procedure :: mpi_append_fname
    procedure :: mpi_is_iproc
    procedure :: mpi_addmodel
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
    generic, public :: mpi_global_exchange_sum => mpi_global_exchange_sum1, mpi_global_exchange_sum2
    procedure :: mpi_global_exchange_absmin1
    procedure :: mpi_global_exchange_absmin2
    generic, public :: mpi_global_exchange_absmin => mpi_global_exchange_absmin1, mpi_global_exchange_absmin2
    procedure :: mpi_global_exchange_absmax1
    procedure :: mpi_global_exchange_absmax2
    generic, public :: mpi_global_exchange_absmax => mpi_global_exchange_absmax1, mpi_global_exchange_absmax2
    procedure :: mpi_global_exchange_max => mpi_global_exchange_max_int
    procedure :: mpi_debug
    procedure :: mpi_da
    procedure :: mpi_not_supported
  end type MpiExchangeType
  
  ! -- World communicator
  type(MpiExchangeType), pointer :: MpiWorld => null()

  character(len=LINELENGTH) :: errmsg
  
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
    ! -- dummy
    ! -- local
    character(len=LENORIGIN) :: origin
    integer(I4B) :: memitype, isize
! ------------------------------------------------------------------------------
    !
    ! -- Allocate MpiWorld object
    allocate(MpiWorld)
    !
    ! -- Set name
    MpiWorld%name = 'MPI_WORLD'
    !
    ! -- Allocate scalars
    origin = MpiWorld%name
    call mem_allocate(MpiWorld%comm, 'COMM', origin)
    call mem_allocate(MpiWorld%nrproc, 'NRPROC', origin)
    call mem_allocate(MpiWorld%myrank, 'MYRANK', origin)
    call mem_allocate(MpiWorld%myproc, 'MYPROC', origin)
    !
    MpiWorld%comm   = mpiwrpcommworld()
    MpiWorld%nrproc = mpiwrpcomm_size(MpiWorld%comm)
    MpiWorld%myrank = mpiwrpcomm_rank(MpiWorld%comm)
    MpiWorld%myproc = MpiWorld%myrank + 1
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
    call mpiwrpbarrier(this%comm)
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
  
  subroutine mpi_append_fname(this,f)
! ******************************************************************************
! Append the file name with the process rank ID.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiExchangeType) :: this
    character(len=*), intent(inout) :: f
    ! -- local
! ------------------------------------------------------------------------------
    if (serialrun) then
      return
    end if
    !
    write(f,'(3a)') trim(f),'.',trim(this%partstr)
    !
    ! -- return
    return
  end subroutine mpi_append_fname
  
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
        !write(*,*) '@@@@@ allocating vgbuf for '//trim(vgname),ixp
        allocate(this%lxch(ixp)%vgbuf(MAXNVG))
      end if
      do iex = 1, this%lxch(ixp)%nexchange
        if (.not.associated(this%lxch(ixp)%exchange(iex)%vgvar)) then
          !write(*,*) '@@@@@ allocating vgvar for '//trim(vgname),iex
          allocate(this%lxch(ixp)%exchange(iex)%vgvar(MAXNVG))
        end if
      end do
    end do
    !
    ! -- return
    return
  end subroutine mpi_add_vg

  subroutine mpi_add_vmt(this, vgname, origin, name, nameext, vmttype)
! ******************************************************************************
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiExchangeType) :: this
    character(len=*), intent(in) :: vgname
    character(len=*), intent(in) :: origin
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: nameext
    character(len=*), intent(in) :: vmttype
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
      call store_error('Program error: mpi_add_vmt')
      call ustop()
    end if
    !
    ivg = ifind(this%vg, trim(vgname))
    if (ivg < 0) then
      write(errmsg,'(a)') 'Program error: mpi_add_vmt'
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
          call store_error('Program error: mpi_add_vmt')
          call ustop()
        end if
        vgvar%nvar = iv
        vgvar%var(iv)%origin  = trim(origin)
        vgvar%var(iv)%name    = trim(name)
        vgvar%var(iv)%nameext = trim(nameext)
        select case(trim(vmttype))
          case('SOL')
            vgvar%var(iv)%vmttype = ivmtsol
          case('GWF')
            vgvar%var(iv)%vmttype = ivmtgwf
          case('MVR')        
            vgvar%var(iv)%vmttype = ivmtmvr
          case default
            call store_error('Program error: mpi_add_vmt')
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
    use MemoryTypeModule, only: iaint1d, iadbl1d !@@@@DEBUG
    use MemoryManagerModule, only: mem_get_ptr, mem_setptr, mem_setval,        &
                                   mem_setval_id, mem_check_by_name !@@@@
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
    integer(I4B) :: ixp, iex, nsnd, nrcv, iv, is, i, j, istat
    character(len=LENMODELNAME) :: mname, id, m1_name, m2_name
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
          ! -- Loop over the exchange send variables
          do iv = 1, vgvar%nvar
            var => vgvar%var(iv)
            if (var%lsnd) then
              is = is + 1
              ! -- check
              if (is > vgbuf%nsnd) then
                !write(*,*) '@@@ packing '//trim(vgname)//' '//trim(var%name),is
                call store_error('Program error: mpi_local_exchange')
                call ustop()
              end if
              ! -- Set the origin and offset
              write(mod_origin,'(a,1x,a)') trim(m1_name), trim(var%nameext)
              select case(var%vmttype)
                case(ivmtgwf)
                  src_origin = mod_origin
                  moffset = 0
                case(ivmtsol)
                  read(solname,*) sol_origin
                  write(sol_origin,'(a,1x,a)') trim(sol_origin), trim(var%nameext)
                  src_origin = sol_origin
                  call mem_setptr(tmp, 'MOFFSET', trim(m1_name))
                  moffset = tmp
                case(ivmtmvr)
                  src_origin = var%origin
                  tgt_origin = src_origin
                  moffset = 0
              end select
              !  
              call mem_get_ptr(var%name, src_origin, mt)
              !
              select case(var%vmttype)
                case(ivmtgwf,ivmtsol)
                  tgt_origin = mod_origin
                  call mem_setptr(tmp, 'NEXG', trim(ex%name))
                  call mem_setptr(nodem1, 'NODEM1', trim(ex%name))
                  if (.not.associated(tmp).or..not.associated(nodem1)) then
                    call store_error('Program error: mpi_local_exchange')
                    call ustop()
                  end if
                  nexg = tmp
                case(ivmtmvr)
                  nexg = 1
                  nodem1 => var%id
                end select
                call mpi_pack_mt_nodem1(tgt_origin, nodem1, nexg, moffset, mt,  &
                                       vgbuf%sndmt(is))
            end if
          end do
        end do
      end do
    end if
    !
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
      do irank = 0, this%nrproc-1
        if (irank == this%myrank) then
          write(*,*) '=================myrank',this%myrank
          do ixp = 1, this%nrxp
            vgbuf => this%lxch(ixp)%vgbuf(ivg)
            write(*,*) '---sent to',this%lxch(ixp)%xprnk
            do i = 1, vgbuf%nsnd
              write(*,*) 'name:      ', trim(vgbuf%sndmt(i)%name)
              write(*,*) 'origin:    ', trim(vgbuf%sndmt(i)%origin)
              write(*,*) 'isize:     ', vgbuf%sndmt(i)%isize
              memitype = vgbuf%sndmt(i)%memitype
              !write(*,*) 'memitype:  ', memitype
              select case(memitype)
                case(iadbl1d)
!                  do j = 1, vgbuf%sndmt(i)%isize
                  do j = 1, 1
                    write(*,*) 'dval:      ', vgbuf%sndmt(i)%adbl1d(j)
                  end do
                case(iaint1d)
!                  do j = 1, vgbuf%sndmt(i)%isize
                  do j = 1, 1
                    write(*,*) 'ival:      ', vgbuf%sndmt(i)%aint1d(j)
                  end do
              end select
            end do
            write(*,*) '---received from',this%lxch(ixp)%xprnk
            do i = 1,vgbuf%nrcv
              write(*,*) 'name:      ', trim(vgbuf%rcvmt(i)%name)
              write(*,*) 'origin:    ', trim(vgbuf%rcvmt(i)%origin)
              write(*,*) 'isize:     ', vgbuf%rcvmt(i)%isize
              !write(*,*) 'memitype:  ', vgbuf%rcvmt(i)%memitype
              memitype = vgbuf%rcvmt(i)%memitype
              !write(*,*) 'memitype:  ', memitype
              select case(memitype)
                case(iadbl1d)
!                  do j = 1, vgbuf%sndmt(i)%isize
                  do j = 1, 1
                    write(*,*) 'dval:      ', vgbuf%rcvmt(i)%adbl1d(j)
                  end do
                case(iaint1d)
!                  do j = 1, vgbuf%rcvmt(i)%isize
                  do j = 1, 1
                    write(*,*) 'ival:      ', vgbuf%rcvmt(i)%aint1d(j)
                  end do
               end select
            end do
          end do
        end if
        call mpiwrpbarrier(this%comm)
      end do
    end if
    !
    ! -- Set the received data for the halo (m2) models
    if (lunpack) then
      do ixp = 1, this%nrxp
        vgbuf => this%lxch(ixp)%vgbuf(ivg)
        is = 0
        do iex = 1, this%lxch(ixp)%nexchange
          ex => this%lxch(ixp)%exchange(iex)
          ! -- halo model name
          m2_name = ex%m2_name
          vgvar => ex%vgvar(ivg)
          !write(*,*) '@@@@@@@ vgvar%nvar =',vgvar%nvar, this%myrank
          do iv = 1, vgvar%nvar
            var => vgvar%var(iv)
            if (var%lrcv) then
              is = is + 1
              ! check
              if (vgbuf%rcvmt(is)%name /= var%name) then
                call store_error('Program error: mpi_local_exchange')
                call ustop()
              end if
              select case(var%vmttype)
                case(ivmtgwf,ivmtsol)
                  src_origin = vgbuf%rcvmt(is)%origin
                  read(src_origin,*,iostat=istat) mname, id
                  if (istat /= 0) then
                    read(src_origin,*,iostat=istat) mname
                    id = ''
                  end if
                  if (index(id,'HALO') > 0) then
                    id = ''
                  end if
                  write(tgt_origin, '(a,1x,a)') trim(m2_name), trim(id)
                  vgbuf%rcvmt(is)%origin = tgt_origin
                  !write(*,*) '@@@ setting '//'"'//trim(vgbuf%rcvmt(is)%origin),'" "',trim(vgbuf%rcvmt(is)%name)
                  call mem_setval(vgbuf%rcvmt(is))
                case(ivmtmvr)
                  tgt_origin = var%origin
                  !if (this%myrank == 0) then
                  !  write(*,*) '@@@ setting MPI '//'"'//trim(tgt_origin), &
                  !    '" "',trim(vgbuf%rcvmt(is)%name),'": ',var%id, vgbuf%rcvmt(is)%adbl1d(1)
                  !end if
                  vgbuf%rcvmt(is)%origin = tgt_origin
                  call mem_setval_id(vgbuf%rcvmt(is), var%id, 1)
              end select
            end if
          end do
        end do
      end do   
    end if
    !
    ! -- return
    return
  end subroutine mpi_local_exchange
  
  subroutine mpi_mv_halo(this, solname, vgname, x)
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
    ! -- local
    type(MpiGwfBuf), pointer :: vgbuf
    character(len=LENORIGIN) :: sol_origin
    integer(I4B) :: ivg, ixp, ix, iv, iexg, n, i
    integer(I4B), pointer :: nexg
    integer(I4B), dimension(:), pointer :: nodem1
    integer(I4B), dimension(:), pointer :: active
    integer(I4B), pointer :: moffset
    real(DP), dimension(:), pointer :: cond
    real(DP), dimension(:), pointer :: newtonterm
    real(DP) :: v
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    end if
    !
    ! -- find the variable group name
    ivg = ifind(this%vg, trim(vgname))
    if (ivg < 0) then
      call store_error('Program error: mpi_local_exchange for '//trim(vgname))
      call ustop()
    end if
    !
    read(solname,*) sol_origin
    call mem_setptr(active, 'IACTIVE', trim(sol_origin))
    !
    do ixp = 1, this%nrxp
      vgbuf => this%lxch(ixp)%vgbuf(ivg)
      do ix = 1, this%lxch(ixp)%nexchange
        ! -- get exchange nodes and conductance
        call mem_setptr(nexg, 'NEXG', trim(this%lxch(ixp)%exchange(ix)%name))
        call mem_setptr(nodem1, 'NODEM1', trim(this%lxch(ixp)%exchange(ix)%name))
        call mem_setptr(cond, 'COND', trim(this%lxch(ixp)%exchange(ix)%name))
        call mem_setptr(newtonterm, 'NEWTONTERM', trim(this%lxch(ixp)%exchange(ix)%name))
        call mem_setptr(moffset, 'MOFFSET', trim(this%lxch(ixp)%exchange(ix)%m1_name))
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
            x(n) = x(n) + (cond(iexg) - newtonterm(iexg))*v
          end if
        end do
      end do 
    end do
    !
    ! -- return
    return
  end subroutine mpi_mv_halo
  
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
    call mpiwrpallreduce(dwrk, 1, 'mpi_sum', this%comm)
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
    call mpiwrpallreduce(dwrk, 2, 'mpi_sum', this%comm)
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
    call mpiwrpallreduce(dwrk, 1, 'mpi_max', this%comm)
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
    call mpiwrpallreduce(dwrk, 2, 'mpi_max', this%comm)
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
    call mpiwrpallreduce(iwrk, 1, 'mpi_max', this%comm)
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
    call mpiwrpallreduce(dwrk, 1, 'mpi_min', this%comm)
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
    call mpiwrpallreduce(dwrk, 2, 'mpi_min', this%comm)
    dval1 = dwrk(1)
    dval2 = dwrk(2)
    !
    ! -- return
    return
  end subroutine mpi_global_exchange_absmin2
  
  subroutine mpi_pack_mt_nodem1(origin, node, nexg, moffset, mti, mto)
! ******************************************************************************
! Pack memory type for point-to-point communication.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryTypeModule, only: iaint1d, iadbl1d
    ! -- dummy
    character(len=*), intent(in) :: origin
    integer(I4B), intent(in) :: nexg
    integer(I4B), intent(in) :: moffset
    integer(I4B), dimension(nexg), intent(in) :: node
    type(MemoryType), intent(in) :: mti
    type(MemoryType), intent(out) :: mto
    ! -- local
    integer(I4B) :: i, n
! ------------------------------------------------------------------------------
    !
    write(errmsg,'(a)') 'Program error in mpi_pack_mt_nodem1.'
    !
    mto%name     = mti%name
    mto%origin   = trim(origin)
    mto%memitype = mti%memitype
    mto%isize    = nexg
    !
    select case(mti%memitype)
      case(iaint1d)
        allocate(mto%aint1d(nexg))
        do i = 1, nexg
          n = node(i) + moffset
          if (n < 0 .or. n > size(mti%aint1d)) then
            call store_error(errmsg)
            call ustop()
          end if
          mto%aint1d(i) = mti%aint1d(n)
        end do
        !write(*,*) '# int n, isize',trim(mti%name)//' '//trim(mti%origin), n,size(mti%aint1d)
      case(iadbl1d)
        allocate(mto%adbl1d(nexg))
        do i = 1, nexg
          n = node(i) + moffset
          if (n < 0 .or. n > size(mti%adbl1d)) then
            call store_error(errmsg)
            call ustop()
          end if
          mto%adbl1d(i) = mti%adbl1d(n)
        end do
        !write(*,*) '# dbl n, isize',trim(mti%name)//' '//trim(mti%origin), n,size(mti%adbl1d)
      case default
        call store_error(errmsg)
        call ustop()
    end select
    !
    ! -- return
    return
  end subroutine mpi_pack_mt_nodem1
  
  subroutine mpi_pack_mt(origin, mti, mto)
! ******************************************************************************
! Pack memory type for point-to-point communication.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryTypeModule, only: iintsclr, idblsclr
    ! -- dummy
    character(len=*), intent(in) :: origin
    type(MemoryType), intent(in) :: mti
    type(MemoryType), intent(out) :: mto
    ! -- local
! ------------------------------------------------------------------------------
    !
    write(errmsg,'(a)') 'Program error in mpi_pack_mt.'
    !
    mto%name     = mti%name
    mto%origin   = trim(origin)
    mto%memitype = mti%memitype
    mto%isize    = mti%isize
    !
    select case(mti%memitype)
      case(iintsclr)
        allocate(mto%intsclr)
        mto%intsclr = mti%intsclr
      case(idblsclr)
        allocate(mto%dblsclr)
        mto%dblsclr = mti%dblsclr
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
      case default
        write(errmsg,'(a)') 'Program error in mpi_addmodel.'
        call store_error(errmsg)
        call ustop()
    end select
    !
    ! -- return
    return
  end subroutine mpi_addmodel
  
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
    use MpiExchangeGenModule, only: nhalo, modelname_halo,                      &
                                    mpi_create_modelname_halo
    ! -- dummy
    class(MpiExchangeType) :: this
    character(len=*), intent(inout) :: mname
    logical, intent(out) :: lok
    ! -- local
    integer(I4B) :: m
! ------------------------------------------------------------------------------
    !
    lok = .true.
    if (serialrun) then
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
    if (m <= 0) then
      call mpi_create_modelname_halo(1, mname)
      m = ifind(modelname_halo, mname)
      if (m <= 0) then
        lok = .false.
      end if
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
      
  subroutine mpi_add_halo_model(im, modelname)
! ******************************************************************************
! This subroutine sets the list of halo (m2) models
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    use MpiExchangeGenModule, only: nhalo, modelname_halo,                      &
                                    mpi_create_modelname_halo
    ! -- dummy
    integer, intent(in) :: im
    character(len=*), intent(inout) :: modelname
    ! -- local
    integer(I4B) :: m
! ------------------------------------------------------------------------------
    call mpi_create_modelname_halo(im, modelname)
    m = ifind(modelname_halo, modelname)
    if (m < 0) then
      nhalo = nhalo + 1
      call ExpandArray(modelname_halo)
      modelname_halo(nhalo) = modelname
    end if
    !
    ! -- return
    return
  end subroutine mpi_add_halo_model
  
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
      cmt(is)%aint1d = mt%aint1d
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
    type(VarGroupType), pointer :: vgvar
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