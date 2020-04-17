module GwfGwfExchangeModule

  use KindModule,              only: DP, I4B
  use ConstantsModule,         only: DZERO !HALO2
  use ArrayHandlersModule,     only: ExpandArray
  use BaseModelModule,         only: GetBaseModelFromList
  use BaseExchangeModule,      only: BaseExchangeType, AddBaseExchangeToList
  use ConstantsModule,         only: LENBOUNDNAME, NAMEDBOUNDFLAG, LINELENGTH, &
                                     TABCENTER, TABLEFT,                       &
                                     LENMODELNAME !HALO2
  use ListsModule,             only: basemodellist
  use NumericalExchangeModule, only: NumericalExchangeType
  use NumericalModelModule,    only: NumericalModelType
  use GwfModule,               only: GwfModelType
  use GwfHaloModule,           only: GwfHaloModelType !HALO2
  use GhostNodeModule,         only: GhostNodeType
  use GwfMvrModule,            only: GwfMvrType
  use ObserveModule,           only: ObserveType
  use ObsModule,               only: ObsType
  use SimModule,               only: count_errors, store_error,                &
                                     store_error_unit, ustop
  use BlockParserModule,       only: BlockParserType
  use TableModule,             only: TableType, table_cr

  implicit none

  private
  public :: gwfexchange_create
  public :: GwfExchangeType !HALO2
  public :: gwf_mpi_halo_init !HALO2

  type, extends(NumericalExchangeType) :: GwfExchangeType
    type(GwfHaloModelType), pointer                  :: gwfhalo     => null()    ! pointer to halo model !HALO2
    integer(I4B), pointer                            :: m1id        => null()    !HALO2
    integer(I4B), pointer                            :: m2id        => null()    !HALO2
    type(GwfModelType), pointer                      :: gwfmodel1   => null()    ! pointer to GWF Model 1
    type(GwfModelType), pointer                      :: gwfmodel2   => null()    ! pointer to GWF Model 2
    integer(I4B), pointer                            :: inewton     => null()    ! newton flag (1 newton is on)
    integer(I4B), pointer                            :: icellavg    => null()    ! cell averaging
    integer(I4B), pointer                            :: ivarcv      => null()    ! variable cv
    integer(I4B), pointer                            :: idewatcv    => null()    ! dewatered cv
    integer(I4B), pointer                            :: ianglex     => null()    ! flag indicating anglex was read, if read, ianglex is index in auxvar
    integer(I4B), pointer                            :: icdist      => null()    ! flag indicating cdist was read, if read, icdist is index in auxvar
    integer(I4B), pointer                            :: inamedbound => null()    ! flag to read boundnames
    real(DP), pointer                                :: satomega    => null()    ! saturation smoothing
    integer(I4B), dimension(:), pointer, contiguous  :: ihc         => null()    ! horizontal connection indicator array
    real(DP), dimension(:), pointer, contiguous      :: condsat     => null()    ! saturated conductance
    real(DP), dimension(:), pointer, contiguous      :: cl1         => null()    ! connection length 1
    real(DP), dimension(:), pointer, contiguous      :: cl2         => null()    ! connection length 2
    real(DP), dimension(:), pointer, contiguous      :: hwva        => null()    ! horizontal widths, vertical flow areas
    integer(I4B), pointer                            :: ingnc       => null()    ! unit number for gnc (0 if off)
    type(GhostNodeType), pointer                     :: gnc         => null()    ! gnc object
    integer(I4B), pointer                            :: inmvr       => null()    ! unit number for mover (0 if off)
    type(GwfMvrType), pointer                        :: mvr         => null()    ! water mover object
    integer(I4B), pointer                            :: inobs       => null()    ! unit number for GWF-GWF observations
    type(ObsType), pointer                           :: obs         => null()    ! observation object
    character(len=LENBOUNDNAME), dimension(:),                                  &
                                 pointer, contiguous :: boundname   => null()    ! boundnames
    !
    ! -- table objects
    type(TableType), pointer :: outputtab1 => null()
    type(TableType), pointer :: outputtab2 => null()

  contains

    procedure          :: exg_df      => gwf_gwf_df
    procedure          :: exg_ac      => gwf_gwf_ac
    procedure          :: exg_mc      => gwf_gwf_mc
    procedure          :: exg_ar      => gwf_gwf_ar
    procedure          :: exg_rp      => gwf_gwf_rp
    procedure          :: exg_ad      => gwf_gwf_ad
    procedure          :: exg_cf      => gwf_gwf_cf
    procedure          :: exg_fc      => gwf_gwf_fc
    procedure          :: exg_fn      => gwf_gwf_fn
    procedure          :: exg_cq      => gwf_gwf_cq
    procedure          :: exg_bd      => gwf_gwf_bd
    procedure          :: exg_ot      => gwf_gwf_ot
    procedure          :: exg_da      => gwf_gwf_da
    procedure          :: exg_fp      => gwf_gwf_fp
    procedure          :: get_iasym   => gwf_gwf_get_iasym
    procedure          :: get_m1m2    => gwf_gwf_get_m1m2 !PAR
    procedure          :: allocate_scalars
    procedure          :: allocate_arrays
    procedure          :: read_options
    procedure          :: read_data
    procedure          :: read_gnc
    procedure          :: read_mvr
    procedure, private :: rewet
    procedure, private :: qcalc
    procedure, private :: gwf_gwf_df_obs
    procedure, private :: gwf_gwf_rp_obs
    procedure, public  :: gwf_gwf_save_simvals
  end type GwfExchangeType

contains

  subroutine gwfexchange_create(filename, id, m1i, m2i, mname1i, mname2i, m2_bympi) !PAR
! ******************************************************************************
! Create a new GWF to GWF exchange object.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ConstantsModule, only: LINELENGTH
    use BaseModelModule, only: BaseModelType, AddBaseModelToList !HALO2
    use ListsModule, only: baseexchangelist, halomodellist !HALO2
    use ObsModule, only: obs_cr
    use GwfHaloModule, only: gwfhalo_cr !HALO2
    use MpiExchangeModule, only: MpiWorld !PAR
    ! -- dummy
    character(len=*),intent(in) :: filename
    integer(I4B), intent(in) :: id, m1i, m2i !PAR
    character(len=*), intent(in) :: mname1i, mname2i !PAR
    logical, intent(in) :: m2_bympi !PAR
    ! -- local
    type(GwfExchangeType), pointer :: exchange
    class(BaseModelType), pointer :: mb
    class(BaseExchangeType), pointer :: baseexchange
    character(len=20) :: cint
    integer(I4B) :: m1, m2 !PAR
    character(len=LINELENGTH) :: hmodel, m1name, m2name !PAR
! ------------------------------------------------------------------------------
    !
    ! -- Return in case this exchange does not have connected models
    if (m1i < 0 .and. m2i < 0) then !PAR
      return !PAR
    endif !PAR
    !
    ! -- Create a new exchange and add it to the baseexchangelist container
    allocate(exchange)
    baseexchange => exchange
    call AddBaseExchangeToList(baseexchangelist, baseexchange)
    !
    ! -- Assign id and name
    exchange%id = id
    write(cint, '(i0)') id
    exchange%name = 'GWF-GWF_' // trim(adjustl(cint))
    !
    ! -- allocate scalars and set defaults
    call exchange%allocate_scalars()
    exchange%filename = filename
    exchange%typename = 'GWF-GWF'
    exchange%implicit = .true.
    exchange%m2_bympi = m2_bympi !PAR
    !
    if (m1i > 0 .or. .not.m2_bympi) then !PAR
      m1 = m1i !PAR
      m2 = m2i !PAR
      m1name = mname1i !PAR
      m2name = mname2i !PAR
      exchange%m1m2_swap = .false. !PAR
    else !PAR
      m1 = m2i !PAR
      m2 = m1i !PAR
      m1name = mname2i !PAR
      m2name = mname1i !PAR
      exchange%m1m2_swap = .true. !PAR
    endif !PAR
    !
    ! -- Create the halo model
    !write(*,*) 'Myrank:',MpiWorld%myrank,' creating ','GWFHALO_'//trim(adjustl(cint)),' for "',trim(m1name),'" and "',trim(m2name),'"'
    hmodel = 'GWFHALO_'//trim(adjustl(cint)) !HALO2
    call gwfhalo_cr(exchange%gwfhalo, id, hmodel, m1name, m2name) !HALO2
    !
    ! -- Add to halo model list
    mb => exchange%gwfhalo !HALO2
    call AddBaseModelToList(halomodellist, mb) !HALO2
    !
    ! -- Add halo model to MPI world communicator
    if (m2_bympi) then !PAR
      call MpiWorld%mpi_addhmodel(hmodel, m1name, m2name) !PAR
    endif !PAR
    !
    ! -- set exchange%m1
    mb => GetBaseModelFromList(basemodellist, m1)
    select type (mb)
    class is (NumericalModelType)
      exchange%m1=>mb
    end select
    !
    ! -- set exchange%m2
    if (.not.m2_bympi) then !PAR
      mb => GetBaseModelFromList(basemodellist, m2)
      select type (mb)
      class is (NumericalModelType)
        exchange%m2=>mb
      end select
    endif !PAR
    !
    ! -- set gwfmodel1
    allocate(exchange%m1id) !HALO2
    exchange%m1id = m1 !HALO2
    mb => GetBaseModelFromList(basemodellist, m1)
    select type (mb)
    type is (GwfModelType)
      exchange%gwfmodel1 => mb
    end select
    !
    ! -- set gwfmodel2
    allocate(exchange%m2id) !HALO2
    exchange%m2id = m2 !HALO2
    if (.not.m2_bympi) then
      mb => GetBaseModelFromList(basemodellist, m2)
      select type (mb)
      type is (GwfModelType)
        exchange%gwfmodel2 => mb
      end select
    endif
    !
    ! -- Create the obs package
    call obs_cr(exchange%obs, exchange%inobs)
    !
    ! -- return
    return
  end subroutine gwfexchange_create

  subroutine gwf_gwf_df(this, iopt) !HALO2
! ******************************************************************************
! gwf_gwf_df -- Define GWF to GWF exchange object.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use SimVariablesModule, only: iout
    use InputOutputModule, only: getunit, openfile
    use GhostNodeModule, only: gnc_cr
    use BaseModelModule, only: BaseModelType !HALO2
    ! -- dummy
    class(GwfExchangeType) :: this
    integer(I4B), intent(in) :: iopt !HALO2
    ! -- local
    integer(I4B) :: inunit
! ------------------------------------------------------------------------------
    !
    if (iopt == 2) then !HALO2
      ! -- Define the halo model part 3
      call this%gwfhalo%gwfhalo_df3(this%nodem2, this%ihc, this%cl1,                & !HALO2 
        this%cl2, this%hwva) !HALO2
      return !HALO2
    end if !HALO2
    
    ! -- Define the halo model part 1
    call this%gwfhalo%gwfhalo_df1(this%m1id, this%m2id) !HALO2
    !
    ! -- open the file
    inunit = getunit()
    write(iout,'(/a,a)') ' Creating exchange: ', this%name
    call openfile(inunit, iout, this%filename, 'GWF-GWF')
    !
    call this%parser%Initialize(inunit, iout)
    !
    ! -- Ensure models are in same solution
!TODO    if(this%gwfmodel1%idsoln /= this%gwfmodel2%idsoln) then
!TODO      call store_error('ERROR.  TWO MODELS ARE CONNECTED ' //                  &
!TODO        'IN A GWF EXCHANGE BUT THEY ARE IN DIFFERENT SOLUTIONS. ' //           &
!TODO        'GWF MODELS MUST BE IN SAME SOLUTION: ' //                             &
!TODO        trim(this%gwfmodel1%name) // ' ' // trim(this%gwfmodel2%name) )
!TODO      call this%parser%StoreErrorUnit()
!TODO      call ustop()
!TODO    endif
    !
    ! -- read options
    call this%read_options(iout)
    !
    ! -- read dimensions
    call this%read_dimensions(iout)
    !
    ! -- allocate arrays
    call this%allocate_arrays()
    !
    ! -- read exchange data
    call this%read_data(iout)
    !
    ! -- call each model and increase the edge count
    call this%gwfmodel1%npf%increase_edge_count(this%nexg)
    if(associated(this%gwfmodel2)) then !HALO2
      call this%gwfmodel2%npf%increase_edge_count(this%nexg)
    endif !HALO2
    !
    ! -- Create and read ghost node information
    if(this%ingnc > 0) then
      call gnc_cr(this%gnc, this%name, this%ingnc, iout)
      call this%read_gnc(iout)
    endif
    !
    ! -- Read mover information
    if(this%inmvr > 0) then
      call this%read_mvr(iout)
    endif
    !
    ! -- close the file
    close(inunit)
    !
    ! -- Store obs
    call this%gwf_gwf_df_obs()
    call this%obs%obs_df(iout, this%name, 'GWF-GWF', this%gwfhalo%dis) !HALO2
    !
    ! -- Define the halo model part 2
    call this%gwfhalo%gwfhalo_df2(this%nodem1, this%inewton, this%m2_bympi) !HALO2
    !
    ! -- return
    return
  end subroutine gwf_gwf_df

  subroutine gwf_gwf_ac(this, sparse)
! ******************************************************************************
! gwf_gwf_ac -- override parent exg_ac so that gnc can add
!   connections here.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use SparseModule, only:sparsematrix
    ! -- dummy
    class(GwfExchangeType) :: this
    type(sparsematrix), intent(inout) :: sparse
    ! -- local
! ------------------------------------------------------------------------------
    !
    ! -- call parent model to add exchange connections
    call this%NumericalExchangeType%exg_ac(sparse)
    !
    ! -- add gnc connections
    if(this%ingnc > 0) then
      call this%gnc%gnc_ac(sparse)
    endif
    !
    ! -- Return
    return
  end subroutine gwf_gwf_ac

  subroutine gwf_gwf_mc(this, iasln, jasln)
! ******************************************************************************
! gwf_gwf_mc -- Map the connections in the global matrix
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use SparseModule, only:sparsematrix
    ! -- dummy
    class(GwfExchangeType) :: this
    integer(I4B), dimension(:), intent(in) :: iasln
    integer(I4B), dimension(:), intent(in) :: jasln
    ! -- local
! ------------------------------------------------------------------------------
    !
    ! -- call parent model to map exchange connections
    call this%NumericalExchangeType%exg_mc(iasln, jasln)
    !
    ! -- map gnc connections
    if(this%ingnc > 0) then
      call this%gnc%gnc_mc(iasln, jasln)
    endif
    !
    ! -- Return
    return
  end subroutine gwf_gwf_mc

  subroutine gwf_gwf_ar(this)
! ******************************************************************************
! gwf_gwf_ar -- Calculate the saturated conductance.  Must be called after
!               npf_ar for both GWF models.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ConstantsModule, only: LINELENGTH, DZERO, DHALF, DONE, DPIO180
    use SimModule, only: store_error, ustop
    use GwfNpfModule, only: condmean, vcond, hcond
    ! -- dummy
    class(GwfExchangeType) :: this
    ! -- local
    integer(I4B) :: iexg
    integer(I4B) :: n, m, ihc
    real(DP) :: topn, topm
    real(DP) :: botn, botm
    real(DP) :: satn, satm
    real(DP) :: thickn, thickm
    real(DP) :: angle, hyn, hym
    real(DP) :: csat
    real(DP) :: fawidth
    real(DP), dimension(3) :: vg
    character(len=LINELENGTH) :: errmsg
! ------------------------------------------------------------------------------
    !
    ! -- halo ar
    call this%gwfhalo%gwfhalo_ar() !HALO2
    !
    ! -- If mover is active, then call ar routine
    if(this%inmvr > 0) call this%mvr%mvr_ar()
    !
    ! -- Check to see if horizontal anisotropy is in either model1 or model2.
    !    If so, then ANGLDEGX must be provided as an auxiliary variable for this
    !    GWF-GWF exchange (this%ianglex > 0).
    if(this%gwfhalo%npf%ik22 /= 0) then !HALO2
      if(this%ianglex == 0) then
        write(errmsg, '(a)') 'Error.  GWF-GWF requires that ANGLDEGX be ' //   &
                             'specified as an auxiliary variable because ' //  &
                             'K22 was specified in one or both ' // &
                             'groundwater models.'
        call store_error(errmsg)
        call ustop()
      endif
    endif
    !
    ! -- Check to see if specific discharge is needed for model1 or model2.
    !    If so, then ANGLDEGX must be provided as an auxiliary variable for this
    !    GWF-GWF exchange (this%ianglex > 0).
    if(this%gwfhalo%npf%icalcspdis /= 0) then !HALO2
      if(this%ianglex == 0) then
        write(errmsg, '(a)') 'Error.  GWF-GWF requires that ANGLDEGX be ' //   &
                             'specified as an auxiliary variable because ' //  &
                             'specific discharge is being calculated in' // &
                             ' one or both groundwater models.'
        call store_error(errmsg)
        call ustop()
      endif
      if(this%icdist == 0) then
        write(errmsg, '(a)') 'Error.  GWF-GWF requires that CDIST be ' //   &
                             'specified as an auxiliary variable because ' //  &
                             'specific discharge is being calculated in' // &
                             ' one or both groundwater models.'
        call store_error(errmsg)
        call ustop()
      endif
    endif
    !
    ! -- Go through each connection and calculate the saturated conductance
    do iexg = 1, this%nexg
      !
      ihc = this%ihc(iexg)
      n = this%gwfhalo%imapnodem1tohalo(iexg) !HALO2
      m = this%gwfhalo%imapnodem2tohalo(iexg) + this%gwfhalo%offset !HALO2
      topn = this%gwfhalo%dis%top(n) !HALO2
      topm = this%gwfhalo%dis%top(m) !HALO2
      botn = this%gwfhalo%dis%bot(n) !HALO2
      botm = this%gwfhalo%dis%bot(m) !HALO2
      satn = this%gwfhalo%npf%sat(n) !HALO2
      satm = this%gwfhalo%npf%sat(m) !HALO2
      thickn = (topn - botn) * satn
      thickm = (topm - botm) * satm
      !
      ! -- Calculate conductance depending on connection orientation
      if(ihc == 0) then
        !
        ! -- Vertical conductance for fully saturated conditions
        vg(1) = DZERO
        vg(2) = DZERO
        vg(3) = DONE
        hyn = this%gwfhalo%npf%hy_eff(n, 0, ihc, vg=vg) !HALO2
        hym = this%gwfhalo%npf%hy_eff(m, 0, ihc, vg=vg) !HALO2
        csat = vcond(1, 1, 1, 1, 0, 1, 1, DONE,                                &
                      botn, botm,                                              &
                      hyn, hym,                                                &
                      satn, satm,                                              &
                      topn, topm,                                              &
                      botn, botm,                                              &
                      this%hwva(iexg))
      else
        !
        ! -- Calculate horizontal conductance
        hyn = this%gwfhalo%npf%k11(n) !HALO2
        hym = this%gwfhalo%npf%k11(m) !HALO2
        !
        ! -- Check for anisotropy in models, and recalculate hyn and hym
        if(this%ianglex > 0) then
          angle = this%auxvar(this%ianglex, iexg) * DPIO180
          vg(1) = abs(cos(angle))
          vg(2) = abs(sin(angle))
          vg(3) = DZERO
          !
          ! -- anisotropy in model 1 !TODO
          if(this%gwfhalo%npf%ik22 /= 0) then !HALO2
            hyn = this%gwfhalo%npf%hy_eff(n, 0, ihc, vg=vg) !HALO2
            hym = this%gwfhalo%npf%hy_eff(m, 0, ihc, vg=vg) !HALO2
          endif
        endif
        !
        fawidth = this%hwva(iexg)
        csat = hcond(1, 1, 1, 1, this%inewton, 0, ihc,                        &
                      this%icellavg, 0, 0, DONE,                              &
                      topn, topm, satn, satm, hyn, hym,                       &
                      topn, topm,                                             &
                      botn, botm,                                             &
                      this%cl1(iexg), this%cl2(iexg),                         &
                      fawidth, this%satomega)
      endif
      !
      ! -- store csat in condsat
      this%condsat(iexg) = csat
    enddo
    !
    ! -- Observation AR
    call this%obs%obs_ar()
    !
    ! -- Return
    return
  end subroutine gwf_gwf_ar

  subroutine gwf_gwf_rp(this)
! ******************************************************************************
! gwf_gwf_rp -- Read and prepare
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use TdisModule, only: readnewdata
    ! -- dummy
    class(GwfExchangeType) :: this
! ------------------------------------------------------------------------------
    !
    ! -- Check with TDIS on whether or not it is time to RP
    if (.not. readnewdata) return
    !
    ! -- Read and prepare for mover
    if(this%inmvr > 0) call this%mvr%mvr_rp()
    !
    ! -- Read and prepare for observations
    call this%gwf_gwf_rp_obs()
    !
    ! -- Return
    return
  end subroutine gwf_gwf_rp

  subroutine gwf_gwf_ad(this, isolnid, kpicard, isubtime)
! ******************************************************************************
! gwf_gwf_ad -- Initialize package x values to zero for explicit exchanges
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(GwfExchangeType) :: this
    integer(I4B), intent(in) :: isolnid
    integer(I4B), intent(in) :: kpicard
    integer(I4B), intent(in) :: isubtime
    ! -- local
! ------------------------------------------------------------------------------
    !
    ! -- Advance mover
    if(this%inmvr > 0) call this%mvr%mvr_ad()
    !
    ! -- Push simulated values to preceding time/subtime step
    call this%obs%obs_ad()
    !
    ! -- advance halo model
    call this%gwfhalo%model_ad(kpicard, isubtime) !HALO2
    !
    ! -- Return
    return
  end subroutine gwf_gwf_ad

  subroutine gwf_gwf_cf(this, kiter)
! ******************************************************************************
! gwf_gwf_cf -- Calculate the conductance term.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(GwfExchangeType) :: this
    integer(I4B), intent(in) :: kiter
    ! -- local
! ------------------------------------------------------------------------------
    !
    ! -- Set x and ibound for halo model
    call this%gwfhalo%gwfhalo_cf1(kiter) !HALO2
    !
    ! -- Rewet cells across models using the wetdry parameters in each model's
    !    npf package, and the head in the connected model.
    call this%rewet(kiter)
    !
    ! -- handle cf routines for halo model
    call this%gwfhalo%gwfhalo_cf2(kiter) !HALO2
    !
    ! -- Return
    return
  end subroutine gwf_gwf_cf

  subroutine gwf_gwf_fc(this, kiter, iasln, amatsln, inwtflag)
! ******************************************************************************
! gwf_gwf_fc -- Fill the matrix
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ConstantsModule, only: DHALF
    use GwfNpfModule, only: hcond, vcond
    ! -- dummy
    class(GwfExchangeType) :: this
    integer(I4B), intent(in) :: kiter
    integer(I4B), dimension(:), intent(in) :: iasln
    real(DP), dimension(:), intent(inout) :: amatsln
    integer(I4B), optional, intent(in) :: inwtflag
    ! -- local
    integer(I4B) :: inwt, iexg
    integer(I4B) :: njasln
    integer(I4B) :: n, m, nodensln, nodemsln, idiagsln !HALO2
    real(DP), dimension(3, 2) :: terms !HALO2
! ------------------------------------------------------------------------------
    !
    !
    ! -- Call fill method of parent to put this%cond into amatsln
    !call this%NumericalExchangeType%exg_fc(kiter, iasln, amatsln)
    !
    ! -- Use halo model to calculate amat terms
    !call this%gwfhalo%gwfhalo_fc(kiter, this%cond) !HALO2
    !
    if (this%m2_bympi) then !PAR
      do iexg = 1, this%nexg
        !
        ! -- initialize terms to zero
        terms(:, :) = DZERO !HALO2
        !
        n = this%nodem1(iexg) !HALO2
        nodensln = this%nodem1(iexg) + this%m1%moffset !HALO2
        call this%gwfhalo%gwfhalo_fc_calc(iexg, terms) !HALO2
        !
        ! -- row n
        idiagsln = iasln(nodensln) !HALO2
        amatsln(idiagsln) = amatsln(idiagsln) + terms(1, 1) !- this%cond(i) !HALO2
        this%cond(iexg) = terms(2, 1) !HALO2
        this%gwfmodel1%rhs(n) = this%gwfmodel1%rhs(n) + terms(3, 1) !HALO2
      end do 
    else !PAR
      do iexg = 1, this%nexg !HALO2
        !
        ! -- initialize terms to zero
        terms(:, :) = DZERO !HALO2
        !
        n = this%nodem1(iexg) !HALO2
        m = this%nodem2(iexg) !HALO2
        nodensln = this%nodem1(iexg) + this%m1%moffset !HALO2
        nodemsln = this%nodem2(iexg) + this%m2%moffset !HALO2
        call this%gwfhalo%gwfhalo_fc_calc(iexg, terms) !HALO2
        this%cond(iexg) = terms(2, 1) !HALO2
        !
        ! -- row n
        idiagsln = iasln(nodensln) !HALO2
        amatsln(idiagsln) = amatsln(idiagsln) + terms(1, 1) !- this%cond(i) !HALO2
        amatsln(this%idxglo(iexg)) = amatsln(this%idxglo(iexg)) + terms(2, 1) !this%cond(i) !HALO2
        this%gwfmodel1%rhs(n) = this%gwfmodel1%rhs(n) + terms(3, 1)
        !
        ! -- row m
        idiagsln = iasln(nodemsln) !HALO2
        amatsln(idiagsln) = amatsln(idiagsln) + terms(1, 2) !- this%cond(i) !HALO2
        amatsln(this%idxsymglo(iexg)) = amatsln(this%idxsymglo(iexg)) + terms(2, 2) !this%cond(i)!HALO2 
        this%gwfmodel2%rhs(m) = this%gwfmodel2%rhs(m) + terms(3, 2) !HALO2
      enddo !HALO2
    endif !HALO2
    !
    ! -- if gnc is active, then copy cond into gnc cond (might consider a
    !    pointer here in the future)
    if(this%ingnc > 0) then !HALO2
      do iexg = 1, this%nexg !HALO2
        this%gnc%cond(iexg) = this%cond(iexg) !HALO2
      enddo !HALO2
    endif !HALO2
    !
    ! -- Fill the gnc terms in the solution matrix
    if(this%ingnc > 0) then
      call this%gnc%gnc_fc(kiter, amatsln)
    endif
    !
    ! -- Call mvr fc routine
    if(this%inmvr > 0) call this%mvr%mvr_fc()
    !
    ! -- Set inwt to exchange newton, but shut off if requested by caller
    this%newtonterm = DZERO !PAR
    inwt = this%inewton
    if(present(inwtflag)) then
      if (inwtflag == 0) inwt = 0
    endif
    if (inwt /= 0) then
      call this%exg_fn(kiter, iasln, amatsln)
    endif
    !
    ! -- Ghost node Newton-Raphson
    if (this%ingnc > 0) then
      if (inwt /= 0) then
        njasln = size(amatsln)
        call this%gnc%gnc_fn(kiter, njasln, amatsln, this%condsat,             &
          ihc_opt=this%ihc, ivarcv_opt=this%ivarcv,                            &
          ictm1_opt=this%gwfhalo%npf%icelltype,                                & !HALO2
          ictm2_opt=this%gwfhalo%npf%icelltype) !HALO2
      endif
    endif
    !
    ! -- Return
    return
  end subroutine gwf_gwf_fc

  subroutine gwf_gwf_fn(this, kiter, iasln, amatsln)
! ******************************************************************************
! gwf_gwf_fn -- Fill amatsln with Newton terms
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use SmoothingModule, only: sQuadraticSaturationDerivative
    ! -- dummy
    class(GwfExchangeType) :: this
    integer(I4B), intent(in) :: kiter
    integer(I4B), dimension(:), intent(in) :: iasln
    real(DP), dimension(:), intent(inout) :: amatsln
    ! -- local
    integer(I4B) :: iexg, n, m, nodensln, nodemsln, idiagsln !HALO2 
    real(DP), dimension(3, 2) :: terms !HALO2
! ------------------------------------------------------------------------------
    !
    if (this%m2_bympi) then !PAR
      do iexg = 1, this%nexg !HALO2
        ! -- initialize terms to zero
        terms(:, :) = DZERO !HALO2
        !
        n = this%gwfhalo%nodem1(iexg) !HALO2
        call this%gwfhalo%gwfhalo_fn_calc(iexg, terms) !HALO2
        !
        nodensln = this%nodem1(iexg) + this%m1%moffset !HALO2
        idiagsln = iasln(nodensln) !HALO2
        amatsln(idiagsln) = amatsln(idiagsln) + terms(1, 1) !- this%cond(i) !HALO2
        this%newtonterm(iexg) = terms(2, 1) !HALO2
        this%gwfmodel1%rhs(n) = this%gwfmodel1%rhs(n) + terms(3, 1) !HALO2
      enddo
    else
      do iexg = 1, this%nexg !HALO2
        !
        ! -- initialize terms to zero
        terms(:, :) = DZERO !HALO2
        !
        n = this%gwfhalo%nodem1(iexg) !HALO2
        m = this%gwfhalo%nodem2(iexg) !HALO2
        call this%gwfhalo%gwfhalo_fn_calc(iexg, terms) !HALO2
        !
        nodensln = this%nodem1(iexg) + this%m1%moffset
        nodemsln = this%nodem2(iexg) + this%m2%moffset
        !
        ! -- row n
        idiagsln = iasln(nodensln) !HALO2
        amatsln(idiagsln) = amatsln(idiagsln) + terms(1, 1) !- this%cond(i) !HALO2
        amatsln(this%idxglo(iexg)) = amatsln(this%idxglo(iexg)) + terms(2, 1) !this%cond(i) !HALO2
        this%gwfmodel1%rhs(n) = this%gwfmodel1%rhs(n) + terms(3, 1) !HALO2
        !
        ! -- row m
        idiagsln = iasln(nodemsln) !HALO2
        amatsln(idiagsln) = amatsln(idiagsln) + terms(1, 2) !- this%cond(i) !HALO2
        amatsln(this%idxsymglo(iexg)) = amatsln(this%idxsymglo(iexg)) + terms(2, 2) !this%cond(i) !HALO2 
        this%gwfmodel2%rhs(m) = this%gwfmodel2%rhs(m) + terms(3, 2) !HALO2
      enddo !HALO2
    endif
    !
    ! -- Return
    return
  end subroutine gwf_gwf_fn

  subroutine gwf_gwf_cq(this, icnvg, isuppress_output, isolnid)
! ******************************************************************************
! gwf_gwf_cq -- Calculate flow between two cells
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ConstantsModule, only: DZERO, DPIO180
    use GwfNpfModule, only: thksatnm
    ! -- dummy
    class(GwfExchangeType) :: this
    integer(I4B), intent(inout) :: icnvg
    integer(I4B), intent(in) :: isuppress_output
    integer(I4B), intent(in) :: isolnid
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: n1
    integer(I4B) :: n2
    integer(I4B) :: ihc
    integer(I4B) :: ibdn1
    integer(I4B) :: ibdn2
    integer(I4B) :: ictn1
    integer(I4B) :: ictn2
    integer(I4B) :: iusg
    real(DP) :: topn1
    real(DP) :: topn2
    real(DP) :: botn1
    real(DP) :: botn2
    real(DP) :: satn1
    real(DP) :: satn2
    real(DP) :: hn1
    real(DP) :: hn2
    real(DP) :: rrate
    real(DP) :: thksat
    real(DP) :: angle
    real(DP) :: nx
    real(DP) :: ny
    real(DP) :: distance
    real(DP) :: dltot
    real(DP) :: hwva
    real(DP) :: area
! ------------------------------------------------------------------------------
    !
    ! -- Return if there neither model needs to calculate specific discharge
    if (this%gwfhalo%npf%icalcspdis == 0) return !HALO2
    !
    ! -- initialize
    iusg = 0
    !
    ! -- Loop through all exchanges
    do i = 1, this%nexg
      rrate = DZERO
      n1 = this%gwfhalo%imapnodem1tohalo(i) !HALO2
      n2 = this%gwfhalo%imapnodem2tohalo(i) + this%gwfhalo%offset !HALO2
      ihc = this%ihc(i)
      hwva = this%hwva(i)
      ibdn1 = this%gwfhalo%ibound(n1) !HALO2
      ibdn2 = this%gwfhalo%ibound(n2) !HALO2
      ictn1 = this%gwfhalo%npf%icelltype(n1) !HALO2
      ictn2 = this%gwfhalo%npf%icelltype(n2) !HALO2
      topn1 = this%gwfhalo%dis%top(n1) !HALO2
      topn2 = this%gwfhalo%dis%top(n2) !HALO2
      botn1 = this%gwfhalo%dis%bot(n1) !HALO2
      botn2 = this%gwfhalo%dis%bot(n2) !HALO2
      satn1 = this%gwfhalo%npf%sat(n1) !HALO2
      satn2 = this%gwfhalo%npf%sat(n2) !HALO2
      hn1 = this%gwfhalo%x(n1) !HALO2
      hn2 = this%gwfhalo%x(n2) !HALO2
      !
      ! -- If both cells are active then calculate flow rate, and add ghost
      !    node contribution
      if(ibdn1 /= 0 .and. ibdn2 /= 0) then
        rrate = this%qcalc(i, n1, n2)
        if(this%ingnc > 0) then
          rrate = rrate + this%gnc%deltaqgnc(i)
        endif
      endif
      !
      ! -- Calculate face normal components
      if(ihc == 0) then
        nx = DZERO
        ny = DZERO
        area = hwva
        if (botn1 < botn2) then
          ! -- n1 is beneath n2, so rate is positive downward.  Flip rate
          !    upward so that points in positive z direction
          rrate = - rrate
        endif
      else
        if(this%ianglex > 0) then
          angle = this%auxvar(this%ianglex, i) * DPIO180
          nx = cos(angle)
          ny = sin(angle)
        else
          ! error?
          call ustop('error in gwf_gwf_cq')
        endif
        !
        ! -- Calculate the saturated thickness at interface between n1 and n2
        thksat = thksatnm(ibdn1, ibdn2, ictn1, ictn2, this%inewton, ihc,       & 
                          iusg, hn1, hn2, satn1, satn2,                        &
                          topn1, topn2, botn1, botn2, this%satomega)
        area = hwva * thksat
      endif
      !
      ! -- Submit this connection and flow information to the npf
      !    package of gwfmodel1
      if(this%icdist > 0) then
        dltot = this%auxvar(this%icdist, i)
      else
        call ustop('error in gwf_gwf_cq')
      endif
      distance = dltot * this%cl1(i) / (this%cl1(i) + this%cl2(i))
!TODO      if (this%gwfmodel1%npf%icalcspdis == 1) then
!TODO        call this%gwfmodel1%npf%set_edge_properties(n1, ihc, rrate, area,      &
!TODO                                                    nx, ny, distance)
!TODO      endif
      !
      ! -- Submit this connection and flow information to the npf
      !    package of gwfmodel2
      if(this%icdist > 0) then
        dltot = this%auxvar(this%icdist, i)
      else
        call ustop('error in gwf_gwf_cq')
      endif
!TODO      if (this%gwfmodel2%npf%icalcspdis == 1) then
!TODO        distance = dltot * this%cl2(i) / (this%cl1(i) + this%cl2(i))
!TODO        if (ihc /= 0) rrate = -rrate
!TODO        call this%gwfmodel2%npf%set_edge_properties(n2, ihc, rrate, area,     &
!TODO                                                    nx, ny, distance)
!TODO      endif
      !
    enddo
    !
    ! -- return
    return
  end subroutine gwf_gwf_cq
  
  subroutine gwf_gwf_bd(this, icnvg, isuppress_output, isolnid)
! ******************************************************************************
! gwf_gwf_bd -- Budget for implicit gwf to gwf exchange; the budget for the
!               explicit exchange connections is handled for each model by
!               the exchange boundary package.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ConstantsModule, only: DZERO, LENBUDTXT, LENPACKAGENAME
    !use TdisModule, only: kstp, kper
    ! -- dummy
    class(GwfExchangeType) :: this
    integer(I4B), intent(inout) :: icnvg
    integer(I4B), intent(in) :: isuppress_output
    integer(I4B), intent(in) :: isolnid
    ! -- local
    character(len=LENBOUNDNAME) :: bname
    character(len=LENPACKAGENAME+4) :: packname1
    character(len=LENPACKAGENAME+4) :: packname2
    character(len=LENBUDTXT), dimension(1) :: budtxt
    character(len=20) :: nodestr
    integer(I4B) :: ntabrows
    integer(I4B) :: nodeu
    real(DP), dimension(2, 1) :: budterm
    integer(I4B) :: i, n1, n2, n1u, n2u, n1h, n2h !HALO2
    integer(I4B) :: ibinun1, ibinun2
    integer(I4B) :: icbcfl, ibudfl
    real(DP) :: ratin, ratout, rrate, deltaqgnc
    ! -- formats
! ------------------------------------------------------------------------------
    !
    ! -- initialize local variables
    budtxt(1) = '    FLOW-JA-FACE'
    packname1 = 'EXG '//this%name
    packname1 = adjustr(packname1)
    packname2 = 'EXG '//this%name
    packname2 = adjustr(packname2)
    !
    ! -- update output tables
    if (this%iprflow /= 0) then
      !
      ! -- update titles
      if (this%gwfmodel1%oc%oc_save('BUDGET')) then
        call this%outputtab1%set_title(packname1)
      end if
      if (this%gwfmodel2%oc%oc_save('BUDGET')) then 
        call this%outputtab2%set_title(packname2)
      end if
      !
      ! -- update maxbound of tables
      ntabrows = 0
      do i = 1, this%nexg
        n1 = this%nodem1(i)
        n2 = this%nodem2(i)
        !
        ! -- If both cells are active then calculate flow rate
        if (this%gwfmodel1%ibound(n1) /= 0 .and.                                  &
            this%gwfmodel2%ibound(n2) /= 0) then
          ntabrows = ntabrows + 1
        end if
      end do
      if (ntabrows > 0) then
        call this%outputtab1%set_maxbound(ntabrows)
        call this%outputtab2%set_maxbound(ntabrows)
      end if
    end if
    !
    ! -- Print and write budget terms for model 1
    !
    ! -- Set binary unit numbers for saving flows
    if(this%ipakcb /= 0) then
      ibinun1 = this%gwfmodel1%oc%oc_save_unit('BUDGET')
    else
      ibinun1 = 0
    endif
    !
    ! -- If save budget flag is zero for this stress period, then
    !    shut off saving
    if(.not. this%gwfmodel1%oc%oc_save('BUDGET')) ibinun1 = 0
    if(isuppress_output /= 0) then
      ibinun1 = 0
    endif
    !
    ! -- If cell-by-cell flows will be saved as a list, write header.
    if(ibinun1 /= 0) then
      call this%gwfmodel1%dis%record_srcdst_list_header(budtxt(1),             &
                                       this%gwfhalo%m1name, this%name,         & !HALO2
                                       this%gwfhalo%m2name, this%name,         & !HALO2
                                       this%naux, this%auxname,                &
                                       ibinun1, this%nexg, this%gwfmodel1%iout)
    endif
    !
    ! Initialize accumulators
    ratin = DZERO
    ratout = DZERO
    !
    ! -- Loop through all exchanges
    do i = 1, this%nexg
      !
      ! -- Assign boundary name
      if (this%inamedbound>0) then
        bname = this%boundname(i)
      else
        bname = ''
      endif
      !
      ! -- Calculate the flow rate between n1 and n2
      rrate = DZERO
      n1 = this%nodem1(i)
      n2 = this%nodem2(i)
      n1h = this%gwfhalo%imapnodem1tohalo(i) !HALO2
      n2h = this%gwfhalo%imapnodem2tohalo(i) + this%gwfhalo%offset !HALO2
      !
      ! -- If both cells are active then calculate flow rate
      if(this%gwfhalo%ibound(n1h) /= 0 .and. & !HALO2
         this%gwfhalo%ibound(n2h) /= 0) then !HALO2
        rrate = this%qcalc(i, n1h, n2h) !HALO2
        !
        ! -- add ghost node contribution
        if(this%ingnc > 0) then
          deltaqgnc = this%gnc%deltaqgnc(i)
          rrate = rrate + deltaqgnc
        endif
        !
        ! -- Print the individual rates to model list files if requested
        if(this%iprflow /= 0) then
          if(this%gwfmodel1%oc%oc_save('BUDGET')) then
            !
            ! -- set nodestr and write outputtab table
            nodeu = this%gwfmodel1%dis%get_nodeuser(n1)
            call this%gwfmodel1%dis%nodeu_to_string(nodeu, nodestr)
            call this%outputtab1%print_list_entry(i, trim(adjustl(nodestr)),     &
                                                  rrate, bname)
          end if
        endif
        if(rrate < DZERO) then
          ratout = ratout - rrate
        else
          ratin = ratin + rrate
        endif
      endif
      !
      ! -- If saving cell-by-cell flows in list, write flow
      n1u = this%gwfmodel1%dis%get_nodeuser(n1)
      if (this%m2_bympi) then !PAR
        n2u = this%nodeum2(i) !PAR
      else !PAR
        n2u = this%gwfmodel2%dis%get_nodeuser(n2)
      endif !PAR
      !
      if(ibinun1 /= 0)                                                         &
        call this%gwfmodel1%dis%record_mf6_list_entry(                         &
          ibinun1, n1u, n2u, rrate, this%naux, this%auxvar(:, i),              &
          .false., .false.)
      !
    enddo
    !
    ! -- Add the budget terms to model 1
    budterm(1, 1) = ratin
    budterm(2, 1) = ratout
    call this%m1%model_bdentry(budterm, budtxt, this%name)
    !
    if(.not.this%m2_bympi) then !PAR
    !
    ! -- Print and write budget terms for model 2
    !
    ! -- Set binary unit numbers for saving flows
    if(this%ipakcb /= 0) then
      ibinun2 = this%gwfmodel2%oc%oc_save_unit('BUDGET')
    else
      ibinun2 = 0
    endif
    !
    ! -- If save budget flag is zero for this stress period, then
    !    shut off saving
    if(.not. this%gwfmodel2%oc%oc_save('BUDGET')) ibinun2 = 0
    if(isuppress_output /= 0) then
      ibinun2 = 0
    endif
    !
    ! -- If cell-by-cell flows will be saved as a list, write header.
    if(ibinun2 /= 0) then
      call this%gwfmodel2%dis%record_srcdst_list_header(budtxt(1),             &
                                       this%gwfhalo%m2name, this%name,         & !HALO2
                                       this%gwfhalo%m1name, this%name,         & !HALO2
                                       this%naux, this%auxname,                &
                                       ibinun2, this%nexg, this%gwfmodel2%iout)
    endif
    !
    ! Initialize accumulators
    ratin = DZERO
    ratout = DZERO
    !
    ! -- Loop through all exchanges
    do i = 1, this%nexg
      !
      ! -- Assign boundary name
      if (this%inamedbound>0) then
        bname = this%boundname(i)
      else
        bname = ''
      endif
      !
      ! -- Calculate the flow rate between n1 and n2
      rrate = DZERO
      n1 = this%nodem1(i)
      n2 = this%nodem2(i)
      !
      ! -- If both cells are active then calculate flow rate
      if(this%gwfmodel1%ibound(n1) /= 0 .and. &
          this%gwfmodel2%ibound(n2) /= 0) then
        rrate = this%cond(i) * this%m2%x(n2) - this%cond(i) * this%m1%x(n1)
        !
        ! -- add ghost node contribution
        if(this%ingnc > 0) then
          deltaqgnc = this%gnc%deltaqgnc(i)
          rrate = rrate + deltaqgnc
        endif
        !
        ! -- Print the individual rates to model list files if requested
        if(this%iprflow /= 0) then
          if(this%gwfmodel2%oc%oc_save('BUDGET')) then
            !
            ! -- set nodestr and write outputtab table
            nodeu = this%gwfmodel2%dis%get_nodeuser(n2)
            call this%gwfmodel2%dis%nodeu_to_string(nodeu, nodestr)
            call this%outputtab2%print_list_entry(i, trim(adjustl(nodestr)),     &
                                                  -rrate, bname)
          end if
        endif
        if(rrate < DZERO) then
          ratout = ratout - rrate
        else
          ratin = ratin + rrate
        endif
      endif
      !
      ! -- If saving cell-by-cell flows in list, write flow
      n1u = this%gwfmodel1%dis%get_nodeuser(n1)
      n2u = this%gwfmodel2%dis%get_nodeuser(n2)
      if(ibinun2 /= 0)                                                         &
        call this%gwfmodel2%dis%record_mf6_list_entry(                         &
          ibinun2, n2u, n1u, -rrate, this%naux, this%auxvar(:, i),             &
          .false., .false.)
      !
    enddo
    !
    ! -- Add the budget terms to model 2
    budterm(1, 1) = ratout
    budterm(2, 1) = ratin
    call this%m2%model_bdentry(budterm, budtxt, this%name)
    !
    endif !PAR
    !
    ! -- Set icbcfl, ibudfl to zero so that flows will be printed and
    !    saved, if the options were set in the MVR package
    icbcfl = 1
    ibudfl = 1
    !
    ! -- Call mvr bd routine
    if(this%inmvr > 0) call this%mvr%mvr_bd(icbcfl, ibudfl, isuppress_output)
    !
    ! -- Calculate and write simulated values for observations
    if(this%inobs /= 0) then
      call this%gwf_gwf_save_simvals()
    endif
    !
    ! -- return
    return
  end subroutine gwf_gwf_bd

  subroutine gwf_gwf_ot(this)
! ******************************************************************************
! gwf_gwf_ot
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use SimVariablesModule, only: iout
    use ConstantsModule, only: DZERO, LINELENGTH
    ! -- dummy
    class(GwfExchangeType) :: this
    ! -- local
    integer(I4B) :: iexg, n1, n2
    real(DP) :: flow, deltaqgnc
    character(len=LINELENGTH) :: node1str, node2str
    ! -- format
    character(len=*), parameter :: fmtheader =                                 &
     "(/1x, 'SUMMARY OF EXCHANGE RATES FOR EXCHANGE ', a, ' WITH ID ', i0, /,  &
       &2a16, 5a16, /, 112('-'))"
    character(len=*), parameter :: fmtheader2 =                                &
     "(/1x, 'SUMMARY OF EXCHANGE RATES FOR EXCHANGE ', a, ' WITH ID ', i0, /,  &
       &2a16, 4a16, /, 96('-'))"
    character(len=*), parameter :: fmtdata =                                   &
     "(2a16, 5(1pg16.6))"
! ------------------------------------------------------------------------------
    !
    ! -- Initialize
    deltaqgnc = DZERO
    !
    ! -- Write a table of exchanges
    if(this%iprflow /= 0) then
      if(this%ingnc > 0) then
        write(iout, fmtheader) trim(adjustl(this%name)), this%id, 'NODEM1',    &
                             'NODEM2', 'COND', 'X_M1', 'X_M2', 'DELTAQGNC',    &
                             'FLOW'
      else
        write(iout, fmtheader2) trim(adjustl(this%name)), this%id, 'NODEM1',   &
                             'NODEM2', 'COND', 'X_M1', 'X_M2', 'FLOW'
      endif
      do iexg = 1, this%nexg
        n1 = this%gwfhalo%imapnodem1tohalo(iexg) !HALO2 TODO
        n2 = this%gwfhalo%imapnodem2tohalo(iexg) + this%gwfhalo%offset !HALO2 TODO
        flow = this%cond(iexg) * (this%gwfhalo%x(n2) - this%gwfhalo%x(n1)) !HALO2
        call this%gwfhalo%dis%noder_to_string(n1, node1str)!HALO2 TODO
        call this%gwfhalo%dis%noder_to_string(n2, node2str)!HALO2 TODO
        if(this%ingnc > 0) then
          deltaqgnc = this%gnc%deltaqgnc(iexg)
          write(iout, fmtdata) trim(adjustl(node1str)),                        &
                               trim(adjustl(node2str)),                        &
                               this%cond(iexg), this%gwfhalo%x(n1),            & !HALO2
                               this%gwfhalo%x(n2),                             & !HALO2
                               deltaqgnc, flow + deltaqgnc
        else
          write(iout, fmtdata) trim(adjustl(node1str)),                        &
                               trim(adjustl(node2str)),                        &
                               this%cond(iexg), this%gwfhalo%x(n1),            & !HALO2
                               this%gwfhalo%x(n2),                             & !HALO2
                               flow
        endif
      enddo
    endif
    !
    ! -- Mover budget output
    if(this%inmvr > 0) call this%mvr%mvr_ot()
    !
    ! -- OBS output
    call this%obs%obs_ot()
    !
    ! -- return
    return
  end subroutine gwf_gwf_ot

  subroutine read_options(this, iout)
! ******************************************************************************
! read_options -- Read Options
! Subroutine: (1) read options from input file
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ArrayHandlersModule, only: ifind
    use ConstantsModule, only: LINELENGTH, DEM6
    use InputOutputModule, only: getunit, openfile, urdaux
    use SimModule, only: store_error, store_error_unit, ustop
    ! -- dummy
    class(GwfExchangeType) :: this
    integer(I4B), intent(in) :: iout
    ! -- local
    character(len=LINELENGTH) :: line, errmsg, keyword, fname
    integer(I4B) :: istart,istop,lloc,ierr,ival
    integer(I4B) :: inobs
    logical :: isfound, endOfBlock
! ------------------------------------------------------------------------------
    !
    ! -- get options block
    call this%parser%GetBlock('OPTIONS', isfound, ierr,                        &
      supportOpenClose=.true., blockRequired=.false.)
    !
    ! -- parse options block if detected
    if (isfound) then
      write(iout,'(1x,a)')'PROCESSING GWF EXCHANGE OPTIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
          case('AUXILIARY')
            call this%parser%GetRemainingLine(line)
            lloc = 1
            call urdaux(this%naux, this%parser%iuactive, iout, lloc, istart,   &
                        istop, this%auxname, line, 'GWF_GWF_Exchange')
            !
            ! -- If ANGLDEGX is an auxiliary variable, then anisotropy can be
            !    used in either model.  Store ANGLDEGX position in this%ianglex
            ival = ifind(this%auxname, 'ANGLDEGX')
            if(ival > 0) this%ianglex = ival
            ival = ifind(this%auxname, 'CDIST')
            if(ival > 0) this%icdist = ival
          case ('PRINT_INPUT')
            this%iprpak = 1
            write(iout,'(4x,a)') &
              'THE LIST OF EXCHANGES WILL BE PRINTED.'
          case ('PRINT_FLOWS')
            this%iprflow = 1
            write(iout,'(4x,a)') &
              'EXCHANGE FLOWS WILL BE PRINTED TO LIST FILES.'
          case ('SAVE_FLOWS')
            this%ipakcb = -1
            write(iout,'(4x,a)') &
              'EXCHANGE FLOWS WILL BE SAVED TO BINARY BUDGET FILES.'
          case ('ALTERNATIVE_CELL_AVERAGING')
            call this%parser%GetStringCaps(keyword)
            select case(keyword)
            case('LOGARITHMIC')
              this%icellavg = 1
            case('AMT-LMK')
              this%icellavg = 2
            case default
              write(errmsg,'(4x,a,a)')'UNKNOWN CELL AVERAGING METHOD: ',       &
                                       trim(keyword)
              call store_error(errmsg)
              call this%parser%StoreErrorUnit()
              call ustop()
            end select
            write(iout,'(4x,a,a)')                                             &
              'CELL AVERAGING METHOD HAS BEEN SET TO: ', trim(keyword)
          case ('VARIABLECV')
            this%ivarcv = 1
            write(iout,'(4x,a)')                                               &
              'VERTICAL CONDUCTANCE VARIES WITH WATER TABLE.'
            call this%parser%GetStringCaps(keyword)
            if(keyword == 'DEWATERED') then
              this%idewatcv = 1
              write(iout,'(4x,a)')                                             &
                'VERTICAL CONDUCTANCE ACCOUNTS FOR DEWATERED PORTION OF   ' // &
                'AN UNDERLYING CELL.'
            endif
          case ('NEWTON')
            this%inewton = 1
            write(iout, '(4x,a)')                                              &
                             'NEWTON-RAPHSON method used for unconfined cells'
          case ('GNC6')
            call this%parser%GetStringCaps(keyword)
            if(keyword /= 'FILEIN') then
              call store_error('GNC6 KEYWORD MUST BE FOLLOWED BY ' //          &
                '"FILEIN" then by filename.')
              call this%parser%StoreErrorUnit()
              call ustop()
            endif
            call this%parser%GetString(fname)
            if(fname == '') then
              call store_error('NO GNC6 FILE SPECIFIED.')
              call this%parser%StoreErrorUnit()
              call ustop()
            endif
            this%ingnc = getunit()
            call openfile(this%ingnc, iout, fname, 'GNC')
            write(iout,'(4x,a)')                                               &
              'GHOST NODES WILL BE READ FROM ', trim(fname)
          case ('MVR6')
            call this%parser%GetStringCaps(keyword)
            if(keyword /= 'FILEIN') then
              call store_error('MVR6 KEYWORD MUST BE FOLLOWED BY ' //          &
                '"FILEIN" then by filename.')
              call this%parser%StoreErrorUnit()
              call ustop()
            endif
            call this%parser%GetString(fname)
            if(fname == '') then
              call store_error('NO MVR6 FILE SPECIFIED.')
              call this%parser%StoreErrorUnit()
              call ustop()
            endif
            this%inmvr = getunit()
            call openfile(this%inmvr, iout, fname, 'MVR')
            write(iout,'(4x,a)')                                               &
              'WATER MOVER INFORMATION WILL BE READ FROM ', trim(fname)
          case ('BOUNDNAMES')
            this%inamedbound = 1
            write(iout,'(4x,a)') 'EXCHANGE BOUNDARIES HAVE NAMES' // &
                                      ' IN LAST COLUMN.'
          case ('OBS6')
            call this%parser%GetStringCaps(keyword)
            if(keyword /= 'FILEIN') then
              call store_error('OBS8 KEYWORD MUST BE FOLLOWED BY ' //         &
                '"FILEIN" then by filename.')
              call this%parser%StoreErrorUnit()
              call ustop()
            endif
            this%obs%active = .true.
            call this%parser%GetString(this%obs%inputFilename)
            inobs = GetUnit()
            call openfile(inobs, iout, this%obs%inputFilename, 'OBS')
            this%obs%inUnitObs = inobs
          case default
            write(errmsg,'(4x,a,a)')'***ERROR. UNKNOWN GWF EXCHANGE OPTION: ', &
                                     trim(keyword)
            call store_error(errmsg)
            call this%parser%StoreErrorUnit()
            call ustop()
        end select
      end do
      write(iout,'(1x,a)')'END OF GWF EXCHANGE OPTIONS'
    end if
    !
    ! -- set omega value used for saturation calculations
    if (this%inewton > 0) then
      this%satomega = DEM6
    end if
    !
    ! -- return
    return
  end subroutine read_options

  subroutine read_data(this, iout)
! ******************************************************************************
! read_data -- Read EXGDATA block
! Subroutine: (1) read list of EXGs from input file
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ConstantsModule, only: LINELENGTH
    use SimModule, only: ustop, store_error, store_error_unit, count_errors
    use BaseModelModule, only: BaseModelType
    use BaseDisModule, only: DisBaseType !PAR
    use MpiExchangeGenModule, only: mpi_is_halo !PAR
    use MpiExchangeGwfModule, only: mpi_set_gwfhalo_world_dis !PAR
    use MpiExchangeModule, only: MpiWorld !PAR !@@@@
    ! -- dummy
    class(gwfExchangeType) :: this
    integer(I4B), intent(in) :: iout
    ! -- local
    character(len=LINELENGTH) :: errmsg, nodestr, node1str, node2str, cellid
    character(len=2) :: cnfloat
    integer(I4B) :: lloc, ierr, nerr, iaux
    integer(I4B) :: iexg, nodem1, nodem2, nodeum1, nodeum2
    logical :: isfound, endOfBlock
    class(NumericalModelType), pointer :: m2 !PAR
    ! -- format
    character(len=*), parameter :: fmtexglabel = "(5x, 3a10, 50(a16))"
    character(len=*), parameter :: fmtexgdata  =                               &
      "(5x, a, 1x, a ,I10, 50(1pg16.6))"
    character(len=40) :: fmtexgdata2
! ------------------------------------------------------------------------------
    !
    if (this%m2_bympi) then !PAR
      call mpi_set_gwfhalo_world_dis(this%gwfhalo%m2name, m2) !PAR
!      call m2tmp%dis_da() !PAR
!      deallocate(m2tmp) !PAR
    else !PAR
      m2 => this%m2
    endif !PAR
    !
    ! -- get ExchangeData block
    call this%parser%GetBlock('EXCHANGEDATA', isfound, ierr,                   &
                              supportOpenClose=.true.)
    !
    ! -- parse ExchangeData block if detected
    if (isfound) then
      write(iout,'(1x,a)')'PROCESSING EXCHANGEDATA'
      if(this%iprpak /= 0) then
        if (this%inamedbound==0) then
          write(iout, fmtexglabel) 'NODEM1', 'NODEM2', 'IHC',                  &
              'CL1', 'CL2', 'HWVA', (adjustr(this%auxname(iaux)),              &
              iaux = 1, this%naux)
        else
          write(iout, fmtexglabel) 'NODEM1', 'NODEM2', 'IHC', 'CL1', 'CL2',    &
              'HWVA', (adjustr(this%auxname(iaux)),iaux=1,this%naux),          &
              ' BOUNDNAME      '
          ! Define format suitable for writing input data,
          ! any auxiliary variables, and boundname.
          write(cnfloat,'(i0)') 3+this%naux
          fmtexgdata2 = '(5x, a, 1x, a, i10, ' // trim(cnfloat) //             &
            '(1pg16.6), 1x, a)'
        endif
      endif
      !
      do iexg = 1, this%nexg
        call this%parser%GetNextLine(endOfBlock)
        lloc = 1
        !
        if (.not.this%m1m2_swap) then !PAR
          ! -- Read and check node 1
          call this%parser%GetCellid(this%m1%dis%ndim, cellid, flag_string=.true.)
          nodem1 = this%m1%dis%noder_from_cellid(cellid, this%parser%iuactive,   &
                                                 iout, flag_string=.true.)
          this%nodem1(iexg) = nodem1
          !
          ! -- Read and check node 2
          call this%parser%GetCellid(m2%dis%ndim, cellid, flag_string=.true.)
          nodem2 = m2%dis%noder_from_cellid(cellid, this%parser%iuactive,   &
                                            iout, flag_string=.true.)
          this%nodem2(iexg) = nodem2
          if (this%m2_bympi) then !PAR
            this%nodeum2(iexg) = nodem2 !PAR
            this%nodem2(iexg) = iexg !PAR
          endif !PAR
        else
          ! -- Read and check node 2
          call this%parser%GetCellid(m2%dis%ndim, cellid, flag_string=.true.) !HALO2
          nodem2 = m2%dis%noder_from_cellid(cellid, this%parser%iuactive,   & !HALO2
                                            iout, flag_string=.true.) !HALO2
          this%nodem2(iexg) = nodem2
          if (this%m2_bympi) then !PAR
            this%nodeum2(iexg) = nodem2 !PAR
            this%nodem2(iexg) = iexg !PAR
          endif !PAR
          !
          ! -- Read and check node 1
          call this%parser%GetCellid(this%m1%dis%ndim, cellid, flag_string=.true.) !HALO2
          nodem1 = this%m1%dis%noder_from_cellid(cellid, this%parser%iuactive,   & !HALO2
                                                 iout, flag_string=.true.) !HALO2
          this%nodem1(iexg) = nodem1
        endif
        !
        ! -- Read rest of input line
        this%ihc(iexg) = this%parser%GetInteger()
        if (.not.this%m1m2_swap) then !PAR
          this%cl1(iexg) = this%parser%GetDouble()
          this%cl2(iexg) = this%parser%GetDouble()
        else
          this%cl2(iexg) = this%parser%GetDouble() !PAR
          this%cl1(iexg) = this%parser%GetDouble() !PAR
        endif !PAR
        this%hwva(iexg) = this%parser%GetDouble()
        do iaux = 1, this%naux
          this%auxvar(iaux, iexg) = this%parser%GetDouble()
        enddo
        if (this%inamedbound==1) then
          call this%parser%GetStringCaps(this%boundname(iexg))
        endif
        !
        ! -- Write the data to listing file if requested
        !HALO2 TODO
        !if(this%iprpak /= 0) then
        !  nodeum1 = this%m1%dis%get_nodeuser(nodem1)
        !  call this%m1%dis%nodeu_to_string(nodeum1, node1str)
        !  if (.not. mpi_is_halo(this%m2%name)) then !PAR
        !    nodeum2 = this%m2%dis%get_nodeuser(nodem2)
        !  else !PAR
        !    nodeum2 = nodem2 !PAR
        !  endif !PAR 
        !  call this%m2%dis%nodeu_to_string(nodeum2, node2str)
        !  if (this%inamedbound == 0) then
        !    write(iout, fmtexgdata) trim(node1str), trim(node2str),            &
        !                this%ihc(iexg), this%cl1(iexg), this%cl2(iexg),        &
        !                this%hwva(iexg),                                       &
        !                (this%auxvar(iaux, iexg), iaux=1,this%naux)
        !  else
        !    write(iout, fmtexgdata2) trim(node1str), trim(node2str),           &
        !                this%ihc(iexg), this%cl1(iexg), this%cl2(iexg),        &
        !                this%hwva(iexg),                                       &
        !                (this%auxvar(iaux, iexg), iaux=1,this%naux),           &
        !                trim(this%boundname(iexg))
        !  endif
        !endif
        !
        ! -- Check to see if nodem1 is outside of active domain
        if(nodem1 <= 0) then
          call this%gwfmodel1%dis%nodeu_to_string(nodeum1, nodestr)
          write(errmsg, *)                                                     &
                  trim(adjustl(this%gwfmodel1%name)) //                        &
                  ' Cell is outside active grid domain: ' //                   &
                  trim(adjustl(nodestr))
          call store_error(errmsg)
        endif
        !
        ! -- Check to see if nodem2 is outside of active domain
        if(nodem2 <= 0) then
          call this%gwfmodel2%dis%nodeu_to_string(nodeum2, nodestr)
          write(errmsg, *)                                                     &
                  trim(adjustl(this%gwfmodel2%name)) //                        &
                  ' Cell is outside active grid domain: ' //                   &
                  trim(adjustl(nodestr))
          call store_error(errmsg)
        endif
      enddo
      !
      ! -- Stop if errors
      nerr = count_errors()
      if(nerr > 0) then
        call store_error('Errors encountered in exchange input file.')
        call this%parser%StoreErrorUnit()
        call ustop()
      endif
      !
      write(iout,'(1x,a)')'END OF EXCHANGEDATA'
    else
      write(errmsg, '(1x,a)')'ERROR.  REQUIRED EXCHANGEDATA BLOCK NOT FOUND.'
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
      call ustop()
    end if
    !
    if (this%m2_bympi) then !PAR
      call m2%dis%dis_da() !PAR
      deallocate(m2) !PAR
    endif !PAR
    !
    ! -- return
    return
  end subroutine read_data

  subroutine read_gnc(this, iout)
! ******************************************************************************
! read_gnc -- Read ghost node information.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use SimModule, only: store_error, store_error_unit, count_errors, ustop
    use ConstantsModule, only: LINELENGTH
    ! -- dummy
    class(GwfExchangeType) :: this
    integer(I4B), intent(in) :: iout
    ! -- local
    integer(I4B) :: i, nm1, nm2, nmgnc1, nmgnc2
    character(len=LINELENGTH) :: errmsg
    character(len=*), parameter :: fmterr = &
      "('EXCHANGE NODES ', i0, ' AND ', i0,"  // &
      "' NOT CONSISTENT WITH GNC NODES ', i0, ' AND ', i0)"
! ------------------------------------------------------------------------------
    !
    ! -- If exchange has ghost nodes, then initialize ghost node object
    !    This will read the ghost node blocks from the gnc input file.
    call this%gnc%gnc_df(this%m1, m2=this%m2)
    !
    ! -- Verify gnc is implicit if exchange has Newton Terms
    if(.not. this%gnc%implicit .and. this%inewton /= 0) then
      call store_error('GNC IS EXPLICIT, BUT GWF EXCHANGE HAS ACTIVE NEWTON.')
      call store_error('ADD IMPLICIT OPTION TO GNC OR REMOVE NEWTON FROM ' // &
        'GWF EXCHANGE.')
      call store_error_unit(this%ingnc)
      call ustop()
    endif
    !
    ! -- Perform checks to ensure GNCs match with GWF-GWF nodes
    if(this%nexg /= this%gnc%nexg) then
      call store_error('NUMBER OF EXCHANGES DOES NOT MATCH NUMBER OF GNCs')
      call store_error_unit(this%ingnc)
      call ustop()
    endif
    !
    ! -- Go through each entry and confirm
    do i = 1, this%nexg
      if(this%nodem1(i) /= this%gnc%nodem1(i) .or.                             &
          this%nodem2(i) /= this%gnc%nodem2(i) ) then
        nm1 = this%gwfmodel1%dis%get_nodeuser(this%nodem1(i))
        nm2 = this%gwfmodel2%dis%get_nodeuser(this%nodem2(i))
        nmgnc1 = this%gwfmodel1%dis%get_nodeuser(this%gnc%nodem1(i))
        nmgnc2 = this%gwfmodel2%dis%get_nodeuser(this%gnc%nodem2(i))
        write(errmsg, fmterr) nm1, nm2, nmgnc1, nmgnc2
        call store_error(errmsg)
      endif
    enddo
    if(count_errors() > 0) then
      call store_error_unit(this%ingnc)
      call ustop()
    endif
    !
    ! -- close the file
    close(this%ingnc)
    !
    ! -- return
    return
  end subroutine read_gnc

  subroutine read_mvr(this, iout)
! ******************************************************************************
! read_mvr -- Read water mover information.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use GwfMvrModule, only: mvr_cr
    ! -- dummy
    class(GwfExchangeType) :: this
    integer(I4B), intent(in) :: iout
    ! -- local
! ------------------------------------------------------------------------------
    !
    ! -- Create and initialize the mover object
    call mvr_cr(this%mvr, this%name, this%inmvr, iout, iexgmvr=1)
    !
    ! -- Return
    return
  end subroutine read_mvr

  subroutine rewet(this, kiter)
! ******************************************************************************
! rewet -- Check for rewetting across models
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use TdisModule, only: kper, kstp
    ! -- dummy
    class(GwfExchangeType) :: this
    integer(I4B), intent(in) :: kiter
    ! -- local
    integer(I4B) :: iexg
    integer(I4B) :: n, m, nn, mm, nu, mu !HALO2
    integer(I4B) :: ibdn, ibdm
    integer(I4B) :: ihc
    real(DP) :: hn, hm
    integer(I4B) :: irewet
    character(len=30) :: nodestrn, nodestrm
    character(len=*),parameter :: fmtrwt =                                     &
      "(1x, 'CELL ',A,' REWET FROM GWF MODEL ',A,' CELL ',A,                   &
       &' FOR ITER. ',I0, ' STEP ',I0, ' PERIOD ', I0)"
! ------------------------------------------------------------------------------
    !
    ! -- Use model 1 to rewet model 2 and vice versa
    do iexg = 1, this%nexg
      n = this%gwfhalo%imapnodem1tohalo(iexg) !HALO2
      m = this%gwfhalo%imapnodem2tohalo(iexg) + this%gwfhalo%offset !HALO2
      hn = this%gwfhalo%x(n) !HALO2 
      hm = this%gwfhalo%x(m) !HALO2
      ibdn = this%gwfhalo%ibound(n) !HALO2
      ibdm = this%gwfhalo%ibound(m) !HALO2
      ihc = this%ihc(iexg)
      this%gwfhalo%npf%irewet = this%gwfhalo%m1irewet !HALO2
      call this%gwfhalo%npf%rewet_check(kiter, n, hm, ibdm, ihc,               & !HALO2
        this%gwfhalo%x, irewet) !HALO2
      if(irewet == 1) then !HALO2 TODO
        nn = this%gwfhalo%nodem1(iexg) !HALO2
        this%gwfhalo%gwf1%ibound(nn) = this%gwfhalo%ibound(n) !HALO2
        this%gwfhalo%gwf1%x(nn) = this%gwfhalo%x(n) !HALO2
        if(this%gwfhalo%ibound(n)==30000) then !HALO2
          this%gwfhalo%ibound(n)=1 !HALO2
        endif !HALO2
        call this%gwfhalo%dis%noder_to_string(n, nodestrn) !HALO2
        call this%gwfhalo%dis%noder_to_string(m, nodestrm) !HALO2
        !TODO write(this%gwfmodel1%iout, fmtrwt) trim(nodestrn),               & !HALO2
        !TODO  trim(this%gwfmodel2%name), trim(nodestrm), kiter, kstp, kper !HALO2
      endif !HALO2
      this%gwfhalo%npf%irewet = this%gwfhalo%m2irewet !HALO2
      call this%gwfhalo%npf%rewet_check(kiter, m, hn, ibdn, ihc,               & !HALO2
        this%gwfhalo%x, irewet) !HALO2
      if(irewet == 1) then !HALO2 TODO
        if (associated(this%gwfhalo%gwf2)) then !HALO2
          mm = this%gwfhalo%nodem2(iexg) !HALO2
          this%gwfhalo%gwf2%ibound(mm) = this%gwfhalo%ibound(m) !HALO2
          this%gwfhalo%gwf2%x(mm) = this%gwfhalo%x(m) !HALO2
          if(this%gwfhalo%ibound(m)==30000) then !HALO2
            this%gwfhalo%ibound(m)=1 !HALO2
          endif !HALO2
        endif
        call this%gwfhalo%dis%noder_to_string(n, nodestrm) !HALO2
        call this%gwfhalo%dis%noder_to_string(m, nodestrn) !HALO2
        !TODOwrite(this%gwfmodel2%iout, fmtrwt) trim(nodestrn),                & !HALO2
        !TODO  trim(this%gwfmodel1%name), trim(nodestrm), kiter, kstp, kper !HALO2
      endif
      !
    enddo
    !
    ! -- Return
    return
  end subroutine rewet

  subroutine allocate_scalars(this)
! ******************************************************************************
! allocate_scalars
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    use ConstantsModule, only: LENORIGIN, DZERO
    ! -- dummy
    class(GwfExchangeType) :: this
    ! -- local
    character(len=LENORIGIN) :: origin
! ------------------------------------------------------------------------------
    !
    ! -- create the origin name
    origin = trim(this%name)
    !
    ! -- Call parent type allocate_scalars
    call this%NumericalExchangeType%allocate_scalars()
    !
    call mem_allocate(this%icellavg, 'ICELLAVG', origin)
    call mem_allocate(this%ivarcv, 'IVARCV', origin)
    call mem_allocate(this%idewatcv, 'IDEWATCV', origin)
    call mem_allocate(this%inewton, 'INEWTON', origin)
    call mem_allocate(this%ianglex, 'IANGLEX', origin)
    call mem_allocate(this%icdist, 'ICDIST', origin)
    call mem_allocate(this%ingnc, 'INGNC', origin)
    call mem_allocate(this%inmvr, 'INMVR', origin)
    call mem_allocate(this%inobs, 'INOBS', origin)
    call mem_allocate(this%inamedbound, 'INAMEDBOUND', origin)
    call mem_allocate(this%satomega, 'SATOMEGA', origin)
    this%icellavg = 0
    this%ivarcv = 0
    this%idewatcv = 0
    this%inewton = 0
    this%ianglex = 0
    this%icdist = 0
    this%ingnc = 0
    this%inmvr = 0
    this%inobs = 0
    this%inamedbound = 0
    this%satomega = DZERO
    this%m2_bympi = .false. !PAR
    this%m1m2_swap = .false. !PAR
    !
    ! -- return
    return
  end subroutine allocate_scalars

  subroutine gwf_gwf_da(this)
! ******************************************************************************
! gwf_gwf_da
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(GwfExchangeType) :: this
    ! -- local
! ------------------------------------------------------------------------------
    !
    ! -- Call parent type allocate_scalars
    call this%NumericalExchangeType%exg_da()
    !
    ! -- objects
    if(this%ingnc > 0) then
      call this%gnc%gnc_da()
      deallocate(this%gnc)
    endif
    if (this%inmvr > 0) then
      call this%mvr%mvr_da()
      deallocate(this%mvr)
    endif
    call this%obs%obs_da()
    deallocate(this%obs)
    !
    ! -- arrays
    call mem_deallocate(this%ihc)
    call mem_deallocate(this%cl1)
    call mem_deallocate(this%cl2)
    call mem_deallocate(this%hwva)
    call mem_deallocate(this%condsat)
    deallocate(this%boundname)
    !
    ! -- output table objects
    if (associated(this%outputtab1)) then
      call this%outputtab1%table_da()
      deallocate(this%outputtab1)
      nullify(this%outputtab1)
    end if
    if (associated(this%outputtab2)) then
      call this%outputtab2%table_da()
      deallocate(this%outputtab2)
      nullify(this%outputtab2)
    end if
    !
    ! -- scalars
    call mem_deallocate(this%icellavg)
    call mem_deallocate(this%ivarcv)
    call mem_deallocate(this%idewatcv)
    call mem_deallocate(this%inewton)
    call mem_deallocate(this%ianglex)
    call mem_deallocate(this%icdist)
    call mem_deallocate(this%ingnc)
    call mem_deallocate(this%inmvr)
    call mem_deallocate(this%inobs)
    call mem_deallocate(this%inamedbound)
    call mem_deallocate(this%satomega)
    !
    ! -- return
    return
  end subroutine gwf_gwf_da

  subroutine allocate_arrays(this)
! ******************************************************************************
! allocate_scalars
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    use ConstantsModule, only: LENORIGIN
    ! -- dummy
    class(GwfExchangeType) :: this
    ! -- local
    character(len=LINELENGTH) :: text
    character(len=LENORIGIN) :: origin
    integer(I4B) :: ntabcol
! ------------------------------------------------------------------------------
    !
    ! -- create the origin name
    origin = trim(this%name)
    !
    ! -- Call parent type allocate_scalars
    call this%NumericalExchangeType%allocate_arrays()
    !
    call mem_allocate(this%ihc, this%nexg, 'IHC', origin)
    call mem_allocate(this%cl1, this%nexg, 'CL1', origin)
    call mem_allocate(this%cl2, this%nexg, 'CL2', origin)
    call mem_allocate(this%hwva, this%nexg, 'HWVA', origin)
    call mem_allocate(this%condsat, this%nexg, 'CONDSAT', origin)
    !
    ! -- Allocate boundname
    if(this%inamedbound==1) then
      allocate(this%boundname(this%nexg))
    else
      allocate(this%boundname(1))
    endif
    this%boundname(:) = ''
    !
    ! -- allocate and initialize the output table
    if (this%iprflow /= 0) then
      !
      ! -- dimension table
      ntabcol = 3
      if (this%inamedbound > 0) then
        ntabcol = ntabcol + 1
      end if
      !
      ! -- initialize the output table objects
      !    outouttab1
      call table_cr(this%outputtab1, this%name, '    ')
      call this%outputtab1%table_df(this%nexg, ntabcol, this%gwfmodel1%iout,     &
                                    transient=.TRUE.)
      text = 'NUMBER'
      call this%outputtab1%initialize_column(text, 10, alignment=TABCENTER)
      text = 'CELLID'
      call this%outputtab1%initialize_column(text, 20, alignment=TABLEFT)
      text = 'RATE'
      call this%outputtab1%initialize_column(text, 15, alignment=TABCENTER)
      if (this%inamedbound > 0) then
        text = 'NAME'
        call this%outputtab1%initialize_column(text, 20, alignment=TABLEFT)
      end if
      !    outouttab2
      call table_cr(this%outputtab2, this%name, '    ')
      call this%outputtab2%table_df(this%nexg, ntabcol, this%gwfmodel2%iout,     &
                                    transient=.TRUE.)
      text = 'NUMBER'
      call this%outputtab2%initialize_column(text, 10, alignment=TABCENTER)
      text = 'CELLID'
      call this%outputtab2%initialize_column(text, 20, alignment=TABLEFT)
      text = 'RATE'
      call this%outputtab2%initialize_column(text, 15, alignment=TABCENTER)
      if (this%inamedbound > 0) then
        text = 'NAME'
        call this%outputtab2%initialize_column(text, 20, alignment=TABLEFT)
      end if
    end if
    !
    ! -- return
    return
  end subroutine allocate_arrays

  subroutine gwf_gwf_df_obs(this)
! ******************************************************************************
! gwf_gwf_df_obs
!   -- Store observation type supported by GWF-GWF exchange.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(GwfExchangeType) :: this
    ! -- local
    integer(I4B) :: indx
! ------------------------------------------------------------------------------
    !
    ! -- Store obs type and assign procedure pointer
    !    for gwf-gwf observation type.
    call this%obs%StoreObsType('flow-ja-face', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => gwf_gwf_process_obsID
    !
    ! -- return
    return
  end subroutine gwf_gwf_df_obs

  subroutine gwf_gwf_rp_obs(this)
! ******************************************************************************
! gwf_gwf_rp_obs
!   -- Handle observation IDs that are exchange-boundary names.
!      Store exchange numbers included in observation.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy
    class(GwfExchangeType) :: this
    ! -- local
    integer(I4B) :: i, j, n
    class(ObserveType), pointer :: obsrv => null()
    character(len=LENBOUNDNAME) :: bname
    character(len=1000) :: ermsg
    logical :: jfound
    ! -- formats
10  format('Error: Boundary "',a,'" for observation "',a,               &
           '" is invalid in package "',a,'"')
! ------------------------------------------------------------------------------
    !
    do i=1,this%obs%npakobs
      obsrv => this%obs%pakobs(i)%obsrv
      !
      ! -- indxbnds needs to be deallocated and reallocated (using
      !    ExpandArray) each stress period because list of boundaries
      !    can change each stress period.
      ! -- Not true for exchanges, but leave this in for now anyway.
      if (allocated(obsrv%indxbnds)) then
        deallocate(obsrv%indxbnds)
      endif
      obsrv%BndFound = .false.
      !
      bname = obsrv%FeatureName
      if (bname /= '') then
        ! -- Observation location(s) is(are) based on a boundary name.
        !    Iterate through all boundaries to identify and store
        !    corresponding index(indices) in bound array.
        jfound = .false.
        do j=1,this%nexg
          if (this%boundname(j) == bname) then
            jfound = .true.
            obsrv%BndFound = .true.
            obsrv%CurrentTimeStepEndValue = DZERO
            call ExpandArray(obsrv%indxbnds)
            n = size(obsrv%indxbnds)
            obsrv%indxbnds(n) = j
          endif
        enddo
        if (.not. jfound) then
          write(ermsg,10)trim(bname)
          call store_error(ermsg)
        endif
      else
        ! -- Observation location is a single exchange number
        if (obsrv%intPak1 <= this%nexg) then
          jfound = .true.
          obsrv%BndFound = .true.
          obsrv%CurrentTimeStepEndValue = DZERO
          call ExpandArray(obsrv%indxbnds)
          n = size(obsrv%indxbnds)
          obsrv%indxbnds(n) = obsrv%intPak1
        else
          jfound = .false.
        endif
      endif
    enddo
    !
    if (count_errors() > 0) then
      call store_error_unit(this%inobs)
      call ustop()
    endif
    !
    ! -- Return
    return
  end subroutine gwf_gwf_rp_obs

  subroutine gwf_gwf_fp(this)
! ******************************************************************************
! gwf_gwf_fp
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(GwfExchangeType) :: this
! ------------------------------------------------------------------------------
    !
    return
  end subroutine gwf_gwf_fp
  
  function qcalc(this, iexg, n1, n2)
! ******************************************************************************
! qcalc -- calculate flow between two cells, positive into n1
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- return
    real(DP) :: qcalc
    ! -- dummy
    class(GwfExchangeType) :: this
    integer(I4B), intent(in) :: iexg
    integer(I4B), intent(in) :: n1
    integer(I4B), intent(in) :: n2
    ! -- local
! ------------------------------------------------------------------------------
    !
    ! -- Calculate flow between nodes in the two models
    qcalc = this%cond(iexg) * (this%gwfhalo%x(n2) - this%gwfhalo%x(n1)) !HALO2
    !
    ! -- return
    return
  end function qcalc

  function gwf_gwf_get_iasym(this) result (iasym)
! ******************************************************************************
! gwf_gwf_get_iasym -- return 1 if any option causes the matrix to be asymmetric.
!   Otherwise return 0.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(GwfExchangeType) :: this
    ! -- local
    integer(I4B) :: iasym
! ------------------------------------------------------------------------------
    !
    ! -- Start by setting iasym to zero
    iasym = 0
    !
    ! -- Groundwater flow
    if (this%inewton /= 0) iasym = 1
    !
    ! -- GNC
    if (this%ingnc > 0) then
      if (this%gnc%iasym /= 0) iasym = 1
    endif
    !
    ! -- return
    return
  end function gwf_gwf_get_iasym

  subroutine gwf_gwf_get_m1m2(this, m1, m2) !PAR
! ******************************************************************************
! gwfp_gwfp_get_m1m2
!   -- Get the pointer to gwfmodel1 and gwfmodel2.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    class(gwfExchangeType) :: this
    type(GwfModelType), pointer, intent(out) :: m1, m2
    ! -- local
! ------------------------------------------------------------------------------
    !
    m1 => this%gwfmodel1
    m2 => this%gwfmodel2
    !
    return
  end subroutine gwf_gwf_get_m1m2
  
  subroutine gwf_gwf_save_simvals(this)
! ******************************************************************************
! gwf_gwf_save_simvals
!   -- Calculate observations this time step and call
!      ObsType%SaveOneSimval for each GWF-GWF Type observation.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- dummy
    use SimModule, only: store_error, store_error_unit, ustop
    use ConstantsModule, only: DZERO
    use ObserveModule, only: ObserveType
    class(GwfExchangeType), intent(inout) :: this
    ! -- local
    integer(I4B) :: i, j, n1, n2, nbndobs
    integer(I4B) :: iexg
    real(DP) :: v
    character(len=100) :: msg
    type(ObserveType), pointer :: obsrv => null()
! ------------------------------------------------------------------------------
    !
    ! -- Write simulated values for all gwf-gwf observations
    if (this%obs%npakobs > 0) then
      call this%obs%obs_bd_clear()
      do i = 1, this%obs%npakobs
        obsrv => this%obs%pakobs(i)%obsrv
        nbndobs = size(obsrv%indxbnds)
        do j = 1,  nbndobs
          iexg = obsrv%indxbnds(j)
          v = DZERO
          select case (obsrv%ObsTypeId)
          case ('FLOW-JA-FACE')
            n1 = this%gwfhalo%imapnodem1tohalo(iexg) !HALO2
            n2 = this%gwfhalo%imapnodem2tohalo(iexg) + this%gwfhalo%offset !HALO2
            v = this%cond(iexg) * (this%gwfhalo%x(n2) - this%gwfhalo%x(n1)) !HALO2
            if(this%ingnc > 0) then
              v = v + this%gnc%deltaqgnc(iexg)
            endif
          case default
            msg = 'Error: Unrecognized observation type: ' //                  &
                  trim(obsrv%ObsTypeId)
            call store_error(msg)
            call store_error_unit(this%inobs)
            call ustop()
          end select
          call this%obs%SaveOneSimval(obsrv, v)
        enddo
      enddo
    endif
    !
    return
  end subroutine gwf_gwf_save_simvals

  subroutine gwf_gwf_process_obsID(obsrv, dis, inunitobs, iout)
! ******************************************************************************
! -- This procedure is pointed to by ObsDataType%ProcesssIdPtr. It processes
!    the ID string of an observation definition for GWF-GWF-package observations
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ConstantsModule, only: LINELENGTH
    use InputOutputModule, only: urword
    use ObserveModule, only: ObserveType
    use BaseDisModule, only: DisBaseType
    ! -- dummy
    type(ObserveType),      intent(inout) :: obsrv
    class(DisBaseType), intent(in)    :: dis
    integer(I4B),            intent(in)    :: inunitobs
    integer(I4B),            intent(in)    :: iout
    ! -- local
    integer(I4B) :: n, iexg, istat
    integer(I4B) :: icol, istart, istop
    real(DP) :: r
    character(len=LINELENGTH) :: strng
! ------------------------------------------------------------------------------
    !
    strng = obsrv%IDstring
    icol = 1
    ! -- get exchange index
    call urword(strng, icol, istart, istop, 0, n, r, iout, inunitobs)
    read (strng(istart:istop), '(i10)', iostat=istat) iexg
    if (istat == 0) then
      obsrv%intPak1 = iexg
    else
      ! Integer can't be read from strng; it's presumed to be an exchange
      ! boundary name (already converted to uppercase)
      obsrv%FeatureName = strng(istart:istop)
      ! -- Observation may require summing rates from multiple exchange
      !    boundaries, so assign intPak1 as a value that indicates observation
      !    is for a named exchange boundary or group of exchange boundaries.
      obsrv%intPak1 = NAMEDBOUNDFLAG
    endif
    !
    return
  end subroutine gwf_gwf_process_obsID

  subroutine gwf_cr_halo(this, id, modelname) !HALO2
! ******************************************************************************
! gwf_cr_halo -- Create a new groundwater flow model object for exchange
! Subroutine: (1) creates model object and add to exchange modellist
!             (2) assign values
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use GwfDisuModule,              only: disu_cr
    use GwfNpfModule,               only: npf_cr
    use Xt3dModule,                 only: xt3d_cr
    use GhostNodeModule,            only: gnc_cr
    use GwfMvrModule,               only: mvr_cr
    use GwfIcModule,                only: ic_cr
    ! -- dummy
    type(GwfModelType), pointer   :: this
    integer(I4B), intent(in)      :: id
    character(len=*), intent(in)  :: modelname
    ! -- local
    integer :: in_dum, iout_dum
    ! -- format
! ------------------------------------------------------------------------------
    !
    ! -- Allocate a new GWF Model (this)
    allocate(this)
    call this%allocate_scalars(modelname)
    !
    ! -- Assign values
    this%name = modelname
    this%macronym = 'GWF'
    this%id = id
    !
    ! -- TODO: to be replaced by optional arguments
    in_dum = -1
    iout_dum = -1
    !
    ! -- Create discretization object
    call disu_cr(this%dis, this%name, in_dum, iout_dum)
    !
    ! -- Create packages that are tied directly to model
    call npf_cr(this%npf, this%name, in_dum, iout_dum)
    call xt3d_cr(this%xt3d, this%name, in_dum, iout_dum)
    call gnc_cr(this%gnc, this%name, in_dum, iout_dum)
    call ic_cr(this%ic, this%name, in_dum, iout_dum, this%dis)
    call mvr_cr(this%mvr, this%name, in_dum, iout_dum)
    !
    ! -- return
    return
  end subroutine gwf_cr_halo
  
subroutine gwf_mpi_halo_init() !PAR
! ******************************************************************************
! Initialize the MpiWorld communicator local exchange for gwfhalo
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MpiExchangeGenModule, only: serialrun
    use MpiExchangeModule, only: MpiWorld
    use BaseExchangeModule, only: GetBaseExchangeFromList
    use ListsModule, only: baseexchangelist
    use ArrayHandlersModule, only: ifind
    use ConstantsModule, only: LINELENGTH
    use MemoryManagerModule, only: mem_allocate
    use MpiExchangeColModule, only: mpi_get_distype_str
    ! -- dummy
    ! -- local
    class(BaseExchangeType), pointer :: bp
    type(GwfExchangeType), pointer :: gp
    integer(I4B) :: ic, i, isub, ip, ixp, nex, neq
    integer(I4B), dimension(:), allocatable :: iwrk
    character(len=LENMODELNAME) :: m2name
    character(len=LINELENGTH) :: errmsg
! ------------------------------------------------------------------------------
    if (serialrun) then
      return
    end if
    !
    MpiWorld%nrxp = 0
    allocate(iwrk(MpiWorld%nrproc))
    iwrk = 0
    do ic=1,baseexchangelist%Count()
      bp => GetBaseExchangeFromList(baseexchangelist, ic)
      select type (bp)
      class is (GwfExchangeType)
        gp => bp
      end select
      if (gp%m2_bympi) then
        m2name = gp%gwfhalo%m2name
        i = ifind(MpiWorld%gmodelnames, m2name)
        isub = 0
        if (i > 0) then
          isub = MpiWorld%gsubs(i)
        endif
        if (isub > 0) then
          iwrk(isub) = iwrk(isub) + 1
        else
          write(errmsg,'(a)') 'Program error in gwf_mpi_halo_init.'
          call store_error(errmsg)
          call ustop()
        end if
      end if
    end do
    !
    MpiWorld%nrxp = 0
    do ip = 1, MpiWorld%nrproc
      if (iwrk(ip) > 0) then
        MpiWorld%nrxp = MpiWorld%nrxp + 1
      endif
    end do
    !
    ! -- Allocate local communication data structure
    if (MpiWorld%nrxp > 0) then
      allocate(MpiWorld%lxch(MpiWorld%nrxp))
    endif
    !
    ixp = 0
    do ip = 1, MpiWorld%nrproc
      nex = iwrk(ip)
      if (nex > 0) then
        ixp = ixp + 1
        allocate(MpiWorld%lxch(ixp)%nexchange)
        allocate(MpiWorld%lxch(ixp)%exchange(nex))
        allocate(MpiWorld%lxch(ixp)%xprnk)
        MpiWorld%lxch(ixp)%nexchange = 0
        MpiWorld%lxch(ixp)%xprnk = ip-1
      endif
      ! -- Set mapping to exchange partner index
      iwrk(ip) = ixp
    enddo
    !
    ! -- loop over exchanges and initialize
    do ic=1,baseexchangelist%Count()
      bp => GetBaseExchangeFromList(baseexchangelist, ic)
      select type (bp)
      class is (GwfExchangeType)
        gp => bp
      end select
      if (gp%m2_bympi) then
        m2name = gp%gwfhalo%m2name
        i = ifind(MpiWorld%gmodelnames, m2name)
        isub = MpiWorld%gsubs(i)
        ixp = iwrk(isub)
        nex = MpiWorld%lxch(ixp)%nexchange
        nex = nex + 1
        ! -- set pointer to exchange
        MpiWorld%lxch(ixp)%exchange(nex)%name = gp%name
        MpiWorld%lxch(ixp)%exchange(nex)%halo_name = gp%gwfhalo%name
        MpiWorld%lxch(ixp)%exchange(nex)%m1_name = gp%gwfhalo%m1name
        MpiWorld%lxch(ixp)%exchange(nex)%m2_name = gp%gwfhalo%m2name
        call mpi_get_distype_str(gp%gwfhalo%m1name, MpiWorld%lxch(ixp)%exchange(nex)%m1_dis)
        call mpi_get_distype_str(gp%gwfhalo%m2name, MpiWorld%lxch(ixp)%exchange(nex)%m2_dis)
        MpiWorld%lxch(ixp)%nexchange = nex
      end if
    end do
    !
    MpiWorld%linit = .true.
    !
    ! -- add variables
    call MpiWorld%mpi_add_vg('HALO_INIT_CON_S')
    call MpiWorld%mpi_add_vmt('HALO_INIT_CON_S', 'NODES', 'CON', 'HLL', 'ALL', 'HLL', 'MEM')
    call MpiWorld%mpi_add_vmt('HALO_INIT_CON_S', 'NJA', 'CON', 'HLL', 'ALL', 'HLL', 'MEM')
    call MpiWorld%mpi_add_vmt('HALO_INIT_CON_S', 'NJAS', 'CON', 'HLL', 'ALL', 'HLL', 'MEM')
    call MpiWorld%mpi_init_vg('HALO_INIT_CON_S')
    
    call MpiWorld%mpi_add_vg('HALO_INIT_CON_A')
    call MpiWorld%mpi_add_vmt('HALO_INIT_CON_A', 'NBNODES', '','HLL', 'ALL', 'HLL', 'MEM')
    call MpiWorld%mpi_add_vmt('HALO_INIT_CON_A', 'IMAPNODEMTOHALO', '',  'HLL', 'ALL', 'HLL', 'MEM')
    call MpiWorld%mpi_add_vmt('HALO_INIT_CON_A', 'IMAPMTOHALO', '',  'HLL', 'ALL', 'HLL', 'MEM')
    call MpiWorld%mpi_add_vmt('HALO_INIT_CON_A', 'IA','CON', 'HLL', 'ALL', 'HLL', 'MEM')
    call MpiWorld%mpi_add_vmt('HALO_INIT_CON_A', 'JA','CON', 'HLL', 'ALL', 'HLL', 'MEM')
    call MpiWorld%mpi_add_vmt('HALO_INIT_CON_A', 'JAS','CON', 'HLL', 'ALL', 'HLL', 'MEM')
    call MpiWorld%mpi_add_vmt('HALO_INIT_CON_A', 'IHC','CON', 'HLL', 'ALL', 'HLL', 'MEM')
    call MpiWorld%mpi_add_vmt('HALO_INIT_CON_A', 'CL1','CON', 'HLL', 'ALL', 'HLL', 'MEM')
    call MpiWorld%mpi_add_vmt('HALO_INIT_CON_A', 'CL2','CON', 'HLL', 'ALL', 'HLL', 'MEM')
    call MpiWorld%mpi_add_vmt('HALO_INIT_CON_A', 'HWVA','CON', 'HLL', 'ALL', 'HLL', 'MEM')
    call MpiWorld%mpi_init_vg('HALO_INIT_CON_A')
    !
    call MpiWorld%mpi_add_vg('HALO_INIT_DIS')
    call MpiWorld%mpi_add_vmt('HALO_INIT_DIS', 'TOP',  'DIS', 'GWF', 'HM1', 'HAL', 'HAL')
    call MpiWorld%mpi_add_vmt('HALO_INIT_DIS', 'BOT',  'DIS', 'GWF', 'HM1', 'HAL', 'HAL')
    call MpiWorld%mpi_add_vmt('HALO_INIT_DIS', 'AREA', 'DIS', 'GWF', 'HM1', 'HAL', 'HAL')
    call MpiWorld%mpi_init_vg('HALO_INIT_DIS')
    
    call MpiWorld%mpi_add_vg('HALO_INIT_NPFIC')
    call MpiWorld%mpi_add_vmt('HALO_INIT_NPFIC', 'ICELLTYPE', 'NPF', 'GWF', 'HM1', 'HAL', 'HAL')
    call MpiWorld%mpi_add_vmt('HALO_INIT_NPFIC', 'K11',       'NPF', 'GWF', 'HM1', 'HAL', 'HAL')
    call MpiWorld%mpi_add_vmt('HALO_INIT_NPFIC', 'K22',       'NPF', 'GWF', 'HM1', 'HAL', 'HAL')
    call MpiWorld%mpi_add_vmt('HALO_INIT_NPFIC', 'K33',       'NPF', 'GWF', 'HM1', 'HAL', 'HAL')
    call MpiWorld%mpi_add_vmt('HALO_INIT_NPFIC', 'WETDRY',    'NPF', 'GWF', 'HM1', 'HAL', 'HAL')
    call MpiWorld%mpi_add_vmt('HALO_INIT_NPFIC', 'ANGLE1',    'NPF', 'GWF', 'HM1', 'HAL', 'HAL')
    call MpiWorld%mpi_add_vmt('HALO_INIT_NPFIC', 'ANGLE2',    'NPF', 'GWF', 'HM1', 'HAL', 'HAL')
    call MpiWorld%mpi_add_vmt('HALO_INIT_NPFIC', 'ANGLE3',    'NPF', 'GWF', 'HM1', 'HAL', 'HAL')
    call MpiWorld%mpi_add_vmt('HALO_INIT_NPFIC', 'SAT',       'NPF', 'GWF', 'HM1', 'HAL', 'HAL')
    call MpiWorld%mpi_add_vmt('HALO_INIT_NPFIC', 'STRT',      'IC',  'GWF', 'HM1', 'HAL', 'HAL')
    call MpiWorld%mpi_init_vg('HALO_INIT_NPFIC')
    !
    ! -- local exchange of connection data scalars
    call MpiWorld%mpi_local_exchange('', 'HALO_INIT_CON_S', .true.) !PAR
    !
    ! Allocate connection array for M2
    do ic=1,baseexchangelist%Count()
      bp => GetBaseExchangeFromList(baseexchangelist, ic)
      select type (bp)
      class is (GwfExchangeType)
        gp => bp
      end select
      if (gp%m2_bympi) then
        call gp%gwfhalo%m2con%allocate_arrays()
        call mem_allocate(gp%gwfhalo%m2nbnod, gp%gwfhalo%nband, 'NBNODES', trim(gp%gwfhalo%name)//'_M2')
        call mem_allocate(gp%gwfhalo%imapnodem2tohalo, gp%gwfhalo%nexg, 'IMAPNODEMTOHALO', trim(gp%gwfhalo%name)//'_M2')
        neq = gp%gwfhalo%m2con%nodes
        call mem_allocate(gp%gwfhalo%imapm2tohalo, neq, 'IMAPMTOHALO', trim(gp%gwfhalo%name)//'_M2')
        call mem_allocate(gp%gwfhalo%m2nodes, neq, 'MNODES', trim(gp%gwfhalo%name)//'_M2')
      end if
    end do
    !
    ! -- set arrays
    call MpiWorld%mpi_local_exchange('', 'HALO_INIT_CON_A', .true.) !PAR
    
    deallocate(iwrk)
    !
    ! -- return
    return
end subroutine gwf_mpi_halo_init
  
end module GwfGwfExchangeModule