  MODULE IMSLinearModule
  
  use KindModule, only: DP, I4B
  use ConstantsModule, only: LINELENGTH, LENSOLUTIONNAME, LENMEMPATH,          &
                             IZERO, DZERO, DPREC, DSAME,                       &
                             DEM8, DEM6, DEM5, DEM4, DEM3, DEM2, DEM1,         &
                             DHALF, DONE, DTWO,                                &
                             VDEBUG
  use GenericUtilitiesModule, only: sim_message, IS_SAME
  use IMSReorderingModule, only: ims_genrcm, ims_odrv, ims_dperm, ims_vperm
  use BlockParserModule, only: BlockParserType
  use MpiExchangeModule,  only: MpiExchangeType !PAR

  IMPLICIT NONE
  private
  
  TYPE, PUBLIC :: ImsLinearDataType
    character(len=LENMEMPATH), pointer :: memoryPath !PAR                !< the path for storing variables in the memory manager
    integer(I4B), POINTER :: iout => NULL()
    integer(I4B), POINTER :: IPRIMS => NULL()
    integer(I4B), POINTER :: ILINMETH => NULL()
    integer(I4B), POINTER :: ITER1 => NULL()
    integer(I4B), POINTER :: IPC => NULL()
    integer(I4B), POINTER :: ISCL => NULL()
    integer(I4B), POINTER :: IORD => NULL()
    integer(I4B), POINTER :: NORTH => NULL()
    integer(I4B), POINTER :: ICNVGOPT => NULL()
    integer(I4B), POINTER :: IACPC => NULL()
    integer(I4B), POINTER :: NITERC => NULL()
    integer(I4B), POINTER :: NIABCGS => NULL()
    integer(I4B), POINTER :: NIAPC => NULL()
    integer(I4B), POINTER :: NJAPC => NULL()
    real(DP), POINTER :: DVCLOSE => NULL()
    real(DP), POINTER :: RCLOSE => NULL()
    real(DP), POINTER :: RELAX => NULL()
    real(DP), POINTER :: EPFACT => NULL()
    real(DP), POINTER :: L2NORM0 => NULL()
    real(DP), POINTER :: RHOTOL => NULL() !SOL
    real(DP), POINTER :: ALPHATOL => NULL() !SOL
    real(DP), POINTER :: OMEGATOL => NULL() !SOL
    integer(I4B), POINTER :: IRCLOSEPRECHK => NULL() !SOL
    ! ILUT VARIABLES
    integer(I4B), POINTER :: LEVEL => NULL()
    real(DP), POINTER :: DROPTOL => NULL()
    integer(I4B), POINTER :: NJLU => NULL()
    integer(I4B), POINTER :: NJW => NULL()
    integer(I4B), POINTER :: NWLU => NULL()
    ! POINTERS TO SOLUTION VARIABLES
    integer(I4B), POINTER :: NEQ => NULL()
    integer(I4B), POINTER :: NJA => NULL()
    integer(I4B), dimension(:), pointer, contiguous :: IA => NULL()
    integer(I4B), dimension(:), pointer, contiguous :: JA => NULL()
    real(DP), dimension(:), pointer, contiguous :: AMAT => NULL()
    real(DP), dimension(:), pointer, contiguous :: RHS => NULL()
    real(DP), dimension(:), pointer, contiguous :: X => NULL()
    integer(I4B), pointer :: IBJFLAG => NULL() !BJ
    ! VECTORS
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: DSCALE => NULL()
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: DSCALE2 => NULL()
    integer(I4B), POINTER,DIMENSION(:),CONTIGUOUS :: IAPC => NULL()
    integer(I4B), POINTER,DIMENSION(:),CONTIGUOUS :: JAPC => NULL()
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: APC => NULL()
    integer(I4B), POINTER, DIMENSION(:), CONTIGUOUS :: LORDER => NULL()
    integer(I4B), POINTER, DIMENSION(:), CONTIGUOUS :: IORDER => NULL()
    integer(I4B), POINTER, DIMENSION(:), CONTIGUOUS :: IARO => NULL()
    integer(I4B), POINTER, DIMENSION(:), CONTIGUOUS :: JARO => NULL()
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: ARO => NULL()
    ! WORKING ARRAYS
    integer(I4B), POINTER, DIMENSION(:), CONTIGUOUS :: IW => NULL()
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: W => NULL()
    integer(I4B), POINTER, DIMENSION(:), CONTIGUOUS :: ID => NULL()
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: D => NULL()
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: P => NULL()
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: Q => NULL()
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: Z => NULL()
    ! BICGSTAB WORKING ARRAYS
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: T => NULL()
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: V => NULL()
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: DHAT => NULL()
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: PHAT => NULL()
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: QHAT => NULL()
    ! POINTERS FOR USE WITH BOTH ORIGINAL AND RCM ORDERINGS
    integer(I4B), POINTER, DIMENSION(:), CONTIGUOUS :: IA0 => NULL()
    integer(I4B), POINTER, DIMENSION(:), CONTIGUOUS :: JA0 => NULL()
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: A0 => NULL()
    ! ILUT WORKING ARRAYS
    integer(I4B),POINTER,DIMENSION(:),CONTIGUOUS :: JLU => NULL()
    integer(I4B),POINTER,DIMENSION(:),CONTIGUOUS :: JW => NULL()
    real(DP),POINTER,DIMENSION(:),CONTIGUOUS::WLU => NULL()
    type(MpiExchangeType), pointer:: MpiSol => NULL() !PAR
    ! COARSE GRID CORRECTION ARRAYS
    integer(I4B), pointer                                :: icgc     => NULL() !CGC
    integer(I4B), pointer                                :: cgcsol   => NULL() !CGC
    !
    integer(I4B), pointer                                :: cgcgneq  => NULL() !CGC
    integer(I4B), pointer                                :: cgcgnia  => NULL() !CGC
    integer(I4B), pointer                                :: cgcgnja  => NULL() !CGC
    integer(I4B), dimension(:), pointer, contiguous      :: cgcgia   => NULL() !CGC
    integer(I4B), dimension(:), pointer, contiguous      :: cgcgja   => NULL() !CGC
    integer(I4B), dimension(:), pointer, contiguous      :: cgcgiapc => NULL() !CGC
    integer(I4B), dimension(:), pointer, contiguous      :: cgcgjapc => NULL() !CGC
    !
    integer(I4B), pointer                                :: cgclnia  => NULL() !CGC
    integer(I4B), pointer                                :: cgclnja  => NULL() !CGC
    integer(I4B), dimension(:), pointer, contiguous      :: cgclia   => NULL() !CGC
    integer(I4B), dimension(:), pointer, contiguous      :: cgclja   => NULL() !CGC
    integer(I4B), dimension(:), pointer, contiguous      :: cgclgja  => NULL() !CGC
    !
    integer(I4B), dimension(:), pointer, contiguous      :: cgcl2gid => NULL() !CGC: local  --> global mapping
    integer(I4B), dimension(:), pointer, contiguous      :: cgcg2lid => NULL() !CGC: global --> local  mapping
    integer(I4B), dimension(:), pointer, contiguous      :: cgcg2cid => NULL() !CGC: global --> coarse mapping
    integer(I4B), dimension(:), pointer, contiguous      :: cgcc2gid => NULL() !CGC: coarse --> global mapping
    integer(I4B), dimension(:), pointer, contiguous      :: cgccizc  => NULL() !CGC
    !
    real(DP), pointer, dimension(:,:), contiguous        :: acpc  => NULL() !CGC
    real(DP), pointer, dimension(:), contiguous          :: acpci => NULL() !CGC
    real(DP), pointer, dimension(:), contiguous          :: xc    => NULL() !CGC
    real(DP), pointer, dimension(:), contiguous          :: rhsc  => NULL() !CGC
    real(DP), pointer, dimension(:), contiguous          :: wc    => NULL() !CGC
    real(DP), pointer, dimension(:), contiguous          :: zc    => NULL() !CGC
    real(DP), dimension(:,:), pointer, contiguous        :: a0zc  => NULL() !CGC
    integer(I4B), pointer                                :: acpcb => NULL() !CGC
    ! PROCEDURES (METHODS)
    CONTAINS
      PROCEDURE :: IMSLINEAR_ALLOCATE => IMSLINEAR_AR
      procedure :: imslinear_summary
      PROCEDURE :: IMSLINEAR_ALLOCATE_CGC => IMSLINEAR_AR_CGC ! CGC
      PROCEDURE :: IMSLINEAR_APPLY => IMSLINEAR_AP
      procedure :: IMSLINEAR_DA
      procedure, private :: allocate_scalars
      ! -- PRIVATE PROCEDURES
      PROCEDURE, PRIVATE :: SET_IMSLINEAR_INPUT
      procedure, private :: IMSLINEARSUB_CG !CGC
      procedure, private :: IMSLINEARSUB_BCGS !CGC
      procedure, private :: IMSLINEARSUB_PCC !CGC
      procedure, private :: IMSLINEARSUB_PCC_SOL !CGC
  END TYPE ImsLinearDataType
  
  type(BlockParserType), private :: parser

  
  CONTAINS
    SUBROUTINE IMSLINEAR_AR(THIS, NAME, IN, IOUT, IPRIMS, MXITER, IFDPARAM, &
                            IMSLINEARM, NEQ, NJA, IA, JA, AMAT, RHS, X,     &
                            NINNER, &
                            IBJFLAG, & !BJ
                            LFINDBLOCK)
!     ******************************************************************
!     ALLOCATE STORAGE FOR PCG ARRAYS AND READ IMSLINEAR DATA
!     ******************************************************************
!
!        SPECIFICATIONS:
!     ------------------------------------------------------------------
      use MemoryManagerModule, only: mem_allocate
      use MemoryHelperModule,  only: create_mem_path
      use SimModule, only: ustop, store_error, count_errors,            &
                           deprecation_warning
      use MpiExchangeGenModule, only: parallelrun !PAR
      !IMPLICIT NONE
!     + + + DUMMY VARIABLES + + +
      CLASS(ImsLinearDataType), INTENT(INOUT) :: THIS
      CHARACTER (LEN=LENSOLUTIONNAME), INTENT(IN) :: NAME
      integer(I4B), INTENT(IN) :: IN
      integer(I4B), INTENT(IN) :: IOUT
      integer(I4B), TARGET, INTENT(IN) :: IPRIMS
      integer(I4B), INTENT(IN) :: MXITER
      integer(I4B), INTENT(IN) :: IFDPARAM
      integer(I4B), INTENT(INOUT) :: IMSLINEARM
      integer(I4B), TARGET, INTENT(IN) :: NEQ
      integer(I4B), TARGET, INTENT(IN) :: NJA
      integer(I4B), DIMENSION(NEQ+1), TARGET, INTENT(IN) :: IA
      integer(I4B), DIMENSION(NJA), TARGET, INTENT(IN) :: JA
      real(DP), DIMENSION(NJA), TARGET, INTENT(IN) :: AMAT
      real(DP), DIMENSION(NEQ), TARGET, INTENT(INOUT) :: RHS
      real(DP), DIMENSION(NEQ), TARGET, INTENT(INOUT) :: X
      integer(I4B), TARGET, INTENT(INOUT) :: NINNER
      integer(I4B), INTENT(IN), TARGET :: IBJFLAG !BJ
      integer(I4B), INTENT(IN), OPTIONAL :: LFINDBLOCK
!     + + + LOCAL VARIABLES + + +
      LOGICAL :: lreaddata
      character(len=LINELENGTH) :: errmsg
      character(len=LINELENGTH) :: warnmsg
      character(len=LINELENGTH) :: keyword
      integer(I4B) :: i, n
      integer(I4B) :: i0
      integer(I4B) :: iscllen, iolen
      integer(I4B) :: ierr
      real(DP) :: r
      logical :: isfound, endOfBlock
      integer(I4B) :: ijlu
      integer(I4B) :: ijw
      integer(I4B) :: iwlu
      integer(I4B) :: iwk
!     + + + PARAMETERS + + +
!     + + + OUTPUT FORMATS + + +
!------------------------------------------------------------------
!
!-------SET LREADDATA
      IF (PRESENT(LFINDBLOCK)) THEN
        IF (LFINDBLOCK < 1) THEN
          lreaddata = .FALSE.
        ELSE
          lreaddata = .TRUE.
        END IF
      ELSE
        lreaddata = .TRUE.
      END IF
!
!-------DEFINE NAME      
      ALLOCATE(THIS%memoryPath) !CGC
      this%memoryPath = create_mem_path(name, 'IMSLINEAR') !PAR
!
!-------SET POINTERS TO SOLUTION STORAGE
      THIS%IPRIMS => IPRIMS
      THIS%NEQ => NEQ
      THIS%NJA => NJA
      THIS%IA => IA
      THIS%JA => JA
      THIS%AMAT => AMAT
      THIS%RHS => RHS
      THIS%X => X
      THIS%IBJFLAG => IBJFLAG !BJ
!-------ALLOCATE SCALAR VARIABLES
      call this%allocate_scalars()
!
!-------initialize iout
      this%iout = iout
!
!-------DEFAULT VALUES
      THIS%IORD = 0
      THIS%ISCL = 0
      THIS%IPC = 0
      THIS%LEVEL = 0
      this%cgcsol = 1 !CGC
      
      !clevel = ''
      !cdroptol = ''
!
!-------TRANSFER COMMON VARIABLES FROM IMS TO IMSLINEAR
      THIS%ILINMETH = 0
      
      THIS%IACPC = 0
      THIS%RELAX = DZERO !0.97
      
      THIS%DROPTOL = DZERO
      
      THIS%NORTH = 0

      THIS%ICNVGOPT = 0
!
!-------PRINT A MESSAGE IDENTIFYING IMSLINEAR SOLVER PACKAGE
      write(iout,2000)
02000 FORMAT (1X,/1X,'IMSLINEAR -- UNSTRUCTURED LINEAR SOLUTION',               &
     &        ' PACKAGE, VERSION 8, 04/28/2017')
!
!-------SET DEFAULT IMSLINEAR PARAMETERS
      CALL THIS%SET_IMSLINEAR_INPUT(IFDPARAM)
      NINNER = this%iter1
!
!-------Initialize block parser
      call parser%Initialize(in, iout)
!
! -- get IMSLINEAR block
      if (lreaddata) then
        call parser%GetBlock('LINEAR', isfound, ierr, &
          supportOpenClose=.true., blockRequired=.FALSE.)
      else
        isfound = .FALSE.
      end if
!
! -- parse IMSLINEAR block if detected
      if (isfound) then
        write(iout,'(/1x,a)')'PROCESSING LINEAR DATA'
        do
          call parser%GetNextLine(endOfBlock)
          if (endOfBlock) exit
          call parser%GetStringCaps(keyword)
          ! -- parse keyword
          select case (keyword)
            case ('INNER_DVCLOSE')
              this%DVCLOSE = parser%GetDouble()
            case ('INNER_RCLOSE')
              this%rclose = parser%GetDouble()
              ! -- look for additional key words
              call parser%GetStringCaps(keyword)
              if (keyword == 'STRICT') then
                THIS%ICNVGOPT = 1
              else if (keyword == 'L2NORM_RCLOSE') then
                THIS%ICNVGOPT = 2
              else if (keyword == 'RELATIVE_RCLOSE') then
                THIS%ICNVGOPT = 3
              else if (keyword == 'L2NORM_RELATIVE_RCLOSE') then
                THIS%ICNVGOPT = 4
              end if
            case ('INNER_RCLOSE_PRE_CHECK') !SOL
              this%ircloseprechk = parser%GetInteger() !SOL
            case ('INNER_MAXIMUM')
              i = parser%GetInteger()
              this%iter1 = i
              NINNER = i
            case ('LINEAR_ACCELERATION')
              call parser%GetStringCaps(keyword)
              if (keyword.eq.'CG') then
                THIS%ILINMETH = 1
              else if (keyword.eq.'BICGSTAB') then
                THIS%ILINMETH = 2
              else
                THIS%ILINMETH = 0
                write(errmsg,'(3a)')                                             &
                  'UNKNOWN IMSLINEAR LINEAR_ACCELERATION METHOD (',              &
                  trim(keyword), ').'
                call store_error(errmsg)
              end if
            case ('SCALING_METHOD')
              call parser%GetStringCaps(keyword)
              i = 0
              if (keyword.eq.'NONE') then
                i = 0
              else if (keyword.eq.'DIAGONAL') then
                i = 1
              else if (keyword.eq.'L2NORM') then
                i = 2
              else
                write(errmsg,'(3a)')                                             &
                  'UNKNOWN IMSLINEAR SCALING_METHOD (', trim(keyword), ').'
                call store_error(errmsg)
              end if
              THIS%ISCL = i
            case ('RED_BLACK_ORDERING')
              i = 0
            case ('REORDERING_METHOD')
              call parser%GetStringCaps(keyword)
              i = 0
              if (keyword == 'NONE') then
                i = 0
              else if (keyword == 'RCM') then
                i = 1
              else if (keyword == 'MD') then
                i = 2
              else
                write(errmsg,'(3a)')                                             &
                  'UNKNOWN IMSLINEAR REORDERING_METHOD (', trim(keyword), ').'
                call store_error(errmsg)
              end if
              THIS%IORD = i
            case ('NUMBER_ORTHOGONALIZATIONS')
              this%north = parser%GetInteger()
            case ('RELAXATION_FACTOR')
              this%relax = parser%GetDouble()
            case ('PRECONDITIONER_LEVELS')
              i = parser%GetInteger()
              this%level = i
              if (i < 0) then
                write(errmsg,'(a,1x,a)')                                         &
                  'IMSLINEAR PRECONDITIONER_LEVELS MUST BE GREATER THAN',        &
                  'OR EQUAL TO ZERO'
                call store_error(errmsg)
              end if
              !write(clevel, '(i15)') i
            case ('PRECONDITIONER_DROP_TOLERANCE')
              r = parser%GetDouble()
              THIS%DROPTOL = r
              if (r < DZERO) then
                write(errmsg,'(a,1x,a)')                                         &
                  'IMSLINEAR PRECONDITIONER_DROP_TOLERANCE',                     &
                  'MUST BE GREATER THAN OR EQUAL TO ZERO'
                call store_error(errmsg)
              end if
            !
            ! -- deprecated variables
            case ('INNER_HCLOSE')
              this%DVCLOSE = parser%GetDouble()
              !
              ! -- create warning message
              write(warnmsg,'(a)')                                               &
                'SETTING INNER_DVCLOSE TO INNER_HCLOSE VALUE'
              !
              ! -- create deprecation warning
              call deprecation_warning('LINEAR', 'INNER_HCLOSE', '6.1.1',        &
                                       warnmsg, parser%GetUnit())
            !
           case ('INNER_RHO_TOL') !SOL
              r = parser%GetDouble() !SOL
              this%rhotol = r !SOL
              if (r < DZERO) then !SOL
                write(errmsg,'(4x,a,a)') & !SOL
     &            '****ERROR. INNER_RHO_TOL: ', & !SOL
     &            'MUST BE GREATER THAN OR EQUAL TO ZERO' !SOL
                call store_error(errmsg) !SOL
              end if !SOL
            case ('INNER_ALPHA_TOL') !SOL
              r = parser%GetDouble() !SOL
              this%alphatol = r !SOL
              if (r < DZERO) then !SOL
                write(errmsg,'(4x,a,a)') & !SOL
     &            '****ERROR. INNER_ALPHA_TOL: ', & !SOL
     &            'MUST BE GREATER THAN OR EQUAL TO ZERO' !SOL
                call store_error(errmsg) !SOL
              end if !SOL
            case ('INNER_OMEGA_TOL') !SOL
              r = parser%GetDouble() !SOL
              this%omegatol = r !SOL
              if (r < DZERO) then !SOL
                write(errmsg,'(4x,a,a)') & !SOL
     &            '****ERROR. INNER_OMEGA_TOL: ', & !SOL
     &            'MUST BE GREATER THAN OR EQUAL TO ZERO' !SOL
                call store_error(errmsg) !SOL
              end if !SOL
            case ('CGC_SOLVER') !CGC
              i = parser%GetInteger() !CGC
              this%cgcsol = i !CGC
              if ((i < 0).or.(i > 2)) then !CGC
                write(errmsg,'(4x,a,a)') & !CGC
     &            '****ERROR. CGC_SOLVER: ', & !CGC
     &            'MUST BE 1 FOR BAND LU OR 2 FOR ILU(0)' !CGC
                call store_error(errmsg) !CGC
              end if !CGC
            ! -- default
            case default
                write(errmsg,'(3a)')                                             &
                  'UNKNOWN IMSLINEAR KEYWORD (', trim(keyword), ').'
              call store_error(errmsg)
          end select
        end do
        write(iout,'(1x,a)') 'END OF LINEAR DATA'
      else
        if (IFDPARAM ==  0) THEN
          write(errmsg,'(a)') 'NO LINEAR BLOCK DETECTED.'
          call store_error(errmsg)
        end if
      end if
      
      IMSLINEARM = THIS%ILINMETH
!
!-------DETERMINE PRECONDITIONER
      IF (THIS%LEVEL > 0 .OR. THIS%DROPTOL > DZERO) THEN
        THIS%IPC = 3
      ELSE
        THIS%IPC = 1
      END IF
      IF (THIS%RELAX > DZERO) THEN
        THIS%IPC = THIS%IPC + 1
      END IF
!
!-------ERROR CHECKING FOR OPTIONS
      IF (THIS%ISCL < 0 ) THIS%ISCL = 0
      IF (THIS%ISCL > 2 ) THEN
        WRITE(errmsg,'(A)') 'IMSLINEAR7AR ISCL MUST BE <= 2'
        call store_error(errmsg)
      END IF
      IF (THIS%IORD < 0 ) THIS%IORD = 0
      IF (THIS%IORD > 2) THEN
        WRITE(errmsg,'(A)') 'IMSLINEAR7AR IORD MUST BE <= 2'
        call store_error(errmsg)
      END IF
      IF (THIS%NORTH < 0) THEN
        WRITE(errmsg,'(A)') 'IMSLINEAR7AR NORTH MUST >= 0'
        call store_error(errmsg)
      END IF
      IF (THIS%RCLOSE == DZERO) THEN
        IF (THIS%ICNVGOPT /= 3) THEN
          WRITE(errmsg,'(A)') 'IMSLINEAR7AR RCLOSE MUST > 0.0'
          call store_error(errmsg)
        END IF
      END IF
      IF (THIS%RELAX < DZERO) THEN
        WRITE(errmsg,'(A)') 'IMSLINEAR7AR RELAX MUST BE >= 0.0'
        call store_error(errmsg)
      END IF
      IF (THIS%RELAX > DONE) THEN
        WRITE(errmsg,'(A)') 'IMSLINEAR7AR RELAX MUST BE <= 1.0'
        call store_error(errmsg)
      END IF
      IF (THIS%IBJFLAG.EQ.1) THEN !BJ
        IF (THIS%IPC > 2) THEN !BJ
          WRITE( errmsg,'(A)' ) 'IMSLINEAR7AR: BLOCK JACOBI '//     & !BJ
          'OPTION NOT YET SUPPORTED FOR USED PRECONDITIONER' !BJ
          call store_error(errmsg) !BJ
        END IF !BJ
      END IF !BJ
!
!-------CHECK FOR ERRORS IN IMSLINEAR      
      if (count_errors() > 0) then
        call parser%StoreErrorUnit()
        call ustop()
      endif
!
!-------INITIALIZE IMSLINEAR VARIABLES
      THIS%NITERC = 0
!
!-------ALLOCATE AND INITIALIZE MEMORY FOR IMSLINEAR
      iscllen  = 1
      IF (THIS%ISCL.NE.0 ) iscllen  = NEQ
      CALL mem_allocate(THIS%DSCALE, iscllen, 'DSCALE', TRIM(THIS%memoryPath))
      CALL mem_allocate(THIS%DSCALE2, iscllen, 'DSCALE2', TRIM(THIS%memoryPath))
      
!-------ALLOCATE MEMORY FOR PRECONDITIONING MATRIX
      ijlu      = 1
      ijw       = 1
      iwlu      = 1
      ! -- ILU0 AND MILU0
      THIS%NIAPC = THIS%NEQ
      THIS%NJAPC = THIS%NJA
      ! -- ILUT AND MILUT
      IF (THIS%IPC ==  3 .OR. THIS%IPC ==  4) THEN
        THIS%NIAPC = THIS%NEQ
        IF (THIS%LEVEL > 0) THEN
          iwk      = THIS%NEQ * (THIS%LEVEL * 2 + 1)
        ELSE
          iwk = 0
          DO n = 1, NEQ
            i = IA(n+1) - IA(n)
            IF (i > iwk) THEN
              iwk = i
            END IF
          END DO
          iwk      = THIS%NEQ * iwk
        END IF
        THIS%NJAPC = iwk
        ijlu       = iwk
        ijw        = 2 * THIS%NEQ
        iwlu       = THIS%NEQ + 1
      END IF
      THIS%NJLU = ijlu
      THIS%NJW  = ijw
      THIS%NWLU = iwlu
!-------ALLOCATE BASE PRECONDITIONER VECTORS
      CALL mem_allocate(THIS%IAPC, THIS%NIAPC+1, 'IAPC', TRIM(THIS%memoryPath))
      CALL mem_allocate(THIS%JAPC, THIS%NJAPC, 'JAPC', TRIM(THIS%memoryPath))
      CALL mem_allocate(THIS%APC, THIS%NJAPC, 'APC', TRIM(THIS%memoryPath))
!-------ALLOCATE MEMORY FOR ILU0 AND MILU0 NON-ZERO ROW ENTRY VECTOR
      CALL mem_allocate(THIS%IW, THIS%NIAPC, 'IW', TRIM(THIS%memoryPath))
      CALL mem_allocate(THIS%W, THIS%NIAPC, 'W', TRIM(THIS%memoryPath))
!-------ALLOCATE MEMORY FOR ILUT VECTORS
      CALL mem_allocate(THIS%JLU, ijlu, 'JLU', TRIM(THIS%memoryPath))
      CALL mem_allocate(THIS%JW, ijw, 'JW', TRIM(THIS%memoryPath))
      CALL mem_allocate(THIS%WLU, iwlu, 'WLU', TRIM(THIS%memoryPath))
!-------GENERATE IAPC AND JAPC FOR ILU0 AND MILU0
      IF (THIS%IPC ==  1 .OR. THIS%IPC ==  2) THEN
        CALL IMSLINEARSUB_PCCRS(THIS%NEQ,THIS%NJA,THIS%IA,THIS%JA,              &
                                THIS%IAPC,THIS%JAPC)
      END IF
!-------ALLOCATE SPACE FOR PERMUTATION VECTOR
      i0     = 1
      iolen  = 1
      IF (THIS%IORD.NE.0) THEN
        i0     = THIS%NEQ
        iolen  = THIS%NJA
      END IF
      CALL mem_allocate(THIS%LORDER, i0, 'LORDER', TRIM(THIS%memoryPath))
      CALL mem_allocate(THIS%IORDER, i0, 'IORDER', TRIM(THIS%memoryPath))
      CALL mem_allocate(THIS%IARO, i0+1, 'IARO', TRIM(THIS%memoryPath))
      CALL mem_allocate(THIS%JARO, iolen, 'JARO', TRIM(THIS%memoryPath))
      CALL mem_allocate(THIS%ARO, iolen, 'ARO', TRIM(THIS%memoryPath))
!-------ALLOCATE WORKING VECTORS FOR IMSLINEAR SOLVER
      CALL mem_allocate(THIS%ID, THIS%NEQ, 'ID', TRIM(THIS%memoryPath))
      CALL mem_allocate(THIS%D, THIS%NEQ, 'D', TRIM(THIS%memoryPath))
      CALL mem_allocate(THIS%P, THIS%NEQ, 'P', TRIM(THIS%memoryPath))
      CALL mem_allocate(THIS%Q, THIS%NEQ, 'Q', TRIM(THIS%memoryPath))
      CALL mem_allocate(THIS%Z, THIS%NEQ, 'Z', TRIM(THIS%memoryPath))
!-------ALLOCATE MEMORY FOR BCGS WORKING ARRAYS
      THIS%NIABCGS = 1
      IF (THIS%ILINMETH ==  2) THEN
        THIS%NIABCGS = THIS%NEQ
      END IF
      CALL mem_allocate(THIS%T, THIS%NIABCGS, 'T', TRIM(THIS%memoryPath))
      CALL mem_allocate(THIS%V, THIS%NIABCGS, 'V', TRIM(THIS%memoryPath))
      CALL mem_allocate(THIS%DHAT, THIS%NIABCGS, 'DHAT', TRIM(THIS%memoryPath))
      CALL mem_allocate(THIS%PHAT, THIS%NIABCGS, 'PHAT', TRIM(THIS%memoryPath))
      CALL mem_allocate(THIS%QHAT, THIS%NIABCGS, 'QHAT', TRIM(THIS%memoryPath))
      ! -- add variables for MPI point-to-point
      if (parallelrun) then !PAR
        call this%MpiSol%mpi_add_vg('IMS-X') !PAR
        call this%MpiSol%mpi_add_vmt('IMS-X','X','', 'SOL', 'MM1', 'SOL', 'MEM') !PAR
        call this%MpiSol%mpi_init_vg('IMS-X') !PAR
        IF (THIS%ILINMETH.EQ.1) THEN !PAR
          call this%MpiSol%mpi_add_vg('IMS-P') !PAR
          call this%MpiSol%mpi_add_vmt('IMS-P','P','IMSLINEAR', &
            'SOL', 'MM1', 'SOL', 'MEM') !PAR
          call this%MpiSol%mpi_init_vg('IMS-P') !PAR
        ELSEIF (THIS%ILINMETH.EQ.2) THEN !PAR
          call this%MpiSol%mpi_add_vg('IMS-PHAT') !PAR
          call this%MpiSol%mpi_add_vmt('IMS-PHAT', 'PHAT','IMSLINEAR', &
            'SOL', 'MM1', 'SOL', 'MEM') !PAR
          call this%MpiSol%mpi_init_vg('IMS-PHAT') !PAR
          call this%MpiSol%mpi_add_vg('IMS-QHAT') !PAR
          call this%MpiSol%mpi_add_vmt('IMS-QHAT','QHAT','IMSLINEAR', &
            'SOL', 'MM1', 'SOL', 'MEM') !PAR
          call this%MpiSol%mpi_init_vg('IMS-QHAT') !PAR
        END IF !PAR
      end if !PAR
!-------INITIALIZE IMSLINEAR VECTORS
      DO n = 1, iscllen
        THIS%DSCALE(n)  = DONE
        THIS%DSCALE2(n) = DONE
      END DO
      DO n = 1, THIS%NJAPC
        THIS%APC(n)  = DZERO
      END DO
!-------WORKING VECTORS
      DO n = 1, THIS%NEQ
        THIS%ID(n)   = IZERO
        THIS%D(n)    = DZERO
        THIS%P(n)    = DZERO
        THIS%Q(n)    = DZERO
        THIS%Z(n)    = DZERO
      END DO
      DO n = 1, THIS%NIAPC
        THIS%IW(n)   = IZERO
        THIS%W(n)    = DZERO
      END DO
!-------BCGS WORKING VECTORS
      DO n = 1, THIS%NIABCGS
        THIS%T(n)    = DZERO
        THIS%V(n)    = DZERO
        THIS%DHAT(n) = DZERO
        THIS%PHAT(n) = DZERO
        THIS%QHAT(n) = DZERO
      END DO
!-------ILUT AND MILUT WORKING VECTORS      
      DO n = 1, ijlu
        THIS%JLU(n)   = DZERO
      END DO
      DO n = 1, ijw
        THIS%JW(n)   = DZERO
      END DO
      DO n = 1, iwlu
        THIS%WLU(n)  = DZERO
      END DO
!-------REORDERING VECTORS
      DO n = 1, i0 + 1
        THIS%IARO(n) = IZERO
      END DO
      DO n = 1, iolen
        THIS%JARO(n) = IZERO
        THIS%ARO(n)  = DZERO
      END DO
!
!-------REVERSE CUTHILL MCKEE AND MINIMUM DEGREE ORDERING
      IF (THIS%IORD.NE.0) THEN
        CALL IMSLINEARSUB_CALC_ORDER(IOUT, THIS%IPRIMS, THIS%IORD,THIS%NEQ,     &
                                     THIS%NJA,THIS%IA,THIS%JA,                  &
                                     THIS%LORDER,THIS%IORDER)
      END IF
!      
!-------ALLOCATE MEMORY FOR STORING ITERATION CONVERGENCE DATA
      
!
!-------RETURN
      RETURN
    END SUBROUTINE IMSLINEAR_AR

    subroutine imslinear_summary(this, mxiter)
      class(ImsLinearDataType), intent(inout) :: this
      integer(I4B), intent(in) :: mxiter
!     + + + LOCAL VARIABLES + + +
      CHARACTER (LEN= 10) :: clin(0:2)
      CHARACTER (LEN= 31) :: clintit(0:2)
      CHARACTER (LEN= 20) :: cipc(0:4)
      CHARACTER (LEN= 20) :: cscale(0:2)
      CHARACTER (LEN= 25) :: corder(0:2)
      CHARACTER (LEN= 16), DIMENSION(0:4) :: ccnvgopt
      CHARACTER (LEN= 15) :: clevel
      CHARACTER (LEN= 15) :: cdroptol 
      integer(I4B) :: i, j
!     + + + PARAMETERS + + +
!       DATA
      DATA clin  /'UNKNOWN   ', &
                  'CG        ', &
     &            'BCGS      '/
      DATA clintit  /'             UNKNOWN           ', &
                     '       CONJUGATE-GRADIENT      ', &
     &               'BICONJUGATE-GRADIENT STABILIZED'/
      DATA cipc  /'UNKNOWN             ', &
     &            'INCOMPLETE LU       ', &
     &            'MOD. INCOMPLETE LU  ', &
     &            'INCOMPLETE LUT      ', &
     &            'MOD. INCOMPLETE LUT '/
      DATA cscale/'NO SCALING          ', &
     &            'SYMMETRIC SCALING   ', &
     &            'L2 NORM SCALING     '/
      DATA corder/'ORIGINAL ORDERING        ', &
     &            'RCM ORDERING             ', &
     &            'MINIMUM DEGREE ORDERING  '/
      DATA ccnvgopt  /'INFINITY NORM   ', &
     &                'INFINITY NORM S ', &
     &                'L2 NORM         ', &
     &                'RELATIVE L2NORM ', &
                      'L2 NORM W. REL. '/
!       OUTPUT FORMATS
02010 FORMAT (1X,/,7X,'SOLUTION BY THE',1X,A31,1X,'METHOD', &
     &        /,1X,66('-'),/, &
     &        ' MAXIMUM OF ',I0,' CALLS OF SOLUTION ROUTINE',/, &
     &        ' MAXIMUM OF ',I0, &
     &        ' INTERNAL ITERATIONS PER CALL TO SOLUTION ROUTINE',/, &
     &        ' LINEAR ACCELERATION METHOD            =',1X,A,/, &
     &        ' MATRIX PRECONDITIONING TYPE           =',1X,A,/, &
     &        ' MATRIX SCALING APPROACH               =',1X,A,/, &
     &        ' MATRIX REORDERING APPROACH            =',1X,A,/, &
     &        ' NUMBER OF ORTHOGONALIZATIONS          =',1X,I0,/, &
     &        ' HEAD CHANGE CRITERION FOR CLOSURE     =',E15.5,/, &
     &        ' RESIDUAL CHANGE CRITERION FOR CLOSURE =',E15.5,/, &
     &        ' RESIDUAL CONVERGENCE OPTION           =',1X,I0,/, &
     &        ' RESIDUAL CONVERGENCE NORM             =',1X,A,/, &
     &        ' RELAXATION FACTOR                     =',E15.5,/, &
     &        ' RELATIVE TOLERANCE FOR RHO            =',E15.5,/, & !SOL
     &        ' RELATIVE TOLERANCE FOR ALPHA          =',E15.5,/, & !SOL
     &        ' RELATIVE TOLERANCE FOR OMEGA          =',E15.5,/, & !SOL
     &        ' RCLOSE PRE CHECK                      =',I9) !SOL
02015 FORMAT (' NUMBER OF LEVELS                      =',A15,/, &
     &        ' DROP TOLERANCE                        =',A15,//)
2030  FORMAT(1X,A20,1X,6(I6,1X))
2040  FORMAT(1X,20('-'),1X,6(6('-'),1X))
2050  FORMAT(1X,62('-'),/)      !
!------------------------------------------------------------------
      !
      ! -- initialize clevel and cdroptol
      clevel = ''
      cdroptol = ''
      !
      ! -- PRINT MXITER,ITER1,IPC,ISCL,IORD,DVCLOSE,RCLOSE
      write(this%iout,2010)                                         &
                        clintit(THIS%ILINMETH), MXITER, THIS%ITER1, &
                        clin(THIS%ILINMETH), cipc(THIS%IPC),        &
                        cscale(THIS%ISCL), corder(THIS%IORD),       &
                        THIS%NORTH, THIS%DVCLOSE, THIS%RCLOSE,      &
                        THIS%ICNVGOPT, ccnvgopt(THIS%ICNVGOPT),     &
                        THIS%RELAX,                                 &
                        THIS%RHOTOL, THIS%ALPHATOL,THIS%OMEGATOL,   & !SOL
                        THIS%IRCLOSEPRECHK !SOL
      if (this%level > 0) then
        write(clevel, '(i15)') this%level
      end if
      if (this%droptol > DZERO) then
        write(cdroptol, '(e15.5)') this%droptol
      end if
      IF (this%level > 0 .or. this%droptol > DZERO) THEN
        write(this%iout,2015) trim(adjustl(clevel)),               &
                               trim(adjustl(cdroptol))
      ELSE
         write(this%iout,'(//)')
      END IF
      
      if (this%iord /= 0) then
        !                                                                       
        ! -- WRITE SUMMARY OF REORDERING INFORMATION TO LIST FILE                                                  
        if (this%iprims ==  2) then 
          DO i = 1, this%neq, 6 
            write(this%iout,2030) 'ORIGINAL NODE      :',                      &
                              (j,j=i,MIN(i+5,this%neq))                      
            write(this%iout,2040) 
            write(this%iout,2030) 'REORDERED INDEX    :',                      &
                              (this%lorder(j),j=i,MIN(i+5,this%neq))              
            write(this%iout,2030) 'REORDERED NODE     :',                      &
                              (this%iorder(j),j=i,MIN(i+5,this%neq))              
            write(this%iout,2050) 
          END DO 
        END IF 
      end if
      !
      ! -- return
    return
    end subroutine imslinear_summary 
    
    subroutine imslinear_ar_cgc(this, icgc,                             & !CGC
      cgcgneq, cgcgnia, cgcgnja, cgcgia, cgcgja,                        &
      cgclnia, cgclnja, cgclia, cgclja, cgclgja,                        &
      cgcl2gid, cgcg2lid, cgcg2cid, cgcc2gid, cgccizc) !CGC
!     ******************************************************************
!     ALLOCATE STORAGE FOR COARSE GRID CORRECTION ARRAYS
!     ******************************************************************
!
!        SPECIFICATIONS:
!     ------------------------------------------------------------------
      use MemoryManagerModule, only: mem_allocate
      use SimModule, only: ustop, store_error, count_errors
      use MpiExchangeGenModule, only: parallelrun !PAR
      
      implicit none
!     + + + dummy variables + + +
      class(ImsLinearDataType), intent(inout) :: this
      integer(I4B), target                            :: icgc
      !
      integer(I4B), target                            :: cgcgneq
      integer(I4B), target                            :: cgcgnia
      integer(I4B), target                            :: cgcgnja
      integer(I4B), dimension(:), target, contiguous  :: cgcgia
      integer(I4B), dimension(:), target, contiguous  :: cgcgja
      !
      integer(I4B), target                            :: cgclnia
      integer(I4B), target                            :: cgclnja
      integer(I4B), dimension(:), target, contiguous  :: cgclia
      integer(I4B), dimension(:), target, contiguous  :: cgclja
      integer(I4B), dimension(:), target, contiguous  :: cgclgja
      !
      integer(I4B), dimension(:), target, contiguous  :: cgcl2gid
      integer(I4B), dimension(:), target, contiguous  :: cgcg2lid
      integer(I4B), dimension(:), target, contiguous  :: cgcg2cid
      integer(I4B), dimension(:), target, contiguous  :: cgcc2gid
      integer(I4B), dimension(:), target, contiguous  :: cgccizc
!     + + + local variables + + +
      integer(I4B) :: il, ig, j, m, n, nc, nr, ic0, ic1, nlblk, ngblk
      integer(I4B), dimension(:), allocatable :: iwrk
!------------------------------------------------------------------
      ! -- set pointers
      this%icgc     => icgc
      !
      this%cgcgneq  => cgcgneq
      this%cgcgnia  => cgcgnia
      this%cgcgnja  => cgcgnja
      this%cgcgia   => cgcgia
      this%cgcgja   => cgcgja
      !
      this%cgclnia  => cgclnia
      this%cgclnja  => cgclnja
      this%cgclia   => cgclia
      this%cgclja   => cgclja
      this%cgclgja   => cgclgja
      !
      this%cgcl2gid => cgcl2gid
      this%cgcg2lid => cgcg2lid
      this%cgcg2cid => cgcg2cid
      this%cgcc2gid => cgcc2gid
      this%cgccizc  => cgccizc
      !
      call mem_allocate(this%acpc, this%cgcgneq, this%cgcgneq,  &
        'ACPC', trim(this%memoryPath))
      call mem_allocate(this%acpci, this%cgcgnja, 'ACPCI', trim(this%memoryPath))
      call mem_allocate(this%acpcb, 'ACPCB', trim(this%memoryPath))
      call mem_allocate(this%xc, this%cgcgneq, 'XC', trim(this%memoryPath))
      call mem_allocate(this%rhsc, this%cgcgneq, 'RHSC', trim(this%memoryPath))
      call mem_allocate(this%wc, this%cgcgneq, 'WC', trim(this%memoryPath))
      call mem_allocate(this%cgcgiapc, this%cgcgnia, 'CGCGIAPC', trim(this%memoryPath))
      call mem_allocate(this%cgcgjapc, this%cgcgnja, 'CGCGJAPC', trim(this%memoryPath))
      call mem_allocate(this%zc, this%neq, 'ZC', trim(this%memoryPath))
      !
      this%acpcb = 0
      do n = 1, this%neq
        this%zc(n) = DONE
      end do
      !
      ! -- allocate A0ZC
      nr = this%neq
      nc = 0
      nlblk = size(this%cgcl2gid)
      ngblk =  this%cgcgneq
      allocate(iwrk(ngblk))
      iwrk = 0
      do il = 1, nlblk ! my local blocks
        ig = this%cgcl2gid(il)
        ic0 = this%cgcgia(ig)
        ic1 = this%cgcgia(ig+1)-1
        do m = ic0, ic1
          j = this%cgcgja(m)
          iwrk(j) = 1
        end do
      end do ! my local blocks
      nc = sum(iwrk)
      deallocate(iwrk)
      !
      call mem_allocate(this%a0zc, nc, nr, 'A0ZC', trim(this%memoryPath))
      !
      ! -- add variables for MPI point-to-point
      if (parallelrun) then
        call this%MpiSol%mpi_add_vg('IMS-ZC') !PAR
        call this%MpiSol%mpi_add_vmt('IMS-ZC','ZC','IMSLINEAR', &
          'SOL', 'MM1', 'SOL', 'MEM') !PAR
        call this%MpiSol%mpi_init_vg('IMS-ZC') !PAR
      end if
      !
!-------return
      return
    end subroutine imslinear_ar_cgc
      
    subroutine allocate_scalars(this)
      use MemoryManagerModule, only: mem_allocate
      class(ImsLinearDataType), intent(inout) :: this
      !
      ! -- scalars
      call mem_allocate(this%iout, 'IOUT', this%memoryPath)
      call mem_allocate(this%ilinmeth, 'ILINMETH', this%memoryPath)
      call mem_allocate(this%iter1, 'ITER1', this%memoryPath)
      call mem_allocate(this%ipc, 'IPC', this%memoryPath)
      call mem_allocate(this%iscl, 'ISCL', this%memoryPath)
      call mem_allocate(this%iord, 'IORD', this%memoryPath)
      call mem_allocate(this%north, 'NORTH', this%memoryPath)
      call mem_allocate(this%icnvgopt, 'ICNVGOPT', this%memoryPath)
      call mem_allocate(this%iacpc, 'IACPC', this%memoryPath)
      call mem_allocate(this%niterc, 'NITERC', this%memoryPath)
      call mem_allocate(this%niabcgs, 'NIABCGS', this%memoryPath)
      call mem_allocate(this%niapc, 'NIAPC', this%memoryPath)
      call mem_allocate(this%njapc, 'NJAPC', this%memoryPath)
      call mem_allocate(this%dvclose, 'DVCLOSE', this%memoryPath)
      call mem_allocate(this%rclose, 'RCLOSE', this%memoryPath)
      call mem_allocate(this%relax, 'RELAX', this%memoryPath)
      call mem_allocate(this%epfact, 'EPFACT', this%memoryPath)
      call mem_allocate(this%l2norm0, 'L2NORM0', this%memoryPath)
      call mem_allocate(this%droptol, 'DROPTOL', this%memoryPath)
      call mem_allocate(this%level, 'LEVEL', this%memoryPath)
      call mem_allocate(this%njlu, 'NJLU', this%memoryPath)
      call mem_allocate(this%njw, 'NJW', this%memoryPath)
      call mem_allocate(this%nwlu, 'NWLU', this%memoryPath)
      call mem_allocate(this%rhotol, 'RHOTOL', this%memoryPath) !SOL
      call mem_allocate(this%alphatol, 'ALPHATOL', this%memoryPath) !SOL
      call mem_allocate(this%omegatol, 'OMEGATOL', this%memoryPath) !SOL
      call mem_allocate(this%ircloseprechk, 'IRCLOSEPRECHK', this%memoryPath) !SOL
      call mem_allocate(this%icgc, 'ICGC', this%memoryPath) !CGC
      call mem_allocate(this%cgcsol, 'CGCSOL', this%memoryPath) !CGC
      !
      ! -- initialize
      this%iout = 0
      this%ilinmeth = 0
      this%iter1 = 0
      this%ipc = 0
      this%iscl = 0
      this%iord = 0
      this%north = 0
      this%icnvgopt = 0
      this%iacpc = 0
      this%niterc = 0
      this%niabcgs = 0
      this%niapc = 0
      this%njapc = 0
      this%dvclose = DZERO
      this%rclose = DZERO
      this%relax = DZERO
      this%epfact = DZERO
      this%l2norm0 = 0
      this%droptol = DZERO
      this%level = 0
      this%njlu = 0
      this%njw = 0
      this%nwlu = 0
      this%rhotol = DZERO !SOL
      this%alphatol = DZERO !SOL
      this%omegatol = DZERO !SOL
      this%ircloseprechk = 0 !SOL
      this%icgc = 0 !CGC
      this%cgcsol = 0 !CGC
      !
      ! --Return
      return
    end subroutine allocate_scalars
    
    subroutine IMSLINEAR_DA(this)
      use MemoryManagerModule, only: mem_deallocate
      class(ImsLinearDataType), intent(inout) :: this
      !
      ! -- arrays
      call mem_deallocate(this%dscale)
      call mem_deallocate(this%dscale2)
      call mem_deallocate(this%iapc)
      call mem_deallocate(this%japc)
      call mem_deallocate(this%apc)
      call mem_deallocate(this%iw)
      call mem_deallocate(this%w)
      call mem_deallocate(this%jlu)
      call mem_deallocate(this%jw)
      call mem_deallocate(this%wlu)
      call mem_deallocate(this%lorder)
      call mem_deallocate(this%iorder)
      call mem_deallocate(this%iaro)
      call mem_deallocate(this%jaro)
      call mem_deallocate(this%aro)
      call mem_deallocate(this%id)
      call mem_deallocate(this%d)
      call mem_deallocate(this%p)
      call mem_deallocate(this%q)
      call mem_deallocate(this%z)
      call mem_deallocate(this%t)
      call mem_deallocate(this%v)
      call mem_deallocate(this%dhat)
      call mem_deallocate(this%phat)
      call mem_deallocate(this%qhat)
      !
      ! -- scalars
      call mem_deallocate(this%iout)
      call mem_deallocate(this%ilinmeth)
      call mem_deallocate(this%iter1)
      call mem_deallocate(this%ipc)
      call mem_deallocate(this%iscl)
      call mem_deallocate(this%iord)
      call mem_deallocate(this%north)
      call mem_deallocate(this%icnvgopt)
      call mem_deallocate(this%iacpc)
      call mem_deallocate(this%niterc)
      call mem_deallocate(this%niabcgs)
      call mem_deallocate(this%niapc)
      call mem_deallocate(this%njapc)
      call mem_deallocate(this%dvclose)
      call mem_deallocate(this%rclose)
      call mem_deallocate(this%relax)
      call mem_deallocate(this%epfact)
      call mem_deallocate(this%l2norm0)
      call mem_deallocate(this%droptol)
      call mem_deallocate(this%level)
      call mem_deallocate(this%njlu)
      call mem_deallocate(this%njw)
      call mem_deallocate(this%nwlu)
      call mem_deallocate(this%rhotol) !SOL
      call mem_deallocate(this%alphatol) !SOL
      call mem_deallocate(this%omegatol) !SOL
      call mem_deallocate(this%ircloseprechk) !SOL
      call mem_deallocate(this%icgc) !CGC
      call mem_deallocate(this%cgcsol) !CGC
      !
      ! -- nullify pointers
      nullify(this%iprims)
      nullify(this%neq)
      nullify(this%nja)
      nullify(this%ia)
      nullify(this%ja)
      nullify(this%amat)
      nullify(this%rhs)
      nullify(this%x)
      !
      ! --Return
      return
    end subroutine IMSLINEAR_DA
    
    SUBROUTINE SET_IMSLINEAR_INPUT(THIS, IFDPARAM)
      IMPLICIT NONE
!     + + + DUMMY ARGUMENTS + + +
      CLASS(ImsLinearDataType), INTENT(INOUT) :: THIS
      integer(I4B), INTENT(IN) :: IFDPARAM
!     + + + LOCAL DEFINITIONS + + +
!     + + + PARAMETERS + + +
!     + + + FUNCTIONS + + +
!
!     + + + CODE + + +
      SELECT CASE ( IFDPARAM )
        ! Simple option
        CASE(1)
          THIS%ITER1 = 50
          THIS%ILINMETH=1
          THIS%IPC = 1
          THIS%ISCL = 0
          THIS%IORD = 0
          THIS%DVCLOSE = DEM3
          THIS%RCLOSE = DEM1
          THIS%RELAX = DZERO
          THIS%LEVEL = 0
          THIS%DROPTOL = DZERO
          THIS%NORTH = 0
          THIS%RHOTOL = DSAME !SOL
          THIS%ALPHATOL = DSAME !SOL
          THIS%OMEGATOL = DSAME !SOL
          THIS%IRCLOSEPRECHK = 0 !SOL
        ! Moderate
        CASE(2)
          THIS%ITER1 = 100
          THIS%ILINMETH=2
          THIS%IPC = 2
          THIS%ISCL = 0
          THIS%IORD = 0
          THIS%DVCLOSE = DEM2
          THIS%RCLOSE = DEM1
          THIS%RELAX = 0.97D0
          THIS%LEVEL = 0
          THIS%DROPTOL = DZERO
          THIS%NORTH = 0
          THIS%RHOTOL = DSAME !SOL
          THIS%ALPHATOL = DSAME !SOL
          THIS%OMEGATOL = DSAME !SOL
          THIS%IRCLOSEPRECHK = 0 !SOL
        ! Complex
        CASE(3)
          THIS%ITER1 = 500
          THIS%ILINMETH=2
          THIS%IPC = 3
          THIS%ISCL = 0
          THIS%IORD = 0
          THIS%DVCLOSE = DEM1
          THIS%RCLOSE = DEM1
          THIS%RELAX = DZERO
          THIS%LEVEL = 5
          THIS%DROPTOL = DEM4
          THIS%NORTH = 2
          THIS%RHOTOL = DSAME !SOL
          THIS%ALPHATOL = DSAME !SOL
          THIS%OMEGATOL = DSAME !SOL
          THIS%IRCLOSEPRECHK = 0 !SOL
      END SELECT
      RETURN
    END SUBROUTINE SET_IMSLINEAR_INPUT

      SUBROUTINE IMSLINEAR_AP(THIS,ICNVG,KSTP,KITER,IN_ITER,                  &
                              NCONV, CONVNMOD, CONVMODSTART, LOCDV, LOCDR,    &
                              CACCEL, ITINNER, CONVLOCDV, CONVLOCDR,          &
                              DVMAX, DRMAX, CONVDVMAX, CONVDRMAX,             &
                              ICNVGPREV) !SOL
!
!     ******************************************************************
!     SOLUTION BY THE CONJUGATE GRADIENT METHOD -
!                                          UP TO ITER1 ITERATIONS
!     ******************************************************************
!
!        SPECIFICATIONS:
!     ------------------------------------------------------------------
      USE SimModule
      IMPLICIT NONE
!     + + + DUMMY ARGUMENTS + + +
      CLASS(ImsLinearDataType), INTENT(INOUT) :: THIS
      integer(I4B), INTENT(INOUT)                          :: ICNVG
      integer(I4B), INTENT(INOUT)                          :: ICNVGPREV !SOL
      integer(I4B), INTENT(IN)                             :: KSTP
      integer(I4B), INTENT(IN)                             :: KITER
      integer(I4B), INTENT(INOUT)                          :: IN_ITER
      ! CONVERGENCE INFORMATION
      integer(I4B), INTENT(IN) :: NCONV
      integer(I4B), INTENT(IN) :: CONVNMOD
      integer(I4B), DIMENSION(CONVNMOD+1), INTENT(INOUT) ::CONVMODSTART
      integer(I4B), DIMENSION(CONVNMOD), INTENT(INOUT) :: LOCDV
      integer(I4B), DIMENSION(CONVNMOD), INTENT(INOUT) :: LOCDR
      character(len=31), DIMENSION(NCONV), INTENT(INOUT) :: CACCEL
      integer(I4B), DIMENSION(NCONV), INTENT(INOUT) :: ITINNER
      integer(I4B), DIMENSION(CONVNMOD, NCONV), INTENT(INOUT) :: CONVLOCDV
      integer(I4B), DIMENSION(CONVNMOD, NCONV), INTENT(INOUT) :: CONVLOCDR
      real(DP), DIMENSION(CONVNMOD), INTENT(INOUT) :: DVMAX
      real(DP), DIMENSION(CONVNMOD), INTENT(INOUT) :: DRMAX
      real(DP), DIMENSION(CONVNMOD, NCONV), INTENT(INOUT) :: CONVDVMAX
      real(DP), DIMENSION(CONVNMOD, NCONV), INTENT(INOUT) :: CONVDRMAX
      
!     + + + LOCAL DEFINITIONS + + +
      integer(I4B) :: n
      integer(I4B) :: innerit
      integer(I4B) :: irc
      integer(I4B) :: itmax
      real(DP) :: tv
      real(DP) :: rmax
      logical :: lsolve !SOL
      real(DP) :: rcnvg !SOL
!     + + + PARAMETERS + + +
!     + + + FUNCTIONS + + +
!
!     + + + CODE + + +
!
!-------SET EPFACT BASED ON MFUSG TIMESTEP
      IF (THIS%ICNVGOPT ==  2) THEN
        IF (KSTP ==  1) THEN
          THIS%EPFACT = 0.01
        ELSE
          THIS%EPFACT = 0.10
        END IF
      ELSE IF (THIS%ICNVGOPT ==  4) THEN
        THIS%EPFACT = DEM4
      ELSE
        THIS%EPFACT = DONE
      END IF

!-------SCALE PROBLEM
      IF (THIS%ISCL.NE.0) THEN
        CALL IMSLINEARSUB_SCALE(0,THIS%ISCL,                                    &
                         THIS%NEQ,THIS%NJA,THIS%IA,THIS%JA,                     &
                         THIS%AMAT,THIS%X,THIS%RHS,                             &
                         THIS%DSCALE,THIS%DSCALE2)
      END IF
!
!-------PERMUTE ROWS, COLUMNS, AND RHS
      IF (THIS%IORD.NE.0) THEN
        CALL ims_dperm(THIS%NEQ, THIS%NJA, THIS%AMAT,THIS%JA,THIS%IA, &
     &                 THIS%ARO,THIS%JARO,THIS%IARO,THIS%LORDER,THIS%ID,1)
        CALL ims_vperm(THIS%NEQ, THIS%X, THIS%LORDER)
        CALL ims_vperm(THIS%NEQ, THIS%RHS, THIS%LORDER)
        THIS%IA0 => THIS%IARO
        THIS%JA0 => THIS%JARO
        THIS%A0  => THIS%ARO
      ELSE
        THIS%IA0 => THIS%IA
        THIS%JA0 => THIS%JA
        THIS%A0  => THIS%AMAT
      END IF
!
!-------UPDATE PRECONDITIONER
!      CALL IMSLINEARSUB_PCU(this%iout,THIS%NJA,THIS%NEQ,THIS%NIAPC,THIS%NJAPC,  & !SOL
!                            THIS%IPC, THIS%RELAX, THIS%A0, THIS%IA0, THIS%JA0,  & !SOL
!                            THIS%APC,THIS%IAPC,THIS%JAPC,THIS%IW,THIS%W,        & !SOL
!                            THIS%LEVEL, THIS%DROPTOL, THIS%NJLU, THIS%NJW,      & !SOL
!                            THIS%NWLU, THIS%JLU, THIS%JW, THIS%WLU) !SOL
!-------INITIALIZE SOLUTION VARIABLE AND ARRAYS
      IF (KITER ==  1 ) THIS%NITERC = 0
      irc    = 1
      ICNVG  = 0
      DO n = 1, THIS%NEQ
        THIS%D(n) = DZERO
        THIS%P(n) = DZERO
        THIS%Q(n) = DZERO
        THIS%Z(n) = DZERO
      END DO
!-------CALCULATE INITIAL RESIDUAL
      CALL IMSLINEARSUB_MV(THIS%NJA,THIS%NEQ,THIS%A0,THIS%X,THIS%D,             &
                           THIS%IA0,THIS%JA0)
      ! -- MPI parallel: point-to-point comm. of THIS%X and update of THIS%D
      call this%MpiSol%mpi_local_exchange(this%memoryPath, 'IMS-X', .false.) !PAR
      call this%MpiSol%mpi_mv_halo(this%memoryPath, 'IMS-X', this%d) !PAR
      rmax = DZERO
      THIS%L2NORM0 = DZERO
      DO n = 1, THIS%NEQ
        tv   = THIS%D(n)
        THIS%D(n) = THIS%RHS(n) - tv
        IF (ABS( THIS%D(n) ) > rmax ) rmax = ABS( THIS%D(n) )
        THIS%L2NORM0 = THIS%L2NORM0 + THIS%D(n) * THIS%D(n)
      END DO
      ! -- MPI parallel: collective sum of THIS%L2NORM0
      call this%MpiSol%mpi_global_exchange_all_sum(THIS%L2NORM0) !PAR
      ! -- MPI parallel: collective sum of rmax
      call this%MpiSol%mpi_global_exchange_all_absmax(rmax) !PAR
      THIS%L2NORM0 = SQRT(THIS%L2NORM0)
!-------CHECK FOR EXACT SOLUTION
      itmax = THIS%ITER1
      IF (rmax ==  DZERO) THEN
        itmax = 0
        ICNVG = 1
      END IF
!-------RCLOSE PRE CHECK BEFORE SOLVING
      LSOLVE = .TRUE. !SOL
      IF ((ICNVGPREV /= 0) .AND. (THIS%IRCLOSEPRECHK == 1)) THEN !SOL
        IF (THIS%ICNVGOPT.EQ.2 .OR. THIS%ICNVGOPT.EQ.3 .OR.                     & !SOL
            THIS%ICNVGOPT.EQ.4) THEN !SOL
          RCNVG = THIS%L2NORM0 !SOL
        ELSE !SOL
          RCNVG = RMAX !SOL
        END IF !SOL
        CALL IMSLINEARSUB_TESTCNVG(THIS%ICNVGOPT, ICNVG, 1,                     & !SOL
                                   DZERO, RCNVG,                                & !SOL
                                   THIS%L2NORM0,                                & !SOL
                                   THIS%EPFACT, THIS%DVCLOSE, THIS%RCLOSE) !SOL
        IF (ICNVG /= 0) THEN !SOL
          LSOLVE = .FALSE. !SOL
          INNERIT = 0 !SOL
        END IF !SOL
      END IF !SOL
! 
      IF (LSOLVE) THEN !SOL
!-------UPDATE PRECONDITIONER
      CALL IMSLINEARSUB_PCU(this%iout,THIS%NJA,THIS%NEQ,THIS%NIAPC,THIS%NJAPC,  & !SOL
                            THIS%IPC, THIS%RELAX, THIS%A0, THIS%IA0, THIS%JA0,  & !SOL
                            THIS%APC,THIS%IAPC,THIS%JAPC,THIS%IW,THIS%W,        & !SOL
                            THIS%LEVEL, THIS%DROPTOL, THIS%NJLU, THIS%NJW,      & !SOL
                            THIS%NWLU, THIS%JLU, THIS%JW, THIS%WLU,             & !SOL
                            THIS%IBJFLAG, CONVNMOD, CONVMODSTART) !BJ
!-------UPDATE COARSE GRID PRECONDITIONER
      if (this%icgc == 1) then !CGC
        call this%imslinearsub_pcc(convnmod, convmodstart) !CGC
      end if !CGC
!      
!-------SOLUTION BY THE CONJUGATE GRADIENT METHOD
      IF (THIS%ILINMETH ==  1) THEN
        CALL THIS%IMSLINEARSUB_CG(ICNVG, itmax, innerit,                        &
                             NCONV, CONVNMOD, CONVMODSTART, LOCDV, LOCDR,       &
                             CACCEL, ITINNER, CONVLOCDV, CONVLOCDR,             &
                             DVMAX, DRMAX, CONVDVMAX, CONVDRMAX)
!-------SOLUTION BY THE BICONJUGATE GRADIENT STABILIZED METHOD
      ELSE IF (THIS%ILINMETH ==  2) THEN
        CALL THIS%IMSLINEARSUB_BCGS(ICNVG, itmax, innerit,                      &
                               NCONV, CONVNMOD, CONVMODSTART, LOCDV, LOCDR,     &
                               CACCEL, ITINNER, CONVLOCDV, CONVLOCDR,           &
                               DVMAX, DRMAX, CONVDVMAX, CONVDRMAX)
      END IF
      END IF !SOL
!
!-------BACK PERMUTE AMAT, SOLUTION, AND RHS
      IF (THIS%IORD.NE.0) THEN
        CALL ims_dperm(THIS%NEQ, THIS%NJA, THIS%A0, THIS%JA0, THIS%IA0,         &
     &                 THIS%AMAT, THIS%JA,THIS%IA,THIS%IORDER,THIS%ID,1)
        CALL ims_vperm(THIS%NEQ, THIS%X, THIS%IORDER)
        CALL ims_vperm(THIS%NEQ, THIS%RHS, THIS%IORDER)
      END IF
!
!-------UNSCALE PROBLEM
      IF (THIS%ISCL.NE.0) THEN
        CALL IMSLINEARSUB_SCALE(1, THIS%ISCL,                                   &
                                THIS%NEQ, THIS%NJA, THIS%IA, THIS%JA,           &
                                THIS%AMAT, THIS%X, THIS%RHS,                    &
                                THIS%DSCALE, THIS%DSCALE2)
      END IF
!
!-------SET IMS INNER ITERATION NUMBER (IN_ITER) TO NUMBER OF
!       IMSLINEAR INNER ITERATIONS (innerit)
      IN_ITER = innerit
!
!     -- MPI parallel: update the solution for the halo models
      call this%MpiSol%mpi_local_exchange(this%memoryPath, 'IMS-X', .false.) !PAR
!
!-------RETURN
      RETURN
!
      END SUBROUTINE IMSLINEAR_AP

    
! -- IMSLinearModule subroutines that do not depend on data stored in the IMSLinearModule class
!    all data is passed through subroutine calls
!                                                                       
!-------ROUTINE TO CALCULATE LORDER AND IORDER FOR REORDERING           
      SUBROUTINE IMSLINEARSUB_CALC_ORDER(IOUT, IPRIMS, IORD, NEQ, NJA, IA, JA,  &
     &                                   LORDER, IORDER)                     
      use SimModule, only: ustop, store_error, count_errors
      IMPLICIT NONE 
!       + + + DUMMY ARGUMENTS + + +                                       
        integer(I4B), INTENT(IN) :: IOUT 
        integer(I4B), INTENT(IN) :: IPRIMS 
        integer(I4B), INTENT(IN) :: IORD 
        integer(I4B), INTENT(IN) :: NEQ 
        integer(I4B), INTENT(IN) :: NJA 
        integer(I4B), DIMENSION(NEQ+1), INTENT(IN)  :: IA 
        integer(I4B), DIMENSION(NJA),   INTENT(IN)  :: JA 
        integer(I4B), DIMENSION(NEQ), INTENT(INOUT) :: LORDER 
        integer(I4B), DIMENSION(NEQ), INTENT(INOUT) :: IORDER 
!       + + + LOCAL DEFINITIONS + + +                                     
        character (len=LINELENGTH) :: errmsg 
        integer(I4B) :: n 
        integer(I4B) :: nsp 
        integer(I4B), DIMENSION(:), ALLOCATABLE :: iwork0, iwork1 
        integer(I4B) :: iflag 
!       + + + PARAMETERS + + +                                            
!       + + + FUNCTIONS + + +                                             
!       + + + FORMATS + + +                                               
!       + + + CODE + + +                                                  
        DO n = 1, NEQ 
          LORDER(n) = IZERO 
          IORDER(n) = IZERO 
        END DO 
        ALLOCATE ( iwork0(NEQ)  ) 
        SELECT CASE ( IORD ) 
          CASE ( 1 ) 
            ALLOCATE ( iwork1(NEQ) ) 
            CALL ims_genrcm(NEQ, NJA, IA, JA,                           &
     &                      LORDER, iwork0, iwork1 )                        
          CASE ( 2 ) 
            nsp = 3 * NEQ + 4 * NJA 
            ALLOCATE ( iwork1(nsp)  ) 
            CALL ims_odrv(NEQ, NJA, nsp, IA, JA, LORDER, iwork0,                 &
                          iwork1, iflag)                           
            IF (iflag.NE.0) THEN 
              write(errmsg,'(A,1X,A)')                                           &
                'IMSLINEARSUB_CALC_ORDER ERROR CREATING MINIMUM DEGREE ',        &
                'ORDER PERMUTATION '                           
              call store_error(errmsg) 
              !call ustop()                                             
            END IF 
        END SELECT 
!                                                                       
!         GENERATE INVERSE OF LORDER                                    
        DO n = 1, NEQ 
          IORDER( LORDER(n) ) = n 
        END DO 
!
!         DEALLOCATE TEMPORARY STORAGE                                  
        DEALLOCATE ( iwork0, iwork1 ) 
!
        if (count_errors() > 0) then
          call parser%StoreErrorUnit()
          call ustop()
        endif
!
!---------RETURN                                                        
        RETURN 
      END SUBROUTINE IMSLINEARSUB_CALC_ORDER 
!                                                                       
!-------ROUTINE TO SCALE THE COEFFICIENT MATRIX (AMAT),                 
!       THE RHS (B), AND THE ESTIMATE OF X (X)                          
      SUBROUTINE IMSLINEARSUB_SCALE(IOPT, ISCL, NEQ, NJA, IA, JA, AMAT, X, B,   &
     &                              DSCALE, DSCALE2)                            
        IMPLICIT NONE 
!       + + + DUMMY ARGUMENTS + + +                                       
        integer(I4B), INTENT(IN) :: IOPT 
        integer(I4B), INTENT(IN) :: ISCL 
        integer(I4B), INTENT(IN) :: NEQ 
        integer(I4B), INTENT(IN) :: NJA 
        integer(I4B), DIMENSION(NEQ+1), INTENT(IN) :: IA 
        integer(I4B), DIMENSION(NJA),   INTENT(IN) :: JA 
        real(DP), DIMENSION(NJA),  INTENT(INOUT) :: AMAT 
        real(DP), DIMENSION(NEQ),  INTENT(INOUT) :: X 
        real(DP), DIMENSION(NEQ),  INTENT(INOUT) :: B 
        real(DP), DIMENSION(NEQ),  INTENT(INOUT) :: DSCALE 
        real(DP), DIMENSION(NEQ), INTENT(INOUT)  :: DSCALE2 
!       + + + LOCAL DEFINITIONS + + +                                     
        integer(I4B) :: i, n 
        integer(I4B) :: id, jc 
        integer(I4B) :: i0, i1 
        real(DP) :: v, c1, c2 
!       + + + FUNCTIONS + + +                                             
!       + + + CODE + + +                                                  
!                                                                       
!---------SCALE SCALE AMAT, X, AND B                                    
        IF (IOPT ==  0) THEN 
!-----------SYMMETRIC SCALING                                           
          SELECT CASE ( ISCL ) 
            CASE ( 1 ) 
              DO n = 1, NEQ 
                id   = IA(n) 
                v    = AMAT(id) 
                c1   = DONE / SQRT( ABS( v ) ) 
                DSCALE(n)  = c1 
                DSCALE2(n) = c1 
              END DO 
!               SCALE AMAT -- AMAT = DSCALE(row) * AMAT(i) * DSCALE2(col)
              DO n = 1, NEQ 
                c1 = DSCALE(n) 
                i0 = IA(n) 
                i1 = IA(n+1) - 1 
                DO i = i0, i1 
                  jc = JA(i) 
                  c2 = DSCALE2(jc) 
                  AMAT(i) = c1 * AMAT(i) * c2 
                END DO 
              END DO 
!-----------L-2 NORM SCALING                                            
            CASE ( 2 ) 
!               SCALE EACH ROW SO THAT THE L-2 NORM IS 1                
              DO n = 1, NEQ 
                c1 = DZERO 
                i0 = IA(n) 
                i1 = IA(n+1) - 1 
                DO i = i0, i1 
                  c1 = c1 + AMAT(i) * AMAT(i) 
                END DO 
                c1 = SQRT( c1 ) 
                IF (c1 ==  DZERO) THEN 
                  c1 = DONE 
                ELSE 
                  c1 = DONE / c1 
                END IF 
                DSCALE(n) = c1 
!                 INITIAL SCALING OF AMAT -- AMAT = DSCALE(row) * AMAT(i)
                DO i = i0, i1 
                  AMAT(i) = c1 * AMAT(i) 
                END DO 
              END DO 
!               SCALE EACH COLUMN SO THAT THE L-2 NORM IS 1             
              DO n = 1, NEQ 
                DSCALE2(n) = DZERO 
              END DO 
              c2 = DZERO 
              DO n = 1, NEQ 
                i0 = IA(n) 
                i1 = IA(n+1) - 1 
                DO i = i0, i1 
                  jc = JA(i) 
                  c2 = AMAT(i) 
                  DSCALE2(jc) = DSCALE2(jc) + c2 * c2 
                END DO 
              END DO 
              DO n = 1, NEQ 
                c2 = DSCALE2(n) 
                IF (c2 ==  DZERO) THEN 
                  c2 = DONE 
                ELSE 
                  c2 = DONE / SQRT( c2 ) 
                END IF 
                DSCALE2(n) = c2 
              END DO 
!               FINAL SCALING OF AMAT -- AMAT = DSCALE2(col) * AMAT(i)  
              DO n = 1, NEQ 
                i0 = IA(n) 
                i1 = IA(n+1) - 1 
                DO i = i0, i1 
                  jc = JA(i) 
                  c2 = DSCALE2(jc) 
                  AMAT(i) = c2 * AMAT(i) 
                END DO 
              END DO 
          END SELECT 
!-----------SCALE X AND B                                               
          DO n = 1, NEQ 
            c1    = DSCALE(n) 
            c2    = DSCALE2(n) 
            X(n)  = X(n) / c2 
            B(n)  = B(n) * c1 
          END DO 
!---------UNSCALE SCALE AMAT, X, AND B                                  
        ELSE 
          DO n = 1, NEQ 
            c1 = DSCALE(n) 
            i0 = IA(n) 
            i1 = IA(n+1) - 1 
!             UNSCALE AMAT                                              
            DO i = i0, i1 
              jc = JA(i) 
              c2 = DSCALE2(jc) 
              AMAT(i) = ( DONE / c1 ) * AMAT(i) * ( DONE / c2 ) 
            END DO 
!             UNSCALE X AND B                                           
            c2   = DSCALE2(n) 
            X(n) = X(n) * c2 
            B(n) = B(n) / c1 
          END DO 
        END IF 
!---------RETURN                                                        
        RETURN 
      END SUBROUTINE IMSLINEARSUB_SCALE 
!                                                                       
!-------ROUTINE TO UPDATE THE PRECONDITIONER                            
      SUBROUTINE IMSLINEARSUB_PCU(IOUT, NJA, NEQ, NIAPC, NJAPC, IPC, RELAX,     &
                                  AMAT, IA, JA, APC, IAPC, JAPC, IW, W,         &
                                  LEVEL, DROPTOL, NJLU, NJW, NWLU, JLU, JW, WLU,&
                                  IBJFLAG, CONVNMOD, CONVMODSTART) !BJ
      use SimModule, only: ustop, store_error, count_errors
!       + + + DUMMY ARGUMENTS + + +                                       
        integer(I4B), INTENT(IN) :: IOUT 
        integer(I4B), INTENT(IN) :: NJA 
        integer(I4B), INTENT(IN) :: NEQ 
        integer(I4B), INTENT(IN) :: NIAPC 
        integer(I4B), INTENT(IN) :: NJAPC 
        integer(I4B), INTENT(IN) :: IPC 
        real(DP), INTENT(IN) :: RELAX 
        real(DP), DIMENSION(NJA),  INTENT(IN)     :: AMAT 
        integer(I4B), DIMENSION(NEQ+1), INTENT(IN)    :: IA 
        integer(I4B), DIMENSION(NJA), INTENT(IN)      :: JA 
        real(DP), DIMENSION(NJAPC), INTENT(INOUT) :: APC 
        integer(I4B), DIMENSION(NIAPC+1), INTENT(INOUT) :: IAPC 
        integer(I4B), DIMENSION(NJAPC), INTENT(INOUT)   :: JAPC 
        integer(I4B), DIMENSION(NIAPC), INTENT(INOUT)   :: IW 
        real(DP), DIMENSION(NIAPC), INTENT(INOUT) :: W 
        ! ILUT
        integer(I4B), INTENT(IN) :: LEVEL
        real(DP), INTENT(IN) :: DROPTOL
        integer(I4B), INTENT(IN) :: NJLU
        integer(I4B), INTENT(IN) :: NJW
        integer(I4B), INTENT(IN) :: NWLU
        integer(I4B), DIMENSION(NJLU), INTENT(INOUT) :: JLU
        integer(I4B), DIMENSION(NJW),  INTENT(INOUT) :: JW
        real(DP), DIMENSION(NWLU),  INTENT(INOUT) :: WLU
        integer(I4B), INTENT(IN) :: IBJFLAG !BJ
        integer(I4B), INTENT(IN) :: CONVNMOD !BJ
        integer(I4B), DIMENSION(CONVNMOD+1), INTENT(INOUT) ::CONVMODSTART !BJ
!       + + + LOCAL DEFINITIONS + + +                                     
        character(len=LINELENGTH) :: errmsg
        character(len=80), dimension(3) :: cerr
        integer(I4B) :: izero 
        integer(I4B) :: ierr
        real(DP) :: delta 
!       + + + FUNCTIONS + + +  
!       + + + DATA + + +
        DATA cerr  /'INCOMPREHENSIBLE ERROR - MATRIX MUST BE WRONG.              ', &
                    'INSUFFICIENT STORAGE IN ARRAYS ALU, JLU TO STORE FACTORS.   ', &
                    'ZERO ROW ENCOUNTERED.                                       '/
        
!       + + + FORMATS + + +                                               
 2000   FORMAT (/,' MATRIX IS SEVERELY NON-DIAGONALLY DOMINANT.',               &
     &          /,' ADDING SMALL VALUE TO PIVOT (IMSLINEARSUB_PCU)')           
!       + + + CODE + + +                                                  
        izero = 0 
        delta = DZERO 
        PCSCALE: DO
          SELECT CASE(IPC) 
!             ILU0 AND MILU0                                              
            CASE (1,2) 
              CALL IMSLINEARSUB_PCILU0(NJA, NEQ, AMAT, IA, JA,                   &
                                       APC, IAPC, JAPC, IW, W,                   &
                                       RELAX, izero, delta,                      &
                                       IBJFLAG, CONVNMOD, CONVMODSTART) !BJ
!             ILUT AND MILUT                                             
            CASE (3,4) 
              ierr = 0
              CALL IMSLINEARSUB_PCMILUT(NEQ, AMAT, JA, IA,                      &
                                        LEVEL, DROPTOL, RELAX,                  &
                                        APC, JLU, IW, NJAPC, WLU, JW, ierr,     &
                                        izero, delta)
              IF (ierr.NE.0) THEN
                write(errmsg,'(a,1x,a)') 'ILUT ERROR: ', cerr(-ierr)
                call store_error(errmsg)
                call parser%StoreErrorUnit()
                call ustop()
              END IF
!           ADDITIONAL PRECONDITIONERS                     
            CASE DEFAULT
              izero = 0
          END SELECT
          IF (izero < 1) THEN 
            EXIT PCSCALE 
          END IF 
          delta = 1.5D0 * delta + 0.001 
          izero = 0 
          IF (delta > DHALF) THEN 
            WRITE(IOUT,2000) 
            delta = DHALF 
            izero = 2 
          END IF
        END DO PCSCALE
!---------RETURN                                                        
        RETURN 
      END SUBROUTINE IMSLINEARSUB_PCU 

!-------COARSE GRID PRECONDITIONER
      subroutine imslinearsub_pcc(this, nlblk, lblkeq)
        use MpiExchangeGenModule, only: parallelrun !PAR
        use SimModule, only: ustop
!       + + + dummy arguments + + +                                       
        class(ImsLinearDataType), intent(inout) :: this
        integer(I4B), intent(in) :: nlblk
        integer(I4B), dimension(:), intent(in) :: lblkeq
!       + + + local definitions + + +                                     
        integer(I4B) :: i, j, k, ngblk, lna1d
        integer(I4B) :: icbl, icbg
        integer(I4B) :: n, m, icol, ic, ic0, ic1, jc, jc0, jc1
        integer(I4B) ::  ireq0, ireq1, iceq0, iceq1, ieq
        integer(I4B) :: p, q
        integer(I4B) :: izero
        integer(I4B), dimension(2) :: eqbnd
        real(DP) :: delta, relax
        real(DP), dimension(:), allocatable :: la1d
        real(DP) :: tv
        real(DP), dimension(:), allocatable :: dwrk
        
!       + + + code + + +                                                  
        !
        !-- correct ZC for constant head cells
        do i = 1, this%neq
          ic0 = this%ia(i)
          ic1 = this%ia(i+1)-1
          tv = DZERO
          do m = ic0+1, ic1
            j = this%ja(m)
            tv = tv + abs(this%a0(m))
          end do
          if (tv == DZERO) then
            j = this%ja(ic0)
            this%zc(j) = DZERO
          end if
        end do
        !
        ngblk = this%cgcgneq
        lna1d = 0
        do i = 1, nlblk ! my local blocks
          j = this%cgcl2gid(i)
          j = this%cgcg2cid(j)
          jc0 = this%cgcgia(j)
          jc1 = this%cgcgia(j+1)-1
          lna1d = lna1d + (jc1 - jc0 + 1)
        end do ! my local blocks
        !
        allocate(la1d(lna1d))
        !
        ! initalize and fill local A*ZC
        do i = 1, size(this%a0zc,2)
          do j = 1, size(this%a0zc,1)
            this%a0zc(j,i) = DZERO
          end do  
        end do
        !
        ! -- allocate work array for interface contributions
        if (parallelrun) then
          allocate(dwrk(this%neq))
        end if
        !
        ! -- point-to-point communication of ZC
        call this%MpiSol%mpi_local_exchange(this%memoryPath, 'IMS-ZC', .false.) !PAR
        !
        ! -- loop over local topology and fill A0*ZC
        do i = 1, this%cgclnia-1 ! my local blocks
          ic0 = this%cgclia(i)
          ic1 = this%cgclia(i+1)-1 
          !
          ireq0 = lblkeq(i)
          ireq1 = lblkeq(i+1)-1
          !
          do ic = ic0, ic1
            j    = this%cgclja(ic) ! generated local index
            icbg = this%cgclgja(ic) ! global block number
            if (this%cgcg2lid(icbg) > 0) then ! block belongs to me
              icbl = this%cgcg2lid(icbg) ! local block number
              iceq0 = lblkeq(icbl)
              iceq1 = lblkeq(icbl+1)-1
              do ieq = ireq0, ireq1
                tv = DZERO
                jc0 = this%ia(ieq)
                jc1 = this%ia(ieq+1)-1
                do jc = jc0, jc1
                  icol = this%ja(jc)
                  if ((icol >= iceq0).and.(icol <= iceq1)) then
                    tv  = tv + this%a0(jc) * this%zc(icol)
                  end if
                end do
                this%a0zc(j,ieq) = tv
              end do
            else ! block does not belong to me (parallel only)
              do ieq = 1, this%neq
                dwrk(ieq) = DZERO
              end do
              call this%MpiSol%mpi_mv_halo(this%memoryPath, 'IMS-ZC', dwrk, icbg) !PAR
              do ieq = ireq0, ireq1
                this%a0zc(j,ieq) = dwrk(ieq)
              end do
            end if
          end do
          !
        end do
        !
        ! -- clean up
        if (allocated(dwrk)) then
          deallocate(dwrk)
        end if
        !
        ! -- loop over local topology and fill AC coeffients
        m = 0
        do i = 1, this%cgclnia-1 ! my local blocks
          ic0 = this%cgclia(i)
          ic1 = this%cgclia(i+1)-1 
          !
          ireq0 = lblkeq(i)
          ireq1 = lblkeq(i+1)-1
          !
          do ic = ic0, ic1
            j = this%cgclja(ic) ! generated local index
            tv = DZERO
            do ieq = ireq0, ireq1
              tv = tv + this%zc(ieq)*this%a0zc(j,ieq)
            end do
            m = m + 1
            la1d(m) = tv
          end do
        end do
        !
        ! -- MPI all-to-all of matrix coefficients
        call this%MpiSol%mpi_global_exchange_all_cgc(lna1d, this%cgcgnja, &
          la1d, this%acpci, this%cgcgnia, this%cgcgnja, this%cgcgia, this%cgcgja, &
          1)
        !
        do i = 1, ngblk
          do j = 1, ngblk
            this%acpc(j,i) = DZERO
          end do
        end do
        !
        n = 0
        do i = 1, ngblk
          ic0 = this%cgcgia(i)
          ic1 = this%cgcgia(i+1)-1
          do m = ic0, ic1
            j = this%cgcgja(m)
            n = n + 1
            this%acpc(j,i) = this%acpci(n)
          end do
        end do  
        if (this%cgcsol == 1) then
          ! -- band LU decomposition
          call imslinearsub_slu_band(this%acpc, ngblk, p, q)
          this%acpcb = p
          call imslinearsub_slu(this%acpc, ngblk, p, q)
        !
        else
          allocate(dwrk(this%cgcgnja))
          do i = 1, this%cgcgnja
            dwrk(i) = this%acpci(i)
          end do
          ! -- ILU(0) decomposition
          call imslinearsub_pccrs(ngblk, this%cgcgnja, this%cgcgia, this%cgcgja, &
                                  this%cgcgiapc, this%cgcgjapc)
          izero = 0
          delta = dzero
          relax = done
          eqbnd(1) = 1
          eqbnd(2) = ngblk + 1
          call imslinearsub_pcilu0(this%cgcgnja, ngblk, dwrk, this%cgcgia,       &
                                   this%cgcgja, this%acpci, this%cgcgiapc,       &
                                   this%cgcgjapc, this%iw, this%w,               &
                                   relax, izero, delta, &
                                   0, 1, eqbnd)
          deallocate(dwrk)
        end if
        !
        !-- clean up
        deallocate(la1d)
        !
        !-- return
        return
      end subroutine imslinearsub_pcc
      
      subroutine imslinearsub_slu_band(a, n, p, q)
!         ******************************************************************
!         Determine the bands of the matrix.
!         ******************************************************************
!     
!            SPECIFICATIONS:
!         ------------------------------------------------------------------
        implicit none
!         + + + DUMMY ARGUMENTS + + +                                       
        integer(I4b), intent(in) :: n
        integer(I4b), intent(out) :: p
        integer(I4b), intent(out) :: q
        real(DP), dimension(n,n), intent(inout) :: a
!         + + + LOCAL DEFINITIONS + + +                                     
        integer(I4b) :: i, j
!         + + + CODE + + +                                                  
      
        !  check for pivoting and determine the band length
        p = 0 ! lower bandwidth
        q = 0 ! upper bandwidth
        do i = 1, n
          ! check for pivoting
          if (a(i,i) == dzero) then
            write(*,*) 'Error, pivoting required.'
            stop 1
          end if
          do j = 1, i ! lower band
            if (a(j,i) /= dzero) then
              p = max(p, i - j)
            end if
          end do
          do j = i, n ! upper band
            if (a(j,i) /= dzero) then
              q = max(q, j - i)
            end if
          end do
        end do
      
!---------return                                                        
        return 
      end subroutine imslinearsub_slu_band
      
      subroutine imslinearsub_slu(a, n, p, q)
!         ******************************************************************
!         Band LU-decomposition according to Golub & van Loan, p152.
!         ******************************************************************
!     
!            SPECIFICATIONS:
!         ------------------------------------------------------------------
        implicit none
!         + + + DUMMY ARGUMENTS + + +                                       
        integer(I4b), intent(in) :: n
        real(DP), dimension(n,n), intent(inout) :: a
        integer(I4b), intent(in) :: p, q
!         + + + LOCAL DEFINITIONS + + +                                     
        integer(i4b) :: i, j, k, min_p, min_q
!         + + + CODE + + +                                                  
        !
        do k = 1, n - 1
          min_p = min(k + p, n)
          min_q = min(k + q, n)
          do i = k + 1, min_p
            a(k,i) = a(k,i) / a(k,k) ! L
          end do
          do j = k + 1, min_q
            do i = k + 1, min_p 
              a(j,i) = a(j,i) - a(k,i)*a(j,k) ! U
            end do
          end do
        end do
        
!---------return                                                        
        return 
      end subroutine imslinearsub_slu
      
      subroutine imslinearsub_slu_solve(a, n, p, q, x)
!         ******************************************************************
!         Band LU-decomposition according to Golub & van Loan, p152.
!         ******************************************************************
!     
!            SPECIFICATIONS:
!         ------------------------------------------------------------------
        implicit none
!         + + + DUMMY ARGUMENTS + + +                                       
        integer(I4b), intent(in) :: n, p, q
        real(DP), dimension(n,n), intent(in) :: a
        real(DP), dimension(n), intent(inout) :: x
!         + + + LOCAL DEFINITIONS + + +                                     
        integer(i4b) :: i, j, min_p, max_q
!         + + + CODE + + +                                                  
        
        forward: do j = 1, n
          min_p = min(j + p, n)
          do i = j + 1, min_p
            x(i) = x(i) - a(j,i)*x(j)
          end do
        end do forward
        
        backward: do j = n, 1, -1
          x(j) = x(j)/a(j,j)
          max_q = max(1, j - q)
          do i = max_q, j - 1
            x(i) = x(i) - a(j,i)*x(j)
          end do
        end do backward
        
!---------return                                                        
        return 
      end subroutine imslinearsub_slu_solve
!
!-------JACOBI PRECONDITIONER - INVERSE OF DIAGONAL                     
      SUBROUTINE IMSLINEARSUB_PCJ(NJA, NEQ, AMAT, APC, IA, JA) 
!       + + + DUMMY ARGUMENTS + + +                                       
        integer(I4B), INTENT(IN) :: NJA 
        integer(I4B), INTENT(IN) :: NEQ 
        real(DP), DIMENSION(NJA),  INTENT(IN)      :: AMAT 
        real(DP), DIMENSION(NEQ),  INTENT(INOUT)   :: APC 
        integer(I4B), DIMENSION(NEQ+1), INTENT(IN) :: IA 
        integer(I4B), DIMENSION(NJA),   INTENT(IN) :: JA 
!       + + + LOCAL DEFINITIONS + + +                                     
        integer(I4B) :: i, n 
        integer(I4B) :: ic0, ic1 
        integer(I4B) :: id 
        real(DP) :: tv 
!       + + + PARAMETERS + + +                                            
!       + + + FUNCTIONS + + +                                             
!       + + + CODE + + +                                                  
        DO n = 1, NEQ 
            ic0 = IA(n) 
            ic1 = IA(n+1) - 1 
            id = IA(n) 
            DO i = ic0, ic1 
              IF (JA(i) ==  n) THEN 
                id = i 
                EXIT 
              END IF 
            END DO 
            tv  = AMAT(id) 
            IF (ABS( tv ) > DZERO ) tv = DONE / tv 
            APC(n) = tv 
        END DO 
!---------RETURN                                                        
        RETURN 
      END SUBROUTINE IMSLINEARSUB_PCJ 
                                                                        
      SUBROUTINE IMSLINEARSUB_JACA(NEQ, A, D1, D2) 
        IMPLICIT NONE 
!       + + + DUMMY ARGUMENTS + + +                                       
        integer(I4B), INTENT(IN) :: NEQ 
        real(DP), DIMENSION(NEQ),  INTENT(IN)    :: A 
        real(DP), DIMENSION(NEQ),  INTENT(IN)    :: D1 
        real(DP), DIMENSION(NEQ),  INTENT(INOUT) :: D2 
!       + + + LOCAL DEFINITIONS + + +                                     
        integer(I4B) :: n 
        real(DP) :: tv 
!       + + + PARAMETERS + + +                                            
!       + + + FUNCTIONS + + +                                             
!       + + + CODE + + +                                                  
        DO n = 1, NEQ 
          tv     = A(n) * D1(n) 
          D2(n) = tv 
        END DO 
!---------RETURN                                                        
        RETURN 
      END SUBROUTINE IMSLINEARSUB_JACA 
                                                                        
      SUBROUTINE IMSLINEARSUB_PCILU0(NJA, NEQ, AMAT, IA, JA,                    &
                                     APC, IAPC, JAPC, IW, W,                    &
                                     RELAX, IZERO, DELTA,                       &
                                     IBJFLAG, CONVNMOD, CONVMODSTART) !BJ
        IMPLICIT NONE 
!       + + + DUMMY ARGUMENTS + + +                                       
        integer(I4B), INTENT(IN) :: NJA 
        integer(I4B), INTENT(IN) :: NEQ 
        real(DP), DIMENSION(NJA),  INTENT(IN)     :: AMAT 
        integer(I4B), DIMENSION(NEQ+1), INTENT(IN)    :: IA 
        integer(I4B), DIMENSION(NJA), INTENT(IN)      :: JA 
        real(DP), DIMENSION(NJA), INTENT(INOUT)   :: APC 
        integer(I4B), DIMENSION(NEQ+1), INTENT(INOUT) :: IAPC 
        integer(I4B), DIMENSION(NJA), INTENT(INOUT)   :: JAPC 
        integer(I4B), DIMENSION(NEQ), INTENT(INOUT)   :: IW 
        real(DP), DIMENSION(NEQ), INTENT(INOUT)   :: W 
        real(DP), INTENT(IN) :: RELAX 
        integer(I4B), INTENT(INOUT) :: IZERO 
        real(DP), INTENT(IN) :: DELTA 
        integer(I4B), INTENT(IN) :: IBJFLAG !BJ
        integer(I4B), INTENT(IN) :: CONVNMOD !BJ
        integer(I4B), DIMENSION(CONVNMOD+1), INTENT(INOUT) ::CONVMODSTART !BJ
!       + + + LOCAL DEFINITIONS + + +                                     
        integer(I4B) :: ic0, ic1 
        integer(I4B) :: iic0, iic1 
        integer(I4B) :: iu, iiu 
        integer(I4B) :: j, n 
        integer(I4B) :: jj 
        integer(I4B) :: jcol, jw 
        integer(I4B) :: jjcol 
        real(DP) :: drelax 
        real(DP) :: sd1 
        real(DP) :: tl 
        real(DP) :: rs 
        real(DP) :: d 
        integer(I4B) :: im, ieq0, ieq1 !BJ
!       + + + PARAMETERS + + +                                            
!       + + + FUNCTIONS + + +                                             
!       + + + CODE + + +                                                  
        drelax = RELAX 
        DO n = 1, NEQ 
          IW(n)  = 0 
          W(n)   = DZERO 
        END DO 
        
        MAIN: DO im = 1, CONVNMOD !BJ
          ieq0 = CONVMODSTART(im) !BJ
          ieq1 = CONVMODSTART(im+1)-1 !BJ
          DO n = ieq0, ieq1 !BJ
          ic0 = IA(n) 
          ic1 = IA(n+1) - 1 
          DO j = ic0, ic1 
            jcol      = JA(j) 
            IF (IBJFLAG.EQ.1) THEN !BJ
              IF (jcol >= ieq0 .and. jcol <= ieq1) THEN !BJ
                IW(jcol) = 1 !BJ
                W(jcol) = W(jcol) + AMAT(j) !BJ
              END IF !BJ
            ELSE !BJ
            IW(jcol) = 1 
            W(jcol) = W(jcol) + AMAT(j) 
            END IF !BJ
          END DO !BJ
          ic0 = IAPC(n) 
          ic1 = IAPC(n+1) - 1 
          iu  = JAPC(n) 
          rs   = DZERO 
          LOWER: DO j = ic0, iu-1 
            jcol     = JAPC(j) 
            iic0     = IAPC(jcol) 
            iic1     = IAPC(jcol+1) - 1 
            iiu      = JAPC(jcol) 
            tl       = W(jcol) * APC(jcol) 
            W(jcol) = tl 
            DO jj = iiu, iic1 
              jjcol = JAPC(jj) 
              jw    = IW(jjcol) 
              IF (jw.NE.0) THEN 
                W(jjcol) = W(jjcol) - tl * APC(jj) 
              ELSE 
                rs = rs + tl * APC(jj) 
              END IF 
            END DO 
          END DO LOWER 
!           DIAGONAL - CALCULATE INVERSE OF DIAGONAL FOR SOLUTION       
          d   = W(n) 
          tl  = ( DONE + DELTA ) * d - ( drelax * rs ) 
!-----------ENSURE THAT THE SIGN OF THE DIAGONAL HAS NOT CHANGED AND IS 
          sd1 = SIGN(d,tl) 
          IF (sd1.NE.d) THEN 
!             USE SMALL VALUE IF DIAGONAL SCALING IS NOT EFFECTIVE FOR
!             PIVOTS THAT CHANGE THE SIGN OF THE DIAGONAL               
            IF (IZERO > 1) THEN 
              tl = SIGN(DEM6,d) 
!             DIAGONAL SCALING CONTINUES TO BE EFFECTIVE                
            ELSE 
              IZERO = 1 
              EXIT MAIN 
            END IF 
          END IF 
          IF (ABS(tl) ==  DZERO) THEN 
!             USE SMALL VALUE IF DIAGONAL SCALING IS NOT EFFECTIVE FOR
!             ZERO PIVOTS                                               
            IF (IZERO > 1) THEN 
              tl = SIGN(DEM6,d) 
!             DIAGONAL SCALING CONTINUES TO BE EFFECTIVE FOR ELIMINATING 
            ELSE 
              IZERO = 1 
              EXIT MAIN 
            END IF 
          END IF 
          APC(n) = DONE / tl 
!           RESET POINTER FOR IW TO ZERO                                
          IW(n) = 0 
          W(n)  = DZERO 
          DO j = ic0, ic1 
            jcol = JAPC(j) 
            APC(j) = W(jcol) 
            IW(jcol) = 0 
            W(jcol) = DZERO 
          END DO 
        END DO
      END DO MAIN !BJ
!                                                                       
!---------RESET IZERO IF SUCCESSFUL COMPLETION OF MAIN                  
        IZERO = 0 
!                                                                       
!---------RETURN                                                        
        RETURN 
      END SUBROUTINE IMSLINEARSUB_PCILU0 
                                                                        
      SUBROUTINE IMSLINEARSUB_ILU0A(NJA, NEQ, APC, IAPC, JAPC, R, D) 
        IMPLICIT NONE 
!       + + + DUMMY ARGUMENTS + + +                                       
        integer(I4B), INTENT(IN) :: NJA 
        integer(I4B), INTENT(IN) :: NEQ 
        real(DP), DIMENSION(NJA),  INTENT(IN)  :: APC 
        integer(I4B), DIMENSION(NEQ+1), INTENT(IN) :: IAPC 
        integer(I4B), DIMENSION(NJA), INTENT(IN)   :: JAPC 
        real(DP), DIMENSION(NEQ),  INTENT(IN)     :: R 
        real(DP), DIMENSION(NEQ),  INTENT(INOUT)  :: D 
!       + + + LOCAL DEFINITIONS + + +                                     
        integer(I4B) :: ic0, ic1 
        integer(I4B) :: iu 
        integer(I4B) :: jcol 
        integer(I4B) :: j, n 
        real(DP) :: tv 
!       + + + FUNCTIONS + + +                                             
!       + + + CODE + + +                                                  
!         FORWARD SOLVE - APC * D = R                                   
        FORWARD: DO n = 1, NEQ 
          tv   = R(n) 
          ic0 = IAPC(n) 
          ic1 = IAPC(n+1) - 1 
          iu  = JAPC(n) - 1 
          LOWER: DO j = ic0, iu 
            jcol = JAPC(j) 
            tv    = tv - APC(j) * D(jcol) 
          END DO LOWER 
          D(n) = tv 
        END DO FORWARD 
!         BACKWARD SOLVE - D = D / U                                    
        BACKWARD: DO n = NEQ, 1, -1 
          ic0 = IAPC(n) 
          ic1 = IAPC(n+1) - 1 
          iu  = JAPC(n) 
          tv   = D(n) 
          UPPER: DO j = iu, ic1 
            jcol = JAPC(j) 
            tv    = tv - APC(j) * D(jcol) 
          END DO UPPER 
!           COMPUTE D FOR DIAGONAL - D = D / U                          
          D(n) =  tv * APC(n) 
        END DO BACKWARD 
!---------RETURN                                                        
        RETURN 
      END SUBROUTINE IMSLINEARSUB_ILU0A 
      SUBROUTINE IMSLINEARSUB_CG(THIS, ICNVG, ITMAX, INNERIT,                   &
                                 NCONV, CONVNMOD, CONVMODSTART, LOCDV, LOCDR,   &
                                 CACCEL, ITINNER, CONVLOCDV, CONVLOCDR,         &
                                 DVMAX, DRMAX, CONVDVMAX, CONVDRMAX)                                        
        use MpiExchangeGenModule, only: writestd !PAR
        IMPLICIT NONE 
!       + + + DUMMY ARGUMENTS + + +                                       
        CLASS(ImsLinearDataType), INTENT(INOUT) :: THIS !CGC
        integer(I4B), INTENT(INOUT) :: ICNVG 
        integer(I4B), INTENT(IN)    :: ITMAX 
        integer(I4B), INTENT(INOUT) :: INNERIT 
        ! CONVERGENCE INFORMATION
        integer(I4B), INTENT(IN) :: NCONV
        integer(I4B), INTENT(IN) :: CONVNMOD
        integer(I4B), DIMENSION(CONVNMOD+1), INTENT(INOUT) ::CONVMODSTART
        integer(I4B), DIMENSION(CONVNMOD), INTENT(INOUT) :: LOCDV
        integer(I4B), DIMENSION(CONVNMOD), INTENT(INOUT) :: LOCDR
        character(len=31), DIMENSION(NCONV), INTENT(INOUT) :: CACCEL
        integer(I4B), DIMENSION(NCONV), INTENT(INOUT) :: ITINNER
        integer(I4B), DIMENSION(CONVNMOD, NCONV), INTENT(INOUT) :: CONVLOCDV
        integer(I4B), DIMENSION(CONVNMOD, NCONV), INTENT(INOUT) :: CONVLOCDR
        real(DP), DIMENSION(CONVNMOD), INTENT(INOUT) :: DVMAX
        real(DP), DIMENSION(CONVNMOD), INTENT(INOUT) :: DRMAX
        real(DP), DIMENSION(CONVNMOD, NCONV), INTENT(INOUT) :: CONVDVMAX
        real(DP), DIMENSION(CONVNMOD, NCONV), INTENT(INOUT) :: CONVDRMAX
!       + + + LOCAL DEFINITIONS + + + 
        integer(I4B),               POINTER :: NEQ      => NULL()
        integer(I4B),               POINTER :: NJA      => NULL()
        integer(I4B),               POINTER :: NIAPC    => NULL()
        integer(I4B),               POINTER :: NJAPC    => NULL()
        integer(I4B),               POINTER :: IPC      => NULL()
        integer(I4B),               POINTER :: NITERC   => NULL()
        integer(I4B),               POINTER :: ICNVGOPT => NULL() 
        integer(I4B),               POINTER :: NORTH    => NULL()
        real(DP),                   POINTER :: DVCLOSE  => NULL()
        real(DP),                   POINTER :: RCLOSE   => NULL()
        real(DP),                   POINTER :: L2NORM0  => NULL()
        real(DP),                   POINTER :: EPFACT   => NULL()
        integer(I4B), DIMENSION(:), POINTER :: IA0      => NULL() 
        integer(I4B), DIMENSION(:), POINTER :: JA0      => NULL()
        real(DP),     DIMENSION(:), POINTER :: A0       => NULL()
        integer(I4B), DIMENSION(:), POINTER :: IAPC     => NULL()
        integer(I4B), DIMENSION(:), POINTER :: JAPC     => NULL()
        real(DP),     DIMENSION(:), POINTER :: APC      => NULL()
        real(DP),     DIMENSION(:), POINTER :: X        => NULL()
        real(DP),     DIMENSION(:), POINTER :: B        => NULL()
        real(DP),     DIMENSION(:), POINTER :: D        => NULL()
        real(DP),     DIMENSION(:), POINTER :: P        => NULL()
        real(DP),     DIMENSION(:), POINTER :: Q        => NULL()
        real(DP),     DIMENSION(:), POINTER :: Z        => NULL()
        integer(I4B),               POINTER :: NJLU     => NULL()
        integer(I4B), DIMENSION(:), POINTER :: IW       => NULL()
        integer(I4B), DIMENSION(:), POINTER :: JLU      => NULL()
        integer(I4B),               pointer :: IPRIMS   => NULL() !PAR
        type(MpiExchangeType),      pointer :: MpiSol   => NULL() !PAR
        real(DP),                   POINTER :: RHOTOL   => NULL() !SOL
        !
        LOGICAL :: LORTH
        logical :: lsame 
        character(len=31) :: cval
        integer(I4B) :: n 
        integer(I4B) :: iiter 
        integer(I4B) :: xloc, rloc
        integer(I4B) :: im, im0, im1
        real(DP) :: tv 
        real(DP) :: deltax 
        real(DP) :: rmax 
        real(DP) :: l2norm 
        real(DP) :: rcnvg 
        real(DP) :: denom
        real(DP) :: alpha, beta 
        real(DP) :: rho, rho0 
        CHARACTER(LEN=100), DIMENSION(4) :: SA !PAR
!       + + + PARAMETERS + + +                                            
!       + + + FUNCTIONS + + +                                         
!                                                                       
!         + + + CODE + + +                                              
        NEQ      => THIS%NEQ
        NJA      => THIS%NJA
        NIAPC    => THIS%NIAPC
        NJAPC    => THIS%NJAPC
        IPC      => THIS%IPC
        NITERC   => THIS%NITERC
        ICNVGOPT => THIS%ICNVGOPT
        NORTH    => THIS%NORTH
        DVCLOSE  => THIS%DVCLOSE
        RCLOSE   => THIS%RCLOSE
        L2NORM0  => THIS%L2NORM0
        EPFACT   => THIS%EPFACT
        IA0      => THIS%IA0
        JA0      => THIS%JA0
        A0       => THIS%A0
        IAPC     => THIS%IAPC
        JAPC     => THIS%JAPC
        APC      => THIS%APC
        X        => THIS%X
        B        => THIS%RHS
        D        => THIS%D
        P        => THIS%P
        Q        => THIS%Q
        Z        => THIS%Z
        NJLU     => THIS%NJLU
        IW       => THIS%IW
        JLU      => THIS%JLU
        IPRIMS   => THIS%IPRIMS !PAR
        MpiSol   => THIS%MPISOL !PAR
        RHOTOL   => THIS%RHOTOL !SOL
! 
        rho0 = DZERO 
        rho = DZERO 
        INNERIT  = 0 
!                                                                       
!-------INNER ITERATION                                                 
        INNER: DO iiter = 1, itmax 
           INNERIT = INNERIT + 1 
           NITERC  = NITERC  + 1 
!----------APPLY PRECONDITIONER                                         
          SELECT CASE (IPC) 
!             ILU0 AND MILU0              
            CASE (1,2) 
              CALL IMSLINEARSUB_ILU0A(NJA, NEQ, APC, IAPC, JAPC, D, Z) 
!             ILUT AND MILUT
            CASE (3,4)
              CALL IMSLINEARSUB_PCMILUT_LUSOL(NEQ, D, Z, APC, JLU, IW) 
          END SELECT 
          IF (THIS%ICGC == 1) THEN !CGC
!----------APPLY COARSE GRID PRECONDITIONER
            CALL THIS%IMSLINEARSUB_PCC_SOL(NEQ, CONVNMOD, CONVMODSTART, D, Z) !CGC
          END IF !CGC
          rho = IMSLINEARSUB_DP(NEQ, D, Z) 
          ! -- MPI parallel: collective comm. of rho (sum)
          call MpiSol%mpi_global_exchange_all_sum(rho) !PAR
!-----------COMPUTE DIRECTIONAL VECTORS                                 
          IF (IITER ==  1) THEN 
            DO n = 1, NEQ 
              P(n) = Z(n) 
            END DO 
          ELSE
            !denom = rho0 + SIGN(DPREC,rho0)
            !beta = rho / denom
            beta = rho / rho0 
            DO n = 1, NEQ 
              P(n) = Z(n) + beta * P(n) 
            END DO 
          END IF 
!-----------COMPUTE ITERATES                                            
!           UPDATE Q                                                   
          CALL IMSLINEARSUB_MV(NJA, NEQ, A0, P, Q, IA0, JA0) 
          ! -- MPI parallel: point-to-point comm. of THIS%X and update of THIS%D
          call MpiSol%mpi_local_exchange(this%memoryPath, 'IMS-P', .false.) !PAR
          call MpiSol%mpi_mv_halo(this%memoryPath, 'IMS-P', q) !PAR
          denom =  IMSLINEARSUB_DP(NEQ, P, Q)
          ! -- MPI parallel: collective comm. of denom (sum)
          call MpiSol%mpi_global_exchange_all_sum(denom) !PAR
          denom = denom + SIGN(DPREC, denom) 
          alpha = rho / denom
!-----------UPDATE X AND RESIDUAL                                       
          deltax = DZERO 
          rmax   = DZERO 
          l2norm = DZERO 
          DO im = 1, CONVNMOD
            DVMAX(im) = DZERO
            DRMAX(im) = DZERO
          END DO
          im  = 1
          im0 = CONVMODSTART(1)
          im1 = CONVMODSTART(2)
          DO n = 1, NEQ
            ! -- determine current model index
            if (n == im1) then
              im = im + 1
              im0 = CONVMODSTART(im)
              im1 = CONVMODSTART(im+1)
            end if
            ! -- identify deltax and rmax
            tv      = alpha * P(n) 
            X(n)  = X(n) + tv 
            IF (ABS(tv) > ABS(deltax)) THEN
              deltax = tv
              xloc = n
            END IF
            IF (ABS(tv) > ABS(DVMAX(im))) THEN
              DVMAX(im) = tv
              LOCDV(im) = n
            END IF
            tv      = D(n) 
            tv      = tv - alpha * Q(n) 
            D(n)  = tv 
            IF (ABS(tv) > ABS(rmax)) THEN
              rmax = tv
              rloc = n
            END IF
            IF (ABS(tv) > ABS(DRMAX(im))) THEN
              DRMAX(im) = tv
              LOCDR(im) = n
            END IF
            l2norm = l2norm + tv * tv 
          END DO 
          ! -- MPI parallel: collective comm. of deltax and rmax (max)
          call MpiSol%mpi_global_exchange_all_absmax(deltax, rmax) !PAR
          ! -- MPI parallel: collective comm. of l2norm (sum)
          call MpiSol%mpi_global_exchange_all_sum(l2norm) !PAR
          l2norm = SQRT(l2norm) 
!-----------SAVE SOLVER CONVERGENCE INFORMATION
          IF (NCONV > 1) THEN
            n = NITERC
            WRITE(cval, '(g15.7)') alpha
            CACCEL(n) = cval
            ITINNER(n)    = iiter
            DO im = 1, CONVNMOD
              CONVLOCDV(im, n) = LOCDV(im)
              CONVLOCDR(im, n) = LOCDR(im)
              CONVDVMAX(im, n) = DVMAX(im)
              CONVDRMAX(im, n) = DRMAX(im)
            END DO
          END IF
!-----------TEST FOR SOLVER CONVERGENCE                                 
          IF (ICNVGOPT ==  2 .OR. ICNVGOPT ==  3 .OR. ICNVGOPT ==  4) THEN 
            rcnvg = l2norm 
          ELSE 
            rcnvg = rmax 
          END IF 
          CALL IMSLINEARSUB_TESTCNVG(ICNVGOPT, ICNVG, INNERIT,                  &
                                     deltax, rcnvg,                             &
                                     L2NORM0, EPFACT, DVCLOSE, RCLOSE) 
!
!-----------CHECK FOR EXACT SOLUTION                                    
          IF (rcnvg == DZERO) ICNVG = 1 
!          
!-----------CHECK FOR STANDARD CONVERGENCE 
          IF (ICNVG.NE.0) EXIT INNER 
!
!-----------CHECK THAT CURRENT AND PREVIOUS rho ARE DIFFERENT           
          lsame = IS_SAME(rho, rho0, RHOTOL) !SOL
          IF (lsame) THEN 
            EXIT INNER 
          END IF 
!
!-----------RECALCULATE THE RESIDUAL
          IF (NORTH > 0) THEN
            LORTH = mod(iiter+1,NORTH) == 0
            IF (LORTH) THEN
              CALL IMSLINEARSUB_MV(NJA, NEQ, A0, X, D, IA0, JA0)
             ! -- MPI parallel: point-to-point comm. of X and update of D
              call MpiSol%mpi_local_exchange(this%memoryPath, 'IMS-X', .false.) !PAR
              call MpiSol%mpi_mv_halo(this%memoryPath, 'IMS-X', d) !PAR
 
            END IF
          END IF
!-----------SAVE CURRENT INNER ITERATES                                 
          rho0 = rho 
        END DO INNER 

        if (iprims == 3 .and. writestd) then
          write(sa(1),*) real(deltax)
          write(sa(2),*) real(rmax)
          write(*,'(1x,a,1x,i5.5,a,1x,i5.5,a,1x,2a,1x,a)')                      &
          'cg: oit, iit, bih, bir =',                                           &
          niterc,',',innerit,',',trim(adjustl(sa(1))),                          &
          ',',trim(adjustl(sa(2)))
!
        endif
!---------RESET ICNVG        
        IF (ICNVG < 0) ICNVG = 0
!---------RETURN                                                        
        RETURN 
      END SUBROUTINE IMSLINEARSUB_CG 
                                                                        
      SUBROUTINE IMSLINEARSUB_BCGS(THIS, ICNVG, ITMAX, INNERIT,                 &
                                   NCONV, CONVNMOD, CONVMODSTART, LOCDV, LOCDR, &
                                   CACCEL, ITINNER, CONVLOCDV, CONVLOCDR,       &
                                   DVMAX, DRMAX, CONVDVMAX, CONVDRMAX)                                
        use MpiExchangeGenModule, only: writestd !PAR

        IMPLICIT NONE 
!       + + + DUMMY ARGUMENTS + + +                                       
        CLASS(ImsLinearDataType), INTENT(INOUT) :: THIS !CGC
        integer(I4B), INTENT(INOUT) :: ICNVG 
        integer(I4B), INTENT(IN)    :: ITMAX 
        integer(I4B), INTENT(INOUT) :: INNERIT 
        ! CONVERGENCE INFORMATION
        integer(I4B), INTENT(IN) :: NCONV
        integer(I4B), INTENT(IN) :: CONVNMOD
        integer(I4B), DIMENSION(CONVNMOD+1), INTENT(INOUT) ::CONVMODSTART
        integer(I4B), DIMENSION(CONVNMOD), INTENT(INOUT) :: LOCDV
        integer(I4B), DIMENSION(CONVNMOD), INTENT(INOUT) :: LOCDR
        character(len=31), DIMENSION(NCONV), INTENT(INOUT) :: CACCEL
        integer(I4B), DIMENSION(NCONV), INTENT(INOUT) :: ITINNER
        integer(I4B), DIMENSION(CONVNMOD, NCONV), INTENT(INOUT) :: CONVLOCDV
        integer(I4B), DIMENSION(CONVNMOD, NCONV), INTENT(INOUT) :: CONVLOCDR
        real(DP), DIMENSION(CONVNMOD), INTENT(INOUT) :: DVMAX
        real(DP), DIMENSION(CONVNMOD), INTENT(INOUT) :: DRMAX
        real(DP), DIMENSION(CONVNMOD, NCONV), INTENT(INOUT) :: CONVDVMAX
        real(DP), DIMENSION(CONVNMOD, NCONV), INTENT(INOUT) :: CONVDRMAX
!       + + + LOCAL DEFINITIONS + + +  
        integer(I4B),               POINTER :: NEQ
        integer(I4B),               POINTER :: NJA
        integer(I4B),               POINTER :: NIAPC
        integer(I4B),               POINTER :: NJAPC
        integer(I4B),               POINTER :: IPC
        integer(I4B),               POINTER :: NITERC
        integer(I4B),               POINTER :: ICNVGOPT
        integer(I4B),               POINTER :: NORTH
        integer(I4B),               POINTER :: ISCL
        real(DP),     DIMENSION(:), POINTER :: DSCALE
        real(DP),                   POINTER :: DVCLOSE
        real(DP),                   POINTER :: RCLOSE
        real(DP),                   POINTER :: L2NORM0
        real(DP),                   POINTER :: EPFACT
        integer(I4B), DIMENSION(:), POINTER :: IA0
        integer(I4B), DIMENSION(:), POINTER :: JA0
        real(DP),     DIMENSION(:), POINTER :: A0
        integer(I4B), DIMENSION(:), POINTER :: IAPC
        integer(I4B), DIMENSION(:), POINTER :: JAPC
        real(DP),     DIMENSION(:), POINTER :: APC
        real(DP),     DIMENSION(:), POINTER :: X
        real(DP),     DIMENSION(:), POINTER :: B
        real(DP),     DIMENSION(:), POINTER :: D
        real(DP),     DIMENSION(:), POINTER :: P
        real(DP),     DIMENSION(:), POINTER :: Q
        real(DP),     DIMENSION(:), POINTER :: T
        real(DP),     DIMENSION(:), POINTER :: V
        real(DP),     DIMENSION(:), POINTER :: DHAT
        real(DP),     DIMENSION(:), POINTER :: PHAT
        real(DP),     DIMENSION(:), POINTER :: QHAT
        ! ILUT
        integer(I4B),               POINTER :: NJLU
        integer(I4B), DIMENSION(:), POINTER :: IW
        integer(I4B), DIMENSION(:), POINTER :: JLU
        integer(I4B),               POINTER :: IPRIMS !PAR
        type(MpiExchangeType),      POINTER :: MpiSol !PAR
        REAL(DP),                   POINTER :: RHOTOL !SOL
        REAL(DP),                   POINTER :: ALPHATOL !SOL
        REAL(DP),                   POINTER :: OMEGATOL !SOL
        !
        LOGICAL :: LORTH
        logical :: lsame 
        character(len=15) :: cval1, cval2
        integer(I4B) :: n 
        integer(I4B) :: iiter 
        integer(I4B) :: xloc, rloc
        integer(I4B) :: im, im0, im1
        real(DP) :: tv 
        real(DP) :: deltax 
        real(DP) :: rmax
        real(DP) :: l2norm 
        real(DP) :: rcnvg 
        real(DP) :: alpha, alpha0 
        real(DP) :: beta 
        real(DP) :: rho, rho0 
        real(DP) :: omega, omega0 
        real(DP) :: numer, denom
        CHARACTER(LEN=100), DIMENSION(4) :: SA !PAR
        
!       + + + PARAMETERS + + +                                            
!       + + + FUNCTIONS + + +                                         
!                                                                       
!         + + + CODE + + +                                              
        NEQ      => THIS%NEQ
        NJA      => THIS%NJA
        NIAPC    => THIS%NIAPC
        NJAPC    => THIS%NJAPC
        IPC      => THIS%IPC
        NITERC   => THIS%NITERC
        ICNVGOPT => THIS%ICNVGOPT
        NORTH    => THIS%NORTH
        ISCL     => THIS%ISCL
        DSCALE   => THIS%DSCALE
        DVCLOSE  => THIS%DVCLOSE
        RCLOSE   => THIS%RCLOSE
        L2NORM0  => THIS%L2NORM0
        EPFACT   => THIS%EPFACT
        IA0      => THIS%IA0
        JA0      => THIS%JA0
        A0       => THIS%A0
        IAPC     => THIS%IAPC
        JAPC     => THIS%JAPC
        APC      => THIS%APC
        X        => THIS%X
        B        => THIS%RHS
        D        => THIS%D
        P        => THIS%P
        Q        => THIS%Q
        T        => THIS%T
        V        => THIS%V
        DHAT     => THIS%DHAT
        PHAT     => THIS%PHAT
        QHAT     => THIS%QHAT
        NJLU     => THIS%NJLU
        IW       => THIS%IW
        JLU      => THIS%JLU
        IPRIMS   => THIS%IPRIMS !PAR
        MpiSol   => THIS%MpiSol !PAR
        RHOTOL   => THIS%RHOTOL !SOL
        ALPHATOL => THIS%ALPHATOL !SOL
        OMEGATOL => THIS%OMEGATOL !SOL
        
        INNERIT  = 0 
                                                                        
        alpha  = DZERO
        alpha0 = DZERO
        beta   = DZERO 
        rho    = DZERO 
        rho0   = DZERO
        omega  = DZERO
        omega0 = DZERO
        deltax = DZERO !PAR
        rmax   = DZERO !PAR
!                                                                       
!-------SAVE INITIAL RESIDUAL                                           
        DO n = 1, NEQ 
          DHAT(n) = D(n)
        END DO
!                                                                       
!-------INNER ITERATION                                                 
        INNER: DO iiter = 1, itmax 
           INNERIT = INNERIT + 1 
           NITERC = NITERC + 1 
!----------CALCULATE rho                                                
          rho = IMSLINEARSUB_DP(NEQ, DHAT, D) 
          ! -- MPI parallel: collective comm. of rho (sum)
          call MpiSol%mpi_global_exchange_all_sum(rho) !PAR
!-----------COMPUTE DIRECTIONAL VECTORS                                 
          IF (IITER ==  1) THEN 
            DO n = 1, NEQ 
              P(n) = D(n) 
            END DO 
          ELSE 
            rho0   = rho0 + SIGN(DPREC, rho0) !PAR
            omega0 = omega0 + SIGN(DPREC, omega0) !PAR
            beta = ( rho / rho0 ) * ( alpha0 / omega0 ) 
            DO n = 1, NEQ 
              P(n) = D(n) + beta * ( P(n) - omega0 * V(n) ) 
            END DO 
          END IF 
!----------APPLY PRECONDITIONER TO UPDATE PHAT                          
          SELECT CASE (IPC) 
!             ILU0 AND MILU0
            CASE (1,2) 
              CALL IMSLINEARSUB_ILU0A(NJA, NEQ, APC, IAPC, JAPC, P, PHAT) 
!             ILUT AND MILUT
            CASE (3,4)
              CALL IMSLINEARSUB_PCMILUT_LUSOL(NEQ, P, PHAT, APC, JLU, IW) 
          END SELECT 
          IF (THIS%ICGC == 1) THEN !CGC
!----------APPLY COARSE GRID PRECONDITIONER
            CALL THIS%IMSLINEARSUB_PCC_SOL(NEQ, CONVNMOD, CONVMODSTART, P, PHAT) !CGC
          END IF !CGC
!-----------COMPUTE ITERATES                                            
!           UPDATE V WITH A AND PHAT                                    
          CALL IMSLINEARSUB_MV(NJA, NEQ, A0, PHAT, V, IA0, JA0) 
          ! -- MPI parallel: point-to-point comm. of PHAT and update of V
          call MpiSol%mpi_local_exchange(this%memoryPath, 'IMS-PHAT', .false.) !PAR
          call MpiSol%mpi_mv_halo(this%memoryPath, 'IMS-PHAT', v) !PAR
!           UPDATE alpha WITH DHAT AND V                                
          denom = IMSLINEARSUB_DP(NEQ, DHAT, V) 
          ! -- MPI parallel: collective comm. of denom (sum)
          call MpiSol%mpi_global_exchange_all_sum(denom) !PAR
          denom = denom + SIGN(DPREC, denom) 
          alpha = rho / denom 
!-----------UPDATE Q                                                    
          DO n = 1, NEQ 
            Q(n) = D(n) - alpha * V(n)  
          END DO 
!!-----------CALCULATE INFINITY NORM OF Q - TEST FOR TERMINATION         
!!           TERMINATE IF rmax IS LESS THAN MACHINE PRECISION (DPREC) 
!          rmax = DZERO 
!          DO n = 1, NEQ 
!              tv = Q(n) 
!              IF (ISCL.NE.0 ) tv = tv / DSCALE(n) 
!              IF (ABS(tv) > ABS(rmax) ) rmax = tv 
!          END DO 
!          IF (ABS(rmax).LE.DPREC) THEN 
!            deltax = DZERO 
!            DO n = 1, NEQ 
!              tv      = alpha * PHAT(n) 
!              IF (ISCL.NE.0) THEN 
!                tv = tv * DSCALE(n) 
!              END IF 
!              X(n)  = X(n) + tv 
!              IF (ABS(tv) > ABS(deltax) ) deltax = tv 
!            END DO 
!            CALL IMSLINEARSUB_TESTCNVG(ICNVGOPT, ICNVG, INNERIT,                &
!                                       deltax, rmax,                            &
!                                       rmax, EPFACT, DVCLOSE, RCLOSE )          
!            IF (ICNVG.NE.0 ) EXIT INNER 
!          END IF 
!-----------APPLY PRECONDITIONER TO UPDATE QHAT                         
          SELECT CASE (IPC) 
!            ILU0 AND MILU0            
            CASE (1,2) 
              CALL IMSLINEARSUB_ILU0A(NJA, NEQ, APC, IAPC, JAPC, Q, QHAT) 
!             ILUT AND MILUT
            CASE (3,4)
              CALL IMSLINEARSUB_PCMILUT_LUSOL(NEQ, Q, QHAT, APC, JLU, IW)
          END SELECT
!----------APPLY COARSE GRID PRECONDITIONER
          IF (THIS%ICGC == 1) THEN !CGC
            CALL THIS%IMSLINEARSUB_PCC_SOL(NEQ, CONVNMOD, CONVMODSTART, Q, QHAT) !CGC
          END IF !CGC
!           UPDATE T WITH A AND QHAT                                    
          CALL IMSLINEARSUB_MV(NJA, NEQ, A0, QHAT, T, IA0, JA0) 
          ! -- MPI parallel: point-to-point comm. of QHAT and update of T
          call MpiSol%mpi_local_exchange(this%memoryPath, 'IMS-QHAT', .false.) !PAR
          call MpiSol%mpi_mv_halo(this%memoryPath, 'IMS-QHAT', t) !PAR
!-----------UPDATE omega                                                
          numer = IMSLINEARSUB_DP(NEQ, T, Q) 
          denom = IMSLINEARSUB_DP(NEQ, T, T)
          ! -- MPI parallel: collective comm. of numer and denom (sum)
          call MpiSol%mpi_global_exchange_all_sum(numer, denom) !PAR
          denom = denom + SIGN(DPREC,denom) 
          omega = numer / denom 
!-----------UPDATE X AND RESIDUAL                                       
          deltax = DZERO 
          rmax   = DZERO 
          l2norm = DZERO 
          DO im = 1, CONVNMOD
            DVMAX(im) = DZERO
            DRMAX(im) = DZERO
          END DO
          im  = 1
          im0 = CONVMODSTART(1)
          im1 = CONVMODSTART(2)
          DO n = 1, NEQ 
            ! -- determine current model index
            if (n == im1) then
              im = im + 1
              im0 = CONVMODSTART(im)
              im1 = CONVMODSTART(im+1)
            end if
!-------------X AND DX                                                  
            tv = alpha * PHAT(n) + omega * QHAT(n) 
            X(n)  = X(n) + tv 
            IF (ISCL.NE.0) THEN 
              tv = tv * DSCALE(n) 
            END IF 
            IF (ABS(tv) > ABS(deltax)) THEN
              deltax = tv
              xloc = n
            END IF
            IF (ABS(tv) > ABS(DVMAX(im))) THEN
              DVMAX(im) = tv
              LOCDV(im) = n
            END IF
!-------------RESIDUAL                                                  
            tv      = Q(n) - omega * T(n) 
            D(n)  = tv 
            IF (ISCL.NE.0) THEN 
              tv = tv / DSCALE(n) 
            END IF 
            IF (ABS(tv) > ABS(rmax)) THEN
              rmax = tv
              rloc = n
            END IF
            IF (ABS(tv) > ABS(DRMAX(im))) THEN
              DRMAX(im) = tv
              LOCDR(im) = n
            END IF
            l2norm = l2norm + tv * tv 
          END DO
          ! -- MPI parallel: collective comm. of l2norm (sum) and rmax (max)
          call MpiSol%mpi_global_exchange_all_sum(l2norm) !PAR
          call MpiSol%mpi_global_exchange_all_absmax(deltax, rmax) !PAR
          l2norm = sqrt(l2norm)
!-----------SAVE SOLVER CONVERGENCE INFORMATION
          IF (NCONV > 1) THEN
            n = NITERC
            WRITE(cval1,'(g15.7)') alpha
            WRITE(cval2,'(g15.7)') omega
            CACCEL(n) = trim(adjustl(cval1)) // ',' // trim(adjustl(cval2))
            ITINNER(n)    = iiter
            DO im = 1, CONVNMOD
              CONVLOCDV(im, n) = LOCDV(im)
              CONVLOCDR(im, n) = LOCDR(im)
              CONVDVMAX(im, n) = DVMAX(im)
              CONVDRMAX(im, n) = DRMAX(im)
            END DO
          END IF
!-----------TEST FOR SOLVER CONVERGENCE                                 
          IF (ICNVGOPT ==  2 .OR. ICNVGOPT ==  3 .OR. ICNVGOPT ==  4) THEN 
            rcnvg = l2norm 
          ELSE 
            rcnvg = rmax 
          END IF 
          CALL IMSLINEARSUB_TESTCNVG(ICNVGOPT, ICNVG, INNERIT,                  &
                                     deltax, rcnvg,                             &
                                     L2NORM0, EPFACT, DVCLOSE, RCLOSE)         
!          
!-----------CHECK FOR EXACT SOLUTION                                    
          IF (rcnvg == DZERO) ICNVG = 1 
!          
!-----------CHECK FOR STANDARD CONVERGENCE 
          IF (ICNVG.NE.0) EXIT INNER
!
!-----------CHECK THAT CURRENT AND PREVIOUS rho, alpha, AND omega ARE 
!           DIFFERENT
          lsame = IS_SAME(rho, rho0, RHOTOL) !SOL
          IF (lsame) THEN 
            EXIT INNER 
          END IF 
          lsame = IS_SAME(alpha, alpha0, ALPHATOL) !SOL
          IF (lsame) THEN 
            EXIT INNER 
          END IF 
          lsame = IS_SAME(omega, omega0, OMEGATOL) !SOL
          IF (lsame) THEN 
            EXIT INNER 
          END IF 
!-----------RECALCULATE THE RESIDUAL
          IF (NORTH > 0) THEN
            LORTH = mod(iiter+1,NORTH) == 0
            IF (LORTH) THEN
              CALL IMSLINEARSUB_MV(NJA, NEQ, A0,X , D, IA0, JA0)
              ! -- MPI parallel: point-to-point comm. of X and update of D
              call MpiSol%mpi_local_exchange(this%memoryPath, 'IMS-X', .false.) !PAR
              call MpiSol%mpi_mv_halo(this%memoryPath, 'IMS-X', d) !PAR
              CALL IMSLINEARSUB_AXPY(NEQ, B, -DONE, D, D)
              !DO n = 1, NEQ
              !  tv   = D(n)
              !  D(n) = B(n) - tv
              !END DO
            END IF
          END IF
!-----------SAVE CURRENT INNER ITERATES                                 
          rho0   = rho
          alpha0 = alpha
          omega0 = omega
        END DO INNER
!                                                                       
        if (iprims == 3 .and. writestd) then !PAR
          write(sa(1),*) real(deltax) !PAR
          write(sa(2),*) real(rmax) !PAR
          write(*,'(1x,a,1x,i5.5,a,1x,i5.5,a,1x,2a,1x,a)')                  & !PAR
          'bicgstab: oit, iit, bih, bir =',                                 & !PAR
          niterc,',',innerit,',',trim(adjustl(sa(1))),                      & !PAR
          ',',trim(adjustl(sa(2))) !PAR
        endif !PAR
!
!---------RESET ICNVG
        IF (ICNVG < 0) ICNVG = 0
!---------RETURN                                                        
        RETURN 
      END SUBROUTINE IMSLINEARSUB_BCGS 
!                                                                       
!---------TEST FOR SOLVER CONVERGENCE                                   
        SUBROUTINE IMSLINEARSUB_TESTCNVG(Icnvgopt, Icnvg, Iiter,                &
                                         Hmax, Rmax,                            &
                                         Rmax0, Epfact, Dvclose, Rclose )         
        IMPLICIT NONE 
!       + + + DUMMY ARGUMENTS + + +                                       
        integer(I4B), INTENT(IN)         :: Icnvgopt 
        integer(I4B), INTENT(INOUT)      :: Icnvg 
        integer(I4B), INTENT(IN)         :: Iiter
        real(DP), INTENT(IN) :: Hmax 
        real(DP), INTENT(IN) :: Rmax 
        real(DP), INTENT(IN) :: Rmax0 
        real(DP), INTENT(IN) :: Epfact 
        real(DP), INTENT(IN) :: Dvclose 
        real(DP), INTENT(IN) :: Rclose 
!       + + + LOCAL DEFINITIONS + + +                                     
!       + + + FUNCTIONS + + +                                             
!       + + + CODE + + +                                                  
        IF (Icnvgopt ==  0) THEN 
          IF (ABS(Hmax) <= Dvclose .AND. ABS(Rmax) <= Rclose) THEN 
            Icnvg = 1 
          END IF 
        ELSE IF (Icnvgopt == 1) THEN 
          IF (ABS(Hmax) <= Dvclose .AND. ABS(Rmax) <= Rclose) THEN
            IF (iiter == 1) THEN 
              Icnvg = 1
            ELSE
              Icnvg = -1
            END IF
          END IF 
        ELSE IF (Icnvgopt == 2) THEN 
          IF (ABS(Hmax) <= Dvclose .OR. Rmax <= Rclose) THEN 
            Icnvg = 1 
          ELSE IF (Rmax <= Rmax0*Epfact) THEN 
            Icnvg = -1 
          END IF
        ELSE IF (Icnvgopt == 3) THEN 
          IF (ABS(Hmax) <= Dvclose) THEN
            Icnvg = 1 
          ELSE IF (Rmax <= Rmax0*Rclose) THEN  
            Icnvg = -1 
          END IF
        ELSE IF (Icnvgopt == 4) THEN 
          IF (ABS(Hmax) <= Dvclose .AND. Rmax <= Rclose) THEN
            Icnvg = 1 
          ELSE IF (Rmax <=  Rmax0*Epfact) THEN  
            Icnvg = -1 
          END IF
       END IF 
!---------RETURN                                                        
        RETURN 
      END SUBROUTINE IMSLINEARSUB_TESTCNVG 
!                                                                       
!---------GENERATE IAPC AND JAPC FROM IA AND JA                         
!         JAPC(1:NEQ) HAS THE POSITION OF THE UPPER ENTRY FOR A ROW     
!         JAPC(NEQ+1:NJA) IS THE COLUMN POSITION FOR ENTRY              
!         APC(1:NEQ) PRECONDITIONED INVERSE OF THE DIAGONAL             
!         APC(NEQ+1:NJA) PRECONDITIONED ENTRIES FOR OFF DIAGONALS       
        SUBROUTINE IMSLINEARSUB_PCCRS(NEQ, NJA, IA, JA,                         &
                                      IAPC,JAPC)                               
        IMPLICIT NONE 
!       + + + DUMMY ARGUMENTS + + +                                       
        integer(I4B), INTENT(IN)         :: NEQ 
        integer(I4B), INTENT(IN)         :: NJA 
        integer(I4B), DIMENSION(NEQ+1), INTENT(IN)    :: IA 
        integer(I4B), DIMENSION(NJA), INTENT(IN)      :: JA 
        integer(I4B), DIMENSION(NEQ+1), INTENT(INOUT) :: IAPC 
        integer(I4B), DIMENSION(NJA), INTENT(INOUT)   :: JAPC 
!       + + + LOCAL DEFINITIONS + + +                                     
        integer(I4B) :: n, j 
        integer(I4B) :: i0, i1 
        integer(I4B) :: nlen 
        integer(I4B) :: ic,ip 
        integer(I4B) :: jcol 
        integer(I4B), DIMENSION(:), ALLOCATABLE :: iarr 
!       + + + FUNCTIONS + + +                                             
!       + + + CODE + + +                                                  
        ip = NEQ + 1 
        DO n = 1, NEQ 
          i0 = IA(n) 
          i1 = IA(n+1) - 1 
          nlen = i1 - i0 
          ALLOCATE( iarr(nlen) ) 
          ic = 0 
          DO j = i0, i1 
            jcol = JA(j) 
            IF (jcol ==  n) CYCLE 
            ic = ic + 1 
            iarr(ic) = jcol 
          END DO 
          CALL IMSLINEARSUB_ISORT(nlen,iarr) 
          IAPC(n) = ip 
          DO j = 1, nlen 
            jcol = iarr(j) 
            JAPC(ip) = jcol 
            ip = ip + 1 
          END DO 
          DEALLOCATE(iarr) 
        END DO 
        IAPC(NEQ+1) = NJA + 1 
!---------POSITION OF THE FIRST UPPER ENTRY FOR ROW                     
        DO n = 1, NEQ 
          i0 = IAPC(n) 
          i1 = IAPC(n+1) - 1 
          JAPC(n) = IAPC(n+1) 
          DO j = i0, i1 
            jcol = JAPC(j) 
            IF (jcol > n) THEN 
              JAPC(n) = j 
              EXIT 
            END IF 
          END DO 
        END DO 
!---------RETURN                                                        
        RETURN 
      END SUBROUTINE IMSLINEARSUB_PCCRS 
!                                                                       
!-------SIMPLE IN-PLACE SORTING ROUTINE FOR AN INTEGER ARRAY             
      SUBROUTINE IMSLINEARSUB_ISORT(NVAL, IARRAY) 
        IMPLICIT NONE 
!       + + + DUMMY ARGUMENTS + + +                                       
        integer(I4B),INTENT(IN) :: NVAL 
        integer(I4B),DIMENSION(NVAL),INTENT(INOUT) :: IARRAY 
!       + + + LOCAL DEFINITIONS + + +                                     
        integer(I4B) :: i, j, itemp 
!       + + + FUNCTIONS + + +                                             
!       + + + CODE + + +                                                  
        DO i = 1, NVAL-1 
            DO j = i+1, NVAL 
                if(IARRAY(i) > IARRAY(j)) then 
                    itemp = IARRAY(j) 
                    IARRAY(j) = IARRAY(i) 
                    IARRAY(i) = itemp 
                END IF 
            END DO                                                      
        END DO                                                          
!---------RETURN                                                        
        RETURN                                                          
      END SUBROUTINE IMSLINEARSUB_ISORT
!      
!-------INITIALIZE REAL VECTOR      
      SUBROUTINE IMSLINEARSUB_SETX(NR, D1, C)
        IMPLICIT NONE
!     + + + DUMMY ARGUMENTS + + +
        integer(I4B), INTENT(IN) :: NR
        real(DP), DIMENSION(NR),  INTENT(INOUT) :: D1
        real(DP),  INTENT(IN)                   :: C
!     + + + LOCAL DEFINITIONS + + +
        INTEGER :: n
!     + + + FUNCTIONS + + +
!     + + + CODE + + +
!
        DO n = 1, NR
          D1(n) = C
        END DO
!---------RETURN
        RETURN
      END SUBROUTINE IMSLINEARSUB_SETX
                                                                        
!       COPY ONE real(DP) VECTOR TO ANOTHER                      
      SUBROUTINE IMSLINEARSUB_DCOPY(NR, V, R)                                 
        IMPLICIT NONE                                                   
!       + + + DUMMY ARGUMENTS + + +                                       
        integer(I4B), INTENT(IN) :: NR                                       
        real(DP), DIMENSION(NR), INTENT(IN)    :: V              
        real(DP), DIMENSION(NR), INTENT(INOUT) :: R              
!       + + + LOCAL DEFINITIONS + + +                                     
        integer(I4B) :: n                                                    
!       + + + FUNCTIONS + + +                                             
!       + + + CODE + + +                                                  
        DO n = 1, NR                                                    
          R(n) = V(n)                                                   
        END DO                                                          
!---------RETURN                                                        
        RETURN                                                          
      END SUBROUTINE IMSLINEARSUB_DCOPY                                       
                                                                        
!       COPY ONE INTEGER VECTOR TO ANOTHER                              
      SUBROUTINE IMSLINEARSUB_ICOPY(NR, V, R)                                 
        IMPLICIT NONE                                                   
!       + + + DUMMY ARGUMENTS + + +                                       
        integer(I4B), INTENT(IN) :: NR                                       
        integer(I4B), DIMENSION(NR), INTENT(IN)    :: V                      
        integer(I4B), DIMENSION(NR), INTENT(INOUT) :: R                      
!       + + + LOCAL DEFINITIONS + + +                                     
        integer(I4B) :: n                                                    
!       + + + FUNCTIONS + + +                                             
!       + + + CODE + + +                                                  
        DO n = 1, NR                                                    
          R(n) = V(n)                                                   
        END DO                                                          
!---------RETURN                                                        
        RETURN                                                          
      END SUBROUTINE IMSLINEARSUB_ICOPY                                       
!      
!-------SCALE A REAL VECTOR WITH A CONSTANT      
      SUBROUTINE IMSLINEARSUB_RSCAL(NR, C, D1)
        IMPLICIT NONE
!     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NR
        real(DP), INTENT(IN) :: C
        real(DP), DIMENSION(NR),  INTENT(INOUT) :: D1
!     + + + LOCAL DEFINITIONS + + +
        INTEGER :: n
!     + + + FUNCTIONS + + +
!     + + + CODE + + +
        DO n = 1, NR
          D1(n) = C * D1(n)
        END DO
!---------RETURN
        RETURN
      END SUBROUTINE IMSLINEARSUB_RSCAL

      
      SUBROUTINE IMSLINEARSUB_MV(NJA, NEQ, A, D1, D2, IA, JA)                        
        IMPLICIT NONE                                                   
!       + + + DUMMY ARGUMENTS + + +                                       
        integer(I4B), INTENT(IN) :: NJA                                      
        integer(I4B), INTENT(IN) :: NEQ                                      
        real(DP), DIMENSION(NJA),  INTENT(IN)    :: A            
        real(DP), DIMENSION(NEQ),  INTENT(IN)    :: D1           
        real(DP), DIMENSION(NEQ),  INTENT(INOUT) :: D2           
        integer(I4B), DIMENSION(NEQ+1), INTENT(IN) :: IA                     
        integer(I4B), DIMENSION(NJA), INTENT(IN)   :: JA                     
!       + + + LOCAL DEFINITIONS + + +                                     
        integer(I4B) :: ic0, ic1                                             
        integer(I4B) :: icol                                                 
        integer(I4B) :: m, n                                                 
        real(DP) :: tv                                           
!       + + + PARAMETERS + + +                                            
!       + + + FUNCTIONS + + +                                             
!       + + + CODE + + +                                                  
        DO n = 1, NEQ                                                   
!           ADD DIAGONAL AND OFF-DIAGONAL TERMS                         
          tv     = DZERO                                                
          ic0   = IA(n)                                                 
          ic1   = IA(n+1)-1                                             
          DO m = ic0, ic1                                               
            icol = JA(m)                                                
            tv  = tv + A(m) * D1(icol)                                  
          END DO                                                        
          D2(n) = tv                                                    
        END DO                                                          
!---------RETURN                                                        
        RETURN                                                          
      END SUBROUTINE IMSLINEARSUB_MV  
      
      SUBROUTINE IMSLINEARSUB_AXPY(NEQ, D1, DC, D2, DR)
        IMPLICIT NONE
!     + + + DUMMY ARGUMENTS + + +
        integer(I4B), INTENT(IN) :: NEQ
        real(DP), DIMENSION(NEQ), INTENT(IN)    :: D1
        real(DP), INTENT(IN) :: DC
        real(DP), DIMENSION(NEQ), INTENT(IN)    :: D2
        real(DP), DIMENSION(NEQ), INTENT(INOUT) :: DR
!     + + + LOCAL DEFINITIONS + + +
        integer(I4B) :: n
!     + + + FUNCTIONS + + +
!     + + + CODE + + +
         DO n = 1, NEQ
          DR(n) = D1(n) + DC * D2(n)
         END DO
!---------RETURN
        RETURN
      END SUBROUTINE IMSLINEARSUB_AXPY

      
    FUNCTION IMSLINEARSUB_DP(neq, a, b) RESULT(c)
      ! -- return variable
      real(DP) :: c
!     + + + dummy arguments + + +
      integer(I4B), intent(in) :: neq
      real(DP), dimension(neq),  intent(in) :: a
      real(DP), dimension(neq),  intent(in) :: b
!     + + + local definitions + + +
      integer(I4B) :: n
!     + + + parameters + + +
!     + + + functions + + +
!     + + + code + + +
      c = DZERO
      do n = 1, neq
        c = c + a(n) * b(n)
      end do
      !---------return
      return
    END FUNCTION IMSLINEARSUB_DP

      
    FUNCTION IMSLINEARSUB_RNRM2(neq, a) RESULT(c)
      ! -- return variable
      real(DP) :: c
!     + + + dummy arguments + + +
      integer(I4B), intent(in) :: neq
      real(DP), dimension(neq),  intent(in) :: a
!     + + + local definitions + + +
      integer(I4B) :: n
      real(DP) :: ssq
      real(DP) :: scale
      real(DP) :: norm
      real(DP) :: absan
!     + + + parameters + + +
!     + + + functions + + +
!     + + + code + + +
      if (neq < 1) then
        norm = DZERO
      else if (neq == 1) then
        norm = ABS(a(1))
      else
        scale = DZERO
        ssq = DONE
        do n = 1, neq
          if (a(n) /= DZERO) then
            absan = abs(a(n))
            if (scale < absan) then
              ssq = DONE + ssq * (scale/absan)**2
              scale = absan
            else
              ssq = ssq + (absan/scale)**2
            end if
          end if
        end do
        norm = scale * sqrt(ssq)
      END IF
      c = norm
      !---------return
      return
    END FUNCTION IMSLINEARSUB_RNRM2
!
!    
!-------BEGINNING OF SUBROUTINES FROM OTHER LIBRARIES                   
                                                                        
      !       SUBSET OF SPARSKIT VERSION 2 SOURCE CODE
      !
      !  SPARSKIT VERSION 2 SUBROUTINES INCLUDED INCLUDE:
      !
      !    1 - IMSLINEARSUB_PCMILUT
      !    2 - IMSLINEARSUB_PCMILUT_LUSOL
      !    3 - IMSLINEARSUB_PCMILUT_QSPLIT
      !
      !-----------------------------------------------------------------------
      !                   S P A R S K I T   V E R S I O N  2.
      !-----------------------------------------------------------------------
      !
      !Latest update : Tue Mar  8 11:01:12 CST 2005
      !
      !-----------------------------------------------------------------------
      !
      !Welcome  to SPARSKIT  VERSION  2.  SPARSKIT is  a  package of  FORTRAN
      !subroutines  for working  with  sparse matrices.  It includes  general
      !sparse  matrix  manipulation  routines  as  well as  a  few  iterative
      !solvers, see detailed description of contents below.
      !
      ! Copyright (C) 2005, the Regents of the University of Minnesota
      !
      !SPARSKIT is  free software; you  can redistribute it and/or  modify it
      !under the terms of the  GNU Lesser General Public License as published
      !by the  Free Software Foundation [version  2.1 of the  License, or any
      !later version.]
      !
      !A copy of  the licencing agreement is attached in  the file LGPL.  For
      !additional information  contact the Free Software  Foundation Inc., 59
      !Temple Place - Suite 330, Boston, MA 02111, USA or visit the web-site
      !
      ! http://www.gnu.org/copyleft/lesser.html
      !
      !
      !DISCLAIMER
      !----------
      !
      !SPARSKIT  is distributed  in  the hope  that  it will  be useful,  but
      !WITHOUT   ANY  WARRANTY;   without  even   the  implied   warranty  of
      !MERCHANTABILITY  or FITNESS  FOR A  PARTICULAR PURPOSE.   See  the GNU
      !Lesser General Public License for more details.
      !
      !For more information contact saad@cs.umn.edu
      !
      !

      SUBROUTINE IMSLINEARSUB_PCMILUT(n, a, ja, ia, lfil, droptol, relax,       &
                                      alu, jlu, ju, iwk, w, jw, ierr,           &
                                      izero, delta)
      !-----------------------------------------------------------------------
        integer(I4B) :: n
        real(DP) :: a(*),alu(*),w(n+1),droptol,relax
        integer(I4B) :: ja(*),ia(n+1),jlu(*),ju(n),jw(2*n),lfil,iwk,ierr
        integer(I4B) :: izero
        real(DP) :: delta
        !----------------------------------------------------------------------*
        !                      *** ILUT preconditioner ***                     *
        !      incomplete LU factorization with dual truncation mechanism      *
        !----------------------------------------------------------------------*
        !     Author: Yousef Saad *May, 5, 1990, Latest revision, August 1996  *
        !----------------------------------------------------------------------*
        ! PARAMETERS
        !-----------
        !
        ! on entry:
        !==========
        ! n       = integer. The row dimension of the matrix A. The matrix
        !
        ! a,ja,ia = matrix stored in Compressed Sparse Row format.
        !
        ! lfil    = integer. The fill-in parameter. Each row of L and each row
        !           of U will have a maximum of lfil elements (excluding the
        !           diagonal element). lfil must be .ge. 0.
        !           ** WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
        !           EARLIER VERSIONS.
        !
        ! droptol = real. Sets the threshold for dropping small terms 
        !           in the factorization. See below for details on dropping 
        !           strategy.
        !
        !
        ! iwk     = integer. The lengths of arrays alu and jlu. If the arrays
        !           are not big enough to store the ILU factorizations, ilut
        !           will stop with an error message.
        !
        ! On return:
        !===========
        !
        ! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
        !           the L and U factors together. The diagonal (stored in
        !           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
        !           contains the i-th row of L (excluding the diagonal entry=1)
        !           followed by the i-th row of U.
        !
        ! ju      = integer array of length n containing the pointers to
        !           the beginning of each row of U in the matrix alu,jlu.
        !
        ! ierr    = integer. Error message with the following meaning.
        !           ierr  = 0    --> successful return.
        !           ierr .gt. 0  --> zero pivot encountered at step number ierr.
        !           ierr  = -1   --> Error. input matrix may be wrong.
        !                            (The elimination process has generated a
        !                            row in L or U whose length is .gt.  n.)
        !           ierr  = -2   --> The matrix L overflows the array al.
        !           ierr  = -3   --> The matrix U overflows the array alu.
        !           ierr  = -4   --> Illegal value for lfil.
        !           ierr  = -5   --> zero row encountered.
        !
        ! work arrays:
        !=============
        ! jw      = integer work array of length 2*n.
        ! w       = real work array of length n+1.
        !
        !----------------------------------------------------------------------
        ! w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u]
        ! jw(n+1:2n)  stores nonzero indicators
        !
        ! Notes:
        ! ------
        ! The diagonal elements of the input matrix must be  nonzero (at least
        ! 'structurally').
        !
        !----------------------------------------------------------------------*
        !---- Dual drop strategy works as follows.                             *
        !                                                                      *
        !     1) Thresholding in L and U as set by droptol. Any element whose  *
        !        magnitude is less than some tolerance (relative to the abs    *
        !        value of diagonal element in u) is dropped.                   *
        !                                                                      *
        !     2) Keeping only the largest lfil elements in the i-th row of L   *
        !        and the largest lfil elements in the i-th row of U (excluding *
        !        diagonal elements).                                           *
        !                                                                      *
        ! Flexibility: one  can use  droptol=0  to get  a strategy  based on   *
        ! keeping  the largest  elements in  each row  of L  and U.   Taking   *
        ! droptol .ne.  0 but lfil=n will give  the usual threshold strategy   *
        ! (however, fill-in is then unpredictable).                            *
        !----------------------------------------------------------------------*
        !     locals
        character(len=LINELENGTH) :: line
        integer(I4B) :: ju0,k,j1,j2,j,ii,i,lenl,lenu,jj,jrow,jpos,ilen
        real(DP) :: tnorm, t, abs, s, fact
        real(DP) :: rs, d, sd1, tl
        !     format
        character(len=*), parameter :: fmterr = "(//,1x,a)"
        !     code
        if (lfil .lt. 0) goto 998
        !-----------------------------------------------------------------------
        !     initialize ju0 (points to next element to be added to alu,jlu)
        !     and pointer array.
        !-----------------------------------------------------------------------
        ju0 = n+2
        jlu(1) = ju0
        !
        !     initialize nonzero indicator array.
        !
        do j = 1, n
          jw(n+j)  = 0
        end do
        !-----------------------------------------------------------------------
        !     beginning of main loop.
        !-----------------------------------------------------------------------
        main: do ii = 1, n
          j1 = ia(ii)
          j2 = ia(ii+1) - 1
          rs = DZERO
          tnorm = DZERO
          do k = j1, j2
            tnorm = tnorm+abs(a(k))
          end do
          if (tnorm .eq. DZERO) goto 999
          tnorm = tnorm/real(j2-j1+1)
          !
          !     unpack L-part and U-part of row of A in arrays w
          !
          lenu = 1
          lenl = 0
          jw(ii) = ii
          w(ii) = DZERO
          jw(n+ii) = ii
          !
          do j = j1, j2
            k = ja(j)
            t = a(j)
            if (k .lt. ii) then
              lenl = lenl+1
              jw(lenl) = k
              w(lenl) = t
              jw(n+k) = lenl
            else if (k .eq. ii) then
              w(ii) = t
            else
              lenu = lenu+1
              jpos = ii+lenu-1
              jw(jpos) = k
              w(jpos) = t
              jw(n+k) = jpos
            end if
          end do
          jj = 0
          ilen = 0
          !
          !     eliminate previous rows
          !
150       jj = jj+1
          if (jj .gt. lenl) goto 160
          !-----------------------------------------------------------------------
          !     in order to do the elimination in the correct order we must select
          !     the smallest column index among jw(k), k=jj+1, ..., lenl.
          !-----------------------------------------------------------------------
          jrow = jw(jj)
          k = jj
          !
          !     determine smallest column index
          !
          do j = jj+1, lenl
            if (jw(j) .lt. jrow) then
              jrow = jw(j)
              k = j
            end if
          end do
          !
          if (k .ne. jj) then
            !     exchange in jw
            j = jw(jj)
            jw(jj) = jw(k)
            jw(k) = j
            !     exchange in jr
            jw(n+jrow) = jj
            jw(n+j) = k
            !     exchange in w
            s = w(jj)
            w(jj) = w(k)
            w(k) = s
          end if
          !
          !     zero out element in row by setting jw(n+jrow) to zero.
          !
          jw(n+jrow) = 0
          !
          !     get the multiplier for row to be eliminated (jrow).
          !
          fact = w(jj)*alu(jrow)
          if (abs(fact) .le. droptol) then
            rs = rs + w(jj)
            goto 150
          end if
          !
          !     combine current row and row jrow
          !
          do k = ju(jrow), jlu(jrow+1)-1
            s = fact*alu(k)
            j = jlu(k)
            jpos = jw(n+j)
            if (j .ge. ii) then
              !
              !     dealing with upper part.
              !
              if (jpos .eq. 0) then
                !
                !     this is a fill-in element
                !
                lenu = lenu+1
                if (lenu .gt. n) goto 995
                i = ii+lenu-1
                jw(i) = j
                jw(n+j) = i
                w(i) = - s
              else
                !
                !     this is not a fill-in element
                !
                w(jpos) = w(jpos) - s

              end if
            else
              !
              !     dealing  with lower part.
              !
              if (jpos .eq. 0) then
                !
                !     this is a fill-in element
                !
                lenl = lenl+1
                if (lenl .gt. n) goto 995
                jw(lenl) = j
                jw(n+j) = lenl
                w(lenl) = - s
              else
                !
                !     this is not a fill-in element
                !
                w(jpos) = w(jpos) - s
              end if
            end if
          end do
          !
          !     store this pivot element -- (from left to right -- no danger of
          !     overlap with the working elements in L (pivots).
          !
          ilen = ilen+1
          w(ilen) = fact
          jw(ilen) = jrow
          goto 150
160       continue
          !
          !     reset double-pointer to zero (U-part)
          !
          do k = 1, lenu
            jw(n+jw(ii+k-1)) = 0
          end do
          !
          !     update L-matrix
          !
          lenl = ilen
          ilen = min0(lenl,lfil)
          !
          !     sort by quick-split
          !
          call IMSLINEARSUB_PCMILUT_QSPLIT(lenl, w, jw, ilen)
          !
          !     store L-part
          !
          do k = 1, ilen
            !            if (ju0 .gt. iwk) goto 996
            if (ju0 .gt. iwk) then
              write(line, '(2i10)') ju0, iwk
              call sim_message(line, fmt=fmterr, level=VDEBUG)
              goto 996
            end if
            alu(ju0) =  w(k)
            jlu(ju0) =  jw(k)
            ju0 = ju0+1
          end do
          !
          !     save pointer to beginning of row ii of U
          !
          ju(ii) = ju0
          !
          !     update U-matrix -- first apply dropping strategy
          !
          ilen = 0
          do k = 1, lenu-1
            if (abs(w(ii+k)) .gt. droptol*tnorm) then
              ilen = ilen+1
              w(ii+ilen) = w(ii+k)
              jw(ii+ilen) = jw(ii+k)
            else
              rs = rs + w(ii+k)
            end if
          end do
          lenu = ilen+1
          ilen = min0(lenu,lfil)
          !
          call IMSLINEARSUB_PCMILUT_QSPLIT(lenu-1, w(ii+1), jw(ii+1), ilen)
          !
          !     copy
          !
          t = abs(w(ii))
          !         if (ilen + ju0 .gt. iwk) goto 997
          if (ilen + ju0 .gt. iwk) then
            write(line, '(2i10)') (ilen + ju0), iwk
            call sim_message(line, fmt=fmterr, level=VDEBUG)
            goto 997
          end if
          do k = ii+1, ii+ilen-1
            jlu(ju0) = jw(k)
            alu(ju0) = w(k)
            t = t + abs(w(k) )
            ju0 = ju0+1
          end do
          !!
          !!     add dropped terms to diagonal element
          !!
          !IF (relax > DZERO) THEN
          !  w(ii) = w(ii) + relax * rs
          !END IF
          !!
          !!     store inverse of diagonal element of u
          !!
          !if (w(ii) == DZERO) w(ii) = (DEM4 + droptol)*tnorm
          !!
          !alu(ii) = DONE / w(ii)

          !    diagonal - calculate inverse of diagonal for solution       
          d   = w(ii) 
          tl  = ( DONE + delta ) * d + ( relax * rs ) 

          !    ensure that the sign of the diagonal has not changed
          sd1 = SIGN(d,tl) 
          IF (sd1.NE.d) THEN 
            !  use small value if diagonal scaling is not effective for 
            !    pivots that change the sign of the diagonal               
            IF (izero > 1) THEN 
              tl = SIGN(DONE,d) * (DEM4 + droptol) * tnorm 
            !  diagonal scaling continues to be effective                
            ELSE 
              izero = 1 
              exit main 
            END IF 
          END IF
          !    ensure that the diagonal is not zero
          IF (ABS(tl) == DZERO) THEN 
            !  use small value if diagonal scaling is not effective
            !    zero pivots                                               
            IF (izero > 1) THEN 
              tl = SIGN(DONE,d) * (DEM4 + droptol) * tnorm 
            !  diagonal scaling continues to be effective 
            ELSE 
              izero = 1 
              exit main 
            END IF 
          END IF
          w(ii) = tl
          alu(ii) = DONE / w(ii)           
          !
          !     update pointer to beginning of next row of U.
          !
          jlu(ii+1) = ju0
          !-----------------------------------------------------------------------
          !     end main loop
          !-----------------------------------------------------------------------
        end do main
        ierr = 0
        return
        !
        !     incomprehensible error. Matrix must be wrong.
        !
995     ierr = -1
        return
        !
        !     insufficient storage in L.
        !
996     ierr = -2
        return
        !
        !     insufficient storage in U.
        !
997     ierr = -3
        return
        !
        !     illegal lfil entered.
        !
998     ierr = -4
        return
        !
        !     zero row encountered
        !
999     ierr = -5
        return
      !----------------end-of-ilut--------------------------------------------
      !-----------------------------------------------------------------------
      END SUBROUTINE IMSLINEARSUB_PCMILUT

      !-----------------------------------------------------------------------
      SUBROUTINE IMSLINEARSUB_PCMILUT_LUSOL(n, y, x, alu, jlu, ju)
        integer(I4B) :: n
        real(DP) :: x(n), y(n), alu(*)
        integer(I4B) :: jlu(*), ju(*)
        !-----------------------------------------------------------------------
        !
        ! This routine solves the system (LU) x = y,
        ! given an LU decomposition of a matrix stored in (alu, jlu, ju)
        ! modified sparse row format
        !
        !-----------------------------------------------------------------------
        ! on entry:
        ! n   = dimension of system
        ! y   = the right-hand-side vector
        ! alu, jlu, ju
        !     = the LU matrix as provided from the ILU routines.
        !
        ! on return
        ! x   = solution of LU x = y.
        !-----------------------------------------------------------------------
        !
        ! Note: routine is in place: call IMSLINEARSUB_PCMILUT_LUSOL (n, x, x, alu, jlu, ju)
        !       will solve the system with rhs x and overwrite the result on x .
        !
        !-----------------------------------------------------------------------
        ! -- local
        !
        integer(I4B) :: i, k
        !
        ! forward solve
        !
        do i = 1, n
          x(i) = y(i)
          do k = jlu(i), ju(i)-1
            x(i) = x(i) - alu(k)* x(jlu(k))
          end do
        end do
        !
        !     backward solve.
        !
        do i = n, 1, -1
          do k = ju(i), jlu(i+1)-1
            x(i) = x(i) - alu(k)*x(jlu(k))
          end do
          x(i) = alu(i)*x(i)
        end do
        !
        return
      !----------------end of IMSLINEARSUB_PCMILUT_LUSOL ------------------------------------------
      !-----------------------------------------------------------------------
      END SUBROUTINE IMSLINEARSUB_PCMILUT_LUSOL
                                           
      subroutine imslinearsub_pcc_sol(this, neq, nlblk, lblkeq, r, s) !CGC
!       ******************************************************************
!       APPLY COARSE GRID CORRECTION
!       ******************************************************************
!    
!          SPECIFICATIONS:
!       ------------------------------------------------------------------
        use SimModule, only: ustop
!       + + + dummy variables + + +
        class(ImsLinearDataType), intent(inout) :: this
        integer(I4B), intent(in) :: neq
        integer(I4B), intent(in) :: nlblk
        integer(I4B), dimension(:), intent(in) :: lblkeq
        real(DP), dimension(neq), intent(in) :: r
        real(DP), dimension(neq), intent(inout) :: s
!       + + + local variables + + +
        integer(I4B) :: irbl, irbg, irbc, ieq, ireq0, ireq1
        real(DP) :: tv
!--  ----------------------------------------------------------------
        !
        !-- construct RHS
        do irbl = 1, nlblk ! my local blocks
          irbg = this%cgcl2gid(irbl)
          irbc = this%cgcg2cid(irbg) ! row block
          !
          ireq0 = lblkeq(irbl)
          ireq1 = lblkeq(irbl+1)-1
          !
          tv = DZERO
          do ieq = ireq0, ireq1
            tv = tv + this%zc(ieq)*r(ieq)
          end do
          this%wc(irbl) = tv
        end do
        !
        ! -- MPI all-to-all for global RHS
        call this%MpiSol%mpi_global_exchange_all_cgc(nlblk, this%cgcgnia-1,  &
          this%wc, this%rhsc, this%cgcgnia, this%cgcgnja, this%cgcgia, this%cgcgja, &
          2)
        !
        !-- solve the coarse system
        this%xc = this%rhsc
        if (this%cgcsol == 1) then
          !-- Band LU solve of the coarse system
          call imslinearsub_slu_solve(this%acpc, this%cgcgneq, this%acpcb,       &
                                      this%acpcb, this%xc)
        else
          ! -- ILU(0) solve of the coarse system
          call imslinearsub_ilu0a(this%cgcgnja, this%cgcgneq, this%acpci,        &
                                  this%cgcgiapc, this%cgcgjapc,                  &
                                  this%xc, this%xc)
        end if
        !
        !-- Correct search direction s
        do irbl = 1, nlblk ! my local blocks
          irbg = this%cgcl2gid(irbl)
          !
          ireq0 = lblkeq(irbl)
          ireq1 = lblkeq(irbl+1)-1
          !
          do ieq = ireq0, ireq1
            s(ieq) = s(ieq) + this%xc(irbg)*this%zc(ieq)
          end do
        end do
        !
        !-- return
        return
      end subroutine imslinearsub_pcc_sol
      
      !-----------------------------------------------------------------------
      SUBROUTINE IMSLINEARSUB_PCMILUT_QSPLIT(n, a, ind, ncut)
        integer(I4B) :: n
        real(DP) :: a(n)
        integer(I4B) :: ind(n), ncut
        !-----------------------------------------------------------------------
        !     does a quick-sort split of a real array.
        !     on input a(1:n). is a real array
        !     on output a(1:n) is permuted such that its elements satisfy:
        !
        !     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
        !     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
        !
        !     ind(1:n) is an integer array which permuted in the same way as a(*
        !-----------------------------------------------------------------------
        real(DP) :: tmp, abskey
        integer(I4B) :: itmp, first, last
        integer(I4B) :: mid
        integer(I4B) :: j
        !-----
        first = 1
        last = n
        if (ncut .lt. first .or. ncut .gt. last) return
        !
        !     outer loop -- while mid .ne. ncut do
        !
00001   mid = first
        abskey = abs(a(mid))
        do j = first+1, last
          if (abs(a(j)) .gt. abskey) then
            mid = mid+1
            !     interchange
            tmp = a(mid)
            itmp = ind(mid)
            a(mid) = a(j)
            ind(mid) = ind(j)
            a(j)  = tmp
            ind(j) = itmp
          end if
        end do
        !
        !     interchange
        !
        tmp = a(mid)
        a(mid) = a(first)
        a(first)  = tmp
        !
        itmp = ind(mid)
        ind(mid) = ind(first)
        ind(first) = itmp
        !
        !     test for while loop
        !
        if (mid .eq. ncut) return
        if (mid .gt. ncut) then
          last = mid-1
        else
          first = mid+1
        end if
        goto 1
        !----------------end-of-IMSLINEARSUB_PCMILUT_QSPLIT------------------------------------------
        !-----------------------------------------------------------------------
      END SUBROUTINE IMSLINEARSUB_PCMILUT_QSPLIT

END MODULE IMSLinearModule
