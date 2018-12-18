!
program mf6
! ******************************************************************************
! Main MODFLOW Version 6 program.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  ! -- modules
  use KindModule,             only: DP, I4B
  use ConstantsModule,        only: ISTDOUT
  use VersionModule,          only: VERSION, MFVNAM, MFTITLE, FMTDISCLAIMER,   &
                                    IDEVELOPMODE
  use CompilerVersion
  use CommandArguments,       only: GetCommandLineArguments
  use InputOutputModule,      only: write_centered
  use SimulationCreateModule, only: simulation_cr, simulation_da
  use TimerModule,            only: start_time, elapsed_time
  use MemoryManagerModule,    only: mem_usage, mem_da, &
                                    mem_timing !TIM
  use BaseModelModule,        only: BaseModelType, GetBaseModelFromList
  use BaseExchangeModule,     only: BaseExchangeType, GetBaseExchangeFromList
  use BaseSolutionModule,     only: BaseSolutionType, GetBaseSolutionFromList
  use SolutionGroupModule,    only: SolutionGroupType, GetSolutionGroupFromList
  use ListsModule,            only: basesolutionlist, solutiongrouplist,       &
                                    basemodellist, baseexchangelist,           &
                                    lists_da
  use SimVariablesModule,     only: iout, isimdd !PAR
  use SimModule,              only: converge_reset, converge_check,            &
                                    final_message
  use TdisModule,             only: tdis_tu, tdis_da,                          &
                                    endofsimulation
  use MpiExchangeGenModule,   only: serialrun, writestd !PAR
  use MpiExchangeModule,      only: mpi_initialize_world, mpi_world_da,        & !PAR
                                    MpiWorld !PAR
  use MpiExchangeGwfModule,   only: mpi_gwfhalo_world,                        &
                                    mpi_set_gwfhalo_world_mvr !PAR
  use NumericalSolutionModule, only: NumericalSolutionType !PAR
  use GwfGwfExchangeModule, only: gwf_mpi_halo_init
  implicit none
  ! -- local
  class(SolutionGroupType), pointer :: sgp => null()
  class(BaseSolutionType),  pointer :: sp => null()
  class(NumericalSolutionType), pointer :: nsp !PAR
  class(BaseModelType),     pointer :: mp => null()
  class(BaseExchangeType),  pointer :: ep => null()
  integer(I4B) :: im, ic, is, isg, iprtim
  logical :: exit_tsloop
  character(len=80) :: compiler
  ! -- formats
! ------------------------------------------------------------------------------
  !
  ! -- Initialize MPI if required
  call mpi_initialize_world() !PAR
  !
  ! -- parse any command line arguments
  call GetCommandLineArguments()
  !
  ! -- Write banner to screen (unit 6) and start timer
  call write_centered('MODFLOW'//MFVNAM, ISTDOUT, 80)
  call write_centered(MFTITLE, ISTDOUT, 80)
  call write_centered('VERSION '//VERSION, ISTDOUT, 80)
  if (.not.serialrun) call write_centered('***RUNNING IN PARALLEL MODE WITH '& !PAR
    //TRIM(MpiWorld%nrprocstr)//' MPI PROCESSES***',ISTDOUT, 80) !PAR
  !
  ! -- Write if develop mode
  if (IDEVELOPMODE == 1) call write_centered('***DEVELOP MODE***', ISTDOUT, 80)
  !
  ! -- Write compiler version
  call get_compiler(compiler)
  call write_centered(' ', ISTDOUT, 80)
  call write_centered(trim(adjustl(compiler)), ISTDOUT, 80)
  !
  ! -- Write disclaimer
  if (writestd) write(ISTDOUT, FMTDISCLAIMER) !PAR
  ! -- get start time
  call MpiWorld%mpi_barrier() !PAR 
  call start_time(writestd) !PAR
  !
  !
  ! -- CREATE (CR)
  call simulation_cr()
  !
  !
  ! -- DEFINE (DF)
  ! -- Define each model
  do im = 1, basemodellist%Count()
    mp => GetBaseModelFromList(basemodellist, im)
    call mp%model_df()
  enddo
  !
  ! -- Collective MPI communication scalars DIS
  if (isimdd == 1) then !PAR
    call mpi_gwfhalo_world(2) !PAR
    call mpi_set_gwfhalo_world_mvr() !PAR
  endif !PAR
  !
  ! -- Define each exchange phase 1
  do ic = 1, baseexchangelist%Count()
    ep => GetBaseExchangeFromList(baseexchangelist, ic)
    call ep%exg_df(1)
  enddo
  !
  ! -- MPI local exchange 
  if (isimdd == 1) then !PAR
    call gwf_mpi_halo_init()
    call MpiWorld%mpi_local_exchange('', 'HALO_INIT_DIS', .false.)
  endif
  !  
  ! -- Define each exchange phase 2
  do ic = 1, baseexchangelist%Count()
    ep => GetBaseExchangeFromList(baseexchangelist, ic)
    call ep%exg_df(2)
  enddo
  !
  ! -- Define each solution
  do is = 1, basesolutionlist%Count()
    sp => GetBaseSolutionFromList(basesolutionlist, is)
    call sp%sln_df()
  enddo
  !
  !
  ! -- ALLOCATE AND READ (AR)
  ! -- Allocate and read each model
  do im = 1, basemodellist%Count()
    mp => GetBaseModelFromList(basemodellist, im)
    call mp%model_ar()
  enddo
  !
  if (isimdd == 1) then !PAR
    ! -- Local MPI exchange of NPF flags
    call mpi_gwfhalo_world(3)
    !-- Local MPI exchange of NPF arrays
    call MpiWorld%mpi_local_exchange('', 'HALO_INIT_NPFIC', .false.) !PAR
    do is=1,basesolutionlist%Count() !PAR
      sp => GetBaseSolutionFromList(basesolutionlist, is) !PAR
      select type (sp) !PAR
      class is (NumericalSolutionType) !PAR
        nsp => sp !PAR
      end select !PAR
      call nsp%MpiSol%mpi_local_exchange(nsp%name, 'X_IACTIVE', .false.) !PAR
    enddo !PAR
  endif !PAR
  !
  ! -- Allocate and read each exchange
  do ic = 1, baseexchangelist%Count()
    ep => GetBaseExchangeFromList(baseexchangelist, ic)
    call ep%exg_ar()
  enddo
  !
  ! -- Allocate and read each solution
  do is=1,basesolutionlist%Count()
    sp => GetBaseSolutionFromList(basesolutionlist, is)
    call sp%sln_ar()
  enddo
  !
  !
  ! -- TIME STEP LOOP
  tsloop: do
    !
    ! -- TIME UPDATE (TU)
    call tdis_tu()
    !
    !
    ! -- READ AND PREPARE (RP)
    ! -- Read and prepare each model
    do im = 1, basemodellist%Count()
      mp => GetBaseModelFromList(basemodellist, im)
      call mp%model_rp()
    enddo
    !
    ! -- Read and prepare each exchange
    do ic = 1, baseexchangelist%Count()
      ep => GetBaseExchangeFromList(baseexchangelist, ic)
      call ep%exg_rp()
    enddo
    !
    ! -- MPI parallel: initialize point-to-point mover
    if (isimdd == 1) then !PAR
      do is=1,basesolutionlist%Count() !PAR
        sp => GetBaseSolutionFromList(basesolutionlist, is) !PAR
        select type (sp) !PAR
        class is (NumericalSolutionType) !PAR
          nsp => sp !PAR
        end select !PAR
        call nsp%slnmpimvrinit(nsp%name) !PAR
      enddo !PAR
    endif !PAR
    !
    ! -- Read and prepare each solution
    do is=1,basesolutionlist%Count()
      sp => GetBaseSolutionFromList(basesolutionlist, is)
      call sp%sln_rp()
    enddo
    !
    !
    ! -- CALCULATE (CA)
    call converge_reset()
    do isg = 1, solutiongrouplist%Count()
      sgp => GetSolutionGroupFromList(solutiongrouplist, isg)
      call sgp%sgp_ca()
    enddo
    !
    !
    ! -- OUTPUT (OT)
    ! -- Write output for each model
    do im = 1, basemodellist%Count()
      mp => GetBaseModelFromList(basemodellist, im)
      call mp%model_ot()
    enddo
    !
    ! -- Write output for each exchange
    do ic = 1, baseexchangelist%Count()
      ep => GetBaseExchangeFromList(baseexchangelist, ic)
      call ep%exg_ot()
    enddo
    !
    ! -- Write output for each solution
    do is=1,basesolutionlist%Count()
      sp => GetBaseSolutionFromList(basesolutionlist, is)
      call sp%sln_ot()
    enddo
    !
    ! -- Time step exit conditions
    call converge_check(exit_tsloop)
    if(exit_tsloop) exit tsloop
    if(endofsimulation) exit tsloop
    !
  enddo tsloop
  !
  !
  ! -- FINAL PROCESSING (FP)
  ! -- Final processing for each model
  do im = 1, basemodellist%Count()
    mp => GetBaseModelFromList(basemodellist, im)
    call mp%model_fp()
  enddo
  !
  ! -- Final processing for each exchange
  do ic = 1, baseexchangelist%Count()
    ep => GetBaseExchangeFromList(baseexchangelist, ic)
    call ep%exg_fp()
  enddo
  !
  ! -- Final processing for each solution
  do is=1,basesolutionlist%Count()
    sp => GetBaseSolutionFromList(basesolutionlist, is)
    call sp%sln_fp()
  enddo
  !
  ! -- DEALLOCATE (DA)
  ! -- Deallocate tdis
  call tdis_da()
  !
  ! -- Deallocate for each model
  do im = 1, basemodellist%Count()
    mp => GetBaseModelFromList(basemodellist, im)
    call mp%model_da()
    deallocate(mp)
  enddo
  !
  ! -- Deallocate for each exchange
  do ic = 1, baseexchangelist%Count()
    ep => GetBaseExchangeFromList(baseexchangelist, ic)
    call ep%exg_da()
    deallocate(ep)
  enddo
  !
  ! -- Deallocate for each solution
  do is=1,basesolutionlist%Count()
    sp => GetBaseSolutionFromList(basesolutionlist, is)
    call sp%sln_da()
    deallocate(sp)
  enddo
  !
  ! -- Deallocate solution group and simulation variables
  do isg = 1, solutiongrouplist%Count()
    sgp => GetSolutionGroupFromList(solutiongrouplist, isg)
    call sgp%sgp_da()
    deallocate(sgp)
  enddo
  call simulation_da()
  call lists_da()
  !
  ! -- Deallocate MPI world
  call MpiWorld%mpi_barrier() !PAR
  call mpi_world_da() !PAR
  !
  ! -- Calculate memory usage, elapsed time and terminate
  call mem_usage(iout)
  call mem_timing(iout) !TIM
  call mem_da()
  call elapsed_time(iout, 1, writestd)
  call final_message()
  !
end program mf6

