module Mf6CoreModule
  use KindModule,             only: I4B
  use ListsModule,            only: basesolutionlist, solutiongrouplist, basemodellist, baseexchangelist
  use BaseModelModule,        only: BaseModelType, GetBaseModelFromList
  use BaseExchangeModule,     only: BaseExchangeType, GetBaseExchangeFromList
  use BaseSolutionModule,     only: BaseSolutionType, GetBaseSolutionFromList
  use SolutionGroupModule,    only: SolutionGroupType, GetSolutionGroupFromList
  use NumericalSolutionModule, only: NumericalSolutionType !PAR
  use SimVariablesModule,      only: isimdd !PAR
  implicit none
contains
  
  subroutine Mf6Run
  ! ******************************************************************************
  ! Main MODFLOW Version 6 program.
  ! ******************************************************************************
  !
  !    SPECIFICATIONS:
  ! ------------------------------------------------------------------------------
    ! -- modules 
    
    use CommandArguments, only: GetCommandLineArguments
    use TdisModule, only: totim, totalsimtime
    use KindModule, only: DP
    use MpiExchangeModule, only: mpi_initialize_world !PAR
    logical :: hasConverged
    !
    ! -- Initialize MPI if required
    call mpi_initialize_world() !PAR
    !
    ! -- parse any command line arguments
    call GetCommandLineArguments()
    !
    ! initialize simulation
    call Mf6Initialize()
    !
    ! -- time loop
    tsloop: do while (totim < totalsimtime)
      
      ! perform a time step
      hasConverged = Mf6Update()
      
      ! if not converged, break
      if(.not. hasConverged) exit tsloop   
      
    enddo tsloop
    !
    ! -- finalize simulation
    call Mf6Finalize()
    
  end subroutine Mf6Run
  
  subroutine Mf6Initialize()
    use SimulationCreateModule, only: simulation_cr
    
    ! -- print banner and info to screen
    call printInfo()
    
    ! -- create
    call simulation_cr()
    
    ! -- define
    call simulation_df()
       
    ! -- allocate and read
    call simulation_ar()
    
  end subroutine Mf6Initialize
  
  function Mf6Update() result(hasConverged)
    logical :: hasConverged
    !
    ! -- prepare timestep
    call Mf6PrepareTimestep()
    !
    ! -- do timestep
    call Mf6DoTimestep()      
    !
    ! -- after timestep
    hasConverged = Mf6FinalizeTimestep()
    !
  end function Mf6Update
  
  subroutine Mf6Finalize()
    use MpiExchangeGenModule, only: parallelrun, writestd !PAR
    use MpiExchangeModule, only: MpiWorld !PAR
    use, intrinsic :: iso_fortran_env, only: output_unit
    use ListsModule,            only: lists_da
    use MemoryManagerModule,    only: mem_usage, mem_da,                         &
                                      mem_timing,                                & !TIM
                                      bytes !PAR
    use TimerModule,            only: elapsed_time
    use SimVariablesModule,     only: iout
    use SimulationCreateModule, only: simulation_cr, simulation_da  
    use TdisModule,             only: tdis_tu, tdis_da
    use SimModule,              only: final_message
    use MpiExchangeModule,      only: mpi_world_da !PAR
    use MpiExchangeGenModule,   only: writestd !PAR

    integer(I4B) :: im, ic, is, isg
    class(SolutionGroupType), pointer :: sgp => null()
    class(BaseSolutionType), pointer :: sp => null()
    class(BaseModelType), pointer :: mp => null()
    class(BaseExchangeType), pointer :: ep => null()
    
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
    ! -- Calculate memory usage, elapsed time and terminate
    call mem_usage(iout)
    call MpiWorld%mpi_total_memory(bytes) !PAR
    !
    call mem_timing(iout) !TIM
    call MpiWorld%mpi_barrier() !PAR
    call elapsed_time(iout, 1, writestd) !PAR
    !
    ! -- Deallocate MPI world
    call mpi_world_da() !PAR
    !
    call mem_da()
    !
    call final_message()
    !        
  end subroutine Mf6Finalize
  
  subroutine printInfo()         
    use CompilerVersion
    use VersionModule,          only: VERSION, MFVNAM, MFTITLE, FMTDISCLAIMER,   &
                                      IDEVELOPMODE
    use TimerModule,            only: start_time    
    use GenericUtilitiesModule, only: write_centered, sim_message
    use MpiExchangeGenModule,   only: serialrun, writestd !PAR
    use MpiExchangeModule,      only: MpiWorld !PAR
    use SimVariablesModule,     only: istdout !PAR
    character(len=80) :: compiler
    
    ! -- Write banner to screen (unit stdout) and start timer
    call write_centered('MODFLOW'//MFVNAM, 80)
    call write_centered(MFTITLE, 80)
    call write_centered('VERSION '//VERSION, 80)
    if (.not.serialrun) call write_centered('***RUNNING IN PARALLEL MODE WITH '  & !PAR
      //TRIM(MpiWorld%nrprocstr)//' MPI PROCESSES***', 80) !PAR
    !
    ! -- Write if develop mode
    if (IDEVELOPMODE == 1) then
      call write_centered('***DEVELOP MODE***', 80)
    end if
    !
    ! -- Write compiler version
    call get_compiler(compiler)
    call write_centered(' ', 80)
    call write_centered(trim(adjustl(compiler)), 80)
    !
    ! -- Write disclaimer
    if (writestd) call sim_message('', fmt=FMTDISCLAIMER) !PAR
    ! -- get start time
    call MpiWorld%mpi_barrier() !PAR
    call start_time(writestd) !PAR

  end subroutine printInfo
  
  subroutine simulation_df()
    use MpiExchangeGwfModule, only: mpi_gwfhalo_world,                           & !PAR
                                    mpi_set_gwfhalo_world_mvr !PAR
    use GwfGwfExchangeModule, only: gwf_mpi_halo_init !PAR
    use MpiExchangeModule,    only: MpiWorld !PAR
    integer(I4B) :: im, ic, is
    class(BaseSolutionType), pointer :: sp => null()
    class(BaseModelType), pointer :: mp => null()
    class(BaseExchangeType), pointer :: ep => null()

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
    do ic = 1, baseexchangelist%Count() !HALO2
      ep => GetBaseExchangeFromList(baseexchangelist, ic) !HALO2
      call ep%exg_df(1) !HALO2
    enddo !HALO2
    !
    ! -- MPI local exchange
    if (isimdd == 1) then !PAR
      call gwf_mpi_halo_init()
      call MpiWorld%mpi_local_exchange('', 'HALO_INIT_DIS', .false.)
    endif
    !
    ! -- Define each exchange phase 2
    do ic = 1, baseexchangelist%Count() !HALO2
      ep => GetBaseExchangeFromList(baseexchangelist, ic) !HALO2
      call ep%exg_df(2) !HALO2
    enddo !HALO2
    !
    ! -- Define each solution
    do is = 1, basesolutionlist%Count()
      sp => GetBaseSolutionFromList(basesolutionlist, is)
      call sp%sln_df()
    enddo
  
  end subroutine simulation_df
  
  subroutine simulation_ar()
    use MpiExchangeModule, only: MpiWorld !PAR
    use MpiExchangeGwfModule, only: mpi_gwfhalo_world !PAR
    integer(I4B) :: im, ic, is
    class(BaseSolutionType), pointer :: sp => null()
    class(BaseModelType), pointer :: mp => null()
    class(BaseExchangeType), pointer :: ep => null()
    class(NumericalSolutionType), pointer :: nsp => null() !PAR

    ! -- Allocate and read each model
    do im = 1, basemodellist%Count()
      mp => GetBaseModelFromList(basemodellist, im)
      call mp%model_ar()
    enddo
    !
    if (isimdd == 1) then !PAR
      ! -- Local MPI exchange of NPF flags
      call mpi_gwfhalo_world(3)
      !
      !-- Local MPI exchange of NPF arrays
      call MpiWorld%mpi_local_exchange('', 'HALO_INIT_NPFIC', .false.) !PAR
      !
      do is=1,basesolutionlist%Count() !PAR
        sp => GetBaseSolutionFromList(basesolutionlist, is) !PAR
        select type (sp) !PAR
        class is (NumericalSolutionType) !PAR
          nsp => sp !PAR
        end select !PAR
        call nsp%MpiSol%mpi_local_exchange(nsp%name, 'X_IACTIVE', .false.) !PAR
      enddo !PAR
    endif !PAR

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
  end subroutine simulation_ar
  
  subroutine Mf6PrepareTimestep()
    use KindModule,             only: I4B
    use TdisModule,             only: tdis_tu
    use ListsModule,            only: basesolutionlist, basemodellist, baseexchangelist
    use BaseModelModule,        only: BaseModelType, GetBaseModelFromList
    use BaseExchangeModule,     only: BaseExchangeType, GetBaseExchangeFromList
    use BaseSolutionModule,     only: BaseSolutionType, GetBaseSolutionFromList
    use SimModule,              only: converge_reset
    
    integer(I4B) :: im, ic, is
    class(BaseSolutionType), pointer :: sp => null()
    class(BaseModelType), pointer :: mp => null()
    class(BaseExchangeType), pointer :: ep => null()
    class(NumericalSolutionType), pointer :: nsp !PAR
    
    ! -- time update
    call tdis_tu()
    
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
    call converge_reset()
    
  end subroutine Mf6PrepareTimestep
  
  subroutine Mf6DoTimestep()
    use KindModule,           only: I4B
    use ListsModule,          only: solutiongrouplist
    use SolutionGroupModule,  only: SolutionGroupType, GetSolutionGroupFromList
    class(SolutionGroupType), pointer :: sgp => null()
    integer(I4B) :: isg
    
    do isg = 1, solutiongrouplist%Count()
      sgp => GetSolutionGroupFromList(solutiongrouplist, isg)
      call sgp%sgp_ca()
    enddo
      
  end subroutine Mf6DoTimestep
  
  function Mf6FinalizeTimestep() result(hasConverged)
    use KindModule,             only: I4B
    use ListsModule,            only: basesolutionlist, basemodellist, baseexchangelist    
    use BaseModelModule,        only: BaseModelType, GetBaseModelFromList
    use BaseExchangeModule,     only: BaseExchangeType, GetBaseExchangeFromList
    use BaseSolutionModule,     only: BaseSolutionType, GetBaseSolutionFromList
    use SimModule,              only: converge_check
    logical :: hasConverged    
    integer(I4B) :: im, ic, is
    class(BaseSolutionType), pointer :: sp => null()
    class(BaseModelType), pointer :: mp => null()
    class(BaseExchangeType), pointer :: ep => null()
    
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
    ! -- Check if we're done
    call converge_check(hasConverged)
    
  end function Mf6FinalizeTimestep
  
end module Mf6CoreModule
