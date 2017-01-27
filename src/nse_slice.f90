Program nse_slice
  Use controls
  Use nse_module
  Use nuclear_data
  Use timers
  Use, intrinsic :: iso_fortran_env, only: stdin=>input_unit, stdout=>output_unit
  Character (LEN=80) :: data_desc,data_dir
  Character (LEN=80) :: diag_file_base='nse_diag'
  Character (LEN=80) :: diag_file
  Character (LEN=80) :: timer_file_base='nse_timer'
  Character (LEN=80) :: timer_file
  Real(8) :: rho,ye,t9,rhoo,yeo
  Integer :: lun_control,ierr,knse,i
  Real(8) :: atot,ztot,ytot,abar,zbar,benuc,ubind

  Open(NEWUNIT=lun_control,FILE='control')

! Read Job Controls
  call find_controls_block(lun_control,'Job Controls',ierr)
  Read(lun_control,*) iscrn        ! controls the treatment of nuclear screening
  Read(lun_control,*) iprocess     ! controls the runtime pre-processing of the network data

! Read Output Controls
  call find_controls_block(lun_control,'Output Controls',ierr)
  Read(lun_control,*) idiag        ! sets diagnostic output level
  Read(lun_control,*) itsout       ! sets per timestep output level

! Read Input Controls
  call find_controls_block(lun_control,'Input Controls',ierr)
  Read(lun_control,"(72x)")
  Read(lun_control,"(a80)") data_dir

  Close(lun_control)

! Open diagnositic output file, per thread if OMP
  If(idiag>=0) Then
    diag_file=trim(diag_file_base)
    Open(NEWUNIT=lun_diag,file=diag_file)
    
    timer_file=trim(timer_file_base)
    Open(NEWUNIT=lun_timer,file=timer_file)
  EndIf

! In requested, pre-process the nuclear and reaction data.
  If(iprocess>0) Call net_preprocess(stdout,data_dir,data_dir)

! Initialize EoS for screening
  If(iscrn>0) Call eos_initialize

! Read nuclear dataset
  Call read_nuclear_data(data_dir,data_desc)

! Allocate and initialize NSE arrays
  Call nse_init

  knse = 0
  iconverge(:) = 0

  DO

    Write(stdout,'(a)') 'Rho, Ye, T9?'
    Read(stdin,*,iostat=ierr) rho, ye, t9
    IF(ierr/=0) Exit

    knse = knse + 1
    IF(knse>1) THEN
      IF(rhoo/=rho) THEN
        Write(lun_timer,*)
        Write(lun_diag,*)
      END IF
      IF(yeo/=ye) THEN
        Write(lun_timer,*)
        Write(lun_diag,*)
!       Write(stdout,"(31i5)") (iconverge(i),i=-15,15)
        iconverge(:) = 0
      END IF
    END IF

! Calculate NSE abundances
    Call nse_solve( rho, t9, ye )
    atot =sum(aa*ynse)
    ztot =sum(zz*ynse)
    ytot =sum(ynse)
    abar =atot/ytot
    zbar =ztot/ytot
    benuc=sum(be*ynse)
    ubind=-9.616221376e17*benuc

    IF ( itsout >= 0 ) THEN
      Write(stdout,'(a)')             'NSE composition written to '//TRIM(diag_file)
      Write(stdout,'(a,1x,3i15)')     'NR iters, LS iters, F evals =',knrtot(1),knrtot(2),knrtot(3)
      Write(stdout,'(a,2es23.15)')    'zbar, zbar                  =',zbar,abar
      Write(stdout,'(a,4es23.15)')    'xn, xp, alpha, xni56        =',xnse(i_nn),xnse(i_pp),xnse(i_he4),xnse(i_ni56)

      ii = MAXLOC( xnse(:), 1 )
      Write(stdout,'(a,a23,es23.15)') 'largest mass fraction       =',nname(ii),xnse(ii)
    END IF
    IF ( itsout >= 1 ) THEN
      Write(stdout,'(a)')             'NSE composition (mass fractions):' 
      Write(stdout,'(5(a6,1es10.3))') (nname(i),xnse(i),i=1,ny)
    END IF

    Write(lun_diag,'(a12,3es12.5,5es13.6,3i5)') 'NSE solved',rho,ye,t9,ynse(1:2),zbar,abar,-ubind,knrtot(1:3)
    Write(lun_diag,"(5(a6,1es10.3))") (nname(i),xnse(i),i=1,ny)

    Write(lun_timer,'(7es12.5,3i5)') rho,ye,t9,timer_nsesolv,timer_nsenrap,timer_nseeval,timer_nsescrn,knrtot(1:3)

    rhoo=rho
    yeo=ye

  END DO

! Write(stdout,"(31i5)") (iconverge(i),i=-15,15)

  Close(lun_timer)
  Close(lun_diag)

End Program nse_slice
