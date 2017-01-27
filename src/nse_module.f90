MODULE nse_module
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: stdin=>input_unit, stdout=>output_unit
  USE constants
  USE controls, ONLY: lun_diag, idiag, itsout, iscrn
  USE part_funct_data
  USE nuclear_data
  USE timers
  IMPLICIT NONE
  PUBLIC

  REAL(8), PARAMETER :: onethird  = 1.0d0/3.0d0
  REAL(8), PARAMETER :: twothird  = 2.0d0*onethird
  REAL(8), PARAMETER :: fivethird = 5.0d0*onethird

  INTEGER, PARAMETER :: lnorm    = 0
  REAL(8), PARAMETER :: tolu     = 1.0d-14
  REAL(8), PARAMETER :: tolf     = 1.0d-8
  REAL(8), PARAMETER :: tolmin   = 1.0d-14
  REAL(8), PARAMETER :: tolbeta  = 3.0d-15
  REAL(8), PARAMETER :: maxstep0 = 2.0d0
  REAL(8), PARAMETER :: typf     = 1.0d0

  REAL(8), PARAMETER :: alpha  = 1.0d-8
  REAL(8), PARAMETER :: gamma1 = 0.1d0
  REAL(8), PARAMETER :: gamma2 = 0.5d0

  INTEGER, PARAMETER :: nritmax = 200
  INTEGER, PARAMETER :: lsitmax = 30

  ! NSE variables
  REAL(8) :: rhonse, t9nse, yense
  REAL(8), ALLOCATABLE :: xnse(:), ynse(:), unse(:)
  REAL(8), ALLOCATABLE :: mm(:), mm52(:)
!$OMP THREADPRIVATE(rhonse,t9nse,yense,xnse,ynse,unse,mm,mm52)

  REAL(8) :: fvec(2), fjac(2,2)
!$OMP THREADPRIVATE(fvec,fjac)

  ! Scaling variables
  REAL(8) :: typu(2), typfvec(2)
  REAL(8) :: scaleu(2), scalefvec(2)
!$OMP THREADPRIVATE(typu,typfvec,scaleu,scalefvec)

  ! Counters
  INTEGER :: knrtot(3), knr(3)
  INTEGER :: itry, iconverge(-15:15)
  LOGICAL :: use_CP98
!$OMP THREADPRIVATE(knrtot,knr,itry,iconverge,use_CP98)

  ! Screening variables
  INTEGER :: zmax
  INTEGER, ALLOCATABLE :: intz(:)
  REAL(8), ALLOCATABLE :: hnse(:)
!$OMP THREADPRIVATE(zmax,intz,hnse)

  CONTAINS

  FUNCTION safe_exp( x ) RESULT( y )

    ! Input variables
    REAL(8), INTENT(IN) :: x(:)

    ! Function variable
    REAL(8) :: y(SIZE(x))

    ! Local variables
    REAL(8), PARAMETER :: exp_max = MAXEXPONENT(1.0d0)*LOG(2.0d0)*0.99d0
    REAL(8), PARAMETER :: exp_min = MINEXPONENT(1.0d0)*LOG(2.0d0)*0.99d0

    y(:) = EXP( MIN( exp_max, MAX( exp_min, x(:) ) ) )

    RETURN
  END FUNCTION safe_exp

  FUNCTION l2norm( x ) RESULT( y )

    ! Input variables
    REAL(8), INTENT(IN) :: x(:)

    ! Function variable
    REAL(8) :: y

    ! Local variables
    INTEGER :: n, incx

    ! External function
    REAL(8), EXTERNAL :: dnrm2

    n = SIZE(x)
    incx = 1

    y = dnrm2( n, x, incx )

    RETURN
  END FUNCTION l2norm

  FUNCTION dot_prod( x, y ) RESULT( z )

    ! Input variables
    REAL(8), INTENT(IN) :: x(:)
    REAL(8), INTENT(IN) :: y(:)

    ! Function variable
    REAL(8) :: z

    ! Local variables
    INTEGER :: n, incx, incy

    ! External function
    REAL(8), EXTERNAL :: ddot

    n = SIZE(x)
    incx = 1
    incy = 1

    z = ddot( n, x, incx, y, incy )

    RETURN
  END FUNCTION dot_prod

  SUBROUTINE nse_init

    ALLOCATE(xnse(ny),ynse(ny),unse(ny),mm(ny),mm52(ny))

    mm(:) = aa(:) / avn! + mex(:)*epmev/(clt*clt)
!   mm(:) = zz(:)*(mp+me) + nn(:)*mn - be(:)*epmev/(clt*clt)
    mm52(:) = mm(:) * mm(:) * SQRT(mm(:))

    typu(:)    = (/ 1.0d0, 1.0d0 /)
    typfvec(:) = (/ 1.0d0, 1.0d0 /)

    scaleu(:)    = 1.0d0 / typu(:)
    scalefvec(:) = 1.0d0 / typfvec(:)

    zmax = NINT( MAXVAL( zz(:), 1 ) )
    ALLOCATE(hnse(zmax),intz(ny))
    WHERE ( zz(:) > 1.0d0 )
      intz(:) = NINT(zz(:))
    ELSE WHERE
      intz(:) = 1
    END WHERE

    RETURN
  END SUBROUTINE nse_init

  SUBROUTINE nse_solve(rho,t9,ye)

    ! Input variables
    REAL(8), INTENT(IN) :: rho, t9, ye

    ! Local variables
    REAL(8) :: uvec(2), uvec0(2)
    INTEGER :: i, ii
    LOGICAL :: check
    INTEGER :: info, info0, itry0

    start_timer = xnet_wtime()
    timer_nsesolv = timer_nsesolv - start_timer

    ! Reset counters
    knrtot(:) = 0

    ! Initialize NSE variables
    rhonse = rho
    t9nse = t9
    yense = ye
    CALL nse_pf

    ! Start with easy screening to get a good guess
    use_CP98 = .true.
    CALL nse_screen_CP98
    DO itry = 1, 6
      unse(:) = 0.0d0
      xnse(:) = 0.0d0
      ynse(:) = 0.0d0
      CALL nse_guess(uvec)
      CALL nse_nr(uvec,check,info)
      IF ( info > 0 ) EXIT
    END DO
    IF ( itsout >=2 ) WRITE(stdout,'(a5,2i2,8es23.15)') 'CP98:', &
    & info,itry,uvec(1),uvec(2),xnse(i_nn),xnse(i_pp),xnse(i_he4),xnse(i_ni56),fvec(1),fvec(2)

    ! Save the CP98 result
    itry0 = itry
    info0 = info
    uvec0(:) = uvec(:)

    IF ( iscrn > 0 ) THEN
      ! Now use XNet's screening to match network
      use_CP98 = .false.
      IF ( info > 0 ) THEN
        itry = 0
        ! Iterate a few times to make sure screening and composition are consistent
        DO i = 1, 3
          CALL nse_screen
          CALL nse_nr(uvec,check,info)
          IF ( itsout >=2 ) WRITE(stdout,'(a5,2i2,8es23.15)') 'XNET:', &
          & info,itry,uvec(1),uvec(2),xnse(i_nn),xnse(i_pp),xnse(i_he4),xnse(i_ni56),fvec(1),fvec(2)
        END DO
      END IF

      ! Try different guesses if both screening approaches fail
      IF ( info <= 0 ) THEN
        DO itry = 1, 6
          unse(:) = 0.0d0
          xnse(:) = 0.0d0
          ynse(:) = 0.0d0
          hnse(:) = 0.0d0
          CALL nse_guess(uvec)
          CALL nse_screen
          CALL nse_nr(uvec,check,info)
          IF ( info > 0 ) EXIT
        END DO
        IF ( itsout >=2 ) WRITE(stdout,'(a5,2i2,8es23.15)') 'XNET:', &
        & info,itry,uvec(1),uvec(2),xnse(i_nn),xnse(i_pp),xnse(i_he4),xnse(i_ni56),fvec(1),fvec(2)
      END IF
    END IF

    IF ( info <= 0 ) THEN
      IF ( itsout >= 2 ) THEN
        IF ( info == 0 ) THEN
          WRITE(stdout,'(a)') 'Exceeded max NR iterations'
        ELSE IF ( info == -1 ) THEN
          WRITE(stdout,'(a)') 'Exceeded max LS iterations'
        ELSE IF ( info == -2 ) THEN
          WRITE(stdout,'(a)') 'Failed to find sufficient decrease in LS'
        ELSE IF ( info == -3 ) THEN
          WRITE(stdout,'(a)') 'Could not converge to solution'
        ELSE IF ( info == -4 ) THEN
          WRITE(stdout,'(a)') 'Non-negative slope in LS'
        ELSE IF ( info == -6 ) THEN
          WRITE(stdout,'(a)') 'Minimizer is not a root'
        ELSE IF ( info == -7 ) THEN
          WRITE(stdout,'(a)') 'Singular matrix'
        END IF
      END IF
      IF ( info0 > 0 ) THEN
        IF ( itsout >= 2 ) WRITE(stdout,'(a)') 'XNet screening failed; using CP98 screening result'
        itry = 7+itry0
        iconverge(itry) = iconverge(itry) + 1
      ELSE
        WRITE(stdout,'(a,i2,a)') 'Check convergence of root finder (',info,')'
        STOP 'ERROR NR Failed'
      END IF
    ELSE
      iconverge(itry) = iconverge(itry) + 1
    END IF
    
    ! Evaluate NSE composition with converged solution
    IF ( info > 0 ) THEN
      use_CP98 = .false.
      CALL nse_eval(uvec)
      knrtot(3) = knrtot(3) + 1
    ELSE IF ( info0 > 0 ) THEN
      ! If XNet screening failed to converge, use the CP98 result
      use_CP98 = .true.
      uvec(:) = uvec0(:)
      CALL nse_screen_CP98
      CALL nse_eval(uvec)
      knrtot(3) = knrtot(3) + 1
    END IF

    stop_timer = xnet_wtime()
    timer_nsesolv = timer_nsesolv + stop_timer

    IF ( itsout >= 0 ) THEN
      Write(stdout,'(a)')             'NSE solved'
      IF ( itsout >= 1 ) THEN
        Write(stdout,'(a,2es23.15)')  'roots                       =',uvec(1),uvec(2)
        Write(stdout,'(a,2es23.15)')  'residuals                   =',fvec(1),fvec(2)
      END IF
    END IF

    RETURN
  END SUBROUTINE nse_solve

  SUBROUTINE nse_guess(uvec)

    ! Output variables
    REAL(8), INTENT(OUT) :: uvec(2)

    ! Local variables
    REAL(8) :: xmax, xmin
    REAL(8) :: bkt, c1, c2, det, lhs(2,2), rhs(2)
    REAL(8) :: zatst, za_ye(2), za_min, za_max
    INTEGER :: i, ii, iye(2), ibe
    INTEGER :: iye_za, iye_za_l, iye_za_u, iza_min, iza_max
    INTEGER :: iye_hvy, iye_hvy_l, iye_hvy_u, ihvy_min, ihvy_max
    INTEGER :: iye_fe, iye_fe_l, iye_fe_u, ife_min, ife_max

    xmax = 1.0d0
    xmin = 0.0d0

    ! Find maximally bound nuclei
    ibe = MAXLOC( be(:)/aa(:), 1 )

    ! Find min/max Z/A in network
    iza_min = 3 + MINLOC( zz(3:ny)/aa(3:ny), 1 ) - 1
    iza_max = 3 + MAXLOC( zz(3:ny)/aa(3:ny), 1 ) - 1
    za_min = zz(iza_min)/aa(iza_min)
    za_max = zz(iza_max)/aa(iza_max)

    ! Find A>1 nuclei with Z/A closest to Ye from above and below 
    iye_za_u = ny + 1 - MINLOC( ABS(zz(ny:1:-1)/aa(ny:1:-1) - yense), 1, &
    &                           MASK = zz(ny:1:-1)/aa(ny:1:-1) >= yense .and. &
    &                                  aa(ny:1:-1)>1.0d0 )
    IF ( iye_za_u > ny ) THEN
      iye_za_u = iza_max
    END IF
    iye_za_l = MINLOC( ABS(zz(:)/aa(:) - yense), 1, &
    &                  MASK = zz(:)/aa(:) <= yense .and. &
    &                         aa(:)>1.0d0 )
    IF ( iye_za_l < 1 ) THEN
      iye_za_l = iza_min
    END IF

    zatst = zz(iye_za_l)/aa(iye_za_l)
    IF ( ABS(zatst-yense) < ABS(zz(iye_za_u)/aa(iye_za_u)-yense) ) THEN
      iye_za = iye_za_l
    ELSE
      iye_za = iye_za_u
    END IF

    ! Find Z>=20 nuclei with Z/A closest to Ye from above and below 
    ihvy_min = MINLOC( zz(:)/aa(:), 1, zz(:)>=20.0d0 )
    ihvy_max = ny + 1 - MAXLOC( zz(ny:1:-1)/aa(ny:1:-1), 1, zz(ny:1:-1)>=20.0d0 )

    iye_hvy_u = ny + 1 - MINLOC( ABS(zz(ny:1:-1)/aa(ny:1:-1) - yense), 1, &
    &                           MASK = zz(ny:1:-1)/aa(ny:1:-1) >= yense .and. &
    &                                  zz(ny:1:-1)>=20.0d0 )
    IF ( iye_hvy_u > ny ) THEN
      iye_hvy_u = ihvy_max
    END IF
    iye_hvy_l = MINLOC( ABS(zz(:)/aa(:) - yense), 1, &
    &                  MASK = zz(:)/aa(:) <= yense .and. &
    &                         zz(:)>=20.0d0 )
    IF ( iye_hvy_l < 1 ) THEN
      iye_hvy_l = ihvy_min
    END IF

    zatst = zz(iye_hvy_l)/aa(iye_hvy_l)
    IF ( ABS(zatst-yense) < ABS(zz(iye_hvy_u)/aa(iye_hvy_u)-yense) ) THEN
      iye_hvy = iye_hvy_l
    ELSE
      iye_hvy = iye_hvy_u
    END IF

    ! Find Fe/Ni nuclei with Z/A closest to Ye from above and below 
    ife_min = MINLOC( zz(:)/aa(:), 1, (zz(:)==26.0d0 .or. zz(:)==28.0d0) .and. &
    &                                 (aa(:)==52.0d0 .or. aa(:)==54.0d0 .or. &
    &                                  aa(:)==56.0d0 .or. aa(:)==58.0d0) )
    ife_max = ny + 1 - MAXLOC( zz(ny:1:-1)/aa(ny:1:-1), 1, (zz(ny:1:-1)==26.0d0 .or. zz(ny:1:-1)==28.0d0) .and. &
    &                                                      (aa(ny:1:-1)==52.0d0 .or. aa(ny:1:-1)==54.0d0 .or. &
    &                                                       aa(ny:1:-1)==56.0d0 .or. aa(ny:1:-1)==58.0d0) )

    iye_fe_u = ny + 1 - MINLOC( ABS(zz(ny:1:-1)/aa(ny:1:-1) - yense), 1, &
    &                           MASK = zz(ny:1:-1)/aa(ny:1:-1) >= yense .and. &
    &                                  ( zz(ny:1:-1)==26.0d0 .or. zz(ny:1:-1)==28.0d0 ) .and. &
    &                                  ( aa(ny:1:-1)==52.0d0 .or. aa(ny:1:-1)==54.0d0 .or. &
    &                                    aa(ny:1:-1)==56.0d0 .or. aa(ny:1:-1)==58.0d0 ) )
    IF ( iye_fe_u > ny ) THEN
      iye_fe_u = ife_max
    END IF
    iye_fe_l = MINLOC( ABS(zz(:)/aa(:) - yense), 1, &
    &                  MASK = zz(:)/aa(:) <= yense .and. &
    &                         ( zz(:)==26.0d0 .or. zz(:)==28.0d0 ) .and. &
    &                         ( aa(:)==52.0d0 .or. aa(:)==54.0d0 .or. &
    &                           aa(:)==56.0d0 .or. aa(:)==58.0d0 ) )
    IF ( iye_fe_l < 1 ) THEN
      iye_fe_l = ihvy_min
    END IF

    zatst = zz(iye_fe_l)/aa(iye_fe_l)
    IF ( ABS(zatst-yense) < ABS(zz(iye_fe_u)/aa(iye_fe_u)-yense) ) THEN
      iye_fe = iye_fe_l
    ELSE
      iye_fe = iye_fe_u
    END IF

    IF ( itry <= 1 .and. ( (rhonse <= 1.0d5 .and. t9nse >= 7.0d0) .or. &
    &                      (rhonse <= 1.0d6 .and. t9nse >= 8.0d0) .or. &
    &                      (rhonse <= 1.0d7 .and. t9nse >= 9.0d0) .or. &
    &                      (rhonse <= 1.0d8 .and. t9nse >= 11.0d0) .or. &
    &                      (rhonse <= 1.0d9 .and. t9nse >= 14.0d0) ) ) THEN
      iye(1) = i_nn
      iye(2) = i_pp
    ELSE IF ( yense < za_min ) THEN
      IF ( itry <= 1 ) THEN
        iye(1) = iye_hvy_l
        iye(2) = i_nn
      ELSE IF ( itry == 2 ) THEN
        iye(1) = i_ni56
        iye(2) = i_nn
      ELSE IF ( itry == 3 ) THEN
        iye(:) = i_ni56
      ELSE IF ( itry == 4 ) THEN
        iye(:) = iye_hvy_l
      ELSE IF ( itry == 5 ) THEN
        iye(1) = ibe
        iye(2) = i_nn
      ELSE IF ( itry == 6 ) THEN
        iye(:) = ibe
      ELSE
        iye(1) = i_nn
        iye(2) = i_pp
      END IF
    ELSE IF ( yense >= za_min .and. yense < 0.499d0 ) THEN
      IF ( itry <= 1 ) THEN
        iye(1) = iye_fe_l
        iye(2) = iye_fe_u
      ELSE IF ( itry == 2 ) THEN
        iye(1) = i_ni56
        iye(2) = iye_fe_l
      ELSE IF ( itry == 3 ) THEN
        iye(1) = i_he4
        iye(2) = iye_fe_l
      ELSE IF ( itry == 4 ) THEN
        iye(:) = i_ni56
      ELSE IF ( itry == 5 ) THEN
        iye(1) = iye_hvy_l
        iye(2) = iye_hvy_u
      ELSE IF ( itry == 6 ) THEN
        iye(:) = ibe
      ELSE
        iye(1) = i_nn
        iye(2) = i_pp
      END IF
    ELSE IF ( yense >= 0.501d0 ) THEN
      IF ( itry <= 1 ) THEN
        iye(1) = i_ni56
        iye(2) = i_pp
      ELSE IF ( itry == 2 ) THEN
        iye(:) = i_ni56
      ELSE IF ( itry == 3 ) THEN
        iye(1) = i_he4
        iye(2) = iye_za_u
      ELSE IF ( itry == 4 ) THEN
        iye(1) = i_ni56
        iye(2) = iye_za_u
      ELSE IF ( itry == 5 ) THEN
        iye(:) = ibe
      ELSE IF ( itry == 6 ) THEN
        iye(1) = ibe
        iye(2) = i_pp
      ELSE
        iye(1) = i_nn
        iye(2) = i_pp
      END IF
    ELSE
      IF ( itry <= 1 ) THEN
        iye(:) = i_ni56
      ELSE IF ( itry == 2 ) THEN
        iye(:) = iye_fe_l
      ELSE IF ( itry == 3 ) THEN
        iye(:) = ibe
      ELSE IF ( itry == 4 ) THEN
        iye(:) = iye_hvy_u
      ELSE IF ( itry == 5 ) THEN
        iye(:) = i_si28
      ELSE IF ( itry == 6 ) THEN
        iye(:) = i_he4
      ELSE
        iye(1) = i_nn
        iye(2) = i_pp
      END IF
    END IF
    za_ye(:) = zz(iye(:))/aa(iye(:))

    IF ( itry /= 0 ) THEN
      IF ( iye(1) == iye(2) ) THEN
        xnse(iye(1)) = xmax
      ELSE IF ( za_ye(1) == za_ye(2) ) THEN
        xnse(iye(1)) = xmax
        xnse(iye(2)) = xmin
      ELSE
        xnse(iye(1)) = MAX( MIN( ( yense - za_ye(2) ) / ( za_ye(1) - za_ye(2) ), xmax ), xmin )
        xnse(iye(2)) = MAX( MIN( 1.0d0 - xnse(iye(1)), xmax ), xmin )
      END IF
      ynse(iye) = xnse(iye) / (mm(iye)*avn)
    ELSE
      iye(:) = MINLOC( (zz(:)*xnse(:)/aa(:) - yense)**2 + (xnse(:)-1.0d0)**2, 1 )
!     iye(1) = i_nn
!     iye(2) = i_pp
    END IF
    
    ! Update screening corrections
    CALL nse_screen

    bkt = t9nse*bok*epmev
    c1 = bkt/(2.0*pi*hbar*hbar*epmev*epmev)
    IF ( iye(1) == iye(2) .or. ANY( xnse(iye(:)) <= 0.0d0 ) ) THEN
      ! Assume un = up
      ii = MAXLOC( xnse(iye(:)), 1 )
      i = iye(ii)
      c2 = mm52(i) * angm(i)*gg(i) * c1 * SQRT(c1) / rhonse
      uvec(1) = (bkt*(LOG(xnse(i)/c2) - hnse(intz(i))) - be(i)*epmev)/aa(i)
      uvec(2) = uvec(1)
    ELSE 
      ! Solve for un and up
      DO ii = 1, 2
        i = iye(ii)
        c2 = mm52(i) * angm(i)*gg(i) * c1 * SQRT(c1) / rhonse
        lhs(ii,1) = nn(i)
        lhs(ii,2) = zz(i)
        rhs(ii) = bkt*(LOG(xnse(i)/c2) - hnse(intz(i))) - be(i)*epmev
      END DO
      det = lhs(1,1)*lhs(2,2)-lhs(1,2)*lhs(2,1)
      uvec(1) = rhs(1)*lhs(2,2) - rhs(2)*lhs(1,2)
      uvec(2) = rhs(2)*lhs(1,1) - rhs(1)*lhs(2,1)
      uvec(:) = uvec(:) / det
    END IF
    IF ( itsout >= 2 ) THEN
      WRITE(stdout,'(a,2es23.15)') 'Initial guess               =', uvec(1), uvec(2)
      WRITE(stdout,'(a5,f11.7,es23.15)') nname(iye(1)), za_ye(1), xnse(iye(1))
      WRITE(stdout,'(a5,f11.7,es23.15)') nname(iye(2)), za_ye(2), xnse(iye(2))
    END IF
    IF ( itsout >= 2 ) THEN
      WRITE(stdout,'(a,2a5)') 'iye_za_l,  iye_za_u  =', nname(iye_za_l),  nname(iye_za_u)
      WRITE(stdout,'(a,2a5)') 'iye_hvy_l, iye_hvy_u =', nname(iye_hvy_l), nname(iye_hvy_u)
      WRITE(stdout,'(a,2a5)') 'iye_fe_l,  iye_fe_u  =', nname(iye_fe_l),  nname(iye_fe_u)
    END IF

    RETURN
  END SUBROUTINE nse_guess

  SUBROUTINE nse_composition(uvec) 

    ! Input variables
    REAL(8), INTENT(IN) :: uvec(2)         ! chemical potentials for n & p

    ! Local variables
    REAL(8) :: c1,temp(ny)
    REAL(8) :: bkt,bktinv,rhoinv
    INTEGER :: ii,jj

    start_timer = xnet_wtime()
    timer_nseeval = timer_nseeval - start_timer

    ! Useful scalars
    bkt = t9nse*bok*epmev
    bktinv = 1.0/bkt
    rhoinv = 1.0/rhonse
    c1 = bkt/(2.0*pi*hbar*hbar*epmev*epmev)
    c1 = c1*SQRT(c1)

    ! Evaluate chemical potentials and mass fractions
    unse(:) = nn(:)*uvec(1) + zz(:)*uvec(2)
    xnse(:) = c1 * rhoinv * angm(1:ny)*gg(1:ny) * mm52(:)
    temp(:) = (unse(:)+be(:)*epmev)*bktinv + hnse(intz(:))
    IF ( itsout >= 5 ) THEN
      ii = MAXLOC( temp, 1 )
      jj = MINLOC( temp, 1 )
      WRITE(stdout,'(a,a5,4es13.5)') 'max exponent: ',nname(ii), temp(ii), unse(ii)*bktinv, be(ii)*epmev*bktinv, hnse(intz(ii))
      WRITE(stdout,'(a,a5,4es13.5)') 'min exponent: ',nname(jj), temp(jj), unse(jj)*bktinv, be(jj)*epmev*bktinv, hnse(intz(jj))
    END IF
    xnse(:) = xnse(:) * safe_exp( temp(:) )
    ynse(:) = xnse(:) / (mm(:)*avn)

    stop_timer = xnet_wtime()
    timer_nseeval = timer_nseeval + stop_timer

    RETURN
  END SUBROUTINE nse_composition

  SUBROUTINE nse_eval(uvec) 

    ! Input variables
    REAL(8), INTENT(IN) :: uvec(2)         ! chemical potentials for n & p

    ! Local variables
    REAL(8) :: temp(ny)

    ! Update the NSE composition variables, unse, xnse, ynse
    CALL nse_composition(uvec)

    start_timer = xnet_wtime()
    timer_nseeval = timer_nseeval - start_timer

    ! Mass and charge conservation
    fvec(1) = SUM(xnse(:)) - 1.0d0

    temp(:) = (zz(:)-aa(:)*yense)/(avn*mm(:))
    fvec(2) = SUM( temp(:)*xnse(:) )
!   temp(:) = zz(:)/(avn*mm(:))
!   fvec(2) = SUM( temp(:)*xnse(:) ) - yense

    knr(3) = knr(3) + 1

    stop_timer = xnet_wtime()
    timer_nseeval = timer_nseeval + stop_timer

    RETURN
  END SUBROUTINE nse_eval

  SUBROUTINE nse_jac(uvec,feval)

    ! Input variables
    REAL(8), INTENT(IN) :: uvec(2)
    LOGICAL, INTENT(IN) :: feval

    ! Local variables
    REAL(8) :: dxdn(ny),dxdp(ny)
    REAL(8) :: bkt,bktinv,temp(ny)

    ! Useful scalars
    bkt = t9nse*bok*epmev
    bktinv = 1.0d0/bkt

    IF ( feval ) THEN
      CALL nse_eval(uvec)
    ELSE
      CALL nse_composition(uvec)
    END IF

    temp(:) = (zz(:)-aa(:)*yense)/(avn*mm(:))
!   temp(:) = zz(:)/(avn*mm(:))

    ! Jacobian
    dxdn(:)   = bktinv * nn(:) * xnse(:)
    fjac(1,1) = SUM( dxdn(:) )
    fjac(2,1) = SUM( temp(:)*dxdn(:) )

    dxdp(:)   = bktinv * zz(:) * xnse(:)
    fjac(1,2) = SUM( dxdp(:) )
    fjac(2,2) = SUM( temp(:)*dxdp(:) )

    RETURN
  END SUBROUTINE nse_jac

  SUBROUTINE nse_nr(uvec,check,info)

    ! Input/Output variables
    REAL(8), INTENT(INOUT) :: uvec(2)

    ! Output variables
    LOGICAL, INTENT(OUT) :: check
    INTEGER, INTENT(OUT) :: info

    ! Local variables
    REAL(8) :: uvec0(2), gvec(2), pvec(2)
    REAL(8) :: det, jcond, fjaci(2,2)
    REAL(8) :: f, f0, maxstep, unorm
    REAL(8) :: testf, testu, testg

    INTEGER :: nrit, ipiv(2)
    INTEGER :: neval
    INTEGER :: i, n, info2

    ! Initialize return values
    info = 0
    info2 = -1
    check = .false.
    knr(:) = 0

    n = SIZE(uvec)
    testu = 0.0d0
    testg = 0.0d0

    ! Check for initial guess being correct, but with stricter tolerance
    f = nse_func(uvec)
    knrtot(:) = knrtot(:) + knr(:)
    testf = MAXVAL( ABS(scalefvec(:)*fvec(:)), 1 )
    IF ( testf < 0.01d0*tolf ) THEN
      check = .false.
      info = 1
      RETURN
    END IF

    start_timer = xnet_wtime()
    timer_nsenrap = timer_nsenrap - start_timer

    ! Calculate maximum relative step length
    unorm = l2norm( scaleu(:)*uvec(:) )
    maxstep = maxstep0 * MAX( unorm, l2norm( scaleu(:) ) )

    IF ( itsout >= 2 ) WRITE(stdout,'(3i3,7es23.15)') 0, knr(2:3), uvec(1), uvec(2), fvec(1), fvec(2), f, testf, testu
    IF ( itsout >= 3 ) WRITE(stdout,'(9x,2es23.15)') xnse(1), xnse(2)

    DO nrit = 1, nritmax

      ! Reset per iteration counters
      knr(:) = 0

      ! Build jacobian
!     CALL nse_jac(uvec,.true.)

      ! Get estimate of condition number
      det = fjac(1,1)*fjac(2,2) - fjac(1,2)*fjac(2,1)
      IF ( det > 2.0d0*EPSILON(det) ) THEN
        fjaci(1,1) = fjac(2,2)
        fjaci(2,1) = -fjac(2,1)
        fjaci(1,2) = -fjac(1,2)
        fjaci(2,2) = fjac(1,1)
        fjaci(:,:) = fjaci(:,:) / det
        jcond = MAXVAL( SUM( ABS(fjac(:,:)), 1 ) ) * MAXVAL( SUM( ABS(fjaci(:,:)), 1 ) )
      ELSE
        jcond = 0.0d0
      END IF

      IF ( det <= 2.0d0*EPSILON(det) .or. jcond*tolmin > 1.0d0 ) THEN
        IF ( itsout >=1 ) WRITE(stdout,'(a,2es23.15)') 'Singular Matrix: det,jcond = ',det,jcond
        info = -7
        EXIT
      END IF

      ! Calculate search direction, pvec
      pvec(1) = -fvec(1)*fjaci(1,1) - fvec(2)*fjaci(1,2)
      pvec(2) = -fvec(1)*fjaci(2,1) - fvec(2)*fjaci(2,2)

      ! Calculate gradient of fvec, gvec
      CALL dgemv('T',n,n,1.0d0,fjac,n,scalefvec(:)*scalefvec(:)*fvec(:),1,0.0d0,gvec,1)

!     pvec(:) = -fvec(:)
!     CALL dgesv(n,1,fjac,n,ipiv,pvec,n,info)

      ! Save old values
      uvec0(:) = uvec(:)
      f0 = f

      ! Backtracking line search to get step length and update uvec, f, fvec
      CALL nse_lnsrch(uvec0,f0,uvec,f,gvec,pvec,maxstep,check,info2)

      knrtot(:) = knrtot(:) + knr(:)

      ! test for convergence
      IF ( lnorm == 0 ) THEN
        testf = MAXVAL( ABS(scalefvec(:)*fvec(:)), 1 )
        testu = MAXVAL( ABS(uvec(:)-uvec0(:)) / MERGE( ABS(uvec(:)), typu(:), ABS(uvec(:))>typu(:) ), 1 )
        testg = MAXVAL( ABS(gvec(:)) * MERGE( ABS(uvec(:)), typu(:), ABS(uvec(:))>typu(:) ) / MAX(f, 0.5d0*n), 1 )
      ELSE IF ( lnorm == 1 ) THEN
        testf = SUM( ABS(scalefvec(:)*fvec(:)) )
        testu = SUM( ABS(uvec(:)-uvec0(:)) / MERGE( ABS(uvec(:)), typu(:), ABS(uvec(:))>typu(:) ) )
        testg = SUM( ABS(gvec(:)) * MERGE( ABS(uvec(:)), typu(:), ABS(uvec(:))>typu(:) ) / MAX(f, 0.5d0*n) )
      ELSE IF ( lnorm == 2 ) THEN
        testf = l2norm( ABS(scalefvec(:)*fvec(:)) )
        testu = l2norm( ABS(uvec(:)-uvec0(:)) / MERGE( ABS(uvec(:)), typu(:), ABS(uvec(:))>typu(:) ) )
        testg = l2norm( ABS(gvec(:)) * MERGE( ABS(uvec(:)), typu(:), ABS(uvec(:))>typu(:) ) / MAX(f, 0.5d0*n) )
      END IF

      IF ( itsout >= 2 ) WRITE(stdout,'(3i3,8es23.15)') nrit, knr(2:3), uvec(1), uvec(2), fvec(1), fvec(2), f, testf, testu, testg
      IF ( itsout >= 3 ) WRITE(stdout,'(7x,i2,6es23.15)') info2, xnse(1), xnse(2), pvec(1), pvec(2), gvec(1), gvec(2)

      IF ( info2 < 0 ) THEN
        info = info2
        EXIT
      ELSE IF ( testf < tolf ) THEN
        info = 1
        EXIT
      ELSE IF ( check ) THEN
        IF ( testg < tolmin ) THEN
          info = -6
          EXIT
        ELSE
          info = -3
          EXIT
        END IF
      ELSE IF ( testu < tolu ) THEN
        info = 2
        EXIT
      ELSE
        info = 0
      END IF

    END DO

    IF ( info <= 0 ) THEN
      check = .true.
    END IF
    knr(1) = nrit
    knrtot(1) = knrtot(1) + MIN( nrit, nritmax )

    stop_timer = xnet_wtime()
    timer_nsenrap = timer_nsenrap + stop_timer

    RETURN
  END SUBROUTINE nse_nr

  SUBROUTINE nse_lnsrch(uvec0,f0,uvec,f,gvec,pvec,maxstep,check,info)

    ! Input variables
    REAL(8), INTENT(IN) :: uvec0(2),f0,maxstep

    ! Input/Output variables
    REAL(8), INTENT(INOUT) :: pvec(2), gvec(2)

    ! Output variables
    REAL(8), INTENT(OUT) :: uvec(2),f
    LOGICAL, INTENT(OUT) :: check
    INTEGER, INTENT(OUT) :: info

    ! Local variables
    INTEGER :: i, lsit, ipiv(4), ierr, j

    REAL(8) :: gveck(2)
    REAL(8) :: betak(0:lsitmax+1), fk(0:lsitmax+1), slopek(0:lsitmax+1)
    REAL(8) :: newtlen, slope, plength, beta_min
    REAL(8) :: rbsq, rb0sq, rbeta, rbeta0, rdbeta
    REAL(8) :: rhs(4), c1, c2, c3, disc, vmat(4,4)

    check = .false.
    info = -1

    ! Scale step length if necessary
    newtlen = l2norm( scaleu(:)*pvec(:) )
    IF ( newtlen > maxstep ) THEN
      pvec(:) = pvec(:) * maxstep/newtlen
    END IF

    ! Calculate "slope" of step direction
    slope = dot_prod( gvec(:), pvec(:) )
!   slope = -dot_prod( fvec(:), fvec(:) )
    IF ( slope >= 0.0d0 ) THEN
      info = -4
      RETURN
    END IF

    ! Determine minimum allowable step length
    plength = MAXVAL( ABS(pvec(:)) / MERGE( ABS(uvec0(:)), typu(:), ABS(uvec0(:))>typu(:) ), 1 )
    beta_min = tolbeta/plength

    ! Start with full step length
    fk(:) = f0
    betak(0) = 0.0d0
    betak(1) = 1.0d0
    slopek(:) = slope
    IF ( itsout >= 4 ) WRITE(stdout,'(3x,i3,3x,4es23.15)') 0, beta_min, f0, f0+alpha*betak(1)*slope, slope

    DO lsit = 1, lsitmax

      ! Update uvec, f, fvec
      uvec(:) = uvec0(:) + betak(lsit)*pvec(:)
      fk(lsit) = nse_func(uvec)

      IF ( itsout >= 4 ) WRITE(stdout,'(3x,i3,3x,5es23.15)') &
      &         lsit, betak(lsit), fk(lsit), fk(lsit-1) + alpha*betak(lsit)*slopek(lsit-1), slopek(lsit-1), &
      &         dot_prod( (gveck(:)-gvec(:)), pvec )

      ! Test if sufficient decrease condition has been met
      IF ( fk(lsit) <= f0 + alpha*betak(lsit)*slope ) THEN
        check = .false.
        info = 0
        f = fk(lsit)
        EXIT
      ELSE IF ( betak(lsit) < beta_min ) THEN
        check = .true.
        info = 3
        uvec(:) = uvec0(:)
        f = f0
        EXIT
      ELSE IF ( lsit > 200 ) THEN
        CALL dgemv('T',2,2,1.0d0,fjac,2,scalefvec(:)*scalefvec(:)*fvec(:),1,0.0d0,gveck,1)
        slopek(lsit) = dot_prod( gveck(:), pvec(:) )
        rhs(:) = fk(lsit-3:lsit)
        vmat(:,1) = betak(lsit-3:lsit)**3
        vmat(:,2) = betak(lsit-3:lsit)**2
        vmat(:,3) = betak(lsit-3:lsit)
        vmat(:,4) = 1.0d0
      ELSE IF ( lsit > 100 ) THEN
        CALL dgemv('T',2,2,1.0d0,fjac,2,scalefvec(:)*scalefvec(:)*fvec(:),1,0.0d0,gveck,1)
        slopek(lsit) = dot_prod( gveck(:), pvec(:) )
        rhs(1:2) = fk(lsit-1:lsit)
        rhs(3:4) = slopek(lsit-1:lsit)
        vmat(1:2,1) = betak(lsit-1:lsit)**3
        vmat(3:4,1) = 3.0d0*betak(lsit-1:lsit)**2
        vmat(1:2,2) = betak(lsit-1:lsit)**2
        vmat(3:4,2) = 2.0d0*betak(lsit-1:lsit)
        vmat(1:2,3) = betak(lsit-1:lsit)
        vmat(3:4,3) = 1.0d0
        vmat(1:2,4) = 1.0d0
        vmat(3:4,4) = 0.0d0
      ELSE IF ( lsit >= 1 ) THEN
!       ! Choose next step length that minimizes the quadratic interpolating polynomial:
!       !     P = c1*beta(k+1)^2 + c2*beta(k+1) + c3
!       c2 = (fk(lsit) - f0 - slope*betak(lsit))*rbsq
!       c3 = slope
!       betak(lsit+1) = -0.5d0 * c2 / c1

        ! We know fk(0), slope(0), fk(k), and slope(k), so...
        ! Find betak(k+1) that minimizes the cubic interpolating polynomial: 
        !     P = c1*beta(k+1)^3 + c2*beta(k+1)^2 + c3*beta(k+1) + c4
        !     Get the coefficients from cn = V^-1 * y_n, for Vandermonde matrix V
        CALL dgemv('T',2,2,1.0d0,fjac,2,scalefvec(:)*scalefvec(:)*fvec(:),1,0.0d0,gveck,1)
        slopek(lsit) = dot_prod( gveck(:), pvec(:) )
        rhs(1:2) = fk(lsit-1:lsit)
        rhs(3:4) = slopek(lsit-1:lsit)
        vmat(1:2,1) = betak(lsit-1:lsit)**3
        vmat(3:4,1) = 3.0d0*betak(lsit-1:lsit)**2
        vmat(1:2,2) = betak(lsit-1:lsit)**2
        vmat(3:4,2) = 2.0d0*betak(lsit-1:lsit)
        vmat(1:2,3) = betak(lsit-1:lsit)
        vmat(3:4,3) = 1.0d0
        vmat(1:2,4) = 1.0d0
        vmat(3:4,4) = 0.0d0
      ELSE
        CALL dgemv('T',2,2,1.0d0,fjac,2,scalefvec(:)*scalefvec(:)*fvec(:),1,0.0d0,gveck,1)
        slopek(lsit) = dot_prod( gveck(:), pvec(:) )
        ! We know fk(0), slope(0), fk(k-1), and fk(k), so...
        rhs(1:3) = fk(lsit-2:lsit)
        rhs(4) = slopek(lsit-2)
        vmat(1:3,1) = betak(lsit-2:lsit)**3
        vmat(4,1) = 3.0d0*betak(lsit-2)**2
        vmat(1:3,2) = betak(lsit-2:lsit)**2
        vmat(4,2) = 2.0d0*betak(lsit-2)
        vmat(1:3,3) = betak(lsit-2:lsit)
        vmat(4,3) = 1.0d0
        vmat(1:3,4) = 1.0d0
        vmat(4,4) = 0.0d0
      END IF

      CALL dgesv(4,1,vmat,4,ipiv,rhs,4,ierr)
      c1 = rhs(1)
      c2 = rhs(2)
      c3 = rhs(3)
      IF ( c1 == 0.0d0 ) THEN
        ! Minimum of quadratic interpolating polynomial
        betak(lsit+1) = -0.5d0*c3/c2
      ELSE
        ! Choose the root which corresponds to the local minimum
        disc = c2*c2 - 3.0d0*c1*c3
        IF ( disc < 0.0d0 ) THEN
          ! No roots exist, so set next step length to its upper-bound
          betak(lsit+1) = gamma2*betak(lsit)
        ELSE IF ( c2 < 0.0d0 ) THEN
          betak(lsit+1) = (-c2 + SQRT(disc)) / (3.0d0*c1)
        ELSE
          betak(lsit+1) = -c3 / (c2 + SQRT(disc))
        END IF
      END IF

      ! Limit change in step length
      betak(lsit+1) = MIN( betak(lsit+1), gamma2*betak(lsit) )
      betak(lsit+1) = MAX( betak(lsit+1), gamma1*betak(lsit) )
    END DO

    IF ( lsit > lsitmax ) THEN
      check = .true.
      f = fk(lsitmax)
    END IF
    knr(2) = lsit

    RETURN
  END SUBROUTINE nse_lnsrch

  REAL(8) FUNCTION nse_func(uvec)

    REAL(8), INTENT(IN) :: uvec(2)

    CALL nse_jac(uvec,.true.)

    nse_func = 0.5d0 * dot_prod( scalefvec(:)*fvec(:), scalefvec(:)*fvec(:) )

    RETURN
  END FUNCTION nse_func

  SUBROUTINE nse_pf

    ! Local variables
    REAL(8) :: rdt9
    INTEGER :: i, ii

    DO i = 1, 24
      IF ( t9i(i) >= t9nse ) EXIT
    END DO
    ii = i

    SELECT CASE (ii)
      CASE(2:24)
        rdt9 = (t9nse-t9i(ii-1)) / (t9i(ii)-t9i(ii-1))
        gg(1:ny) = EXP( rdt9*LOG(g(ii,1:ny)) + (1-rdt9)*LOG(g(ii-1,1:ny)) )
      CASE (1)
        gg(1:ny) = g(1,1:ny)
      CASE (25)
        gg(1:ny) = g(24,1:ny)
    END SELECT

    RETURN
  END SUBROUTINE nse_pf

  SUBROUTINE nse_screen

    ! Local variables
    INTEGER, PARAMETER :: iz1 = 1
    REAL(8), PARAMETER :: z1 = 1.0d0
    REAL(8) :: cds(5) = (/ -0.899172d0, 0.602249d0, -0.274823d0, -1.401915d0, 0.3230064d0 /)

    INTEGER :: j,iz2(zmax-1)
    REAL(8) :: z2(zmax-1), h0(zmax-1), hw(zmax-1), hi(zmax-1), hs(zmax-1), lambda12(zmax-1)
!   REAL(8) :: gamma12(zmax-1)
    REAL(8) :: fhs(0:zmax+1), fhi(0:zmax+1)
    REAL(8) :: ztot, ztilde, zinter, lambda0, gammae, dztildedt9
    REAL(8) :: gammaz, z

    start_timer = xnet_wtime()
    timer_nsescrn = timer_nsescrn - start_timer

    IF ( iscrn <= 0 ) THEN
      hnse = 0.0d0
    ELSE IF ( .not. use_CP98 ) THEN
      ! Call EOS to get plasma quantities
      CALL eos_interface(t9nse,rhonse,ynse,ztot,ztilde,zinter,lambda0,gammae,dztildedt9)

      !-----------------------------------------------------------------------------
      ! Calculate screening energies as a function of Z, for prescriptions that 
      ! follow this approach 
      !-----------------------------------------------------------------------------
      fhi(0) = 0.0d0
      fhs(0) = 0.0d0
      DO j = 1, zmax+1
        z = REAL( j, KIND(1.0d0) )
        fhi(j) = 0.38d0*zinter*lambda0**0.86d0 * z**1.86d0
        gammaz = z**fivethird*gammae                            
        fhs(j) = cds(1)*gammaz + cds(2)*gammaz**cds(5)/cds(5) + cds(3)*LOG(gammaz) + cds(4)
      END DO

      DO j = 1, zmax-1
        iz2(j) = j
      END DO
      z2(:) = REAL( iz2, KIND(1.0d0) )

      ! Weak and intermediate screening factors, Table 4 of Graboske et al.         
      lambda12 = z1*z2*ztilde*lambda0    
      hw = lambda12
!     hw = lambda12 * ( 1.0d0 + lambda12 * ( LOG(lambda12) + 0.8364d0 ) )

!     hi = 0.38d0*zinter*lambda0**0.86d0 * ( (z1+z2)**1.86d0 - z1**1.86d0 - z2**1.86d0 )
      hi = fhi(iz1+iz2) - fhi(iz1) - fhi(iz2)                                                     

      ! Strong screening from Dewitt & Slattery using linear mixing.
      hs = fhs(iz1) + fhs(iz2) - fhs(iz1+iz2)

      ! Strong Screening from Itoh, non-radial component
!     gamma12 = 2.0d0*z1*z2*gammae / ( z1**onethird + z2**onethird )
!     hs = 1.25d0*gamma12

      ! Select Screening factor
      WHERE ( iz2 == 0 )
        h0 = 0.0
      ELSE WHERE ( lambda12 < 0.1d0 )
        h0 = hw
      ELSE WHERE ( lambda12 < 2.0d0 )
        h0 = hi
      ELSE WHERE ( lambda12 > 5.0d0 )
        h0 = hs
      ELSE WHERE
        h0 = MERGE( hi, hs, hi<hs )
      END WHERE

      ! Add succeeding screening factors
      hnse(1) = 0.0d0
      DO j = 2, zmax
!       hnse(j) = hnse(j-1) + h0(j-1)
        hnse(j) = SUM( h0(1:(j-1)) )
      END DO

      IF ( idiag >= 4 ) THEN
        DO j = 2, zmax
          WRITE(lun_diag,'(a5,3i6,7es13.5)') 'HNSE',iz1,iz2(j-1),j,lambda12(j-1),h0(j-1),hw(j-1),hi(j-1),hs(j-1),hnse(j)
        END DO
      END IF

    END IF

    stop_timer = xnet_wtime()
    timer_nsescrn = timer_nsescrn + stop_timer

    RETURN
  END SUBROUTINE nse_screen

  SUBROUTINE nse_screen_CP98

    ! Local variables
    INTEGER, PARAMETER :: iz1 = 1
    REAL(8), PARAMETER :: z1 = 1.0d0
    REAL(8), PARAMETER :: a1 = -0.9052d0, a2 = 0.6322d0
    REAL(8), PARAMETER :: a2inv = 1.0d0 / a2
    REAL(8), PARAMETER :: a3 = -SQRT(3.0d0)*2.0d0 - a1*SQRT(a2inv)

    INTEGER :: j,iz2(zmax)
    REAL(8) :: z2(zmax),h0(zmax)
    REAL(8) :: fh0(zmax+1)
    REAL(8) :: ztot, ztilde, zinter, lambda0, gammae, dztildedt9
    REAL(8) :: gz, z, sqrtgz, ae

    start_timer = xnet_wtime()
    timer_nsescrn = timer_nsescrn - start_timer

    IF ( iscrn <= 0 ) THEN
      hnse = 0.0d0
    ELSE
      ae=(3./(4.*pi*avn*rhonse*yense))**onethird ! electron-sphere radius
      gammae=e2/(ae*t9nse*bok) ! electron Coulomb coupling parameter 

      !-----------------------------------------------------------------------------
      ! Calculate screening energies as a function of Z, for prescriptions that 
      ! follow this approach 
      !-----------------------------------------------------------------------------
      DO j = 1, zmax
        z = REAL( j, KIND(1.0d0) )
        gz = z**fivethird*gammae                            
        sqrtgz = SQRT(gz)
        fh0(j) = a1*( SQRT(gz*(a2+gz)) - a2*LOG(SQRT(gz*a2inv)+SQRT(1.0d0+gz*a2inv)) ) &
        &        + 2.0d0*a3*( sqrtgz - ATAN(sqrtgz) )
      END DO

      DO j = 1, zmax
        iz2(j) = j
      END DO
      z2(:) = REAL( iz2, KIND(1.0d0) )

!     h0(:) = z2(:)*fh0(1) - fh0(:)
      h0(:) = fh0(iz1) + fh0(iz2(:)) - fh0(iz1+iz2(:))

      hnse(1) = 0.0d0
      hnse(2:zmax) = z2(2:zmax)*fh0(1) - fh0(2:zmax)
      ! Add succeeding screening factors
!     hnse(1) = 0.0d0
!     DO j = 2, zmax
!       hnse(j) = SUM( h0(1:(j-1)) )
!     END DO

      IF ( idiag >= 4 ) THEN
        DO j = 2, zmax
          WRITE(lun_diag,'(a5,3i6,13x,es13.5,39x,es13.5)') 'HNSE',iz1,iz2(j-1),j,SUM(h0(1:(j-1))),hnse(j)
        END DO
      END IF

    END IF

    stop_timer = xnet_wtime()
    timer_nsescrn = timer_nsescrn + stop_timer

    RETURN
  END SUBROUTINE nse_screen_CP98

END MODULE nse_module