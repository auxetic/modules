module mo_fire
    use mo_syst
    use mo_config
    use mo_list
    use mo_network
    use mo_force
    implicit none

!-- fire var
    real(8), private, parameter :: fmax    = 1.d-12
    integer, private, parameter :: stepmax = 1e6
    integer, private :: step
!--
    real(8), private, parameter :: dt0   = 3.d-2
    real(8), private, parameter :: dtmax = 3.d-1
    real(8), private, parameter :: beta0 = 1.d-1
    real(8), private, parameter :: finc  = 1.1d0
    real(8), private, parameter :: fdec  = 0.5d0
    real(8), private, parameter :: fbeta = 0.99d0
    integer, private, parameter :: nmin  = 5
!--
    real(8), private :: dt, dt2, dt22, beta, power, fn, vn
    integer, private :: count


    type(tpcon)  :: confire, confire2
    type(tplist) :: nbfire

contains

    subroutine init_confire( tconfire, tcon, tnet )
        implicit none

        ! para list
        type(tpcon),     intent(inout)        :: tconfire
        type(tpcon),     intent(in)           :: tcon
        type(tpnetwork), intent(in), optional :: tnet

        ! copy con to confire
        tconfire = tcon

        ! if net don't exist, allocate list
        if ( .not. present( tnet ) ) then
            call init_list( nbfire, tconfire )
        end if

    end subroutine init_confire


    subroutine check_system_force( tcon, tnet )
        implicit none

        ! para list
        type(tpcon),     intent(inout)           :: tcon
        type(tpnetwork), intent(inout), optional :: tnet

        if ( .not. present( tnet ) ) then
            call make_list( nbfire, tcon )
            call calc_force( tcon, nbfire )
        else
           !call make_network( tnet, tcon )
            call calc_force_spring( tcon, tnet )
        end if

    end subroutine check_system_force

    ! main 1
    subroutine mini_fire_cv( tcon, tnet )
        implicit none

        ! para list
        type(tpcon),     intent(inout)        :: tcon
        type(tpnetwork), intent(in), optional :: tnet

        ! local
        logical :: nonnetwork_flag
        real(8) :: dt, beta, temp
        real(8) :: onembeta, betavndfn
        integer :: cumn

        ! network system or not
        nonnetwork_flag = .true.
        if ( present( tnet ) ) nonnetwork_flag = .false.

        associate(              &
            ra    => tcon.ra,   &
            va    => tcon.va,   &
            fa    => tcon.fa,   &
            natom => tcon.natom &
            )

            ! initial sets
            fa   = 0.d0 ; va  = 0.d0   ;
            dt   = dt0  ; beta = beta0 ;
            cumn = 0

            ! calc fortran before iteration
            if ( nonnetwork_flag ) then
                call make_list( nbfire, tcon )
                call calc_force( tcon, nbfire )
            else
                call calc_force_spring( tcon, tnet )
            end if

            ! main
            do step=1, stepmax

                if ( nonnetwork_flag ) then
                    if( check_list( nbfire, tcon ) ) then
                        call make_list( nbfire, tcon )
                    end if
                end if

                dt2  = 0.5d0 * dt
                dt22 = 0.5d0 * dt**2

                ! velocity verlet method / move 1
                ra = ra + va * dt + fa * dt22
                va = va + fa * dt2

                ! velocity verlet method / force
                if ( nonnetwork_flag ) then
                    call calc_force( tcon, nbfire )
                else
                    call calc_force_spring( tcon, tnet )
                end if

                ! velocity verlet method / move 2
                va = va + fa * dt2

                ! fire
                cumn  = cumn + 1
                power = sum( fa * va )

                fn = norm2( fa )
                vn = norm2( va )

                onembeta = 1.d0 - beta
                betavndfn = beta*vn/fn
                va = onembeta * va + betavndfn * fa
               !va = ( 1.d0 - beta ) * va + beta * (vn/fn) * fa

                if ( power > 0.d0 .and. cumn > nmin ) then
                    dt   = min( dt*finc, dtmax )
                    beta = beta * fbeta
                end if

                if ( power < 0.d0 ) then
                    cumn = 0
                    dt   = dt * fdec
                    beta = beta0
                    va   = 0.d0
                end if

                temp = maxval( abs(fa) )
                if ( temp < fmax ) exit

               !if( step == 1 ) write(*,*) 'step', '    dt         ', '              fmax   '

               !if( mod(step,1) == 0 ) then
               !    write(*,'(i6,2e25.15)') step, dt, tcon.stress
               !end if

            end do

        end associate

    end subroutine mini_fire_cv

    ! main 2
    subroutine mini_fire_cp( tcon, tnet, opboxp_set, opxyp_set, opxp_set, opyp_set, opstress_set )
        implicit none

        type(tpcon),     intent(inout)        :: tcon
        type(tpnetwork), intent(in), optional :: tnet
        real(8),         intent(in), optional :: opboxp_set    ! target press
        real(8),         intent(in), optional :: opxyp_set
        real(8),         intent(in), optional :: opxp_set, opyp_set
        real(8),         intent(in), optional :: opstress_set

        real(8) :: dstrain
        logical :: nonnetwork_flag
        real(8) :: boxp_set, xyp_set, xp_set, yp_set, stress_set
        logical :: cp_flag, boxcp_flag, xycp_flag, xcp_flag, ycp_flag
        logical :: cs_flag

        real(8) :: dt, beta, temp
        real(8) :: onembeta, betavndfn
        integer :: cumn

        ! 1. network or not
        nonnetwork_flag = .true.
        if ( present( tnet ) ) nonnetwork_flag = .false.

        ! 2.0 model of constant pressure
        cp_flag    = .false.
        boxcp_flag = .false.
        xycp_flag  = .false.    ! unused now
        xcp_flag   = .false.
        ycp_flag   = .false.
        ! 2.1 constant box press, change box xy length simutaneously
        if ( present( opboxp_set ) ) then
            if ( present(opxp_set) .or. present(opyp_set) .or. present(opxyp_set) ) then
                print*, "Error set in cp"
                stop
            end if
            boxp_set   = opboxp_set
            boxcp_flag = .true.
        end if
        ! 2.2 constant box press, change box xy length indenpently
        if ( present(opxyp_set) ) then
            if ( present(opxp_set) .or. present(opyp_set) ) then
                print*, "Error set in xyp_set"
                stop
            end if
            xp_set   = opxyp_set ; yp_set   = opxyp_set
            xcp_flag = .true.    ; ycp_flag = .true.
        else
            ! 2.3 constant press of x side of box
            if( present( opxp_set ) ) then
                xp_set = opxp_set
                xcp_flag = .true.
            end if
            ! 2.4 constant press of x side of box
            if ( present( opyp_set ) ) then
                yp_set = opyp_set
                ycp_flag = .true.
            end if
        end if

        ! cp or not
        if ( boxcp_flag .or. xycp_flag .or. xcp_flag .or. ycp_flag ) cp_flag = .true.

        ! 3. cs or not
        cs_flag = .false.
        if ( present( opstress_set ) ) then
            stress_set = opstress_set
            cs_flag    = .true.
        end if

        associate(                   &
            natom   => tcon.natom,   &
            ra      => tcon.ra,      &
            va      => tcon.va,      &
            fa      => tcon.fa,      &
            press   => tcon.press,   &
            pressx  => tcon.pressx,  &
            pressy  => tcon.pressy,  &
            strain  => tcon.strain,  &
            strainv => tcon.strainv, &
            strainf => tcon.strainf, &
            stress  => tcon.stress,  &
            la      => tcon.la,      &
            lainv   => tcon.lainv,   &
            lav     => tcon.lav,     &
            laf     => tcon.laf      &
            )

            ! initial sets
            fa      = 0.d0 ; va      = 0.d0
            lav     = 0.d0 ; laf     = 0.d0
            strainv = 0.d0 ; strainf = 0.d0
            dt      = dt0  ; beta    = beta0
            cumn = 0

            ! pre
            if ( nonnetwork_flag ) then
                call make_list( nbfire, tcon )
                call calc_force( tcon, nbfire )
            else
                call calc_force_spring( tcon, tnet )
            end if
            if ( boxcp_flag ) laf     = press       - boxp_set
            if ( xcp_flag   ) laf(1)  = pressx      - xp_set
            if ( ycp_flag   ) laf(2)  = pressy      - yp_set
            if ( cs_flag    ) strainf = stress_set - stress

            do step=1, stepmax

                if ( nonnetwork_flag .and. check_list( nbfire, tcon ) ) then
                    call make_list( nbfire, tcon )
                end if

                dt2  = dt    * 0.5d0
                dt22 = dt**2 * 0.5d0

                ! velocity verlet method / move 1 / config
                ra = ra + va * dt + fa * dt22
                va = va + fa * dt2

                ! velocity verlet method / move 1 / box
                if ( cp_flag ) then
                    la  = la  + lav * dt + laf * dt22
                    lav = lav + laf * dt2
                    !v affine deformation
                    if ( boxcp_flag ) ra      = ra      * ( lainv(1) * la(1) )
                    if ( xcp_flag   ) ra(1,:) = ra(1,:) * ( lainv(1) * la(1) )
                    if ( ycp_flag   ) ra(2,:) = ra(2,:) * ( lainv(2) * la(2) )
                    !^
                    lainv = 1.d0 / la
                end if

                ! velocity verlet method / move 1 / shear
                if ( cs_flag ) then
                    dstrain = strainv * dt + strainf * dt22
                    strain  = strain  + dstrain
                    strainv = strainv + strainf * dt2

                    ra(1,:) = ra(1,:) + dstrain * ra(2,:)

                end if

                ! velocity verlet method / force
                if ( nonnetwork_flag ) then
                    call calc_force( tcon, nbfire )
                else
                    call calc_force_spring( tcon, tnet )
                end if
                if ( cp_flag ) then
                    if ( boxcp_flag ) laf     = press  - boxp_set
                    if ( xcp_flag   ) laf(1)  = pressx - xp_set
                    if ( ycp_flag   ) laf(2)  = pressy - yp_set
                end if
                if ( cs_flag ) strainf = stress_set - stress

                ! velocity verlet method / move 2
                va = va + fa * dt2
                if ( cp_flag ) lav = lav + laf * dt2
                if ( cs_flag ) strainv = strainv + strainf * dt2

                ! fire
                cumn  = cumn + 1
                power = sum( fa * va ) + sum( laf*lav ) + strainv * strainf

                fn = sqrt( sum( fa**2 ) + sum(laf**2) + strainf**2 )
                vn = sqrt( sum( va**2 ) + sum(lav**2) + strainv**2 )

                onembeta  = 1.d0 - beta
                betavndfn = beta*vn/fn

                va      = onembeta * va      + betavndfn * fa
                lav     = onembeta * lav     + betavndfn * laf
                strainv = onembeta * strainv + betavndfn * strainf

                if ( power > 0.d0 .and. cumn > nmin ) then
                    dt   = min( dt*finc, dtmax )
                    beta = beta * fbeta
                end if

                if ( power < 0.d0 ) then
                    cumn    = 0
                    dt      = dt * fdec
                    beta    = beta0
                    va      = 0.d0
                    lav     = 0.d0
                    strainv = 0.d0
                end if

                temp = maxval( abs(fa) )
                if ( boxcp_flag ) temp = max( temp , abs( press  - boxp_set  ) )
                if ( xcp_flag   ) temp = max( temp , abs( pressx - xp_set ) )
                if ( ycp_flag   ) temp = max( temp , abs( pressy - yp_set ) )
                if ( cs_flag    ) temp = max( temp , abs( stress - stress_set ) )

                if ( temp < fmax .or. step == (stepmax-1)  ) then
                    !write(*,'(5es16.6)') 1.0, tcon.press, tcon.pressx, tcon.pressy, tcon.stress
                    exit
                end if

               !if( step == 1 ) write(*,*) 'step', '    dt         ', '              fmax   '

               !if( mod(step,10) == 0 ) then
               !    write(*,'(i6,3e16.6)') step, dt, temp, press
               !end if

            end do

        end associate

    end subroutine mini_fire_cp

end module
