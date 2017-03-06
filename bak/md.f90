module mo_md
    use mo_syst
    use mo_config
    use mo_list
    use mo_network
    use mo_force
    implicit none

    type tpmdargs
        real(8) :: temper
        real(8) :: press
    end type

    type(tpmdargs) :: mdargs

    real(8), private, parameter :: dt   = 1.d-2
    real(8), private, parameter :: dt2  = 1.d-4
    real(8), private, parameter :: hdt  = 5.d-3
    real(8), private, parameter :: hdt2 = 5.d-5

    integer, private, save      :: cumstep = 0
    integer, private, parameter :: nscale = 33

contains
    
    subroutine init_md( ttemper, oppress )
        implicit none

        ! para list
        real(8), intent(in) :: ttemper
        real(8), intent(in), optional :: oppress

        mdargs.temper = ttemper
        
        if ( present( oppress ) ) then
            mdargs.press  = oppress
        end if

    end subroutine init_md

    subroutine md_nvt( tcon, tnb )
        implicit none

        ! para list
        type(tpcon),    intent(inout)        :: tcon
        type(tplist),   intent(in), optional :: tnb

        ! local
        real(8) :: chipxi

        associate(               &
            ra  => tcon.ra,      &
            va  => tcon.va,      &
            fa  => tcon.fa,      &
            Tk  => mdargs.temper &
            )

            ! velocity verlet : move 1
            va = va + fa * hdt
            ra = ra + va * dt

            ! velocity verlet : force
            call calc_force( tcon, tnb )
            
            ! constraint methond to control Temper
            chipxi = sum( fa*va ) / sum(va**2)
            fa = fa - chipxi * va

            ! velocity verlet : move 2
            va = va + fa * hdt

            cumstep = cumstep + 1
            if ( cumstep == nscale ) then
                cumstep = 0
                call scale_temper( tcon, Tk )
            end if

        end associate

    end subroutine

    subroutine scale_temp( tcon, tTk )
        implicit none

        ! para list
        type(tpcon), intent(inout) :: tcon
        real(8), intent(in) :: tTk

        ! local
        real(8) :: temper_now


        temper_now = 0.5d0 * sum( tcon.va**2 ) / ( tcon.natom * free )

        tcon.va = sqrt( tTk / temper_now ) * tcon.va

    end subroutine

end module

