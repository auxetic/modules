program main
    use mo_syst
    use mo_var
    use mo_config
    use mo_list
    use mo_fire
    use mo_md
    use mo_dynamic
    implicit none

    ! vars
    call testvar

    ! system
    call init_system( con, sets%natom, sets%phi )
    call gen_rand_config( con, sets%seed )

    ! list
    call init_list( nb, con )
    call make_list( nb, con )

    ! md
    call init_md( 1.d-3 )
    call pre_nvt( con )

    ! msd
    !call init_msd( con, msd1, 100, 200 )
    call init_vcorr( con, vcorr, 10, 200 )

    ! main
    do step=1, 10000
        if ( check_list( nb, con ) ) call make_list( nb, con )
        call md_nvt( con, nb )
    end do

    do step=1, 200000

        if ( check_list( nb, con ) ) call make_list( nb, con )

        call md_nvt( con, nb )

        if ( mod( step, 10000 ) == 0 ) print*, step
        !call calc_msd( con, msd1 )
        call calc_vcorr( con, vcorr )

    end do

    !call save_config_to( con, "./con.dat" )

    call endof_vcorr( vcorr )

    !do step=1, msd1%nhistmax
    !    print*, 1.d-2 * (step-1) * msd1%ndt, msd1%msd(step), msd1%alpha2(step)
    !end do
    do step=1, vcorr%nhistmax
        print*, 1.d-2 * (step-1) * vcorr%ndt, vcorr%vcorr(step)
    end do
contains

    subroutine testvar
        implicit none

        sets%natom = 128
        sets%phi = 0.84d0
        sets%seed = 202
    end subroutine testvar

end program main
