program main
    use mo_syst
    use mo_var
    use mo_config
    use mo_list
    use mo_fire
    use mo_network
    use mo_disorder
    implicit none

    ! vars
    call testvar

    ! system
    call init_system( con, sets.natom, sets.phi )
    call gen_lattice_triangle( con )

    ! network
    call init_network( net, con )
    call make_network( net, con )
    !print*, maxval(net.sps(:).l0)

    ! fire
    call init_confire( confire, con )
    call check_system_force( confire, tnet = net )

    ! disorder
    call init_disorder( radisorder, con, sets.np )

   !do step=5, 48
    do step=1, 48
    do step2=1, 1, 2000

        eta = 0.01 * step
        temp1 = 1.d-5 * step2

        confire = con

        call make_disorder( confire, radisorder, eta )
       !call remake_network( net, confire )
        call change_k_spring( net, confire, opktan = 6.5/eta )
       !call change_k_spring_2( net, confire, opktan = 10.d0/eta )

        confire = con
       !call remake_network( net, confire )

        ! Ex, Mu_xy
        testp = 1.d-8
        confire2 = confire
        call mini_fire_cp( confire2, tnet=net, txp_set = testp, typ_set = 0.d0 )
        epl_xx = ( confire2.la(1) - confire.la(1) ) / confire.la(1)
        epl_yy = ( confire2.la(2) - confire.la(2) ) / confire.la(2)
        Mk_x = - testp / epl_xx
        mu_xy = - epl_yy / epl_xx

        ! Ey, Mu_xy
        testp = 1.d-8
        confire2 = confire
        call mini_fire_cp( confire2, tnet=net, txp_set = 0.d0, typ_set = testp )
        epl_xx = ( confire2.la(1) - confire.la(1) ) / confire.la(1)
        epl_yy = ( confire2.la(2) - confire.la(2) ) / confire.la(2)
        Mk_y = - testp / epl_yy
        mu_yx = - epl_xx / epl_yy

        ! Gs
        call check_system_force( confire, tnet = net )
        confire2 = confire
        confire2.strain = confire2.strain + testp
        call mini_fire_cv( confire2, tnet=net )
       !MG_s = ( confire2.stress - confire.stress ) / testp
        MG_s = confire2.Ea / ( 0.5d0 * product(confire2.la(1:free)) * testp**2 )

        ! Gxy
        confire2 = confire
        associate( la => confire2.la, lainv => confire2.lainv )
            la(1) = ( 1.d0 - testp ) * la(1)
            la(2) = ( 1.d0 + testp ) * la(2)
            lainv = 1.d0 / la
        end associate
        call mini_fire_cv( confire2, tnet=net )
        MG_xy = ( confire2.pressx - confire2.pressy ) / testp / 4.d0
       !MG_xy = confire2.Ea / ( 0.5d0 * product(confire2.la(1:free)) * testp**2 )

        write(*,'(9es16.6)') eta, temp1, Mk_x, Mk_y, mu_xy, mu_yx, MG_s, MG_xy, 0.5d0*(MG_s + MG_xy)

    end do
    end do


contains

    subroutine testvar
        implicit none

        sets.natom = 256
        sets.phi = 0.91d0
        sets.np = 202

    end subroutine testvar

end program main
