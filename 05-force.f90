module mo_force
    use mo_syst
    use mo_config
    use mo_list
    use mo_network
    implicit none

contains

    subroutine calc_force( tcon, tnb )
        implicit none

        type(tpcon)  :: tcon
        type(tplist) :: tnb

        real(8), dimension(free) :: rai, raj, dra
        real(8) :: ri, rj, rij2, rij, dij, fr, wij, wili, wilix, wiliy
        integer :: iround(free), cory
        integer :: i, j, k, jj

        associate(                 &
            natom  => tcon.natom,  &
            radius => tcon.r,      &
            ra     => tcon.ra,     &
            fa     => tcon.fa,     &
            Ea     => tcon.Ea,     &
            la     => tcon.la,     &
            lainv  => tcon.lainv,  &
            strain => tcon.strain, &
            stress => tcon.stress, &
            press  => tcon.press,  &
            pressx => tcon.pressx, &
            pressy => tcon.pressy, &
            list   => tnb.list     &
            )

            Ea     = 0.d0
            fa     = 0.d0
            press  = 0.d0; pressx = 0.d0; pressy = 0.d0
            stress = 0.d0
            wili   = 0.d0; wilix  = 0.d0; wiliy  = 0.d0

            do i=1, natom

                rai = ra(:,i)
                ri  = radius(i)

                do jj=1, list(i).nbsum

                    j = list(i).nblist(jj)
                    iround = list(i).iround(:,jj)
                    cory = list(i).cory(jj)

                    raj = ra(:,j)
                    rj  = radius(j)

                    dra = raj - rai
                    dra(1) = dra(1) - cory * strain * la(free)

                    do k=1, free
                        dra(k) = dra(k) - iround(k) * la(k)
                    end do

                    rij2 = sum( dra**2 )
                    dij = ri + rj

                    if ( rij2 > dij**2 ) cycle

                    rij = sqrt( rij2 )

                    Ea = Ea + ( 1.d0 - rij/dij )**alpha/alpha

                    wij = (1.d0 - rij/dij)**(alpha-1) * rij / dij
                    wili = wili + wij

                    fr = wij / rij2

                    wilix = wilix + fr * dra(1)**2
                    wiliy = wiliy + fr * dra(2)**2

                    fa(:,j) = fa(:,j) + fr * dra
                    fa(:,i) = fa(:,i) - fr * dra

                    stress = stress - 2 * dra(1) * dra(2) * fr

                end do

            end do

            stress = stress * product(lainv) / free
            press  = wili   * product(lainv) / free
            pressx = wilix  * product(lainv)
            pressy = wiliy  * product(lainv)

        end associate

    end subroutine calc_force

    subroutine calc_force_withoutlist( tcon )
        implicit none

        ! para list
        type(tpcon)  :: tcon

        ! local
        real(8), dimension(free) :: rai, raj, dra
        real(8) :: ri, rj, rij2, rij, dij, fr, wij, wili, wilix, wiliy
        integer :: iround(free), cory
        integer :: i, j, k, jj

        associate(                 &
            natom  => tcon.natom,  &
            radius => tcon.r,      &
            ra     => tcon.ra,     &
            fa     => tcon.fa,     &
            Ea     => tcon.Ea,     &
            la     => tcon.la,     &
            lainv  => tcon.lainv,  &
            strain => tcon.strain, &
            stress => tcon.stress, &
            press  => tcon.press,  &
            pressx => tcon.pressx, &
            pressy => tcon.pressy  &
            )

            Ea     = 0.d0
            fa     = 0.d0
            press  = 0.d0; pressx = 0.d0; pressy = 0.d0
            stress = 0.d0
            wili   = 0.d0; wilix  = 0.d0; wiliy  = 0.d0

            do i=1, natom

                rai = ra(:,i)
                ri  = radius(i)

                do j=i+1, natom

                    raj = ra(:,j)
                    rj  = radius(j)

                    dra = raj - rai

                    cory = nint( dra(free) * lainv(free) )
                    dra(1) = dra(1) - strain * la(free) * cory

                    do k=1, free-1
                        iround(k) = nint( dra(k) * lainv(k) )
                    end do
                    iround(free) = cory

                    do k=1, free
                        dra(k) = dra(k) - iround(k) * la(k)
                    end do

                    rij2 = sum( dra**2 )
                    dij = ri + rj

                    if ( rij2 > dij**2 ) cycle

                    rij = sqrt( rij2 )

                    Ea = Ea + ( 1.d0 - rij/dij )**alpha/alpha

                    wij = (1.d0 - rij/dij)**(alpha-1) * rij / dij
                    wili = wili + wij

                    fr = wij / rij2

                    wilix = wilix + fr * dra(1)**2
                    wiliy = wiliy + fr * dra(2)**2

                    fa(:,j) = fa(:,j) + fr * dra
                    fa(:,i) = fa(:,i) - fr * dra

                    stress = stress - 2 * dra(1) * dra(2) * fr

                end do

            end do

            stress = stress * product(lainv) / free
            press  = wili   * product(lainv) / free
            pressx = wilix  * product(lainv)
            pressy = wiliy  * product(lainv)

        end associate

    end subroutine calc_force_withoutlist

    subroutine calc_force_spring( tcon, tnet )
        implicit none

        ! para list
        type(tpcon),     intent(inout) :: tcon
        type(tpnetwork), intent(in)    :: tnet

        ! local
        real(8), dimension(free) :: rai, raj, dra
        real(8) :: ri, rij2, rij, dij, l0, fr, wij, wili, wilix, wiliy, ks
        integer :: iround(free), cory
        integer :: i, j, k, ii, jj

        associate(                 &
            natom  => tcon.natom,  &
            radius => tcon.r,      &
            ra     => tcon.ra,     &
            fa     => tcon.fa,     &
            Ea     => tcon.Ea,     &
            la     => tcon.la,     &
            lainv  => tcon.lainv,  &
            strain => tcon.strain, &
            stress => tcon.stress, &
            press  => tcon.press,  &
            pressx => tcon.pressx, &
            pressy => tcon.pressy, &
            list   => tnet.sps,    &
            nlist  => tnet.nsps    &
            )

            Ea     = 0.d0
            fa     = 0.d0
            press  = 0.d0; pressx = 0.d0; pressy = 0.d0
            stress = 0.d0
            wili   = 0.d0; wilix  = 0.d0; wiliy  = 0.d0

            do ii=1, nlist

                i      = list(ii).i
                j      = list(ii).j
                cory   = list(ii).cory
                iround = list(ii).iround
                l0     = list(ii).l0
                ks     = list(ii).ks

                rai = ra(:,i)
                raj = ra(:,j)

                dra = raj - rai
                dra(1) = dra(1) - cory * strain * la(free)

                do k=1, free
                    dra(k) = dra(k) - iround(k) * la(k)
                end do

                rij2 = sum( dra**2 )
                rij  = sqrt( rij2 )

                Ea = Ea + 0.5d0 * ks * ( l0 - rij )**2

                wij  = ks * ( l0 - rij ) * rij
                wili = wili + wij

                fr = wij / rij2

                wilix = wilix + fr * dra(1)**2
                wiliy = wiliy + fr * dra(2)**2

                fa(:,j) = fa(:,j) + fr * dra
                fa(:,i) = fa(:,i) - fr * dra

                stress = stress - 2 * dra(1) * dra(2) * fr  ! 3d ?

            end do

            stress = stress * product(lainv) / free
            press  = wili   * product(lainv) / free
            pressx = wilix  * product(lainv)
            pressy = wiliy  * product(lainv)

        end associate
    end subroutine

end module
