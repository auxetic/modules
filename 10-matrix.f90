module mo_mode
    use mo_syst
    use mo_config
    implicit none

    type tpmatrix
        real(8), allocatable, dimension(:,:) :: dymatrix, dymatrix0
        real(8), allocatable, dimension(:)   :: egdymatrix, pw
        real(8), allocatable, dimension(:,:,:) :: trimatrix
        real(8), allocatable, dimension(:)   :: modez
        integer :: ndim, mdim, natom
        logical :: boxflag   = .false.      ! can box length change
        logical :: xyflag    = .false.      ! can box length change seperately
        logical :: shearflag = .false.      ! can box be changed by shear
    contains
        procedure :: solve   => solve_mode
        procedure :: calc_pw => calc_pw
    end type

    type(tpmatrix) :: mode, mode1, mode2

contains

    subroutine init_mode( tmode, tcon, opboxflag, opxyflag, opshearflag )
        implicit none

        ! para list
        type(tpmatrix), intent(inout)        :: tmode
        type(tpcon),    intent(in)           :: tcon
        logical,        intent(in), optional :: opboxflag, opxyflag, opshearflag

        ! local

        associate(                       &
            natom     => tmode%natom,    &
            mdim      => tmode%mdim,     &
            ndim      => tmode%ndim,     &
            boxflag   => tmode%boxflag,  &
            xyflag    => tmode%xyflag,   &
            shearflag => tmode%shearflag &
            )

            boxflag   = .false.
            xyflag    = .false.
            shearflag = .false.
            if ( present(opboxflag) .and. present(opxyflag) ) then
                stop "boxflag and xyflag can not exist simultaneously"
            end if
            if ( present(opboxflag) )   boxflag   = opboxflag
            if ( present(opxyflag) )    xyflag    = opxyflag
            if ( present(opshearflag) ) shearflag = opshearflag

            natom = tcon%natom
            ndim  = free * tcon%natom
            mdim  = ndim

            if ( boxflag   ) mdim = mdim + 1
            if ( xyflag    ) mdim = mdim + free
            if ( shearflag ) mdim = mdim + 1

            if ( .not. allocated(tmode%dymatrix) ) then
                allocate( tmode%dymatrix(mdim,mdim), tmode%egdymatrix(mdim) )
               !allocate( tmode%dymatrix0(mdim,mdim) )
               !allocate( tmode%trimatrix(mdim,mdim,mdim) )
               !allocate( tmode%modez(mdim) )
            end if

        end associate
    end subroutine

    subroutine make_dymatrix( tmode, tcon )
        implicit none

        ! para list
        type(tpmatrix), intent(inout) :: tmode
        type(tpcon),    intent(in)    :: tcon

        ! local
        integer :: i, j
        real(8) :: dra(free), rij, rij2, rij3, dij, xij, yij
        real(8) :: vr, vrr, rx, ry, rxx, ryy, rxy
        real(8) :: mij(free,free)   ! \p^2v / ( \pxi * \pxj ) ....

        associate(                      &
            dymatrix => tmode%dymatrix, &
            radius   => tcon%r,         &
            natom    => tcon%natom      &
            )

            dymatrix = 0.d0
        
            do i=1, natom-1
                do j=i+1, natom

                    dra  = calc_dra( tcon, i, j )
                    rij2 = sum(dra**2)

                    dij = radius(i) + radius(j)

                    if ( rij2 > dij**2 ) cycle

                    rij  = sqrt(rij2)
                    rij3 = rij2*rij

                    ! \partial V_ij / \partial r_ij
                    vr = - ( 1.d0 - rij/dij )**(alpha-1) / dij
                    ! \partial^2 V_ij / \partial r_ij^2
                    vrr = ( alpha-1) * ( 1.d0 - rij/dij )**(alpha-2) / dij**2

                    xij = dra(1)
                    yij = dra(2)

                    ! \partial r_ij / \partial x_j
                    rx = xij / rij
                    ry = yij / rij

                    ! \partial r_ij / [ \partial x_j \partial x_j ]
                    rxx = 1.d0/rij - xij**2 /rij3
                    ryy = 1.d0/rij - yij**2 /rij3
                    rxy =          - xij*yij/rij3
                    
                    ! \partial^2 V_ij / [ \partial x_i \partial x_j ]
                    ! = - \partial V_ij / [ \partial x_j \partial x_j ]
                    ! =   - [ Vrr * rx**2   + vr * rxx ]
                    ! yy: - [ Vrr * rx**2 + Vr * ryy ]
                    ! xy: - [ Vrr * rx*ry + Vr * rxy ]
                    mij(1,1) = - ( vrr*rx**2 + vr*rxx )
                    mij(2,2) = - ( vrr*ry**2 + vr*ryy )
                    mij(1,2) = - ( vrr*rx*ry + vr*rxy )
                    mij(2,1) = mij(1,2)

                    dymatrix(free*(i-1)+1:free*i, free*(i-1)+1:free*i ) = dymatrix(free*(i-1)+1:free*i, free*(i-1)+1:free*i ) - mij
                    dymatrix(free*(j-1)+1:free*j, free*(j-1)+1:free*j ) = dymatrix(free*(j-1)+1:free*j, free*(j-1)+1:free*j ) - mij
                    dymatrix(free*(i-1)+1:free*i, free*(j-1)+1:free*j ) = mij
                    dymatrix(free*(j-1)+1:free*j, free*(i-1)+1:free*i ) = mij

                end do
            end do

        end associate
    end subroutine

    subroutine make_dymatrix_net( tmode, tcon, tnet )
        use mo_network
        implicit none

        ! para list
        type(tpmatrix), intent(inout) :: tmode
        type(tpcon), intent(in) :: tcon
        type(tpnetwork), intent(in) :: tnet

        ! local
        integer :: i, j, ii
        real(8) :: dra(free), rij, rij2, rij3, dij, xj, yj
        real(8) :: l0, ks
        real(8) :: vr, vrr, rx, ry, rxx, ryy, rxy
        real(8) :: mij(free,free)   ! \p^2v / ( \pxi * \pxj ) ....

        associate(                      &
            dymatrix => tmode%dymatrix, &
            radius   => tcon%r,         &
            natom    => tcon%natom,     &
            list     => tnet%sps,       &
            nlist    => tnet%nsps       &
            )

            dymatrix = 0.d0
        
            do ii=1, nlist

                i  = list(ii)%i
                j  = list(ii)%j
                l0 = list(ii)%l0
                ks = list(ii)%ks

                dra = calc_dra( tcon, i, j )
                rij2 = sum(dra**2)

                dij = radius(i) + radius(j)

                rij  = sqrt(rij2)
                rij3 = rij2*rij

                ! \partial V_ij / \partial r_ij
               !vr = - ( 1.d0 - rij/dij )**(alpha-1) / dij
                vr = 0.d0
                ! \partial^2 V_ij / \partial r_ij^2
               !vrr = ( alpha-1) * ( 1.d0 - rij/dij )**(alpha-2) / dij**2
                vrr = ks

                xj = dra(1)
                yj = dra(2)

                ! \partial r_ij / \partial x_j
                rx = xj / rij
                ry = yj / rij

                ! \partial r_ij / [ \partial x_j \partial x_j ]
                rxx = 1.d0/rij - xj**2/rij3
                ryy = 1.d0/rij - yj**2/rij3
                rxy =          - xj*yj/rij3
                
                ! \partial^2 V_ij / [ \partial x_i \partial x_j ]
                ! = - \partial V_ij / [ \partial x_j \partial x_j ]
                ! =   - [ Vrr * rx**2   + vr * rxx ]
                ! yy: - [ Vrr * rx**2 + Vr * ryy ]
                ! xy: - [ Vrr * rx*ry + Vr * rxy ]
                mij(1,1) = - ( vrr*rx**2 + vr*rxx )
                mij(2,2) = - ( vrr*ry**2 + vr*ryy )
                mij(1,2) = - ( vrr*rx*ry + vr*rxy )
                mij(2,1) = mij(1,2)

                dymatrix(free*(i-1)+1:free*i, free*(i-1)+1:free*i ) = dymatrix(free*(i-1)+1:free*i, free*(i-1)+1:free*i ) - mij
                dymatrix(free*(j-1)+1:free*j, free*(j-1)+1:free*j ) = dymatrix(free*(j-1)+1:free*j, free*(j-1)+1:free*j ) - mij
                dymatrix(free*(i-1)+1:free*i, free*(j-1)+1:free*j ) = mij
                dymatrix(free*(j-1)+1:free*j, free*(i-1)+1:free*i ) = mij

            end do

        end associate
    end subroutine

    subroutine make_dymatrix_free_boundary( tmode, tcon )
        implicit none

        ! para list
        type(tpmatrix), intent(inout) :: tmode
        type(tpcon), intent(in) :: tcon

        ! local
        integer :: i, j
        real(8) :: dra(free), rij, rij2, rij3, dij, xj, yj
        real(8) :: vr, vrr, rx, ry, rxx, ryy, rxy
        real(8) :: mij(free,free)   ! \p^2v / ( \pxi * \pxj ) ....

        associate(                      &
            dymatrix => tmode%dymatrix, &
            radius   => tcon%r,         &
            natom    => tcon%natom      &
            )

            dymatrix = 0.d0
        
            do i=1, natom-1
                do j=i+1, natom

                    dra = calc_dra( tcon, i, j )
                    rij2 = sum(dra**2)

                    dij = radius(i) + radius(j)

                    if ( rij2 > dij**2 ) cycle

                    rij  = sqrt(rij2)
                    rij3 = rij2*rij

                    ! \partial V_ij / \partial r_ij
                    vr = - ( 1.d0 - rij/dij )**(alpha-1) / dij
                    ! \partial^2 V_ij / \partial r_ij^2
                    vrr = ( alpha-1) * ( 1.d0 - rij/dij )**(alpha-2) / dij**2

                    xj = dra(1)
                    yj = dra(2)

                    ! \partial r_ij / \partial x_j
                    rx = xj / rij
                    ry = yj / rij

                    ! \partial r_ij / [ \partial x_j \partial x_j ]
                    rxx = 1.d0/rij - xj**2/rij3
                    ryy = 1.d0/rij - yj**2/rij3
                    rxy =          - xj*yj/rij3
                    
                    ! \partial^2 V_ij / [ \partial x_i \partial x_j ]
                    ! = - \partial V_ij / [ \partial x_j \partial x_j ]
                    ! =   - [ Vrr * rx**2   + vr * rxx ]
                    ! yy: - [ Vrr * rx**2 + Vr * ryy ]
                    ! xy: - [ Vrr * rx*ry + Vr * rxy ]
                    mij(1,1) = - ( vrr*rx**2 + vr*rxx )
                    mij(2,2) = - ( vrr*ry**2 + vr*ryy )
                    mij(1,2) = - ( vrr*rx*ry + vr*rxy )
                    mij(2,1) = mij(1,2)

                    ! ij
                    dymatrix(free*(i-1)+1:free*i, free*(i-1)+1:free*i ) = dymatrix(free*(i-1)+1:free*i, free*(i-1)+1:free*i ) - mij
                    dymatrix(free*(j-1)+1:free*j, free*(j-1)+1:free*j ) = dymatrix(free*(j-1)+1:free*j, free*(j-1)+1:free*j ) - mij
                    dymatrix(free*(i-1)+1:free*i, free*(j-1)+1:free*j ) = mij
                    dymatrix(free*(j-1)+1:free*j, free*(i-1)+1:free*i ) = mij

                    ! compress 4+4
                    ! i ex(y)
                    dymatrix(free*(i-1)+1:free*i, free*natom+1) = dymatrix(free*(i-1)+1:free*i, free*natom+1) + mij(:,1)*dra(1)
                    dymatrix(free*natom+1, free*(i-1)+1:free*i) = dymatrix(free*natom+1, free*(i-1)+1:free*i) + mij(1,:)*dra(1)
                    dymatrix(free*(i-1)+1:free*i, free*natom+2) = dymatrix(free*(i-1)+1:free*i, free*natom+2) + mij(:,2)*dra(2)
                    dymatrix(free*natom+2, free*(i-1)+1:free*i) = dymatrix(free*natom+2, free*(i-1)+1:free*i) + mij(2,:)*dra(2)
                    ! j ex(y)
                    dymatrix(free*(j-1)+1:free*j, free*natom+1) = dymatrix(free*(j-1)+1:free*j, free*natom+1) + mij(:,1)*dra(1)
                    dymatrix(free*natom+1, free*(j-1)+1:free*j) = dymatrix(free*natom+1, free*(j-1)+1:free*j) + mij(1,:)*dra(1)
                    dymatrix(free*(j-1)+1:free*j, free*natom+2) = dymatrix(free*(j-1)+1:free*j, free*natom+2) + mij(:,2)*dra(2)
                    dymatrix(free*natom+2, free*(j-1)+1:free*j) = dymatrix(free*natom+2, free*(j-1)+1:free*j) + mij(2,:)*dra(2)
                    ! ex ey
                    dymatrix(free*natom+1,free*natom+1) = dymatrix(free*natom+1,free*natom+1) - mij(1,1)*dra(1)**2
                    dymatrix(free*natom+2,free*natom+2) = dymatrix(free*natom+2,free*natom+2) - mij(2,2)*dra(2)**2
                    dymatrix(free*natom+1,free*natom+2) = dymatrix(free*natom+1,free*natom+2) - mij(1,2)*dra(1)*dra(2)
                    dymatrix(free*natom+2,free*natom+1) = dymatrix(free*natom+2,free*natom+1) - mij(2,1)*dra(1)*dra(2)

                    ! shear 2+1 + 4
                    ! i gamma
                    dymatrix(free*(i-1)+1:free*i, free*natom+3) = dymatrix(free*(i-1)+1:free*i, free*natom+3) + mij(:,1)*dra(2)
                    dymatrix(free*natom+3, free*(i-1)+1:free*i) = dymatrix(free*natom+3, free*(i-1)+1:free*i) + mij(1,:)*dra(2)
                    ! j gamma
                    dymatrix(free*(j-1)+1:free*j, free*natom+3) = dymatrix(free*(j-1)+1:free*j, free*natom+3) + mij(:,1)*dra(2)
                    dymatrix(free*natom+3, free*(j-1)+1:free*j) = dymatrix(free*natom+3, free*(j-1)+1:free*j) + mij(1,:)*dra(2)
                    ! gamma gamma
                    dymatrix(free*natom+3,free*natom+3) = dymatrix(free*natom+3,free*natom+3) - mij(1,1)*dra(2)**2

                    ! ex gamma
                    dymatrix(free*natom+1,free*natom+3) = dymatrix(free*natom+1,free*natom+3) - mij(1,1)*dra(1)*dra(2)
                    dymatrix(free*natom+3,free*natom+1) = dymatrix(free*natom+3,free*natom+1) - mij(1,1)*dra(1)*dra(2)
                    ! ey gamma
                    dymatrix(free*natom+2,free*natom+3) = dymatrix(free*natom+2,free*natom+3) - mij(1,2)*dra(2)**2
                    dymatrix(free*natom+3,free*natom+2) = dymatrix(free*natom+3,free*natom+2) - mij(2,1)*dra(2)**2

                end do
            end do

        end associate
    end subroutine

    subroutine make_dymatrix_free_boundary_net_3( tmode, tcon, tnet )
        use mo_network
        implicit none

        ! para list
        type(tpmatrix), intent(inout) :: tmode
        type(tpcon), intent(in) :: tcon
        type(tpnetwork), intent(in) :: tnet

        ! local
        integer :: i, j, ii
        real(8) :: dra(free), rij, rij2, rij3, dij, xj, yj
        real(8) :: vr, vrr, rx, ry, rxx, ryy, rxy
        real(8) :: mij(free,free)   ! \p^2v / ( \pxi * \pxj ) ....
        real(8) :: l0, ks

        associate(                      &
            dymatrix => tmode%dymatrix, &
            radius   => tcon%r,         &
            natom    => tcon%natom,     &
            list     => tnet%sps,       &
            nlist    => tnet%nsps       &
            )

            dymatrix = 0.d0
        
            do ii=1, nlist

                i  = list(ii)%i
                j  = list(ii)%j
                l0 = list(ii)%l0
                ks = list(ii)%ks

                dra = calc_dra( tcon, i, j )
                rij2 = sum(dra**2)

                rij  = sqrt(rij2)
                rij3 = rij2*rij

                ! \partial V_ij / \partial r_ij
                vr = 0.d0
                ! \partial^2 V_ij / \partial r_ij^2
                vrr = ks

                xj = dra(1)
                yj = dra(2)

                ! \partial r_ij / \partial x_j
                rx = xj / rij
                ry = yj / rij

                ! \partial r_ij / [ \partial x_j \partial x_j ]
                rxx = 1.d0/rij - xj**2/rij3
                ryy = 1.d0/rij - yj**2/rij3
                rxy =          - xj*yj/rij3
                
                ! \partial^2 V_ij / [ \partial x_i \partial x_j ]
                ! = - \partial V_ij / [ \partial x_j \partial x_j ]
                ! =   - [ Vrr * rx**2   + vr * rxx ]
                ! yy: - [ Vrr * rx**2 + Vr * ryy ]
                ! xy: - [ Vrr * rx*ry + Vr * rxy ]
                mij(1,1) = - ( vrr*rx**2 + vr*rxx )
                mij(2,2) = - ( vrr*ry**2 + vr*ryy )
                mij(1,2) = - ( vrr*rx*ry + vr*rxy )
                mij(2,1) = mij(1,2)

                ! ij
                dymatrix(free*(i-1)+1:free*i, free*(i-1)+1:free*i ) = dymatrix(free*(i-1)+1:free*i, free*(i-1)+1:free*i ) - mij
                dymatrix(free*(j-1)+1:free*j, free*(j-1)+1:free*j ) = dymatrix(free*(j-1)+1:free*j, free*(j-1)+1:free*j ) - mij
                dymatrix(free*(i-1)+1:free*i, free*(j-1)+1:free*j ) = mij
                dymatrix(free*(j-1)+1:free*j, free*(i-1)+1:free*i ) = mij

                ! compress 4+4
                dymatrix(free*(i-1)+1:free*i, free*natom+1) = dymatrix(free*(i-1)+1:free*i, free*natom+1) + mij(:,1)*dra(1)
                dymatrix(free*natom+1, free*(i-1)+1:free*i) = dymatrix(free*natom+1, free*(i-1)+1:free*i) + mij(1,:)*dra(1)
                dymatrix(free*(i-1)+1:free*i, free*natom+2) = dymatrix(free*(i-1)+1:free*i, free*natom+2) + mij(:,2)*dra(2)
                dymatrix(free*natom+2, free*(i-1)+1:free*i) = dymatrix(free*natom+2, free*(i-1)+1:free*i) + mij(2,:)*dra(2)
                !
                dymatrix(free*(j-1)+1:free*j, free*natom+1) = dymatrix(free*(j-1)+1:free*j, free*natom+1) - mij(:,1)*dra(1)
                dymatrix(free*natom+1, free*(j-1)+1:free*j) = dymatrix(free*natom+1, free*(j-1)+1:free*j) - mij(1,:)*dra(1)
                dymatrix(free*(j-1)+1:free*j, free*natom+2) = dymatrix(free*(j-1)+1:free*j, free*natom+2) - mij(:,2)*dra(2)
                dymatrix(free*natom+2, free*(j-1)+1:free*j) = dymatrix(free*natom+2, free*(j-1)+1:free*j) - mij(2,:)*dra(2)
                !
                dymatrix(free*natom+1,free*natom+1) = dymatrix(free*natom+1,free*natom+1) - mij(1,1)*dra(1)**2
                dymatrix(free*natom+2,free*natom+2) = dymatrix(free*natom+2,free*natom+2) - mij(2,2)*dra(2)**2
                dymatrix(free*natom+1,free*natom+2) = dymatrix(free*natom+1,free*natom+2) - mij(1,2)*dra(1)*dra(2)
                dymatrix(free*natom+2,free*natom+1) = dymatrix(free*natom+2,free*natom+1) - mij(2,1)*dra(1)*dra(2)

                ! shear 2+1 + 4
                dymatrix(free*(i-1)+1:free*i, free*natom+3) = dymatrix(free*(i-1)+1:free*i, free*natom+3) + mij(:,1)*dra(2)
                dymatrix(free*natom+3, free*(i-1)+1:free*i) = dymatrix(free*natom+3, free*(i-1)+1:free*i) + mij(1,:)*dra(2)
                !
                dymatrix(free*(j-1)+1:free*j, free*natom+3) = dymatrix(free*(j-1)+1:free*j, free*natom+3) - mij(:,1)*dra(2)
                dymatrix(free*natom+3, free*(j-1)+1:free*j) = dymatrix(free*natom+3, free*(j-1)+1:free*j) - mij(1,:)*dra(2)
                !
                dymatrix(free*natom+3,free*natom+3) = dymatrix(free*natom+3,free*natom+3) - mij(1,1)*dra(2)**2

                dymatrix(free*natom+1,free*natom+3) = dymatrix(free*natom+1,free*natom+3) - mij(1,1)*dra(1)*dra(2)
                dymatrix(free*natom+3,free*natom+1) = dymatrix(free*natom+3,free*natom+1) - mij(1,1)*dra(1)*dra(2)

                dymatrix(free*natom+2,free*natom+3) = dymatrix(free*natom+2,free*natom+3) - mij(1,2)*dra(2)**2
                dymatrix(free*natom+3,free*natom+2) = dymatrix(free*natom+3,free*natom+2) - mij(2,1)*dra(2)**2

            end do

        end associate
    end subroutine

    subroutine make_dymatrix_free_boundary_net_2_box( tmode, tcon, tnet )
        use mo_network
        implicit none

        ! para list
        type(tpmatrix), intent(inout) :: tmode
        type(tpcon), intent(in) :: tcon
        type(tpnetwork), intent(in) :: tnet

        ! local
        integer :: i, j, ii
        real(8) :: dra(free), rij, rij2, rij3, dij, xj, yj
        real(8) :: l0, ks
        real(8) :: vr, vrr, rx, ry, rxx, ryy, rxy
        real(8) :: rex, rey, rexx, rexy, reyx, reyy, rexex, reyey, rexey
        real(8) :: mij(free,free)   ! \p^2v / ( \pxi * \pxj ) ....
        real(8) :: miex(free), miey(free), mee(free,free)

        associate(                      &
            dymatrix => tmode%dymatrix, &
            radius   => tcon%r,         &
            natom    => tcon%natom,     &
            list     => tnet%sps,       &
            nlist    => tnet%nsps       &
            )

            dymatrix = 0.d0
        
            do ii=1, nlist

                i  = list(ii)%i
                j  = list(ii)%j
                l0 = list(ii)%l0
                ks = list(ii)%ks

                dra = calc_dra( tcon, i, j )
                rij2 = sum(dra**2)

                dij = radius(i) + radius(j)

                rij  = sqrt(rij2)
                rij3 = rij2*rij

                ! \partial V_ij / \partial r_ij
                vr = 0.d0
               !vr = - ( 1.d0 - rij/dij )**(alpha-1) / dij
                ! \partial^2 V_ij / \partial r_ij^2
               !vrr = ( alpha-1) * ( 1.d0 - rij/dij )**(alpha-2) / dij**2
                vrr = ks

                xj = dra(1)
                yj = dra(2)

                ! \partial r_ij / \partial x_j
                rx = xj / rij
                ry = yj / rij
                ! \partial r_ij / \partial \epsilon_x = rx * x_ij
                rex = rx * xj
                rey = ry * yj

                ! \partial^2 r_ij / [ \partial x_j \partial x_j ]
                rxx = 1.d0/rij - xj**2/rij3
                ryy = 1.d0/rij - yj**2/rij3
                rxy =          - xj*yj/rij3
                ! \partial^2 r_ij / [ \partial x_j \partial \epsilon_x ]
                rexx = rxx * xj
                rexy = rxy * xj
                reyx = rxy * yj
                reyy = ryy * yj
                !
                rexex = rxx * xj**2
                reyey = ryy * yj**2
                rexey = rxy * xj*yj

                ! \partial^2 V_ij / [ \partial x_i \partial x_j ]
                ! = - \partial V_ij / [ \partial x_j \partial x_j ]
                ! =   - [ Vrr * rx**2   + vr * rxx ]
                ! yy: - [ Vrr * rx**2 + Vr * ryy ]
                ! xy: - [ Vrr * rx*ry + Vr * rxy ]
                mij(1,1) = - ( vrr*rx**2 + vr*rxx )
                mij(2,2) = - ( vrr*ry**2 + vr*ryy )
                mij(1,2) = - ( vrr*rx*ry + vr*rxy )
                mij(2,1) = mij(1,2)

                miex(1) = - ( vrr*rx*rex + vr*rexx )
                miex(2) = - ( vrr*ry*rex + vr*rexy )
                miey(1) = - ( vrr*rx*rey + vr*reyx )
                miey(2) = - ( vrr*ry*rey + vr*reyy )
                
                mee(1,1) = vrr*rex*rex + vr*rexex
                mee(2,2) = vrr*rey*rey + vr*reyey
                mee(1,2) = vrr*rex*rey + vr*rexey
                mee(2,1) = mee(1,2)

                ! ij
                dymatrix(free*(i-1)+1:free*i, free*(i-1)+1:free*i ) = dymatrix(free*(i-1)+1:free*i, free*(i-1)+1:free*i ) - mij
                dymatrix(free*(j-1)+1:free*j, free*(j-1)+1:free*j ) = dymatrix(free*(j-1)+1:free*j, free*(j-1)+1:free*j ) - mij
                dymatrix(free*(i-1)+1:free*i, free*(j-1)+1:free*j ) = mij
                dymatrix(free*(j-1)+1:free*j, free*(i-1)+1:free*i ) = mij

                dymatrix(free*(i-1)+1:free*i, free*natom+1) = dymatrix(free*(i-1)+1:free*i, free*natom+1) + miex
                dymatrix(free*natom+1, free*(i-1)+1:free*i) = dymatrix(free*natom+1, free*(i-1)+1:free*i) + miex
                dymatrix(free*(i-1)+1:free*i, free*natom+2) = dymatrix(free*(i-1)+1:free*i, free*natom+2) + miey
                dymatrix(free*natom+2, free*(i-1)+1:free*i) = dymatrix(free*natom+2, free*(i-1)+1:free*i) + miey

                dymatrix(free*(j-1)+1:free*j, free*natom+1) = dymatrix(free*(j-1)+1:free*j, free*natom+1) - miex
                dymatrix(free*natom+1, free*(j-1)+1:free*j) = dymatrix(free*natom+1, free*(j-1)+1:free*j) - miex
                dymatrix(free*(j-1)+1:free*j, free*natom+2) = dymatrix(free*(j-1)+1:free*j, free*natom+2) - miey
                dymatrix(free*natom+2, free*(j-1)+1:free*j) = dymatrix(free*natom+2, free*(j-1)+1:free*j) - miey

                dymatrix(free*natom+1:free*natom+2,free*natom+1:free*natom+2) = dymatrix(free*natom+1:free*natom+2,free*natom+1:free*natom+2) + mee

            end do

        end associate
    end subroutine

    subroutine make_dymatrix_free_boundary_net_2( tmode, tcon, tnet )
        use mo_network
        implicit none

        ! para list
        type(tpmatrix), intent(inout) :: tmode
        type(tpcon), intent(in) :: tcon
        type(tpnetwork), intent(in) :: tnet

        ! local
        integer :: i, j, ii
        real(8) :: dra(free), rij, rij2, rij3, dij, xj, yj
        real(8) :: l0, ks
        real(8) :: vr, vrr, rx, ry, rxx, ryy, rxy
        real(8) :: rex, rey, rexx, rexy, reyx, reyy, rexex, reyey, rexey
        real(8) :: mij(free,free)   ! \p^2v / ( \pxi * \pxj ) ....
        real(8) :: miex(free), miey(free), mee(free,free)

        associate(                      &
            dymatrix => tmode%dymatrix, &
            radius   => tcon%r,         &
            natom    => tcon%natom,     &
            list     => tnet%sps,       &
            nlist    => tnet%nsps       &
            )

            dymatrix = 0.d0
        
            do ii=1, nlist

                i  = list(ii)%i
                j  = list(ii)%j
                l0 = list(ii)%l0
                ks = list(ii)%ks

                dra = calc_dra( tcon, i, j )
                rij2 = sum(dra**2)

                dij = radius(i) + radius(j)

                rij  = sqrt(rij2)
                rij3 = rij2*rij

                ! \partial V_ij / \partial r_ij
                vr = 0.d0
               !vr = - ( 1.d0 - rij/dij )**(alpha-1) / dij
                ! \partial^2 V_ij / \partial r_ij^2
               !vrr = ( alpha-1) * ( 1.d0 - rij/dij )**(alpha-2) / dij**2
                vrr = ks

                xj = dra(1)
                yj = dra(2)

                ! \partial r_ij / \partial x_j
                rx = xj / rij
                ry = yj / rij
                ! \partial r_ij / \partial \epsilon_x = rx * x_ij
                rex = rx * xj
                rey = ry * yj

                ! \partial^2 r_ij / [ \partial x_j \partial x_j ]
                rxx = 1.d0/rij - xj**2/rij3
                ryy = 1.d0/rij - yj**2/rij3
                rxy =          - xj*yj/rij3
                ! \partial^2 r_ij / [ \partial x_j \partial \epsilon_x ]
                rexx = rxx * xj + rx
                rexy = rxy * xj
                reyx = rxy * yj
                reyy = ryy * yj + ry
                !
                rexex = rxx * xj**2 + rx * xj
                reyey = ryy * yj**2 + ry * yj

                ! \partial^2 V_ij / [ \partial x_i \partial x_j ]
                ! = - \partial V_ij / [ \partial x_j \partial x_j ]
                ! =   - [ Vrr * rx**2   + vr * rxx ]
                ! yy: - [ Vrr * rx**2 + Vr * ryy ]
                ! xy: - [ Vrr * rx*ry + Vr * rxy ]
                mij(1,1) = - ( vrr*rx**2 + vr*rxx )
                mij(2,2) = - ( vrr*ry**2 + vr*ryy )
                mij(1,2) = - ( vrr*rx*ry + vr*rxy )
                mij(2,1) = mij(1,2)

                miex(1) = - ( vrr*rx*rex + vr*rexx )
                miex(2) = - ( vrr*ry*rex + vr*rexy )
                miey(1) = - ( vrr*rx*rey + vr*reyx )
                miey(2) = - ( vrr*ry*rey + vr*reyy )
                
                mee(1,1) = vrr*rex*rex + vr*rexex
                mee(2,2) = vrr*rey*rey + vr*reyey
                mee(1,2) = vrr*rex*rey + vr*rexey
                mee(2,1) = mee(1,2)

                ! ij
                dymatrix(free*(i-1)+1:free*i, free*(i-1)+1:free*i ) = dymatrix(free*(i-1)+1:free*i, free*(i-1)+1:free*i ) - mij
                dymatrix(free*(j-1)+1:free*j, free*(j-1)+1:free*j ) = dymatrix(free*(j-1)+1:free*j, free*(j-1)+1:free*j ) - mij
                dymatrix(free*(i-1)+1:free*i, free*(j-1)+1:free*j ) = mij
                dymatrix(free*(j-1)+1:free*j, free*(i-1)+1:free*i ) = mij

                dymatrix(free*(i-1)+1:free*i, free*natom+1) = dymatrix(free*(i-1)+1:free*i, free*natom+1) + miex
                dymatrix(free*natom+1, free*(i-1)+1:free*i) = dymatrix(free*natom+1, free*(i-1)+1:free*i) + miex
                dymatrix(free*(i-1)+1:free*i, free*natom+2) = dymatrix(free*(i-1)+1:free*i, free*natom+2) + miey
                dymatrix(free*natom+2, free*(i-1)+1:free*i) = dymatrix(free*natom+2, free*(i-1)+1:free*i) + miey

                dymatrix(free*(j-1)+1:free*j, free*natom+1) = dymatrix(free*(j-1)+1:free*j, free*natom+1) - miex
                dymatrix(free*natom+1, free*(j-1)+1:free*j) = dymatrix(free*natom+1, free*(j-1)+1:free*j) - miex
                dymatrix(free*(j-1)+1:free*j, free*natom+2) = dymatrix(free*(j-1)+1:free*j, free*natom+2) - miey
                dymatrix(free*natom+2, free*(j-1)+1:free*j) = dymatrix(free*natom+2, free*(j-1)+1:free*j) - miey

                dymatrix(free*natom+1:free*natom+2,free*natom+1:free*natom+2) = dymatrix(free*natom+1:free*natom+2,free*natom+1:free*natom+2) + mee

            end do

        end associate
    end subroutine

    subroutine make_trimatrix( tmode, tcon )
        implicit none

        ! para list
        type(tpmatrix), intent(inout) :: tmode
        type(tpcon), intent(in) :: tcon

        ! local
        integer :: i, j
        real(8) :: dra(free), rij, rij2, rij3, rij5, dij, xi, xj, yi, yj
        real(8) :: vr, vrr, vrrr
        real(8) :: rx, ry 
        real(8) :: rxx, ryy, rxy
        real(8) :: rxxx, rxxy, rxyy, ryyy
        real(8) :: mijk(free,free,free)   ! \p^2v / ( \pxi * \pxj ) ....

        associate(                        &
            trimatrix => tmode%trimatrix, &
            radius    => tcon%r,          &
            natom     => tcon%natom       &
            )

            trimatrix = 0.d0
        
            do i=1, natom-1
                do j=i+1, natom

                    dra = calc_dra( tcon, i, j )
                    rij2 = sum(dra**2)

                    dij = radius(i) + radius(j)

                    if ( rij2 > dij**2 ) cycle

                    rij  = sqrt(rij2)
                    rij3 = rij2*rij
                    rij5 = rij3*rij2

                    !vr = - ( 1.d0 - rij/dij )**(alpha-1) / dij
                    !vrr = ( alpha-1) * ( 1.d0 - rij/dij )**(alpha-2) / dij**2
                    ! alpha = 2
                    vr   = - ( 1.d0 - rij/dij ) / dij
                    vrr  = 1.d0 / dij**2
                    !vrrr = 0.d0

                    xj = dra(1); yj = dra(2)
                    xi =-dra(1); yi =-dra(2)

                    rx = xj / rij
                    ry = yj / rij

                    rxx = 1.d0/rij - rx**2/rij3
                    ryy = 1.d0/rij - ry**2/rij3
                    rxy =          - rx*ry/rij3

                    ! \pr / ( \pxj \pxj \pxj )
                    rxxx = 3.d0*xj/rij3 + xj**3/rij5
                    ryyy = 3.d0*yj/rij3 + yj**3/rij5
                    rxxy = -yj/rij3 + 3.d0*xj**2*yj/rij5
                    rxyy = -xj/rij3 + 3.d0*xj*yj**2/rij5

                    ! \pv / ( \pxi \pxj \pxj )
                    mijk(1,1,1) = -( vrr*3.d0*rxx*rx - vr*rxxx )
                    mijk(1,1,2) = -( vrr*(2.d0*rxx*ry+rx*rxy) - vr*rxxy )
                    mijk(1,2,1) = mijk(1,1,2)
                    mijk(1,2,2) = -( vrr*(2.d0*ryy*rx+ry*rxy) - vr*rxyy )
                    mijk(2,1,1) = -mijk(1,1,2)
                    mijk(2,1,2) = -mijk(1,2,2)
                    mijk(2,2,1) = -mijk(1,2,2)
                    mijk(2,2,2) = -( vrr*3.d0*ryy*ry - vr*ryyy )

                    ! i
                    trimatrix(free*(i-1)+1:free*i, free*(j-1)+1:free*j, free*(j-1)+1:free*j ) =  mijk
                    trimatrix(free*(i-1)+1:free*i, free*(j-1)+1:free*j, free*(i-1)+1:free*i ) = -mijk
                    trimatrix(free*(i-1)+1:free*i, free*(i-1)+1:free*i, free*(j-1)+1:free*j ) = -mijk
                    ! j
                    trimatrix(free*(j-1)+1:free*j, free*(j-1)+1:free*j, free*(i-1)+1:free*i ) =  mijk
                    trimatrix(free*(j-1)+1:free*j, free*(i-1)+1:free*i, free*(i-1)+1:free*i ) = -mijk
                    trimatrix(free*(j-1)+1:free*j, free*(i-1)+1:free*i, free*(j-1)+1:free*j ) =  mijk
                    ! iii, jjj
                    trimatrix(free*(i-1)+1:free*i, free*(i-1)+1:free*i, free*(i-1)+1:free*i ) = &
                        trimatrix(free*(i-1)+1:free*i, free*(i-1)+1:free*i, free*(i-1)+1:free*i ) + mijk
                    trimatrix(free*(j-1)+1:free*j, free*(j-1)+1:free*j, free*(j-1)+1:free*j ) = &
                        trimatrix(free*(j-1)+1:free*j, free*(j-1)+1:free*j, free*(j-1)+1:free*j ) - mijk

                end do
            end do

        end associate
    end subroutine

    subroutine clac_modez( tmode, tb )
        implicit none

        ! para list
        type(tpmatrix), intent(inout) :: tmode
        real(8), intent(out) :: tb

        ! local
        integer :: i, j, k
        real(8) :: kappa, tau
        real(8) :: pkpz(tmode%mdim), ptpz(tmode%mdim), pbpz(tmode%mdim)

        associate( &
            mdim => tmode%mdim, &
            dymatrix => tmode%dymatrix0, &
            trimatrix => tmode%trimatrix, &
            modez => tmode%modez &
            )
        
            kappa = 0.d0
            do i=1, mdim
                do j=1, mdim
                    kappa = kappa + dymatrix(j,i) * modez(i) * modez(j)
                end do
            end do

            tau = 0.d0
            do i=1, mdim
                do j=1, mdim
                    do k=1, mdim
                        tau = tau + trimatrix(k,j,i) * modez(i) * modez(j) * modez(k)
                    end do
                end do
            end do

            pkpz = 0.d0
            do i=1, mdim
                do j=1, mdim
                    pkpz(i) = pkpz(i) + dymatrix(j,i) * modez(j)
                end do
            end do

            ptpz = 0.d0
            do i=1, mdim
                do j=1, mdim
                    do k=1, mdim
                        ptpz(i) = ptpz(i) + trimatrix(k,j,i) * modez(k) * modez(j)
                    end do
                end do
            end do

            pbpz = 0.5d0*kappa**2/tau**2 * pkpz - 4.d0*kappa**3/3.d0/tau**3 * ptpz

            modez = modez - pbpz*1.d-2

            !modez = modez / norm2(modez)
            modez = modez / sqrt(sum(modez**2))

            tb = 2.d0*kappa**3/3.d0/tau**2
            !print*, norm2(pbpz), tb, dot_product(pbpz,modez)/norm2(pbpz)/norm2(modez)
        end associate
    end subroutine

    subroutine solve_mode( this, oprange )
        implicit none

        class(tpmatrix) :: this
        integer, optional :: oprange

        ! local
        integer :: rangevar

        rangevar = 0
        if ( present( oprange ) ) rangevar = oprange

        associate(                        &
            mdim       => this%mdim,      &
            dymatrix   => this%dymatrix,  &
            egdymatrix => this%egdymatrix &
            )

            call solve_matrix( dymatrix, mdim, egdymatrix, rangevar )

        end associate
    end subroutine

    subroutine calc_pw( this )
        implicit none

        class(tpmatrix) :: this

        ! local
        integer :: i, j 
        real(8) :: sum4

        if ( .not. allocated(this%pw) ) allocate( this%pw(this%mdim) )

        associate(                    &
            mdim   => this%mdim,      &
            natom  => this%natom,     &
            pw     => this%pw,        &
            dymatrix => this%dymatrix &
            )
            
            do i=1, mdim
                sum4 = 0.d0
                do j=1, natom
                    sum4 = sum4 + sum( dymatrix(free*(j-1)+1:free*j,i)**2 )**2
                end do
                if ( mdim > natom * free ) then
                    sum4 = sum4 + sum( dymatrix(free*natom+1:mdim,i)**4 )
                end if
                pw(i) = 1.d0 / dble(natom) / sum4
            end do

        end associate
    end subroutine

    subroutine solve_matrix(a,order,b, rangevar)
        implicit none
        
        !-- variables input and output    
        integer :: order
        real(8),dimension(1:order,1:order) :: a
        real(8),dimension(1:order) :: b
        integer :: rangevar

        !-- local variables    
        character :: jobz = 'V'
        character :: range = 'A'
        character :: uplo = 'U'
        
        integer :: n
        !real, dimension  :: a
        integer :: lda
        integer :: vl, vu, il, iu   ! will not referenced when range=A 
        real(8) :: abstol           ! important variable
        integer :: m
        !real, dimension :: w        ! use b above ! output eigenvalues
        real(8), allocatable, dimension(:,:) :: z
        integer :: ldz
        integer, allocatable, dimension(:) :: isuppz
        real(8), allocatable, dimension(:) :: work
        integer :: lwork
        integer, allocatable, dimension(:) :: iwork
        integer :: liwork
        integer :: info
       
        n = order
        lda = order
        ldz= order

        if( rangevar == 0 ) then
            range = 'A'
        elseif( rangevar > 0 .and. rangevar <= order ) then
            range = 'I'
            il = 1
            iu = rangevar
        else
            print*, "error rangevar set"
            stop
        end if
      
        !-- 
        abstol = -1 
        !abstol = dlamch('S')        ! high precision
       
        !allocate(a(order,order)); allocate(b(order))
        allocate(z(order,order))
        allocate(isuppz(2*order))
        

        !- query the optimal workspace    
        allocate(work(1000000)); allocate(iwork(1000000))
        lwork  = -1
        liwork = -1
        
        call dsyevr(jobz,range,uplo,n,a,lda,vl,vu,il,iu,abstol,m,b,z,       &
                   & ldz, isuppz, work, lwork, iwork, liwork, info)

        lwork  = int(work(1))
        liwork = int(iwork(1))
        
        deallocate(work); deallocate(iwork)
        
        allocate(work(lwork)); allocate(iwork(liwork))
        
        !vvv- main

        call dsyevr(jobz,range,uplo,n,a,lda,vl,vu,il,iu,abstol,m,b,z,       &
                   & ldz, isuppz, work, lwork, iwork, liwork, info)

        if(info .ne. 0) then
            write(*,*) '** Find error in subroutine diago'
            stop
        end if
        !^^^- done    
        
        deallocate(work); deallocate(iwork)
        
        !- output results via a    
        a = z
        
        !deallocate(a); deallocate(b)
        deallocate(z); deallocate(isuppz)
    end subroutine

end module
