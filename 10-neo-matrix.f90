module mo_mode
    use mo_syst
    use mo_config
    implicit none

    type tpmatrix
        real(8), allocatable, dimension(:,:)   :: dymatrix, dymatrix0
        real(8), allocatable, dimension(:)     :: egdymatrix, pw
        real(8), allocatable, dimension(:,:,:) :: trimatrix
        real(8), allocatable, dimension(:)     :: modez
        integer :: ndim, mdim, natom
        logical :: boxflag   = .false.      ! can box length change
        logical :: xyflag    = .false.      ! can box length change seperately
        logical :: shearflag = .false.      ! can box be changed by shear
    contains
!       procedure :: solve   => solve_mode
!       procedure :: calc_pw => calc_pw
    end type

!   type(tpneomatrix) :: mode, mode0, mode1, mode2

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

    subroutine kernel_matrix_fix( mdim, dymatrix, i, j, vr, vrr, xij, yij, rij )
        implicit none

        ! para list
        integer, intent(in)    :: mdim
        real(8), intent(inout) :: dymatrix(mdim,mdim)
        integer, intent(in)    :: i, j
        real(8), intent(in)    :: vr, vrr
        real(8), intent(in)    :: xij, yij, rij

        ! local
        real(8) :: rx, ry, rxx, rxy, ryy
        real(8) :: mij(free,free)

        ! \partial r_ij / \partial x_j
        rx = xij / rij
        ry = yij / rij

        ! \partial r_ij / [ \partial x_j \partial x_j ]
        rxx = 1.d0/rij - xij**2 /rij**3
        ryy = 1.d0/rij - yij**2 /rij**3
        rxy =          - xij*yij/rij**3

        ! \partial^2 V_ij / [ \partial x_i \partial x_j ]
        ! = - \partial V_ij / [ \partial x_j \partial x_j ]
        ! =   - [ Vrr * rx**2   + vr * rxx ]
        ! yy: - [ Vrr * rx**2 + Vr * ryy ]
        ! xy: - [ Vrr * rx*ry + Vr * rxy ]
        mij(1,1) = - ( vrr*rx**2 + vr*rxx )
        mij(2,2) = - ( vrr*ry**2 + vr*ryy )
        mij(1,2) = - ( vrr*rx*ry + vr*rxy )
        mij(2,1) = mij(1,2)

        dymatrix( free*(i-1)+1:free*i, free*(i-1)+1:free*i ) = &
      & dymatrix( free*(i-1)+1:free*i, free*(i-1)+1:free*i ) - mij
        dymatrix( free*(j-1)+1:free*j, free*(j-1)+1:free*j ) = &
      & dymatrix( free*(j-1)+1:free*j, free*(j-1)+1:free*j ) - mij

        dymatrix( free*(i-1)+1:free*i, free*(j-1)+1:free*j ) = mij
        dymatrix( free*(j-1)+1:free*j, free*(i-1)+1:free*i ) = mij
    end subroutine 

    subroutine kernel_matrix_compress( mdim, dymatrix, i, j, natom, vr, vrr, xij, yij, rij )
        implicit none

        ! para list
        integer, intent(in)    :: mdim
        real(8), intent(inout) :: dymatrix(mdim,mdim)
        integer, intent(in)    :: i, j
        integer, intent(in)    :: natom
        real(8), intent(in)    :: vr, vrr
        real(8), intent(in)    :: xij, yij, rij

        ! local
        real(8) :: rx, ry, rxx, rxy, ryy
        real(8) :: rex, rey
        real(8) :: rexx, rexy, reyx, reyy
        real(8) :: rexex, rexey, reyey

        real(8) :: mij(free,free)
        real(8) :: miex(free), miey(free)
        real(8) :: mee(2,2)

        ! \partial r_ij / \partial x_ij
        rx = xij / rij
        ry = yij / rij
        ! \partial r_ij / \partial \epsilon_x = rx * x_ij
        rex = rx * xij
        rey = ry * yij

        ! \partial r_ij / [ \partial x_j \partial x_j ]
        rxx = 1.d0/rij - xij**2 /rij**3
        ryy = 1.d0/rij - yij**2 /rij**3
        rxy =          - xij*yij/rij**3
        ! \partial^2 r_ij / [ \partial x_j \partial x_j ]
        rxx = 1.d0/rij - xij**2 /rij**3
        ryy = 1.d0/rij - yij**2 /rij**3
        rxy =          - xij*yij/rij**3
        ! \partial^2 r_ij / [ \partial x_j \partial \epsilon_x ]
        rexx = rxx * xij
        rexy = rxy * xij
        reyx = rxy * yij
        reyy = ryy * yij
        !
        rexex = rxx * xij**2 
        reyey = ryy * yij**2
        rexey = ryy * xij*yij

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

        ! con ij
        dymatrix( free*(i-1)+1:free*i, free*(i-1)+1:free*i ) = &
      & dymatrix( free*(i-1)+1:free*i, free*(i-1)+1:free*i ) - mij
        dymatrix( free*(j-1)+1:free*j, free*(j-1)+1:free*j ) = &
      & dymatrix( free*(j-1)+1:free*j, free*(j-1)+1:free*j ) - mij

        dymatrix( free*(i-1)+1:free*i, free*(j-1)+1:free*j ) = mij
        dymatrix( free*(j-1)+1:free*j, free*(i-1)+1:free*i ) = mij
        
        ! ex ey
        dymatrix( free*(i-1)+1:free*i, free*natom+1 ) = &
      & dymatrix( free*(i-1)+1:free*i, free*natom+1 ) + miex
        dymatrix( free*natom+1, free*(i-1)+1:free*i ) = &
      & dymatrix( free*natom+1, free*(i-1)+1:free*i ) + miex
        dymatrix( free*(i-1)+1:free*i, free*natom+2 ) = &
      & dymatrix( free*(i-1)+1:free*i, free*natom+2 ) + miey
        dymatrix( free*natom+2, free*(i-1)+1:free*i ) = &
      & dymatrix( free*natom+2, free*(i-1)+1:free*i ) + miey

        dymatrix( free*(j-1)+1:free*j, free*natom+1 ) = &
      & dymatrix( free*(j-1)+1:free*j, free*natom+1 ) - miex
        dymatrix( free*natom+1, free*(j-1)+1:free*j ) = &
      & dymatrix( free*natom+1, free*(j-1)+1:free*j ) - miex
        dymatrix( free*(j-1)+1:free*j, free*natom+2 ) = &
      & dymatrix( free*(j-1)+1:free*j, free*natom+2 ) - miey
        dymatrix( free*natom+2, free*(j-1)+1:free*j ) = &
      & dymatrix( free*natom+2, free*(j-1)+1:free*j ) - miey

        dymatrix(free*natom+1:free*natom+2,free*natom+1:free*natom+2) = &
      & dymatrix(free*natom+1:free*natom+2,free*natom+1:free*natom+2) + mee
        
    end subroutine

    subroutine kernel_matrix_compress_shear( mdim, dymatrix, i, j, natom, vr, vrr, xij, yij, rij )
        implicit none

        ! para list
        integer, intent(in)    :: mdim
        real(8), intent(inout) :: dymatrix(mdim,mdim)
        integer, intent(in)    :: i, j
        integer, intent(in)    :: natom
        real(8), intent(in)    :: vr, vrr
        real(8), intent(in)    :: xij, yij, rij

        ! local
        real(8) :: rx, ry, rxx, rxy, ryy
        real(8) :: rex, rey, res
        real(8) :: rexx, rexy, reyx, reyy, resx, resy
        real(8) :: rexex, rexey, rexes, reyey, reyes, reses

        real(8) :: mij(free,free)
        real(8) :: miex(free), miey(free), mies(free)
        real(8) :: mee(3,3)

        ! \partial r_ij / \partial x_ij
        rx = xij / rij
        ry = yij / rij
        ! \partial r_ij / \partial \epsilon_x = rx * x_ij
        rex = rx * xij
        rey = ry * yij
        res = rx * yij

        ! \partial r_ij / [ \partial x_j \partial x_j ]
        rxx = 1.d0/rij - xij**2 /rij**3
        ryy = 1.d0/rij - yij**2 /rij**3
        rxy =          - xij*yij/rij**3
        ! \partial^2 r_ij / [ \partial x_j \partial x_j ]
        rxx = 1.d0/rij - xij**2 /rij**3
        ryy = 1.d0/rij - yij**2 /rij**3
        rxy =          - xij*yij/rij**3
        ! \partial^2 r_ij / [ \partial x_j \partial \epsilon_x ]
        rexx = rxx * xij
        rexy = rxy * xij
        reyx = rxy * yij
        reyy = ryy * yij
        resx = rxx * yij
        resy = rxy * yij
        !
        rexex = rxx * xij**2 
        rexey = ryy * xij*yij
        rexes = rxx * xij*yij
        reyey = ryy * yij**2
        reyes = rxy * yij**2
        reses = rxx * yij**2

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
        mies(1) = - ( vrr*rx*res + vr*resx )
        mies(2) = - ( vrr*ry*res + vr*resy )

        mee(1,1) = vrr*rex*rex + vr*rexex
        mee(2,2) = vrr*rey*rey + vr*reyey
        mee(3,3) = vrr*res*res + vr*reses
        mee(1,2) = vrr*rex*rey + vr*rexey
        mee(2,1) = mee(1,2)

        ! con ij
        dymatrix( free*(i-1)+1:free*i, free*(i-1)+1:free*i ) = &
      & dymatrix( free*(i-1)+1:free*i, free*(i-1)+1:free*i ) - mij
        dymatrix( free*(j-1)+1:free*j, free*(j-1)+1:free*j ) = &
      & dymatrix( free*(j-1)+1:free*j, free*(j-1)+1:free*j ) - mij

        dymatrix( free*(i-1)+1:free*i, free*(j-1)+1:free*j ) = mij
        dymatrix( free*(j-1)+1:free*j, free*(i-1)+1:free*i ) = mij
        
        ! ex ey
        dymatrix( free*(i-1)+1:free*i, free*natom+1 ) = &
      & dymatrix( free*(i-1)+1:free*i, free*natom+1 ) + miex
        dymatrix( free*natom+1, free*(i-1)+1:free*i ) = &
      & dymatrix( free*natom+1, free*(i-1)+1:free*i ) + miex
        dymatrix( free*(i-1)+1:free*i, free*natom+2 ) = &
      & dymatrix( free*(i-1)+1:free*i, free*natom+2 ) + miey
        dymatrix( free*natom+2, free*(i-1)+1:free*i ) = &
      & dymatrix( free*natom+2, free*(i-1)+1:free*i ) + miey

        dymatrix( free*(j-1)+1:free*j, free*natom+1 ) = &
      & dymatrix( free*(j-1)+1:free*j, free*natom+1 ) - miex
        dymatrix( free*natom+1, free*(j-1)+1:free*j ) = &
      & dymatrix( free*natom+1, free*(j-1)+1:free*j ) - miex
        dymatrix( free*(j-1)+1:free*j, free*natom+2 ) = &
      & dymatrix( free*(j-1)+1:free*j, free*natom+2 ) - miey
        dymatrix( free*natom+2, free*(j-1)+1:free*j ) = &
      & dymatrix( free*natom+2, free*(j-1)+1:free*j ) - miey

        dymatrix(free*natom+1:free*natom+2,free*natom+1:free*natom+2) = &
      & dymatrix(free*natom+1:free*natom+2,free*natom+1:free*natom+2) + mee
        
    end subroutine
end module
end module
