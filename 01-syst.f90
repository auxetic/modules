module mo_syst
    implicit none
        
    type tpset
        integer :: natom
        real(8) :: phi
        integer :: np
    end type
    type(tpset) :: sets

    real(8), parameter :: pi    = 3.1415926535897932d0
    real(8), parameter :: ratio = 1.4d0
    real(8), parameter :: alpha = 2.d0
    integer, parameter :: free  = 2

contains

!    subroutine testvar
!        implicit none
!            
!        sets.natom = 64
!        sets.phi = 0.91d0
!        sets.np = 201
!
!    end subroutine testvar
    
end module mo_syst

module mo_var
    implicit none
        
    integer :: i, j, k, ii, jj, kk, itemp, step, step1, step2, nstep
    real(8) :: temp1, temp2, temp3

    ! main
    real(8) :: testp
    real(8) :: Mk_x, Mk_y, mu_xy, mu_yx, MB, MG_s, MG_xy
    real(8) :: epl_xx, epl_yy
    
    logical :: exist_flag, eof_flag

    character(250) :: filename, chtemp
    character(250) :: faverage, fresults, fdump, ftemp
end module
