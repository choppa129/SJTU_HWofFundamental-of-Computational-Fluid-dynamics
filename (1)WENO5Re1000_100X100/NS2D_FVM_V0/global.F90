module global
!------------------------------------------------------------------------   
    IMPLICIT REAL*8(A-H,O-Z) 
    
    integer      IM, JM, & 
                 NX, NY, nstep,nstart
    integer      neib(4), mype, npes 

    real*8       time, dt, dt_max, dt_old, DTx, viscoeff, delta, U_ref, dx, dy, cfl, cfl_c

!       ..1: water; ..2:vapor
    real*8       Gravity, visc1, den1
! B C
    integer      U_type(4), V_type(4), P_type(4)  
    real*8       U_valu(4), V_valu(4), P_valu(4) 
!   real*8	 cube_xmin, cube_xmax, cube_ymin, cube_ymax
!   Frequently used constant variables    
    real*8       tttw,     &
                 fourth,   &
                 third,    &
                 sevsix,   & 
                 elvsix,   &
                 sixth,    &
                 fivsix,   &
                 tenth,    &
                 sixten,   &
                 treten,   & 
                 smallnum, PI
!------------------------------------------------------------------------
end module 
