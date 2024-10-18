

    
!========================================================================
    subroutine Momentum_AB2
!========================================================================
    use global
    use common_grid
    
    real*8  ::  AB_coeff, inv_J_
    real*8	Px(-2:NX+3,-2:NY+3), &
		Py(-2:NX+3,-2:NY+3)  

    write(*,*) 'step <1> : u v w mom equation' 
! ini
! center pt
    dUdx = 0.0  
    dVdy = 0.0  
    Ux   = 0.0  
    Vy   = 0.0  
! cross pt
    dVdx = 0.0 
    dUdy = 0.0
    Uy   = 0.0 
    Vx   = 0.0

! * - Gradient_UV and UV_face 
    call Gradient_UVflux_face 

! * - Reconstruct face val - 
    call Adv_UV_WENO5  !Hybrid !!Adv_UV_QUICK ! 
 
! * - Advection - Diffusion -
    FT = 0.
    GT = 0.

    call Convection_NS
    call Diffusion_NS 

! * - update U, V-
    do j=1, ny 
    do i=1, nx-1
        Uflux(i,j) = Uflux(i,j) - dt*FT(i,j)  
    enddo
    enddo

    do j=1, ny-1 
    do i=1, nx
        Vflux(i,j) = Vflux(i,j) - dt*GT(i,j) 
    enddo
    enddo

! b. c.
    call Velocity_bc
    call ib_object

! * - Projection Method -
    write(*,*) 'step <2> : projection for pressure'
    call Projection_solver  

!------------------------------------------------------------------------
    end subroutine

    
