


!========================================================================
    subroutine Convection_NS
!========================================================================
    use global
    use common_grid 
!                     Uy,Vy (i,j)
!c                |----------------| 
!c                |                |
!c                |                |
!c  Ux,Vx (i-1,j) |      (i,j)     | Ux,Vx (i,j)
!c                |                |
!c                |                |
!c                |----------------|
!                    Uy,Vy (i,j-1)
    write(*,*) 'Convection_NS...' 
! -U-
    do j=1, ny 
    do i=1, nx  
	Uadv1 = (Uflux(i-1,j)+Uflux(i,j))*0.5
	Uadv2 = (Uflux(i,j)  +Uflux(i+1,j))*0.5
	Vadv1 = (Vflux(i,j-1)+Vflux(i+1,j-1))*0.5  
	Vadv2 = (Vflux(i,j)  +Vflux(i+1,j))*0.5  
	U_W   = Ux(i,j)
	U_E   = Ux(i+1,j)
	U_S   = Uy(i,j-1)
	U_N   = Uy(i,j)

	SW = 0.
	SE = 0.
	SS = 0.
	SN = 0.
	SP = 0.
	! normal cell 
	SW = delta
	SE = delta
	SS = delta
	SN = delta
	SP = 0.
	Vol = delta*delta
	FT(i,j) =(Uadv2*U_E*SE - Uadv1*U_W*SW + Vadv2*U_N*SN - Vadv1*U_S*SS )  /Vol 
    enddo
    enddo

! -V-
    do j=1, ny 
    do i=1, nx  
	Uadv1 = (Uflux(i-1,j)+Uflux(i-1,j+1) )*0.5
	Uadv2 = (Uflux(i,j)  +Uflux(i,j+1) )*0.5
	Vadv1 = (Vflux(i,j)  +Vflux(i,j-1))*0.5  
	Vadv2 = (Vflux(i,j+1)+Vflux(i,j))*0.5
	V_W   = Vx(i-1,j)
	V_E   = Vx(i,j)
	V_S   = Vy(i,j)
	V_N   = Vy(i,j+1) 

	SW = 0.
	SE = 0.
	SS = 0.
	SN = 0.
	SP = 0.
	! normal cell  
	SW = delta
	SE = delta
	SS = delta
	SN = delta
	SP = 0.
	Vol = delta*delta 
	GT(i,j) =(Uadv2*V_E*SE - Uadv1*V_W*SW + Vadv2*V_N*SN - Vadv1*V_S*SS ) /Vol  
    enddo
    enddo
!stop
    end subroutine
!------------------------------------------------------------------------
    



!========================================================================
    subroutine Diffusion_NS
!========================================================================
    use global
    use common_grid

!c|----------------| 
!c|4     -4-      3|
!c|                |
!c|-1-          -3-|
!c|                |
!c|1     -2-      2|
!c|----------------| 

    write(*,*) 'Diffusion_NS...' 
! -U-
    do j=1, ny 
    do i=1, nx
	dUdx_W   = dUdx(i,j) 
	dUdx_E   = dUdx(i+1,j)
	dUdy_S   = dUdy(i,j-1) 
	dUdy_N   = dUdy(i,j) 

	SW = 0.
	SE = 0.
	SS = 0.
	SN = 0.
	SP = 0.
	! normal cell 	
	SW = delta
	SE = delta
	SS = delta
	SN = delta
	SP = 0.
	Vol = delta*delta

	den_ = den1
	visc_W = visc1
	visc_E = visc1
	visc_S = visc1
	visc_N = visc1

	FT(i,j) = FT(i,j) - (dUdx_E*SE*visc_E - dUdx_W*SW*visc_W + &
			     dUdy_N*SN*visc_N - dUdy_S*SS*visc_S ) /den_/Vol 
    enddo
    enddo

! -V-
    do j=1, ny 
    do i=1, nx
	dVdx_W   = dVdx(i-1,j) 
	dVdx_E   = dVdx(i,j) 
	dVdy_S   = dVdy(i,j) 
	dVdy_N   = dVdy(i,j+1) 

	SW = 0.
	SE = 0.
	SS = 0.
	SN = 0.
	SP = 0.
	! normal cell 	
	SW = delta
	SE = delta
	SS = delta
	SN = delta
	SP = 0.
	Vol = delta*delta

	den_ = den1 
	visc_W = visc1
	visc_E = visc1
	visc_S = visc1
	visc_N = visc1

	GT(i,j) = GT(i,j) - (dVdx_E*SE*visc_E - dVdx_W*SW*visc_W + &
			     dVdy_N*SN*visc_N - dVdy_S*SS*visc_S ) /den_/Vol  
    enddo
    enddo

    end subroutine
!------------------------------------------------------------------------




!========================================================================
    subroutine Convection_P(Px, Py)
!========================================================================
    use global
    use common_grid
    real*8	Px(-2:NX+3,-2:NY+3), &
		Py(-2:NX+3,-2:NY+3) 
! IND(i,j) = 1: inerface cell
! IND(i,j) =-2: solid cell
! IND(i,j) = 0: normal cell
! -U-
    do j=1, ny 
    do i=1, nx  
	Uadv1 = Uflux(i-1,j) 
	Uadv2 = Uflux(i,j)  
	Vadv1 = Vflux(i,j-1) 
	Vadv2 = Vflux(i,j) 
	P_W   = Px(i,j)
	P_E   = Px(i+1,j)
	P_S   = Py(i,j)
	P_N   = Py(i,j+1)

	SW = 0.
	SE = 0.
	SS = 0.
	SN = 0.
	SP = 0.
	! normal cell 	 
	SW = delta
	SE = delta
	SS = delta
	SN = delta
	SP = 0.
	Vol = delta*delta 
	HT_ij = (Uadv2*P_E*SE - Uadv1*P_W*SW + Vadv2*P_N*SN - Vadv1*P_S*SS ) /Vol 
	KT_ij = (Uadv2*SE - Uadv1*SW + Vadv2*SN - Vadv1*SS ) /Vol*P(i,j)
	Pg(i,j) = P(i,j) - dt * (HT_ij-KT_ij)  

    enddo
    enddo

    end subroutine
!------------------------------------------------------------------------





!========================================================================
    subroutine Gradient_UVflux_face
!========================================================================
    use global
    use common_grid 
    write(*,*) 'Gradient_UVflux_face...'

    do j=1-2, ny+2 
    do i=1-2, nx+2
	dUdx(i,j) = (Uflux(i,j) - Uflux(i-1,j))/delta
	dUdy(i,j) = (Uflux(i,j+1) - Uflux(i,j))/delta 
    enddo
    enddo

    do j=1-2, ny+2 
    do i=1-2, nx+2
	dVdx(i,j) = (Vflux(i+1,j) - Vflux(i,j))/delta
	dVdy(i,j) = (Vflux(i,j) - Vflux(i,j-1))/delta
    enddo
    enddo 

    end subroutine
!------------------------------------------------------------------------
    
