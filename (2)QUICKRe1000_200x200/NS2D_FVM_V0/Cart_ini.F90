
!========================================================================
    subroutine  Cart_ini
!========================================================================
    use global
    use common_grid 
! constant 
	tttw     = (13.0/12.0)
	fourth   = (1.0/4.0)
	third    = (1.0/3.0)
	sevsix   = (7.0/6.0)
	elvsix   = (11.0/6.0)
	sixth    = (1.0/6.0)
	fivsix   = (5.0/6.0)
	tenth    = (1.0/10.0)
	sixten   = (6.0/10.0)
	treten   = (3.0/10.0)
	smallnum = (1.0D-20)
	PI = 3.141592657589793d0
! parameter
	nx = 200
	ny = 200
	
	time = 0.0 
	delta = 1.D0/real(ny-1)

	cfl = 0.2  
	dt = 0.001
	dt_old = dt  
	dt_max = dt

	visc1 = 0.001
	den1 =  1. 
	Gravity = 0.
 
        dx = delta
        dy = delta 
! Type,1:Dirichlet; 2:0-Neumman 
!c|----------------| 
!c|      -4-       |
!c|                |
!c|-1-          -3-|
!c|                |
!c|      -2-       |
!c|----------------| 
	U_type = (/ 1, 1, 1, 1 /)
	V_type = (/ 1, 1, 1, 1 /)
	P_type = (/ 0, 0, 0, 0 /)

	U_valu = (/0.,0.,0., 1./)
	V_valu = (/0.,0.,0., 0./)
	P_valu = (/0.,0.,0., 0./)

! immersed object 
!	cube_xmin = 0.4
!	cube_xmax = 0.6
!	cube_ymin =-0.1 
!	cube_ymax = 0.1

! allocate space    
	include 'Cart_alloc.F90'

	do i = -2, nx+3 
	    grid_x(i,:) = (1.0D0*real(i)-1.D0)*dx 
	enddo
	do j = -2, ny+3
	    grid_y(:,j) = (1.0D0*real(j)-1.D0)*dy-0.5D0

	enddo 

! ini vel
	do j = -2, ny+3 
	do i = -2, nx+3 
		Uflux(i,j) = U_ref 
	enddo
	enddo

	do j = -2, ny+3 
	do i = -2, nx+3  
		Vflux(i,j) = 0.0  
	enddo
	enddo 
end subroutine 
!------------------------------------------------------------------------
