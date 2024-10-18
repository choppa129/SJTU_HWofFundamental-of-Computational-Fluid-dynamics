

!========================================================================
    subroutine ib_object
!========================================================================
    use global
    use common_grid
! U
    do j=-2, ny+3 
    do i=-2, nx+3  
	xi = grid_x(i,j) + delta 
	yi = grid_y(i,j) 
!	if  ( xi>=cube_xmin .and. xi<=cube_xmax .and. &
!		yi>=cube_ymin .and. yi<=cube_ymax) then
!	Uflux(i,j) = 0.
!	endif
    enddo
    enddo
! V
    do j=-2, ny+3 
    do i=-2, nx+3  
	xi = grid_x(i,j)
	yi = grid_y(i,j) + delta  
!	if  ( xi>=cube_xmin .and. xi<=cube_xmax .and. &
!		yi>=cube_ymin .and. yi<=cube_ymax) then
!	Vflux(i,j) = 0.
!	endif
    enddo
    enddo
    end subroutine
!------------------------------------------------------------------------
    



!========================================================================
    subroutine  Scalar_bc(F)
!========================================================================
    use global     
    real*8   ::   F   (-2:NX+3,-2:NY+3) 
! 
! Bottom bc.
    do m=-2,0; F(:, m) = F(:, 1); enddo

! Top bc.
    do m=NY+1,NY+3; F(:,m) = F(:,NY); enddo 

! Left bc.
    do m=-2,0; F(m, :) = F(1, :); enddo 

! Right bc.
    do m=NX+1,NX+3; F(m,:) = F(NX,:); enddo 
    end subroutine
!------------------------------------------------------------------------
    



!========================================================================
    subroutine  Velocity_bc
!========================================================================
    use global
    use common_grid 
! Default
! Type,1:Dirichlet; 0:0-Neumman 
!c|----------------| 
!c|      -4-       |
!c|                |
!c|-1-          -3-|
!c|                |
!c|      -2-       |
!c|----------------| 

! Left bc. -x
    if    ( U_type(1) == 0) then
	  do m=-2,0; Uflux(m,:) = Uflux(1,:); enddo
    elseif( U_type(1) == 1) then
	  do m=-2,0; Uflux(m,:) = U_valu(1); enddo 
    endif 

    if    ( V_type(1) == 0) then
	  do m=-2,0; Vflux(m,:) = Vflux(1,:); enddo 
    elseif( V_type(1)== 1) then
	  do m=-2,0; Vflux(m,:) = V_valu(1); enddo  
    endif

! Right bc. +x
    if    ( U_type(3) == 0) then
	  do m=NX,NX+3; Uflux(m,:) = Uflux(NX-1,:); enddo   
    elseif( U_type(3) == 1) then
	  do m=NX,NX+3; Uflux(m,:) = U_valu(3); enddo 
    endif

    if    ( V_type(3) == 0) then
	  do m=NX+1,NX+3; Vflux(m,:) = Vflux(NX,:); enddo   
    elseif( V_type(3) == 1) then
	  do m=NX+1,NX+3; Vflux(m,:) = V_valu(3); enddo
    endif

! Bottom bc. -y
    if    ( U_type(2) == 0) then
	  do m=-2,0;Uflux(:,m) = Uflux(:,1); enddo
    elseif( U_type(2) == 1) then
	  do m=-2,0;Uflux(:,m) = U_valu(2); enddo 
    endif 

    if    ( V_type(2) == 0) then
	  do m=-2,0;Vflux(:,m) = Vflux(:,1); enddo
    elseif( V_type(2)== 1) then
	  do m=-2,0;Vflux(:,m) = V_valu(2); enddo 
    endif

! Top bc. +y
    if    ( U_type(4) == 0) then
	  do m=NY+1,NY+3;Uflux(:,m) = Uflux(:,NY); enddo 
    elseif( U_type(4) == 1) then
	  do m=NY+1,NY+3;Uflux(:,m) = U_valu(4); enddo 
    endif 

    if    ( V_type(4) == 0) then
	  do m=NY,NY+3;Vflux(:,m) = Vflux(:,NY-1); enddo 
    elseif( V_type(4) == 1) then
	  do m=NY,NY+3;Vflux(:,m) = V_valu(4); enddo 
    endif 
    end subroutine
!------------------------------------------------------------------------
     



!========================================================================
    subroutine  Press_bc
!========================================================================
    use global
    use common_grid 
    real*8 :: P_valu_func
! Default
! Type,1:Dirichlet; 0:0-Neumman 
!c|----------------| 
!c|      -4-       |
!c|                |
!c|-1-          -3-|
!c|                |
!c|      -2-       |
!c|----------------| 

! Left bc. -x
    if    ( P_type(1) == 0) then
	  do m=-2,0; P(m,:) = P(1,:); enddo
    elseif( P_type(1) == 1) then
	  do m=-2,0;do j=-2,NY+3
		xi = grid_x(m,j)
		yi = grid_y(m,j) 
		P(m,j) = P_valu_func(xi,yi) 
	  enddo;enddo 
    endif 

! Right bc. +x
    if    ( P_type(3) == 0) then
	  do m=NX+1,NX+3; P(m,:) = P(NX,:); enddo   
    elseif( P_type(3) == 1) then
	  do m=NX+1,NX+3;do j=-2,NY+3
		xi = grid_x(m,j)
		yi = grid_y(m,j) 
		P(m,j) = P_valu_func(xi,yi) 
	  enddo;enddo 
    endif

! Bottom bc. -y
    if    ( P_type(2) == 0) then
	  do m=-2,0;P(:,m) = P(:,1); enddo
    elseif( P_type(2) == 1) then
	  do m=-2,0;do i=-2,NX+3
		xi = grid_x(i,m)
		yi = grid_y(i,m) 
		P(i,m) = P_valu_func(xi,yi) 
	  enddo;enddo 
    endif 

! Top bc. +y
    if    ( P_type(4) == 0) then
	  do m=NY+1,NY+3;P(:,m) = P(:,NY); enddo 
    elseif( P_type(4) == 1) then
	  do m=NY+1,NY+3;do i=-2,NX+3
		xi = grid_x(i,m)
		yi = grid_y(i,m) 
		P(i,m) = P_valu_func(xi,yi) 
	  enddo;enddo 
    endif 
    end subroutine
!------------------------------------------------------------------------
     



real*8 function P_valu_func(xi,yi) 
    P_valu_func = 0.0
end function
