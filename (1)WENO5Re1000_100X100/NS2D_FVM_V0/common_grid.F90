
module common_grid
!------------------------------------------------------------------------    
   
	real*8, allocatable :: &
		 P(:,:)  & 
		,Pg(:,:)  & 
		,S(:,:)  & 
		,Ux(:,:)  &
		,Uy(:,:)  &
		,Vx(:,:)  & 
		,Vy(:,:)  & 
		,Uflux(:,:)  &
		,Vflux(:,:)  & 
		,dUdx(:,:)  & ! center
		,dVdy(:,:)  & ! center
		,dVdx(:,:)  & ! cross
		,dUdy(:,:)  & ! cross
		,FT(:,:)  &
		,GT(:,:)  &  
		,HAB(:,:) &
		,UAB(:,:) &
		,VAB(:,:) & 
		,DIV(:,:)

        real*8, allocatable :: &
		 GRID_x(:,:) &
		,GRID_y(:,:) 

!------------------------------------------------------------------------
end module 
