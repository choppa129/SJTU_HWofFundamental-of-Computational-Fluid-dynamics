    
    

    allocate(P (-2:NX+3,-2:NY+3)  &
            ,Pg(-2:NX+3,-2:NY+3)  & 
            ,S (-2:NX+3,-2:NY+3)  &  
            ,Ux(-2:NX+3,-2:NY+3)  &
            ,Vx(-2:NX+3,-2:NY+3)  & 
            ,Uy(-2:NX+3,-2:NY+3)  &
            ,Vy(-2:NX+3,-2:NY+3)  & 
            ,Uflux(-2:NX+3,-2:NY+3)  &
            ,Vflux(-2:NX+3,-2:NY+3)  & 
	    ,dUdx(-2:NX+3,-2:NY+3)  & 
	    ,dVdy(-2:NX+3,-2:NY+3)  & 
	    ,dVdx(-2:NX+3,-2:NY+3)  & 
	    ,dUdy(-2:NX+3,-2:NY+3)  & 
            ,FT (-2:NX+3,-2:NY+3)  &
            ,GT (-2:NX+3,-2:NY+3)  &
	    ,HAB(-2:NX+3,-2:NY+3) &
	    ,UAB(-2:NX+3,-2:NY+3) &
	    ,VAB(-2:NX+3,-2:NY+3) & 
            ,DIV(-2:NX+3,-2:NY+3)  &
            ,GRID_x(-2:NX+3,-2:NY+3)  &
            ,GRID_y(-2:NX+3,-2:NY+3) ) 

	P = 0.0  
	Pg= 0.0  
	S = 0.0  
 
	Ux = 0.0
	Vx = 0.0
	Uy = 0.0
	Vy = 0.0
	Uflux = 0.0
	Vflux = 0.0
	dUdx = 0.0
	dVdy = 0.0
	dVdx = 0.0
	dUdy = 0.0

	FT = 0.0
	GT = 0.0
	HAB = 0.0
	UAB = 0.0
	VAB = 0.0
	DIV = 0.0 

	GRID_x = 0.0
	GRID_y = 0.0 

