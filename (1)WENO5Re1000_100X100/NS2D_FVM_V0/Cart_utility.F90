


!========================================================================
    subroutine  Check_CFL
!========================================================================
    use global
    use common_grid 
    cfl_c = 0. 
    cell_u = 0.0d0
    cell_v = 0.0d0

    do j=1,ny 
    do i=1,nx 
	cell_u = MAX(cell_u, abs(Uflux(i,j)))  
	cell_v = MAX(cell_v, abs(Vflux(i,j)))  
	cfl_c = MAX(cell_u*dt/dx, cell_v*dt/dy)
    enddo 
    enddo
	write(*,* ) ' cell_u =', cell_u,' cell_v =', cell_v

    if( cell_u==0. .and. cell_v==0. ) then
	dt = 0.001
    else
	dt = min(cfl*dx/cell_u, cfl*dy/cell_v)
    endif

    dt = dt*0.7 + dt_old*0.3
    dt = min(dt_max, dt)

    end subroutine
!------------------------------------------------------------------------
    
        


!------------------------------------------------------------------------
    real*8 function sign_(a)
!------------------------------------------------------------------------
	if(a>0) then
		sign_ = 1
	elseif(a<0) then
		sign_ = -1
	else
		sign_ = 1 !0
	endif
    end function
!------------------------------------------------------------------------
    


    
!========================================================================
    subroutine  PARAVIEW(nt)
!========================================================================
    use global
    use common_grid
    
    character*100 tec, tim, ch
    data tim       /'0000000'/
    data tec       /'FXXXXXXX.vtk'/
    
    write(tim,'(I7)') nt
    call modify_array(tim,7)
    tec(2:8)   = tim
    
    open (unit=1, file=tec)
    !paraview vtk output
    write(1,'(a)') '# vtk DataFile Version 3.0'
    write(1,'(a)') 'Non-uniform Rectilinear - Rectilinear Grid' 
    write(1,'(a)') 'ASCII' 
    write(1,'(a)') 'DATASET RECTILINEAR_GRID' 
    write(1,'(a,3i10)') 'DIMENSIONS',nx+2,ny+2,1
    write(1,'(a,i10,a)') 'X_COORDINATES',nx+2,'  float'
    write(1,'(f18.9)') (grid_x(i,1),i=1-1,nx+1)
    write(1,'(a,i10,a)') 'Y_COORDINATES',ny+2,'  float'
    write(1,'(f18.9)') (grid_y(1,j),j=1-1,ny+1)
    write(1,'(a,i10,a)') 'Z_COORDINATES',1,'  float'
    write(1,'(f18.9)') (0.0,i=1,1)
    write(1,'(a,i10)') 'POINT_DATA',(nx+2)*(ny+2)*1 
    write(1,'(a)') 'SCALARS P float 1'
    write(1,'(a)') 'LOOKUP_TABLE default'
    do j = 1-1,ny+1
    do i = 1-1,nx+1 
        write(1,'(F18.9)') P(i,j) 
    enddo
    enddo
    write(1,'(a)') 'VECTORS VELOCITY float'
    do j = 1-1,ny+1
    do i = 1-1,nx+1 
    write(1,'(3F18.9)') (Uflux(i,j)+Uflux(i-1,j))*0.5, (Vflux(i,j)+Vflux(i,j-1))*0.5, 0.0
    enddo
    enddo
    close(1)
    end subroutine
!------------------------------------------------------------------------
    
    
    
    
    
    subroutine modify_array(array,na)
    character  ::  array(na)
    integer    ::  i
    
    do i = 1, na
        if( array(i)==' ' ) then
        array(i)='0'
        endif
    enddo
    end subroutine
    
    
     

!------------------------------------------------------------------------
