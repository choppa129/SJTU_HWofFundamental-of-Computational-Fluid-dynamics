


    
!========================================================================
    subroutine  save_data
!========================================================================
    use global
    use common_grid
    use elliptic_solv

    character*100    :: tec1 
    data tec1           /'Z0000000.0000X000.dat'/

    open (unit=100, file=tec1, form='unformatted')
    write(100) nstep
    write(100) time
    write(100) dt

    write(*, *) 'write to t= ', time, 'with dt= ', dt

    do j = 0, ny+1
    do i = 0, nx+1
        write(100) P (i,j), S (i,j), Ux(i,j), Vx(i,j)  & 
            ,Uy(i,j), Vy(i,j), Uflux(i,j), Vflux(i,j), dUdx(i,j)  & 
	    ,dVdy(i,j), dVdx(i,j), dUdy(i,j), FT (i,j),GT (i,j)  &
            ,DIV(i,j),GRID_x(i,j),GRID_y(i,j)  
    enddo
    enddo

    close(100)
    end subroutine
!------------------------------------------------------------------------
    
  



    
!========================================================================
    subroutine  read_data
!========================================================================
    use global
    use common_grid
    use elliptic_solv

    character*100    :: tec1 
    data tec1           /'Z0000000.0000X000.dat'/

    open (unit=100, file=tec1, form='unformatted')
    read(100) nstep
    read(100) time
    read(100) dt

    write(*, *) 'read from t= ', time, 'with dt= ', dt
    nstart = nstep

    do j = 0, ny+1
    do i = 0, nx+1
        read(100)P (i,j), S (i,j), Ux(i,j), Vx(i,j)  & 
            ,Uy(i,j), Vy(i,j), Uflux(i,j), Vflux(i,j), dUdx(i,j)  & 
	    ,dVdy(i,j), dVdx(i,j), dUdy(i,j), FT (i,j),GT (i,j)  &
            ,DIV(i,j),GRID_x(i,j),GRID_y(i,j)  
    enddo
    enddo

    close(100)
    end subroutine
!------------------------------------------------------------------------
    
  
