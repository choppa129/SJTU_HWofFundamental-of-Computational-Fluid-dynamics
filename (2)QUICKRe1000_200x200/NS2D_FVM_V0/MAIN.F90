
    
    
!------------------------------------------------------------------------
    program main  
!------------------------------------------------------------------------
    use global 
    use common_grid  
    use elliptic_solv 
    
    call Cart_ini  

    icontinue = 0 
    nstart = 1 
!......... MAIN LOOP .......... 
    call Velocity_bc
    call ib_object

    do nstep = nstart, 10000000   
	dt_old = dt 
	call Check_CFL 

        time = time + dt
        write(*,*) 
        write(*,'(a,I13)'  ) ' Nstep =', nstep
        write(*,'(a,F13.7)') ' time  =', time 
        write(*,'(a,F13.7)') ' Dt    =', dt
        write(*,'(a,2F13.7)') ' CFL,CFLC   =', cfl, cfl_c 
	
        call Momentum_AB2

	if(mod(nstep,200)==0) then
		call PARAVIEW(nstep)
		!pause
	endif
	if( nstep==5000000 ) then
		call save_data
		stop
	endif
    enddo
!------------------------------------------------------------------------
    end program

