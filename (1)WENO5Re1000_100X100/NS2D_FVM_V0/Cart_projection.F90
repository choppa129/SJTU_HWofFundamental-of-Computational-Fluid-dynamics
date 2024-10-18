


!========================================================================
    subroutine Projection_solver
!========================================================================
    use global
    use common_grid
    use elliptic_solv

! * make coefficient matrix and r.h.s vector *  
    DT1 = 1.0d0 / DT 

    do j=1, ny 
    do i=1, nx 
        DIV(i,j) = DT1 * ((Uflux(I,J) - Uflux(I-1,J)) + & 
                     (Vflux(I,J) - Vflux(I,J-1)) ) * delta
    enddo
    enddo

! * Start the Elliptic solver
    call SOR_solver

! * pressure bc.
    call Scalar_bc(P) 

! * Face vel. corrector 
    do j=1,ny
    do i=1,nx-1  
	rr1 = 1.
        Uflux(i,j) = Uflux(i,j) - DT * & 
		( P(i+1,j) - P(i,j) )/delta*rr1 
    enddo
    enddo

    do j=1,ny-1
    do i=1,nx  
	rr2 = 1. 
        Vflux(i,j) = Vflux(i,j) - DT * &  
          	( P(i,j+1) - P(i,j) )/delta*rr2  
    enddo
    enddo

! * B. C.
    call Velocity_bc 
 
! debug, validate div.U = 0
    do j = 1, NY 
    do i = 1, NX  
	DIV(i,j) = DT1 * ((Uflux(I,J) - Uflux(I-1,J)) + & 
	     (Vflux(I,J) - Vflux(I,J-1)) ) * delta 
    enddo
    enddo 
    return
!---------------------------------------------------------------------
    end subroutine




!========================================================================
    subroutine JACOBI_solver
!========================================================================
    use global
    use common_grid 
    
! * SOR  
    itermin=100
    itermax=5000
    epsi=1d-3  

    do l=1,itermax
	emaxp=0. 
	do j=1,ny
	do i=1,nx
	aw1=1.
	ae1=1.
	ab1=1.
	at1=1.
	aa= -(aw1 + ae1 + ab1 + at1) 
	ff= DIV(i,j)
! jacobi solver
	P(i,j)=(ff - aw1*Pg(i-1,j) - ae1*Pg(i+1,j) - ab1*Pg(i,j-1) - at1*Pg(i,j+1))/aa 
	enddo
	enddo
! bc
	call Press_bc
! jacobi solver
	do j=-2,ny+3
	do i=-2,nx+3
		emx=abs(P(i,j)-Pg(i,j))
		if(emaxp<emx) emaxp=emx 
		Pg(i,j) = P(i,j)
	enddo
	enddo
	if (l<itermin) cycle
	if (emaxp<epsi) exit 
    enddo !l
    write(*,*) 'JACOBI iter',l,emaxp 
!---------------------------------------------------------------------
    end subroutine




!========================================================================
    subroutine SOR_solver
!========================================================================
    use global
    use common_grid 
    
! * SOR 
    alpha=1.2d0
    itermin=100
    itermax=5000
    epsi=1d-3 

    do l=1,itermax
	emaxp=0. 
	do j=1,ny
	do i=1,nx
	aw1=1. 
	ae1=1. 
	ab1=1. 
	at1=1. 
	aa= -(aw1 + ae1 + ab1 + at1) 
	ff= DIV(i,j) 
! sor
	ww=(ff - aw1*P(i-1,j) - ae1*P(i+1,j) &
		- ab1*P(i,j-1) - at1*P(i,j+1))/aa &
		- P(i,j)
	P(i,j)=P(i,j)+alpha*ww

	emx=abs(ww)
	if(emaxp < emx)emaxp=emx 
	enddo
	enddo
! bc
	call Press_bc

	if (l<itermin) cycle
	if (emaxp<epsi) exit 
    enddo !l
	
    write(*,*) 'SOR iter',l,emaxp
!---------------------------------------------------------------------
    end subroutine


