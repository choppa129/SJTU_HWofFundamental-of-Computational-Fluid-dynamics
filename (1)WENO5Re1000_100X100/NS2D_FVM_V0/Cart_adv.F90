



!========================================================================
    subroutine Adv_UV_QUICK
!========================================================================
    use global
    use common_grid

    write(*,*) 'Adv_UV...'

    do j=1, ny  
    do i=1, nx  
	! X-face for U
	Uf_ = (Uflux(i-1,j) + Uflux(i,j))*0.5
	if(Uf_>0) then
		q2 = Uflux(i-2, j) 
		q3 = Uflux(i-1, j) 
		q4 = Uflux(i  , j)  
	elseif(Uf_<0) then
		q2 = Uflux(i+1, j) 
		q3 = Uflux(i  , j) 
		q4 = Uflux(i-1, j) 
	else
		Ux(i,j) = Uf_; cycle
	endif
	Ux(i,j) = -1./8.*q2 + 3./4.*q3 + 3./8.*q4
    enddo
    enddo

    call Scalar_bc(Ux)

    do j=1, ny  
    do i=1, nx 
	! X-face for V 
	Uf_ = Uflux(i,j) + Uflux(i,j+1)
	if(Uf_>0) then
		q2 = Vflux(i-1, j) 
		q3 = Vflux(i, j) 
		q4 = Vflux(i+1, j)  
	elseif(Uf_<0) then 
		q2 = Vflux(i+2, j)  
		q3 = Vflux(i+1, j) 
		q4 = Vflux(i, j) 
	else
		Vx(i,j) = (Vflux(i, j) + Vflux(i+1, j))*0.5; cycle
	endif
	Vx(i,j) = -1./8.*q2 + 3./4.*q3 + 3./8.*q4  
    enddo
    enddo


    call Scalar_bc(Vx)

    do j=1, ny
    do i=1, nx   
	! Y-face for U 
	Vf_ = Vflux(i,j)+Vflux(i+1,j)
	if(Vf_>0) then
		q2 = Uflux(i, j-1) 
		q3 = Uflux(i, j) 
		q4 = Uflux(i, j+1)  
	elseif(Vf_<0) then
		q2 = Uflux(i, j+2)  
		q3 = Uflux(i, j+1) 
		q4 = Uflux(i, j) 
	else
		Uy(i,j) = (Uflux(i, j+1)+Uflux(i, j))*0.5; cycle
	endif
	Uy(i,j) = -1./8.*q2 + 3./4.*q3 + 3./8.*q4 
    enddo
    enddo

    call Scalar_bc(Uy)

    do j=1, ny
    do i=1, nx   
	! Y-face for V 
	Vf_ = (Vflux(i, j-1) + Vflux(i, j))*0.5 
	if(Vf_>0) then
		q2 = Vflux(i, j-2) 
		q3 = Vflux(i, j-1) 
		q4 = Vflux(i, j)  
	elseif(Vf_<0) then
		q2 = Vflux(i, j+1) 
		q3 = Vflux(i, j) 
		q4 = Vflux(i, j-1) 
	else
		Vy(i,j) = Vf_; cycle
	endif
	Vy(i,j) = -1./8.*q2 + 3./4.*q3 + 3./8.*q4 
    enddo
    enddo

    call Scalar_bc(Vy)

    end subroutine
!------------------------------------------------------------------------




!========================================================================
    subroutine Adv_UV_WENO5
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

    write(*,*) 'Adv_CellFlux...' 
    do j=1, ny  
    do i=1, nx  
	! X-face for U 
	Uf_ = (Uflux(i-1,j) + Uflux(i,j))*0.5
	if(Uf_==0.) then 
		Ux(i,j) = Uf_
		cycle 
	elseif(Uf_>0.) then 
		q1 = Uflux(i-3, j) 
		q2 = Uflux(i-2, j) 
		q3 = Uflux(i-1, j) 
		q4 = Uflux(i  , j) 
		q5 = Uflux(i+1, j) 
	elseif(Uf_<0.) then 
		q1 = Uflux(i+2, j) 
		q2 = Uflux(i+1, j) 
		q3 = Uflux(i  , j) 
		q4 = Uflux(i-1, j) 
		q5 = Uflux(i-2, j)
	endif
	call alpha_weight(q1,q2,q3,q4,q5,w1,w2,w3) 
	Ux(i,j) =  (w1*( q1*third - q2*sevsix + q3*elvsix) &
		+ w2*(-q2*sixth + q3*fivsix + q4*third)  &
		+ w3*( q3*third + q4*fivsix - q5*sixth)) 
    enddo
    enddo

    call Scalar_bc(Ux)

    do j=1, ny  
    do i=1, nx 
	! X-face for V 
	Uf_ = Uflux(i,j) + Uflux(i,j+1)
	if(Uf_==0.0) then
		Vx(i,j) = (Vflux(i, j) + Vflux(i+1, j))*0.5
		cycle
	elseif(Uf_>0.) then
		q1 = Vflux(i-2, j) 
		q2 = Vflux(i-1, j) 
		q3 = Vflux(i  , j) 
		q4 = Vflux(i+1, j) 
		q5 = Vflux(i+2, j) 
	elseif(Uf_<0.) then
		q1 = Vflux(i+3, j) 
		q2 = Vflux(i+2, j) 
		q3 = Vflux(i+1, j) 
		q4 = Vflux(i  , j) 
		q5 = Vflux(i-1, j) 
	endif 
	call alpha_weight(q1,q2,q3,q4,q5,w1,w2,w3) 
	Vx(i,j) =  (w1*( q1*third - q2*sevsix + q3*elvsix) &
		+ w2*(-q2*sixth + q3*fivsix + q4*third)  &
		+ w3*( q3*third + q4*fivsix - q5*sixth)) 
    enddo
    enddo 

    call Scalar_bc(Vx)

    do j=1, ny
    do i=1, nx   
	! Y-face for U 
	Vf_ = Vflux(i,j)+Vflux(i+1,j)
	if(Vf_==0.0) then
		Uy(i,j) = (Uflux(i, j+1)+Uflux(i, j))*0.5
		cycle
	elseif(Vf_>0.) then 
		q1 = Uflux(i, j-2) 
		q2 = Uflux(i, j-1) 
		q3 = Uflux(i, j  ) 
		q4 = Uflux(i, j+1) 
		q5 = Uflux(i, j+2) 
	elseif(Vf_<0.) then 
		q1 = Uflux(i, j+3) 
		q2 = Uflux(i, j+2) 
		q3 = Uflux(i, j+1) 
		q4 = Uflux(i, j  ) 
		q5 = Uflux(i, j-1)
	endif 
	call alpha_weight(q1,q2,q3,q4,q5,w1,w2,w3) 
	Uy(i,j) =  (w1*( q1*third - q2*sevsix + q3*elvsix) &
		+ w2*(-q2*sixth + q3*fivsix + q4*third)  &
		+ w3*( q3*third + q4*fivsix - q5*sixth)) 
    enddo
    enddo

    call Scalar_bc(Uy)

    do j=1, ny
    do i=1, nx   
	! Y-face for V 
	Vf_ = (Vflux(i, j-1) + Vflux(i, j))*0.5 
	if(Vf_==0.0) then
		Vy(i,j) = Vf_
		cycle 
	elseif(Vf_>0.) then
		q1 = Vflux(i, j-3) 
		q2 = Vflux(i, j-2) 
		q3 = Vflux(i, j-1) 
		q4 = Vflux(i, j  ) 
		q5 = Vflux(i, j+1) 
	elseif(Vf_<0.) then 
		q1 = Vflux(i, j+2) 
		q2 = Vflux(i, j+1) 
		q3 = Vflux(i, j  ) 
		q4 = Vflux(i, j-1) 
		q5 = Vflux(i, j-2)
	endif
	call alpha_weight(q1,q2,q3,q4,q5,w1,w2,w3) 
	Vy(i,j) =  (w1*( q1*third - q2*sevsix + q3*elvsix) &
		+ w2*(-q2*sixth + q3*fivsix + q4*third)  &
		+ w3*( q3*third + q4*fivsix - q5*sixth)) 
    enddo
    enddo

    call Scalar_bc(Vy)

    end subroutine
!------------------------------------------------------------------------




!========================================================================
    subroutine Adv_P_WENO5(Px, Py)
!========================================================================
    use global
    use common_grid

    real*8	Px(-2:NX+3,-2:NY+3), &
		Py(-2:NX+3,-2:NY+3)  

    write(*,*) 'Adv_P...'

    do j=1, ny  
    do i=1, nx  
	! X-face for P  
	Uf_ = (Uflux(i-1,j) + Uflux(i,j))*0.5 
	if(Uf_==0.) then 
		Px(i,j) = 0.0
		cycle 
	elseif(Uf_>0.) then
		q1 = P(i-3, j)  !Uflux(i-3, j) 
		q2 = P(i-2, j)  !Uflux(i-2, j) 
		q3 = P(i-1, j)  !Uflux(i-1, j) 
		q4 = P(i  , j)  !Uflux(i  , j) 
		q5 = P(i+1, j)  !Uflux(i+1, j) 
	elseif(Uf_<0.) then
		q1 = P(i+2, j)  !Uflux(i+2, j) 
		q2 = P(i+1, j)  !Uflux(i+1, j) 
		q3 = P(i  , j)  !Uflux(i  , j) 
		q4 = P(i-1, j)  !Uflux(i-1, j) 
		q5 = P(i-2, j)  !Uflux(i-2, j)
	endif
	call alpha_weight(q1,q2,q3,q4,q5,w1,w2,w3) 
	Px(i,j) =(w1*( q1*third - q2*sevsix + q3*elvsix) &
		+ w2*(-q2*sixth + q3*fivsix + q4*third)  &
		+ w3*( q3*third + q4*fivsix - q5*sixth)) 
    enddo
    enddo

    call Scalar_bc(Px)

    do j=1, ny
    do i=1, nx   
	! Y-face for P 
	Vf_ = (Vflux(i, j-1) + Vflux(i, j))*0.5 
	if(Vf_==0.0) then
		Py(i,j) = 0.0
		cycle 
	elseif(Vf_>0.) then
		q1 = P(i, j-3) !Vflux(i, j-3) 
		q2 = P(i, j-2) !Vflux(i, j-2) 
		q3 = P(i, j-1) !Vflux(i, j-1) 
		q4 = P(i, j  ) !Vflux(i, j  ) 
		q5 = P(i, j+1) !Vflux(i, j+1) 
	elseif(Vf_<0.) then 
		q1 = P(i, j+2) !Vflux(i, j+2) 
		q2 = P(i, j+1) !Vflux(i, j+1) 
		q3 = P(i, j  ) !Vflux(i, j  ) 
		q4 = P(i, j-1) !Vflux(i, j-1) 
		q5 = P(i, j-2) !Vflux(i, j-2)
	endif
	call alpha_weight(q1,q2,q3,q4,q5,w1,w2,w3) 
	Py(i,j) =  (w1*( q1*third - q2*sevsix + q3*elvsix) &
		+ w2*(-q2*sixth + q3*fivsix + q4*third)  &
		+ w3*( q3*third + q4*fivsix - q5*sixth)) 
    enddo
    enddo

    call Scalar_bc(Py)

    end subroutine
!------------------------------------------------------------------------





!========================================================================
    subroutine alpha_weight(q1,q2,q3,q4,q5,w1,w2,w3)! 
!========================================================================
    use global
    
    epsilon  =  0.000001D0
    
	si1 = tttw*(q1-2.0*q2+q3)**( 2.0) + fourth*(q1-4.0*q2+3.0*q3)**( 2.0) 
	si2 = tttw*(q2-2.0*q3+q4)**( 2.0) + fourth*(q2-q4)**( 2.0) 
	si3 = tttw*(q3-2.0*q4+q5)**( 2.0) + fourth*(3.0*q3-4.0*q4+q5)**( 2.0) 
	
	alpha1 = tenth*(epsilon+si1)**(-2.0)
	alpha2 = sixten*(epsilon+si2)**(-2.0)
	alpha3 = treten*(epsilon+si3)**(-2.0)
	
	w1 = alpha1/(alpha1+alpha2+alpha3)
	w2 = alpha2/(alpha1+alpha2+alpha3)
	w3 = alpha3/(alpha1+alpha2+alpha3)
!------------------------------------------------------------------------
    endsubroutine
    
