module elliptic_solv
!------------------------------------------------------------------------   
    IMPLICIT REAL*8(A-H,O-Z) 
    
    integer      MEM, &                ! mem size for ADJ0, COF0
                 BWH, &                ! stencil size of GG_stl, GG_cof 
                 idx_num, &            ! matrix size
                 xbc_type(2), &        ! 1: Dirichlet BC. 0: Neumman BC.
                 ybc_type(2), & 
                 
                 nibstart, nibend, &   ! start and end block NO. in each processor / no parallel
                 nunknowns, &          ! total number of unknowns
                 noffrow               ! offset of row for each processor 

    
    integer, allocatable   ::   ADJ(:), ADJrow(:)
    real*8, allocatable    ::   COF(:), rhs(:)
!------------------------------------------------------------------------
endmodule 
