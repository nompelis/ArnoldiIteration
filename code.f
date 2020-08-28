       Program main
       Implicit none
       Real*8,allocatable,dimension(:,:) :: amat
       Real*8,allocatable,dimension(:) :: bvec
       Real*8,allocatable,dimension(:,:) :: q,h
       Integer*4 :: size, num, i,k

       size = 10
       num = 4

       allocate( amat(size,size), bvec(size), q(size,num), h(num,num) )
       amat(:,:) = 0.0d0
       bvec(:) = 0.0d0

       call setup( size, amat, bvec )

       call arnoldi_basis( size, amat, bvec, num, q, h )

       PRINT*,'Vectors'
       do k = 1,num
          PRINT*,'k=',k
          do i = 1,size
             PRINT*,'i=',i,'   q=',q(i,k)
          enddo
       enddo

       call check_orthogonality ( size, num, q )

       End

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       Subroutine arnoldi_basis( size, amat, bvec, num, qvecs, hmat )
       Implicit none
       Integer*4, intent(IN) :: size, num
       Real*8,dimension( size, size ), intent(IN) :: amat
       Real*8,dimension( size ), intent(IN) :: bvec
       Real*8,dimension( size, num ), intent(OUT) :: qvecs
       Real*8,dimension( num, num ), intent(OUT) :: hmat
       Integer*4 :: i,j,k,kk
       Real*8 :: r

       !--- form the first vector assuming it is not normalized
       r = 0.0d0
       do i = 1,size
          r = r + bvec(i)*bvec(i)
       enddo
       r = sqrt(r)
       !--- magnitude of the initial vector; assume it is not a unit vector
       hmat(1,1) = r
       !--- store the first vector normalized
       r = 1.0d0/r
       do i = 1,size
          qvecs(i,1) = bvec(i) * r
       enddo

       !--- form the subsequent vectors (so that there are "num" in total)
       do k = 2,num
          !--- map to a new vector using the matrix
          do i = 1,size
             r = 0.0d0
             do j = 1,size
                r = r + amat(i,j) * qvecs(j,k-1)
             enddo
             qvecs(i,k) = r
          enddo

          !--- make orthogonal to all previous vectors (up to k-1)
          do kk = 1,k-1
             r = 0.0d0
             do i = 1,size
                r = r + qvecs(i,k) * qvecs(i,kk)
             enddo
             !--- store projetion to the Hezenberg matrix
             hmat(k,kk) = r
             !--- subtract the projection vectorwise
             do i = 1,size
                qvecs(i,k) = qvecs(i,k) - r * qvecs(i,kk)
             enddo
          enddo
          !--- calculate the magnitude
          r = 0.0d0
          do i = 1,size
             r = r + qvecs(i,k) * qvecs(i,k)
          enddo
          r = sqrt(r)
          hmat(k,k) = r
          r = 1.0d0/r
          do i = 1,size
             qvecs(i,k) = qvecs(i,k) * r
          enddo
       enddo

       End

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       Subroutine setup( size, amat, bvec )
       Implicit none
       Integer*4, intent(IN) :: size
       Real*8,dimension( size, size ), intent(INOUT) :: amat
       Real*8,dimension( size ), intent(INOUT) :: bvec
       Integer*4 :: i,j

       do i = 1,size
       do j = 1,size
          if( i.eq.j ) amat(i,j) = -2.0d0
          if( i.eq.j-1) amat(i,j) = 1.0d0
          if( i.eq.j+1) amat(i,j) = 1.0d0
       enddo
       enddo
       amat(1,:) = 0.0d0
       amat(1,1) = 1.0d0
       amat(size,:) = 0.0d0
       amat(size,size) = 1.0d0

       PRINT*,'Matrix'
       do i = 1,size
          PRINT*,'Row',i
       do j = 1,size
          PRINT*,amat(i,j)
       enddo
       enddo

       do i = 1,size
          bvec(i) = 1.0d0
       enddo

       PRINT*,'Vector'
       do i = 1,size
          PRINT*,'Row',i, bvec(i)
       enddo

       End

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       Subroutine check_orthogonality ( size, num, q )
       Implicit none
       Integer*4, intent(IN) :: size, num
       Real*8,dimension( size, num ), intent(INOUT) :: q
       Integer*4 :: m,n,i
       Real*8 :: r

       do m = 1,num
       do n = 1,num
          r = 0.0d0
          do i = 1,size
             r = r + q(i,m) * q(i,n)
          enddo
          PRINT*,'Test ortho: ',m,n,r
       enddo
       enddo

       End

