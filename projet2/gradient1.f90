		module gradientconjugue
			use mpi
			use matrices
			implicit none
			contains


						subroutine gc(A,b,NxNy,Nx,Ny,xn)
							implicit none
							real , dimension(:),allocatable ::xn, p, Ap, r,b,r1, produitmatriciel



							REAL :: alpha, rr, rr1, pAp,produitvectoriel
							integer :: j,i,k,statinfo,Np,me, NxNy,Nx,Ny,i1,iN,he
							real,dimension(:,:),allocatable::A
							allocate (r(NxNy))
							allocate (r1(NxNy))

							allocate (p(NxNy))
							allocate(xn(NxNy))
	
							xn=0
							r = b
							p = b
							call produitvec(Nx,Ny,r,r,produitvectoriel,statinfo,me,i1,iN)
							rr = produitvectoriel
							DO k = 1, NxNy
								call produitvec(Nx,Ny,r,r,produitvectoriel,statinfo,me,i1,iN)
								rr = produitvectoriel
								if (sqrt(rr)<=1)then
									exit
								end if
								call produitmat(A,p,Nx,Ny,produitmatriciel,i1,iN,statinfo)
								Ap = produitmatriciel
								deallocate(produitmatriciel)


								call produitvec(Nx,Ny,p,Ap,produitvectoriel,statinfo,me,i1,iN)
								pAp=produitvectoriel
								alpha = rr / pAp
								xn =xn + alpha * p

								r1 = r-alpha * Ap

								call produitvec(Nx,Ny,r1,r1,produitvectoriel,statinfo,me,i1,iN)
								rr1=produitvectoriel
								p=r1+p*(rr1/rr)
								r=r1
							end do



					end subroutine gc


				subroutine produitmat(A,xx,Nx,Ny,produitmatriciel,i1,iN,statinfo)
					implicit none
					integer::Nx,Ny,i,i1,iN,he,me,statinfo,Np,j1,jN
					real,dimension(:,:),allocatable::A
					real,dimension(:),allocatable::xx,pm,produitmatriciel
					allocate(produitmatriciel(Nx*Ny))
					do i=1,Nx*Ny
						if (i==1) then
							produitmatriciel(i)=A(3,i)*xx(i)+A(4,i)*xx(i+1)+A(5,i)*xx(i+Nx)

						else if ((i>1) .and. (i<=Nx)) then
							produitmatriciel(i)=A(2,i)*xx(i-1)+A(3,i)*xx(i)+A(4,i)*xx(i+1)+A(5,i)*xx(i+Nx)

						else if ((i>Nx*Ny-Nx) .and. (i<Nx*Ny)) then
							produitmatriciel(i)= A(1,i)*xx(i-Nx)+A(2,i)*xx(i-1)+A(3,i)*xx(i)+A(4,i)*xx(i+1)

						else if (i==Nx*Ny) then
							produitmatriciel(i)=A(1,i)*xx(i-Nx)+A(2,i)*xx(i-1)+A(3,i)*xx(i)

						else
							produitmatriciel(i)=A(1,i)*xx(i-Nx)+A(2,i)*xx(i-1)+A(3,i)*xx(i)+A(4,i)*xx(i+1)+A(5,i)*xx(i+Nx)
						end if
					end do


				end subroutine produitmat

				subroutine produitvec(Nx,Ny,v,vv,produitvectoriel,statinfo,me,i1,iN)
					implicit none
					integer::Nx,Ny,	i1,iN,i,me,statinfo,Np
					real,dimension(:),allocatable::v,vv
					real:: pv,produitvectoriel

					pv=0
					do i=i1,iN
						pv= pv + v(i)*vv(i)
					end do
					call MPI_ALLREDUCE(pv,produitvectoriel,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,statinfo)
				end subroutine produitvec




		end module gradientconjugue
