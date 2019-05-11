		module gradientconjugue
			use mpi
			use matrices
			implicit none
			contains


						subroutine gc(nnz,colonnes,ic,b,NxNy,Nx,Ny,xn,i1,iN)
						use mpi
						use matrices
							implicit none
							double precision , dimension(:),allocatable ::b,xn, p, Ap, r,r1, produitmatriciel,nnz,colonnes,ic


							double precision :: alpha, rr, rr1, pAp,produitscalaire
							integer :: j,i,k,statinfo,Np,me, NxNy,Nx,Ny,i1,iN

							call MPI_COMM_RANK(MPI_COMM_WORLD,me,statinfo)
							call MPI_COMM_SIZE(MPI_COMM_WORLD,Np,statinfo)
							allocate (r(i1:iN))
							allocate (r1(i1:iN))
							allocate (p(i1:iN))
							!allocate(xn(i1:iN))
							allocate(Ap(i1:iN))


							xn=0.d0
							r(i1:iN) = b(i1:iN)
							p(i1:iN) = b(i1:iN)
							call produitsca(Nx,Ny,r,r,produitscalaire,statinfo,me,i1,iN)
							rr = produitscalaire
							DO k = 1, NxNy+1
								call produitsca(Nx,Ny,r,r,produitscalaire,statinfo,me,i1,iN)
								rr = produitscalaire
								if (sqrt(rr)<=0.0000000000001)then
									exit
								end if
								call produitmat(nnz,colonnes,ic,p,Nx,Ny,produitmatriciel,i1,iN,statinfo)
								Ap(i1:iN) = produitmatriciel
								deallocate(produitmatriciel)


								call produitsca(Nx,Ny,p,Ap,produitscalaire,statinfo,me,i1,iN)
								pAp=produitscalaire
								alpha = rr / pAp

								xn(i1:iN) =xn(i1:iN) + alpha * p(i1:iN)

								r1(i1:iN) = r(i1:iN)-alpha * Ap(i1:iN)

								call produitsca(Nx,Ny,r1,r1,produitscalaire,statinfo,me,i1,iN)
								rr1=produitscalaire
								p(i1:iN)=r1(i1:iN)+(p(i1:iN))*(rr1/rr)
								r(i1:iN)=r1(i1:iN)
							end do



					end subroutine gc


				subroutine produitmat(nnz,colonnes,ic,xx,Nx,Ny,produitmatriciel,i1,iN,statinfo)
					use mpi
					use matrices
					implicit none
					integer::Nx,Ny,i,i1,iN,me,statinfo,Np,j1,jN,k1,k2,k3,j,tag4,tag6,request
					double precision,dimension(:,:),allocatable::A
					double precision,dimension(:),allocatable::xx,produitmatriciel,nnz,colonnes,ic,pp
					integer,dimension(MPI_STATUS_SIZE):: status
					call MPI_COMM_RANK(MPI_COMM_WORLD,me,statinfo)
					call MPI_COMM_SIZE(MPI_COMM_WORLD,Np,statinfo)
					allocate(produitmatriciel((iN-i1+1)))
					produitmatriciel=0
					tag4=10;tag6=120
					allocate(pp((Nx+2)*(Ny+2)))
					pp=0
					pp(i1:iN)=xx(i1:iN)

				if(Np>1) then
					if (me==Np-1) then
						call MPI_ISEND(xx(i1:Nx+2),Nx+2,MPI_REAL,me-1,tag4,MPI_COMM_WORLD,request,statinfo )
						call MPI_IRECV(pp(i1-Nx-2:i1-1),Nx+2,MPI_REAL,me-1,tag6,MPI_COMM_WORLD,request,statinfo )
						call MPI_WAIT(request,status,statinfo)


					else if(me==0) then
						call MPI_ISEND(xx(iN-Nx-2+1:iN),Nx+2,MPI_REAL,me+1,tag6,MPI_COMM_WORLD,request,statinfo )
						call MPI_IRECV(pp(iN+1:iN+Nx+2),Nx+2,MPI_REAL,me+1,tag4,MPI_COMM_WORLD,request,statinfo)
						call MPI_WAIT(request,status,statinfo)


					else

						call MPI_ISEND(xx(iN-Nx-2+1:iN),Nx+2,MPI_REAL,me+1,tag6,MPI_COMM_WORLD,request,statinfo )
						call MPI_IRECV(pp(i1-Nx-2:i1-1),Nx+2,MPI_REAL,me-1,tag6,MPI_COMM_WORLD,request,statinfo )
						call MPI_WAIT(request,status,statinfo)
						call MPI_ISEND(xx(i1:Nx+2),Nx+2,MPI_REAL,me-1,tag4,MPI_COMM_WORLD,request,statinfo )
						call MPI_IRECV(pp(iN+1:iN+Nx+2),Nx+2,MPI_REAL,me+1,tag4,MPI_COMM_WORLD,request,statinfo )
						call MPI_WAIT(request,status,statinfo)
					end if
				end if
					call MPI_BARRIER(MPI_COMM_WORLD,statinfo)
					do i=i1,iN
						k1=ic( i-i1+1)
						k2=ic(i-i1+2)-1
						do j=k1,k2
							k3=(colonnes(j))
							produitmatriciel((i-i1+1))=produitmatriciel(i-i1+1) +nnz(j)*pp(k3)
						end do
					end do
					deallocate(pp)
				end subroutine produitmat

				subroutine produitsca(Nx,Ny,v,vv,produitscalaire,statinfo,me,i1,iN)
					use mpi
					use matrices
					implicit none
					integer::Nx,Ny,	i1,iN,i,me,statinfo,Np
					double precision,dimension(:),allocatable::v,vv
					double precision:: produitscalaire,pv
					call MPI_COMM_RANK(MPI_COMM_WORLD,me,statinfo)
					call MPI_COMM_SIZE(MPI_COMM_WORLD,Np,statinfo)
					pv=0.d0
					do i=i1,iN
						pv= pv + v(i)*vv(i)
					end do
					call MPI_ALLREDUCE(pv,produitscalaire,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,statinfo)
				end subroutine produitsca




		end module gradientconjugue
