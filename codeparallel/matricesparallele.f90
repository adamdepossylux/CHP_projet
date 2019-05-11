module matrices
		use mpi
		use fonctions
		implicit none

contains
			subroutine matrice(Nx,Ny,hx,hy,nnz,colonnes,ic,D,i1,iN,dt)
			use mpi
			implicit none
			integer,intent(in)::Nx,Ny,D
			integer::i,j,cont,diim,Np,me,statinfo,i1,iN,n4
			double precision,intent(in) :: hx,hy,dt
			double precision,parameter::delta=0.01d0
			double precision,dimension(:),allocatable::nnz,colonnes,ic
			double precision,dimension(:,:),allocatable::A,AA

			double precision,dimension(Nx+2)::x
			double precision, dimension(Ny+2)::y
			double precision, dimension(100)::temps




				call MPI_COMM_RANK(MPI_COMM_WORLD,me,statinfo)
				call MPI_COMM_SIZE(MPI_COMM_WORLD,Np,statinfo)

				allocate (AA(5,(Nx+2)*(Ny+2)))
				allocate (A((Nx+2)*(Ny+2),(Nx+2)*(Ny+2)))
				call laplacien((1.0/hx),(1.0/hy),dt,A,Nx,Ny)

! axe des abscisses
				do i=1,Nx+2
					x(i)=(i-1)*hx
				end do
! axe des ordonnées
				do i=1,Ny+2
					y(i)=(i-1)*hy
				end do
! échelle de temps
				do i=1,100
					temps(i)=(i-1)*delta
				end do




! nombre d'éléments non nuls dans la matrice


				diim=0
				do i=i1,iN
					do j=1,(Nx+2)*(Ny+2)
						if (A(i,j)/=0) then
							diim=diim+1
						end if
					end do
				end do
! construction de la matrice csr
				!allocate (nnz(diim))
				!allocate (colonnes(diim))
				cont=0
				!allocate (ic( iN-i1+2))

				do i=i1,iN
					ic(i-i1+1)=(cont+1)
					do j=1,(Nx+2)*(Ny+2)
						if (A(i,j)/=0) then
							cont=cont+1
							nnz(cont)=A(i,j)
							colonnes(cont)=j
						end if
					end do
				end do
				ic(iN-i1+2)=cont+1

				deallocate(A)

			end subroutine matrice


			subroutine laplaciencsrp(hx,hy,dt,ic1,colonnes,nnz,Nx,Ny,i1,iN)
			  implicit none
			  integer :: i,n,Nx,Ny,j,Njic,Nj,Njc,i1,iN
			  double precision :: dx,dy,dt,hx,hy
				double precision, dimension(5) ::fg
			  double precision, dimension(:), allocatable :: colonnes,nnz,ic,ic1
				allocate(ic(1:(nx+2)*2+nx*ny+ny*2))
				!hx=1.0/nx;
				!hy=1.0/ny;
				!print*,"pq"
				!print*,"vr"
				ic1=0
				nnz=0
				ic=0
				colonnes=0
				do i=1,nx+1
					nnz(i)=1
					ic(i)=i
					colonnes(i)=i
				end do
				if(i1<nx+1) then
				do i=1,min(nx+2,iN-i1)
					nnz(i)=1
					ic(i)=i
					colonnes(i)=i
				end do
				endif
				if(iN)
				Njc=0
				Njc=0
				Nj=0
				!print*,"salit"
				fg=[-dt*(1/(hy**2)),-dt*(1/(hx**2)),1+2*dt*(1/(hx**2))+dt*2*(1/(hy**2)),-dt*(1/(hx**2)),-1*dt*(1/(hy**2))]
				do j=max(1,i1),min(ny-1,int(iN/5.))
				    Njic=(j-1)*(nx+1)+1
				    Nj=(j-1)*5*(nx-1)
				    Njc=Njc+2
				    !ic(nx+2+Njic)=ic(nx+1+Njic)+1
						!ic(nx+2+Njic)=ic(nx+1+Njic)+1
				    colonnes(nx+1+Nj+Njc-1)=j*(nx+1)+1
				    nnz(nx+1+Nj+Njc-1)=1
						!print*,"hey"
				    do i=1,nx-1

				        ic(nx+2+Njic+i)=ic(nx+1+Njic+i)+5
								colonnes(nx+1+Nj+Njc+5*(i-1):nx+Nj+Njc+i*5)=[j*(nx+1)+i-nx,j*(nx+1)+i,j*(nx+1)+i+1,j*(nx+1)+i+2,j*(nx+1)+i+nx+2]
								nnz(nx+1+Nj+Njc+5*(i-1):nx+Nj+Njc+i*5)=fg
						end do
				    ic(nx+2+Njic+nx)=ic(nx+2+Njic+nx-1)+1;
				    colonnes(nx+Nj+Njc+(nx-1)*5+1)=j*(nx+1)+nx+1;
				    nnz(nx+Nj+Njc+(nx-1)*5+1)=1;
				end do
				!do i=0,nx+1
			!		ic((nx+1)*(ny+1)-nx-1+i)=5*(nx-1)*(ny-1)+(ny-1)*2+nx+1+i
				!	nnz(5*(nx-1)*(ny-1)+(ny-1)*2+nx+2+i)=1
				!	colonnes(5*(nx-1)*(ny-1)+(ny-1)*2+nx+2 + i)=(nx+1)*(ny+1)-nx+i
				!end do
				do i=0,nx
					ic((nx+1)*(ny+1)-nx-1+i)=5*(nx-1)*(ny-1)+(ny-1)*2+nx+1+i
					nnz(5*(nx-1)*(ny-1)+(ny-1)*2+nx+2+i)=1
					colonnes(5*(nx-1)*(ny-1)+(ny-1)*2+nx+2 + i)=(nx+1)*(ny+1)-nx+i
				end do
				do i=0,iN-i1
					ic1(i+1)=ic(i1+i)
				end do
			return
			end subroutine



			subroutine laplaciencsr(hx,hy,dt,ic1,colonnes,nnz,Nx,Ny,i1,iN)
			  implicit none
			  integer :: i,n,Nx,Ny,j,Njic,Nj,Njc,i1,iN
			  double precision :: dx,dy,dt,hx,hy
								double precision, dimension(5) ::fg
			  double precision, dimension(:), allocatable :: colonnes,nnz,ic,ic1
				allocate(ic(1:(nx+2)*2+nx*ny+ny*2))
				!hx=1.0/nx;
				!hy=1.0/ny;
				!print*,"pq"
				!print*,"vr"
				ic1=0
				nnz=0
				ic=0
				colonnes=0
				do i=1,nx+1
					nnz(i)=1
					ic(i)=i
					colonnes(i)=i
				end do
				do i=1,nx+2
					nnz(i)=1
					ic(i)=i
					colonnes(i)=i
				end do
				Njc=0
				Njc=0
				Nj=0
				!print*,"salit"
				fg=[-dt*(1/(hy**2)),-dt*(1/(hx**2)),1+2*dt*(1/(hx**2))+dt*2*(1/(hy**2)),-dt*(1/(hx**2)),-1*dt*(1/(hy**2))]
				do j=1,ny-1
				    Njic=(j-1)*(nx+1)+1
				    Nj=(j-1)*5*(nx-1)
				    Njc=Njc+2
				    !ic(nx+2+Njic)=ic(nx+1+Njic)+1
						!ic(nx+2+Njic)=ic(nx+1+Njic)+1
				    colonnes(nx+1+Nj+Njc-1)=j*(nx+1)+1
				    nnz(nx+1+Nj+Njc-1)=1
						!print*,"hey"
				    do i=1,nx-1

				        ic(nx+2+Njic+i)=ic(nx+1+Njic+i)+5
								colonnes(nx+1+Nj+Njc+5*(i-1):nx+Nj+Njc+i*5)=[j*(nx+1)+i-nx,j*(nx+1)+i,j*(nx+1)+i+1,j*(nx+1)+i+2,j*(nx+1)+i+nx+2]
								nnz(nx+1+Nj+Njc+5*(i-1):nx+Nj+Njc+i*5)=fg
						end do
				    ic(nx+2+Njic+nx)=ic(nx+2+Njic+nx-1)+1;
				    colonnes(nx+Nj+Njc+(nx-1)*5+1)=j*(nx+1)+nx+1;
				    nnz(nx+Nj+Njc+(nx-1)*5+1)=1;
				end do
				!do i=0,nx+1
			!		ic((nx+1)*(ny+1)-nx-1+i)=5*(nx-1)*(ny-1)+(ny-1)*2+nx+1+i
				!	nnz(5*(nx-1)*(ny-1)+(ny-1)*2+nx+2+i)=1
				!	colonnes(5*(nx-1)*(ny-1)+(ny-1)*2+nx+2 + i)=(nx+1)*(ny+1)-nx+i
				!end do
				do i=0,nx
					ic((nx+1)*(ny+1)-nx-1+i)=5*(nx-1)*(ny-1)+(ny-1)*2+nx+1+i
					nnz(5*(nx-1)*(ny-1)+(ny-1)*2+nx+2+i)=1
					colonnes(5*(nx-1)*(ny-1)+(ny-1)*2+nx+2 + i)=(nx+1)*(ny+1)-nx+i
				end do
				do i=0,iN-i1
					ic1(i+1)=ic(i1+i)
				end do
			return
			end subroutine

					subroutine laplaciencsr2(dx,dy,dt,ic1,colonnes,nnz,Nx,Ny,i1,iN)
					  implicit none
					  integer :: i,n,Nx,Ny,j,Njic,Nj,Njc,i1,iN
					  double precision :: dx,dy,dt,hx,hy
					  double precision, dimension(:), allocatable :: colonnes,nnz,ic,ic1
						allocate(ic(1:(nx+2)*2+nx*ny+ny*2))
						hx=1.0/nx;
						hy=1.0/ny;
						!print*,"pq"
						!print*,"vr"
						nnz=0
						ic=0
						colonnes=0
						do i=1,nx+1
							nnz(i)=1
							ic(i)=i
							colonnes(i)=i
						end do

						Njc=0
						Njc=0
						Nj=0
						!print*,"salit"
						do j=1,ny-1
						    Njic=(j-1)*(nx+1)
						    Nj=(j-1)*5*(nx-1)
						    Njc=Njc+2
						    ic(nx+2+Njic)=ic(nx+1+Njic)+1
						    colonnes(nx+1+Nj+Njc-1)=j*(nx+1)+1
						    nnz(nx+1+Nj+Njc-1)=1
								!print*,"hey"
						    do i=1,nx-1

						        ic(nx+2+Njic+i)=ic(nx+1+Njic+i)+5

						        colonnes(nx+1+Nj+Njc+5*(i-1))=j*(nx+1)+i-nx
										colonnes(nx+1+Nj+Njc+5*(i-1)+1)=j*(nx+1)+i
										colonnes(nx+1+Nj+Njc+5*(i-1)+2)=j*(nx+1)+i+1
										colonnes(nx+1+Nj+Njc+5*(i-1)+3)=j*(nx+1)+i+2
										colonnes(nx+1+Nj+Njc+5*(i-1)+4)=j*(nx+1)+i+nx+2

									  nnz(nx+1+Nj+Njc+5*(i-1))=-dt*(1/(hy**2))
										nnz(nx+1+Nj+Njc+5*(i-1)+1)=-dt*(1/(hx**2))
										nnz(nx+1+Nj+Njc+5*(i-1)+2)=1+dt*2*(1/(hx**2))+dt*2/(1/(hy**2))
										nnz(nx+1+Nj+Njc+5*(i-1)+3)=-dt*(1/(hx**2))
										nnz(nx+1+Nj+Njc+5*(i-1)+4)=-dt*(1/(hy**2))
								end do
						    ic(nx+2+Njic+nx)=ic(nx+2+Njic+nx-1)+1;
						    colonnes(nx+Nj+Njc+(nx-1)*5+1)=j*(nx+1)+nx+1;
						    nnz(nx+Nj+Njc+(nx-1)*5+1)=1;
						end do
						do i=0,nx+1
							ic((nx+1)*(ny+1)-nx-1+i)=5*(nx-1)*(ny-1)+(ny-1)*2+nx+1+i
							nnz(5*(nx-1)*(ny-1)+(ny-1)*2+nx+2+i)=1
							colonnes(5*(nx-1)*(ny-1)+(ny-1)*2+nx+2 + i)=(nx+1)*(ny+1)-nx+i
						end do
						!print*,"i1 = ",i1,"iN = ",iN
						do i=0,iN-i1
							ic1(i+1)=ic(i1+i)
						end do
					return
					end subroutine

			subroutine sme(x,y,Nx,Ny,i1,iN,Uinitial)
				use mpi
				implicit none
				integer::Nx,Ny,i,j,k,i1,iN,d
				double precision, dimension(Nx+2)::x
				double precision, dimension(Ny+2)::y
				double precision, dimension(:),allocatable::Uinitial
				allocate(Uinitial(i1:iN))
				do k=i1,iN
					d=k/(Nx+2)
					i= k-d*(Nx+2)
					if (i==0) then
						j=d
						i=Nx+2
					else
						j=d+1
					end if

					if ((j==1) .OR. (j==Ny+2))  then
						Uinitial(k)=g(x(i),y(j))
					else if ((i==1) .OR. (i==Nx+2)) then
						Uinitial(k)=h(x(i),y(j))
					else
						Uinitial(k)=0.d0
					end if
				end do
			end subroutine sme


			subroutine ajoutf(x,y,Nx,Ny,i1,iN,smf,t,Ly,Lx)
				use mpi
				implicit none
				integer::Nx,Ny,i,j,k,i1,iN,d
				double precision::t,Ly,Lx
				double precision, dimension(Nx+2)::x
				double precision, dimension(Ny+2)::y
				double precision, dimension(:),allocatable::smf
				do k=i1,iN
					d=k/(Nx+2)
					i= k-d*(Nx+2)
					if (i==0) then
						j=d
						i=Nx+2
					else
						j=d+1
					end if
					if ((j==1) .OR. (j==Ny+2))  then
						smf(k)=0.d0
					else if ((i==1) .OR. (i==Nx+2)) then
						smf(k)=0.d0
					else
						smf(k)=ff(x(i),y(j),t,Ly,Lx)
					end if
					!print*,k
				end do

			end subroutine ajoutf

			subroutine laplacien(dx,dy,dt,M,Nx,Ny)
			  implicit none
			  integer :: i,Np,n,Nx,Ny
			  double precision :: dx,dy,dt
			  double precision, dimension(:,:), allocatable :: A,B,K,M
			  dx=1/real(Nx)
			  dy=1/real(Ny)
			  allocate(A(1:Ny+1,1:Nx+1));allocate(B(1:Ny+1,1:Nx+1));allocate(K(1:Ny+1,1:3*(Nx+1)));!allocate(M(1:Np**2,1:Np**2));
			  M=0
			  A=0
			  K=0
			  A(1,1)=1;
			  do i=2,Nx!proche x et y loin
			    A(i,i-1)=-dt*(1/(dx**2))
			    A(i,i)=1+dt*2/(dx**2)+dt*2/(dy**2)
			    A(i,i+1)=-dt*(1/(dx**2))
			  enddo
			  B=0;
			  A(Nx+1,Nx+1)=1
			  do i=2,Nx
			    B(i,i)=-dt/(dy**2)
			  enddo
			  do i=1,Nx+1
			      M(i,i)=1
			      M(i+Nx*(Nx+1),i+Ny*(Ny+1))=1
			  enddo
			  K(:,1:(Ny+1))=B
			  K(:,Ny+2:2*(Ny+1))=A
			  K(:,2*(Ny+1)+1:3*(Ny+1))=B
			  do i=1,Nx-1
			    M(i*(Nx+1)+1:(i+1)*(Nx+1),(i-1)*(Ny+1)+1:(i+2)*(Ny+1))=K
			  enddo
			  deallocate(A)
			  deallocate(B)
			  deallocate(K)
			return
			end subroutine




			subroutine charge(me,Np,N,i1,iN)
				implicit none
				integer ::i1,iN,Me,Np,N,r,q
				q=N/np
				r=N-q*np
				if (me<r) then
					i1=me*(q+1) +1
					iN=(me+1)*(q+1)
				else
					i1= 1 + r + me*q
					iN= i1 +q -1
				endif
			end subroutine charge


end module  matrices
