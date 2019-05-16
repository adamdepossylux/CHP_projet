module matrices
		use mpi
		use fonctions
		implicit none

contains
			subroutine matrice(Nx,Ny,hx,hy,nnz,colonnes,ic,D,i1,iN,dt)
			use mpi
			implicit none
			integer,intent(in)::Nx,Ny
			integer::i,j,cont,diim,Np,me,statinfo,i1,iN
			double precision,intent(in) :: hx,hy,dt,D
			double precision,dimension(:),allocatable::nnz,colonnes,ic
			double precision,dimension(:,:),allocatable::A,AA

			double precision,dimension(Nx+2)::x
			double precision, dimension(Ny+2)::y




				call MPI_COMM_RANK(MPI_COMM_WORLD,me,statinfo)
				call MPI_COMM_SIZE(MPI_COMM_WORLD,Np,statinfo)

				allocate (AA(5,(Nx+2)*(Ny+2)))
				allocate (A((Nx+2)*(Ny+2),(Nx+2)*(Ny+2)))




!print*,D
! construction des 5 diagonales
				AA=0.d0
				AA(3,:)=1.d0
				do j=1,Ny
					do i=1,Nx
						AA(2,i+j*(Nx+2))=-dt*D*(1.d0/(hx**2))
					end do
					do i=3,Nx+2
						AA(4,i+j*(Nx+2))=-dt*D*(1.d0/(hx**2))
					end do
					do i=2,Nx+1
						AA(1,i+(j-1)*(Nx+2))=-dt*D*(1.d0/(hy**2))
						AA(3,i+j*(Nx+2))=1.d0-dt*D*(-2.d0/(hx**2) -2.d0/(hy**2))
						AA(5,i+(j+1)*(Nx+2))=-dt*D*(1.d0/(hy**2))
					end do
				end do


! stockage de ces 5 diagonales dans la matrice de taille (Nx+2)*(Ny+2)
				A=0.d0
				A(1,1)=AA(3,1)
				A((Nx+2)*(Ny+2),(Nx+2)*(Ny+2))=AA(3,(Nx+2)*(Ny+2))
				DO i=2,(Nx+2)*(Ny+2)-1
					A(i,i)=AA(3,i)	;A(i,i-1)=AA(2,i-1) ; A(i,i+1)=AA(4,i+1)
					if (i>=(Nx+2+1)) then
						A(i,i-Nx-2)=AA(1,i-Nx-2)
					end if
					if (i<=(Nx+2)*(Ny+2)-Nx-2-1) then
						A(i,i+Nx+2)=AA(5,i+Nx+2)
					end if

				ENDDO



				deallocate(AA)


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
				allocate (nnz(diim))
				allocate (colonnes(diim))
				cont=0
				allocate (ic( iN-i1+2))

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

			subroutine laplaciencsr(hx,hy,dt,ic1,colonnes,nnz,Nx,Ny,i1,iN)
				  implicit none
				  integer :: i,n,Nx,Ny,j,Njic,Nj,Njc,i1,iN,d
				  double precision :: dx,dy,dt,hx,hy
					double precision, dimension(5) ::fg
				  double precision, dimension(:), allocatable :: colonnes,nnz,ic,ic1
					allocate(ic1(1:iN-i1+2))! +1 ou 2
					allocate(ic(1:nx*ny+1))
					allocate(nnz(1:(2*(2*3+4*(Nx-2))+(Ny-2)*(2*4+5*(Nx-2)))))
					allocate(colonnes(1:(2*(2*3+4*(Nx-2))+(Ny-2)*(2*4+5*(Nx-2)))))
					nnz=0.d0
					ic=0.d0
					colonnes=0.d0

					!disparition de g10 et de g01
					ic(1)=1
					nnz(1:3)=[1+2*dt*(1/(hx**2))+dt*2*(1/(hy**2)),-dt*(1/(hx**2)),-1*dt*(1/(hy**2))]
					colonnes(1:3)=[1,2,2+nx-1]
					ic(2)=ic(1)+3
					!k=i*nx+j
					Nj=4
					Njic=2
					do i=1,nx-2
						!on se balade de l'indice 2 à l'indice nx-1, donc on est en i+1 et non en i
						nnz((i-1)*4+Nj:i*4+nj)=[-dt*(1/(hx**2)),1+2*dt*(1/(hx**2))+dt*2*(1/(hy**2)),-dt*(1/(hx**2)),-1*dt*(1/(hy**2))]
						ic(1+Njic)=ic(Njic)+4
						colonnes((i-1)*4+Nj:i*4+nj)=[i+1-1,i+1,i+1+1,i+1+nx]
						Njic=Njic+1
					end do
					Nj=Nj+4*(nx-2)
					nnz(Nj:Nj+2)=[-dt*(1/(hx**2)),1+2*dt*(1/(hx**2))+dt*2*(1/(hy**2)),-1*dt*(1/(hy**2))]
					colonnes(Nj:Nj+2)=[nx-1,nx,nx+nx]
					ic(Njic+1)=ic(Njic)+3 !i=i+1, à la fin i=nx-1
					Nj=Nj+3
					Njic=Njic+1
					!print*,"njic = ",njic
					fg=[-1*dt*(1/(hy**2)),-dt*(1/(hx**2)),1+2*dt*(1/(hx**2))+dt*2*(1/(hy**2)),-dt*(1/(hx**2)),-1*dt*(1/(hy**2))]
					do j=1,ny-2
							nnz(Nj:Nj+3)=[-1*dt*(1/(hy**2)),1+2*dt*(1/(hx**2))+dt*2*(1/(hy**2)),-dt*(1/(hx**2)),-1*dt*(1/(hy**2))]
							colonnes(Nj:Nj+3)=[j*nx+1-nx,j*nx+1,j*nx+1+1,j*nx+1+nx]
							ic(Njic+1)=ic(Njic)+4
							!k=i*nx+j
							Nj=Nj+4
							Njic=Njic+1
							do i=1,nx-2
								!on se balade de l'indice 2 à l'indice nx-1, donc on est en i+1 et non en i
								nnz((i-1)*5+Nj:i*5+nj)=fg
								ic(Njic+1)=ic(Njic)+5
								colonnes((i-1)*5+Nj:i*5+nj)=[j*nx+i+1-nx,j*nx+i+1-1,j*nx+i+1,j*nx+i+1+1,j*nx+i+1+nx]
								Njic=Njic+1
							end do
							Nj=Nj+5*(nx-2)
							nnz(Nj:Nj+3)=[-1*dt*(1/(hy**2)),-dt*(1/(hx**2)),1+2*dt*(1/(hx**2))+dt*2*(1/(hy**2)),-1*dt*(1/(hy**2))]
							colonnes(Nj:Nj+3)=[j*nx+1+nx-1-nx,j*nx+1+nx-1-1,j*nx+1+nx-1,j*nx+1+nx-1+nx]
							ic(Njic+1)=ic(Njic)+4 !i=i+1, à la fin i=nx-1

							Njic=Njic+1 !décalage pour i
							Nj=Nj+4 !
							Njc=Njc+2 ! valeurs aux bords
					end do

					nnz(Nj:Nj+2)=[-1*dt*(1/(hy**2)),1+2*dt*(1/(hx**2))+dt*2*(1/(hy**2)),-dt*(1/(hx**2))]
					colonnes(Nj:Nj+2)=[(ny-1)*nx+1-nx,(ny-1)*nx+1,(ny-1)*nx+1+1]
					ic(Njic+1)=ic(Njic)+3
					!k=i*nx+j
					Nj=Nj+3
					Njic=Njic+1

					do i=1,nx-2
						!on se balade de l'indice 2 à l'indice nx-1, donc on est en i+1 et non en i
						nnz((i-1)*4+Nj:i*4+nj)=[-1*dt*(1/(hy**2)),-dt*(1/(hx**2)),1+2*dt*(1/(hx**2))+dt*2*(1/(hy**2)),-dt*(1/(hx**2))]
						ic(Njic+1)=ic(njic)+4
						colonnes((i-1)*4+Nj:i*4+nj)=[(ny-1)*nx+1+i-nx,(ny-1)*nx+1+i-1,(ny-1)*nx+1+i,(ny-1)*nx+1+i+1]
						Njic=Njic+1
					end do
					Nj=Nj+4*(nx-2)
					nnz(Nj:Nj+2)=[-1*dt*(1/(hy**2)),-dt*(1/(hx**2)),1+2*dt*(1/(hx**2))+dt*2*(1/(hy**2))]
					colonnes(Nj:Nj+2)=[ny*nx-nx,ny*nx-1,ny*nx]
					ic(Njic+1)=ic(Njic)+3 !i=i+1, à la fin i=nx-1
					!ic(Njic+2)=ic(Njic+1)+1
					do i=0,iN-i1+1
						ic1(i+1)=ic(i1+i)
					end do
				return
				end subroutine

			subroutine laplaciencsr2(hx,hy,dt,ic1,colonnes,nnz,Nx,Ny,i1,iN)
							  implicit none
							  integer :: i,n,Nx,Ny,j,Njic,Nj,Njc,i1,iN
							  double precision :: dx,dy,dt,hx,hy
								double precision, dimension(5) ::fg
							  double precision, dimension(:), allocatable :: colonnes,nnz,ic,ic1,colonnes1,nnz1
								allocate(nnz((nx+2)*2+5*nx*ny+ny*2))
								allocate(colonnes((nx+2)*2+5*nx*ny+ny*2))
								allocate (ic1(1:iN-i1+1))
								allocate(ic(1:(nx+2)*2+nx*ny+ny*2))

								ic1=0
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
								fg=[-dt*(1/(hy**2)),-dt*(1/(hx**2)),1+2*dt*(1/(hx**2))+dt*2*(1/(hy**2)),-dt*(1/(hx**2)),-1*dt*(1/(hy**2))]
								do j=1,ny-1
								    Njic=(j-1)*(nx+1)+1
								    Nj=(j-1)*5*(nx-1)
								    Njc=Njc+2
								    colonnes(nx+1+Nj+Njc-1)=j*(nx+1)+1
								    nnz(nx+1+Nj+Njc-1)=1
										ic(nx+2+Njic)=ic(nx+1+Njic)+1
								    do i=1,nx-1
								        ic(nx+2+Njic+i)=ic(nx+1+Njic+i)+5
												colonnes(nx+1+Nj+Njc+5*(i-1):nx+Nj+Njc+i*5)=[j*(nx+1)+i-nx,j*(nx+1)+i,j*(nx+1)+i+1,j*(nx+1)+i+2,j*(nx+1)+i+nx+2]
												nnz(nx+1+Nj+Njc+5*(i-1):nx+Nj+Njc+i*5)=fg
										end do
								    ic(nx+2+Njic+nx)=ic(nx+2+Njic+nx-1)+1;
								    colonnes(nx+Nj+Njc+(nx-1)*5+1)=j*(nx+1)+nx+1;
								    nnz(nx+Nj+Njc+(nx-1)*5+1)=1;
								end do
								do i=1,nx+1
									ic((nx+1)*(ny+1)-nx-1+i)=5*(nx-1)*(ny-1)+(ny-1)*2+nx+1+i
									nnz(5*(nx-1)*(ny-1)+(ny-1)*2+nx+2+i)=1
									colonnes(5*(nx-1)*(ny-1)+(ny-1)*2+nx+2 + i)=(nx+1)*(ny+1)-nx+i
								end do
								do i=0,iN-i1
									ic1(i+1)=ic(i1+i)
								end do
							return
							end subroutine

			subroutine sm3(x,y,Nx,Ny,i1,iN,ssm,u,t,n,nt,dx,dy,dt)
				use mpi
				implicit none
				integer::Nx,Ny,i,j,k,i1,iN,d,Nt,n
				double precision ::dt,Lx,Ly,dx,dy
				double precision, dimension(Nx+2)::x
				double precision, dimension(Ny+2)::y
				double precision, dimension(Nt)::t
				double precision, dimension(:),allocatable::ssm,u
				!allocate(Uinitial(i1:iN))
				!allocate(U(i1:iN))
				Lx=1
				Ly=1
				ssm=0.d0
				do k=i1,iN
					d=k/Nx
					i= k-d*Nx
					if (i==0) then
						j=d
						i=Nx+2
					else
						j=d+1
					end if
					if ((j==1)) then
						ssm(k)=ff(x(i+1),y(j+1),t(n),Ly,Lx)*dt+u(k)+0.d0
					else if (j==Ny) then
						ssm(k)=ff(x(i+1),y(j+1),t(n),Ly,Lx)*dt+u(k)+0.d0
					else if (i==1) then!gouche
						ssm(k)=ff(x(i+1),y(j+1),t(n),Ly,Lx)*dt+u(k)+dt*1.d0/(dx**2)
					else if (i==Nx) then!droite
						ssm(k)=ff(x(i+1),y(j+1),t(n),Ly,Lx)*dt+u(k)+dt*1.d0/(dy**2)
					else
						ssm(k)=ff(x(i+1),y(j+1),t(n),Ly,Lx)*dt+u(k)
						!ssm(k)=1.d0
					end if
				end do
			end subroutine sm3


			subroutine sm1(x,y,Nx,Ny,ssm,dt,u,dx,dy,i1,iN)
				use mpi
				implicit none
				integer::Nx,Ny,i,j,k,i1,iN,d
				double precision :: dt,dx,dy
				double precision, dimension(Nx+2)::x
				double precision, dimension(Ny+2)::y
				double precision, dimension(:),allocatable::u
				double precision, dimension(:),allocatable::ssm
				!allocate(u(i1:iN))
				ssm=0.d0
				do k=i1,iN
					d=k/Nx!d diviseur euclidien de k par Nx
					i= k-d*Nx!i reste de la division euclidienne de k par Nx
					if (i==0) then!divions entière
						j=d
						i=Nx+2
					else
						j=d+1!k=i+j*Nx
					end if
					if ((j==1)) then
						ssm(k)=u(k)+f1(x(i+1),y(j+1))*dt
					else if (j==Ny) then
						ssm(k)=u(k)+f1(x(i+1),y(j+1))*dt
					else if (i==1) then
						ssm(k)=u(k)+f1(x(1),y(j+1))*dt
					else if (i==Nx) then
						ssm(k)=u(k)+f1(x(Nx+2),y(j+1))*dt
					else
						ssm(k)=u(k)+f1(x(i+1),y(j+1))*dt
					end if
				end do
			end subroutine sm1

			subroutine sm2(x,y,Nx,Ny,ssm,dt,u,dx,dy,i1,iN)
				use mpi
				implicit none
				integer::Nx,Ny,i,j,k,i1,iN,d
				double precision :: dt,dx,dy
				double precision, dimension(Nx+2)::x
				double precision, dimension(Ny+2)::y
				double precision, dimension(:),allocatable::u
				double precision, dimension(:),allocatable::ssm
				!allocate(u(i1:iN))
				ssm=0.d0
				do k=i1,iN
					d=k/Nx
					i= k-d*Nx
					if (i==0) then
						j=d
						i=Nx+2
					else
						j=d+1
					end if
					if ((j==1)) then
						ssm(k)=f(x(i+1),y(j+1))*dt+u(k)+dt*f(x(i+1),y(j))/(dy**2)
					else if (j==Ny) then
						ssm(k)=f(x(i+1),y(j+1))*dt+u(k)+dt*f(x(i+1),y(Ny+2))/(dy**2)
					else if (i==1) then!gouche
						ssm(k)=f(x(i+1),y(j+1))*dt+u(k)+dt*f(x(i),y(j+1))/(dx**2)
					else if (i==Nx) then!droite
						ssm(k)=f(x(i+1),y(j+1))*dt+u(k)+dt*f(x(Nx+2),y(j+1))/(dx**2)
					else
						ssm(k)=f(x(i+1),y(j+1))*dt+u(k)
					end if
				end do
			end subroutine sm2

			subroutine uex1m(x,y,Nx,Ny,i1,iN,u)
				use mpi
				implicit none
				integer::Nx,Ny,i,j,k,i1,iN,d
				double precision :: dt,dx,dy
				double precision, dimension(Nx+2)::x
				double precision, dimension(Ny+2)::y
				double precision, dimension(:),allocatable::u
				!allocate(u(i1:iN))
				!ssm=0.d0
				u=0
				do k=i1,iN
					d=k/Nx
					i= k-d*Nx
					if (i==0) then
						j=d
						i=Nx+2
					else
						j=d+1
					end if
					u(k)=uex1(x(i+1),y(j+1))
				end do
			end subroutine uex1m


			subroutine uex2m(x,y,Nx,Ny,i1,iN,u)
				use mpi
				implicit none
				integer::Nx,Ny,i,j,k,i1,iN,d
				double precision :: dt,dx,dy
				double precision, dimension(Nx+2)::x
				double precision, dimension(Ny+2)::y
				double precision, dimension(:),allocatable::u
				!allocate(u(i1:iN))
				!ssm=0.d0
				u=0
				do k=i1,iN
					d=k/Nx
					i= k-d*Nx
					if (i==0) then
						j=d
						i=Nx+2
					else
						j=d+1
					end if
					u(k)=uex2(x(i+1),y(j+1))
				end do
			end subroutine uex2m

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
