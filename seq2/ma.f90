module matrices
		use fonctions
		implicit none



		contains
			subroutine laplaciencsr(hx,hy,dt,ic,colonnes,nnz,Nx,Ny)
			  implicit none
			  integer :: i,n,Nx,Ny,j,Njic,Nj,Njc,i1,iN,d
			  double precision :: dx,dy,dt,hx,hy
				double precision, dimension(5) ::fg
			  double precision, dimension(:), allocatable :: colonnes,nnz,ic,ic1
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

			return
			end subroutine

			subroutine sm1(x,y,Nx,Ny,ssm,dt,u)  !second memebre
				implicit none
				integer::Nx,Ny,i,j,k,me,i1,iN,Np,he,statinfo
				double precision :: dt
				double precision, dimension(Nx+2)::x
				double precision, dimension(Ny+2)::y
				double precision, dimension(Ny*Nx)::u
				!call charge(me, Np, Ny,j1,jN)
				double precision, dimension(:),allocatable::ssm

				allocate(ssm(Nx*Ny))
				ssm=0
				do j=1,Ny
					do i=1,Nx
						k=i+(j-1)*Nx
						if ((j==1) .OR. (j==Ny+2))  then
							ssm(k)=0.d0
						else if ((i==1) .OR. (i==Nx+2)) then
							ssm(k)=0.d0
						else
							ssm(k)=f1(x(i),y(j))*dt+u(k)
						end if
					end do
				end do

			end subroutine sm1

			subroutine sm3(x,y,Nx,Ny,ssm,dt,u,n,Nt)  !second memebre
				implicit none
				integer::Nx,Ny,i,j,k,me,i1,iN,Np,he,statinfo,n,Nt
				double precision :: dt,Lx,Ly
				double precision, dimension(Nx+2)::x
				double precision, dimension(Ny+2)::y
				double precision,dimension(Nt)::t
				double precision, dimension(Ny*Nx)::u
				!call charge(me, Np, Ny,j1,jN)
				double precision, dimension(:),allocatable::ssm

				allocate(ssm(Nx*Ny))
				Lx=1
				Ly=1
				ssm=0
				do j=1,Ny
					do i=1,Nx
						k=i+(j-1)*Nx
						if ((j==1) .OR. (j==Ny+2))  then
							ssm(k)=f(x(i),y(j))
						else if ((i==1) .OR. (i==Nx+2)) then
							ssm(k)=f(x(i),y(j))
						else
							ssm(k)=ff(x(i),y(j),t(n),Ly,Lx)*dt+u(k)
						end if
					end do
				end do

			end subroutine sm3



			subroutine sm2(x,y,Nx,Ny,ssm,dt,u)  !second memebre
				implicit none
				integer::Nx,Ny,i,j,k,me,i1,iN,Np,he,statinfo
				double precision :: dt
				double precision, dimension(Nx+2)::x
				double precision, dimension(Ny+2)::y
				double precision, dimension(Ny*Nx)::u
				!call charge(me, Np, Ny,j1,jN)
				double precision, dimension(:),allocatable::ssm

				allocate(ssm(Nx*Ny))
				ssm=0
				do j=1,Ny
					do i=1,Nx
						k=i+(j-1)*Nx
						if ((j==1) .OR. (j==Ny+2))  then
							ssm(k)=f(x(i),y(j))
						else if ((i==1) .OR. (i==Nx+2)) then
							ssm(k)=f(x(i),y(j))
						else
							ssm(k)=f(x(i),y(j))*dt+u(k)
						end if
					end do
				end do

			end subroutine sm2



end module  matrices
