module matrices
		use fonctions
		implicit none



		contains
		subroutine matrice(Nx,Ny,hx,hy,nnz,colonnes,ic,dt)


			implicit none
			integer,intent(in)::Nx,Ny

			integer::i,i1,iN,statinfo,Np,me,j,cont,diim,n6
			double precision,intent(in) :: hx,hy,dt



			real,dimension(:),allocatable::nnz,colonnes,ic
			double precision,dimension(:,:),allocatable::A,AA


			allocate (AA(5,(Nx+2)*(Ny+2)))
			allocate (A((Nx+2)*(Ny+2),(Nx+2)*(Ny+2)))

			call laplacien(1.0/hx,1.0/hy,dt,A,Nx+1,Ny+1)

! nombre d'éléments non nuls dans la matrice
			diim=0
			do j=1,(Nx+2)*(Ny+2)
				do i=1,(Nx+2)*(Ny+2)
					if (A(i,j)/=0) then
						diim=diim+1
					end if
				end do
			end do
! construction de la matrice csr
			allocate (nnz(diim))! matrice valeurs csr
			allocate (colonnes(diim)) ! matrice coordonnées des colonnes
			allocate (ic((Nx+2)*(Ny+2)+1)) !matrice définie avec la somme depuis la première lignes du nombre d'elements non nulle de la ligne
			cont=0
			do i=1,(Nx+2)*(Ny+2)
				ic(i)=(cont+1)
				do j=1,(Nx+2)*(Ny+2)
					if (A(i,j)/=0) then
						cont=cont+1
						nnz(cont)=A(i,j)
						colonnes(cont)=j
					end if
				end do
			end do
			ic((Nx+2)*(Ny+2)+1)=cont+1


		end subroutine matrice



		subroutine sm(x,y,Nx,Ny,i1,iN,ssm)  !second memebre
			implicit none
			integer::Nx,Ny,i,j,k,me,i1,iN,Np,he,statinfo
			real, dimension(Nx+2)::x
			real, dimension(Ny+2)::y
			!call charge(me, Np, Ny,j1,jN)
			real, dimension(:),allocatable::ssm

			allocate(ssm((Nx+2)*(Ny+2)))

			do j=1,Ny+2
				do i=1,Nx+2
					k=i+(j-1)*(Nx+2)
					if ((j==1) .OR. (j==Ny+2))  then
						ssm(k)=g(x(i),y(j))
					else if ((i==1) .OR. (i==Nx+2)) then
						ssm(k)=h(x(i),y(j))
					else
						ssm(k)=f(x(i),y(j))
					end if
				end do
			end do
		end subroutine sm



		subroutine smt(x,y,Nx,Ny,n,t,Ly,Lx,ssmt) !second membre stationnaire et parallélisé
			implicit none
			integer::Nx,Ny,n,i,j,k,statinfo,me,Np,he,j1,jN
			real,dimension(Nx*Ny)::t
			real::Lx,Ly
			real, dimension(Nx)::x
			real, dimension(Ny)::y
			real, dimension(:),allocatable::ssmt

			call charge(me,Np,Ny,j1,jN)
			allocate(ssmt((jN-j1+1)*Nx))

			do j=j1,jN
				do i=1,Nx
					k=i+(j-1)*(Nx)
					if ((j==1) .OR. (j==Ny))  then
						ssmt(k)=gg(x(i),y(j),t(n))
					else if ((i==1) .OR. (i==Nx)) then
						ssmt(k)=hh(x(i),y(j),t(n))
					else
						ssmt(k)=ff(x(i),y(j),t(n),Ly,Lx)
					end if
				end do
			end do

		end subroutine smt


		subroutine test(i1,iN,N,Np,k,Mef,r)!renvoie la processus Mef ou est i(Me)+k et le position exacte du fichier est i_Mef+r
			implicit none
			integer ::i1,iN,Me,Np,N,r,q,k,mef,t
			t=k/(N/Np)
			r=k-(N/Np)*t
			Mef=me+t
		end subroutine test

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
