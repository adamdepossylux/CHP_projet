module matrices
		use fonctions
		implicit none



		contains
		subroutine matrice(Nx,Ny,hx,hy,A)
			implicit none
			integer::Nx,Ny,i,j,me,Np,j1,jN
			real::hx,hy
			real,dimension(:,:),allocatable::A
			allocate (A(5,Nx*Ny))


			A=0
			A(3,:)=1.
			do j=1,Ny-2
				do i=1,Nx-2
					A(2,i+1+j*Nx)=1./(hx**2)
				end do
				do i=3,Nx
					A(4,i-1+j*Nx)=1./(hx**2)
				end do
				do i=2,Nx-1
					A(1,i+Nx+(j-1)*Nx)=1./(hy**2)
					A(3,i+j*Nx)=-2./(hx**2) -2./(hy**2)
					A(5,i -Nx +(j+1)*Nx)=1./(hy**2)
				end do
			end do
		end subroutine matrice


		subroutine sm(x,y,Nx,Ny,ssm)
			implicit none
			integer::Nx,Ny,i,j,k,me,j1,jN,Np,he,statinfo
			real, dimension(Nx)::x
			real, dimension(Ny)::y

			real, dimension(:),allocatable::ssm
			allocate(ssm(Ny*Nx));

			do j=1,Ny
				do i=1,Nx
					k=i+(j-1)*(Nx)
					if ((j==1) .OR. (j==Ny))  then
						ssm(k)=g(x(i),y(j))
					else if ((i==1) .OR. (i==Nx)) then
						ssm(k)=h(x(i),y(j))
					else
						ssm(k)=f(x(i),y(j))
					end if
				end do
			end do


		end subroutine sm



		subroutine smt(x,y,Nx,Ny,n,t,Ly,Lx,ssmt)
			implicit none
			integer::Nx,Ny,n,i,j,k, statinfo,me,Np,he,j1,jN
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
			integer ::i1,iN,Me,Np,N,r,q
			t=k/(N/Np)
			r=k-(N/Np)*t
			Mef=me+t
		end subroutine test


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
