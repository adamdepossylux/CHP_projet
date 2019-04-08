module matrices
		use fonctions
		implicit none

		
		
		contains	
		subroutine matrice(Nx,Ny,hx,hy,nnz,colonnes,ic)	


			implicit none 
			integer,intent(in)::Nx,Ny
			
			integer::i,i1,iN,statinfo,Np,me,j,cont,diim
			double precision,intent(in) :: hx,hy
			
			

			real,dimension(:),allocatable::nnz,colonnes,ic
			double precision,dimension(:,:),allocatable::A,AA
			
		
			allocate (AA(5,(Nx+2)*(Ny+2)))
			allocate (A((Nx+2)*(Ny+2),(Nx+2)*(Ny+2)))
		
			
		

! construction des 5 diagonales			
			A=0.d0
			
			AA=0.d0
			AA(3,:)=1.d0
			do j=1,Ny
				do i=1,Nx
					AA(2,i+j*(Nx+2))=1.d0/(hx**2)
				end do
				do i=3,Nx+2
					AA(4,i+j*(Nx+2))=1.d0/(hx**2)
				end do	
				do i=2,Nx+1
					AA(1,i+(j-1)*(Nx+2))=1.d0/(hy**2)
					AA(3,i+j*(Nx+2))=-2.d0/(hx**2) -2.d0/(hy**2)
					AA(5,i+(j+1)*(Nx+2))=1.d0/(hy**2)
				end do	
			end do
			A(1,1)=AA(3,1)
			A((Nx+2)*(Ny+2),(Nx+2)*(Ny+2))=AA(3,(Nx+2)*(Ny+2))	


! stockage de ces 5 diagonale dans la matrice de taille Nx*Ny
			DO i=2,(Nx+2)*(Ny+2)-1
				A(i,i)=AA(3,i)	;A(i,i-1)=AA(2,i-1) ; A(i,i+1)=AA(4,i+1)
				if (i>=(Nx+2)) then
					A(i,i-Nx-2)=AA(1,i-Nx-2)
				end if
				if (i<=(Nx+2)*(Ny+2)-Nx-2) then  
					A(i,i+Nx+2)=AA(5,i+Nx+2)
				end if
				
			ENDDO
			
			
			
			deallocate(AA)
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
			allocate (nnz(diim))
			allocate (colonnes(diim))
			allocate (ic((Nx+2)*(Ny+2)+1))
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
		
		
		subroutine sm(x,y,Nx,Ny,i1,iN,ssm)  
			implicit none
			integer::Nx,Ny,i,j,k,me,i1,iN,Np,he,statinfo
			real, dimension(Nx+2)::x
			real, dimension(Ny+2)::y

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
	
  