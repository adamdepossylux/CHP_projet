module matrices
		use mpi
		use fonctions
		implicit none

contains
			subroutine matrice(Nx,Ny,hx,hy,nnz,colonnes,ic,D,i1,iN)
			use mpi
			implicit none 
			integer,intent(in)::Nx,Ny,D
			integer::i,j,cont,diim,Np,me,statinfo,i1,iN
			double precision,intent(in) :: hx,hy
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

! construction des 5 diagonales					
				AA=0.d0
				AA(3,:)=1
				do j=1,Ny
					do i=1,Nx
						AA(2,i+j*(Nx+2))=-delta*(1.d0/(hx**2))
					end do
					do i=3,Nx+2
						AA(4,i+j*(Nx+2))=-delta*(1.d0/(hx**2))
					end do	
					do i=2,Nx+1
						AA(1,i+(j-1)*(Nx+2))=-delta*(1.d0/(hy**2))
						AA(3,i+j*(Nx+2))=1.d0-delta*(-2.d0/(hx**2) -2.d0/(hy**2))
						AA(5,i+(j+1)*(Nx+2))=-delta*(1.d0/(hy**2))
					end do	
				end do
				

! stockage de ces 5 diagonales dans la matrice de taille (Nx+2)*(Ny+2)
				A=0.d0
				A(1,1)=AA(3,1)
				A((Nx+2)*(Ny+2),(Nx+2)*(Ny+2))=AA(3,(Nx+2)*(Ny+2))
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
				allocate(smf(i1:iN))	
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
					
				end do
			end subroutine ajoutf	
			
		
			
	
			
		

		
		
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
	
  