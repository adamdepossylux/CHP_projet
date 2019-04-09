program new22
				use mpi
				implicit none
				integer,parameter :: Nx=2,Ny=2,NxNy=(Nx+2)*(Ny+2)
				
				integer::po, j1,jN,k1,k2,k3,j,tag,tag2,jj,tag3,tag4,tag5,tag6,tag7,tag8,ppp,mevoulu,i,i1,iN,statinfo,Np,me,cont,diim,he
				real,parameter :: D=1.,hx=1./(Nx+1),hy=1./(Ny+1),Lx=1.,Ly=1.,dt=1.
				real, dimension(:),allocatable :: xn,ssm,produitmatriciel,v,vv,nnz,colonnes,ic,xx,pm,pp
				integer, dimension( MPI_STATUS_SIZE) :: status
			double precision,dimension(:,:),allocatable::A,AA

				real,dimension(Nx+2)::x
				real, dimension(Ny+2)::y
				real::produitvectoriel
			
						ppp=5


			
				call MPI_INIT(statinfo)
				call MPI_COMM_RANK(MPI_COMM_WORLD,me,statinfo)
				call MPI_COMM_SIZE(MPI_COMM_WORLD,Np,statinfo)
				
			
			allocate (AA(5,(Nx+2)*(Ny+2)))
			allocate (A((Nx+2)*(Ny+2),(Nx+2)*(Ny+2)))
			call charge(me,Np,(Nx+2)*(Ny+2),i1,iN)
	
		allocate(produitmatriciel((iN-i1+1)))
							
					tag=100
					produitmatriciel=0

	
	
				do i=1,Nx+2
					x(i)=(i-1)*hx
				end do
				do i=1,Ny+2
					y(i)=(i-1)*hy
				end do
				
		

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

				
		tag=100
		tag2=50
		tag3=25
		tag4=10
		tag8=90
		tag7=100
		tag6=120

			 allocate(xx(i1:iN))
			 xx(i1:iN)=1
			allocate(pp(2*(iN-i1+1)))
			pp=0
			mevoulu=0
			he=0
			pp(i1:iN)=xx(i1:iN)

		if (me==Np) then
										
				call MPI_RECV(pp(i1-Nx-2:i1-1),Nx+2,MPI_REAL,me-1,tag3,MPI_COMM_WORLD,status,statinfo )

				call MPI_SEND(xx(i1:Nx+2),Nx+2,MPI_REAL,me-1,tag4,MPI_COMM_WORLD,statinfo )
		
		else if(me==0) then
					call MPI_SEND(xx(iN-Nx-2+1:iN),Nx+2,MPI_REAL,me+1,tag3,MPI_COMM_WORLD,statinfo )



					call MPI_RECV(pp(iN+1:iN+Nx+2),Nx+2,MPI_REAL,me+1,tag4,MPI_COMM_WORLD,status,statinfo)
		
		
		else if (modulo(me,2)==1) then
		
			call MPI_RECV(pp(iN+1:iN+Nx+2),Nx+2,MPI_REAL,me-1,tag7,MPI_COMM_WORLD,status,statinfo )
			call MPI_RECV(pp(i1-Nx-2:i1-1),Nx+2,MPI_REAL,me-1,tag6,MPI_COMM_WORLD,status,statinfo )

				
				call MPI_SEND(xx(i1:Nx+2),Nx+2,MPI_REAL,me-1,tag5,MPI_COMM_WORLD,statinfo )
				call MPI_SEND(xx(iN-Nx-2:iN),Nx+2,MPI_REAL,me+1,tag8,MPI_COMM_WORLD,statinfo )

		else if (modulo(me,2)==0) then
				
				call MPI_SEND(xx(i1:Nx+2),Nx+2,MPI_REAL,me-1,tag7,MPI_COMM_WORLD,statinfo )
				call MPI_SEND(xx(iN-Nx-2:iN),Nx+2,MPI_REAL,me+1,tag6,MPI_COMM_WORLD,statinfo )

				
				call MPI_RECV(pp(i1-Nx-2:i1-1),Nx+2,MPI_REAL,me-1,tag8,MPI_COMM_WORLD,status,statinfo )
				call MPI_RECV(pp(iN+1:iN+Nx+2),Nx+2,MPI_REAL,me-1,tag5,MPI_COMM_WORLD,status,statinfo )

	
		
		
		
		end if
			
			
			call MPI_BARRIER(MPI_COMM_WORLD,statinfo)
			do i=i1,iN
						
				k1=ic( i-i1+1)
				k2=ic(i-i1+2)-1
				do j=k1,k2
					k3=(colonnes(j))
					if (i==6) then
						print*,me,k3
					
					end if
		
				
						produitmatriciel((i-i1+1))=produitmatriciel(i-i1+1) +nnz(j)*pp(k3)

						

					
			end do
			end  do
			
			
			
			call MPI_FINALIZE(statinfo)
			print*,me,'pm' , produitmatriciel
			contains
			
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
		
			
			
			end program
