			program programme
				use mpi
				use fonctions
				use matrices
				use gradientconjugue
				implicit none
				integer::D
				integer,parameter :: Nx=1,Ny=1,NxNy=(Nx+2)*(Ny+2)			
				integer::j1,jN,k1,k2,k3,j,i,i1,iN,statinfo,Np,me,cont,diim
				double precision,parameter :: hx=1.d0/(Nx+1),hy=1.d0/(Ny+1),Lx=1.,Ly=1.,dt=1.
				double precision, dimension(:),allocatable :: b,xn,ssm,ssme,produitmatriciel,v,vv,nnz,colonnes,ic,xx
				double precision, dimension(:),allocatable :: solutionexacte,smf,solution,Uinitial
				integer, dimension( MPI_STATUS_SIZE) :: status
				double precision,dimension(:,:),allocatable::A,AA
				double precision:: produitscalaire,delta,temps_debut,temps_fin,temps_fin_max
				double precision,dimension(Nx+2)::x
				double precision, dimension(Ny+2)::y
				double precision, dimension(100)::temps

				delta=0.01d0


			
				call MPI_INIT(statinfo)
				
				temps_debut= MPI_WTIME()

				call MPI_COMM_RANK(MPI_COMM_WORLD,me,statinfo)
				call MPI_COMM_SIZE(MPI_COMM_WORLD,Np,statinfo)
				call charge(me,Np,(Nx+2)*(Ny+2),i1,iN)
				allocate (b(i1:iN))
				allocate (solution(i1:iN))
	

				do i=1,Nx+2
					x(i)=(i-1)*hx
				end do
				do i=1,Ny+2
					y(i)=(i-1)*hy
				end do				
				do i=1,100
					temps(i)=(i-1)*delta
				end do
				call ajoutf(x,y,Nx,Ny,i1,iN,smf,temps(1),Ly,Lx)
				call sme(x,y,Nx,Ny,i1,iN,Uinitial)

				call matrice(Nx,Ny,hx,hy,nnz,colonnes,ic,D,i1,iN)
		
				b=Uinitial  + delta*smf
				do i=1,10		
					call gc(nnz,colonnes,ic,b,NxNy,Nx,Ny,xn,i1,iN)
					deallocate(smf)
					call ajoutf(x,y,Nx,Ny,i1,iN,smf,temps(i),Ly,Lx)
					b=xn+delta*smf			
					solution=xn
					deallocate(xn)
				end do
				temps_fin= (MPI_WTIME()-temps_debut)
				
				
				call MPI_REDUCE (temps_fin,temps_fin_max,1, MPI_DOUBLE_PRECISION , MPI_MAX ,0,MPI_COMM_WORLD ,statinfo)
				if (me == 0) then 
					print*,"Temps : ",temps_fin_max," secondes"
				end if
				
				
			call MPI_FINALIZE(statinfo)			
				
				
end program programme			