program programme
	use mpi
	use fonctions
	use matrices
	use gradientconjugue
	implicit none
	integer,parameter :: Nx=100,Ny=100,NxNy=Nx*Ny
	integer::j1,jN,k1,k2,k3,j,i,i1,iN,statinfo,Np,me,cont,diim,Nt,n2,n5
	double precision,parameter :: hx=1.d0/(Nx+1),hy=1.d0/(Ny+1),Lx=1.d0,Ly=1.d0,dt=0.01d0,D=1.d0
	double precision, dimension(:),allocatable :: b,xn,ssm,ssme,produitmatriciel,v,vv,nnz,colonnes,ic,xx
	double precision, dimension(:),allocatable :: solutionexacte,smf,U,Uo,Uinitial
	integer, dimension( MPI_STATUS_SIZE) :: status
	double precision,dimension(:,:),allocatable::A,AA
	double precision:: produitscalaire,temps_debut,temps_fin,temps_fin_max
	double precision,dimension(Nx+2)::x
	double precision, dimension(Ny+2)::y
	double precision, dimension(100)::temps



	call MPI_INIT(statinfo)

	temps_debut= MPI_WTIME()

	call MPI_COMM_RANK(MPI_COMM_WORLD,me,statinfo)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,Np,statinfo)
	call charge(me,Np,Nx*Ny,i1,iN)
	allocate (b(i1:iN))
	allocate (U(i1:iN))
	allocate (Uo(i1:iN))
	allocate(ssm(i1:iN))
	Nt=1
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
		temps(i)=(i-1)*dt
	end do
	b=0
	Uo=0
	!print*,x
	!print*,y
	! appel du vecteur smf "interieur",c'est-à-dire celui contenant la fonction f
	! appel du vecteur initial U contenant au bord les fonctions g et h
  call laplaciencsr(hx,hy,dt,ic,colonnes,nnz,Nx,Ny,i1,iN)
			call sm2(x,y,Nx,Ny,i1,iN,ssm,Uo)
			b=ssm
			do i=1,Nt
				call gc(nnz,colonnes,ic,b,NxNy,Nx,Ny,xn,i1,iN)
				U=xn

				call sm2(x,y,Nx,Ny,i1,iN,ssm,U)
				b=ssm
				deallocate(xn)
			end do
			!print*,me,U
			deallocate(ssm)
			temps_fin= (MPI_WTIME()-temps_debut)

			call MPI_REDUCE (temps_fin,temps_fin_max,1, MPI_DOUBLE_PRECISION , MPI_MAX ,0,MPI_COMM_WORLD ,statinfo)
			if (me == 0) then
				print*,"Temps : ",temps_fin_max," secondes"
			end if
			n2=size(b)
			open(unit=2, file="sol.txt",form="formatted",access="sequential")
				do i=1,n2
					write(2,*)U(i)
				enddo
				close(2)


		call MPI_FINALIZE(statinfo)


end program programme
