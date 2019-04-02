			program vr
				use mpi
				USE gradientconjugue
				use matrices
				implicit none
				integer,parameter :: Nx=3,Ny=3,NxNy=Nx*Ny
				integer::i,i1,iN,statinfo,Np,me
				real,parameter :: D=1.,hx=1./(Nx-1),hy=1./(Ny-1),Lx=1.,Ly=1.,dt=1.
				real, dimension(Nx*Ny) :: t
				real, dimension(:),allocatable :: xn,ssm,produitmatriciel,v,vv
				real, dimension(:,:),allocatable :: A
				real,dimension(Nx)::x
				real, dimension(Ny)::y
				integer, dimension(MPI_STATUS_SIZE)::status
				real::produitvectoriel

			
				call MPI_INIT(statinfo)
				call MPI_COMM_RANK(MPI_COMM_WORLD,me,statinfo)
				call MPI_COMM_SIZE(MPI_COMM_WORLD,np,statinfo)
				call charge(me,Np,NxNy,i1,iN)

				
				allocate(v(i1:iN));allocate(vv(i1:iN))
				v=2
				vv=3
				do i=1,Nx
					x(i)=(i-1)*hx
				end do
				do i=1,Ny
					y(i)=(i-1)*hy
				end do
				do i=1,Nx*Ny 
					t(i)=(i-1)*dt
				end do
				call produitvec(Nx,Ny,v,vv,produitvectoriel,statinfo,me,i1,iN)
				
				
				print*,'produitvectoriel',produitvectoriel
				call MPI_FINALIZE
			
		
			end program
		




		