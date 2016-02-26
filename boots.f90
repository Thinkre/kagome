module readin
	implicit none
	character(len=128) fname
	integer::nl,nmax
	! nl--size of system;
	! nmax--maximum number in one site; 
	integer::nequ,nsam
	! nequ,nsam,nbin--numbers of steps of equlirium, sampling
	real(8)::v,t,temp,beta
end module readin

program boots
	use readin
	use bootstrap_mod
	implicit none
	character(len=128)::fname1,dir
	integer::nboot
	real(8),allocatable::en(:),ro(:),super(:)
	real(8)::en_ave,en_err,ro_ave,ro_err,super_ave,super_err
	integer::i,j,k
	logical alive

	nboot=100
	read(*,*)fname,t,v,temp,nl,nmax,nequ,nsam
	beta=1/temp
	allocate(en(nsam),ro(nsam),super(nsam))
	dir=trim(adjustl(fname(scan(fname,'_')+1:len(fname))))
	inquire(directory=trim(adjustl(dir)),exist=alive)
	if(.not.alive)stop "directory dosen't exist"
	inquire(file=trim(adjustl(fname1)),exist=alive)
	if(.not.alive) stop "data file dosen't exist"
    open(8,file='./'//trim(adjustl(dir))//'/'//trim(adjustl(fname))//'.dat')
	do j=1,nsam
		read(8,*)en(j),ro(j),super(j)
	enddo
	close(8)
	call bootstrap(en,nboot,en_ave,en_err)
	call bootstrap(ro,nboot,ro_ave,ro_err)
	call bootstrap(super,nboot,super_ave,super_err)
	write(*,"(10f17.8)")t,v,temp,beta,en_ave,en_err,ro_ave,ro_err,super_ave,super_err
end program boots
