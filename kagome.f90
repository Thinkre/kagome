! Note:
! variables which begin with letter of 'i' are used for cycle, such as ip, ib, il, iv, it and so on
! variables which begin with letter of 'n' are usually integer type number for correponding variable, such as nn,nd,no1,nt,nequ,nsam,nbin
! For globlal variables, we note only in 'var' module;
! for other variables usually used but not globlal:
! ip--address of 12-leg vertex configration
! v12(0:11)--12-leg state

module readin
    ! adjustable parameters which are read in 'input' file
    implicit none
    character(len=128) fname
    ! fname records the name of data file, which includes the executive time.

    integer::nl,nmax
    ! nl--size of system;
    ! nmax--maximum number in one site;

    integer::nequ,nsam
    ! nequ,nsam--number of steps of equlirium, sampling

    real(8)::v,t,temp,beta
end module readin


module loop_mod
    integer,allocatable::link(:),frst(:),last(:)
    ! define the array varialble in subroutine 'loop'. Acatually, these three array are only used in 'loop' subroutine.
end module loop_mod


module var
    use ran_mod 
    ! ran_mod module is external archive file, i.e. libran.a, including two type random function, i.e. ran() and ran(min,max)
    ! ran() generates real random number from 0.0 to 1.0 ,without including 0 and 1;
    ! ran(min,max) generates integer number from min to max, including min and max. 
    ! execute command-- -I$LIB -L$LIB -lran
    use swap_mod
    ! external file, i.e. libswap.a, swaping the same generic type value of two variable
    use readin
    use loop_mod

    integer,parameter::mmax=2000000

    ! integer type variable
    integer::nn,nb
    ! nn--total numbers of particles;
    ! nb--numbers of bonds;
    integer::nt,nd,no1,no2,no3
    !     nd--numbers of diagonal configration type;
    !  2*no1--number of one off-diagonal bonds vertex type;
    !  4*no2--number of two off-diagonal bonds vertex type;
    ! 16*no3--number of three off-diagonal bonds vertex type;
    ! nt=nd+2*no1+4*no2+16*no3
    integer::l,nh,nloop,ig,iga,nha
    ! l--cut-off length;
    ! nh--numbers of updated bonds;
    ! nloop--numbers of loop;
    ! ig--numbers of effective steps in nloop-cycle;
    ! nh--vairiable for adjusting cut-off length;
    ! ig,nh,iga,nha--vairiable for adjusting numbers of loop update steps

    integer::ne(0:11,6)

    ! float type varible
    real(8)::cadd   
    ! cadd--promises diagonal H positive 

    ! new data structure type
    type set
        integer::v(0:11)                                    
        real(8)::w
        real(8)::pb(0:1,0:11,0:11)
        integer::ip(0:1,0:11,0:11)  
    end type set
    ! set--variable set related type configration
    ! v--12-leg state
    ! w--weight
    ! pb(ir,nin,nout)--decrease(ir=0) or increase (ir=1) by one from leg-nin in and leg-nout out
    ! ip(ir,nin,nout)--new after decrease or increase by one from leg-nin in and leg-nout out
    type ops                        
        integer::b                  
        integer::ip                 
    end type ops
    ! ops--variable set for H(a,b)
    ! b--sequence of {a,b};b=even if diagonal bond,or old if off-diagonal bond, b=b/2
    ! ip--address for finding vertex configration, i.e. tp(ip)
    
    ! dynamic integer-type array
    integer,allocatable::lat(:),bond(:,:),ntime(:)
    ! lat(nn)--state on each site;
    ! bond(nb,0:5)--get hexgen; 
    ! ntime(it)--it=0,nh-1, records the effective 'il' in order, i.e. str(il)/=0
    type(set),allocatable::tp(:)        
    ! tp(ip)--configuration type acted by H(a,b)
    type(ops),allocatable::h(:) 
    ! h(il)--H(a,b) sequence


    ! measurement variable
    real(8)::ro,rosq,en,compress,super,qea
    real(8)::ro_err,en_err,compress_err,super_err,qea_err
    ! en--energy;
    ! ro--density;
    ! super--superfluid density;
    ! compress--compression ratio;
    ! err--corresponding tolerance
contains
    function fvd(ip)
        ! convert ip into nmax+1-number system(like binary,decimal), when ip<nd     
        integer::ip,ip1,fvd(0:11)
        if(ip>nd)stop "ip exceed the limits"
        ip1=ip-1
        do i=0,4
            fvd(i)=ip1/(nmax+1)**(5-i)
            ip1=mod(ip1,(nmax+1)**(5-i))
        end do
        fvd(5)=ip1
        fvd(6:11)=fvd(0:5)
    end function fvd

    function fv(ip)
        ! get 12 vertex-legs from ip                           
        implicit none   
        integer::ip,ip1,np,np1,ni,nj,nk,i,j,n1,n2,ino,four(4),v12(0:11),fv(0:11)

        if(ip==0)stop "ip==0"
        if(ip<=nd)then
            fv=fvd(ip)
        elseif(ip<=nd+2*no1)then
            ! ip=np*nmax**2*(nmax+1)**4+(ni-1)*(nmax+1)**4+nj
            ip1=ip-nd-1
            ino=ip1/no1                 
            ip1=mod(ip1,no1)                    
            ! np is the place where off-bond interacts      
            np=ip1/((nmax+1)**4*nmax**2)

            ip1=mod(ip1,(nmax+1)**4*nmax**2)
            ! ni is the sequence of two connected sites.
            ni=ip1/(nmax+1)**4+1                

            ! nj is the sequence of other four sites.
            nj=mod(ip1,(nmax+1)**4)+1
            ! get state of two sites, ni=n1*nmax+n2 when add=1

            n1=(ni-1)/nmax
            n2=mod(ni-1,nmax)+1
            ! get state of other four sites without off-bond interaction
            v12=fvd(nj)

            j=2
            i=0 
            do while(i<6)
                if(i==np)then
                    fv(i)=n1
                    fv(i+6)=n1+1
                    fv(mod(i+1,6))=n2
                    fv(mod(i+1,6)+6)=n2-1
                    if(ino==1)then
                        call swap(fv(i),fv(mod(i+1,6)))
                        call swap(fv(i+6),fv(mod(i+1,6)+6))
                    endif
                    i=i+1
                elseif(i/=mod(np+1,6))then
                    fv(i)=v12(j)
                    fv(i+6)=v12(j)
                    j=j+1
                endif
                i=i+1
            end do

        elseif(ip<=nd+2*no1+4*no2)then
            ! two bonds
            ip1=ip-nd-2*no1-1
            ino=ip1/no2

            ip1=mod(ip1,no2)
            np=ip1/(nmax**4*(nmax+1)**2)
            
            ip1=mod(ip1,nmax**4*(nmax+1)**2)
            ni=ip1/nmax**4

            ip1=mod(ip1,nmax**4)
            nj=ip1/nmax**2

            nk=mod(ip1,nmax**2)

            v12(0)=ni/(nmax+1)
            v12(6)=ni/(nmax+1)
            v12(1)=mod(ni,nmax+1)
            v12(7)=mod(ni,nmax+1)

            v12(2)=nj/nmax
            v12(3)=mod(nj,nmax)+1
            v12(8)=v12(2)+1
            v12(9)=v12(3)-1
            if(ino/2==1)then
                call swap(v12(2),v12(3))
                call swap(v12(8),v12(9))
            endif

            v12(4)=nk/nmax
            v12(5)=mod(nk,nmax)+1
            v12(10)=v12(4)+1
            v12(11)=v12(5)-1
            if(mod(ino,2)==1)then
                call swap(v12(4),v12(5))
                call swap(v12(10),v12(11))
            endif           

            if(np<6)then
                do i=0,5
                    fv(mod(np+i,6))=v12(i)
                    fv(mod(np+i,6)+6)=v12(i+6)
                enddo
            else
                np1=mod(np,6)
                fv(np1)=v12(0)
                fv(np1+3)=v12(1)
                fv(np1+6)=v12(6)
                fv(np1+9)=v12(7)
                do i=1,2
                    fv(np1+i)=v12(i+1)
                    fv(np1+i+6)=v12(i+7)
                    fv(mod(np1+i+3,6))=v12(i+3)
                    fv(mod(np1+i+3,6)+6)=v12(i+9)
                end do
            endif
        elseif(ip<=nd+2*no1+4*no2+16*no3)then
            ! three bonds
            ip1=ip-(nd+2*no1+4*no2)-1
            ino=ip1/no3
            ip1=mod(ip1,no3)
            ni=ip1/nmax**4

            ip1=mod(ip1,nmax**4)
            nj=ip1/nmax**2

            nk=mod(ip1,nmax**2)

            v12(0)=ni/nmax
            v12(1)=mod(ni,nmax)+1
            v12(6)=v12(0)+1
            v12(7)=v12(1)-1

            v12(2)=nj/nmax
            v12(3)=mod(nj,nmax)+1
            v12(8)=v12(2)+1
            v12(9)=v12(3)-1

            v12(4)=nj/nmax
            v12(5)=mod(nj,nmax)+1
            v12(10)=v12(4)+1
            v12(11)=v12(5)-1

            do i=1,4
                four(i)=ino/2**(4-i)
                ino=mod(ino,2**(4-i))
            enddo

            do i=0,2
                if(four(i+2)==1)then
                    call swap(v12(2*i),v12(2*i+1))
                    call swap(v12(2*i+6),v12(2*i+7))
                endif
            end do
            if(four(1)==1)then
                fv(1:5)=v12(0:4)
                fv(0)=v12(5)
                fv(7:11)=v12(6:10)
                fv(6)=v12(11)
            else
                fv=v12
            endif
        endif
    end function fv                
end module var


program main
    use var;implicit none
    integer::i

    read(*,*)fname,t,v,temp,nl,nmax,nequ,nsam

    call config

    call weight

    ! for equlibrium
    do i=1,nequ
        call diagonal
        call loop
        call adjustl
        call adjustnloop
    end do

    ! sampling
    do i=1,nsam
        call diagonal
        call loop
        call measure
    end do
end program


subroutine config
    ! set up configration
    use var;implicit none
    integer,allocatable::site(:,:)
    integer::disp(6,2),i,j,k

    ! integer-type variabel asignment
    nn=3*nl**2
    nb=nl**2
    nd=(nmax+1)**6
    no1=6*nmax**2*(nmax+1)**4
    no2=9*nmax**4*(nmax+1)**2
    no3=nmax**6
    nt=nd+2*no1+4*no2+16*no3

    ! float type
    beta=1.d0/temp

    ! allocate size of dynamic array
    allocate(site(0:2*nl-1,0:2*nl-1))
    ! site(:,:) is local variable
    allocate(lat(nn),bond(nb,0:5))
    allocate(tp(0:nt))
    allocate(frst(nn),last(nn))
    allocate(h(mmax),ntime(0:mmax),link(0:12*mmax-1))

    ! get random state on each site
    call initran    
    lat=[(ran(0,nmax),i=1,nn)]      

    ! map coordinate site(x,y) to each site on kagome lattice
    k=1;site=0      
    do j=0,2*nl-1
        do i=0,2*nl-1
            if(mod(i*j,2)==1) cycle         
            site(i,j)=k
            k=k+1
        end do
    end do
    ! build bond(nb,1:2)    
    disp=reshape([0,1,1,0,-1,-1,-1,-1,0,1,1,0],[6,2])
    bond=reshape([(((site(mod(disp(k,1)+i,2*nl),mod(disp(k,2)+j,2*nl)),k=1,6),i=1,2*nl-1,2),j=1,2*nl-1,2)],[nb,6])

    ! neighbor
    do i=0,11
        ne(i,2)=mod(i,6)
        ne(i,1)=5-mod(6-ne(i,2),6)
        ne(i,3)=mod(ne(i,2)+1,6)
        ne(i,4:6)=ne(i,1:3)+6
    end do

    ! initiates varible
    nh=0
    l=20

    h=ops(0,0)
    tp=set(-1,0.d0,0.d0,0)
    nloop=10
    nha=0
    iga=0   
end subroutine config


subroutine weight
    use var;implicit none
    integer::ip,np,i,inp,fip,ir,nin,six(0:5),v12(0:11)

    ! constant added in dialgonal term of H for promising values positive
    cadd=1.1*v*(6*nmax-3.d0)**2

    ! assign values to tp(1:nt)
    do ip=1,nt
        ! get 12-leg state
        v12=fv(ip)
        tp(ip).v=v12
        six=v12(6:11)-v12(0:5)

        ! calculate weight at ip, i.e. tp(ip).w
        if(ip<=nd)then  
            tp(ip).w=cadd-v*(sum(v12(0:5))-3.d0)**2     
        else
            i=0
            if(ip>nt-8*no3) i=i+1
            ! this is for distinguishing bond-type
            tp(ip).w=0.d0
            do while(i<6)
                if(six(i)==1 .and. six(mod(i+1,6))==-1)then
                    tp(ip).w=tp(ip).w+t*sqrt(dfloat((tp(ip).v(i)+1)*tp(ip).v(mod(i+1,6))))
                    i=i+2
                elseif(six(i)==-1 .and. six(mod(i+1,6))==1)then
                    tp(ip).w=tp(ip).w+t*sqrt(dfloat((tp(ip).v(mod(i+1,6))+1)*tp(ip).v(i)))
                    i=i+2
                else
                    i=i+1
                endif
            end do
        endif
    end do

    ! calculate probability
    do ip=1,nt
        do nin=0,11                      
            do ir=0,1
                ! calculate prbability from nin in and iout out
                call calcpb(ip,ir,nin)
            end do
        end do
    end do
end subroutine weight


function fip(v12,ch)
    ! get ip via v12, ch is for distinguishing bond-type;
    ! when v12 represent 3 off-diagonal bond, we must need another argument to judge where bonds actually acts on;
    ! the reason is that bond my act on 0-1-6-7, 2-3-8-9, 4-5-10-11 as one case and 1-2-7-8, 3-4-9-10, 5-0-11-6 as other case;
    ! So, we define ch=1 for 0-1-6-7, 2-3-8-9, 4-5-10-11 case and ch=2 for 1-2-7-8, 3-4-9-10, 5-0-11-6 case.
    use var;implicit none
    integer::fip,v12(0:11),ch
    integer::i,j,k,icase,n7
    integer::six(0:5),four(0:3)

    fip=0
    do i=0,11
        if(v12(i)<0 .or. v12(i)>nmax)return
    end do

    six=v12(6:11)-v12(0:5)
    do i=0,5
        if(six(i)==0 .and. six(mod(i+1,6))/=0 .and. six(mod(i+2,6))==0)return
    end do
    icase=sum(abs(six))
    if(mod(icase,2)==1)return
    if(icase==6)then
        do i=0,5
            if(six(i)==six(mod(i+1,6)) .and. six(i)==six(mod(i+2,6)))return
        end do
    elseif(icase==4)then
        do i=0,5
            if(six(i)==six(mod(i+1,6)) .and. six(i)==six(mod(i+2,6)))return
        end do
        i=0
        do while(i<6)
            if(six(i)==0 .and. six(mod(i+1,6))==0)then
                four=[(six(mod(i+2+j,6)),j=0,3)]
                exit
            elseif(six(i)==0 .and. six(mod(i+3,6))==0)then
                four(0)=six(mod(i+1,6))
                four(1)=six(mod(i+2,6))
                four(2)=six(mod(i+3,6))
                four(3)=six(mod(i+5,6))
                exit
            endif
            i=i+1
        end do
        if(four(0)==four(1) .or. four(2)==four(3))return
    elseif(icase==2)then
        do i=0,5
            if(six(i)/=0 .and. six(i)==six(mod(i+1,6)))return
        end do
    endif

    if(icase==0)then
        fip=1
        do i=0,5
            fip=fip+v12(i)*(nmax+1)**(5-i)
        end do
    elseif(icase==2)then
        fip=nd+1
        j=4
        i=0
        do while(i<6)
            if(six(i)==0)then
                j=j-1
                fip=fip+v12(i)*((nmax+1)**j)
            elseif(six(i)==1 .and. six(mod(i+1,6))==-1)then
                fip=fip+(v12(i)*nmax+v12(mod(i+1,6))-1)*(nmax+1)**4+i*nmax**2*(nmax+1)**4
                i=i+1
            elseif(six(i)==-1 .and. six(mod(i+1,6))==1)then
                fip=fip+no1+(v12(mod(i+1,6))*nmax+v12(i)-1)*(nmax+1)**4+i*nmax**2*(nmax+1)**4
                i=i+1               
            endif
            i=i+1
        end do
    elseif(icase==4)then
        fip=nd+2*no1+1
        do i=0,5
            if(six(i)==0)then
                if(six(mod(i+1,6))==0)then
                    fip=fip+i*nmax**4*(nmax+1)**2+(v12(i)*(nmax+1)+v12(mod(i+1,6)))*nmax**4
                    if(six(mod(i+2,6))==1 .and. six(mod(i+3,6))==-1)then
                        fip=fip+(v12(mod(i+2,6))*nmax+v12(mod(i+3,6))-1)*nmax**2
                    elseif(six(mod(i+2,6))==-1 .and. six(mod(i+3,6))==1)then
                        fip=fip+(v12(mod(i+2,6))+v12(mod(i+3,6))*nmax-1)*nmax**2+2*no2
                    endif
                    if(six(mod(i+4,6))==1 .and. six(mod(i+5,6))==-1)then
                        fip=fip+(v12(mod(i+4,6))*nmax+v12(mod(i+5,6))-1)
                    elseif(six(mod(i+4,6))==-1 .and. six(mod(i+5,6))==1)then
                        fip=fip+(v12(mod(i+4,6))+v12(mod(i+5,6))*nmax-1)+no2
                    endif
                    exit
                elseif(i<3 .and. six(mod(i+3,6))==0)then
                    fip=fip+(i+6)*nmax**4*(nmax+1)**2+(v12(i)*(nmax+1)+v12(mod(i+3,6)))*nmax**4
                    if(six(mod(i+1,6))==1 .and. six(mod(i+2,6))==-1)then
                        fip=fip+(v12(mod(i+1,6))*nmax+v12(mod(i+2,6))-1)*nmax**2
                    elseif(six(mod(i+1,6))==-1 .and. six(mod(i+2,6))==1)then
                        fip=fip+(v12(mod(i+1,6))+v12(mod(i+2,6))*nmax-1)*nmax**2+2*no2
                    endif
                    if(six(mod(i+4,6))==1 .and. six(mod(i+5,6))==-1)then
                        fip=fip+(v12(mod(i+4,6))*nmax+v12(mod(i+5,6))-1)
                    elseif(six(mod(i+4,6))==-1 .and. six(mod(i+5,6))==1)then
                        fip=fip+(v12(mod(i+4,6))+v12(mod(i+5,6))*nmax-1)+no2
                    endif
                    exit
                endif
            endif
        end do
    elseif(icase==6)then
        if(ch==1)then
            fip=nd+2*no1+4*no2+1
            i=0
        elseif(ch==2)then
            fip=nd+2*no1+4*no2+8*no3+1
            i=1
        endif
        j=2
        k=2
        do while(i<6)
            if(six(i)==1 .and. six(mod(i+1,6))==-1)then
                fip=fip+nmax**j*(v12(i)*nmax+v12(mod(i+1,6))-1)
            else
                fip=fip+nmax**j*(v12(mod(i+1,6))*nmax+v12(i)-1)+2**(2-i/2)
            endif
            i=i+2
            j=j-1
        end do
    endif

    if(fip>nt)stop "fip>nt in fun fip"
end function fip


subroutine calcpb(ip,ir,nin)
    use var;implicit none
    integer::ip,ir,nin
    integer::np,fip,i,j,cnt
    integer::n6(6),ip6(6),nc(6),v12(0:11),six(0:5),loc(1)
    real(8)::w6(6),a(6,6),wt

    if(tp(ip).v(nin)==nmax*ir)return

    n6=ne(nin,1:6)

    do i=1,6
        if(n6(i)==nin)then
            ip6(i)=ip
        else
            v12=tp(ip).v

            v12(nin)=v12(nin)+2*ir-1
            v12(n6(i))=v12(n6(i))+(2*mod(nin/6+n6(i)/6,2)-1)*(2*ir-1)

            ip6(i)=fip(v12,mod(mod(ip-(nd+2*no1+1),no2)/(nmax**4*(nmax+1)**2),2)+1)
        endif
    end do

    if(ip>nt-16*no3)then
        ip6(1+2*mod((ip-(nt-16*no3)-1)/8+nin,2))=0
        ip6(4+2*mod((ip-(nt-16*no3)-1)/8+nin,2))=0
    endif

    w6=[(tp(ip6(i)).w,i=1,6)]

    if(sum(w6)<=0.d0)stop "calcpb: weight<0"
    a=0.d0
    do 
        loc=minloc(w6,w6>=0.d0)
        i=loc(1)
        if(i>6 .or. i<1)exit
        cnt=6
        do j=1,6
            if(w6(j)<0.d0 .or. abs(mod(i-1,3)-mod(j-1,3))==2 .or. i==j)cnt=cnt-1
        end do
        if(cnt==0)then
            a(i,i)=w6(i)
            exit
        endif
        do j=1,6
            if(w6(j)<0.d0 .or. abs(mod(i-1,3)-mod(j-1,3))==2 .or. i==j)cycle
            a(i,j)=w6(i)/cnt
            a(j,i)=a(i,j)
            w6(j)=w6(j)-a(j,i)
        end do
        w6(i)=-1.d0
    end do
    
    nc=[(mod(nin/6+ir,2),i=1,3),(1-mod(nin/6+ir,2),i=1,3)]

    ! calc pb and record new ip
    do i=1,6
        do j=1,6
            tp(ip6(i)).pb(nc(i),n6(i),n6(j))=a(i,j)/tp(ip6(i)).w
            tp(ip6(i)).ip(nc(i),n6(i),n6(j))=ip6(j)
        end do
    end do
end subroutine calcpb


subroutine diagonal
    use var;implicit none
    integer::il,ib,np,i,ip,fip,v12(0:11)

    do il=1,l
        if(h(il).b==0)then
            ! (0,0)-->(1,b), try to update identity bond to dialgonal bond.
            ib=ran(1,nb)
            v12(0:5)=[(lat(bond(ib,i)),i=0,5)]
            v12(6:11)=v12(0:5)
            if((l-nh)<0) stop 'l-nh<0'
            ip=fip(v12,0)
            if(ran()<(beta*nb*(tp(ip).w)/(l-nh)))then 
                ! h(il) is not set up and must be got via nvec
                h(il)=ops(2*ib,ip)
                nh=nh+1
            endif
        elseif(mod(h(il).b,2)==0)then
            ! (1,b)-->(0,0), try to update dialgonal bond to identity bond.
            if(ran()<((l-nh+1)/beta/nb/tp(h(il).ip).w))then
                h(il)=ops(0,0)
                nh=nh-1
            endif
        else
            ! for off-dialgonal term, lattice state needs to be updated and h(il) was set up after loop update  
            lat(bond(h(il).b/2,0:5))=tp(h(il).ip).v(6:11)
        endif
    end do
end subroutine diagonal


subroutine loop
    use var;implicit none
    integer::i,j,k,iv,v0,ip,ir,il,nin,nout,iloop
    
    call linkvrtx     
    ig=0
    if(nh==0) return

    do iloop=1,nloop
        !call link1test

        ! get loop head v0
        call ran_v0(v0,nin,il,ip,ir)
        ! go to next vertex 
        iv=v0
        do
            call whichout(nout)

            ! next-in vertex
            iv=link(12*(iv/12)+nout)

            ! records the numbers of bond-type change
            if(nin/=nout)ig=ig+1            

            ! test whether the loop was closed, exit cycle if yes,otherwise continue cycle 
            if(v0==iv.or.v0==link(iv))exit

            call whichin
        end do
    end do

    ! change state of lattice
    do i=1,nn
        iv=frst(i)
        if(iv/=-1)then
            lat(i)=tp(h(ntime(iv/12)).ip).v(mod(iv,12))
        else
            lat(i)=ran(0,nmax)
        endif
    end do
contains
    subroutine linkvrtx
        integer::it
        link=-1
        frst=-1
        last=-1
        iv=0
        it=0            
        do il=1,l
            if(h(il).b==0)cycle
            ntime(it)=il
            do i=0,5
                j=bond(h(il).b/2,i)
                if(last(j)==-1)then
                    frst(j)=iv+i
                else
                    link(last(j))=iv+i
                    link(iv+i)=last(j)
                endif
                last(j)=iv+i+6
            end do
            iv=iv+12
            it=it+1
        end do
        if(nh/=it)stop "nh/=numbers of bonds, wrong"
        ! due to periodical boundary 
        do i=1,nn
            if(frst(i)==-1)cycle
            link(last(i))=frst(i)
            link(frst(i))=last(i)
        end do
    end subroutine linkvrtx

    subroutine ran_v0(a_v0,a_in,a_il,a_ip,a_ir)
        integer::a_v0,a_in,a_il,a_ip,a_ir
        do
            a_ir=ran(0,1)
            a_v0=ran(0,12*nh-1)
            a_in=mod(a_v0,12)
            a_il=ntime(a_v0/12)
            a_ip=h(a_il).ip
            if(tp(a_ip).v(a_in)/=nmax*a_ir)return
        end do
    end subroutine ran_v0

    subroutine whichout(a_out)
        ! choose leg out
        integer::a_out
        real(8)::pr,pr1,pr2

        pr=ran()
        pr1=0.d0
        do i=1,6
            pr2=pr1
            a_out=ne(nin,i)
            pr1=pr1+tp(ip).pb(ir,nin,a_out)
            if((pr.ge.pr2).and.(pr.lt.pr1))exit
        end do

        ! update H(a,b) at il, i.e. h(il).ip
        h(il).ip=tp(ip).ip(ir,nin,a_out)
        if(h(il).ip<=0)stop "whichout: h(il).ip<=0"
        ! update H(a,b) at il, i.e. h(il).b
        if(h(il).ip>nd)then
            h(il).b=2*(h(il).b/2)+1
        else
            h(il).b=2*(h(il).b/2)
        endif
    end subroutine whichout

    subroutine whichin
        ir=mod(nin/6+nout/6+ir+1,2)
        nin=mod(iv,12)
        il=ntime(iv/12)
        ip=h(il).ip     
    end subroutine whichin
end subroutine loop


subroutine adjustl
    use var;implicit none
    integer::lnew,i,it,nh0,il,inh,iop
    type(ops),allocatable::h1(:)
    integer,allocatable::insert(:)

    allocate(insert(nh))
    allocate(h1(nh))

    nh0=nh/4
    lnew=nh+nh0
    if(lnew<l)then
        return
    elseif(lnew>mmax)then
        stop "L exceeds the maximum"
    endif

    ! h1 stores non-identity bond
    do it=0,nh-1
        h1(it+1)=h(ntime(it))
    end do

    ! random insert identity bond, insert(iop) stands for number of identity bonds
    insert=0
    do i=1,nh0
        iop=ran(1,nh)
        insert(iop)=insert(iop)+1
    end do

    inh=1
    il=1
    do iop=1,nh
        if(insert(iop)/=0)then
            h(il:il+insert(iop)-1)=ops(0,0)
            il=il+insert(iop)
        endif
        h(il)=h1(inh)
        il=il+1
        inh=inh+1
    end do
    if(inh<=nh)stop "Oops, inh<=nh in adjustl"
    l=lnew
end subroutine adjustl


subroutine adjustnloop
    use var;implicit none
    integer::i
    real(8)::ra

    nha=nha+nh
    iga=iga+ig      
    if(mod(i,nequ/10)==0)then
        ra=iga/dfloat(nha)/4
        if(ra<0.8)then
            nloop=nloop*5/3
        elseif(ra>1.5)then
            nloop=max(nloop*3/4,1)
        endif
        if(nloop<1)nloop=1
        nha=0;iga=0
    endif       
end subroutine adjustnloop


subroutine measure
    use var;implicit none
    real(8)::windx,windy
    integer::b,i,ip

    ! energy
    en=cadd/3.d0-dfloat(nh)/dfloat(nn)/beta

    ! density
    ro=dfloat(sum(lat))/dfloat(nn)

    write(*,"(2f16.7)") en,ro
end subroutine measure



! test code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! SUB: calcpb(ir,nin,nout)
! -------------------------------------------------------------------------------------------------
!     do i=1,6
!         do j=1,6
!             if(tp(ip6(i)).w<=0.d0 .or. a(i,j)<=0.d0)cycle
!             tp(ip6(i)).pb(nc(i),n6(i),n6(j))=a(i,j)/tp(ip6(i)).w
!             tp(ip6(i)).ip(nc(i),n6(i),n6(j))=ip6(j)
!             write(7,'(5I7,1f10.4)')ip6(i),nc(i),n6(i),n6(j),ip6(j),a(i,j)/tp(ip6(i)).w
!         end do
!     end do
! -------------------------------------------------------------------------------------------------

! SUB: loop
! -------------------------------------------------------------------------------------------------
! do iv=0,12*nh-1
!   if(link(link(iv))/=iv)then
!       print*,iv,link(iv),link(link(iv))
!       stop "linkvrtx: link error"
!   else
!       if(tp(h(ntime(iv/12)).ip).v(mod(iv,12))/=tp(h(ntime(link(iv)/12)).ip).v(mod(link(iv),12)))then
!           print*,'----------------------------------------'
!           write(*,"(2I6)") iv,mod(iv,12)
!           write(*,"(6I6)") tp(h(ntime(iv/12)).ip).v(0:5)
!           write(*,"(6I6)") tp(h(ntime(iv/12)).ip).v(6:11)
!           print*," "
!           write(*,"(2I6)") link(iv),mod(link(iv),12)
!           write(*,"(6I6)") tp(h(ntime(link(iv)/12)).ip).v(0:5)
!           write(*,"(6I6)") tp(h(ntime(link(iv)/12)).ip).v(6:11)
!           open(8,file='link.dat')
!           do i=0,12*nh-1
!               write(8,'(2I6)')i,link(i)
!           end do
!           close(8)
!           print*," "
!           stop "linkvrtx---------------------------------"
!       endif
!   endif
! end do

! if(abs(mod(nin,6)-mod(nout,6))==2)stop "loop:nin,nout"
! print*,"iv,ip,ir,nin,nout,link-out"
! write(*,"(5I6)") iv,ip,ir,nin,nout
! write(*,"(6I4)") tp(ip).v(0:5)
! write(*,"(6I4)") tp(ip).v(6:11)
! write(*,"(6I4)") tp(h(il).ip).v(0:5)
! write(*,"(6I4)") tp(h(il).ip).v(6:11)
! -------------------------------------------------------------------------------------------------

