module in_data
  double precision, allocatable :: xyz(:,:) 
  double precision ::  dl, dx, dy, xl(2), yl(2), x0, y0, th
  integer :: n0, nx, ny, ndiv
  integer, allocatable :: ain(:,:), num(:,:)
end module

module out_data
  double precision, allocatable :: dat(:,:) 
  integer :: n1
  integer, allocatable :: aout(:,:)
end module


program maketopo
use in_data
use out_data
implicit none
integer :: i, j, k, l, n
double precision :: dum, dx0, dy0

!write(*,*) 'ndiv ?'
!read(*,*) ndiv
ndiv=20 


open(11,file='./input/area4.txt')
read(11,*) dx, dy
read(11,*) x0, y0
read(11,*) dl, th
close(11)

call read_data
dx0 = xl(2) - xl(1)
dy0 = yl(2) - yl(1)

nx = nint(dx)/nint(dl)+1 ; ny = nint(dy)/nint(dl)+1
write(*,'(a,f10.2,a)') 'dx = ', dx, '[m]'
write(*,'(a,f10.2,a)') 'dy = ', dy, '[m]'
write(*,*) 'nx & ny :', nx, ny

n1 = 0 
allocate(dat(nx*ny,3), aout(nx*ny, 2))
dat(:,:) = 0.d0
do j = 1, ny
  do i = 1, nx
    n1 = n1 + 1
    dat(n1,1) = x0 + dble(i-1)*dl*dcos(th) - dble(j-1)*dl*dsin(th)
    dat(n1,2) = y0 + dble(i-1)*dl*dsin(th) + dble(j-1)*dl*dcos(th)
    k = int((dat(n1,1)-xl(1))/dx0*ndiv) + 1
    l = int((dat(n1,2)-yl(1))/dy0*ndiv) + 1
    if(k>ndiv) k = ndiv
    if(l>ndiv) l = ndiv
    aout(n1,1) = k ; aout(n1,2) = l
  end do
end do
write(*,*) 'Number of out_data : ', n1


n = 0
do j = 1,ndiv
  do i = 1,ndiv
    n = n + num(i,j)
    write(*,*) i, j, num(i,j)
  end do
end do
write(*,*) n

do i = 1, n1
  call line_itp(i)
  if(mod(i,10000)==0) write(*,*) i
end do

open(12,file='./input/coordinate.txt',status='replace')
k = 0
do j= 1, ny
  do i = 1, nx
    k = k + 1
    write(12,'(3f10.2)') x0 + dble(i-1)*dl, y0 + dble(j-1)*dl, dat(k,3)
  end do
end do
close(12)

end program

subroutine read_data
use in_data
implicit none
integer :: i, j, nn, nda, k, l
integer, allocatable :: n(:)
double precision :: dum, xx(4), yy(4)
double precision, allocatable :: xyz0(:,:)
character(30) :: cdum
character(30), allocatable :: fname(:)

open(11,file='./input/fname_ndata.txt',status='old')
read(11,*) cdum, nda
allocate(n(nda), fname(nda))

n0 = 0
do i = 1,nda
  read(11,*) fname(i), n(i)
  n0 = n0 + n(i)
end do
close(11)

allocate(xyz0(n0,3))

xx(1) = x0 
xx(2) = x0 + dx*dcos(th) 
xx(3) = x0 - dy*dsin(th)
xx(4) = x0 + dx*dcos(th) - dy*dsin(th)
yy(1) = y0 
yy(2) = y0 + dx*dsin(th)
yy(3) = y0 + dy*dcos(th)
yy(4) = y0 + dx*dsin(th) + dy*dcos(th)

xl(1) = minval(xx) - 4.d0*dl
xl(2) = maxval(xx) + 4.d0*dl
yl(1) = minval(yy) - 4.d0*dl
yl(2) = maxval(yy) + 4.d0*dl

nn = 0
n0 = 0
do j = 1, nda
  open(11,file=fname(j),status='old')
  do i = 1,n(j)
    nn = nn + 1
    read(11,*) dum, xyz0(nn,2), xyz0(nn,1), xyz0(nn,3)
    if(xyz0(nn,1) >= xl(1) .and. xyz0(nn,1) <= xl(2)) then
      if(xyz0(nn,2) >= yl(1) .and. xyz0(nn,2) <= yl(2)) then
        n0 = n0 + 1
      endif
    endif
  end do
  write(*,*) j, fname(j), nn 
  close(11)
end do

write(*,*) "nn: ", nn
write(*,*) "n0: ", n0
allocate(xyz(n0,3), ain(int(n0/ndiv),ndiv*ndiv), num(ndiv,ndiv))

j = 0 ; num(:,:) =0 ; ain(:,:) = 0
do i = 1, nn
  if(xyz0(i,1) >= xl(1) .and. xyz0(i,1) <= xl(2)) then
    if(xyz0(i,2) >= yl(1) .and. xyz0(i,2) <= yl(2)) then
      j = j + 1
      xyz(j,:) = xyz0(i,:)
      k = int((xyz(j,1)-xl(1))/(xl(2)-xl(1))*ndiv) + 1
      l = int((xyz(j,2)-yl(1))/(yl(2)-yl(1))*ndiv) + 1
      num(k,l) = num(k,l) + 1
      ain(num(k,l),(l-1)*ndiv+k) = j
    endif
  endif
end do
deallocate(xyz0)

end subroutine

subroutine line_itp(n)
use in_data 
use out_data
implicit none
integer,intent(in) :: n
integer :: i, j, k, l, ii, pp(4), ij(9,2)
double precision :: d0, dd(4), xx, yy, zz(2)

k=0
do j = -1,1
  do i = -1,1
    k = k + 1
    ij(k,1) = aout(n,1)+i
    ij(k,2) = aout(n,2)+j
  end do
end do

dd(:) = 1000.d0
pp(:) = 0
do k = 1,9
  i = ij(k,1) ; j = ij(k,2)
  if(mod(i,ndiv+1)==0 .or. mod(j,ndiv+1)==0) cycle
  ii = (j-1) * ndiv + i
  do l = 1,num(i,j)
    xx = xyz(ain(l,ii),1) - dat(n,1) ; yy = xyz(ain(l,ii),2) - dat(n,2)
    d0 = dsqrt(xx**2.d0 + yy**2.d0)
    !write(*,*) 'd0 = ', d0
    if(d0==0.d0) then
      dat(n,3) = xyz(ain(l,ii),3)
      exit
    elseif(d0>20.d0) then
      cycle
    else
      if(xx>=0.d0 .and. yy>=0.d0) then
        if(dd(1) > d0) then
          dd(1) = d0 ; pp(1) = ain(l,ii)
        endif
      elseif(xx<=0.d0 .and. yy>=0.d0) then
        if(dd(2) > d0) then
          dd(2) = d0 ; pp(2) = ain(l,ii)
        endif
      elseif(xx<=0.d0 .and. yy<=0.d0) then
        if(dd(3) > d0) then
          dd(3) = d0 ; pp(3) = ain(l,ii)
        endif
      elseif(xx>=0.d0 .and. yy<=0.d0) then
        if(dd(4) > d0) then
          dd(4) = d0 ; pp(4) = ain(l,ii)
        endif
      end if
    end if
  end do
 
  if(d0==0.d0) exit
end do


if(d0 > 0.d0) then
  xx = (xyz(pp(1),2)-xyz(pp(2),2))/(xyz(pp(1),1)-xyz(pp(2),1))* &
       (dat(n,1) - xyz(pp(2),1)) + xyz(pp(2),2)
  yy = (xyz(pp(4),2)-xyz(pp(3),2))/(xyz(pp(4),1)-xyz(pp(3),1))* &
       (dat(n,1) - xyz(pp(3),1)) + xyz(pp(3),2)
  zz(1) = (xyz(pp(1),3)-xyz(pp(2),3))/(xyz(pp(1),1)-xyz(pp(2),1))* &
          (dat(n,1) - xyz(pp(2),1)) + xyz(pp(2),3)
  zz(2) = (xyz(pp(4),3)-xyz(pp(3),3))/(xyz(pp(4),1)-xyz(pp(3),1))* &
          (dat(n,1) - xyz(pp(3),1)) + xyz(pp(3),3)
  dat(n,3) = (zz(1)-zz(2))/(xx-yy)*(dat(n,2)-yy) + zz(2)
end if

if(minval(pp)==0 .and. d0>0) then
  dat(n,3) = -9999
end if
end subroutine


