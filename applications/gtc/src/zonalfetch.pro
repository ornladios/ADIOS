pro zonalfetch
; read run parameters
; XY customize
; this program reads potential data and make zeta fft ( to k_zeta domain) to get powerspectrum in parallel
; for phi^2(k_zata) average over theta and a range in r center at ac
; write parallel+fnstep+.dat to dir ./data

; set up parameters:
;;*******************************

;; file reading control (which time step and how many to read)
ntimes=520 ;how many time steps or files to be read
nstart=5  ;start from which times step
ndelt=5   ;time interval between two ajacent time steps

;;;average r range from m1 to m2
;;need to set m1 and m2 and nmode  int the middle of the program
znmode=16;  # of parallel modes = one half of metza=32



wdir='./data/' ;writing dir for P+fnstep+.dat
dir='./'  ;reading dir for PHI files



frun=dir+'RUNdimen.ncd'
vnmpsi='flux-surface-number'
vnmtheta='poloidal-grids'
vnr='radial-grids'

;;*****************************


print, 'open: ', frun

; open netcdf file
ncid=ncdf_open(frun)

; get variable id;
vidmpsi=ncdf_varid(ncid,vnmpsi)
vidmtheta=ncdf_varid(ncid,vnmtheta)
vidr=ncdf_varid(ncid,vnr)

; read data
ncdf_varget,ncid,vidmpsi,mpsi
mpsi=mpsi-1
mtheta=indgen(mpsi+1)
r=fltarr(mpsi+1)

ncdf_varget,ncid,vidmtheta,mtheta
ncdf_varget,ncid,vidr,r
mtheta=mtheta-1

;grid system
igrid=mtheta
igrid(0)=0
mgrid=0

for i=1,mpsi do begin
    igrid(i)=igrid(i-1)+mtheta(i-1)+1
    mgrid=mgrid+mtheta(i-1)
endfor

mgrid=mgrid+mtheta(mpsi)
print,'mpsi=',mpsi, 'real grid point #=',mgrid
print, 'End of reading parameters'


; compute parallel spectrum
;;average r range from m1 to m2
;;************************
m1=floor(mpsi*1/4) ; consider an annulus with mpsi/4<r<mpsi*3/4
m2=floor(mpsi*3/4)

nmode=znmode ;  # of parallel modes = one half of metza=32

;;*************************



n1=1
n2=ntimes


zonarr=fltarr(m2-m1,ntimes-1)

for it=n1,n2 do begin
nstep=nstart+ndelt*(it-1)
fnstep=strcompress(string(nstep),/remove_all)
fname=dir+'PHI_'+fnstep+'.ncd'
vname='Potential'

; open netcdf file
ncid=ncdf_open(fname)
; get variable id;
varid=ncdf_varid(ncid,vname)
; get (x,y)-coordiates and dimension
ncdf_diminq,ncid,0,xname,xsize
ncdf_diminq,ncid,1,yname,ysize
; read data
data=fltarr(xsize,ysize)
ncdf_varget,ncid,varid,data
mzeta=ysize-1
print,'mzeta=',mzeta


for i=m1,m2 do begin
    zonal=0.0
    for j=1,mtheta(i) do begin
        ij=igrid(i)+j
        for k=1,mzeta do begin
            zonal=zonal+data(ij,k)
        endfor
    endfor
    zonal=zonal/(mzeta*mtheta(i))
    zonarr(i-m1,it-1)=zonal
endfor


print,'done with', fname

;;;plot,sp,charsize=2.0

endfor


fw=wdir+'zonalhist.dat'
openw,1,fw
printf,m1
printf,m2
printf,n1
printf,n2
printf,zonarr
close,1
end
