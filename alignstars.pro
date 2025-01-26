PRO ALIGNSTARS,field_nr, chip

; PURPOSE: Use one (to a few, one is usually sufficient) reference stars per quadrant for a
; pre-alignment with VVV. Then search for more stars and refine the
; alignment.
; 
; This script loops over the four quadrants for each pointing.
; The pointing ('field') must be set by hand.

srej = 2.0 ; sigma rejection for RESISTANT_MEAN


; ------------EDIT HERE---------
; size of VVV reference image
; This is the printed output  the end of prep_refim.pro

xsize_ref = 1800
ysize_ref = 1800

chip_nr = strn(chip)
band = 'H'

ngrid = 60 ; Defines size of star grid for alignment 

; ----------------CAN BE EDITED, but USUALLY NOT NECESSARY ----------


; To create small aligned longexposure images
; So that we can cut away the unnecessary zeros around the images.

; size of dejittered HAWK-I long exposures
xsize_quad = 2600
ysize_quad = 2600

; VVV --> HAWKI
scale = 0.34/0.106

; define paths
VVV       ='/home/data/VVV/'
im_path   = '/home/data/GNS/2021/'+band+'/' + strn(field_nr) + 'HB/ims/'
data_path = '/home/data/GNS/2021/'+band+'/' + strn(field_nr) + 'HB/data/'
vvv_stars =  VVV +'/Fields/J/Field' + strn(field_nr) + '_stars.txt'
tmp_path  = '/home/data/GNS/2021/'+band+'/' + strn(field_nr) + 'HB/tmp/'

; Read list of reference stars in VVV
; and scale to HAWK-I pixel scale
; ------------------------------------
readcol, vvv_stars, x_vvv, y_vvv, a, d, J, Ks, FORMAT='F'
;readcol, vvv_stars, x_vvv, y_vvv, a, d, m, sx, sy, sm, c, FORMAT='F'
; add +1 to x and Y, because this improves the fit
x_vvv = x_vvv + 1
y_vvv = y_vvv + 1
x_vvv_scaled = x_vvv * scale
y_vvv_scaled = y_vvv * scale

; Define size of HAWK-I mosaic
; -----------------------------1
xsize_mosaic = round(xsize_ref * scale) + 50
ysize_mosaic = round(ysize_ref * scale) + 50
print, 'x, y size of HAWK-I mosaic: ' + strn(xsize_mosaic) + ', ' + strn(ysize_mosaic)

;Read list of stars in HAWK-I image
; ----------------------------------
 readcol, data_path + 'stars_' + chip_nr + '.txt', x, y, f
 ; add +1 to x and Y, because this improves the fit
x = x + 1
y = y + 1
n_stars = n_elements(f)

; read lists of common stars found by align_VVV.py
 ; -------------------------------------------------
readcol, data_path + 'aa_stars_hawki_' + chip_nr + '.txt', x_ref, y_ref
readcol, data_path + 'aa_stars_vvv_' + chip_nr + '.txt', x_ref_vvv, y_ref_vvv
x_ref_vvv_scaled = x_ref_vvv * scale
y_ref_vvv_scaled = y_ref_vvv * scale
print, 'Common stars from astroalign: ' + strn(n_elements(x_ref)) + ', ' + strn(n_elements(x_ref_vvv))

; Transform and find common stars
; -------------------------------

 print, 'Initial coarse alignment.'
 degree = 1
 dmax = 2.0
 polywarp, x_ref_vvv_scaled, y_ref_vvv_scaled, x_ref, y_ref, degree, Kx, Ky
 print, Kx
 print, Ky
 xi = Kx[0,0] + Kx[0,1]*x + Kx[1,0]*y + Kx[1,1]*x*y
 yi = Ky[0,0] + Ky[0,1]*x + Ky[1,0]*y + Ky[1,1]*x*y
 compare_lists, x_vvv_scaled, y_vvv_scaled, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
 nc = n_elements(subc1)
 print, 'Found ' + strn(nc) + ' common stars.'

 
; print region file of preliminary alignment stars
openw, out, im_path + 'prelimiinary' + chip_nr + '.reg', /get_lun
printf, out, 'global color=cyan dashlist=8 3 width=2 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
printf, out, 'image'
for s = 0, nc-1 do begin
   printf, out, 'circle('+ strn(x[subc2[s]]+1)+','+strn(y[subc2[s]]+1)+',5)'
endfor
free_lun, out

 ; iterative degree 1 alignment
 ; ------------------------------

  degree = 1
  dmax = 2.0
  for it = 0, 4 do begin
    grid_idx = gridselect(x[subc2], y[subc2],f[subc2],gridsize=ngrid)
    polywarp, x_vvv_scaled[subc1[grid_idx]], y_vvv_scaled[subc1[grid_idx]], x[subc2[grid_idx]], y[subc2[grid_idx]], degree, Kx, Ky
    print, Kx
    print, Ky
    xi = Kx[0,0] + Kx[0,1]*x + Kx[1,0]*y + Kx[1,1]*x*y
    yi = Ky[0,0] + Ky[0,1]*x + Ky[1,0]*y + Ky[1,1]*x*y
    compare_lists, x_vvv_scaled, y_vvv_scaled, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
    nc = n_elements(subc1)
    print, 'Iteration ' + strn(it)
    print, 'Found ' + strn(nc) + ' common stars.'
  endfor

 ; iterative higher degree alignment
 ; ------------------------------

  ; next line serves to switch higher degree fit on/off
  if 1 then begin

   degree = 2
   dmax = 2.0
 
   print, 'Now Degree ' + strn(degree) + ' alignment.'
   for it = 0, 4 do begin
    grid_idx = gridselect(x[subc2], y[subc2],f[subc2],gridsize=ngrid)
    polywarp, x_vvv_scaled[subc1[grid_idx]], y_vvv_scaled[subc1[grid_idx]], x[subc2[grid_idx]], y[subc2[grid_idx]], degree, Kx, Ky

    print, Kx
    print, Ky
    xi = replicate(0.0,n_stars)
    yi = replicate(0.0,n_stars)
    for k = 0, degree do begin
      for m = 0, degree do begin
        xi = xi + Kx[m,k]*x^k*y^m
        yi = yi + Ky[m,k]*x^k*y^m
      endfor
    endfor
    compare_lists, x_vvv_scaled, y_vvv_scaled, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
    nc = n_elements(subc1)
    print, 'Iteration ' + strn(it+1)
    print, 'Found ' + strn(nc) + ' common stars.'
    ; Outlier rejection
    dp = sqrt((x1c-x2c)^2 + (y1c-y2c)^2)
    RESISTANT_Mean, dp, srej, dp_mean, dp_Sigma, Num_RejECTED, GOODVEC=goodvec
    print, 'Rejected ' + strn(Num_RejECTED) + ' outliers, retaining ' + strn(n_elements(goodvec)) + ' stars.'
    subc1 = subc1[goodvec]
    subc2 = subc2[goodvec]
   endfor
endif



  ; Kx and Ky describe second order polynomial transformation into
 ; a rescaled VVV (to HAWK-I pixel scale) reference frame
 ; ------------------------------------------------------------

 grid_idx = gridselect(x[subc2], y[subc2],f[subc2],gridsize=ngrid)
 x_grid =  x[subc2[grid_idx]]
 y_grid =  y[subc2[grid_idx]]
 x_vvv_grid = x_vvv_scaled[subc1[grid_idx]]
 y_vvv_grid = y_vvv_scaled[subc1[grid_idx]]

; Compute transformed HAWK-I star positions
; inside the scaled VVV frame.
; For testing the goodness of the alignment, 
; create an artificial map and save the list of transformed positions.
; -------------------------------------------------------------------
polywarp, x_vvv_grid, y_vvv_grid, x_grid, y_grid, degree,  Kx, Ky
n_stars = n_elements(f)
x_fin = replicate(0.0,n_stars)
y_fin = replicate(0.0,n_stars)
for k = 0, degree do begin
  for m = 0, degree do begin
    x_fin = x_fin + Kx[m,k]*x^k*y^m
    y_fin = y_fin + Ky[m,k]*x^k*y^m
  endfor
endfor
dat = ptr_new({X_size: 10, Y_size: 10, Sigma_x: 0.8, Sigma_y: 0.8, Angle: 0.0})
model = image_model(x_fin,y_fin,f,xsize_mosaic,ysize_mosaic,'gaussian', dat)
writefits, im_path + 'lnx_map' + chip_nr + '.fits', model, /COMPRESS

forprint, TEXTOUT=data_path + 'stars_transformed_' + chip_nr + '_d' + strn(degree) + '.txt', x_fin, y_fin, f, /NOCOMMENT


; Plot distributions of differences of
; astrometric  positions between reference stars in HAWK-I and VVV
; -----------------------------------------------------------------

dx_tmp =  (x_vvv_scaled[subc1] - x_fin[subc2])*0.34
dy_tmp = (y_vvv_scaled[subc1] - y_fin[subc2])*0.34
set_plot,'PS',/interpolate
device, XOFFSET=0, YOFFSET=0, $
      FILENAME =  data_path + 'alignment_uncertainty_' + strn(chip) + '.eps', XSIZE=20., YSIZE=12., $
      /portrait, /color,BITS_PER_PIXEL=8, encapsulated=1
!P.MULTI=[0,2,1]
!P.CHARSIZE=1.2
!P.THICK=4.0
!P.CHARTHICK=4.0

binsize = 0.01
cgHistoplot, dx_tmp, /FILL, HISTDATA=h1, LOCATIONS=loc, BINSIZE=binsize, XTITLE = 'Difference in X-axis (")', YTITLE = 'Number of Stars',xrange = [-0.5,0.5], xticks = 4

binCenters1 = loc + (binsize / 2.0)
yfit = GaussFit(binCenters1, h1, coeff, NTERMS=3, sigma = error)
cgPlot, binCenters1, yfit, COLOR='dodger blue', THICK=2, /OVERPLOT
  
centerfit = String(coeff[1], FORMAT='(F0.3)')
cgText, 0.17, 0.84, /NORMAL, 'Center ' + centerfit + ' arcsec', COLOR='navy'
cgText, 0.17, 0.80, /NORMAL, 'Error ' + strn(coeff[2], FORMAT='(F0.3)') + ' arcsec', COLOR='navy'

cgHistoplot, dy_tmp, /FILL, HISTDATA=h1, LOCATIONS=loc, BINSIZE=binsize, XTITLE = 'Difference in Y-axis (")' , YTITLE = '' ,xrange = [-0.5,0.5], xticks = 4
binCenters1 = loc + (binsize / 2.0)

yfit = GaussFit(binCenters1, h1, coeff, NTERMS=3, sigma = error)
cgPlot, binCenters1, yfit, COLOR='dodger blue', THICK=2, /OVERPLOT
  
centerfit = String(coeff[1], FORMAT='(F0.3)')
cgText, 0.65, 0.84, /NORMAL, 'Center ' + centerfit + ' arcsec', COLOR='navy'
cgText, 0.65, 0.80, /NORMAL, 'Error ' + strn(coeff[2], FORMAT='(F0.3)') + ' arcsec', COLOR='navy'

device, /close


; determine transformation parameters for image and save them
; for later use in alignframes.pro
; The order of VVV and HAWK-I data points must be reversed for the image transforms.
; -------------------------------------------------------------------------------
 polywarp, x_grid, y_grid, x_vvv_grid, y_vvv_grid, degree, Kx, Ky
 SAVE, Kx, Ky, FILENAME=data_path + 'AlignPars_chip' + chip_nr

; print region file of reference stars on HAWK-I image
openw, out, im_path + 'alignstars_' + chip_nr + '.reg', /get_lun
printf, out, 'global color=cyan dashlist=8 3 width=3 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1;'
printf, out, 'image'
n_stars = n_elements(grid_idx)
for s = 0, n_stars-1 do begin
   printf, out, 'circle('+ strn(x_grid[s])+','+strn(y_grid[s])+',5)'
endfor
free_lun, out



 ; transform image and mask
 ; mask pixels that have coverage < 0.9 to avoid 
 ; later problems with division by small numbers
 ; ---------------------------------------------

 im = readfits(im_path + 'lnx_jitter_'+chip_nr + '_wt.fits.gz')
 transim = POLY_2D(im,Kx,Ky,2,xsize_mosaic,ysize_mosaic,CUBIC=-0.5,MISSING=0)
 writefits, im_path + 'lnx_jitter_'+chip_nr+ '_aligned_wt.fits.gz', transim, /COMPRESS

 im = readfits(im_path + 'lnx_jitter_'+chip_nr + '.fits.gz')
 transim = POLY_2D(im,Kx,Ky,2,xsize_mosaic,ysize_mosaic,CUBIC=-0.5,MISSING=0)
 writefits, im_path + 'lnx_jitter_'+chip_nr+ '_aligned.fits.gz', transim, /COMPRESS

 transim_im = transim

 im = readfits(im_path + 'lnx_jitter_'+chip_nr + '_sig.fits.gz')
 transim = POLY_2D(im,Kx,Ky,2,xsize_mosaic,ysize_mosaic,CUBIC=-0.5,MISSING=0 )
 transim_noise = transim
 writefits, im_path + 'lnx_jitter_'+chip_nr+ '_aligned_sig.fits.gz', transim, /COMPRESS

 ; to check alignment with VVV, create
 ; a transformed image inside the VVV frame of reference
 ; ---------------------------------------------
 polywarp,   x[subc2[grid_idx]], y[subc2[grid_idx]], x_vvv[subc1[grid_idx]], y_vvv[subc1[grid_idx]], degree, Kx, Ky
 im = readfits(im_path + 'lnx_jitter_'+chip_nr + '.fits.gz')
 transim = POLY_2D(im,Kx,Ky,2,xsize_ref,ysize_ref,CUBIC=-0.5,MISSING=0)
 writefits, im_path + 'lnx_jitter_'+chip_nr+ '_VVV.fits.gz', transim, /COMPRESS

 im = readfits(im_path + 'lnx_jitter_'+chip_nr + '_wt.fits.gz')
 transim = POLY_2D(im,Kx,Ky,2,xsize_ref,ysize_ref,CUBIC=-0.5,MISSING=0)
 writefits, im_path + 'lnx_jitter_'+chip_nr+ '_VVV_wt.fits.gz', transim, /COMPRESS


; Now cut out the aligned HAWK-I images from the large frames to get
; rid of the edges
; --------------------------------------------------------------------

xlo = x_off[chip-1]
ylo = y_off[chip-1]
xhi = x_off[chip-1] + xsize_quad - 1
yhi = y_off[chip-1] + ysize_quad - 1

lnx = transim_im[xlo:xhi,ylo:yhi]
lnx_noise = transim_noise[xlo:xhi,ylo:yhi]
  
writefits, data_path + 'lnx_aligned_' + chip_nr + '.fits.gz', lnx, /COMPRESS
writefits, data_path + 'lnx_aligned_' + chip_nr + '_sig.fits.gz', lnx_noise, /COMPRESS


;To have an idea of the displacement that we have, we take the initial
;list and the corrected one to see the difference in positions.
; These plots can be used to detect systematic residual distortions.
; ------------------------------------------------------------------------------

; median global offset
x_dif = median(x_fin - x)
y_dif = median(y_fin - y) 
 
; origin of arrows in displaced image
x0_new = x_dif + x
y0_new = y_dif + y
 
; max size of image
x_dis = max(x_fin-x0_new)
y_dis = max(y_fin-y0_new)
 

print,  strn(x_dis*0.106) + 'arcsec    ' + strn(y_dis*0.106) + 'arcsec'

dx = x_dis*0.106
dy = y_dis*0.106

;forprint, TEXTOUT= data_path + 'displ_arcsec_distortion_chip' + chip_nr + '.txt', dx, dy, /NOCOMMENT 



;Drawing the detector


set_plot,'PS', /interpolate

device, XOFFSET=0, YOFFSET=0, $
   FILENAME= data_path + 'Alignment_resid_chip' + chip_nr +'.eps', XSIZE=30., YSIZE=30., $
   /portrait, /color, BITS_PER_PIXEL=8, encapsulated=1
!P.MULTI=[0,1,0]
!P.CHARSIZE=1.8
!X.THICK=4
!Y.THICK=4
!P.THICK=4.0
!P.CHARTHICK=4

!P.COLOR=0

;cgoplot, x1c, y1c, Color='orange', PSYM=3

; compute shift between common re-scaled VVV and HAWK-I positions
; The arrow will point to the VVV position
; -----------------------------------------------------------------
xf = (-xi[subc2]+x1c)*10 
yf = (-yi[subc2]+y1c)*10
; subtract global shift
xf = xf - median(xf)
yf = yf - median(yf)
x0 = xi[subc2] - xlo
y0 = yi[subc2] - ylo
xf = xf+xi[subc2] - xlo
yf = yf+yi[subc2] - ylo

ran = RANDOMU(Seed, n_elements(xf))

num = sort(ran)
max_num = n_elements(ran) < 80
num = num[0:max_num-1]
a = xi[subc2[num]]
b = yi[subc2[num]]

cgplot, a, b, Color='black', XRANGE = [0,xsize_quad], YRANGE = [0,ysize_quad], XTITLE='X-Axis [pixels]', YTITLE='Y-Axis [pixels]', XSTYLE=1, PSYM=3, YSTYLE=1


ARROW, x0[num], y0[num], xf[num], yf[num], /DATA, COLOR=4, HSIZE=100


device, /close




END
