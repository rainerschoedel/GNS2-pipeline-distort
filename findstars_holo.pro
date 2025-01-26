; Auxiliary routine taken almost literally from StarFinder's
; take a close look at selected reference stars
; with the option to reject individual sources
; XPSF_EXTRACT_CONFIRM

PRO confirm_stars, image, wnum, display_opt, x_in, y_in, $
			  psfsize, x_out, y_out

	on_error, 2
	if  n_elements(x_in) eq 0 or n_elements(y_in) eq 0  then  return
   	sub_arrays, image, x_in, y_in, psfsize, stack
   	nstars = n_elements(x_in)
   	for  n = 0L, nstars - 1  do begin
   	   xn = -1  &  yn = -1
   	   opt = default_display_opt(stack[*,*,n])
   	   opt.reverse = display_opt.reverse
   	   opt.stretch = display_opt.stretch
   	   opt.color_table = display_opt.color_table
   	   display_image, stack[*,*,n], wnum, OPTIONS = opt
   	   msg = dialog_message('Confirm this star?', /QUESTION)
   	   if  strlowcase(msg) eq 'no'  then begin
   	      x_in[n] = -1  &  y_in[n] = -1
   	   endif
   	endfor
   	w = where(x_in ge 0 and y_in ge 0, n_confirm)
   	if  n_confirm ne 0  then begin
   	   x_out = x_in[w]  &  y_out = y_in[w]
   	endif
	display_image, image, wnum, OPTIONS = display_opt
	return
END

; ==============================

PRO FINDSTARS_HOLO, field, chip


band = 'H'
basedir = '/home/data/GNS/2021/H/'
in_path = basedir + strn(field) + 'HB/ims/'
data_path = basedir + strn(field) + 'HB/data/'
tmp_path = basedir + strn(field) + 'HB/tmp/'


; General StarFinder settings
; ---------------------------
psf_size = 40.
maskrad = 15
back_box = 0.
deblend = 0
min_correlation = 0.8
weigh_med = 0
unweighted = 1
upper_lev = 3.0e4
n_fwhm_back = 9.0
n_fwhm_fit = 2.0
n_fwhm_match = 1.0
mag_fac = 2L
n_width = 3.0
norm_max = 1
correl_mag = 2.0
niter = 2
compbg = 0
rel_thresh = 1
guide_x = ""
guide_y = ""


chip_nr = strn(chip)
im = readfits(in_path + 'lxp.fits.gz')
im = im[xlo:xhi,ylo:yhi]
sz = size(im)
n1 = sz[1]
n2 = sz[2]
noise = readfits(in_path + 'lxp_sigma.fits.gz')
noise = noise[xlo:xhi,ylo:yhi]
writefits, data_path + 'lxp_' + chip_nr + '.fits', im, /COMPRESS
writefits, data_path + 'lxp_' + chip_nr + '_sigma.fits', noise, /COMPRESS

width = 2400L
if (chip eq 1) then begin
  xlo = 0 & xhi = width -1 
  ylo = 0 & yhi = width -1 
endif
if (chip eq 2) then begin
  xlo = width & xhi = 2*width -1
  ylo = 0 & yhi = width -1 
endif 
if (chip eq 4) then begin
  xlo = width & xhi = 2*width -1
  ylo = width & yhi = 2*width -1
endif 
if (chip eq 3) then begin
  xlo = 0 & xhi = width -1 
  ylo = width & yhi = 2*width -1
endif 

 ; extract PSF
 ; --------------------------

 psf_fwhm = 4.0
 mask = replicate(1,n1,n2)    ; in principle we do not need to mask anything
 saturated = where(im gt upper_lev)
 mask[saturated] = 0
 writefits, 'mask.fits', mask
 psf = extractpsf(im, noise, mask, maskrad, psf_fwhm, tmp_path, /UNWEIGHTED, SF_THRESH = 100., SAT_LEVEL=upper_lev)
 writefits, data_path + 'psf_lxp_'+ strn(chip) + '.fits' , psf
 print, 'Chip ' + chip_nr + '.txt, PSF FWHM: ' + strn(psf_fwhm*0.106) + ' "'

;   STOP
 ; run StarFinder
 ; --------------------------------------------

threshold = [10.,10.]
back_box = maskrad
deblend = 1
deblost = 0
starfinder, im, psf, X_BAD=x_bad, Y_BAD = y_bad, $
        BACKGROUND = background, BACK_BOX = back_box, $
        threshold, REL_THRESHOLD = rel_thresh, /PRE_SMOOTH, $
        NOISE_STD = noise, min_correlation, $
        CORREL_MAG = correl_mag, INTERP_TYPE = 'I', $
	ESTIMATE_BG = compbg, DEBLEND = deblend, DEBLOST=deblost, N_ITER = niter, SILENT=0, $
	GUIDE_X = guide_x, GUIDE_Y = guide_y, $
	SV_SIGMA_R = 0.0085, SV_SIGMA_A= 0.0050, $
      	x, y, f, sx, sy, sf, c, STARS = stars, $
        LOGFILE = logfilename, /CUBIC
forprint, TEXTOUT= data_path + 'stars' + '_' + chip_nr + '_holo.txt', COMMENT = ' x y f c', x, y, f, c
; OPTIONAL OUTPUT FOR DEBUGGING
writefits, data_path + 'lxp_resid' + '_' + chip_nr + '.fits' , im - stars, /COMPRESS
writefits, data_path + 'lxp_stars' + '_' + chip_nr + '.fits' , stars, /COMPRESS
   dat = ptr_new({X_size: 20, Y_size: 20, Sigma_x: 1.5, Sigma_y: 1.5, Angle: 0.0})
map = image_model(x,y,f,n1,n2,'gaussian', dat)
writefits, data_path + 'lxp_map' + '_' + chip_nr + '.fits', map, /COMPRESS

END
