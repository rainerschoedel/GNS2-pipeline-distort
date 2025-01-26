PRO SORT_SUBCUBES, field, chip

; Sort subcubes by FoV so that they can be pased to holo_full.pro

chip_nr = strn(chip)
cube_path = '/home/data/GNS/2021/H/'+strn(field)+'HB/subcubes/'

x_cube = 600 ; xaxis length of sub-cube
y_cube = 600 ; yaxis length of sub-cube
x_large = 4800  ; xaxis length of large cube
y_large = 4800  ; yaxis length of large cube
 
x_sub_shift = x_cube/2
y_sub_shift = y_cube/2
nx = x_large/x_sub_shift - 1
ny = y_large/y_sub_shift - 1

for i_x = 0, nx -1 do begin
  for i_y = 0, ny -1 do begin
   filenam = '_' + strn(i_x) + '_' + strn(i_y)
   indir =  cube_path  + 'chip' + chip_nr + '/'
   innam = 'cube' + filenam + '.fits'
   masknam = 'masks' + filenam + '.fits'
   
   if FILE_TEST(indir + innam) then begin
     cube = readfits(indir + innam)
     masks = readfits(indir + masknam)
     SORT_MASKS, masks, cube
     writefits, indir + innam, cube, /COMPRESS
     writefits, indir + masknam, masks, /COMPRESS
     print, 'Processed cube ' + filenam + '.'
   endif
  endfor
endfor

END
