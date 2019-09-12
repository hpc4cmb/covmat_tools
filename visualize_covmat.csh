#! /bin/csh -f

# This shell script processes a fits format covariance matrix
# and produces gif plots of the correlation of pixel # 0
# (NESTED scheme)
#
# The script requires 
#  - plot_pixel_covariance
#  - map2gif
#
# March 5th 2008 Reijo Keskitalo 

if ($#argv < 3) then
    echo "usage: $0 <covmat_file> <nside> <nstokes> [<pixel>]"
else

set covmat  = $1
set nside   = $2
set nstokes = $3
if ($#argv == 4) then
    set pixel = $4
else
    set pixel = 0
endif

echo " Reading from $covmat"
echo ""
echo " Using   Nside = $nside"
echo "       nstokes = $nstokes"
echo "        #pixel = $pixel"

# Get the fits maps corresponding to the specified column (pixel). 
# Store the corresponding pixel numbers in 'rows' array
set rows = `plot_pixel_covariance $covmat $nside $nstokes $pixel | \
    awk '/row #/ {print $4}'`

set file_root = ` echo $covmat | sed 's/.fits//g' `

# Clean the filenames of obsolete suffixes (.dat, .fits)
foreach srow ( $rows )
    mv -f ${file_root}.fits_row_$srow.fits ${file_root}_row_$srow.fits
end


# plot the maps in gif files
map2gif -inp ${file_root}_row_$rows[1].fits \
    -out \!${file_root}_row_$rows[1].II.gif -bar .true.
map2gif -inp ${file_root}_row_$rows[1].fits \
    -out \!${file_root}_row_$rows[1].IQ.gif -sig 2 -bar .true.
map2gif -inp ${file_root}_row_$rows[1].fits \
    -out \!${file_root}_row_$rows[1].IU.gif -sig 3 -bar .true.

map2gif -inp ${file_root}_row_$rows[2].fits \
    -out \!${file_root}_row_$rows[2].QQ.gif -sig 2 -bar .true.
map2gif -inp ${file_root}_row_$rows[2].fits \
    -out \!${file_root}_row_$rows[2].QU.gif -sig 3 -bar .true.

map2gif -inp ${file_root}_row_$rows[3].fits \
    -out \!${file_root}_row_$rows[3].UU.gif -sig 3 -bar .true.


endif # master loop
