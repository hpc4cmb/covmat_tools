#! /bin/csh -f

if ($#argv < 3) then
  echo "usage: $0 <nside> <nstokes> <matrix1> <matrix2> [<matrix3> [<matrix4> ... ]]"
else

cp $3 full.dat
echo $3

foreach covmat ($argv[4-])
    echo " + $covmat"
    merge_covmat \
	full.dat $1 $2 $covmat next.dat > /dev/null
    mv -f next.dat full.dat
end

echo " == full.dat"

endif
