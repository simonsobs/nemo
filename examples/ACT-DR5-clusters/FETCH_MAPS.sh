mkdir -p maps
cd maps
wget -c https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr5/maps/act_dr5.01_s08s18_AA_f090_daynight_map.fits
wget -c https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr5/maps/act_dr5.01_s08s18_AA_f090_daynight_ivar.fits
wget -c https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr5/maps/act_dr5.01_s08s18_AA_f150_daynight_map.fits
wget -c https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr5/maps/act_dr5.01_s08s18_AA_f150_daynight_ivar.fits
cd ..
