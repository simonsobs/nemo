mkdir -p maps
cd maps
wget -c https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr5/maps/act_dr5.01_s08s18_AA_f090_daynight_map.fits
wget -c https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr5/maps/act_dr5.01_s08s18_AA_f090_daynight_ivar.fits
wget -c https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr5/maps/act_dr5.01_s08s18_AA_f150_daynight_map.fits
wget -c https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr5/maps/act_dr5.01_s08s18_AA_f150_daynight_ivar.fits
cd ..
#wget -c https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr5/maps/act_dr5.01_auxilliary.zip
#unzip act_dr5.01_auxilliary.zip
wget -c https://astro.ukzn.ac.za/~mjh/nemo-dr5-masking.tar.gz
tar -zxvf nemo-dr5-masking.tar.gz
