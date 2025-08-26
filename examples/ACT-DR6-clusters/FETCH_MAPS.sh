# Fetches the DR6 maps used for the ACT DR6 cluster catalog paper
mkdir -p maps
cd maps
curl https://lambda.gsfc.nasa.gov/data/act/maps/published/act-planck_dr4dr6_coadd_AA_daynight_f090_map.fits -o act-planck_dr4dr6_coadd_AA_daynight_f090_map.fits || wget -c https://lambda.gsfc.nasa.gov/data/act/maps/published/act-planck_dr4dr6_coadd_AA_daynight_f090_map.fits
curl https://lambda.gsfc.nasa.gov/data/act/maps/published/act-planck_dr4dr6_coadd_AA_daynight_f150_map.fits -o act-planck_dr4dr6_coadd_AA_daynight_f150_map.fits || wget -c https://lambda.gsfc.nasa.gov/data/act/maps/published/act-planck_dr4dr6_coadd_AA_daynight_f150_map.fits
curl https://lambda.gsfc.nasa.gov/data/act/maps/published/act-planck_dr6.02_coadd_AA_daynight_f220_map.fits -o act-planck_dr6.02_coadd_AA_daynight_f220_map.fits || wget -c https://lambda.gsfc.nasa.gov/data/act/maps/published/act-planck_dr6.02_coadd_AA_daynight_f220_map.fits
curl https://lambda.gsfc.nasa.gov/data/act/maps/published/act-planck_dr4dr6_coadd_AA_daynight_f090_ivar.fits -o act-planck_dr4dr6_coadd_AA_daynight_f090_ivar.fits || wget -c https://lambda.gsfc.nasa.gov/data/act/maps/published/act-planck_dr4dr6_coadd_AA_daynight_f090_ivar.fits
curl https://lambda.gsfc.nasa.gov/data/act/maps/published/act-planck_dr4dr6_coadd_AA_daynight_f150_ivar.fits -o act-planck_dr4dr6_coadd_AA_daynight_f150_ivar.fits || wget -c https://lambda.gsfc.nasa.gov/data/act/maps/published/act-planck_dr4dr6_coadd_AA_daynight_f150_ivar.fits
curl https://lambda.gsfc.nasa.gov/data/act/maps/published/act-planck_dr6.02_coadd_AA_daynight_f220_ivar.fits -o act-planck_dr6.02_coadd_AA_daynight_f220_ivar.fits || wget -c https://lambda.gsfc.nasa.gov/data/act/maps/published/act-planck_dr6.02_coadd_AA_daynight_f220_ivar.fits
cd ..
curl -L "https://dl.dropbox.com/scl/fi/op52a98mpnygpsidxqhvh/ACTDR6ClustersExtras.tar.gz?rlkey=odc6eqyjih01hg27vqxzv90ah&st=8cue3f99&dl=0" -o ACTDR6ClustersExtras.tar.gz || wget -c "https://dl.dropbox.com/scl/fi/op52a98mpnygpsidxqhvh/ACTDR6ClustersExtras.tar.gz?rlkey=odc6eqyjih01hg27vqxzv90ah&st=8cue3f99&dl=0" -O ACTDR6ClustersExtras.tar.gz
tar -zxvf ACTDR6ClustersExtras.tar.gz
