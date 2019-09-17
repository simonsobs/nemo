# Small footprint inside E-D56 - for quick tests
# NOTE: -W switch adds uK/pix white noise to avoid ringing in filtered maps
python3 healpix2CAR.py TOnly_la280.fits template_small.txt -O map -v -W 2 -L small_CAR
python3 healpix2CAR.py TOnly_la225.fits template_small.txt -O map -v -W 2 -L small_CAR
python3 healpix2CAR.py TOnly_la145.fits template_small.txt -O map -v -W 2 -L small_CAR
python3 healpix2CAR.py TOnly_la093.fits template_small.txt -O map -v -W 2 -L small_CAR

# AdvACT footprint - for testing big maps, tiling
python3 healpix2CAR.py TOnly_la280.fits template_AdvACT.txt -O map -v -W 2
python3 healpix2CAR.py TOnly_la225.fits template_AdvACT.txt -O map -v -W 2
python3 healpix2CAR.py TOnly_la145.fits template_AdvACT.txt -O map -v -W 2
python3 healpix2CAR.py TOnly_la093.fits template_AdvACT.txt -O map -v -W 2
