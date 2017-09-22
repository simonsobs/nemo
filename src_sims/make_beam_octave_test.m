addpath("/home/mlungu/work/170212/lib/")

foo  = myload('/data/saiola/depot/Beams/151107_v1/beam_tform_151107_v1_2014_pa2_jitter_deep56.txt');
foo2 = myload('/data/saiola/depot/Beams/151107_v1/beam_profile_151107_v1_2014_pa2_jitmap_deep56.txt');

ell = foo(:,1); bl = foo(:,2);
t   = foo2(:,1); bt = foo2(:,2);

wl = sinc(0.5*(ell/(2.*60.*180.)));

bt2 = transform2profile_c(t/180.*pi,bl*wl);


