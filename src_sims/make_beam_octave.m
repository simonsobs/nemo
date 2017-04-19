addpath("/home/mlungu/work/170212/lib/")

foo  = myload('tform_gauss_1.3arcmin.txt');
foo2 = myload('/data/saiola/depot/Beams/151107_v1/beam_profile_151107_v1_2014_pa2_jitmap_deep56.txt');

ell = foo(:,1); bl = foo(:,2);
t   = foo2(:,1); %I take the same formatting as in Mariu's beam files.

wl = sinc(0.5*(ell/(2.*60.*180.)));

bt2 = transform2profile_c(t/180.*pi,bl.*wl);


fid = fopen('profile_gauss_1.3arcmin_pixwin.txt','w');
for i = 1:length(bt2)
   fprintf(fid,"%0.6f  %1.5e\n", t(i), bt2(i));
endfor
fclose(fid);
