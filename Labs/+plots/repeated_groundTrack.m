function a = repeated_groundTrack(k,m,om_E,mu)

n = om_E*k/m;
a = nthroot(mu/n^2,3);
