%% Lab 3
% Definisci il tempo di riferimento (ad esempio l'attuale data e ora)
r = [1 1 1]'.*1e7;

%background('Milky Way'); hold on;
opts = struct('Position',r.*4)
planet3D('Mars',opts); drawnow; hold on;
opts = struct('Position',r) 
planet3D('Earth',opts); drawnow; hold on;
%% State reconstruction problem
clear; clc;
r1 = [-21800 ; 37900 ; 0]; %km
r2 = [27300 ; 27700 ; 0 ];

%           h                     min              s
d_t =   15 * 3600       +       6 * 60      +     40;

background('Stars'); hold on;
mu = astroConstants(13); Nrev = 0; optionsLMR = 1; orbitType = 0; Ncase = 0;
[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambert.lambertMR(r1,r2,d_t,mu,orbitType,Nrev,Ncase,optionsLMR);
T = 2*pi*sqrt(A^3/mu);
t = linspace(0,T,1000);
[r1_prop,v] = coord.car_time_elapsed(t,r1,VI,mu);
plots.static_orbit(r1_prop);
[r2_prop,v,t] = coord.car_time_elapsed(t,r2,VF,mu);
plots.static_orbit(r2_prop);
plot3(r1(1),r1(2),r1(3),'Marker','o','MarkerSize',20,'LineWidth',10)
plot3(r2(1),r2(2),r2(3),'Marker','o','MarkerSize',20,'LineWidth',10)
%% 
clear; clc;
mu = astroConstants(13);
a1 = 12500; e1 = 0; i1 = 0; OM1 = 0; om1 = 0; th1 = deg2rad(120);
T1 = 2*pi*sqrt(a1^3/mu);
a2 = 9500; e2 = 0.3; i2 = 0; OM2 = 0; om2 = 0; th2 = deg2rad(250);
T2 = 2*pi*sqrt(a2^3/mu);
d_t = 3300;
background('Stars'); hold on;
opts = struct('Position',[0;0;0],'Units','km');
planet3D('Earth',opts); drawnow;

[r1,v1] = coord.kep2car_theta(a1,e1,i1,OM1,om1,th1,mu);
r1_prop = coord.car_time_elapsed([0 T1],r1,v1,mu);
[r2,v2] = coord.kep2car_theta(a2,e2,i2,OM2,om2,th2,mu);
r2_prop = coord.car_time_elapsed([0 T2],r2,v2,mu);

Nrev = 0; optionsLMR = 1; orbitType = 0; Ncase = 0;
[at,P,et,ERROR,vt1,vt2,TPAR,THETA] = lambert.lambertMR(r1,r2,d_t,mu,orbitType,Nrev,Ncase,optionsLMR);
Tt = 2*pi*sqrt(at^3/mu);
[~, ~, it, OMt, omt, th1t] = coord.car2kep_theta(r1,vt1,mu);
[~, ~, ~, ~, ~, th2t] = coord.car2kep_theta(r2,vt2,mu);
%

rt_prop = coord.kep2car_theta(at,et,it,OMt,omt,linspace(th1t,th2t+2*pi,1000),mu);

%
plot3(r1_prop(:,1),r1_prop(:,2),r1_prop(:,3),'LineWidth',2); hold on;
plot3(r2_prop(:,1),r2_prop(:,2),r2_prop(:,3),'LineWidth',2);
plot3(rt_prop(:,1),rt_prop(:,2),rt_prop(:,3),'LineWidth',2);
%%
clear
clc
t1_early = [2003, 4, 1, 0, 0, 0];
t1_late = [2003, 8, 1, 0, 0, 0];
t2_early = [2003, 9, 1, 0, 0, 0];
t2_late = [2004, 3, 1, 0, 0, 0];
[d_v_min,t_dep,t_arr] = lambert.lamb_cost(t1_early,t1_late,t2_early,t2_late,3,4,astroConstants(4));
%%
plots.solar_sys_optim_transfer(3,4,t_dep,t_arr,0.5e-3)
%%
clear
clc
t1_early = [2024, 6, 1, 0, 0, 0];
t1_late = [2026, 11, 1, 0, 0, 0];
t2_early = [2024, 12, 1, 0, 0, 0];
t2_late = [2027, 6, 1, 0, 0, 0];
[d_v_min,t_dep,t_arr] = lambert.lamb_cost(t1_early,t1_late,t2_early,t2_late,3,2,astroConstants(4));

%%
%t_dep = timeConversion.date2mjd2000(t_dep);
%t_arr = timeConversion.date2mjd2000(t_arr);
lambert.compute_delta_v(t_dep,t_arr,3,2,astroConstants(4))

%%
plots.solar_sys_optim_transfer(3,2,t_dep,t_arr,0.5e-3)
%% Mercury
clear
clc
t1_early = [2023, 11, 1, 0, 0, 0];
t1_late = [2025, 1, 1, 0, 0, 0];
t2_early = [2024, 4, 1, 0, 0, 0];
t2_late = [2027, 6, 1, 0, 0, 0];
[d_v_min,t_dep,t_arr] = lambert.lamb_cost(t1_early,t1_late,t2_early,t2_late,3,1,astroConstants(4),30);
%%
plots.solar_sys_optim_transfer(3,1,t_dep,t_arr,0.5e-3)
%% Venus
clear
clc
t1_early = [2024, 6, 1, 0, 0, 0];
t1_late = [2026, 11, 1, 0, 0, 0];
t2_early = [2024, 12, 1, 0, 0, 0];
t2_late = [2027, 6, 1, 0, 0, 0];
[d_v_min,t_dep,t_arr] = lambert.lamb_cost(t1_early,t1_late,t2_early,t2_late,3,2,astroConstants(4),100);
%%
plots.solar_sys_optim_transfer(3,2,t_dep,t_arr,0.5e-3)
%% Mars
clear
clc
t1_early = [2025, 8, 1, 0, 0, 0];
t1_late = [2031, 1, 1, 0, 0, 0];
t2_early = [2026, 1, 1, 0, 0, 0];
t2_late = [2032, 1, 1, 0, 0, 0];
[d_v_min,t_dep,t_arr] = lambert.lamb_cost(t1_early,t1_late,t2_early,t2_late,3,4,astroConstants(4),10);
%%
plots.solar_sys_optim_transfer(3,4,t_dep,t_arr,0.5e-3)
%% Jupiter
clear
clc
t1_early = [2026, 6, 1, 0, 0, 0];
t1_late = [2028, 6, 1, 0, 0, 0];
t2_early = [2028, 6, 1, 0, 0, 0];
t2_late = [2034, 1, 1, 0, 0, 0];
[d_v_min,t_dep,t_arr] = lambert.lamb_cost(t1_early,t1_late,t2_early,t2_late,3,5,astroConstants(4),20);
%%
plots.solar_sys_optim_transfer(3,5,t_dep,t_arr,0.5e-3)
%% Saturn
clear
clc
t1_early = [2027, 9, 1, 0, 0, 0];
t1_late = [2029, 10, 1, 0, 0, 0];
t2_early = [2030, 4, 1, 0, 0, 0];
t2_late = [2036, 3, 1, 0, 0, 0];
[d_v_min,t_dep,t_arr] = lambert.lamb_cost(t1_early,t1_late,t2_early,t2_late,3,6,astroConstants(4),20);
%%
plots.solar_sys_optim_transfer(3,6,t_dep,t_arr,0.5e-3)
%% Uranus
clear
clc
t1_early = [2027, 1, 1, 0, 0, 0];
t1_late = [2029, 1, 1, 0, 0, 0];
t2_early = [2031, 4, 1, 0, 0, 0];
t2_late = [2045, 12, 1, 0, 0, 0];
[d_v_min,t_dep,t_arr] = lambert.lamb_cost(t1_early,t1_late,t2_early,t2_late,3,7,astroConstants(4),20);
%%
plots.solar_sys_optim_transfer(3,7,t_dep,t_arr,0.5e-3)
%% Neptune
clear
clc
t1_early = [2025, 1, 1, 0, 0, 0];
t1_late = [2026, 10, 1, 0, 0, 0];
t2_early = [2036, 4, 1, 0, 0, 0];
t2_late = [2055, 6, 1, 0, 0, 0];
[d_v_min,t_dep,t_arr] = lambert.lamb_cost(t1_early,t1_late,t2_early,t2_late,3,8,astroConstants(4),20);
%%
plots.solar_sys_optim_transfer(3,8,t_dep,t_arr,0.5e-3)

