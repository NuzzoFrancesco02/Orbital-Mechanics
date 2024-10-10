%% Lab_2
% ... Exercise 2.1.b
%% Assigned kep coordinates
clear; clc;
a = 83500; e = 0.1976; i = 60; OM = 270; w = 45; th0 = 230;
mu = astroConstants(13); 
T = 2*pi*sqrt(a^3/mu);
t = linspace(0,2*T,1000);

groundTrack_mine([],[a e deg2rad([i OM w th0])],0,deg2rad(15.05)/3600,t,mu);
%% Molniya orbit
clear; clc;
a = 26600; e =  0.74; i = 63.4; OM = 270; w = 270; th0 = 0;
mu = astroConstants(13); 
T = 2*pi*sqrt(a^3/mu);
t = linspace(0,2*T,1000);

groundTrack_mine([],[a e deg2rad([i OM w th0])],0,deg2rad(15.05)/3600,t,mu);

startTime = datetime(2020,5,11,12,35,38);
stopTime = startTime + days(2);
sampleTime = 60;
sc = satelliteScenario(startTime,stopTime,sampleTime);

sat1 = satellite(sc,a*1e3,e,i,OM,w,th0,"Name",'Sat1');
groundTrack(sat1);
%% Circular LEO i = 0
clear; clc;
a = 300; e =  0; i = 0; OM = 45; w = 45; th0 = 0;
mu = astroConstants(13); 
T = 2*pi*sqrt(a^3/mu);
t = linspace(0,2*T,1000);

groundTrack_mine([],[a e deg2rad([i OM w th0])],0,deg2rad(15.05)/3600,t,mu);
%% Circular LEO i = 98
clear; clc;
a = 300; e =  0; i = 98; OM = 45; w = 45; th0 = 0;
mu = astroConstants(13); 
T = 2*pi*sqrt(a^3/mu);
t = linspace(0,2*T,1000);

groundTrack_mine([],[a e deg2rad([i OM w th0])],0,deg2rad(15.05)/3600,t,mu);
%% Circular LEO i = 30
clear; clc;
a = 300; e =  0; i = 30; OM = 45; w = 45; th0 = 0;
mu = astroConstants(13); 
T = 2*pi*sqrt(a^3/mu);
t = linspace(0,2*T,1000);

groundTrack_mine([],[a e deg2rad([i OM w th0])],0,deg2rad(15.05)/3600,t,mu);
%%
clear; clc;
a = 35786; e =  0; i = 30; OM = 45; w = 45; th0 = 0;
mu = astroConstants(13); 
T = 2*pi*sqrt(a^3/mu);
t = linspace(0,4*T,1000);

groundTrack_mine([],[a e deg2rad([i OM w th0])],0,deg2rad(15.05)/3600,t,mu);
%%

sat1 = satellite(sc,a*1e3,e,i,OM,w,th0,"Name",'Sat1');
groundTrack(sat1)
play(sc)
%%
h = plot3(nan,nan,nan,'or','MarkerFaceColor',"#77AC30",'MarkerEdgeColor',"#77AC30",'MarkerSize',10);
hold on;
traj = plot3(nan,nan,nan,'Color',"#A2142F");
plot3(X,Y,Z,'Color','none');
% Define the step animation
step_animation = 10;
for j = 1 : 3 
    for i = 1:10:length(X)
        set(traj,'XData',X(1:i),'YData',Y(1:i),'ZData',Z(1:i));
        hold on;
        set(h,'XData',X(i),'YData',Y(i),'ZData',Z(i));
        drawnow
        
        pause(0.1);
        
    end
    hold off
end
%% 
clear
clc
e = 0.37255;
mu = astroConstants(13);
t = 1; n = 3.6029; 
a = nthroot(mu/n^2,3)
t0 = 0; th0 = 0; 
[th, E, M] = kepler_equation(t, e, a, mu, t0, th0)
%%
r0 = [7000 -12124 0]'; v0 = [2.6679 4.6210 0]'; t = 3600;
[r,v] = car2ECI([0 t],r0,v0,astroConstants(13));
r(end,:)
v(end,:)