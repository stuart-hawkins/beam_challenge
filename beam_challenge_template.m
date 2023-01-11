% Stuart C. Hawkins - 11 January 2023

% Note: in this code vectors (x,y) are represented by complex numbers 
% x + 1i * y

clear all
close all

%--------------------------------------------
% set main parameters
%--------------------------------------------

% wave speed
c = 1;

% angular frequency
omega = 1;

%--------------------------------------------
% set dependent parameters
%--------------------------------------------

% wavenumber
kwave = omega/c;

% wavelength
lambda = 2*pi/kwave;

% set the domain [D(1) D(2)] x [D(3) D(4)] on which to plot the solution
D = [2-10 2+10 8-10 8+10];

% setup a sensor at the point (0,10) with radius 0.5 to measure the field
s = sensor(10i,0.5);

%--------------------------------------------
% setup object for solving the wave problem
%--------------------------------------------

% setup incident wave to be a beam
inc = 1*beam(pi/4,lambda,kwave);

% setup solver with the incident wave
obj = tmatsolver(inc);

% plot the incident wave in figure(1)
figure(1)
tmplot(D,inc);

% plot the sensor too
hold on
s.plot('k--')
hold off

% make the figure look nice
axis equal
shading interp
view([0 90])
colorbar

% add a title to the figure
title('incident wave and sensor')

%--------------------------------------------
% setup the configuration of particles
%--------------------------------------------

% array of types of scatterer... in this case we will only use sound soft
% circular scatterers with radius=0.2 ie we only have one type of scatterer
types = [sound_soft_disk(kwave,0.2)];

% create a line of closely spaced scatterers
for y = 1:0.401:9.5
    obj.addParticle(particle(types(1),5+1i*y));
end

% plot the configuration in figure(2)
figure(2)
obj.schematic

% plot the sensor too
hold on
s.plot('k--')
hold off

% make the figure look nice
axis equal
view([0 90])

% add a title to the figure
title('configuration and sensor')

%--------------------------------------------
% solve the scattering problem
%--------------------------------------------

obj.solve()

%--------------------------------------------
% plot the total field
%--------------------------------------------

% plot the real part of the total field in figure(3)
figure(3)
obj.plot(D)

% plot the sensor too
hold on
s.plot('k--')
hold off

% make the figure look nice
axis equal
shading interp
view([0 90])
colorbar

% add a title to the figure
title('total field')

%--------------------------------------------
% generate a movie of the total field
%--------------------------------------------

% put the movie in a new figure so it doesn't overwrite
% the other figures
figure(4)

% make movie with domain specified by D
m = obj.movie(D);

% export the movie as a gif
m.export('tmatsolver_movie.gif','gif')

%--------------------------------------------
% plot the far field
%--------------------------------------------

% array of observation directions at which to compute
% the far field
tp = linspace(0,2*pi,1000);

% compute the far field
f = obj.getFarfield(exp(1i*tp));

% plot the far field in figure(5)
figure(5)
polarplot(tp,abs(f),'r-')

% add a title to the figure
title('far field')

%--------------------------------------------
% measure the field at the sensor
%--------------------------------------------

% here obj.soln is a cell array holding the total field...
% we know the incident field is zero at the sensor so we
% don't need to include it
measurement = s.measure(inc,obj.soln{:});

% print out the value
fprintf('Measured field is %0.2f\n',measurement);