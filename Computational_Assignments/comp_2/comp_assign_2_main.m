%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Computational Assignment 2      %%%
%%%  Connor T. O'Reilly              %%%
%%%  email: coor1752@colorado.edu    %%%
%%%  SID: 107054811                  %%%
%%%  Collaborated with: Ethan Fleer  %%%
%%%                     Kevin Yevak  %%%
%%%                     Adam Elsayed %%%
%%%                     
%%%  Date: 09/19/19                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% House keeping 
clc; clear all; close all;
tic
%% Visualization
% Given
% symetric airfoil 
chord = 2; %m
aoa = 0.2; 
v_inf = 68; %m/s
p_inf = 101.3 * 10^3; %Pa
rho_inf = 1.225; %kg/m^3
N = 1200;

Plot_Airfoil_Flow(chord,aoa,v_inf,p_inf,rho_inf,N);

timeElapsed = toc;
fprintf('Total Time elapsed %.2f \n as the Flow Stream Velocity increases, stream lines do not apear to change. \n as the AOA increases, stream lines become more circular. \n as the chord length increases, do not apear to change.  ',timeElapsed)
print('as the Flow Stream Velocity increases, stream lines do not apear to change. ')
print('as the AOA increases, stream lines become more circular. ')
print('as the chord length increases, do not apear to change. ')
function [] = Plot_Airfoil_Flow(c,alpha,V_inf,p_inf,rho_inf,N)
%function will plot stream lines, potential lines, and pressure contours
%using the given inputs.
%define change in x
deltx = c/N;

%define x values and set all yvals along airfoil to zero
xvals = linspace(deltx/2,c-deltx,1200);
yvals = zeros(1,length(xvals));

%define vortex strength
gamma = 2 * alpha * V_inf * sqrt((1-(xvals/c))./(xvals/c));

%circulation
circ = gamma * deltx;

%% define domain

xmin = -2*c;
xmax = 2 * c;
ymin = -c;
ymax = c;

%% Define Number of Grid Points
nx = 100; % steps in the x direction
ny = 100; % steps in the y direction
%% Create mesh over domain using number of grid points specified
[x,y]=meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));
%% Define Flow Parameters
r =@(x1,y1) sqrt((x-x1).^2 + (y-y1).^2);
%% Create Stream Function 
%stream function of 
Stream_Uniform = V_inf * (y * cos(alpha) - x * sin(alpha));
%Stream_VOrtex = (VorStrength(1)*log(r)/(2pi)
%% Define a function which calculates the radius.
Stream_Vortex = 0;
for i = 1:length(xvals)
    Stream_Vortex = Stream_Vortex + (circ(i) * log(r(xvals(i),yvals(i))))/(2*pi);
end
StreamFunction = Stream_Uniform + Stream_Vortex;

%%determin equipotential function where stream functions are equal
%% Determine color levels for contours
levmin = StreamFunction(1,nx); % defines the color levels -> trial and error to find a good representation
levmax = StreamFunction(ny,nx/2);
levels = linspace(levmin,levmax,100)';
%find velocity potential
Pot_Vortex = 0;

levmin2 = StreamFunction(1,nx); % defines the color levels -> trial and error to find a good representation
levmax2 = StreamFunction(ny,nx/50);
levels2 = linspace(levmin2,levmax2,100)';
for i = 1:length(xvals)
    theta = atan2(-y,-x+xvals(i));
    Pot_Vortex = Pot_Vortex + (-circ(i)/(2*pi)*theta);
end
Pot_Uniform = V_inf * (x * cos(alpha) + y * sin(alpha));
Pot = Pot_Uniform + Pot_Vortex;

%find cp
q_inf = .5 * rho_inf * V_inf^2;
Cp = 1-(gradient(Pot)./V_inf).^2;
P = q_inf * Cp + p_inf;
figure(1)
contour(x,y,StreamFunction,levels)
hold on 
plot([0 c],[0 0],'linewidth',2)
axis equal
axis([xmin xmax ymin ymax])
ylabel('y')
xlabel('x')
title('Streamlines with N = 1200 discrete vortices')
figure(2)
contour(x,y,Pot,levels2)
hold on 
plot([0 c],[0 0],'linewidth',2)
axis equal
axis([xmin xmax ymin ymax])
ylabel('y')
xlabel('x')
title('Velocity Potential with N = 1200 discrete vortices')

figure(3)
contour(x,y,StreamFunction,levels)
hold on
contour(x,y,Pot,100)
hold on
plot([0 c],[0 0],'linewidth',2)
axis equal
axis([xmin xmax ymin ymax])
title('Stream lines and Velocity Potential with N = 1200 discrete vortices')
ylabel('y')
xlabel('x')

figure(4)
contourf(x,y,P,100)
hold on
plot([0 c],[0 0],'linewidth',2)
axis equal
title('Pressure Contours with N= 1200')
axis([xmin xmax ymin ymax])
ylabel('y')
xlabel('x')
end
