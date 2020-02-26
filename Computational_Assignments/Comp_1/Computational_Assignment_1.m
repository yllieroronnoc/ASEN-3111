%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Computational Assignment 1      %%%
%%%  Connor T. O'Reilly              %%%
%%%  email: coor1752@colorado.edu    %%%
%%%  SID: 107054811                  %%%
%%%  Collaborated with: Ethan Fleer  %%%
%%%                     Kevin Yevak  %%%
%%%                     Adam Elsayed %%%
%%%                     
%%%  Date: 09/19/19                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%house keeping
clear all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question #1 %%
%%%%%%%%%%%%%%%%%

fprintf('---- Question 1 ----\r \r')


%Given
diameter = 1; %m
v_inf = 30; %m/s
q_inf  = (.5)*(1.225)*(v_inf^2);% dynamic pressure (kg/m^3)
p_inf = 101.3 * 10^3; %freestream pressure(Pa)

%coefficient of pressure
cp = @(p) (p - p_inf)/q_inf;

%deriving an equation for surface pressure distribution
%using two equations provided for the coefficient of pressure
p = @(t) q_inf.*(1 - (4 .* sin(t).^2))+ p_inf;


%initialize vectors for normal and axial force calcs
Norm = zeros(1,100); Axe = zeros(1,100);
for N = 1:100
    %initialize sums
    sumN = 0; sumA = 0;
    %integral will be taken from a to b
    a = 0; %radians
    b = 2*pi; %radians

    %vector for theta  values
    t = linspace(0,2*pi,N+1); %radians

    %define h
    h = (t(N+1)-t(N))/(2*N);
    %simpsons rule
        for i = 1 : (N/2)
            sumN = sumN + ((p(t(2*i - 1)) * cos(t(2*i - 1))) + ((4*p(t(2*i))*cos(t(2*i)))) + (p(t(2*i + 1)) * cos(t(2*i + 1))));
            sumA = sumA + ((p(t(2*i - 1)) * sin(t(2*i - 1))) + ((4*p(t(2*i))*sin(t(2*i)))) + (p(t(2*i + 1)) * sin(t(2*i + 1))));
        end
        %finish simpsons rule calc
        Norm(N) = sumN * (h/3) * (diameter/2);
        Axe(N) = sumA * (h/3) * (diameter/2);
end

%Normal/Axial Force and Lift/Drag are equivelant because every AOA is identical

figure(59)
plot((1:N),Norm,'Color','r')
hold on
plot((1:N),Axe,'Color','g')
legend('Lift','Drag')
xlabel('Number of panels used')

%Find num panels so that solutions are within 0.001
%Newtons of the exact solutions
for i = 2:N
    if abs(Norm(i)) <= .001
        kN = i;
        break
    end
end

for i = 2:N
    
    if abs(Axe(i)) < .001
        jA = i;
        break
    end
end

%exact solutions are zero because ideal flow and pressure
%is equal and opposite for dif side of cylinder.
fprintf('Lift using %i Panels: %f \r',N,Norm(end))
fprintf('Drag using %i Panels: %f \r',N,Axe(end))
%turns out simpons rule only uses 
fprintf('Lift estimation within 0.001 Newtons after: %i Panels \r',kN)
fprintf('Drag estimation within 0.001 Newtons after: %i Panels \r',jA)
fprintf('--------------------------------------------------------------\r\r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% question 2 %%
%%%%%%%%%%%%%%%%
close all;
t = cputime;

%given 
aoa = 9; %degrees
chrd = 2; %m
v_inf = 30; %m/s
rho_inf = 1.225; %kg/m^3
p_inf = 101.3 *10^3;% Pa

%Determine lift and fraag per unit span on NACA 00122 airfoil
load Cp 
%define number of integration points for ideal cn and ca
intPoints = 1000;

%calculate actual and estimated Lift and Drag
[CnLegit,CnActual,CaLegit,CaActual] = naca0012_n_a(Cp_upper,Cp_lower,chrd,intPoints);

ClLegit = CnLegit*cosd(aoa) - CaLegit*sind(aoa);
CdLegit = CnLegit*sind(aoa) + CaLegit*cosd(aoa);

LLegit = ClLegit * (.5) * rho_inf * (v_inf^2) * chrd ;
DLegit = CdLegit * .5 * rho_inf * (v_inf^2) * chrd ;


ClActual =  CnActual.*cosd(aoa) - CaActual.*sind(aoa);
CdActual = CnActual.*sind(aoa) + CaActual.*cosd(aoa);

LActual = ClActual * .5 * rho_inf * v_inf^2 * chrd ;
DActual = CdActual * .5 * rho_inf * v_inf^2 * chrd ;

%figure out relative error for lift and drag
er5 = 0;
er1 = 0;
er110 = 0;
for i = 1:length(LActual)
   curEr = (LLegit - LActual(i))/LLegit;
   if(curEr <= .5)
       er5 = i;
       break
   end
   curEr = 0;
end
for j = 1:length(LActual)
   curEr = (LLegit - LActual(j))/LLegit;
   if(curEr <= .1)
       er1 = j;
       break
   end
   curEr = 0;
end

for k = 1:length(LActual)
   curEr = (LLegit - LActual(k))/LLegit;
   if(curEr <= 0.001)
       er110 = k;
       break
   end
   curEr = 0;
end
for i = 1:length(DActual)
   curEr = (DLegit - DActual(i))/DLegit;
   if(curEr <= .5)
       er5D = i;
       break
   end
   curEr = 0;
end
for j = 1:length(DActual)
   curEr = (DLegit - DActual(j))/DLegit;
   if(curEr <= .1)
       er1D = j;
       break
   end
   curEr = 0;
end

for k = 1:length(LActual)
   curEr = (DLegit - DActual(k))/DLegit;
   if(curEr <= 0.001)
       er110D = k;
       break
   end
   curEr = 0;
end
e = cputime-t;
fprintf('----- Question 2 ----- \r \r')
fprintf('Computation Time: %f \r \r',e)
fprintf('Lift: %f\r',LLegit)
fprintf('Drag: %f\r \r',DLegit)

figure(60)
plot((1:intPoints),LActual,'-*','MarkerSize',3)
hold on
line([er5 er5],[min(LActual),max(LActual)],'LineWidth',3,'Color','g');
line([er1 er1],[min(LActual),max(LActual)],'LineWidth',3,'Color','r');
line([er110 er110],[min(LActual),max(LActual)],'LineWidth',3,'Color','b');
xlabel('Number of Integration Points Used')
ylabel('Lift')
legend('Lift','Relative Error of 5%','Relative Error of 1%','Relative Error of 0.1%')
fprintf('Number of integration points for Lift to be within 5%% relative error:  %i \r',er5)
fprintf('Number of integration points for Lift to be within 1%% relative error:  %i \r',er1)
fprintf('Number of integration points for Lift to be within 0.1%% relative error: %i \r \r',er110)

figure(61)
plot((1:intPoints),DActual,'-*','MarkerSize',3)
hold on
line([er5D er5D],[min(DActual),max(DActual)],'LineWidth',2,'Color','g');
line([er1D er1D],[min(DActual),max(DActual)],'LineWidth',2,'Color','r');
line([er110D er110D],[min(DActual),max(DActual)],'LineWidth',3,'Color','b');
xlabel('Number of Integration Points Used')
ylabel('Drag')
legend('Drag','Relative Error of 5%','Relative Error of 1%','Relative Error of 0.1%')
fprintf('Number of integration points for Drag to be within 5%% relative error:  %i \r',er5D)
fprintf('Number of integration points for Drag to be within 1%% relative error:  %i \r',er1D)
fprintf('Number of integration points for Drag to be within 0.1%% relative error: %i \r',er110D)

fprintf('-------------------------------------------------------------- \r')

%% Function to estimate normal force and axial force coefficients
function [CnLegit,CnActual,CaLegit,CaActual] = naca0012_n_a(cp_upper,cp_lower,c,intPoints)
%Function will calcuate the value of Cn and Ca for a NACA 0012 airfoil
%using trapezodial approximation with 10000000 panels as well as estimating
%Cn and Ca using a set amount of integration points.

% Author: Connor T. O'Reilly
% Collaborators: Davis Peirce, Kevin Yevak
% Date: 09/19/19

%formula for shape of NACA 0012 airfoil
%distance from chord to upper surface
yu = @(x) (.12/0.2).*2 .* ( 0.2969.*sqrt(x./2) -0.1260*(x./2) - 0.3516.*(x./2).^2 + 0.2843.*(x./2).^3 - 0.1036.*(x./2).^4);
%distance from chord to lower surface
%same but negative because airfoil is sym
yl = @(x) -1.*((.12/0.2).*2 .* ( 0.2969.*sqrt(x./2) -0.1260*(x./2) - 0.3516.*(x./2).^2 + 0.2843.*(x./2).^3 - 0.1036.*(x./2).^4));

%find integration points for actual Cn and Ca calculation
xL = linspace(0,c,10000000);

%calculate Cp along upper and lower surfaces
Cpu = fnval(cp_upper,xL/c);
Cpl = fnval(cp_lower,xL/c);

%use trapz to get actual coeffiecients of normal and axial force
CnLegit = (1/c)*(trapz(xL,Cpl) - trapz(xL,Cpu));
CaLegit = (1/c) * (trapz(yu(xL),Cpu) - trapz(yl(xL),Cpl));

%use trapezodial rule to determine Cn and Ca estimate.
CaActual = zeros(1,intPoints);
CnActual = zeros(1,intPoints);
for i = 1:intPoints
    %loop will change the number of integration points for each itteration 
    
    %find the normal force coeff by using equation 1.15 in anderson 
    %assuming friction coeff is negligible 
    xN = linspace(0,c,i);
    x1 = xN(1:i-1);
    x2 = xN(2:i);
    
    %determine infinitesimal distance dx
    dx = x2 - x1;
    
    %find cp values for upper surface
    Cpu1 = fnval(cp_upper,x1/c);
    Cpu2 = fnval(cp_upper,x2/c);
    
    %determine Cp along upper surface by averaging the pressure
    %between pressure points and multiplying by dx
    CpU = (Cpu2+Cpu1)/2  .* dx;
    
    
    Cpl1 = fnval(cp_lower,x1/c);
    Cpl2 = fnval(cp_lower,x2/c);

    CpL = (Cpl2+Cpl1)/2  .* dx;
    
    %find dif between lower and upper surface to determine overall Cp
    %sum and mult by (1/c) to calculate Cn
    
    %store normal coef for num of interation points in vec
    CnActual(i) = (1/c)*sum(CpL-CpU);

    %find the axial force coeff by using equation 1.16 in anderson 
    %assuming friction coeff is negligible 
    
    %find change in y value for upper and lower surface
    y1 = yu(x1);
    y2 = yu(x2);
    y3 = yl(x1);
    y4 = yl(x2);
    dyu = y2-y1;
    dyl = (y4-y3);

    %calculate axial force coeff for upper and lower surface
    CaU = (Cpu2+Cpu1)/2 .* dyu;
    CaL = (Cpl1 + Cpl2)/2 .* dyl;
    
    %store axial coef for num of interation points in vec
    CaActual(i) = (1/c)*sum(CaU - CaL);
end

end



