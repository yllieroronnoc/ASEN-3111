%house keeping 
clc; clear all; close all;
%{
Author: Connor T. O'Reilly
Collaborators: ryan hughes
               kevin yevak
               samuel felice
                
Last Revision: 11/7/2019
 %}
%% variables
chord = 1; %m
vinf = 50; %m/s

%% Question 2
%flow for couple different resolutions
figure('name','cl v N')
k = 1;
for N = 10:2:200
    AOA = 0; %deg
    plot_cp = false;
    %for NACA 0012
    [x,y] = NACA_Airfoils(0,0,12,chord,N);
    c_l(k) = Vortex_Panel(x,y,vinf,AOA,plot_cp);
    scatter(N,c_l(k),'*')
    hold on  
    k = k+1;
end
xlabel('Number of panels used')
ylabel('c_l')
title('c_l calculated vs number of pannels used')
plot([0, 200], [c_l(end),c_l(end)])
%find optimal number of panels to use
fprintf('cl values seem to decrease and level off\n chosing a point relatively close to last cl value to decrease run time')

N = 10:2:200;
mid = floor(length(N)/2);
N = N(mid);
err = ((c_l(mid)-c_l(end))/c_l(end)) * 100;
fprintf('\n optimal number of panels to use is: %i',N)
fprintf('\n %f%% greater when compared to cl at much higher N: %f \n',err-100)

%c_l and c_p for different aoa

AOA = [-5,0,5,10]; %degrees
k = 1;
for i = 1:length(AOA)
     plot_cp = true;
    %for NACA 0012
    [x,y] = NACA_Airfoils(0,0,12,chord,N);
    c_l(k) = Vortex_Panel(x,y,vinf,AOA(i),plot_cp);
    k = k + 1;
end

%% Question 3

%obtain plots of sectional coefficient of lift vs AOA for following
%airfoils

%m,p,and t of airfoils
NACAS = [5,0,12;5,10,12;5,20,12;5,30,12];
AOA = [-5,0,5,10]; %degrees
figure('name','NACA c_l')
w = 1;
plot_cp = false;
for i = 1:4
    k = 1;
    for j = 1:length(AOA)
        [x,y] = NACA_Airfoils(NACAS(i,1),NACAS(i,2),NACAS(i,3),chord,100);
        c_l3(j) = Vortex_Panel(x,y,vinf,AOA(j),plot_cp);
        cl_other(w) = c_l3(j); 
        w = w + 1;
    end
    plot(AOA,c_l3)
    hold on
end
xlabel('AOA in degrees')
ylabel('c_l')
title('cl values for different AOA of different NACA airfoils')
legend('NACA 0012','NACA 2412','NACA 4412','NACA  2424','Location','southeast')

%calculate lift slope

%for NACA 0012, symmetric
p = polyfit(AOA,cl_other(1:length(AOA)) ,1);
y1 = polyval(p,AOA);
fprintf('\n \n estimated lift slope of NACA 0012: %f',p(1))

%for NACA 2412, cambered
p = polyfit(AOA,cl_other(length(AOA)+1:length(AOA)+4) ,1);
y1 = polyval(p,AOA);
fprintf('\n estimated lift slope of NACA 2412: %f',p(1))

%for NACA 4412,  cambered
p = polyfit(AOA,cl_other(length(AOA)+5:length(AOA)+8) ,1);
y1 = polyval(p,AOA);
fprintf('\n estimated lift slope of NACA 4412: %f',p(1))

%for NACA 2424,  cambered
p = polyfit(AOA,cl_other(length(AOA)+9:length(AOA)+12) ,1);
y1 = polyval(p,AOA);
fprintf('\n estimated lift slope of NACA 2424: %f \n',p(1))

fprintf('for thin airfoil theory for symmetric airfoil, lift slope = 2*pi')

%for NACA 0012, symmetric
p = polyfit(linspace(-5,10,1000),linspace(cl_other(1),cl_other(length(AOA)),1000) ,1);
y1 = polyval(p,linspace(-5,10,1000));
yeet = linspace(-5,10,1000);
fprintf('\n estimated zero lift AOA for NACA 0012 is: %f',yeet(289))
fprintf('\n for NACA 0012, lift slope varies drastically from 2*pi but zero lift \n AOA is close to 0 deg')