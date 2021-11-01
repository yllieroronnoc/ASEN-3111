%{
Author: Connor T. O'Reilly
Collaborators: kevin yevak       
Last Revision: 12/12/2019
 %}


%% House Keeping
clc; clear all; close all;

%% Data 

%obtained from ansys fluent at multiple AOA
cfd_data = readtable('CFD_Data_Real.csv');
%delete NAN values from empty rows
cfd_data = cfd_data(~any(ismissing(cfd_data),2),:);

%seperate dat
AOA = cfd_data{:,1};
cd_ansys = cfd_data{:,2};
cl_ansys = cfd_data{:,3};
AOA_rad = cfd_data{:,4};
cos_AOA = cfd_data{:,5};
sin_AOA = cfd_data{:,6};

%Import NASA Experimental Data for NACA 0012 Airfoil
    %Digitized data from Abbott & Von Doenhoff
    
%first column: AOA (deg) Second column: CL
 NASA_Exp_Cl = dlmread('NACA0012cl.dat');
 
%first column: CL Second column: CD
NASA_Exp_ClvCd = dlmread('NASA_Exp_ClCd.dat');
 
%% Problem 1 estimations

%plot cd and cl values estimated using ansys

%cl
figure('name', 'CL ansys')
hold on
%fit line to cfd
AOA_cfd_vals = linspace(0,16,10000);
p = fit(AOA,cl_ansys,'poly5');
cl_cfd_fit = feval(p,AOA_cfd_vals);
%actual plotting
plot(AOA_cfd_vals,cl_cfd_fit)
scatter(AOA,cl_ansys)
title('Sectional Coefficient of Lift vs AOA (RANS model)')
xlabel('Angle of Attack (deg)')
ylabel('C_L')
hold off

figure('name', 'Cd ansys')
hold on 
%actual plotting
plot(cl_ansys,cd_ansys)
scatter(cl_ansys,cd_ansys)
title('C_d vs C_L (RANS model)')
xlabel('C_L')
ylabel('C_d')
hold off

%determine following quantities
    %Lift slope
    %Zero-lift angle of attack
    %Stall angle
    %Maximum sectional coefficient of lift
    
%Lift slope
    %determined from linear portion of lift slope
    %from AOA = 2 deg to AOA = 10 deg
    
p = polyfit(AOA(3:8),cl_ansys(3:8),1);
a_cfd = p(1); %lift slope
a_cfd_rad = p(1) * (180/pi);

%Zero-lift angle of attack
    %when CL is closest to zero
    %use both curve fitted line and cfd data

%take abs of array, grab minimum 
val = min(abs(cl_ansys));
val2 = min(abs(cl_cfd_fit));

if val < val2 
    %if sim values are closer to zero
    idx = find(abs(cl_ansys) == val);
    zero_L_cfd = AOA(idx);
else
    %if curve fitted values are closer
    idx = find(abs(cl_cfd_fit) == val2);
    zero_L_cfd = AOA_cfd_vals(idx);
end

%stall angle and maximum sectional lit coefficient
    %stall angle is at max c_l
    %find max of both fitted and actual values

val = max(cl_ansys);
val2 = max(cl_cfd_fit);
if val > val2
    max_cl_cfd = val;
    idx = find(cl_ansys == val);
    stall_cfd = AOA(idx);
else
    max_cl_cfd = val2;
    idx = find(cl_cfd_fit == val2);
    stall_cfd = AOA_cfd_vals(idx);
end

%weird error correction
if length(stall_cfd) ~= 1
    stall_cfd = max(stall_cfd);
end

    
%% Experimental Data and Numerical Methods

%lift slope for thin airfoil is 2pi

% Plot experimental data obtained from NASA 0012 exp

%cl vs AOA
figure('name', 'NASA Exp CL AOA' )
hold on
%fit line to data
AOA_fitted = linspace(-18,18,10000);
p = fit( NASA_Exp_Cl(:,1), NASA_Exp_Cl(:,2),'poly5');
a_exp = p(1);
a_exp_rad = p(1) * (180/pi);
%smoothing spline is used to find maximum cl values, other curves did not
%work too well
clAOA_exp_fit = feval(p,AOA_fitted);
%actual plotting 
plot(AOA_fitted,clAOA_exp_fit)
scatter(NASA_Exp_Cl(:,1), NASA_Exp_Cl(:,2))
title('C_L vs AOA (Experimental)')
xlabel('AOA (deg)')
ylabel('C_L')
hold off

%characteristics
%zero lift AOA
zero_L_exp = 0;
%stall angle exp
val = max(NASA_Exp_Cl(:,2));
val2 = max(clAOA_exp_fit);
if val > val2
    max_cl_exp = val;
    idx = find(NASA_Exp_Cl(:,2) == val);
    stall_exp =  NASA_Exp_Cl(idx,1);
else
    max_cl_exp = val2;
    idx = find(clAOA_exp_fit == val2);
    stall_exp = AOA_fitted(idx);
end

%weird error correction
if length(stall_exp) ~= 1
    stall_exp = max(stall_exp);
end

%Cd vs CL
figure('name', 'NASA Exp CL CD' )
hold on
%fit line to data
cd_fitted = linspace(-1.5,1.5,10000);
p = fit( NASA_Exp_ClvCd(:,1), NASA_Exp_ClvCd(:,2),'poly2');
clcd_exp_fit = feval(p,cd_fitted);
%actual plotting 
plot(cd_fitted,clcd_exp_fit)
scatter(NASA_Exp_ClvCd(:,1), NASA_Exp_ClvCd(:,2))
title('C_d vs C_L (Experimental)')
xlabel('C_L')
ylabel('C_D')
hold off

%plot for thinairfoil theory

%characteristics of 
a_thin = 2*pi;
AOA_thin = rot90(0:17);
cl_thin = AOA_thin * 2 * (pi^2/180) ;
max_cl_thin = max(cl_thin);

%thin airfoil theory
figure('name', 'thin airfoil theory' )
hold on
%line fitting 
AOA_thin_vals = linspace(0,17,10000);
p = fit(AOA_thin, cl_thin,'poly1');
cl_thin_fit = feval(p,AOA_thin_vals);
%actual plotting
plot(AOA_thin_vals,cl_thin_fit)
scatter(AOA_thin,cl_thin)
title('C_L for AOA (Thin Airfoil Theory)')
xlabel('Angle of attack (deg)')
ylabel('C_L')
hold off

%cahracteristics of vortex panel method
a_vortex = 6.89;
zero_L_vortex = 0.004;
AOA_vortex = rot90(0:17);
cl_vortex = (AOA_vortex + 0.004) * a_vortex * (pi/180);
max_cl_vortex = max(cl_vortex);

%vortex panel method
figure('name', 'vortex panel method' )
hold on
%line fitting 
AOA_vor_vals = linspace(0,17,10000);
p = fit(AOA_vortex, cl_vortex,'poly1');
cl_vor_fit = feval(p,AOA_vor_vals);
%actual plotting
plot(AOA_vor_vals,cl_vor_fit)
scatter(AOA_vortex,cl_vortex)
title('C_L vs AOA (Vortex Panel)')
xlabel('Angle of attack (deg)')
ylabel('C_L')
hold off


%% comparison plots

%plot lines of best fit
figure('name', 'comparison' )
hold on
%line fitting 
plot(AOA_cfd_vals,cl_cfd_fit)
plot(AOA_fitted,clAOA_exp_fit)
plot(AOA_thin_vals,cl_thin_fit)
plot(AOA_vor_vals,cl_vor_fit)
title('Comparison of Lift slopes')
xlabel('Angle of attack (deg)')
ylabel('C_L')
legend('CFD','Experimental','Thin Airfoil','Vortex Menthod','Location','southeast')
grid on 
hold off

%% calculation outputs  
fprintf('RANS Model: \n');
fprintf('\t    Lift Slope: %0.3f (/rad) \n', a_cfd_rad )
fprintf('\t Zero-Lift AOA: %0.3f (degrees)\n', zero_L_cfd)
fprintf('\t   Stall angle: %0.3f (degrees)\n', stall_cfd)
fprintf('\t   Maximum C_l: %0.3f \n', max_cl_cfd)

%% comparison outputs

%thin airfoil 
fprintf('Thin Airfoil Theory: \n');
fprintf('\t    Lift Slope: 2pi (/rad) or %0.3f \n',a_thin)
fprintf('\t Zero-Lift AOA: 0 (degrees)\n')
fprintf('\t   Stall angle: DNE, thin airfoil doesnt\n \t \t \t \t \t \t account for stall \n')
fprintf('\t   Maximum C_l: %0.3f \n', max_cl_thin)

%vortex panel
fprintf('Vortex Panel Method: \n');
fprintf('\t    Lift Slope: %0.4f (/rad) \n', a_vortex )
fprintf('\t Zero-Lift AOA: %0.4f (degrees)\n', zero_L_vortex)
fprintf('\t   Stall angle: DNE \n')
fprintf('\t   Maximum C_l: %0.4f \n', max_cl_vortex)

%experimental panel
fprintf('Experimental: \n');
fprintf('\t    Lift Slope: %0.4f (/rad) \n', a_exp_rad )
fprintf('\t Zero-Lift AOA: %0.4f (degrees)\n', zero_L_exp)
fprintf('\t   Stall angle: %0.4f \n', stall_exp)
fprintf('\t   Maximum C_l: %0.4f \n', max_cl_exp)

%% Error analysis

%Lift slope
rel_a_rans = abs((a_exp_rad - a_cfd_rad)/a_exp_rad) * 100;
rel_a_thin = abs((a_exp_rad - a_thin)/a_exp_rad) * 100;
rel_a_vor = abs((a_exp_rad - a_vortex)/a_exp_rad) * 100;

rel_cl_rans = abs((max_cl_exp - max_cl_cfd)/max_cl_exp) * 100;
rel_cl_thin = abs((max_cl_exp - max_cl_thin)/max_cl_exp) * 100;
rel_cl_vor = abs((max_cl_exp -  max_cl_vortex)/max_cl_exp) * 100;

rel_st_rans = abs((stall_exp -  stall_cfd)/stall_exp) * 100;

fprintf('Error Analysis (values are precentage)')
fprintf('\n \n Lift Slope: \n')
fprintf('\t RANS: %0.3f \n',rel_a_rans)
fprintf('\t Thin Airfoil: %0.3f \n',rel_a_thin)
fprintf('\t Vortex: %0.3f \n',rel_a_vor)

fprintf('\n Zero Lift AOA: \n')
fprintf('\t RANS: %0.3f \n',0)
fprintf('\t Thin Airfoil: %0.3f \n',0)
fprintf('\t Vortex: %0.3f \n', zero_L_vortex * 100)


fprintf('\n Max CL: \n')
fprintf('\t RANS: %0.3f \n',rel_cl_rans)
fprintf('\t Thin Airfoil: %0.3f \n',rel_cl_thin)
fprintf('\t Vortex: %0.3f \n', rel_cl_vor)
fprintf('Thin airfoil and vortex do not account for stall so do not take zeriously')


fprintf('\n Stall Angle: \n')
fprintf('\t RANS: %0.3f \n',rel_st_rans)



    


