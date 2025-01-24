%% Clear previous variables 
clear

%  Add path to local version of SeaFreeze (edit for your personal
%  enviornment)
addpath '/Users/jacobbuffo/Documents/MATLAB/SeaFreeze-master/Matlab'


%% Declare variables and their values
%Declare the initial values and start a timer

% Note to Jacob - clean this section up for future versions!!!

tic
global start_height
global resolution
global final_height

minP=0; % Minimum pressure accepted (MPa)
maxP=2300; % Maximum pressure accepted (MPa)

T_surf=235; % Surface Temp of planet (in K)
Base_Flux=30/1000;  % Heat Flux (in W/m^2)
start_height=0;  % Base of domain (km)
resolution=0.5;  % Domain resolution (km)
final_height=30; % Top of domain (km)
g=24; % Gravity of planet (in m/s^2)
Patmosphere=0.101325; % Atmospheric pressure of planet (in MPa)
Mx = 6.38; % Mass factor of the planet (Earth Mass)   (LHS 1140 - 6.38)
Rx = 1.64; % Radius factor of the planet (Earth Radius) (LHS 1140 - 1.64)
%% Create an initial adiabatic profile for the ice sheet and depth dependent 
% melting temperature using the two subroutines below
[T_Start,dtdz,Pressure,rho,alpha,Cp,K,phase,g_s] = adiabat_profile(Patmosphere,T_surf,(final_height-start_height)*1000,resolution,Mx,Rx);
step=(maxP-minP+1)/numel(K);
[MeltT] = findmeltT2(Pressure);
MeltT = MeltT';
Height_list=linspace(start_height,final_height,numel(T_Start));

%% Plot the initial profiles
%Do you want to plot the intial profiles of parameters?%
%plotpaam = 1 plots the relevant initial profiles
plotparam = 1;
if plotparam == 1
figure('units','normalized','position',[.1 .1 .5 .6])
subplot(2,2,1)
plot(Height_list,T_Start,'k','LineWidth',1.9)
xlabel('Depth (km)')
ylabel('Temperature (K)')
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
xlim([0 final_height])
subplot(2,2,2)
plot(Height_list,rho,'k','LineWidth',1.9)
xlabel('Depth (km)')
ylabel('Density (kg m^{-3})')
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
xlim([0 final_height])
subplot(2,2,3)
plot(Height_list,Cp,'k','LineWidth',1.9)
xlabel('Depth (km)')
ylabel('Cp (J kg^{-1} K^{-1})')
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
xlim([0 final_height])
subplot(2,2,4)
scatter(Height_list,phase,'k','LineWidth',1.9)
xlabel('Depth (km)')
% ylabel('K (W m^{-1} K^{-1})')
ylabel('Phase')
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
box on
xlim([0 final_height])
end

%% The main thermal evolution code starts here
timestep = 10000;                   % Time step in years
dt=86400*365.25*timestep;           % Time step in seconds
tf=50*dt;                           % Final time (s)
k_i=K';                             % Thermal Conductivites of Ices
rho_i=rho;                          % Densities of Ices
rho_w=1000;                         % Density of pur water
c_i=Cp;                             % Specific Heats of Ices
Tm=MeltT;                           % Melting Temps
z = linspace(start_height, final_height, numel(T_Start));
dz=1000*resolution;                 % Resolution in (m)
L=334778;                           % Latent heat of fusion (pure ice)
TTol=0.001;                         % Temp Tolerance
Ra_c=1000;                          % Critical Rayleight number for onset of convection
eta_0_Ih=5*10^13;                   % Reference viscosity of ice (Pa*s) for Arrhenius law (10^13 low viscosity, range of 10^13 to 10^15)
eta_0_II=10^15;                     % Reference viscosity of ice II - Echelmeyer 1986
eta_0_III=10^12;                    % Reference viscosity of ice III - Echelmeyer 1986
eta_0_V=10^16;                      % Reference viscosity of ice V - Sotin & Poirier 1987
eta_0_VI=10^14;                     % Reference viscosity of ice VI - Poirier et al., 1981
eta_0_w=1.8*10^-3;                  % Reference viscosity of water (Pa*s) for Arrhenius law 
A=7.5;                              % Coefficient for Arrhenius law viscosity (see Mitri & Showman 2005 or Cuffery and Patterson 2010)
m_dot=Base_Flux/(rho_w*L);          % Mass flux of basal meltwater through HP ices - Kalousova & Sotin 2018

% plotting mats
hold_temps=[];
hold_phase=[];
hold_conv_top=zeros(1000*final_height/dz,tf/dt);
hold_conv_bottom=zeros(1000*final_height/dz,tf/dt);
hold_delta_thick_top=zeros(1000*final_height/dz,tf/dt);
hold_delta_thick_bottom=zeros(1000*final_height/dz,tf/dt);
hold_Ra=[];
hold_is_convect=zeros(1000*final_height/dz,tf/dt);
layers=[];   % number of different ice/water layers for indexing

% loop over time
%Work on making an ice phase plot showing progressive change with depth and time
PT = [Pressure T_Start];
phasenew = SF_WhichPhase(PT);
test(:,1) = T_Start';
test2(:,1) = phasenew';
skip = 2;
ss = 0;
m_H2O = 0;              % Mass of water flowing through HP ices

for time=[0:dt:tf]
    % melt transport dynamics
    if time>0                           % no transport before first time step
        if phasenew(end)==0 && rho_i(top_hold(end)-1)>rho_w
            m_H2O=m_H2O+m_dot*dt;
        else
        end
    
        if m_H2O>dz && rho_i(top_hold(end)-1)>rho_w
            for i=1:length(phasenew)
                if rho_i(i)<rho_w
                    insert_oc=i+1;
                else
                    break
                end
            end
            T_hold=T_Start;
            T_hold(insert_oc:insert_oc-1+floor(m_H2O/dz))=T_Start(top_hold(end):top_hold(end)-1+floor(m_H2O/dz));
            T_hold(insert_oc+floor(m_H2O/dz):top_hold(end)-1+floor(m_H2O/dz))=T_Start(insert_oc:top_hold(end)-1);
            T_Start=T_hold;
            m_H2O=m_H2O-floor(m_H2O/dz)*dz;
        else
        end
    else
    end
    % end
    
    PT = [Pressure T_Start];
    phasenew = SF_WhichPhase(PT);
    
    ss = ss+1;
    [k_i,rho_i,c_i,alpha] = compute_params_v2(PT,phasenew);
    hold_phase=[hold_phase phasenew];
    hold_temps=[hold_temps T_Start];
        
    % viscosity vector
    for j=1:length(phasenew)
        if phasenew(j)==0
            eta_0(j)=eta_0_w;
        elseif phasenew(j)==1
            eta_0(j)=eta_0_Ih;
        elseif phasenew(j)==2
            eta_0(j)=eta_0_II;
        elseif phasenew(j)==3
            eta_0(j)=eta_0_III;
        elseif phasenew(j)==5
            eta_0(j)=eta_0_V;
        elseif phasenew(j)==6
            eta_0(j)=eta_0_VI;
        else
            eta_0(j)=eta_0_Ih;
        end
    end
    % end
    
    % Conduction code
    [T_new,Conv_top,Conv_bottom,Ra,delta_thick_top,delta_thick_bottom,is_convect,k_bar,top_hold,bottom_hold]...
        =HP_Ice_Evolve_v9(T_Start,k_i,rho_i,c_i,dt,dz,T_surf,Base_Flux,Tm,TTol,Pressure,phasenew,Ra_c,eta_0,alpha,g_s,A,Mx,Rx);
    
    T_Start=T_new;
    T_Start(1)=T_surf;
    PT = [Pressure T_Start];
    phasenew = SF_WhichPhase(PT);
    test(:,skip)  = T_Start';
    test2(:,skip)  = phasenew';
    phasetest(:,skip) = phasenew;
    ktest(:,skip) = k_bar;
    rhotest(:,skip) = rho_i;
    skip = skip+1;
    
    if length(Conv_top)>0
        hold_conv_top(1:length(Conv_top),ss)=Conv_top';
        hold_conv_bottom(1:length(Conv_bottom),ss)=Conv_bottom';
        hold_delta_thick_top(1:length(delta_thick_top),ss)=delta_thick_top';
        hold_delta_thick_bottom(1:length(delta_thick_bottom),ss)=delta_thick_bottom';
    else
    end
    hold_Ra=[hold_Ra Ra'];
    hold_is_convect(1:length(is_convect),ss)=is_convect';
    layers=[layers length(is_convect)];
end

log_Ra=log(hold_Ra);
%%
figure('units','normalized','position',[.1 .1 .3 .6])
ax(1) = subplot(2,1,1);
imagesc([0:dt/(86400*365.25):(time)/(86400*365.25)]./1E6,[Height_list],([test]))
hh = colorbar;
%caxis([T_surf max(T_new)])
ylabel(hh,'Temperature (K)','FontSize',24)
xlabel('Time (Myr)')
ylabel('Depth (km)')
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
colormap hot
% title("T_s = " + Ts + "K, " +"q_s = " + q + "mW m^{-2}") %---------_Change this value
ax(2) = subplot(2,1,2);
imagesc([0:dt/(86400*365.25):(time)/(86400*365.25)]./1E6,[Height_list],([test2]))
cb=colorbar;
xlabel('Time (Myr)')
ylabel('Depth (km)')
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
ylabel(cb,'Ice Phase','FontSize',24)
colormap parula

figure
subplot(2,2,1)
plot(Height_list,T_Start,'k','LineWidth',1.9)
xlabel('Depth (km)')
ylabel('Temperature (K)')
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
xlim([0 final_height])
subplot(2,2,2)
scatter(Height_list,phasenew,'k','LineWidth',1.9)
xlabel('Depth (km)')
% ylabel('K (W m^{-1} K^{-1})')
ylabel('Phase')
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
box on
xlim([0 final_height])
subplot(2,2,3)
hold on
for i=1:time/dt
    k=find(hold_conv_top(:,i));
%     l=find(hold_delta_thick(:,i));
%     if test2(end,i)==0
%         l=[l;l(end)+1];
%         l=[l;l(end)+1];
%     else
%     end
    if isempty(k)
        continue
    else
        scatter(i*dt/(86400*365.25)/1E6,-hold_conv_top(1:k(end),i).*resolution,'b')
        scatter(i*dt/(86400*365.25)/1E6,-hold_conv_bottom(1:k(end),i).*resolution,'r')
        count=0;
        for j=1:layers(i)
            if hold_is_convect(j,i)==1
                count=count+1;
                if hold_delta_thick_top(j,i)+hold_delta_thick_bottom(j,i)<hold_conv_bottom(count,i)-hold_conv_top(count,i)
                    scatter(i*dt/(86400*365.25)/1E6,(-hold_conv_top(count,i)-hold_delta_thick_top(j,i)).*resolution,'k.')
                    scatter(i*dt/(86400*365.25)/1E6,(-hold_conv_bottom(count,i)+hold_delta_thick_bottom(j,i)).*resolution,'k.')
                else
                end
            else
            end
        end            
%         for j=1:length(k)
%             if 2*hold_delta_thick(l(j),i)<hold_conv_bottom(k(j),i)-hold_conv_top(k(j),i)
%                 scatter(i*dt/(86400*365.25)/1E6,(-hold_conv_top(k(j),i)-hold_delta_thick(l(j),i)).*resolution,'k.')
%                 scatter(i*dt/(86400*365.25)/1E6,(-hold_conv_bottom(k(j),i)+hold_delta_thick(l(j),i)).*resolution,'k.')
%             else
%             end
%         end
    end
end
xlabel('Time (Myr)')
ylabel('Depth (km)')
title('Convection Locations')
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
subplot(2,2,4)
imagesc([0:dt/(86400*365.25):(time)/(86400*365.25)]./1E6,[Height_list],([log_Ra]))
cb=colorbar;
xlabel('Time (Myr)')
ylabel('Depth (km)')
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
ylabel(cb,'Log(Ra #)','FontSize',24)
% cmap = [0 0 0; lbmap(1,'RedBlue')];
% colormap(ax(2),flipud(cmap))
%  cb.Ticks = linspace(0, 1, 2) ; %Create 8 ticks from zero to 1
%  cb.TickLabels = num2cell(0:1) ;

toc

