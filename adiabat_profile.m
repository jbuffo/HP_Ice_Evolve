function [T,dT_dz,P,rho,alpha,Cp,K,phase,g_s,z] = adiabat_profile(P_s,T_s,depth,res,Mx,Rx)
%Inputs
% P_s = 0.1;   % Starting Pressure (MPa)
% T_s = 235;   % Starting Temperature (K);
% depth = 80e3;  % thickness of the hydrosphere
% r_b = 6370*1e3;  %Radius rocky core
% rho_core = 5500;    % density of the core (kg/m3);
Gcnst = 6.67E-11; %Gravitational Constant
mass_Earth = 5.97E24;
radius_Earth = 6.371E6;
%Resolution / iterations
res = res*1000;  % depth grid
g_it = 3;  % Gravity convergence iteration (3 is enough)

% %Plots and table
% plot_flg = 1 % plot flag
% tab_flg =1 % Table flag

%Rocky/metal core calculation
% Mass_core = 4/3*pi*r_b^3*rho_core;
Mass = Mx*mass_Earth;
radius = Rx*radius_Earth;
g_s = 6.67430e-11*Mass/radius^2; % Gravity at the Hydrosphere Mantle boundary

    z = linspace(res, depth, depth/res);  % depth grid
    T = zeros(length(z),1);  % Temperature grid
    P = zeros(length(z),1);  % Pressure grid
    rho = zeros(length(z),1);  % Density grid
    alpha = zeros(length(z),1);  % thermal epansivity grid
    K = zeros(length(z),1);
    Cp = zeros(length(z),1);  % Heat capacity grid
    dT_dz = zeros(length(z),1);  % thermal gradient grid
    phase = zeros(length(z),1);  % phase grid
    grav = zeros(length(z),1); % gravity grid
    M_L = zeros(length(z),1); % Mass grid
    grav(:) =g_s; % Constant gravity to start with

%      for k =1:g_it;
        g = flipud(grav); %Unnecessary really since we are keeping g const
        PT = [];
        PT = [P_s T_s];
        phase_s = SF_WhichPhase(PT);
        [out,K] = compute_params(PT,phase_s);
        rho_s = out.rho;  % Density at the surface
        K_s = K;
        alpha_s = out.alpha;  % Thermal expansivity at the surface
        Cp_s = out.Cp;  % Heat capacity at the surface
        dT_dz_s = alpha_s*g(1)*T_s/Cp_s;  % Thermal gradient at the surface
        T(1) = T_s;
        P(1) = P_s;
        rho(1) = rho_s;
        alpha(1) = alpha_s;
        Cp(1) = Cp_s;
        K(1) = K_s;
        dT_dz(1) = dT_dz_s;
        phase(1) = phase_s;
        for i = 1:(size(z,2)-1)%Integration with depth
            T(i+1) = T(i)+dT_dz(i)*(z(i+1)-z(i));
            P(i+1) = P(i)+ rho(i) * g(i) * (z(i+1)-z(i))*1e-6;
            PT = [P(i+1) T(i+1)];
            phase(i+1) = SF_WhichPhase(PT);
            [out,oo] = compute_params(PT,phase(i+1));
            K(i+1) =oo;
            rho(i+1) = out.rho;
            alpha(i+1) = out.alpha;
            Cp(i+1) = out.Cp;
            dT_dz(i+1) = alpha(i+1)*g(i+1)*T(i+1)/Cp(i+1);
        end
     end
% end