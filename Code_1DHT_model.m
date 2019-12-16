%{
This MATLAB-script describes the 1DHT model code developed for the research first presented in Hornum et al. (2020). 
The 1DHT model is a transient one-dimensional heat transfer model suitable for simulating permafrost dynamics. 
The core of the model is an explicit forward-difference time approximation of the one-dimensional heat transfer equation. 
See the supplementary info to Hornum et al. (2020) for model code validation.
This script is tailored to simulate the Holocene ground temperature
development in Adventdalen, Svalbard, but may be modified to fit other purposes.

To save the figures (simulations) produced by this script. Run 'Save_data.m' consecutively. 

Usage must be cited by reference to Hornum et al. (2020). Any
questions or other inquiries should be directed to toftmikkel@gmail.com

For references cited below see:
Hornum, M. T., Hodson, A., Jessen, S., Bense, V. and Senger, K.: Numerical modelling of permafrost spring discharge and open-system pingo formation induced by basal permafrost aggradation, Cryosphere, 2020.
%}

clc, clear, 
close all force
%%      Constants and general model setup
toyr =60*60*24*365;     % to year from s

            %%% Thermal properties of water and ice (from Williams and Smith, 1989) %%%
            
p_w = 1000;             % density [kg/m3]
p_ice = 917;
cp_w = 4180;            % specific heat capacity [J/(kg*K)]
cp_ice = 2100;
k_w = 0.56*toyr;        % thermal conductivity [(J/yr)/(m K)]
k_ice = 2.24*toyr;
L=333.6*1000;           % Latent heat of fusion [J/kg]

            %%% Porosity and thermal properties of earth materials %%%
            
n_Sc=1;                 % Porosity scenario is chosen here: 1=minimum n, 2=intermediate n, 3=maximum n
nQall=[0.3 0.4 0.5];    % [min interm. max]
nSndall=[0.06 0.1 0.15];
nShall=[0.1 0.2 0.3];

                    %%% Material 1 - Quaternary sediments
nQ = nQall(n_Sc);        
p_soilQ = 2400;                      
cp_soilQ = 850;    
k_soilQ = 0.5*toyr;                 
aQ=k_soilQ/(p_soilQ*cp_soilQ);
                    %%% Material 2 - Sandstone
nSnd =nSndall(n_Sc);       
p_soilSnd =2600;                 
cp_soilSnd =900;
k_soilSnd =2.5*toyr;
aSnd=k_soilSnd/(p_soilSnd*cp_soilSnd);

                    %%% Material 3 - Shale
nSh =nShall(n_Sc);
p_soilSh =2600;
cp_soilSh =800;
k_soilSh =1.5*toyr;
aSh=k_soilSh/(p_soilSh*cp_soilSh);

            %%% Model parameters %%%
tstep=0.05;         % time step [yr]
ts_1yr=1/tstep;     % number of time steps in a year
dz=2;               % cell size [m]
grid_depth=300;     % grid depth
z=0:dz:grid_depth;  % cell nodes
nocell=numel(z);    % number of cells

            %%% Temperature %%%
Kv=0;               % set to -273.15 if in Kelvin. 0 if in C
T_gradient=1/40;    % Thermal gradient
T_0=0+Kv;           % Initial surface T. NOT USED if temperature reconstruction is defined
T_end=grid_depth*T_gradient+Kv; % Initial temperature at bottom of grid - NOT USED
T_ini=z*T_gradient+Kv; % Initial temperature distribution
T=zeros(1,nocell);

                    % Holocene temperature reconstruction (11 to 0 ka) based on Mangerud and
                    % Svendsen (2017). See  Sect. 3.2 in Hornum et al. (2020)
TCurve_data=load('HoloceneTemperatureCurve.txt'); % 'HoloceneTemperatureCurve.txt' should be placed in same folder as this script
T_sum=TCurve_data(:,2); % MSAT [degree C or K] at these
tt_yr=TCurve_data(:,1); % times [yr]
T_ann10=T_sum-10;       % MAAT

ttt=0:tstep:11000;      % Converts temperature curve to simulation time step %%
xT_ann10q=interp1(tt_yr,T_ann10,ttt); % Interpolates temperatures for each sim. time step

T_ann10q=xT_ann10q;     % These temperatures are used in the model for-loop

%%      Temperature curve used for simulation
figure(15)    %%%  MAAT during simulation period
plot(ttt,T_ann10q)
xlabel('Time [yr]')
ylabel('MAAT [C^{o}]')
axis([0 (max(ttt)+500) -8 1])
title('Temperature curve used for simulation')
movegui('northeast')
%%      Assigning material IDs to the different model grids and creating empty vectors and matrices.  

            %%% Defining the 1D grids (from hereon 'columns') %%%
% 1=Q, 2=Sndstone & 3=Shale
onegrid=ones(1,nocell);
                % Material IDs ZONEs 2-3 (3), 1-2 (2), 0-1 (1)
Z0_1MatID=ones(1,nocell);
Z0_1MatID(31:80)=2;
Z0_1MatID(81:nocell)=3;
Z1_2MatID=Z0_1MatID;
Z2_3MatID=Z0_1MatID;
                % Cell properties ZONE 5-6 (6), 4-5 (5), 3-4 (4)
Z3_4MatID=ones(1,nocell);
Z3_4MatID(28:65)=2;
Z3_4MatID(66:nocell)=3;
Z4_5MatID=Z3_4MatID;
Z5_6MatID=Z3_4MatID;
                % Cell properties ZONE 8-9 (9), 7-8 (8) & 6-7 (7)
Z6_7MatID=ones(1,nocell);
Z6_7MatID(26:55)=2;
Z6_7MatID(56:nocell)=3;
Z7_8MatID=Z6_7MatID;
Z8_9MatID=Z6_7MatID;
                % Cell properties ZONE 9-10A (10)
Z9_10AMatID=ones(1,nocell);
Z9_10AMatID(6:20)=2;
Z9_10AMatID(21:nocell)=3;
                % Cell properties ZONE 9-10B (11)
Z9_10BMatID=ones(1,nocell);
Z9_10BMatID(6:nocell)=3;
                % Cell properties ZONE 10 (12)
Z10MatID=ones(1,nocell);
Z10MatID(6:nocell)=3;

            %%% Creating one-vectors and empty matrices used in for-loop %%%
n_grid=onegrid;
p_soil=onegrid;
cp_soil=onegrid;
k_soil=onegrid;

F_w(1)=0;
F_w(nocell)=1;
f_w=onegrid;
f_w(1)=0;
F_ice(1)=1;
F_ice(nocell)=0;
f_ice(nocell)=0;

dF_w=zeros(1,(nocell-1));
k_eq=zeros(1,(nocell-1));
C_eq=zeros(1,(nocell));
a_eq=zeros(1,(nocell));

T_matx_0_1=zeros(nocell,500);
T_matx_1_2=zeros(nocell,1500);
T_matx_2_3=zeros(nocell,2500);
T_matx_3_4=zeros(nocell,3500);
T_matx_4_5=zeros(nocell,4500);
T_matx_5_6=zeros(nocell,5500);
T_matx_6_7=zeros(nocell,6500);
T_matx_7_8=zeros(nocell,7500);
T_matx_8_9=zeros(nocell,8500);
T_matx_9_10a=zeros(nocell,9500);
T_matx_9_10b=zeros(nocell,9500);
T_matx_10=zeros(nocell,10000);

F_matx_0_1=zeros(nocell,500);
F_matx_1_2=zeros(nocell,1500);
F_matx_2_3=zeros(nocell,2500);
F_matx_3_4=zeros(nocell,3500);
F_matx_4_5=zeros(nocell,4500);
F_matx_5_6=zeros(nocell,5500);
F_matx_6_7=zeros(nocell,6500);
F_matx_7_8=zeros(nocell,7500);
F_matx_8_9=zeros(nocell,8500);
F_matx_9_10a=zeros(nocell,9500);
F_matx_9_10b=zeros(nocell,9500);
F_matx_10=zeros(nocell,10000);

w=0.957725;     % Correction factor
%%      %%% Numerical model %%%
col_incl=([1:12]);       %%% Specifiy which of the 12 columns (zones) to include in the ground temperature simulation (e.g. [1:12] (all), [4] (only one), [1 5 11] (several specific ones).
for col=col_incl
  if col==1             %%% Defining the column in use
      runtime=500;      % Simulation runtime
      materialid=Z0_1MatID;
  elseif col==2
      runtime=1500;
      materialid=Z1_2MatID;
  elseif col==3
      runtime=2500;
      materialid=Z2_3MatID;
  elseif col==4
      runtime=3500;
      materialid=Z3_4MatID;
  elseif col==5
      runtime=4500;
      materialid=Z4_5MatID;
  elseif col==6
      runtime=5500;
      materialid=Z5_6MatID;
  elseif col==7
      runtime=6500;
      materialid=Z6_7MatID;
  elseif col==8
      runtime=7500;
      materialid=Z7_8MatID;
  elseif col==9
      runtime=8500;
      materialid=Z8_9MatID;
  elseif col==10
      runtime=9500;
      materialid=Z9_10AMatID;
  elseif col==11
      runtime=9500;
      materialid=Z9_10BMatID;
  elseif col==12
      runtime=10000;
      materialid=Z10MatID;
  end
  for ii=1:nocell       %%% Assigning thermal parameters to column %%%
    if materialid(ii)==1                %%% Porosity
        n_grid(ii)=onegrid(ii)*nQ;
    elseif materialid(ii)==2
        n_grid(ii)=onegrid(ii)*nSnd;
    elseif materialid(ii)==3
        n_grid(ii)=onegrid(ii)*nSh; 
    end
    if materialid(ii)==1                %%% Density
        p_soil(ii)=onegrid(ii)*p_soilQ;
    elseif materialid(ii)==2
        p_soil(ii)=onegrid(ii)*p_soilSnd;
    elseif materialid(ii)==3
        p_soil(ii)=onegrid(ii)*p_soilSh; 
    end
    if materialid(ii)==1                %%% Specific heat capacity
        cp_soil(ii)=onegrid(ii)*cp_soilQ;
    elseif materialid(ii)==2
        cp_soil(ii)=onegrid(ii)*cp_soilSnd;
    elseif materialid(ii)==3
        cp_soil(ii)=onegrid(ii)*cp_soilSh; 
    end
    if materialid(ii)==1                %%% Thermal conductivity
        k_soil(ii)=onegrid(ii)*k_soilQ;
    elseif materialid(ii)==2
        k_soil(ii)=onegrid(ii)*k_soilSnd;
    elseif materialid(ii)==3
        k_soil(ii)=onegrid(ii)*k_soilSh; 
    end
    f_w(nocell)=n_grid(nocell)*1;
    f_ice(1)=n_grid(1)*1;
    f_soil=1-n_grid;
  end
  
no_tstep=runtime/tstep; % Total number of time step
T_11=T_ann10q(:,1:(no_tstep+1));   %%% Cutting temperature curve to sim period %%
T_1=fliplr(T_11);       %%% New T_1 %%%

              %%% Stability if <0.5
k_s_ice=(n_grid*sqrt(k_ice)+(1-n_grid).*sqrt(k_soil)).^2;
C_s_ice=n_grid*p_ice*cp_ice+(1-n_grid).*p_soil.*cp_soil;
a_s_ice=k_s_ice./C_s_ice;
stability=(max(a_s_ice)*tstep/dz^2);
                
                 %%%  Heat transfer loops begins here  %%%
    if col==1       %%% 1st column begins %%%
        xT_ini=T_ini;
        k=0;     %%% time step count
        l=0;     %%% year county
        for t=1:no_tstep
            if stability>0.5
                continue
            end
            k=k+1;
                for i=2:(nocell-1)
                              %%% Fraction of w and ice within pore space %%%
                        if xT_ini(i)<-2
                            F_w(i)=0;
                        elseif xT_ini(i)>0
                            F_w(i)=1;
                        else
                            F_w(i)=exp(-(xT_ini(i)/w)^2);   % %% eq. A, Freezing curve
                        end
                            F_ice(i)=1-F_w(i);

                              %%% dF_w/dT %%%
                        if xT_ini(i)<-2
                            dF_w(i)=0;
                        elseif xT_ini(i)>0
                            dF_w(i)=0;
                        else
                            dF_w(i)=-2*xT_ini(i)*exp(-(xT_ini(i)/w)^2); %%% diff of eq. A
                        end
                              %%% Total fractions of w and ice
                        f_w(i)=F_w(i)*n_grid(i);
                        f_ice(i) = F_ice(i)*n_grid(i);

                              %%% Equivalent thermal parameters %%%
                    k_eq(i)=(f_soil(i)*sqrt(k_soil(i))+f_w(i)*sqrt(k_w)+f_ice(i)*sqrt(k_ice))^2;
                    C_eq(i)=f_soil(i)*p_soil(i)*cp_soil(i)+f_w(i)*p_w*(cp_w)+f_ice(i)*p_ice*(cp_ice+L*dF_w(i));
                    a_eq(i)=k_eq(i)/C_eq(i);

                              %%% Forward-difference time approximation of the 1D heat equation %%%
                    if (xT_ini(i+1)-xT_ini(i-1))~=0
                        T(t,i)=(a_eq(i)*tstep/dz^2)*(xT_ini(i+1)+xT_ini(i-1))+((1-(2*a_eq(i)*tstep/dz^2))*xT_ini(i));
                    else
                        T(t,i) = xT_ini(i);  
                    end
                              %%% End-grid temperatures %%%
                    T(t,1)=T_1(t);
                    T(t,nocell)=T(t,(nocell-1))+dz*T_gradient;
                end      
                xT_ini=T(t,:);      %%% Value used in loop
                
                a_eq_0_1=a_eq;      %%% Vectors for figures         <<<<<<<<<<<<<<<--------------------
                C_eq_0_1=C_eq;
                F_w_0_1=F_w;
                F_ice_0_1=F_ice;
                              %%% Creating a matrix with the temperature distribution for
                              %%% each year
            if k==ts_1yr
                l=l+1;

                T_matx_0_1(:,l)=T(t,:);      %%% Vectors for figures         <<<<<<<<<<<<<<<--------------
                F_matx_0_1(:,l)=F_w;

                k=0;
            end
        end
                    %%% 1st column ends %%%    
    elseif col==2   %%% 2nd column begins %%%
        xT_ini=T_ini;
        k=0;     %%% time step count
        l=0;     %%% year count
        for t=1:no_tstep
            if stability>0.5
                continue
            end

            k=k+1;
                for i=2:(nocell-1)
                              %%% Fraction of w and ice within pore space %%%
                        if xT_ini(i)<-2
                            F_w(i)=0;
                        elseif xT_ini(i)>0
                            F_w(i)=1;
                        else
                            F_w(i)=exp(-(xT_ini(i)/w)^2);   %%% eq. A, Freezing curve
                        end
                            F_ice(i)=1-F_w(i);

                              %%% dF_w/dT %%%
                        if xT_ini(i)<-2
                            dF_w(i)=0;
                        elseif xT_ini(i)>0
                            dF_w(i)=0;
                        else
                            dF_w(i)=-2*xT_ini(i)*exp(-(xT_ini(i)/w)^2); %%% diff of eq. A
                        end
                              %%% Total fractions of w and ice
                        f_w(i)=F_w(i)*n_grid(i);
                        f_ice(i) = F_ice(i)*n_grid(i);

                              %%% Equivalent thermal parameters %%%
                    k_eq(i)=(f_soil(i)*sqrt(k_soil(i))+f_w(i)*sqrt(k_w)+f_ice(i)*sqrt(k_ice))^2;
                    C_eq(i)=f_soil(i)*p_soil(i)*cp_soil(i)+f_w(i)*p_w*(cp_w)+f_ice(i)*p_ice*(cp_ice+L*dF_w(i));
                    a_eq(i)=k_eq(i)/C_eq(i);

                              %%% FDM 1D heat equation %%%
                    if (xT_ini(i+1)-xT_ini(i-1))~=0
                        T(t,i)=(a_eq(i)*tstep/dz^2)*(xT_ini(i+1)+xT_ini(i-1))+((1-(2*a_eq(i)*tstep/dz^2))*xT_ini(i));
                    else
                        T(t,i) = xT_ini(i);  
                    end
                              %%% End-grid temperatures %%%
                    T(t,1)=T_1(t);
                    T(t,nocell)=T(t,(nocell-1))+dz*T_gradient;
                end      
                xT_ini=T(t,:);      %%% Value used in loop
                
                a_eq_1_2=a_eq;      %%% Vectors for figures  <<<<<<<<------ adjust names
                C_eq_1_2=C_eq;
                F_w_1_2=F_w;
                F_ice_1_2=F_ice;
                              %%% Creating a matrix with the temperature distribution for
                              %%% each year
            if k==ts_1yr
                l=l+1;

                T_matx_1_2(:,l)=T(t,:);    %%% Vectors for figures  <<<<<<<<------ adjust names
                F_matx_1_2(:,l)=F_w;

                k=0;
            end
        end
                    %%% 2nd column ends %%% 
    elseif col==3   %%% 3rd column begins %%%
        xT_ini=T_ini;
        k=0;     %%% time step count
        l=0;     %%% year count
        for t=1:no_tstep
            if stability>0.5
                continue
            end

            k=k+1;
                for i=2:(nocell-1)
                              %%% Fraction of w and ice within pore space %%%
                        if xT_ini(i)<-2
                            F_w(i)=0;
                        elseif xT_ini(i)>0
                            F_w(i)=1;
                        else
                            F_w(i)=exp(-(xT_ini(i)/w)^2);   %%% eq. A, Freezing curve
                        end
                            F_ice(i)=1-F_w(i);

                              %%% dF_w/dT %%%
                        if xT_ini(i)<-2
                            dF_w(i)=0;
                        elseif xT_ini(i)>0
                            dF_w(i)=0;
                        else
                            dF_w(i)=-2*xT_ini(i)*exp(-(xT_ini(i)/w)^2); %%% diff of eq. A
                        end
                              %%% Total fractions of w and ice
                        f_w(i)=F_w(i)*n_grid(i);
                        f_ice(i) = F_ice(i)*n_grid(i);

                              %%% Equivalent thermal parameters %%%
                    k_eq(i)=(f_soil(i)*sqrt(k_soil(i))+f_w(i)*sqrt(k_w)+f_ice(i)*sqrt(k_ice))^2;
                    C_eq(i)=f_soil(i)*p_soil(i)*cp_soil(i)+f_w(i)*p_w*(cp_w)+f_ice(i)*p_ice*(cp_ice+L*dF_w(i));
                    a_eq(i)=k_eq(i)/C_eq(i);

                              %%% FDM 1D heat equation %%%
                    if (xT_ini(i+1)-xT_ini(i-1))~=0
                        T(t,i)=(a_eq(i)*tstep/dz^2)*(xT_ini(i+1)+xT_ini(i-1))+((1-(2*a_eq(i)*tstep/dz^2))*xT_ini(i));
                    else
                        T(t,i) = xT_ini(i);  
                    end
                              %%% End-grid temperatures %%%
                    T(t,1)=T_1(t);
                    T(t,nocell)=T(t,(nocell-1))+dz*T_gradient;
                end      
                xT_ini=T(t,:);      %%% Value used in loop
                
                a_eq_2_3=a_eq;      %%% Vectors for figures  <<<<<<<<------ adjust names
                C_eq_2_3=C_eq;
                F_w_2_3=F_w;
                F_ice_2_3=F_ice;
                              %%% Creating a matrix with the temperature distribution for
                              %%% each year
            if k==ts_1yr
                l=l+1;

                T_matx_2_3(:,l)=T(t,:);    %%% Vectors for figures  <<<<<<<<------ adjust names
                F_matx_2_3(:,l)=F_w;

                k=0;
            end
        end
                    %%% 3rd column ends %%% 
    elseif col==4   %%% 4th column begins %%%
        xT_ini=T_ini;
        k=0;     %%% time step count
        l=0;     %%% year count
        for t=1:no_tstep
            if stability>0.5
                continue
            end

            k=k+1;
                for i=2:(nocell-1)
                              %%% Fraction of w and ice within pore space %%%
                        if xT_ini(i)<-2
                            F_w(i)=0;
                        elseif xT_ini(i)>0
                            F_w(i)=1;
                        else
                            F_w(i)=exp(-(xT_ini(i)/w)^2);   %%% eq. A, Freezing curve
                        end
                            F_ice(i)=1-F_w(i);

                              %%% dF_w/dT %%%
                        if xT_ini(i)<-2
                            dF_w(i)=0;
                        elseif xT_ini(i)>0
                            dF_w(i)=0;
                        else
                            dF_w(i)=-2*xT_ini(i)*exp(-(xT_ini(i)/w)^2); %%% diff of eq. A
                        end
                              %%% Total fractions of w and ice
                        f_w(i)=F_w(i)*n_grid(i);
                        f_ice(i) = F_ice(i)*n_grid(i);

                              %%% Equivalent thermal parameters %%%
                    k_eq(i)=(f_soil(i)*sqrt(k_soil(i))+f_w(i)*sqrt(k_w)+f_ice(i)*sqrt(k_ice))^2;
                    C_eq(i)=f_soil(i)*p_soil(i)*cp_soil(i)+f_w(i)*p_w*(cp_w)+f_ice(i)*p_ice*(cp_ice+L*dF_w(i));
                    a_eq(i)=k_eq(i)/C_eq(i);

                              %%% FDM 1D heat equation %%%
                    if (xT_ini(i+1)-xT_ini(i-1))~=0
                        T(t,i)=(a_eq(i)*tstep/dz^2)*(xT_ini(i+1)+xT_ini(i-1))+((1-(2*a_eq(i)*tstep/dz^2))*xT_ini(i));
                    else
                        T(t,i) = xT_ini(i);  
                    end
                              %%% End-grid temperatures %%%
                    T(t,1)=T_1(t);
                    T(t,nocell)=T(t,(nocell-1))+dz*T_gradient;
                end      
                xT_ini=T(t,:);      %%% Value used in loop
                
                a_eq_3_4=a_eq;      %%% Vectors for figures  <<<<<<<<------ adjust names
                C_eq_3_4=C_eq;
                F_w_3_4=F_w;
                F_ice_3_4=F_ice;
                              %%% Creating a matrix with the temperature distribution for
                              %%% each year
            if k==ts_1yr
                l=l+1;

                T_matx_3_4(:,l)=T(t,:);    %%% Vectors for figures  <<<<<<<<------ adjust names
                F_matx_3_4(:,l)=F_w;

                k=0;
            end
        end
                    %%% 4th column ends %%% 
    elseif col==5   %%% 5th column begins %%%
        xT_ini=T_ini;
        k=0;     %%% time step count
        l=0;     %%% year count
        for t=1:no_tstep
            if stability>0.5
                continue
            end
            k=k+1;
                for i=2:(nocell-1)
                              %%% Fraction of w and ice within pore space %%%
                        if xT_ini(i)<-2
                            F_w(i)=0;
                        elseif xT_ini(i)>0
                            F_w(i)=1;
                        else
                            F_w(i)=exp(-(xT_ini(i)/w)^2);   %%% eq. A, Freezing curve
                        end
                            F_ice(i)=1-F_w(i);

                              %%% dF_w/dT %%%
                        if xT_ini(i)<-2
                            dF_w(i)=0;
                        elseif xT_ini(i)>0
                            dF_w(i)=0;
                        else
                            dF_w(i)=-2*xT_ini(i)*exp(-(xT_ini(i)/w)^2); %%% diff of eq. A
                        end
                              %%% Total fractions of w and ice
                        f_w(i)=F_w(i)*n_grid(i);
                        f_ice(i) = F_ice(i)*n_grid(i);

                              %%% Equivalent thermal parameters %%%
                    k_eq(i)=(f_soil(i)*sqrt(k_soil(i))+f_w(i)*sqrt(k_w)+f_ice(i)*sqrt(k_ice))^2;
                    C_eq(i)=f_soil(i)*p_soil(i)*cp_soil(i)+f_w(i)*p_w*(cp_w)+f_ice(i)*p_ice*(cp_ice+L*dF_w(i));
                    a_eq(i)=k_eq(i)/C_eq(i);

                              %%% FDM 1D heat equation %%%
                    if (xT_ini(i+1)-xT_ini(i-1))~=0
                        T(t,i)=(a_eq(i)*tstep/dz^2)*(xT_ini(i+1)+xT_ini(i-1))+((1-(2*a_eq(i)*tstep/dz^2))*xT_ini(i));
                    else
                        T(t,i) = xT_ini(i);  
                    end
                              %%% End-grid temperatures %%%
                    T(t,1)=T_1(t);
                    T(t,nocell)=T(t,(nocell-1))+dz*T_gradient;
                end      
                xT_ini=T(t,:);      %%% Value used in loop
                
                a_eq_4_5=a_eq;      %%% Vectors for figures  <<<<<<<<------ adjust names
                C_eq_4_5=C_eq;
                F_w_4_5=F_w;
                F_ice_4_5=F_ice;
                              %%% Creating a matrix with the temperature distribution for
                              %%% each year
            if k==ts_1yr
                l=l+1;

                T_matx_4_5(:,l)=T(t,:);    %%% Vectors for figures  <<<<<<<<------ adjust names
                F_matx_4_5(:,l)=F_w;

                k=0;
            end
        end
                    %%% 5th column ends %%% 
    elseif col==6   %%% 6th column begins %%%
        xT_ini=T_ini;
        k=0;     %%% time step count
        l=0;     %%% year count
        for t=1:no_tstep
            if stability>0.5
                continue
            end

            k=k+1;
                for i=2:(nocell-1)
                              %%% Fraction of w and ice within pore space %%%
                        if xT_ini(i)<-2
                            F_w(i)=0;
                        elseif xT_ini(i)>0
                            F_w(i)=1;
                        else
                            F_w(i)=exp(-(xT_ini(i)/w)^2);   %%% eq. A, Freezing curve
                        end
                            F_ice(i)=1-F_w(i);

                              %%% dF_w/dT %%%
                        if xT_ini(i)<-2
                            dF_w(i)=0;
                        elseif xT_ini(i)>0
                            dF_w(i)=0;
                        else
                            dF_w(i)=-2*xT_ini(i)*exp(-(xT_ini(i)/w)^2); %%% diff of eq. A
                        end
                              %%% Total fractions of w and ice
                        f_w(i)=F_w(i)*n_grid(i);
                        f_ice(i) = F_ice(i)*n_grid(i);

                              %%% Equivalent thermal parameters %%%
                    k_eq(i)=(f_soil(i)*sqrt(k_soil(i))+f_w(i)*sqrt(k_w)+f_ice(i)*sqrt(k_ice))^2;
                    C_eq(i)=f_soil(i)*p_soil(i)*cp_soil(i)+f_w(i)*p_w*(cp_w)+f_ice(i)*p_ice*(cp_ice+L*dF_w(i));
                    a_eq(i)=k_eq(i)/C_eq(i);

                              %%% FDM 1D heat equation %%%
                    if (xT_ini(i+1)-xT_ini(i-1))~=0
                        T(t,i)=(a_eq(i)*tstep/dz^2)*(xT_ini(i+1)+xT_ini(i-1))+((1-(2*a_eq(i)*tstep/dz^2))*xT_ini(i));
                    else
                        T(t,i) = xT_ini(i);  
                    end
                              %%% End-grid temperatures %%%
                    T(t,1)=T_1(t);
                    T(t,nocell)=T(t,(nocell-1))+dz*T_gradient;
                end      
                xT_ini=T(t,:);      %%% Value used in loop
                
                a_eq_5_6=a_eq;      %%% Vectors for figures  <<<<<<<<------ adjust names
                C_eq_5_6=C_eq;
                F_w_5_6=F_w;
                F_ice_5_6=F_ice;
                              %%% Creating a matrix with the temperature distribution for
                              %%% each year
            if k==ts_1yr
                l=l+1;

                T_matx_5_6(:,l)=T(t,:);    %%% Vectors for figures  <<<<<<<<------ adjust names
                F_matx_5_6(:,l)=F_w;

                k=0;
            end
        end
                    %%% 6th column ends %%% 
    elseif col==7   %%% 7th column begins %%%
        xT_ini=T_ini;
        k=0;     %%% time step count
        l=0;     %%% year count
        for t=1:no_tstep
            if stability>0.5
                continue
            end

            k=k+1;
                for i=2:(nocell-1)
                              %%% Fraction of w and ice within pore space %%%
                        if xT_ini(i)<-2
                            F_w(i)=0;
                        elseif xT_ini(i)>0
                            F_w(i)=1;
                        else
                            F_w(i)=exp(-(xT_ini(i)/w)^2);   %%% eq. A, Freezing curve
                        end
                            F_ice(i)=1-F_w(i);

                              %%% dF_w/dT %%%
                        if xT_ini(i)<-2
                            dF_w(i)=0;
                        elseif xT_ini(i)>0
                            dF_w(i)=0;
                        else
                            dF_w(i)=-2*xT_ini(i)*exp(-(xT_ini(i)/w)^2); %%% diff of eq. A
                        end
                              %%% Total fractions of w and ice
                        f_w(i)=F_w(i)*n_grid(i);
                        f_ice(i) = F_ice(i)*n_grid(i);

                              %%% Equivalent thermal parameters %%%
                    k_eq(i)=(f_soil(i)*sqrt(k_soil(i))+f_w(i)*sqrt(k_w)+f_ice(i)*sqrt(k_ice))^2;
                    C_eq(i)=f_soil(i)*p_soil(i)*cp_soil(i)+f_w(i)*p_w*(cp_w)+f_ice(i)*p_ice*(cp_ice+L*dF_w(i));
                    a_eq(i)=k_eq(i)/C_eq(i);

                              %%% FDM 1D heat equation %%%
                    if (xT_ini(i+1)-xT_ini(i-1))~=0
                        T(t,i)=(a_eq(i)*tstep/dz^2)*(xT_ini(i+1)+xT_ini(i-1))+((1-(2*a_eq(i)*tstep/dz^2))*xT_ini(i));
                    else
                        T(t,i) = xT_ini(i);  
                    end
                              %%% End-grid temperatures %%%
                    T(t,1)=T_1(t);
                    T(t,nocell)=T(t,(nocell-1))+dz*T_gradient;
                end      
                xT_ini=T(t,:);      %%% Value used in loop
                
                a_eq_6_7=a_eq;      %%% Vectors for figures  <<<<<<<<------ adjust names
                C_eq_6_7=C_eq;
                F_w_6_7=F_w;
                F_ice_6_7=F_ice;
                              %%% Creating a matrix with the temperature distribution for
                              %%% each year
            if k==ts_1yr
                l=l+1;

                T_matx_6_7(:,l)=T(t,:);    %%% Vectors for figures  <<<<<<<<------ adjust names
                F_matx_6_7(:,l)=F_w;

                k=0;
            end
        end
                    %%% 7th column ends %%% 
    elseif col==8   %%% 8th column begins %%%
        xT_ini=T_ini;
        k=0;     %%% time step count
        l=0;     %%% year count
        for t=1:no_tstep
            if stability>0.5
                continue
            end

            k=k+1;
                for i=2:(nocell-1)
                              %%% Fraction of w and ice within pore space %%%
                        if xT_ini(i)<-2
                            F_w(i)=0;
                        elseif xT_ini(i)>0
                            F_w(i)=1;
                        else
                            F_w(i)=exp(-(xT_ini(i)/w)^2);   %%% eq. A, Freezing curve
                        end
                            F_ice(i)=1-F_w(i);

                              %%% dF_w/dT %%%
                        if xT_ini(i)<-2
                            dF_w(i)=0;
                        elseif xT_ini(i)>0
                            dF_w(i)=0;
                        else
                            dF_w(i)=-2*xT_ini(i)*exp(-(xT_ini(i)/w)^2); %%% diff of eq. A
                        end
                              %%% Total fractions of w and ice
                        f_w(i)=F_w(i)*n_grid(i);
                        f_ice(i) = F_ice(i)*n_grid(i);

                              %%% Equivalent thermal parameters %%%
                    k_eq(i)=(f_soil(i)*sqrt(k_soil(i))+f_w(i)*sqrt(k_w)+f_ice(i)*sqrt(k_ice))^2;
                    C_eq(i)=f_soil(i)*p_soil(i)*cp_soil(i)+f_w(i)*p_w*(cp_w)+f_ice(i)*p_ice*(cp_ice+L*dF_w(i));
                    a_eq(i)=k_eq(i)/C_eq(i);

                              %%% FDM 1D heat equation %%%
                    if (xT_ini(i+1)-xT_ini(i-1))~=0
                        T(t,i)=(a_eq(i)*tstep/dz^2)*(xT_ini(i+1)+xT_ini(i-1))+((1-(2*a_eq(i)*tstep/dz^2))*xT_ini(i));
                    else
                        T(t,i) = xT_ini(i);  
                    end
                              %%% End-grid temperatures %%%
                    T(t,1)=T_1(t);
                    T(t,nocell)=T(t,(nocell-1))+dz*T_gradient;
                end      
                xT_ini=T(t,:);      %%% Value used in loop
                
                a_eq_7_8=a_eq;      %%% Vectors for figures  <<<<<<<<------ adjust names
                C_eq_7_8=C_eq;
                F_w_7_8=F_w;
                F_ice_7_8=F_ice;
                              %%% Creating a matrix with the temperature distribution for
                              %%% each year
            if k==ts_1yr
                l=l+1;

                T_matx_7_8(:,l)=T(t,:);    %%% Vectors for figures  <<<<<<<<------ adjust names
                F_matx_7_8(:,l)=F_w;

                k=0;
            end
        end
                    %%% 8th column ends %%% 
    elseif col==9   %%% 9th column begins %%%
        xT_ini=T_ini;
        k=0;     %%% time step count
        l=0;     %%% year count
        for t=1:no_tstep
            if stability>0.5
                continue
            end

            k=k+1;
                for i=2:(nocell-1)
                              %%% Fraction of w and ice within pore space %%%
                        if xT_ini(i)<-2
                            F_w(i)=0;
                        elseif xT_ini(i)>0
                            F_w(i)=1;
                        else
                            F_w(i)=exp(-(xT_ini(i)/w)^2);   %%% eq. A, Freezing curve
                        end
                            F_ice(i)=1-F_w(i);

                              %%% dF_w/dT %%%
                        if xT_ini(i)<-2
                            dF_w(i)=0;
                        elseif xT_ini(i)>0
                            dF_w(i)=0;
                        else
                            dF_w(i)=-2*xT_ini(i)*exp(-(xT_ini(i)/w)^2); %%% diff of eq. A
                        end
                              %%% Total fractions of w and ice
                        f_w(i)=F_w(i)*n_grid(i);
                        f_ice(i) = F_ice(i)*n_grid(i);

                              %%% Equivalent thermal parameters %%%
                    k_eq(i)=(f_soil(i)*sqrt(k_soil(i))+f_w(i)*sqrt(k_w)+f_ice(i)*sqrt(k_ice))^2;
                    C_eq(i)=f_soil(i)*p_soil(i)*cp_soil(i)+f_w(i)*p_w*(cp_w)+f_ice(i)*p_ice*(cp_ice+L*dF_w(i));
                    a_eq(i)=k_eq(i)/C_eq(i);

                              %%% FDM 1D heat equation %%%
                    if (xT_ini(i+1)-xT_ini(i-1))~=0
                        T(t,i)=(a_eq(i)*tstep/dz^2)*(xT_ini(i+1)+xT_ini(i-1))+((1-(2*a_eq(i)*tstep/dz^2))*xT_ini(i));
                    else
                        T(t,i) = xT_ini(i);  
                    end
                              %%% End-grid temperatures %%%
                    T(t,1)=T_1(t);
                    T(t,nocell)=T(t,(nocell-1))+dz*T_gradient;
                end      
                xT_ini=T(t,:);      %%% Value used in loop
                
                a_eq_8_9=a_eq;      %%% Vectors for figures  <<<<<<<<------ adjust names
                C_eq_8_9=C_eq;
                F_w_8_9=F_w;
                F_ice_8_9=F_ice;
                              %%% Creating a matrix with the temperature distribution for
                              %%% each year
            if k==ts_1yr
                l=l+1;

                T_matx_8_9(:,l)=T(t,:);    %%% Vectors for figures  <<<<<<<<------ adjust names
                F_matx_8_9(:,l)=F_w;

                k=0;
            end
        end
                    %%% 9th column ends %%% 
    elseif col==10   %%% 10th column begins %%%
        xT_ini=T_ini;
        k=0;     %%% time step count
        l=0;     %%% year count
        for t=1:no_tstep
            if stability>0.5
                continue
            end

            k=k+1;
                for i=2:(nocell-1)
                              %%% Fraction of w and ice within pore space %%%
                        if xT_ini(i)<-2
                            F_w(i)=0;
                        elseif xT_ini(i)>0
                            F_w(i)=1;
                        else
                            F_w(i)=exp(-(xT_ini(i)/w)^2);   %%% eq. A, Freezing curve
                        end
                            F_ice(i)=1-F_w(i);

                              %%% dF_w/dT %%%
                        if xT_ini(i)<-2
                            dF_w(i)=0;
                        elseif xT_ini(i)>0
                            dF_w(i)=0;
                        else
                            dF_w(i)=-2*xT_ini(i)*exp(-(xT_ini(i)/w)^2); %%% diff of eq. A
                        end
                              %%% Total fractions of w and ice
                        f_w(i)=F_w(i)*n_grid(i);
                        f_ice(i) = F_ice(i)*n_grid(i);

                              %%% Equivalent thermal parameters %%%
                    k_eq(i)=(f_soil(i)*sqrt(k_soil(i))+f_w(i)*sqrt(k_w)+f_ice(i)*sqrt(k_ice))^2;
                    C_eq(i)=f_soil(i)*p_soil(i)*cp_soil(i)+f_w(i)*p_w*(cp_w)+f_ice(i)*p_ice*(cp_ice+L*dF_w(i));
                    a_eq(i)=k_eq(i)/C_eq(i);

                              %%% FDM 1D heat equation %%%
                    if (xT_ini(i+1)-xT_ini(i-1))~=0
                        T(t,i)=(a_eq(i)*tstep/dz^2)*(xT_ini(i+1)+xT_ini(i-1))+((1-(2*a_eq(i)*tstep/dz^2))*xT_ini(i));
                    else
                        T(t,i) = xT_ini(i);  
                    end
                              %%% End-grid temperatures %%%
                    T(t,1)=T_1(t);
                    T(t,nocell)=T(t,(nocell-1))+dz*T_gradient;
                end      
                xT_ini=T(t,:);      %%% Value used in loop
                
                a_eq_9_10a=a_eq;      %%% Vectors for figures  <<<<<<<<------ adjust names
                C_eq_9_10a=C_eq;
                F_w_9_10a=F_w;
                F_ice_9_10a=F_ice;
                              %%% Creating a matrix with the temperature distribution for
                              %%% each year
            if k==ts_1yr
                l=l+1;

                T_matx_9_10a(:,l)=T(t,:);    %%% Vectors for figures  <<<<<<<<------ adjust names
                F_matx_9_10a(:,l)=F_w;

                k=0;
            end
        end
                    %%% 10th column ends %%% 
    elseif col==11   %%% 11th column begins %%%
        xT_ini=T_ini;
        k=0;     %%% time step count
        l=0;     %%% year count
        for t=1:no_tstep
            if stability>0.5
                continue
            end

            k=k+1;
                for i=2:(nocell-1)
                              %%% Fraction of w and ice within pore space %%%
                        if xT_ini(i)<-2
                            F_w(i)=0;
                        elseif xT_ini(i)>0
                            F_w(i)=1;
                        else
                            F_w(i)=exp(-(xT_ini(i)/w)^2);   %%% eq. A, Freezing curve
                        end
                            F_ice(i)=1-F_w(i);

                              %%% dF_w/dT %%%
                        if xT_ini(i)<-2
                            dF_w(i)=0;
                        elseif xT_ini(i)>0
                            dF_w(i)=0;
                        else
                            dF_w(i)=-2*xT_ini(i)*exp(-(xT_ini(i)/w)^2); %%% diff of eq. A
                        end
                              %%% Total fractions of w and ice
                        f_w(i)=F_w(i)*n_grid(i);
                        f_ice(i) = F_ice(i)*n_grid(i);

                              %%% Equivalent thermal parameters %%%
                    k_eq(i)=(f_soil(i)*sqrt(k_soil(i))+f_w(i)*sqrt(k_w)+f_ice(i)*sqrt(k_ice))^2;
                    C_eq(i)=f_soil(i)*p_soil(i)*cp_soil(i)+f_w(i)*p_w*(cp_w)+f_ice(i)*p_ice*(cp_ice+L*dF_w(i));
                    a_eq(i)=k_eq(i)/C_eq(i);

                              %%% FDM 1D heat equation %%%
                    if (xT_ini(i+1)-xT_ini(i-1))~=0
                        T(t,i)=(a_eq(i)*tstep/dz^2)*(xT_ini(i+1)+xT_ini(i-1))+((1-(2*a_eq(i)*tstep/dz^2))*xT_ini(i));
                    else
                        T(t,i) = xT_ini(i);  
                    end
                              %%% End-grid temperatures %%%
                    T(t,1)=T_1(t);
                    T(t,nocell)=T(t,(nocell-1))+dz*T_gradient;
                end      
                xT_ini=T(t,:);      %%% Value used in loop
                
                a_eq_9_10b=a_eq;      %%% Vectors for figures  <<<<<<<<------ adjust names
                C_eq_9_10b=C_eq;
                F_w_9_10b=F_w;
                F_ice_9_10b=F_ice;
                              %%% Creating a matrix with the temperature distribution for
                              %%% each year
            if k==ts_1yr
                l=l+1;

                T_matx_9_10b(:,l)=T(t,:);    %%% Vectors for figures  <<<<<<<<------ adjust names
                F_matx_9_10b(:,l)=F_w;

                k=0;
            end
        end
                    %%% 11th column ends %%% 
    elseif col==12   %%% 12th column begins %%%
        xT_ini=T_ini;
        k=0;     %%% time step count
        l=0;     %%% year count
        for t=1:no_tstep
            if stability>0.5
                continue
            end

            k=k+1;
                for i=2:(nocell-1)
                              %%% Fraction of w and ice within pore space %%%
                        if xT_ini(i)<-2
                            F_w(i)=0;
                        elseif xT_ini(i)>0
                            F_w(i)=1;
                        else
                            F_w(i)=exp(-(xT_ini(i)/w)^2);   %%% eq. A, Freezing curve
                        end
                            F_ice(i)=1-F_w(i);

                              %%% dF_w/dT %%%
                        if xT_ini(i)<-2
                            dF_w(i)=0;
                        elseif xT_ini(i)>0
                            dF_w(i)=0;
                        else
                            dF_w(i)=-2*xT_ini(i)*exp(-(xT_ini(i)/w)^2); %%% diff of eq. A
                        end
                              %%% Total fractions of w and ice
                        f_w(i)=F_w(i)*n_grid(i);
                        f_ice(i) = F_ice(i)*n_grid(i);

                              %%% Equivalent thermal parameters %%%
                    k_eq(i)=(f_soil(i)*sqrt(k_soil(i))+f_w(i)*sqrt(k_w)+f_ice(i)*sqrt(k_ice))^2;
                    C_eq(i)=f_soil(i)*p_soil(i)*cp_soil(i)+f_w(i)*p_w*(cp_w)+f_ice(i)*p_ice*(cp_ice+L*dF_w(i));
                    a_eq(i)=k_eq(i)/C_eq(i);

                              %%% FDM 1D heat equation %%%
                    if (xT_ini(i+1)-xT_ini(i-1))~=0
                        T(t,i)=(a_eq(i)*tstep/dz^2)*(xT_ini(i+1)+xT_ini(i-1))+((1-(2*a_eq(i)*tstep/dz^2))*xT_ini(i));
                    else
                        T(t,i) = xT_ini(i);  
                    end
                              %%% End-grid temperatures %%%
                    T(t,1)=T_1(t);
                    T(t,nocell)=T(t,(nocell-1))+dz*T_gradient;
                end      
                xT_ini=T(t,:);      %%% Value used in loop
                
                a_eq_10=a_eq;      %%% Vectors for figures  <<<<<<<<------ adjust names
                C_eq_10=C_eq;
                F_w_10=F_w;
                F_ice_10=F_ice;
                              %%% Creating a matrix with the temperature distribution for
                              %%% each year
            if k==ts_1yr
                l=l+1;

                T_matx_10(:,l)=T(t,:);    %%% Vectors for figures  <<<<<<<<------ adjust names
                F_matx_10(:,l)=F_w;

                k=0;
            end
        end
                    %%% 12th column ends %%% 
    end
end
toc
%%                          %%% Figures %%%
%%              %%% Fig - Final temperature distrubution %%%  checked
%%
if any(col_incl==1)
figure(11)  %%%  1st column - Zone 0-1
set(gca,'YDir', 'reverse')
grid on
movegui('northwest')
xlabel('Temperature ({\circ}C)')
ylabel('Depth (m)')
axis([-6 4 0 grid_depth])
title('Temperature profiles - Zone 0-1')
hold on
plot(T_matx_0_1(:,500),z)
hold on 
plot(T_ini,z,'k--')
hleg=legend('0','500','1500','2500' ,'3500','4500','5500', '6500', '7500', '8500', '9500');
title(hleg,'Time [yr BP]')
end
if any(col_incl==2)
figure(21)  %%%  2nd column - Zone 1-2
set(gca,'YDir', 'reverse')
grid on
movegui('northwest')
xlabel('Temperature ({\circ}C)')
ylabel('Depth (m)')
axis([-6 4 0 grid_depth])
title('Temperature profiles - Zone 1-2')
hold on
plot(T_matx_1_2(:,1500),z)
hold on
plot(T_matx_1_2(:,500),z)
hold on 
plot(T_ini,z,'k--')
hleg=legend('0','500','1500','2500' ,'3500','4500','5500', '6500', '7500', '8500', '9500');
title(hleg,'Time [yr BP]')
end
if any(col_incl==3)
figure(31)  %%%  3rd column - Zone 2-3
set(gca,'YDir', 'reverse')
grid on
movegui('northwest')
xlabel('Temperature ({\circ}C)')
ylabel('Depth (m)')
axis([-6 4 0 grid_depth])
title('Temperature profiles - Zone 2-3')
hold on
plot(T_matx_2_3(:,2500),z)
hold on
plot(T_matx_2_3(:,1500),z)
hold on
plot(T_matx_2_3(:,500),z)
hold on 
plot(T_ini,z,'k--')
hleg=legend('0','500','1500','2500' ,'3500','4500','5500', '6500', '7500', '8500', '9500');
title(hleg,'Time [yr BP]')
end
if any(col_incl==4)
figure(41)  %%%  2nd column - Zone 3-4
set(gca,'YDir', 'reverse')
grid on
movegui('northwest')
xlabel('Temperature ({\circ}C)')
ylabel('Depth (m)')
axis([-6 4 0 grid_depth])
title('Temperature profiles - Zone 3-4')
hold on
plot(T_matx_3_4(:,3500),z)
hold on
plot(T_matx_3_4(:,2500),z)
hold on
plot(T_matx_3_4(:,1500),z)
hold on
plot(T_matx_3_4(:,500),z)
hold on 
plot(T_ini,z,'k--')
hleg=legend('0','500','1500','2500' ,'3500','4500','5500', '6500', '7500', '8500', '9500');
title(hleg,'Time [yr BP]')
end
if any(col_incl==5)
figure(51)  %%%  2nd column - Zone 4-5
set(gca,'YDir', 'reverse')
grid on
movegui('northwest')
xlabel('Temperature ({\circ}C)')
ylabel('Depth (m)')
axis([-6 4 0 grid_depth])
title('Temperature profiles - Zone 4-5')
hold on
plot(T_matx_4_5(:,4500),z)
hold on
plot(T_matx_4_5(:,3500),z)
hold on
plot(T_matx_4_5(:,2500),z)
hold on
plot(T_matx_4_5(:,1500),z)
hold on
plot(T_matx_4_5(:,500),z)
hold on 
plot(T_ini,z,'k--')
hleg=legend('0','500','1500','2500' ,'3500','4500','5500', '6500', '7500', '8500', '9500');
title(hleg,'Time [yr BP]')
end
if any(col_incl==6)
figure(61)  %%%  Zone 5-6
set(gca,'YDir', 'reverse')
grid on
movegui('northwest')
xlabel('Temperature ({\circ}C)')
ylabel('Depth (m)')
axis([-6 4 0 grid_depth])
title('Temperature profiles - Zone 5-6')
hold on
plot(T_matx_5_6(:,5500),z)
hold on
plot(T_matx_5_6(:,4500),z)
hold on
plot(T_matx_5_6(:,3500),z)
hold on
plot(T_matx_5_6(:,2500),z)
hold on
plot(T_matx_5_6(:,1500),z)
hold on
plot(T_matx_5_6(:,500),z)
hold on 
plot(T_ini,z,'k--')
hleg=legend('0','500','1500','2500' ,'3500','4500','5500', '6500', '7500', '8500', '9500');
title(hleg,'Time [yr BP]')
end
if any(col_incl==7)
figure(71)  %%%  Zone 6_7
set(gca,'YDir', 'reverse')
grid on
movegui('northwest')
xlabel('Temperature ({\circ}C)')
ylabel('Depth (m)')
axis([-6 4 0 grid_depth])
title('Temperature profiles - Zone 6-7')
hold on
plot(T_matx_6_7(:,6500),z)
hold on
plot(T_matx_6_7(:,5500),z)
hold on
plot(T_matx_6_7(:,4500),z)
hold on
plot(T_matx_6_7(:,3500),z)
hold on
plot(T_matx_6_7(:,2500),z)
hold on
plot(T_matx_6_7(:,1500),z)
hold on
plot(T_matx_6_7(:,500),z)
hold on 
plot(T_ini,z,'k--')
hleg=legend('0','500','1500','2500' ,'3500','4500','5500', '6500', '7500', '8500', '9500');
title(hleg,'Time [yr BP]')
end
if any(col_incl==8)
figure(81)  %%%  2nd column - Zone 7_8
set(gca,'YDir', 'reverse')
grid on
movegui('northwest')
xlabel('Temperature ({\circ}C)')
ylabel('Depth (m)')
axis([-6 4 0 grid_depth])
title('Temperature profiles - Zone 7-8')
hold on
plot(T_matx_7_8(:,7500),z)
hold on
plot(T_matx_7_8(:,6500),z)
hold on
plot(T_matx_7_8(:,5500),z)
hold on
plot(T_matx_7_8(:,4500),z)
hold on
plot(T_matx_7_8(:,3500),z)
hold on
plot(T_matx_7_8(:,2500),z)
hold on
plot(T_matx_7_8(:,1500),z)
hold on
plot(T_matx_7_8(:,500),z)
hold on 
plot(T_ini,z,'k--')
hleg=legend('0','500','1500','2500' ,'3500','4500','5500', '6500', '7500', '8500', '9500');
title(hleg,'Time [yr BP]')
end
if any(col_incl==9)
figure(91)  %%% Zone 8_9
set(gca,'YDir', 'reverse')
grid on
movegui('northwest')
xlabel('Temperature ({\circ}C)')
ylabel('Depth (m)')
axis([-6 4 0 grid_depth])
title('Temperature profiles - Zone 8-9')
hold on
plot(T_matx_8_9(:,8500),z)
hold on
plot(T_matx_8_9(:,7500),z)
hold on
plot(T_matx_8_9(:,6500),z)
hold on
plot(T_matx_8_9(:,5500),z)
hold on
plot(T_matx_8_9(:,4500),z)
hold on
plot(T_matx_8_9(:,3500),z)
hold on
plot(T_matx_8_9(:,2500),z)
hold on
plot(T_matx_8_9(:,1500),z)
hold on
plot(T_matx_8_9(:,500),z)
hold on 
plot(T_ini,z,'k--')
hleg=legend('0','500','1500','2500' ,'3500','4500','5500', '6500', '7500', '8500', '9500');
title(hleg,'Time [yr BP]')
end
if any(col_incl==10)
figure(101)  %%%  Zone 9_10a
set(gca,'YDir', 'reverse')
grid on
movegui('northwest')
xlabel('Temperature ({\circ}C)')
ylabel('Depth (m)')
axis([-6 4 0 grid_depth])
title('Temperature profiles - Zone 9-10a')
hold on
plot(T_matx_9_10a(:,9500),z)
hold on
plot(T_matx_9_10a(:,8500),z)
hold on
plot(T_matx_9_10a(:,7500),z)
hold on
plot(T_matx_9_10a(:,6500),z)
hold on
plot(T_matx_9_10a(:,5500),z)
hold on
plot(T_matx_9_10a(:,4500),z)
hold on
plot(T_matx_9_10a(:,3500),z)
hold on
plot(T_matx_9_10a(:,2500),z)
hold on
plot(T_matx_9_10a(:,1500),z)
hold on
plot(T_matx_9_10a(:,500),z)
hold on 
plot(T_ini,z,'k--')
hleg=legend('0','500','1500','2500' ,'3500','4500','5500', '6500', '7500', '8500', '9500');
title(hleg,'Time [yr BP]')
end
if any(col_incl==11)
figure(111)  %%%  Zone 9_10b
set(gca,'YDir', 'reverse')
grid on
movegui('northwest')
xlabel('Temperature ({\circ}C)')
ylabel('Depth (m)')
axis([-6 4 0 grid_depth])
title('Temperature profiles - Zone 9-10b')
hold on
plot(T_matx_9_10b(:,9500),z)
hold on
plot(T_matx_9_10b(:,8500),z)
hold on
plot(T_matx_9_10b(:,7500),z)
hold on
plot(T_matx_9_10b(:,6500),z)
hold on
plot(T_matx_9_10b(:,5500),z)
hold on
plot(T_matx_9_10b(:,4500),z)
hold on
plot(T_matx_9_10b(:,3500),z)
hold on
plot(T_matx_9_10b(:,2500),z)
hold on
plot(T_matx_9_10b(:,1500),z)
hold on
plot(T_matx_9_10b(:,500),z)
hold on 
plot(T_ini,z,'k--')
hleg=legend('0','500','1500','2500' ,'3500','4500','5500', '6500', '7500', '8500', '9500');
title(hleg,'Time [yr BP]')
end
if any(col_incl==12)
figure(121)  %%%  Zone 10
set(gca,'YDir', 'reverse')
grid on
movegui('northwest')
xlabel('Temperature ({\circ}C)')
ylabel('Depth (m)')
axis([-6 4 0 grid_depth])
title('Temperature profiles - Zone 10')
hold on
plot(T_matx_10(:,10000),z)
hold on
plot(T_matx_10(:,8500),z)
hold on
plot(T_matx_10(:,7500),z)
hold on
plot(T_matx_10(:,6500),z)
hold on
plot(T_matx_10(:,5500),z)
hold on
plot(T_matx_10(:,4500),z)
hold on
plot(T_matx_10(:,3500),z)
hold on
plot(T_matx_10(:,2500),z)
hold on
plot(T_matx_10(:,1500),z)
hold on
plot(T_matx_10(:,500),z)
hold on 
plot(T_ini,z,'k--')
hleg=legend('0','500','1500','2500' ,'3500','4500','5500', '6500', '7500', '8500', '10000');
title(hleg,'Time [yr BP]')
end
%%              %%%  Final equiv. vol. heat capacity distribution  %%%  checked
%%
if any(col_incl==1)
figure(12)   %%% 1st column - Zone 0-1
plot(C_eq_0_1(1,2:100), z(1,2:100),'+');
title('Volumetric heat capacity - Zone 0-1')
set(gca,'YDir', 'reverse')
grid on
movegui('west')
end
if any(col_incl==2)
figure(22)   %%% 2nd column
plot(C_eq_1_2(1,2:100), z(1,2:100),'+');
title('Volumetric heat capacity - Zone 1-2')
set(gca,'YDir', 'reverse')
grid on
movegui('west')
end
if any(col_incl==3)
figure(32)   %%%  Zone 2-3
plot(C_eq_2_3(1,2:100), z(1,2:100),'+');
title('Volumetric heat capacity - Zone 2-3')
set(gca,'YDir', 'reverse')
grid on
movegui('west')
end
if any(col_incl==4)
figure(42)   %%% 4th column - Zone 3-4
plot(C_eq_3_4(1,2:100), z(1,2:100),'+');
title('Volumetric heat capacity - Zone 3-4')
set(gca,'YDir', 'reverse')
grid on
movegui('west')
end
if any(col_incl==5)
figure(52)   %%% Zone 4-5
plot(C_eq_4_5(1,2:100), z(1,2:100),'+');
title('Volumetric heat capacity - Zone 4-5')
set(gca,'YDir', 'reverse')
grid on
movegui('west')
end
if any(col_incl==6)
figure(62)   %%% Zone 5-6
plot(C_eq_5_6(1,2:100), z(1,2:100),'+');
title('Volumetric heat capacity - Zone 5-6')
set(gca,'YDir', 'reverse')
grid on
movegui('west')
end
if any(col_incl==7)
figure(72)   %%% Zone 6-7
plot(C_eq_6_7(1,2:100), z(1,2:100),'+');
title('Volumetric heat capacity - Zone 6-7')
set(gca,'YDir', 'reverse')
grid on
movegui('west')
end
if any(col_incl==8)
figure(82)   %%% Zone 7-8
plot(C_eq_7_8(1,2:100), z(1,2:100),'+');
title('Volumetric heat capacity - Zone 7-8')
set(gca,'YDir', 'reverse')
grid on
movegui('west')
end
if any(col_incl==9)
figure(92)   %%% Zone 8-9
plot(C_eq_8_9(1,2:100), z(1,2:100),'+');
title('Volumetric heat capacity - Zone 8-9')
set(gca,'YDir', 'reverse')
grid on
movegui('west')
end
if any(col_incl==10)
figure(102)   %%% Zone 9-10a
plot(C_eq_9_10a(1,2:100), z(1,2:100),'+');
title('Volumetric heat capacity - Zone 9-10a')
set(gca,'YDir', 'reverse')
grid on
movegui('west')
end
if any(col_incl==11)
figure(112)   %%% Zone 9-10b
plot(C_eq_9_10b(1,2:100), z(1,2:100),'+');
title('Volumetric heat capacity - Zone 9-10b')
set(gca,'YDir', 'reverse')
grid on
movegui('west')
end
if any(col_incl==12)
figure(122)   %%% Zone 10
plot(C_eq_10(1,2:100), z(1,2:100),'+');
title('Volumetric heat capacity - Zone 10')
set(gca,'YDir', 'reverse')
grid on
movegui('west')
end
%%              %%%  Equivalent thermal diffusivity %%%   Checked
%%
if any(col_incl==1)
figure(13)   %%% 1st column - Zone 0-1
hold on
plot(a_eq_0_1(1,2:(nocell-1)), z(1,2:(nocell-1)),':');
title('Equivalent thermal diffusivity - Zone 0-1')
set(gca,'YDir', 'reverse')
axis([0 30 0 200])
ylabel('Depth [m]')
xlabel('\alpha_{eq} [m^{2}/yr]')
legend('Scenario 3', 'Scenario 2')
% grid on
% grid minor
movegui('southwest')
end
if any(col_incl==2)
figure(23)   %%% Zone 1-2
hold on
plot(a_eq_1_2(1,2:(nocell-1)), z(1,2:(nocell-1)),':');
title('Equivalent thermal diffusivity - Zone 1-2')
set(gca,'YDir', 'reverse')
axis([0 30 0 200])
ylabel('Depth [m]')
xlabel('\alpha_{eq} [m^{2}/yr]')
legend('Scenario 3', 'Scenario 2')
% grid on
% grid minor
movegui('southwest')
end
if any(col_incl==3)
figure(33)   %%%  3rd column - Zone 2-3
hold on
plot(a_eq_2_3(1,2:(nocell-1)), z(1,2:(nocell-1)),':');
title('Equivalent thermal diffusivity - Zone 2-3')
set(gca,'YDir', 'reverse')
axis([0 30 0 200])
ylabel('Depth [m]')
xlabel('\alpha_{eq} [m^{2}/yr]')
legend('Scenario 3', 'Scenario 2')
% grid on
% grid minor
movegui('southwest')
end
if any(col_incl==4)
figure(43)   %%% Zone 3-4
hold on
plot(a_eq_3_4(1,2:(nocell-1)), z(1,2:(nocell-1)),':');
title('Equivalent thermal diffusivity - Zone 3-4')
set(gca,'YDir', 'reverse')
axis([0 30 0 200])
ylabel('Depth [m]')
xlabel('\alpha_{eq} [m^{2}/yr]')
legend('Scenario 3', 'Scenario 2')
% grid on
% grid minor
movegui('southwest')
end
if any(col_incl==5)
figure(53)   %%% Zone 4-5
hold on
plot(a_eq_4_5(1,2:(nocell-1)), z(1,2:(nocell-1)),':');
title('Equivalent thermal diffusivity - Zone 4-5')
set(gca,'YDir', 'reverse')
axis([0 30 0 200])
ylabel('Depth [m]')
xlabel('\alpha_{eq} [m^{2}/yr]')
legend('Scenario 3', 'Scenario 2')
% grid on
% grid minor
movegui('southwest')
end
if any(col_incl==6)
figure(63)   %%% Zone 5-6
hold on
plot(a_eq_5_6(1,2:(nocell-1)), z(1,2:(nocell-1)),':');
title('Equivalent thermal diffusivity - Zone 5-6')
set(gca,'YDir', 'reverse')
axis([0 30 0 200])
ylabel('Depth [m]')
xlabel('\alpha_{eq} [m^{2}/yr]')
legend('Scenario 3', 'Scenario 2')
% grid on
% grid minor
movegui('southwest')
end
if any(col_incl==7)
figure(73)   %%% Zone 6-7
hold on
plot(a_eq_6_7(1,2:(nocell-1)), z(1,2:(nocell-1)),':');
title('Equivalent thermal diffusivity - Zone 6-7')
set(gca,'YDir', 'reverse')
axis([0 30 0 200])
ylabel('Depth [m]')
xlabel('\alpha_{eq} [m^{2}/yr]')
legend('Scenario 3', 'Scenario 2')
% grid on
% grid minor
movegui('southwest')
end
if any(col_incl==8)
figure(83)   %%% Zone 7-8
hold on
plot(a_eq_7_8(1,2:(nocell-1)), z(1,2:(nocell-1)),':');
title('Equivalent thermal diffusivity - Zone 7-8')
set(gca,'YDir', 'reverse')
axis([0 30 0 200])
ylabel('Depth [m]')
xlabel('\alpha_{eq} [m^{2}/yr]')
legend('Scenario 3', 'Scenario 2')
% grid on
% grid minor
movegui('southwest')
end
if any(col_incl==9)
figure(93)   %%% Zone 8-9
hold on
plot(a_eq_8_9(1,2:(nocell-1)), z(1,2:(nocell-1)),':');
title('Equivalent thermal diffusivity - Zone 8-9')
set(gca,'YDir', 'reverse')
axis([0 30 0 200])
ylabel('Depth [m]')
xlabel('\alpha_{eq} [m^{2}/yr]')
legend('Scenario 3', 'Scenario 2')
% grid on
% grid minor
movegui('southwest')
end
if any(col_incl==10)
figure(103)   %%% Zone 9-10a
hold on
plot(a_eq_9_10a(1,2:(nocell-1)), z(1,2:(nocell-1)),':');
title('Equivalent thermal diffusivity - Zone 9-10a')
set(gca,'YDir', 'reverse')
axis([0 30 0 200])
ylabel('Depth [m]')
xlabel('\alpha_{eq} [m^{2}/yr]')
legend('Scenario 3', 'Scenario 2')
% grid on
% grid minor
movegui('southwest')
end
if any(col_incl==11)
figure(113)   %%% Zone 9-10b
hold on
plot(a_eq_9_10b(1,2:(nocell-1)), z(1,2:(nocell-1)),':');
title('Equivalent thermal diffusivity - Zone 9-10b')
set(gca,'YDir', 'reverse')
axis([0 30 0 200])
ylabel('Depth [m]')
xlabel('\alpha_{eq} [m^{2}/yr]')
legend('Scenario 3', 'Scenario 2')
% grid on
% grid minor
movegui('southwest')
end
if any(col_incl==12)
figure(123)   %%% Zone 10
hold on
plot(a_eq_10(1,2:(nocell-1)), z(1,2:(nocell-1)),':');
title('Equivalent thermal diffusivity - Zone 10')
set(gca,'YDir', 'reverse')
axis([0 30 0 200])
ylabel('Depth [m]')
xlabel('\alpha_{eq} [m^{2}/yr]')
legend('Scenario 3', 'Scenario 2')
% grid on
% grid minor
movegui('southwest')
end
%%              %%%  Final fraction distribution  Checked
%%
if any(col_incl==1)
figure(14)   %% Zone 0-1
plot(F_w_0_1,z,'r--')
hold on 
plot(F_ice_0_1,z,'b-')
set(gca,'YDir', 'reverse')
movegui('north')
grid on
title('Fractions - Zone 0-1')
end
if any(col_incl==2)
figure(24)      %% Zone 1-2
plot(F_w_1_2,z,'r--')
hold on 
plot(F_ice_1_2,z,'b-')
set(gca,'YDir', 'reverse')
movegui('north')
grid on
title('Fractions - Zone 1-2')
end
if any(col_incl==3)
figure(34)   %%%  Zone 2-3
plot(F_w_2_3,z,'r--')
hold on 
plot(F_ice_2_3,z,'b-')
set(gca,'YDir', 'reverse')
movegui('north')
grid on
title('Fractions - Zone 2-3')
end
if any(col_incl==4)
figure(44)   %% Zone 3-4
plot(F_w_3_4,z,'r--')
hold on 
plot(F_ice_3_4,z,'b-')
set(gca,'YDir', 'reverse')
movegui('north')
grid on
title('Fractions - Zone 3-4')
end
if any(col_incl==5)
figure(54)   %% Zone 4-5
plot(F_w_4_5,z,'r--')
hold on 
plot(F_ice_4_5,z,'b-')
set(gca,'YDir', 'reverse')
movegui('north')
grid on
title('Fractions - Zone 4-5')
end
if any(col_incl==6)
figure(64)   %% Zone 5-6
plot(F_w_5_6,z,'r--')
hold on 
plot(F_ice_5_6,z,'b-')
set(gca,'YDir', 'reverse')
movegui('north')
grid on
title('Fractions - Zone 5-6')
end
if any(col_incl==7)
figure(74)   %% Zone 6_7
plot(F_w_6_7,z,'r--')
hold on 
plot(F_ice_6_7,z,'b-')
set(gca,'YDir', 'reverse')
movegui('north')
grid on
title('Fractions - 6-7')
end
if any(col_incl==8)
figure(84)   %% Zone 7-8
plot(F_w_7_8,z,'r--')
hold on 
plot(F_ice_7_8,z,'b-')
set(gca,'YDir', 'reverse')
movegui('north')
grid on
title('Fractions - Zone 7-8')
end
if any(col_incl==9)
figure(94)   %%   Zone 8-9
plot(F_w_8_9,z,'r--')
hold on 
plot(F_ice_8_9,z,'b-')
set(gca,'YDir', 'reverse')
movegui('north')
grid on
title('Fractions - Zone 8-9')
end
if any(col_incl==10)
figure(104)   %%   Zone 9-10a
plot(F_w_9_10a,z,'r--')
hold on 
plot(F_ice_9_10a,z,'b-')
set(gca,'YDir', 'reverse')
movegui('north')
grid on
title('Fractions - Zone 9-10a')
end
if any(col_incl==11)
figure(114)   %%   Zone 9-10b
plot(F_w_9_10b,z,'r--')
hold on 
plot(F_ice_9_10b,z,'b-')
set(gca,'YDir', 'reverse')
movegui('north')
grid on
title('Fractions - Zone 9-10b')
end
if any(col_incl==12)
figure(124)   %%   Zone 10
plot(F_w_10,z,'r--')
hold on 
plot(F_ice_10,z,'b-')
set(gca,'YDir', 'reverse')
movegui('north')
grid on
title('Fractions - Zone 10')
end
%%              %%% Aggradation rate and depth of PF & FF
%%
TqFF=(-2);
Tq_1=-1;
TqPF=0;
t_intv=1;         %%% Time interval between colomns in T_matx [yr]
nocell_1=nocell-1;

Dcol0=0;
Dcol1=615;
Dcol2=1965;
Dcol3=3375;
Dcol4=4505;
Dcol5=5220;
Dcol6=5670;
Dcol7=6020;
Dcol8=6570;
Dcol9=7337.5;
Dcol10=8712.5;
Dcol11=12305;
Dcol12=15000;
DcolFar=18000;    %%% value not used

presentFFdepth(1)=0;
presentPFdepth(1)=0;
presentAggrRate(1)=0;
DtoDF(1)=0;

colcount=1;

if any(col_incl==1)
    l=500;                  %%% runtime [yr]
    k=0;
    zqFF=zeros(1,l);
    zq_1=zeros(1,l);
    zqPF=zeros(1,l);
    temp=zeros(nocell,l);
    tt=1:t_intv:l;
    tt_shift=1:t_intv:(l-1);
    aggr_rateFF=zeros(1,(l-1));
    aggr_rate_1=aggr_rateFF;
    aggr_ratePF=aggr_rateFF;

    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_0_1(1:nocell,k);             %%% use right T_matx

        zqFF(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),TqFF);
        tt(i)=t_intv*k;
    end

    k=0;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_0_1(1:nocell,k);             %%% use right T_matx

        zq_1(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),Tq_1);
        tt(i)=t_intv*k;
    end

    k=0;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_0_1(1:nocell,k);             %%% use right T_matx

        zqPF(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),TqPF);
    %     ttPF(i)=t_intv*k;
    end
    
    colcount=colcount+1;
    presentFFdepth(colcount)=zqFF(length(zqFF));
    presentPFdepth(colcount)=zqPF(length(zqPF));
    DtoDF(colcount)=Dcol1;

    no_points=numel(zqFF);
    tt_shifted=tt-1/2*t_intv;
    tt_shift(1:(no_points-1))=tt_shifted(2:no_points);
    tt_shift_BP=flip(tt_shift);
    tt_flip=flip(tt);

    for ii=1:(no_points-1)
        aggr_rateFF(ii)=(zqFF(ii+1)-zqFF(ii))/t_intv;
    end

    for ii=1:(no_points-1)
        aggr_rate_1(ii)=(zq_1(ii+1)-zq_1(ii))/t_intv;
    end

    for ii=1:(no_points-1)
        aggr_ratePF(ii)=(zqPF(ii+1)-zqPF(ii))/t_intv;
    end

    presentAggrRate(colcount)=aggr_rate_1(length(aggr_rate_1));
    
    figure(17)       %%% Freezing front and PF depth
    plot(tt_flip,zqFF,'-.')
    hold on
    plot(tt_flip,zqPF,'-.')
    title('Freezing front and permafrost depth - Zone 0-1')
    set(gca,'YDir', 'reverse')
    xlabel('Time [yr BP]')
    ylabel('Depth [m]')
    movegui('center')

    figure(16)         %%% Aggradation rate
    % hold on
    % plot(tt_shift_BP,aggr_rateFF,'b--')
    hold on
    plot(tt_shift_BP,aggr_rate_1,'k')
    % hold on
    % plot(tt_shift_BP,aggr_ratePF,'r--')
    xlabel('Time [yr BP]')
    ylabel('Rate [m/yr]')
    title('Aggradation rate - Zone 0-1')
    movegui('south')
    grid on
    grid minor
end
if any(col_incl==2)
    l=1500;     %%% set time [yr]
    k=0;
    zqFF=zeros(1,l);
    zq_1=zeros(1,l);
    zqPF=zeros(1,l);
    temp=zeros(nocell,l);
    tt=1:t_intv:l;
    tt_shift=1:t_intv:(l-1);
    aggr_rateFF=zeros(1,(l-1));
    aggr_rate_1=aggr_rateFF;
    aggr_ratePF=aggr_rateFF;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_1_2(1:nocell,k);       %%% use right T_matx...     

        zqFF(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),TqFF);
        tt(i)=t_intv*k;
    end

    k=0;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_1_2(1:nocell,k);

        zq_1(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),Tq_1);
        tt(i)=t_intv*k;
    end

    k=0;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_1_2(1:nocell,k);

        zqPF(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),TqPF);
    %     ttPF(i)=t_intv*k;
    end
    
    colcount=colcount+1;
    presentFFdepth(colcount)=zqFF(length(zqFF));
    presentPFdepth(colcount)=zqPF(length(zqPF));
    DtoDF(colcount)=Dcol2;
    
    no_points=numel(zqFF);
    tt_shifted=tt-1/2*t_intv;
    tt_shift(1:(no_points-1))=tt_shifted(2:no_points);
    tt_shift_BP=flip(tt_shift);
    tt_flip=flip(tt);

    for ii=1:(no_points-1)
        aggr_rateFF(ii)=(zqFF(ii+1)-zqFF(ii))/t_intv;
    end

    for ii=1:(no_points-1)
        aggr_rate_1(ii)=(zq_1(ii+1)-zq_1(ii))/t_intv;
    end

    for ii=1:(no_points-1)
        aggr_ratePF(ii)=(zqPF(ii+1)-zqPF(ii))/t_intv;
    end

    presentAggrRate(colcount)=aggr_rate_1(length(aggr_rate_1));
    
    figure(27)       %%% Freezing front and PF depth
    plot(tt_flip,zqFF,'-.')
    hold on
    plot(tt_flip,zqPF,'-.')
    title('Freezing front and permafrost depth - Zone 1-2')
    set(gca,'YDir', 'reverse')
    xlabel('Time [yr BP]')
    ylabel('Depth [m]')
    movegui('center')

    figure(26)         %%% Aggradation rate
    % hold on
    % plot(tt_shift_BP,aggr_rateFF,'b--')
    hold on
    plot(tt_shift_BP,aggr_rate_1,'k')
    % hold on
    % plot(tt_shift_BP,aggr_ratePF,'r--')
    xlabel('Time [yr BP]')
    ylabel('Rate [m/yr]')
    title('Aggradation rate - Zone 1-2')
    movegui('south')
    grid on
    grid minor
end
if any(col_incl==3)
    l=2500;     %%% set time [yr]
    k=0;
    zqFF=zeros(1,l);
    zq_1=zeros(1,l);
    zqPF=zeros(1,l);
    temp=zeros(nocell,l);
    tt=1:t_intv:l;
    tt_shift=1:t_intv:(l-1);
    aggr_rateFF=zeros(1,(l-1));
    aggr_rate_1=aggr_rateFF;
    aggr_ratePF=aggr_rateFF;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_2_3(1:nocell,k);       %%% use right T_matx...     

        zqFF(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),TqFF);
        tt(i)=t_intv*k;
    end

    k=0;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_2_3(1:nocell,k);

        zq_1(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),Tq_1);
        tt(i)=t_intv*k;
    end

    k=0;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_2_3(1:nocell,k);

        zqPF(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),TqPF);
    %     ttPF(i)=t_intv*k;
    end
    
    colcount=colcount+1;
    presentPFdepth(colcount)=zqPF(length(zqPF));
    presentFFdepth(colcount)=zqFF(length(zqFF));
    DtoDF(colcount)=Dcol3;
    
    no_points=numel(zqFF);
    tt_shifted=tt-1/2*t_intv;
    tt_shift(1:(no_points-1))=tt_shifted(2:no_points);
    tt_shift_BP=flip(tt_shift);
    tt_flip=flip(tt);

    for ii=1:(no_points-1)
        aggr_rateFF(ii)=(zqFF(ii+1)-zqFF(ii))/t_intv;
    end

    for ii=1:(no_points-1)
        aggr_rate_1(ii)=(zq_1(ii+1)-zq_1(ii))/t_intv;
    end

    for ii=1:(no_points-1)
        aggr_ratePF(ii)=(zqPF(ii+1)-zqPF(ii))/t_intv;
    end

    presentAggrRate(colcount)=aggr_rate_1(length(aggr_rate_1));
    
    figure(37)       %%% Freezing front and PF depth
    plot(tt_flip,zqFF,'-.')
    hold on
    plot(tt_flip,zqPF,'-.')
    title('Freezing front and permafrost depth - Zone 2-3')
    set(gca,'YDir', 'reverse')
    xlabel('Time [yr BP]')
    ylabel('Depth [m]')
    movegui('center')

    figure(36)         %%% Aggradation rate
    % hold on
    % plot(tt_shift_BP,aggr_rateFF,'b--')
    hold on
    plot(tt_shift_BP,aggr_rate_1,'k')
    % hold on
    % plot(tt_shift_BP,aggr_ratePF,'r--')
    xlabel('Time [yr BP]')
    ylabel('Rate [m/yr]')
    title('Aggradation rate - Zone 2-3')
    movegui('south')
    grid on
    grid minor
end
if any(col_incl==4)
    l=3500;     %%% set time [yr]
    k=0;
    zqFF=zeros(1,l);
    zq_1=zeros(1,l);
    zqPF=zeros(1,l);
    temp=zeros(nocell,l);
    tt=1:t_intv:l;
    tt_shift=1:t_intv:(l-1);
    aggr_rateFF=zeros(1,(l-1));
    aggr_rate_1=aggr_rateFF;
    aggr_ratePF=aggr_rateFF;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_3_4(1:nocell,k);       %%% use right T_matx...     

        zqFF(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),TqFF);
        tt(i)=t_intv*k;
    end

    k=0;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_3_4(1:nocell,k);

        zq_1(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),Tq_1);
        tt(i)=t_intv*k;
    end

    k=0;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_3_4(1:nocell,k);

        zqPF(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),TqPF);
    %     ttPF(i)=t_intv*k;
    end
    
    colcount=colcount+1;
    presentPFdepth(colcount)=zqPF(length(zqPF));
    presentFFdepth(colcount)=zqFF(length(zqFF));
    DtoDF(colcount)=Dcol4;
   
    no_points=numel(zqFF);
    tt_shifted=tt-1/2*t_intv;
    tt_shift(1:(no_points-1))=tt_shifted(2:no_points);
    tt_shift_BP=flip(tt_shift);
    tt_flip=flip(tt);

    for ii=1:(no_points-1)
        aggr_rateFF(ii)=(zqFF(ii+1)-zqFF(ii))/t_intv;
    end

    for ii=1:(no_points-1)
        aggr_rate_1(ii)=(zq_1(ii+1)-zq_1(ii))/t_intv;
    end

    for ii=1:(no_points-1)
        aggr_ratePF(ii)=(zqPF(ii+1)-zqPF(ii))/t_intv;
    end

    presentAggrRate(colcount)=aggr_rate_1(length(aggr_rate_1));
    
    figure(47)       %%% Freezing front and PF depth 3_4
    plot(tt_flip,zqFF,'-.')
    hold on
    plot(tt_flip,zqPF,'-.')
    title('Freezing front and permafrost depth - Zone Zone 3-4')
    set(gca,'YDir', 'reverse')
    xlabel('Time [yr BP]')
    ylabel('Depth [m]')
    movegui('center')

    figure(46)         %%% Aggradation rate 3_4
    % hold on
    % plot(tt_shift_BP,aggr_rateFF,'b--')
    hold on
    plot(tt_shift_BP,aggr_rate_1,'k')
    % hold on
    % plot(tt_shift_BP,aggr_ratePF,'r--')
    xlabel('Time [yr BP]')
    ylabel('Rate [m/yr]')
    title('Aggradation rate - Zone Zone 3-4')
    movegui('south')
    grid on
    grid minor
end
if any(col_incl==5)
    l=4500;     %%% set time [yr]
    k=0;
    zqFF=zeros(1,l);
    zq_1=zeros(1,l);
    zqPF=zeros(1,l);
    temp=zeros(nocell,l);
    tt=1:t_intv:l;
    tt_shift=1:t_intv:(l-1);
    aggr_rateFF=zeros(1,(l-1));
    aggr_rate_1=aggr_rateFF;
    aggr_ratePF=aggr_rateFF;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_4_5(1:nocell,k);       %%% use right T_matx...     

        zqFF(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),TqFF);
        tt(i)=t_intv*k;
    end

    k=0;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_4_5(1:nocell,k);       %%% use right T_matx...     

        zq_1(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),Tq_1);
        tt(i)=t_intv*k;
    end

    k=0;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_4_5(1:nocell,k);       %%% use right T_matx...     

        zqPF(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),TqPF);
    %     ttPF(i)=t_intv*k;
    end
    
    colcount=colcount+1;
    presentPFdepth(colcount)=zqPF(length(zqPF));
    presentFFdepth(colcount)=zqFF(length(zqFF));
    DtoDF(colcount)=Dcol5;
    
    no_points=numel(zqFF);
    tt_shifted=tt-1/2*t_intv;
    tt_shift(1:(no_points-1))=tt_shifted(2:no_points);
    tt_shift_BP=flip(tt_shift);
    tt_flip=flip(tt);

    for ii=1:(no_points-1)
        aggr_rateFF(ii)=(zqFF(ii+1)-zqFF(ii))/t_intv;
    end

    for ii=1:(no_points-1)
        aggr_rate_1(ii)=(zq_1(ii+1)-zq_1(ii))/t_intv;
    end

    for ii=1:(no_points-1)
        aggr_ratePF(ii)=(zqPF(ii+1)-zqPF(ii))/t_intv;
    end

    presentAggrRate(colcount)=aggr_rate_1(length(aggr_rate_1));
    
    figure(57)       %%% Freezing front and PF depth  4_5
    plot(tt_flip,zqFF,'-.')
    hold on
    plot(tt_flip,zqPF,'-.')
    title('Freezing front and permafrost depth - Zone 4-5')
    set(gca,'YDir', 'reverse')
    xlabel('Time [yr BP]')
    ylabel('Depth [m]')
    movegui('center')

    figure(56)         %%% Aggradation rate  4_5
    % hold on
    % plot(tt_shift_BP,aggr_rateFF,'b--')
    hold on
    plot(tt_shift_BP,aggr_rate_1,'k')
    % hold on
    % plot(tt_shift_BP,aggr_ratePF,'r--')
    xlabel('Time [yr BP]')
    ylabel('Rate [m/yr]')
    title('Aggradation rate - Zone 4-5')
    movegui('south')
    grid on
    grid minor
end
if any(col_incl==6)
    l=5500;     %%% set time [yr]
    k=0;
    zqFF=zeros(1,l);
    zq_1=zeros(1,l);
    zqPF=zeros(1,l);
    temp=zeros(nocell,l);
    tt=1:t_intv:l;
    tt_shift=1:t_intv:(l-1);
    aggr_rateFF=zeros(1,(l-1));
    aggr_rate_1=aggr_rateFF;
    aggr_ratePF=aggr_rateFF;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_5_6(1:nocell,k);       %%% use right T_matx...     

        zqFF(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),TqFF);
        tt(i)=t_intv*k;
    end

    k=0;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_5_6(1:nocell,k);

        zq_1(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),Tq_1);
        tt(i)=t_intv*k;
    end

    k=0;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_5_6(1:nocell,k);

        zqPF(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),TqPF);
    %     ttPF(i)=t_intv*k;
    end
    
    colcount=colcount+1;
    presentPFdepth(colcount)=zqPF(length(zqPF));
    presentFFdepth(colcount)=zqFF(length(zqFF));
    DtoDF(colcount)=Dcol6;
   
    no_points=numel(zqFF);
    tt_shifted=tt-1/2*t_intv;
    tt_shift(1:(no_points-1))=tt_shifted(2:no_points);
    tt_shift_BP=flip(tt_shift);
    tt_flip=flip(tt);

    for ii=1:(no_points-1)
        aggr_rateFF(ii)=(zqFF(ii+1)-zqFF(ii))/t_intv;
    end

    for ii=1:(no_points-1)
        aggr_rate_1(ii)=(zq_1(ii+1)-zq_1(ii))/t_intv;
    end

    for ii=1:(no_points-1)
        aggr_ratePF(ii)=(zqPF(ii+1)-zqPF(ii))/t_intv;
    end

    presentAggrRate(colcount)=aggr_rate_1(length(aggr_rate_1));
    
    figure(67)       %%% Freezing front and PF depth
    plot(tt_flip,zqFF,'-.')
    hold on
    plot(tt_flip,zqPF,'-.')
    title('Freezing front and permafrost depth - Zone 5-6')
    set(gca,'YDir', 'reverse')
    xlabel('Time [yr BP]')
    ylabel('Depth [m]')
    movegui('center')

    figure(66)         %%% Aggradation rate
    % hold on
    % plot(tt_shift_BP,aggr_rateFF,'b--')
    hold on
    plot(tt_shift_BP,aggr_rate_1,'k')
    % hold on
    % plot(tt_shift_BP,aggr_ratePF,'r--')
    xlabel('Time [yr BP]')
    ylabel('Rate [m/yr]')
    title('Aggradation rate - Zone 5-6')
    movegui('south')
    grid on
    grid minor
end
if any(col_incl==7)
    l=6500;     %%% set time [yr]
    k=0;
    zqFF=zeros(1,l);
    zq_1=zeros(1,l);
    zqPF=zeros(1,l);
    temp=zeros(nocell,l);
    tt=1:t_intv:l;
    tt_shift=1:t_intv:(l-1);
    aggr_rateFF=zeros(1,(l-1));
    aggr_rate_1=aggr_rateFF;
    aggr_ratePF=aggr_rateFF;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_6_7(1:nocell,k);       %%% use right T_matx...     

        zqFF(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),TqFF);
        tt(i)=t_intv*k;
    end

    k=0;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_6_7(1:nocell,k);

        zq_1(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),Tq_1);
        tt(i)=t_intv*k;
    end

    k=0;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_6_7(1:nocell,k);

        zqPF(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),TqPF);
    %     ttPF(i)=t_intv*k;
    end
    
    colcount=colcount+1;
    presentPFdepth(colcount)=zqPF(length(zqPF));
    presentFFdepth(colcount)=zqFF(length(zqFF));
    DtoDF(colcount)=Dcol7;
   
    no_points=numel(zqFF);
    tt_shifted=tt-1/2*t_intv;
    tt_shift(1:(no_points-1))=tt_shifted(2:no_points);
    tt_shift_BP=flip(tt_shift);
    tt_flip=flip(tt);

    for ii=1:(no_points-1)
        aggr_rateFF(ii)=(zqFF(ii+1)-zqFF(ii))/t_intv;
    end

    for ii=1:(no_points-1)
        aggr_rate_1(ii)=(zq_1(ii+1)-zq_1(ii))/t_intv;
    end

    for ii=1:(no_points-1)
        aggr_ratePF(ii)=(zqPF(ii+1)-zqPF(ii))/t_intv;
    end
    
    presentAggrRate(colcount)=aggr_rate_1(length(aggr_rate_1));
    
    figure(77)       %%% Freezing front and PF depth
    plot(tt_flip,zqFF,'-.')
    hold on
    plot(tt_flip,zqPF,'-.')
    title('Freezing front and permafrost depth - Zone 6-7')
    set(gca,'YDir', 'reverse')
    xlabel('Time [yr BP]')
    ylabel('Depth [m]')
    movegui('center')

    figure(76)         %%% Aggradation rate
    % hold on
    % plot(tt_shift_BP,aggr_rateFF,'b--')
    hold on
    plot(tt_shift_BP,aggr_rate_1,'k')
    % hold on
    % plot(tt_shift_BP,aggr_ratePF,'r--')
    xlabel('Time [yr BP]')
    ylabel('Rate [m/yr]')
    title('Aggradation rate - Zone 6-7')
    movegui('south')
    grid on
    grid minor
end
if any(col_incl==8)
    l=7500;     %%% set time [yr]
    k=0;
    zqFF=zeros(1,l);
    zq_1=zeros(1,l);
    zqPF=zeros(1,l);
    temp=zeros(nocell,l);
    tt=1:t_intv:l;
    tt_shift=1:t_intv:(l-1);
    aggr_rateFF=zeros(1,(l-1));
    aggr_rate_1=aggr_rateFF;
    aggr_ratePF=aggr_rateFF;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_7_8(1:nocell,k);       %%% use right T_matx...     

        zqFF(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),TqFF);
        tt(i)=t_intv*k;
    end

    k=0;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_7_8(1:nocell,k);

        zq_1(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),Tq_1);
        tt(i)=t_intv*k;
    end

    k=0;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_7_8(1:nocell,k);

        zqPF(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),TqPF);
    %     ttPF(i)=t_intv*k;
    end
    
    colcount=colcount+1;
    presentPFdepth(colcount)=zqPF(length(zqPF));
    presentFFdepth(colcount)=zqFF(length(zqFF));
    DtoDF(colcount)=Dcol8;
   
    no_points=numel(zqFF);
    tt_shifted=tt-1/2*t_intv;
    tt_shift(1:(no_points-1))=tt_shifted(2:no_points);
    tt_shift_BP=flip(tt_shift);
    tt_flip=flip(tt);

    for ii=1:(no_points-1)
        aggr_rateFF(ii)=(zqFF(ii+1)-zqFF(ii))/t_intv;
    end

    for ii=1:(no_points-1)
        aggr_rate_1(ii)=(zq_1(ii+1)-zq_1(ii))/t_intv;
    end

    for ii=1:(no_points-1)
        aggr_ratePF(ii)=(zqPF(ii+1)-zqPF(ii))/t_intv;
    end
    
    presentAggrRate(colcount)=aggr_rate_1(length(aggr_rate_1));
    
    figure(87)       %%% Freezing front and PF depth
    plot(tt_flip,zqFF,'-.')
    hold on
    plot(tt_flip,zqPF,'-.')
    title('Freezing front and permafrost depth - Zone 7-8')
    set(gca,'YDir', 'reverse')
    xlabel('Time [yr BP]')
    ylabel('Depth [m]')
    movegui('center')

    figure(86)         %%% Aggradation rate
    % hold on
    % plot(tt_shift_BP,aggr_rateFF,'b--')
    hold on
    plot(tt_shift_BP,aggr_rate_1,'k')
    % hold on
    % plot(tt_shift_BP,aggr_ratePF,'r--')
    xlabel('Time [yr BP]')
    ylabel('Rate [m/yr]')
    title('Aggradation rate - Zone 7-8')
    movegui('south')
    grid on
    grid minor
end
if any(col_incl==9)
    l=8500;     %%% set time [yr]
    k=0;
    zqFF=zeros(1,l);
    zq_1=zeros(1,l);
    zqPF=zeros(1,l);
    temp=zeros(nocell,l);
    tt=1:t_intv:l;
    tt_shift=1:t_intv:(l-1);
    aggr_rateFF=zeros(1,(l-1));
    aggr_rate_1=aggr_rateFF;
    aggr_ratePF=aggr_rateFF;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_8_9(1:nocell,k);       %%% use right T_matx...     

        zqFF(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),TqFF);
        tt(i)=t_intv*k;
    end

    k=0;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_8_9(1:nocell,k);

        zq_1(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),Tq_1);
        tt(i)=t_intv*k;
    end

    k=0;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_8_9(1:nocell,k);

        zqPF(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),TqPF);
    %     ttPF(i)=t_intv*k;
    end
    
    colcount=colcount+1;
    presentPFdepth(colcount)=zqPF(length(zqPF));
    presentFFdepth(colcount)=zqFF(length(zqFF));
    DtoDF(colcount)=Dcol9;
   
    no_points=numel(zqFF);
    tt_shifted=tt-1/2*t_intv;
    tt_shift(1:(no_points-1))=tt_shifted(2:no_points);
    tt_shift_BP=flip(tt_shift);
    tt_flip=flip(tt);

    for ii=1:(no_points-1)
        aggr_rateFF(ii)=(zqFF(ii+1)-zqFF(ii))/t_intv;
    end

    for ii=1:(no_points-1)
        aggr_rate_1(ii)=(zq_1(ii+1)-zq_1(ii))/t_intv;
    end

    for ii=1:(no_points-1)
        aggr_ratePF(ii)=(zqPF(ii+1)-zqPF(ii))/t_intv;
    end

    presentAggrRate(colcount)=aggr_rate_1(length(aggr_rate_1));
    
    figure(97)       %%% Freezing front and PF depth
    plot(tt_flip,zqFF,'-.')
    hold on
    plot(tt_flip,zqPF,'-.')
    title('Freezing front and permafrost depth - Zone 8-9')
    set(gca,'YDir', 'reverse')
    xlabel('Time [yr BP]')
    ylabel('Depth [m]')
    movegui('center')

    figure(96)         %%% Aggradation rate
    % hold on
    % plot(tt_shift_BP,aggr_rateFF,'b--')
    hold on
    plot(tt_shift_BP,aggr_rate_1,'k')
    % hold on
    % plot(tt_shift_BP,aggr_ratePF,'r--')
    xlabel('Time [yr BP]')
    ylabel('Rate [m/yr]')
    title('Aggradation rate - Zone 8-9')
    movegui('south')
    grid on
    grid minor
end
if any(col_incl==10)
    l=9500;     %%% set time [yr]
    k=0;
    zqFF=zeros(1,l);
    zq_1=zeros(1,l);
    zqPF=zeros(1,l);
    temp=zeros(nocell,l);
    tt=1:t_intv:l;
    tt_shift=1:t_intv:(l-1);
    aggr_rateFF=zeros(1,(l-1));
    aggr_rate_1=aggr_rateFF;
    aggr_ratePF=aggr_rateFF;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_9_10a(1:nocell,k);       %%% use right T_matx...     

        zqFF(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),TqFF);
        tt(i)=t_intv*k;
    end

    k=0;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_9_10a(1:nocell,k);

        zq_1(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),Tq_1);
        tt(i)=t_intv*k;
    end

    k=0;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_9_10a(1:nocell,k);

        zqPF(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),TqPF);
    %     ttPF(i)=t_intv*k;
    end
    
    colcount=colcount+1;
    presentPFdepth(colcount)=zqPF(length(zqPF));
    presentFFdepth(colcount)=zqFF(length(zqFF));
    DtoDF(colcount)=Dcol10;
   
    no_points=numel(zqFF);
    tt_shifted=tt-1/2*t_intv;
    tt_shift(1:(no_points-1))=tt_shifted(2:no_points);
    tt_shift_BP=flip(tt_shift);
    tt_flip=flip(tt);

    for ii=1:(no_points-1)
        aggr_rateFF(ii)=(zqFF(ii+1)-zqFF(ii))/t_intv;
    end

    for ii=1:(no_points-1)
        aggr_rate_1(ii)=(zq_1(ii+1)-zq_1(ii))/t_intv;
    end

    for ii=1:(no_points-1)
        aggr_ratePF(ii)=(zqPF(ii+1)-zqPF(ii))/t_intv;
    end
    
    presentAggrRate(colcount)=aggr_rate_1(length(aggr_rate_1));
    
    figure(107)       %%% Freezing front and PF depth
    plot(tt_flip,zqFF,'-.')
    hold on
    plot(tt_flip,zqPF,'-.')
    title('Freezing front and permafrost depth - Zone 9-10a')
    set(gca,'YDir', 'reverse')
    xlabel('Time [yr BP]')
    ylabel('Depth [m]')
    movegui('center')

    figure(106)         %%% Aggradation rate
    % hold on
    % plot(tt_shift_BP,aggr_rateFF,'b--')
    hold on
    plot(tt_shift_BP,aggr_rate_1,'k')
    % hold on
    % plot(tt_shift_BP,aggr_ratePF,'r--')
    xlabel('Time [yr BP]')
    ylabel('Rate [m/yr]')
    title('Aggradation rate - Zone 9-10a')
    movegui('south')
    grid on
    grid minor
end
if any(col_incl==11)
    l=9500;     %%% set time [yr]
    k=0;
    zqFF=zeros(1,l);
    zq_1=zeros(1,l);
    zqPF=zeros(1,l);
    temp=zeros(nocell,l);
    tt=1:t_intv:l;
    tt_shift=1:t_intv:(l-1);
    aggr_rateFF=zeros(1,(l-1));
    aggr_rate_1=aggr_rateFF;
    aggr_ratePF=aggr_rateFF;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_9_10b(1:nocell,k);       %%% use right T_matx...     

        zqFF(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),TqFF);
        tt(i)=t_intv*k;
    end

    k=0;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_9_10b(1:nocell,k);

        zq_1(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),Tq_1);
        tt(i)=t_intv*k;
    end

    k=0;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_9_10b(1:nocell,k);

        zqPF(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),TqPF);
    %     ttPF(i)=t_intv*k;
    end
    
    colcount=colcount+1;
    presentPFdepth(colcount)=zqPF(length(zqPF));
    presentFFdepth(colcount)=zqFF(length(zqFF));
    DtoDF(colcount)=Dcol11;
   
    no_points=numel(zqFF);
    tt_shifted=tt-1/2*t_intv;
    tt_shift(1:(no_points-1))=tt_shifted(2:no_points);
    tt_shift_BP=flip(tt_shift);
    tt_flip=flip(tt);

    for ii=1:(no_points-1)
        aggr_rateFF(ii)=(zqFF(ii+1)-zqFF(ii))/t_intv;
    end

    for ii=1:(no_points-1)
        aggr_rate_1(ii)=(zq_1(ii+1)-zq_1(ii))/t_intv;
    end

    for ii=1:(no_points-1)
        aggr_ratePF(ii)=(zqPF(ii+1)-zqPF(ii))/t_intv;
    end
    
    presentAggrRate(colcount)=aggr_rate_1(length(aggr_rate_1));
    
    figure(117)       %%% Freezing front and PF depth vs time
    plot(tt_flip,zqFF,'-.')
    hold on
    plot(tt_flip,zqPF,'-.')
    title('Freezing front and permafrost depth - Zone 9-10b')
    set(gca,'YDir', 'reverse')
    xlabel('Time [yr BP]')
    ylabel('Depth [m]')
    movegui('center')

    figure(116)         %%% Aggradation rate
    % hold on
    % plot(tt_shift_BP,aggr_rateFF,'b--')
    hold on
    plot(tt_shift_BP,aggr_rate_1,'k')
    % hold on
    % plot(tt_shift_BP,aggr_ratePF,'r--')
    xlabel('Time [yr BP]')
    ylabel('Rate [m/yr]')
    title('Aggradation rate - Zone 9-10b')
    movegui('south')
    grid on
    grid minor
end
if any(col_incl==12)
    l=10000;     %%% set time [yr]
    k=0;
    zqFF=zeros(1,l);
    zq_1=zeros(1,l);
    zqPF=zeros(1,l);
    temp=zeros(nocell,l);
    tt=1:t_intv:l;
    tt_shift=1:t_intv:(l-1);
    aggr_rateFF=zeros(1,(l-1));
    aggr_rate_1=aggr_rateFF;
    aggr_ratePF=aggr_rateFF;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_10(1:nocell,k);       %%% use right T_matx...     

        zqFF(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),TqFF);
        tt(i)=t_intv*k;
    end

    k=0;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_10(1:nocell,k);

        zq_1(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),Tq_1);
        tt(i)=t_intv*k;
    end

    k=0;
    for i=1:l
        k=k+1;
        temp(1:nocell,k)=T_matx_10(1:nocell,k);

        zqPF(i)=interp1(transpose(temp(1:nocell_1,k)),z(1:nocell_1),TqPF);
    %     ttPF(i)=t_intv*k;
    end
    
    colcount=colcount+1;
    presentPFdepth(colcount)=zqPF(length(zqPF));
    presentFFdepth(colcount)=zqFF(length(zqFF));
    DtoDF(colcount)=Dcol12;
   
    no_points=numel(zqFF);
    tt_shifted=tt-1/2*t_intv;
    tt_shift(1:(no_points-1))=tt_shifted(2:no_points);
    tt_shift_BP=flip(tt_shift);
    tt_flip=flip(tt);

    for ii=1:(no_points-1)
        aggr_rateFF(ii)=(zqFF(ii+1)-zqFF(ii))/t_intv;
    end

    for ii=1:(no_points-1)
        aggr_rate_1(ii)=(zq_1(ii+1)-zq_1(ii))/t_intv;
    end

    for ii=1:(no_points-1)
        aggr_ratePF(ii)=(zqPF(ii+1)-zqPF(ii))/t_intv;
    end
    
    presentAggrRate(colcount)=aggr_rate_1(length(aggr_rate_1));
    
    figure(127)       %%% Freezing front and PF depth vs time
    plot(tt_flip,zqFF,'-.')
    hold on
    plot(tt_flip,zqPF,'-.')
    title('Freezing front and permafrost depth - Zone 10')
    set(gca,'YDir', 'reverse')
    xlabel('Time [yr BP]')
    ylabel('Depth [m]')
    movegui('center')

    figure(126)         %%% Aggradation rate
    % hold on
    % plot(tt_shift_BP,aggr_rateFF,'b--')
    hold on
    plot(tt_shift_BP,aggr_rate_1,'k')
    % hold on
    % plot(tt_shift_BP,aggr_ratePF,'r--')
    xlabel('Time [yr BP]')
    ylabel('Rate [m/yr]')
    title('Aggradation rate - Zone 10')
    movegui('south')
    grid on
    grid minor
end
%% Final FF and PF depth

figure(18)
plot(DtoDF,presentFFdepth,'k--+', 'MarkerEdgeColor','blue')
hold on
plot(DtoDF,presentPFdepth,'k--+', 'MarkerEdgeColor', 'red')
xlabel('Distance to delta front [m]')
ylabel('Depth [m]')
title('Permafrost and freezing front depths')
legend('Freezing front','Permafrost')
set(gca,'YDir','reverse')
movegui('southeast')

%% Final aggradation rate and recharge equivalent
figure(19)
plot(DtoDF(2:end),presentAggrRate(2:end),'k+--', 'MarkerEdgeColor','blue')
xlabel('Distance to delta front [m]')
ylabel('Rate [m/yr]')
title('Aggradation rate and recharge equivalent')
legend('Aggradation rate')
movegui('east')
toc