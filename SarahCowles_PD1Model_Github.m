%SCC220322_PD1model
close all
clc

%% Run Model Loop for varying Kds
% Run Model: Loop for varying Kds, Calculate time of PD-1 occupancy above threshold

%For Figure 3b: KdVec = 500, 50, 5, 0.5 nM
% kdVec = [1e5, 5e-2; 1e5, 5e-3; 1e5, 5e-4; 1e5, 5e-5];
% kdVec = kdVec';

%For Figures 4 and 5
%set Kd range
kon = ones(1,100).*1e5; %Adjust vector length to vary plot resolution
koff = logspace(-6,-2,100); %Adjust vector length to vary plot resolution, must match above length
kdVec = [kon; koff];
doseConcVec = [500e-9]; %M %Adjust dose (666 nM = 135 ug, 250 = 50 ug, 500 = 100 ug, 3000 = 600 ug)
%variable = 0.0445/3600; %internalization rate, for sensitivity analysis
variableVec1 = 0.0445/3600; %internalization rate
%kFast = 0.1112/3600; %1/s %curve fitting determined 20% weight towards fast rate
%kSlow = 0.007071/3600; %1/s %curve fitting determined 80% weight towards fast rate
%variable = kFast;
%variable2 = kSlow;
%variableVec1 = [variable*0.8, variable*0.9, variable, variable*1.1, variable*1.2];
%variableVec2 = [variable2*0.8, variable2*0.9, variable2, variable2*1.1, variable2*1.2];
Time95Occ = zeros(length(kdVec),length(doseConcVec)); %initialize vector
for k = 1:length(variableVec1)
    k_int = variableVec1(k);
    %k_fast = variableVec1(k);
    %k_slow = variableVec2(k);
    for j = 1:length(doseConcVec)
        doseConc = doseConcVec(j); %1/s
        kFast = 0.1112/3600; %1/s %curve fitting determined 20% weight towards fast rate
        kSlow = 0.007071/3600; %1/s %curve fitting determined 80% weight towards fast rate
        %k_int = 0.0445/3600; %1/s (Median internalization rate from experiemental data)
        PD1ssConc = 10000*1000*4.3/(6.022*10^23)*10^6; %Expected steady state concentration of receptor in tumor
        Occ = []; %initialize occupancy vector
        for i = 1:length(kdVec) %run model for each Kd value
            doseNum = 1; %Number of doses given
            doseDelay = 20*3600*24; %Time between doses (seconds)
            [Tsol,Xsol] = solveModel(doseConc,doseNum,doseDelay,kdVec(1,i),kdVec(2,i),k_int,PD1ssConc,kFast,kSlow); %run model
            Occ = Xsol(:,3)./(Xsol(:,3)+Xsol(:,4)); %collect occupancy of PD-1 at each timepoint
            idx = find(Occ > 0.99); %find index of Xsol where occupancy > 0.99
            if isempty(idx)
                Time95Occ(i,j,k)=0; %total time = 0 if occupancy never reaches 0.99
            else
                Occ95 = Tsol(idx); %find time points where occupancy > 0.99
                Time95Occ(i,j,k) = max(Occ95)-min(Occ95); %total time when occupancy > 0.99         
            end
            %COPY PLOTTING CODE HERE IF DESIRED
        end
    end
end

%% Plot Time of PD-1 Occupancy Above Threshold for Given Affinity
figure
Time95Occ = Time95Occ./3600./24; %Adjust timescale to days
KD = kdVec(2,:)./kdVec(1,:); %Calculate KD from kon and koff
SensitivityAnalysisColors = [0.114 0.267 0.722
             0.514 0.616 0.922
             0 0 0
             0.926 0.529 0.549
             0.808 0.126 0.161];
colororder(SensitivityAnalysisColors)
plot(KD,Time95Occ(:,1),'linewidth',3) %single curve
%plot(KD,Time95Occ(:,1),KD,Time95Occ(:,2),'linewidth',3) %multiple curves
%plot(KD,Time95Occ(:,1,1),KD,Time95Occ(:,1,2),KD,Time95Occ(:,1,3),KD,Time95Occ(:,1,4),KD,Time95Occ(:,1,5),'linewidth',3) %sensitivity analysis
ylabel({'Time of PD-1 Occupancy';'above 99% in Tumor [Days]'})
xlabel('mAb Kd [nM]')
xlim([10^-11 10^-6])
ylim([0 10])
set(gca,'FontSize', 20,'XScale', 'log')
xticks([10^-11 10^-10 10^-9 10^-8 10^-7 10^-6])
xticklabels({'0.01','0.1','1','10','100','1000'})
%legend('k_{fast}-20%, k_{slow}-20%','k_{fast}-10%, k_{slow}-10%','k_{fast}, k_{slow}','k_{fast}+10%, k_{slow}+10%','k_{fast}+20%, k_{slow}+20%')
%legend('k_{int}-20%','k_{int}-10%','k_{int}','k_{int}+10%','k_{int}+20%')
    

%% Copy this plotting code for specific KD values into if statement aboves
%         newcolors = [.753, .753, 1
%             .451, .451, 1
%             0, 0, .847
%             0, 0, .545];
%         colororder(newcolors)
%         subplot(2,2,1)
%         %Plot Drug in Plasma
%         plot(Tsol./3600./24,Xsol(:,1)./1e-9','Linewidth',1.5)
%         ylabel('[mAb] in Plasma [nM]')
%         xlabel('Time [Days]')
%         set(gca,'FontSize', 16)
%         hold on 
%         subplot(2,2,2)
%         %Plot Free Drug in the Tumor
%         plot(Tsol./3600./24,Xsol(:,2)./1e-9','Linewidth',1.5)
%         ylabel('[mAb] Tumor, Free [nM]')
%         xlabel('Time [Days]')
%         set(gca,'FontSize', 16)
%         hold on 
%         subplot(2,2,3)
%         %Plot Bound Drug Concentration in Tumor (also complex concentration, or bound receptor concentration)
%         plot(Tsol./3600./24,Xsol(:,3)./1e-9','Linewidth',1.5)
%         ylabel('[PD-1/mAb] Tumor, Bound [nM]')
%         xlabel('Time [Days]')
%         set(gca,'FontSize', 16)
%         hold on 
%         subplot(2,2,4)
%         %Plot Free PD-1 (Ag) Receptor Concentration in Tumor
%         plot(Tsol./3600./24,Xsol(:,4)./1e-9,'Linewidth',1.5)
%         ylabel('[PD-1] Tumor, Free [nM]')
%         xlabel('Time [Days]')
%         set(gca,'FontSize', 16)
%         hold on

%% Model Functions

%Shell function to run model
function [Tsolfinal,Xsolfinal] = solveModel(doseConc,doseNum,doseDelay,k_on, k_off,k_int,PD1ssConc,kfast,kslow)
    options = odeset('RelTol',5e-14,'AbsTol',1e-14, 'MaxStep', 1800); %set tolerance a few orders of magnitude below concentrations
    eps = 0.24; %tumor void fraction (also can think about it as accessible soluble volume fraction of tumor)
    PDL1ssConc = 3000*1000*10^-12/33000; %M
    PD1free0 = 0.05*PD1ssConc*0.1; %M %Start with 10% total PD-1 steady state concentration and 95% is bound to PD-L1
    PD1bound0 = 0.95*PD1ssConc*0.1; %M 
    PDL1free0 = PDL1ssConc - PD1bound0; %M
    NewIntVal = [doseConc,0,0,PD1free0,PD1bound0,PDL1free0]; %Set initial values for ODE function 
        %NewIntVal(1) = Conc antibody in plasma at time zero, M
        %NewIntVal(2) = Conc antibody free in tumor plasma at time zero, M
        %NewIntVal(3) = Conc antibody bound in tumor at time zero, M
        %NewIntVal(4) = Conc of free receptor at time zero, M
        %NewIntVal(5) = Conc of PD-L1 bound to PD-1 in the tumor at time zero, M
        %NewIntVal(6) = Conc of PD-L1 free in the tumor at time zero, M
    Xsolfinal = []; %Initialize solution matrix
    Tsolfinal = []; %Initialize solution matrix
    for i = 1:doseNum %iterate through doses
        [Tsol, Xsol] = ode15s(@(t,X)odefunc(t,X,k_on,k_off,k_int,PD1ssConc,kfast,kslow),[((i-1)*doseDelay) (i*doseDelay)],NewIntVal,options); %solve ODE for dose time interval
        NewIntVal = Xsol(end,:);
        NewIntVal(1) = NewIntVal(1)+doseConc; %adjust next starting plasma concentration by new dose
        Xsolfinal = [Xsolfinal;Xsol]; %concatenate solution matrices
        Tsolfinal = [Tsolfinal;Tsol]; %concatenate solution matrices
    end
    Xsolfinal(:,2:4) = Xsolfinal(:,2:4)./eps; %Scale concentration matrix terms in tumor by void fraction to show local tumor concentration
end

%Define Model ODE Function
 function derivX = odefunc(t, X, k_on, k_off, k_int, PD1ssConc,kFast,kSlow)
       Ab_p = X(1); %Conc antibody in plasma, M
       Ab_tf = X(2); %Conc antibody free in tumor plasma, M
       Ab_tb = X(3); %Conc antibody bound in tumor, M
       PD1 = X(4); %Conc of free receptor, M
       PDL1_tb = X(5); %Conc of bound PD-L1, M
       PDL1_tf = X(6); %Conc of free PD-L1, M
       P = 3e-9; %m/s
       Rcap = 8e-6; %capillary radius m
       Rk = 75e-6; %Krogh cylinder radius m
       eps = 0.24; %tumor void fraction
       Vtumor = 100e-6; %L
       Vplasma = 2e-3; %L
       %Biexponential clearance rate: 
       kFast = (0.2)*kFast; %1/s %curve fitting determined 20% weight towards fast rate
       kSlow = (0.8)*kSlow; %1/s %curve fitting determined 80% weight towards fast rate
       %PD-L1 rates
       KD = 4000*10^-9; %M
       k_PDL1off = 0.5; %1/s
       k_PDL1on = k_PDL1off/KD; %1/M*s
       k_PDL1int = 0.0922/3600; %1/s
       PDL1ssConc = 3000*1000*10^-12/33000; %M steady state PD-L1 concentration 
       k_PDL1syn = k_PDL1int*PDL1ssConc; %Pseudo-Steady State (Syn ~ Int)
       
       %ODE Model Set
       derivX = [-kSlow*Ab_p - kFast*Ab_p + 2*P*Rcap/Rk^2 * (Ab_tf/eps - Ab_p) * (Vtumor/Vplasma*eps);...
              2*P*Rcap/Rk^2*(Ab_p*eps - Ab_tf) - k_on*Ab_tf*PD1/eps + k_off*Ab_tb;...
              k_on*Ab_tf*PD1/eps - k_off*Ab_tb - k_int*Ab_tb;...
              -k_on*Ab_tf*PD1/eps + (PD1+Ab_tb)*(1-((PD1+Ab_tb)/(PD1ssConc))) + k_off*Ab_tb - k_int*PD1 - k_PDL1on*PDL1_tf*PD1/eps + k_PDL1off*PDL1_tb;...
              k_PDL1on*PDL1_tf*PD1/eps - k_PDL1off*PDL1_tb;...
              -k_PDL1on*PDL1_tf*PD1/eps + k_PDL1off*PDL1_tb - k_PDL1int*PDL1_tf + k_PDL1syn];
 end	
