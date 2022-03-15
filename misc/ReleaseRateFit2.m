
clear all;
close all; 


%%  Fit data
load CaM.txt
load CaN.txt
load DualSensor_RelRate.txt
load Allosteric_RelRate.txt
load Merged_Data.csv
load time.txt
load ReleaseRate_Goda2.csv
load Merged_Data.csv

Cac = CaM(:);
Can = CaN(:);
tspan = time(:);
tspan_goda = ReleaseRate_Goda2(:,1);
Goda_data = ReleaseRate_Goda2(:,2);



%% Fitting to data


%% Goda Data --- Exponential Fit

a10 = [0.025 6.0 0.00023 160 0.00012];              % Goda Data
lb = [0 0 0 0 0];
ub = [];
[a1,resnorm1,residual1] = lsqcurvefit(@(a1,tspan)myode_for(a1,tspan,"Goda"),a10,tspan_goda,Goda_data,lb,ub);
sprintf('%f %f %f %f %f',a1(1),a1(2),a1(3),a1(4),a1(5))
Goda = myode_for(a1, tspan_goda,"Goda");



%% Allosteric Model --- Exponential Fit


[Allo_max, Ii] = max(DualSensor_RelRate(:));
AllostericModel_Fast = Allosteric_RelRate(Ii:end);
tspan_allo = time(Ii:end);

%AllostericModel_Fast = Merged_Data(:,4);
%tspan_allo = Merged_Data(:,1);
                       
a30 = [0.025 6.0 0.00023 160 0.00012];              % Modified Allosteric Model
lb = [0 0 0 0 0];
ub = [];
[a3,resnorm3,residual3] = lsqcurvefit(@(a3,tspan_allo)myode_for(a3,tspan_allo,"Allosteric"),...
                                      a30,tspan_allo,AllostericModel_Fast,lb,ub);
sprintf('%f %f %f %f %f',a3(1),a3(2),a3(3),a3(4),a3(5))
AllostericModel = myode_for(a3,tspan_allo,"Allosteric");



%% Dual Sensor --- Exponential Fit


[Dual_max, I] = max(DualSensor_RelRate(:));
DualSensorModel_Fast = DualSensor_RelRate(I:end);
tspan_dual = time(I:end);

%DualSensorModel_Fast = Merged_Data(:,3);
%tspan_dual = Merged_Data(:,1);


%a20 = [0.01 0.7 0.0009 7 0.00005 160 0.000015]; % Dual Sensor Model
a20 = [0.025 6.0 0.00023 160 0.00012];     % Dual Sensor Model
lb = [0 0 0 0 0 0 0];
ub = [];
[a2,resnorm2,residual2] = lsqcurvefit(@(a2,tspan_dual)myode_for(a2,tspan_dual,"DualSensor"),...
                                        a20,tspan_dual,DualSensorModel_Fast,lb,ub);
sprintf('%f %f %f %f %f',a2(1),a2(2),a2(3),a2(4),a2(5))
DualSensorModel = myode_for(a2, tspan_dual,"DualSensor");





%% Estimate release rate



%% Plot Figures

figure
semilogy(tspan_goda, Goda_data, 'k.', 'LineWidth', 2, 'MarkerSize', 10)
hold on
semilogy(tspan_goda, Goda, 'LineWidth', 2, 'MarkerSize', 10)
hold on
legend({'Original Data From Goda', 'Double Exponential Fit'},'Location', 'northwest')

%Get Current Figure (GCF) & Set image size before saving image
width = 20;  % cm 
height = 19; % cm
set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])

%Set the resolution of 1000dpi and save the plot in TIFF format 
print -dpng -r1000 Exponential_Fit_To_Goda
hold off


figure
semilogy(tspan_allo, AllostericModel_Fast, 'k.', 'LineWidth', 2, 'MarkerSize', 10)
hold on
semilogy(tspan_allo, AllostericModel, 'LineWidth', 2, 'MarkerSize', 10)
hold on
legend({'Single AP Allosteric Data', 'Double Exponential Fit'},'Location', 'northwest')

%Get Current Figure (GCF) & Set image size before saving image
width = 20;  % cm 
height = 19; % cm
set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])

%Set the resolution of 1000dpi and save the plot in TIFF format 
print -dpng -r1000 Exponential_Fit_To_Allosteric
hold off


figure
semilogy(tspan_dual, DualSensorModel_Fast, 'k.', 'LineWidth', 2, 'MarkerSize', 10)
hold on
semilogy(tspan_dual, DualSensorModel, 'LineWidth', 2, 'MarkerSize', 10)
hold on
legend({'Single AP Dual-Sensor Data','Double Exponential Fit'},'Location', 'northwest')

%Get Current Figure (GCF) & Set image size before saving image
width = 20;  % cm 
height = 19; % cm
set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])

%Set the resolution of 1000dpi and save the plot in TIFF format 
print -dpng -r1000 Exponential_Fit_To_DualSensor
hold off


%% FUNCTIONS


function y_real = myode_for(a, tspan, text)
        
        %% Double Exponential Fit
        
        if text=="text"
            yD = a(1)*exp(-tspan/a(2)) + a(3)*exp(-tspan/a(4)) + a(5)*exp(-tspan/a(6)) + a(7);
        else
            yD = a(1)*exp(-tspan/a(2)) + a(3)*exp(-tspan/a(4)) + a(5);
        end
        y_real= yD; 
    
end 