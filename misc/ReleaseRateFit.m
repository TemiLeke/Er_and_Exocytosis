
close all; clear all;

%%  Fit data

load Merged_Data.csv
load ReleaseRate_VGCC.csv
load ReleaseRate_VGCC.csv
ydata_model = Merged_Data(:,3);
ydata_exp = Merged_Data(:,2);
tspan = Merged_Data(:,1);
Cac = Merged_Data(:,4);
Can = Merged_Data(:,5);

%% Fitting to data


%% Modified Dual Sensor

%{
a0 = [0.061200 2.320000 0.003814 0.014682 0.025008 0.250006 2.000008 0.000002 6.340000 0.000050 0.002200 0.027991 0.005208 0.000150 0.001093]; % Dual Sensor Modified
lb = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
ub = [];
[a,resnorm,residual]= lsqcurvefit(@(a,tspan)myode_for(a,tspan,Cac,Can),a0,tspan,ydata_exp,lb,ub);
sprintf('%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f',a(1),a(2),a(3),a(4),a(5),a(6),a(7),a(8),a(9),a(10),a(11),a(12),a(13),a(14),a(15))
%sprintf('%f %f', resnorm, residual)



%% Estimate release rate


ModifiedDualSensorModel = myode_for(a, tspan, Cac, Can);


%% Plot Figures


figure
plot(tspan, Merged_Data(:,2), 'LineWidth', 3, 'MarkerSize', 20)
hold on
plot(tspan, Merged_Data(:,3), 'LineWidth', 3, 'MarkerSize', 20)
hold on
legend({'Original Data from Goda', 'DualSensor model'},'Location', 'northeast')
title('Release Rates from Dual-Sensor and Goda','FontSize',12,'FontWeight','bold','Color','k')

%Get Current Figure (GCF) & Set image size before saving image
width = 20;  % cm 
height = 19; % cm
set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])

%Set the resolution of 1000dpi and save the plot in TIFF format 
print -dpng -r1000 DualSensorModel_and_Goda
hold off

figure
semilogy(tspan, Merged_Data(:,2), 'LineWidth', 3, 'MarkerSize', 20)
hold on
semilogy(tspan, Merged_Data(:,3), 'LineWidth', 3, 'MarkerSize', 20)
hold on
legend({'Original Data from Goda', 'DualSensor model'},'Location', 'northeast')
title('Log Release Rates from Dual-Sensor and Goda','FontSize',12,'FontWeight','bold','Color','k')

%Get Current Figure (GCF) & Set image size before saving image
width = 20;  % cm 
height = 19; % cm
set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])

%Set the resolution of 1000dpi and save the plot in TIFF format 
print -dpng -r1000 DualSensorModel_and_Goda_log
hold off


figure
semilogy(tspan, ModifiedDualSensorModel(:,1),"r", 'LineWidth', 5, 'MarkerSize', 20)
hold on
semilogy(tspan, ydata_exp, 'k.', 'LineWidth', 5, 'MarkerSize', 20)
hold on
legend({' Fit from Dual-Sensor Model', 'Original Data'},'Location', 'northeast')
title('Fitting Dual-Sensor Model to Release Rate from Goda','FontSize',12,'FontWeight','bold','Color','k')

%Get Current Figure (GCF) & Set image size before saving image
width = 20;  % cm 
height = 19; % cm
set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])

%Set the resolution of 1000dpi and save the plot in TIFF format 
print -dpng -r1000 DualSensorModel_Fit_To_Goda
hold off
%}


%% Modified Allosteric Model

%{
a0 = [1.00e-01 4.00 2.0e-07 31.3 0.5 6.34 5.0e-05 0.0022 0.028 0.0032 0.00015 0.00051];              % Dual Sensor Modified
lb = [0 0 0 0 0 0 0 0 0 0 0 0];
ub = [];
[a,resnorm,residual]= lsqcurvefit(@(a,tspan)myode_for(a,tspan,Cac,Can),a0,tspan,ydata_exp,lb,ub);

sprintf('%f %f %f %f %f %f %f %f %f %f %f %f',a(1),a(2),a(3),a(4),a(5),a(6),a(7),a(8),a(9),a(10),a(11),a(12))
%sprintf('%f %f', resnorm, residual)
%}


%% Estimate release rate
%{
ModifiedAllostericModel = myode_for(a, tspan, Cac, Can);
%}

%% Plot Figures

%{
figure
plot(tspan, Merged_Data(:,2), 'LineWidth', 3, 'MarkerSize', 20)
hold on
plot(tspan, Merged_Data(:,6), 'LineWidth', 3, 'MarkerSize', 20)
hold on
legend({'Original Data from Goda', 'Allosteric model'},'Location', 'northeast')
title('Release Rates from Allosteric and Goda','FontSize',12,'FontWeight','bold','Color','k')

%Get Current Figure (GCF) & Set image size before saving image
width = 20;  % cm 
height = 19; % cm
set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])

%Set the resolution of 1000dpi and save the plot in TIFF format 
print -dpng -r1000 AllostericModel_and_Goda
hold off

figure
semilogy(tspan, Merged_Data(:,2), 'LineWidth', 3, 'MarkerSize', 20)
hold on
semilogy(tspan, Merged_Data(:,6), 'LineWidth', 3, 'MarkerSize', 20)
hold on
legend({'Original Data from Goda', 'Allosteric model'},'Location', 'northeast')
title('Log Release Rates from Allosteric and Goda','FontSize',12,'FontWeight','bold','Color','k')

%Get Current Figure (GCF) & Set image size before saving image
width = 20;  % cm 
height = 19; % cm
set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])

%Set the resolution of 1000dpi and save the plot in TIFF format 
print -dpng -r1000 AllostericModel_and_Goda_log
hold off
%}

%{
figure
semilogy(tspan, ModifiedAllostericModel(:,1),"r", 'LineWidth', 5, 'MarkerSize', 20)
hold on
semilogy(tspan, ydata_exp, 'k.', 'LineWidth', 5, 'MarkerSize', 20)
hold on
legend({' Fit from Allosteric Model', 'Original Data'},'Location', 'northeast')
title('Fitting Allosteric Model to Release Rate from Goda','FontSize',12,'FontWeight','bold','Color','k')

%Get Current Figure (GCF) & Set image size before saving image
width = 20;  % cm 
height = 19; % cm
set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])

%Set the resolution of 1000dpi and save the plot in TIFF format 
print -dpng -r1000 AllostericModel_Fit_To_Goda
hold off
%}

%% Unpriming Model


a0 = [153e-03 3.5 0.25 0.000417e-03 115e-03 0.66e-03 0.6e-03 3.82e-03 60e-03 0.25 0.000417e-03 8e-03 0.16e-03 0.052e-03 0.615e-03 0.75e-03];
lb = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
ub = [];
[a,resnorm3, residual3]= lsqcurvefit(@(a,tspan)myode_for(a,tspan,Cac,Can),a0,tspan,ydata_exp,lb,ub);
sprintf('%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', ...
        a(1),a(2),a(3),a(4),a(5),a(6),a(7),a(8),a(9),a(10),a(11),a(12),a(13),a(14),a(15),a(16))
%sprintf('%f %f', resnorm, residual)



%% Estimate release rate

KR_FFModelresult = myode_for(a, tspan, Cac, Can);


%% Plot Figures


figure
semilogy(tspan, KR_FFModelresult(:,1),"r", 'LineWidth', 5, 'MarkerSize', 20)
hold on
semilogy(tspan, ydata_exp, 'k.', 'LineWidth', 5, 'MarkerSize', 20)
hold on
legend({' Fit from Kiss-and-run and Full Fusion Model', 'Original Data'},'Location', 'northeast')
title('Fitting Unprimed Model to Release Rate from Goda','FontSize',12,'FontWeight','bold','Color','k')

%Get Current Figure (GCF) & Set image size before saving image
width = 20;  % cm 
height = 19; % cm
set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])

%Set the resolution of 1000dpi and save the plot in TIFF format 
print -dpng -r1000 KissAndRun_Fit_To_Goda
hold off


%% FUNCTIONS


function y_real = myode_for(a, tspan, Cac, Can)

    dt = 1e-03;  
    %tspan = 0:dt:1;

    yA = zeros(length(Cac),1);
    yD = zeros(length(Cac),1);
    yU = zeros(length(Cac),1);
    
    yMA = zeros(length(tspan),1); 
    yMA1 = zeros(length(tspan),1);
    yMA2 = zeros(length(tspan),1); 
    yMD = zeros(length(tspan),1); 
    yMD1 = zeros(length(tspan),1);
    yMD2 = zeros(length(tspan),1); 
    yMU = zeros(length(tspan),1); 
    yMU1 = zeros(length(tspan),1); 
    t_real = zeros(length(Cac), 1);
    
    nn = 100;
    
    %% Dual Sensor Model 
    
    y5 = [3 0 0 0 0 0 0 0 0 0 0 0 ...
          0 0 0 0 0 0 3 0 0 0 0 0 0 ...
          0 0 0 0 0 0 0 0 0 0 0 85 10 ...
          0 0];
    
    %% Allosteric Model
    
    y1 = [3 0 0 0 0 0 3 0 0 0 0 0 85 12 0 0];
    
    %% Unpriming Model
    
    y2 = [160 0 0 0 0 0 0 ... 
          40 0 0 0 0 0 0 0 0 0];
    
    %% Loop
    
    for i=2:length(Cac)

        %% Modified Allosteric Model
        
        %{
        I = a(3);           % /ms
        F = a(4);
            
        %[t1, y1] = ode15s(@(t1,y1)ModifiedAllosteric(t1, y1, a, Cac(i), Can(i)),tspan,y0MD);

        yn = y1(i-1,:);   % variable values at previous time step  

        % Runge-Kutta Method
        for j = 1:nn  % this loop is used so that we don't have to save data at each time step

            k1 = ModifiedAllosteric(tspan(i), yn, a, Cac(i), Can(i));                        
            k2 = ModifiedAllosteric(tspan(i), yn + k1*dt/2, a, Cac(i), Can(i));
            k3 = ModifiedAllosteric(tspan(i), yn + k2*dt/2, a, Cac(i), Can(i));
            k4 = ModifiedAllosteric(tspan(i), yn + k3*dt, a, Cac(i), Can(i));
            yn = yn + (k1 + 2*k2+ 2*k3 + k4)* dt/6;
        end

       y1(i,:) = yn;           % update the new value of all variables; 
       %}
        
        %% Modified Dual Sensor Model 
        
        %{
        a_n = a(5);
        gam_2 = a(7);         % /ms Synchronous release rate
        gam_3 = a_n*gam_2;    % /ms Asynchronous release rate
        gam_1 = a(8);         % /ms Spontaneous release rate


        %[t5, y5] = ode15s(@(t5,y5)ModifiedDualSensor(t5, y5, a, Cac(i), Can(i)),tspan,y0MD);

        yn = y5(i-1,:);   % variable values at previous time step  

        % Runge-Kutta Method
        for j = 1:nn  % this loop is used so that we don't have to save data at each time step

            k1 = ModifiedDualSensor(tspan(i), yn, a, Cac(i), Can(i));                        
            k2 = ModifiedDualSensor(tspan(i), yn + k1*dt/2, a, Cac(i), Can(i));
            k3 = ModifiedDualSensor(tspan(i), yn + k2*dt/2, a, Cac(i), Can(i));
            k4 = ModifiedDualSensor(tspan(i), yn + k3*dt, a, Cac(i), Can(i));
            yn = yn + (k1 + 2*k2+ 2*k3 + k4)* dt/6;
        end

       y5(i,:) = yn;           % update the new value of all variables; 
       
       %}
       
        %% Kiss-and-Run and Full-Fusion Model 
    
        
        Kspont4 = params(4); 
        Kevoked4 = params(5);
        Kspont7 = params(11);
        Kevoked7 = params(12);

        yn = y2(i-1,:);       % variable values at previous time step  

        % Runge-Kutta Method
        for j = 1:nn  % this loop is used so that we don't have to save data at each time step

            k1 = KR_FF(tspan(i), yn, a, Cac(i), Can(i));                        
            k2 = KR_FF(tspan(i), yn + k1*dt/2, a, Cac(i), Can(i));
            k3 = KR_FF(tspan(i), yn + k2*dt/2, a, Cac(i), Can(i));
            k4 = KR_FF(tspan(i), yn + k3*dt, a, Cac(i), Can(i));
            yn = yn + (k1 + 2*k2+ 2*k3 + k4)* dt/6;
            
        end

       y2(i,:) = yn;           % update the new value of all variables; 
        
    end    
     
    %% Allosteric Model
    
    %{
    yMA1(:,1) = I.*(y1(:,1) + F.*y1(:,2) + (F^2).*y1(:,3)* + (F^3).*y1(:,4) + (F^4).*y1(:,5) + (F^5).*y1(:,6));          
    
    yMA2(:,1) = I.*(y1(:,7) + F.*y1(:,8) + (F^2).*y1(:,9)* + (F^3).*y1(:,10) + (F^4).*y1(:,11) + (F^5).*y1(:,12));          
                
    yMA(:,1) = yMA2(:,1) + yMA1(:,1);
    %}
       
        
    %% DualSensor Model
    
    %{
    yMD1(:,1) = gam_1.*y5(:,1) + gam_2.*(y5(:,16) + y5(:,17) + y5(:,18))+...
                 gam_3.*(y5(:,3) + y5(:,6) + y5(:,9) + y5(:,12) + y5(:,15) + y5(:,18));          
    
    yMD2(:,1) = gam_1.*y5(:,19) + gam_2.*(y5(:,34) + y5(:,35) + y5(:,35))+...
                 gam_3.*(y5(:,21) + y5(:,24) + y5(:,27) + y5(:,30) + y5(:,33) + y5(:,36));          
                
    yMD(:,1) = yMD2(:,1) + yMD1(:,1);
    %}
    
    
    %% Kiss-and-Run and Full-Fusion Model
    
    yMU1(:,1) = Kspont4.*y(:,1) + Kevoked4.*y(:,3) + Kspont7.*y(:,8) + Kevoked7.*y(:,13);
             
    yMU(:,1) = yMU1(:,1);
    
    y_real(:,1) = yMU(:,1); 
    
end 


function dydt = ModifiedAllosteric(t, y, params, Cac, Can)

%{
Allosteric modulation of the presynaptic Ca21 sensor for vesicle fusion
Xuelin Lou1, Volker Scheuss1? & Ralf Schneggenburger1
Calyx Of Held
%}

    %Can = Cac;
    Kon = params(1);         % /uMms
    Koff = params(2);        % /ms
    I = params(3);           % /ms
    F = params(4);
    b = params(5);  
    Krf = 1/params(6);         % /ms rate of recovery of refractoriness
    Kmob = params(7);        % /uMms
    Kdemob = params(8);      % /ms 
    Kprime = params(9);      % /uMms
    Kupr = params(10);       % /ms
    Kattach = params(11);    % /ms
    Kdetach = params(12);    % /ms

     V0 = y(1);
     V1 = y(2);
     V2 = y(3);
     V3 = y(4);
     V4 = y(5);
     V5 = y(6); 
     
     W0 = y(7);
     W1 = y(8);
     W2 = y(9);
     W3 = y(10);
     W4 = y(11);
     W5 = y(12);
     
    R = y(13);
    U = y(14);
    RFv = y(15);
    RFw = y(16); 

    Vtotal = V0 + V1 + V2 + V3 + V4 + V5; 

    Wtotal = W0 + W1 + W2 + W3 + W4 + W5;

    dRdt = -Kmob*Cac*R + U*Kdemob;
    dUdt =  -Kdemob*U + Kmob*Cac*R - Kprime*U*Cac*(1 - RFv) + Kupr*V0;

    dRFvdt = (1 - RFv)*(I*(V0 + V1*F + V2*F^2 + V3*F^3 + V4*F^4 + V5*F^5))/Vtotal - Krf*RFv;

    dRFwdt = (1 - RFw)*(I*(W0 + W1*F + W2*F^2 + W3*F^3 + W4*F^4 + W5*F^5))/Wtotal - Krf*RFw;

    dV0 = Kprime*U*Cac*(1 - RFv) - Kupr*V0 - Kattach*V0*Can*(1 - RFw) + Kdetach*W0 + Koff*V1 - I*V0 - 5*Kon*Cac*V0; 
    dV1dt = 2*Koff*b*V2 - Koff*V1 + 5*Kon*Cac*V0 - 4*Kon*Cac*V1 - I*F*V1;
    dV2dt = 3*Koff*V3*b^2 - 2*Koff*V2*b + 4*Kon*Cac*V1 - 3*Kon*Cac*V2 - I*V2*F^2;
    dV3dt = 4*Koff*V4*b^3 - 3*Koff*V3*b^2 + 3*Kon*Cac*V2 - 2*Kon*Cac*V3 - I*V3*F^3;
    dV4dt = 5*Koff*V5*b^4 - 4*Koff*V4*b^3 + 2*Kon*Cac*V3 - Kon*Cac*V4 - I*V4*F^4;
    dV5dt = Kon*Cac*V4 - 5*Koff*V5*b^4 - I*V5*F^5;


    dW0 = Kattach*W0*Can*(1 - RFw) - Kdetach*W0 + Koff*W1 - I*W0 - 5*Kon*Can*W0; 
    dW1dt = 2*Koff*b*W2 - Koff*W1 + 5*Kon*Can*W0 - 4*Kon*Can*W1 - I*F*W1;
    dW2dt = 3*Koff*W3*b^2 - 2*Koff*W2*b + 4*Kon*Can*W1 - 3*Kon*Can*W2 - I*W2*F^2;
    dW3dt = 4*Koff*W4*b^3 - 3*Koff*W3*b^2 + 3*Kon*Can*W2 - 2*Kon*Can*W3 - I*W3*F^3;
    dW4dt = 5*Koff*W5*b^4 - 4*Koff*W4*b^3 + 2*Kon*Can*W3 - Kon*Can*W4 - I*W4*F^4;
    dW5dt = Kon*Can*W4 - 5*Koff*W5*b^4 - I*W5*F^5;
       
    
    dydt=[dV0 dV1dt dV2dt dV3dt dV4dt dV5dt dW0 dW1dt dW2dt dW3dt dW4dt dW5dt...
           dRdt dUdt dRFvdt dRFwdt];

end



function dydt_5 = ModifiedDualSensor(t, y, params, Cac, Can)

%{
A dual sensor model for neurotransmitter release in central synapse
Sun et al. 2007
Calyx Of Held
%}
    % Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    alpha = params(1);  % /uMms   Association rate for synchronous release
    beta = params(2);        % /ms   Dissociation= rate for synchronous release
    Chi = params(3);     % /uMms   Association rate for Asynchronous release
    delta = params(4);     % /ms   Dissociation rate for Asynchronous release
    a = params(5);
    bo = params(6);
    gam_2 = params(7);          % /ms Synchronous release rate
    gam_3 = a*gam_2;    % /ms Asynchronous release rate
    gam_1 = params(8);  % /ms Spontaneous release rate
    Krf = params(9);         % /ms rate of recovery of refractoriness
    Kmob = params(10);     % /uMms
    Kdemob = params(11);    % /ms 
    Kprime = params(12);     % /uMms
    Kupr = params(13);      % /ms
    Kattach = params(14);  % /ms
    Kdetach = params(15);  % /ms
    
    
    %Can = Cac;
    V00 = y(1); 
    V01 = y(2);
    V02 = y(3);
    V10 = y(4);
    V11 = y(5);
    V12 = y(6);
    V20 = y(7);
    V21 = y(8);
    V22 = y(9);
    V30 = y(10);
    V31 = y(11);
    V32 = y(12);
    V40 = y(13);
    V41 = y(14);
    V42 = y(15);
    V50 = y(16);
    V51 = y(17);
    V52 = y(18);
    
    W00 = y(19); 
    W01 = y(20);
    W02 = y(21);
    W10 = y(22);
    W11 = y(23);
    W12 = y(24);
    W20 = y(25);
    W21 = y(26);
    W22 = y(27);
    W30 = y(28);
    W31 = y(29);
    W32 = y(30);
    W40 = y(31);
    W41 = y(32);
    W42 = y(33);
    W50 = y(34);
    W51 = y(35);
    W52 = y(36);
    
    R = y(37);
    U = y(38);
    RFv = y(39);
    RFw = y(40); 
    
    Vtotal = V00 + V01 + V02 + ...
             V10 + V11 + V12 + ...
             V20 + V21 + V22 + ...
             V30 + V31 + V32 + ...
             V40 + V41 + V42 + ...
             V50 + V51 + V52; 
         
    Wtotal = W00 + W01 + W02 + ...
             W10 + W11 + W12 + ...
             W20 + W21 + W22 + ...
             W30 + W31 + W32 + ...
             W40 + W41 + W42 + ...
             W50 + W51 + W52;     

    dRdt = -Kmob*Cac*R + U*Kdemob;
    dUdt = -Kprime*U*Cac*(1 - RFv) + Kupr*V00;

    dRFvdt = (1 - RFv)*(gam_3*(V02 + V12 + V22 + V32 + V42 + V52) + ...
                        gam_2*(V50 + V51 + V52) + gam_1*V00)/Vtotal - Krf*RFv;

    dRFwdt = (1 - RFw)*(gam_3*(W02 + W12 + W22 + W32 + W42 + W52) + ...
                        gam_2*(W50 + W51 + W52) + gam_1*W00)/Wtotal - Krf*RFw;

    dV00 = Kprime*U*Cac*(1 - RFv) - Kupr*V00 - Kattach*V00*Can*(1 - RFw) + Kdetach*W00 + ...
             beta*V10 - 5*alpha*V00*Cac + delta*V01 - 2*Chi*V00*Cac - gam_1*V00;
    dV01 = 2*Chi*V00*Cac - delta*V01 + 2*bo*delta*V02 - Chi*Cac*V01 + beta*V11 - 5*alpha*Cac*V01;
    dV02 = Chi*Cac*V01 - 2*bo*delta*V02 - gam_3*V02 - 5*alpha*Cac*V02 + beta*V12;
    dV10 = 5*alpha*Cac*V00 - beta*V10 - 4*alpha*Cac*V10 + 2*bo*beta*V20 - 2*Chi*Cac*V10 + delta*V11;
    dV11 = 5*alpha*Cac*V01 - beta*V11 - 4*alpha*Cac*V11 + 2*bo*beta*V21 + 2*Chi*Cac*V10 - delta*V11 - Chi*Cac*V11 + 2*bo*delta*V12;
    dV12 = 5*alpha*Cac*V02 - beta*V12 - 4*alpha*Cac*V12 + 2*bo*beta*V22 + Chi*Cac*V11 - 2*bo*delta*V12 - gam_3*V12;
    dV20 = 4*alpha*Cac*V10 - 2*bo*beta*V20 - 3*alpha*Cac*V20 + 3*(bo^2)*beta*V30 - 2*Chi*Cac*V20 + delta*V21;
    dV21 = 4*alpha*Cac*V11 - 2*bo*beta*V21 - 3*alpha*Cac*V21 + 3*(bo^2)*beta*V31 + 2*Chi*Cac*V20 - delta*V21 - Chi*Cac*V21 + 2*bo*delta*V22;
    dV22 = 4*alpha*Cac*V12 - 2*bo*beta*V22 - 3*alpha*Cac*V22 + 3*(bo^2)*beta*V32 + Chi*Cac*V21 - 2*bo*delta*V22 - gam_3*V22;
    dV30 = 3*alpha*Cac*V20 - 3*(bo^2)*beta*V30 - 2*alpha*Cac*V30 + 4*(bo^3)*beta*V40 - 2*Chi*Cac*V30 + delta*V31;
    dV31 = 3*alpha*Cac*V21 - 3*(bo^2)*beta*V31 - 2*alpha*Cac*V31 + 4*(bo^3)*beta*V41 + 2*Chi*Cac*V30 - delta*V31 - Chi*Cac*V31 + 2*bo*delta*V32;
    dV32 = 3*alpha*Cac*V22 - 3*(bo^2)*beta*V32 - 2*alpha*Cac*V32 + 4*(bo^2)*beta*V42 + Chi*Cac*V31 - 2*bo*delta*V32 - gam_3*V32;
    dV40 = 2*alpha*Cac*V30 - 4*(bo^3)*beta*V40 - alpha*Cac*V40 + 5*(bo^4)*beta*V50 - 2*Chi*Cac*V40 + delta*V41;
    dV41 = 2*alpha*Cac*V31 - 4*(bo^3)*beta*V41 - alpha*Cac*V41 + 5*(bo^4)*beta*V51 + 2*Chi*Cac*V40 - delta*V41 - Chi*Cac*V41 + 2*bo*delta*V42;
    dV42 = 2*alpha*Cac*V32 - 4*(bo^3)*beta*V42 - alpha*Cac*V42 + 5*(bo^4)*beta*V52 + Chi*Cac*V41 - 2*bo*delta*V42 - gam_3*V42;
    dV50 = alpha*Cac*V40 - 5*(bo^4)*beta*V50 - 2*Chi*Cac*V50 + delta*V51 - gam_2*V50;
    dV51 = alpha*Cac*V41 - 5*(bo^4)*beta*V51 + 2*Chi*Cac*V50 - delta*V51 - Chi*Cac*V51 + 2*bo*delta*V52 - gam_2*V51;
    dV52 = alpha*Cac*V42 - 5*(bo^4)*beta*V52 + Chi*Cac*V51 - 2*bo*delta*V52 - gam_3*V52 - gam_2*V52;
    

    dW00 = Kattach*V00*Can*(1 - RFw) - Kdetach*W00 + ...
             beta*W10 - 5*alpha*W00*Can + delta*W01 - 2*Chi*W00*Can - gam_1*W00;
    dW01 = 2*Chi*W00*Can - delta*W01 + 2*bo*delta*W02 - Chi*Can*W01 + beta*W11 - 5*alpha*Can*W01;
    dW02 = Chi*Can*W01 - 2*bo*delta*W02 - gam_3*W02 - 5*alpha*Can*W02 + beta*W12;
    dW10 = 5*alpha*Can*W00 - beta*W10 - 4*alpha*Can*W10 + 2*bo*beta*W20 - 2*Chi*Can*W10 + delta*W11;
    dW11 = 5*alpha*Can*W01 - beta*W11 - 4*alpha*Can*W11 + 2*bo*beta*W21 + 2*Chi*Can*W10 - delta*W11 - Chi*Can*W11 + 2*bo*delta*W12;
    dW12 = 5*alpha*Can*W02 - beta*W12 - 4*alpha*Can*W12 + 2*bo*beta*W22 + Chi*Can*W11 - 2*bo*delta*W12 - gam_3*W12;
    dW20 = 4*alpha*Can*W10 - 2*bo*beta*W20 - 3*alpha*Can*W20 + 3*(bo^2)*beta*W30 - 2*Chi*Can*W20 + delta*W21;
    dW21 = 4*alpha*Can*W11 - 2*bo*beta*W21 - 3*alpha*Can*W21 + 3*(bo^2)*beta*W31 + 2*Chi*Can*W20 - delta*W21 - Chi*Can*W21 + 2*bo*delta*W22;
    dW22 = 4*alpha*Can*W12 - 2*bo*beta*W22 - 3*alpha*Can*W22 + 3*(bo^2)*beta*W32 + Chi*Can*W21 - 2*bo*delta*W22 - gam_3*W22;
    dW30 = 3*alpha*Can*W20 - 3*(bo^2)*beta*W30 - 2*alpha*Can*W30 + 4*(bo^3)*beta*W40 - 2*Chi*Can*W30 + delta*W31;
    dW31 = 3*alpha*Can*W21 - 3*(bo^2)*beta*W31 - 2*alpha*Can*W31 + 4*(bo^3)*beta*W41 + 2*Chi*Can*W30 - delta*W31 - Chi*Can*W31 + 2*bo*delta*W32;
    dW32 = 3*alpha*Can*W22 - 3*(bo^2)*beta*W32 - 2*alpha*Can*W32 + 4*(bo^2)*beta*W42 + Chi*Can*W31 - 2*bo*delta*W32 - gam_3*W32;
    dW40 = 2*alpha*Can*W30 - 4*(bo^3)*beta*W40 - alpha*Can*W40 + 5*(bo^4)*beta*W50 - 2*Chi*Can*W40 + delta*W41;
    dW41 = 2*alpha*Can*W31 - 4*(bo^3)*beta*W41 - alpha*Can*W41 + 5*(bo^4)*beta*W51 + 2*Chi*Can*W40 - delta*W41 - Chi*Can*W41 + 2*bo*delta*W42;
    dW42 = 2*alpha*Can*W32 - 4*(bo^3)*beta*W42 - alpha*Can*W42 + 5*(bo^4)*beta*W52 + Chi*Can*W41 - 2*bo*delta*W42 - gam_3*W42;
    dW50 = alpha*Can*W40 - 5*(bo^4)*beta*W50 - 2*Chi*Can*W50 + delta*W51 - gam_2*W50;
    dW51 = alpha*Can*W41 - 5*(bo^4)*beta*W51 + 2*Chi*Can*W50 - delta*W51 - Chi*Can*W51 + 2*bo*delta*W52 - gam_2*W51;
    dW52 = alpha*Can*W42 - 5*(bo^4)*beta*W52 + Chi*Can*W51 - 2*bo*delta*W52 - gam_3*W52 - gam_2*W52;

    dydt_5 = [dV00 dV01 dV02 dV10 dV11 dV12 dV20 dV21 dV22 dV30 dV31 dV32...
            dV40 dV41 dV42 dV50 dV51 dV52 dW00 dW01 dW02 dW10 dW11 dW12 dW20...
            dW21 dW22 dW30 dW31 dW32 dW40 dW41 dW42 dW50 dW51 dW52 dRdt dUdt...
            dRFvdt dRFwdt];
        
end


function dydt_6 = KR_FF(t, y, params, Cac, Can)

    Kf4 = params(1);
    Kb4 = params(2);
    b4 = params(3);
    Kspont4 = params(4); 
    Kevoked4 = params(5);
    Kef = params(6);
    Kreacid4 = params(7);
    Kf7 = params(8);
    Kb7 = params(9);
    b7 = params(10);
    Kspont7 = params(11);
    Kevoked7 = params(12);
    Kes = params(13);
    Kreacid7 = params(14);
    Kmob = params(15);
    Kdoc = params(16);

    F0 = y(1); F1 = y(2); F2 = y(3); Evoked4 = y(4); Spont4 = y(5); Endo4 = y(6); Reacid4 = y(7);
    S0 = y(8); S1 = y(9); S2 = y(10); S3 = y(11); S4 = y(12); S5 = y(13); Evoked7 = y(14); Spont7 = y(15); Endo7 = y(16); Reacid7 = y(17);

    dF0dt = Kb4*F1 + Kdoc*S0 - 2*Kf4*Can*F0 - Kspont4*F0;
    dF1dt = 2*Kf4*Can*F0 - Kb4*F1 + 2*b4*Kb4*F2 - Kf4*Can*F1;
    dF2dt = Kf4*Can*F1 - 2*b4*Kb4*F2 - Kevoked4*F2;
    dEvoked4dt = Kevoked4*F2 - Kef*Evoked4;
    dSpont4dt = Kspont4*F0 - Kef*Spont4;
    dEndo4dt = Kef*(Evoked4 + Spont4) - Kreacid4*Endo4;
    dReacid4dt = Kreacid4*Endo4 - Kmob*Reacid4;


    dS0dt = Kb7*S1 - 5*Kf7*Can*S0 + Kmob*(Reacid4 + Reacid7) - Kdoc*S0 - Kspont7*S0;
    dS1dt = 5*Kf7*Can*S0 - Kb7*S1 + 2*b7*Kb7*S2 - 4*Kf7*Can*S1;
    dS2dt = 4*Kf7*Can*S1 - 2*b7*Kb7*S2 + 3*Kb7*(b7^2)*S3 - 3*Kf7*Can*S2;
    dS3dt = 3*Kf7*Can*S2 - 3*Kb7*(b7^2)*S3 + 4*Kb7*(b7^3)*S4 - 2*Kf7*Can*S3;
    dS4dt = 2*Kf7*Can*S3 - 4*Kb7*(b7^3)*S4 + 5*Kb7*(b7^4)*S5 - Kf7*Can*S4;
    dS5dt = Kf7*Can*S4 - 5*Kb7*(b7^4)*S5 - Kevoked7*S5;
    dEvoked7dt = Kevoked7*S5 - Kes*Evoked7;
    dSpont7dt = Kspont7*S0 - Kes*Spont7;
    dEndo7dt = Kes*(Evoked7 + Spont7) - Kreacid7*Endo7;
    dReacid7dt = Kreacid7*Endo7 - Kmob*Reacid7;
    
    dydt_6 = [dF0dt dF1dt dF2dt dEvoked4dt dSpont4dt dEndo4dt dReacid4dt... 
              dS0dt dS1dt dS2dt dS3dt dS4dt dS5dt dEvoked7dt dSpont7dt dEndo7dt dReacid7dt];
end