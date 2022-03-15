clear all;


%%  Fit data


%%  Load Data


load CaM.txt
load CaN.txt
load DualSensor_RelRate.txt
load Allosteric_RelRate.txt
load DualSensor_Dose_Resp.txt
load Allosteric_Dose_Resp.txt


Cac = sort(CaM(:));
Can = sort(CaN(:));
ydata = DualSensor_Dose_Resp(:);
ydata2 = Allosteric_Dose_Resp(:);
DualSensorModel_Fast = sort(DualSensor_RelRate(:));
AllostericModel_Fast = sort(Allosteric_RelRate(:));



%% Fitting to data



%% Modified Dual Sensor

%{
a0 = [0.061200 2.320000 0.003814 0.014682 0.025008 0.250006 2.000008 0.000002...
      6.340000 0.000050 0.0022 0.027991 0.005208 0.000150 0.001093];                  % Modified Dual Sensor Model
  
lb = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
ub = [];

[a,resnorm,residual] = lsqcurvefit(@(a,Can)myode_for(a,Can,Cac,"DualSensor"),a0,Can,ydata,lb,ub);

sprintf('%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f',a(1),a(2),a(3),a(4),...
        a(5),a(6),a(7),a(8),a(9),a(10),a(11),a(12),a(13),a(14),a(15))
%sprintf('%f %f', resnorm, residual)
%}


%% Modified Allosteric Model

% Fitted Parameters = [0.097909 9.316730 0.000000 28.693830 0.558504 6.339942...
                      %0.003858 0.002192 0.028560 0.003124 0.000144 2.413995];
a0 = [0.102953 3.999998 0.0000002 31.300000 0.50000 6.34000 0.003858...
      0.002193 0.028560 0.003186 0.000150 0.000509];   % Modified Allosteric Model
lb = [0 0 0.0000002 0 0 0 0 0 0 0 0 0];
ub = [];

[a,resnorm,residual]= lsqcurvefit(@(a,Can)myode_for(a,Can,Cac,"Allosteric"),a0,Can,ydata,lb,ub);

sprintf('%f %f %f %f %f %f %f %f %f %f %f %f',a(1),a(2),a(3),a(4),a(5),a(6),...
        a(7),a(8),a(9),a(10),a(11),a(12))
%sprintf('%f %f', resnorm, residual)




%% Estimate release rate

%ModifiedDualSensorModel = myode_for(a, Can, Cac);
ModifiedAllostericModel = myode_for(a, Can, Cac, "Allosteric");



%% Plot Figures
%{
figure
loglog(Can, ydata, "k.", 'LineWidth', 1, 'MarkerSize', 5)
%hold on
loglog(Can, ModifiedAllostericModel, 'LineWidth', 2, 'MarkerSize', 10)
hold on
legend({'Log Transoformed Allosteric Model'},'Location', 'northwest')
title("Dose Response Curve From Allosteric Model")
box off

%Get Current Figure (GCF) & Set image size before saving image
width = 20;  % cm 
height = 19; % cm
set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])

%Set the resolution of 1000dpi and save the plot in TIFF format 
print -dpng -r1000 Allosteric_Dose_Resp_Fit
hold off
%}


%% FUNCTIONS


function y_real = myode_for(params, Can, Cac, text)

    dt = 1e-03;  
    tspan = 0:dt:1;
    
    yMA = zeros(length(Cac),1); 
    yMA1 = zeros(length(tspan),1);
    yMA2 = zeros(length(tspan),1); 
    
    yMD = zeros(length(Cac),1); 
    yMD1 = zeros(length(tspan),1);
    yMD2 = zeros(length(tspan),1); 

    y_real = zeros(length(Cac),1); 

    
    for i=1:length(Cac)
        
        
        %% Modified Allosteric Model
        
        if text=="Allosteric"

            I = params(3);           % /ms
            F = params(4);

            y0MD = [3 0 0 0 0 0 3 0 0 0 0 0 85 10 0 0];

            [t1, y1] = ode15s(@(t1,y1)ModifiedAllosteric(t1, y1, params, Cac(i), Can(i)),tspan,y0MD);

            %{
            V = zeros(length(t1),1);
            for iv = 1:14
                V = V + y1(:,iv);
            end
            %}

            yMA1(:,1) = I.*(y1(:,1) + F.*y1(:,2) + (F^2).*y1(:,3)* + (F^3).*y1(:,4) + (F^4).*y1(:,5) + (F^5).*y1(:,6));          

            yMA2(:,1) = I.*(y1(:,7) + F.*y1(:,8) + (F^2).*y1(:,9)* + (F^3).*y1(:,10) + (F^4).*y1(:,11) + (F^5).*y1(:,12));          

            y_real(i,1) = max(yMA2(:,1) + yMA1(:,1));
        
        end 
        
        
        
       %% Modified Dual Sensor Model 
        
        if text=="DualSensor"

            a = params(5);
            gam_2 = params(7);     % /ms Synchronous release rate
            gam_3 = a*gam_2;       % /ms Asynchronous release rate
            gam_1 = params(8);     % /ms Spontaneous release rate


            y0MD = [3 0 0 0 0 0 0 0 0 0 0 0 ...
                    0 0 0 0 0 0 3 0 0 0 0 0 0 ...
                    0 0 0 0 0 0 0 0 0 0 0 85 10 ...
                    0 0];

            [t5, y5] = ode15s(@(t5,y5)ModifiedDualSensor(t5, y5, params, Cac(i), Can(i)),tspan,y0MD);

            %{
            V = zeros(length(t5),1);
            for iv = 1:36 
                V = V + y5(:,iv);
            end
            %}

            yMD1(:,1) = gam_1.*y5(:,1) + gam_2.*(y5(:,16) + y5(:,17) + y5(:,18))+...
                     gam_3.*(y5(:,3) + y5(:,6) + y5(:,9) + y5(:,12) + y5(:,15) + y5(:,18));          

            yMD2(:,1) = gam_1.*y5(:,19) + gam_2.*(y5(:,34) + y5(:,35) + y5(:,35))+...
                     gam_3.*(y5(:,21) + y5(:,24) + y5(:,27) + y5(:,30) + y5(:,33) + y5(:,36));          

            y_real(i,1) = max(yMD2(:,1) + yMD1(:,1));
        
        end
        
              
    end

end 



function dydt = ModifiedAllosteric(t, y, params, Cac, Can)

%{
Allosteric modulation of the presynaptic Ca21 sensor for vesicle fusion
Xuelin Lou1, Volker Scheuss1? & Ralf Schneggenburger1
Calyx Of Held
%}

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
       
    
    dydt=[dV0; dV1dt; dV2dt; dV3dt; dV4dt; dV5dt; dW0; dW1dt; dW2dt; dW3dt; dW4dt; dW5dt;...
           dRdt; dUdt; dRFvdt; dRFwdt];

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

    dydt_5 = [dV00; dV01; dV02; dV10; dV11; dV12; dV20; dV21; dV22; dV30; dV31; dV32;...
            dV40; dV41; dV42; dV50; dV51; dV52; dW00; dW01; dW02; dW10; dW11; dW12; dW20;...
            dW21; dW22; dW30; dW31; dW32; dW40; dW41; dW42; dW50; dW51; dW52; dRdt; dUdt;...
            dRFvdt; dRFwdt];
end

