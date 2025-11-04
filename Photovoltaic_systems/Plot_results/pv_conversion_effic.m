
clc;
clear;
close all; 


%%
set(groot, 'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(0, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [0 0 16 9]);

%%
filenames={      
           'out_lksdcee_80400_111_t30s_n05.mat', 'out_dcee_80400_111_t30s_n05.mat', ...
           'out_po_80400_111_t30s_n05.mat'
          };
 
 for kk=1:numel(filenames)
     load(filenames{kk})
 end

%% 
Ts = 1e-4;

t30 = 0:Ts:30;


%% calculate the conversion efficiency
panel_length = 1.576;
panel_width = 0.825; 
panel_area = panel_length*panel_width;

V_max = 36.6;
I_max = 4.78;

power = V_max*I_max; 
power_area = power/panel_area; 

power_in = 1000; 

conversion_effic = power_area/power_in*100; 


%% conversion efficiency: ratio between output and input energy  
power_lksdcee_80400_111_t30s_n05 = out_lksdcee_80400_111_t30s_n05.power_lksdcee.data;
power_dcee_80400_111_t30s_n05 = out_dcee_80400_111_t30s_n05.power_dcee.data;
power_po_80400_111_t30s_n05 = out_po_80400_111_t30s_n05.power_po.data;


%% power loss 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
power_max_80400_111_t30s_n05 = zeros(size(t30));

power_max_80400_111_t30s_n05(t30>=0 & t30<5) =  139.994; 
power_max_80400_111_t30s_n05(t30>=5 & t30<10) =  139.994; 
power_max_80400_111_t30s_n05(t30>=10 & t30<20) = 172.663; 
power_max_80400_111_t30s_n05(t30>=20 & t30<30) = 204.402; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
powererror_lksdcee_80400_111_t30s_n05 = abs(power_max_80400_111_t30s_n05' - power_lksdcee_80400_111_t30s_n05); 
powererror_dcee_80400_111_t30s_n05 = abs(power_max_80400_111_t30s_n05' - power_dcee_80400_111_t30s_n05);
powererror_po_80400_111_t30s_n05 = abs(power_max_80400_111_t30s_n05 - power_po_80400_111_t30s_n05);

conversion_effic_lksdcee = sum(power_lksdcee_80400_111_t30s_n05*Ts)/sum(power_max_80400_111_t30s_n05'*Ts)*100;
conversion_effic_dcee = sum(power_dcee_80400_111_t30s_n05*Ts)/sum(power_max_80400_111_t30s_n05'*Ts)*100;
conversion_effic_po = sum(power_po_80400_111_t30s_n05*Ts)/sum(power_max_80400_111_t30s_n05'*Ts)*100;


%%
figure(1)
plot(t30, cumtrapz(powererror_lksdcee_80400_111_t30s_n05).*Ts, '-g', 'LineWidth', 1.5);
hold on;
plot(t30, cumtrapz(powererror_dcee_80400_111_t30s_n05)*Ts, '-.b', 'LineWidth', 1.5);
hold on;
plot(t30, cumtrapz(powererror_po_80400_111_t30s_n05)*Ts, '--r', 'LineWidth', 1.5)
hold off;
xlim([0,30]);
my_legend = legend(' LKS-DCEE', ' DCEE', ' P\&O ' );
set(my_legend,'FontSize', 9, 'Location', 'north');
xlabel('$ t [s] $','FontSize', 9);
ylabel(' Energy Loss $ [J] $ ','FontSize', 9);



