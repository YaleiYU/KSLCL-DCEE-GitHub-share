

%% Plot Function 
close all
clear
clc

%%
set(groot, 'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(0, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [0 0 16 9]);


filenames={
           'out_lksdcee_80400_111_t30s_n05.mat', 'out_dcee_80400_111_t30s_n05.mat', ...
           'out_po_80400_111_t30s_n05.mat'
          };
 

 for kk=1:numel(filenames)
     load(filenames{kk})
 end


%% 
power_lksdcee_80400_111_t30s = out_lksdcee_80400_111_t30s_n05.power_lksdcee.data;
voltage_lksdcee_80400_111_t30s = out_lksdcee_80400_111_t30s_n05.voltage_lksdcee.data; 
current_lksdcee_80400_111_t30s = out_lksdcee_80400_111_t30s_n05.current_lksdcee.data; 
theta1_lksdcee_80400_111_t30s = out_lksdcee_80400_111_t30s_n05.theta1_lksdcee.data; 
theta2_lksdcee_80400_111_t30s = out_lksdcee_80400_111_t30s_n05.theta2_lksdcee.data; 
theta3_lksdcee_80400_111_t30s = out_lksdcee_80400_111_t30s_n05.theta3_lksdcee.data;
Vbar_lksdcee_80400_111_t30s = out_lksdcee_80400_111_t30s_n05.Vbar_lksdcee.data;


power_po_80400_111_t30s_n05 = out_po_80400_111_t30s_n05.power_po.data;
voltage_po_80400_111_t30s_n05 = out_po_80400_111_t30s_n05.voltage_po.data; 
current_po_80400_111_t30s_n05 = out_po_80400_111_t30s_n05.current_po.data; 


%%  
T3 = 30;

Ts = 1e-4;
Ts2 = 1e-6;
N = int64(T/Ts + 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

power_max_80400_111_t30s = zeros(length(t3),1);
voltage_max_80400_111_t30s = zeros(length(t3),1);
current_max_80400_111_t30s = zeros(length(t3),1);


power_max_80400_111_t30s(t3>=0 & t3<5) = 139.994; 
power_max_80400_111_t30s(t3>=5 & t3<=10) = 139.994; 
power_max_80400_111_t30s(t3>=10 & t3<20) = 172.663;
power_max_80400_111_t30s(t3>=20 & t3<=30) = 204.402;


voltage_max_80400_111_t30s(t3>=0 & t3<5) =  27.7733; 
voltage_max_80400_111_t30s(t3>=5 & t3<=10) =  27.7733;  
voltage_max_80400_111_t30s(t3>=10 & t3<20) = 34.309;
voltage_max_80400_111_t30s(t3>=20 & t3<=30) = 40.8195; 


current_max_80400_111_t30s(t3>=0 & t3<5) = 5.04061; 
current_max_80400_111_t30s(t3>=5 & t3<=10) = 5.04061;  
current_max_80400_111_t30s(t3>=10 & t3<20) = 5.03257;  
current_max_80400_111_t30s(t3>=20 & t3<=30) = 5.00745;  


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1) 
plot(0:Ts:T3, power_lksdcee_80400_111_t30s, '-. k', ...
     0:Ts:T3, power_dcee_80400_111_t30s, '-. b', ...
     0:Ts2:T3, power_po_80400_111_t30s_n05, '-- r', ... 
     0:Ts:T3, power_max_80400_111_t30s , ': g', ...
    'LineWidth', 1.5)
xlim([0,T3]);
my_legend = legend(' $ P $, LKS-DCEE, $ k=3, N=20, \tau_{d}=400 $', ...
                   ' $ P $, DCEE, $ N=20, \tau_{d}=400 $ ', ...
                   ' $ P $, P \& O',  ...
                   'Optimal Reference');
set(my_legend, 'FontSize', 9, 'Location', 'south');
xlabel('$ t[s]  $','FontSize', 9);
ylabel('Power  $ P [W]$ ', 'FontSize', 9);


figure(2)
plot(0:Ts:T3, voltage_lksdcee_80400_111_t30s, '-. k', ...
     0:Ts:T3, voltage_dcee_80400_111_t30s, '-.b', ...
     0:Ts2:T3, voltage_po_80400_111_t30s_n05, '--r', ...
     0:Ts:T3, voltage_max_80400_111_t30s, ': g', ...
    'LineWidth',1.5)
xlim([0,T3]);
my_legend = legend(' $ V $, LKS-DCEE, $ k=3, N=20, \tau_{d}=400 $', ...
                   ' $ V $, DCEE, $ N=20, \tau_{d}=400 $', ...
                   ' $ V $, P \& O',  ...
                   'Optimal Reference');
set(my_legend,'FontSize', 9, 'Location','south');
xlabel('$t[s]$','FontSize', 9);
ylabel('Voltage  $ V [V] $ ', 'FontSize', 9);


figure(3)
plot(0:Ts:T3, current_lksdcee_80400_111_t30s, '-. k', ... 
     0:Ts:T3, current_dcee_80400_111_t30s, '-.b', ...
     0:Ts2:T3, current_po_80400_111_t30s_n05, '--r', ...
     0:Ts:T3, current_max_80400_111_t30s, ':g', ...
    'LineWidth', 1.5)
xlim([0,T3]);
my_legend = legend(' $ I $, LKS-DCEE, $ k=3, N=20, \tau_{d}=400 $', ...
                   ' $ I $, DCEE, $ N=20, \tau_{d}=400 $', ...
                   ' $ I $, P \& O',  ...
                   'Optimal Reference');
set(my_legend,'FontSize',9, 'Location','south');
xlabel('$t[s]$', 'FontSize', 9);
ylabel('Current $ I [A] $ ', 'FontSize', 9);

