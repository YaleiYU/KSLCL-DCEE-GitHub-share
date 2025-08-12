    clc
    clear
    close all    

    %% Simulation Parameters
    thetaStar = 1;
    
    thetamin = 0.01;
    thetamax = 5;
    
    
    %% Initialise Estimators
    N = 100; 
    theta1 = thetamin + (thetamax - thetamin)*rand(N,1);
    theta = theta1;
    
    % system dynamics and gains
    A = [0 1; 2 1];
    B = [1; 1];
    C = [0 1];
    
    Q = [C zeros(1,1); A-eye(2) B];
    rank(Q)
    X = (Q'*Q)^(-1)*Q'*[1; 0 ; 0];
    Psi = X(1:2);
    G = X(3);
    
    xi_k = [0];
    
    P1= [0.4; 0.7];
    
    K = place(A, B, P1);
    
    xmin = -100;
    xmax = 100;
    
    x0 = [2 1];
    x_k = x0;
    
    y0 = C*x0';
    y_k = y0;

    xi_k_in = 0.1; 

    
    % Learning rates
    eta0 = 0.00018;
    d1 = 0.01;
    delta0 = 0.42;
    d2 = 0.4;
    alpha = 0.5; 
    
    theta_store_lksdcee_n8 = [];
    y_k_store_lksdcee_n8 =[];
    xi_k_store_lksdcee_n8 =[];
    x_k_store_lksdcee_n8 =[];
    Reward_store_lksdcee_n8 = [];
    Jbelief_store_lksdcee_n8 = []; 
    uk_store_lksdcee_n8 = []; 
    ref_store_lksdcee_n8 = []; 
    
    thetaMean_store = [];
    eta_store = []; 
    dual_store = []; 
    v_kN_store = []; 

    para_conv_theta1_taud = [];
    para_conv_theta2_taud = [];

    para_conv_theta1_store = [];
    para_conv_theta2_store = [];

    para_conv_theta1_in_taud = [];
    para_conv_theta2_in_taud = [];

    para_conv_theta1_in_store = [];
    para_conv_theta2_in_store = [];


    %% 
    T = 80;
    n = 8;
    
    %% Iteration
    for i = 1: T   
        eta = eta0; 
        delta = delta0;
        
        %%
        RewardSim = reward(thetaStar, y_k);
        noise = 0.5*randn(size(RewardSim));
        RewardSim = RewardSim + noise;
    
        Reward_store_lksdcee_n8 = [Reward_store_lksdcee_n8; RewardSim];
    
        % update estimators
        Jbelief = reward(theta, y_k);
        Jbelief_store_lksdcee_n8 = [Jbelief_store_lksdcee_n8; Jbelief']; 
    
        x_k_store_lksdcee_n8 = [x_k_store_lksdcee_n8; x_k];
        y_k_store_lksdcee_n8 = [y_k_store_lksdcee_n8; y_k];
    
        %% concurrent learning 
        tau_d = 80;
        
        theta1_sum_cl_store = [];
        theta2_sum_cl_store = []; 

        if  i < tau_d 


            eta = eta*(tau_d)/i*0.1; 
 

            for j = 1:i
                theta1_sum_cl = (reward(theta, y_k_store_lksdcee_n8(j)) - Reward_store_lksdcee_n8(j)*ones(N,1)).*(-y_k_store_lksdcee_n8(j).^2)./(1+y_k_store_lksdcee_n8(j).^4.*eta); 
                 
                theta1_sum_cl_store = [theta1_sum_cl_store; theta1_sum_cl'];
           

                para_conv_theta1 =  y_k_store_lksdcee_n8(j).^2.*y_k_store_lksdcee_n8(j).^2; 
                

                para_conv_theta1_taud = [para_conv_theta1_taud para_conv_theta1];
             
            end
            
        else
            for j = (i-tau_d+1):i
                theta1_sum_cl = (reward(theta, y_k_store_lksdcee_n8(j)) - Reward_store_lksdcee_n8(j)*ones(N,1)).*(-y_k_store_lksdcee_n8(j).^2)./(1+y_k_store_lksdcee_n8(j).^4.*eta);
                 
                theta1_sum_cl_store = [theta1_sum_cl_store; theta1_sum_cl'];
                 

                para_conv_theta1 =  y_k_store_lksdcee_n8(j).^2.*y_k_store_lksdcee_n8(j).^2; 
               

                para_conv_theta1_taud = [para_conv_theta1_taud para_conv_theta1];
                
            end
        end
    
        if  tau_d == 0
%            
            theta1 = theta1 - eta.*(Jbelief - RewardSim*ones(N,1)).*(-y_k_store_lksdcee_n8(j).^2)./(1+y_k.^4.*eta);
            
        else

            theta1_sum_cl = sum(theta1_sum_cl_store', 2);
           
            
            theta1 = theta1 - eta.*theta1_sum_cl; 
          

            para_conv_theta1_taud_sum = abs(1 - eta.*sum(para_conv_theta1_taud));
            
        end

        %%
        theta = [theta1 ];

        theta_store_lksdcee_n8 = [theta_store_lksdcee_n8 theta];
    
        theta_size = size(theta);
    

        %% the innner loop updates in the lkscl-dcee algorithm 
        Jbelief_in_store = [];
        Reward_in_store = [];
        y_k_in_store = [];
        dual_store = []; 

        theta1_in = theta1; 
        theta_in = [theta1_in ];

        y_k_in = y_k; 
        x_k_in = x_k; 
        

        tic;


        for j = 1 : n

            RewardSim_in = reward(thetaStar, y_k_in);
            noise = 0*randn(size(RewardSim_in));
            RewardSim_in = RewardSim_in + noise;

            Reward_in_store = [Reward_in_store; RewardSim_in];

            Jbelief = reward(theta, y_k_in);
            Jbelief_in_store = [Jbelief_in_store; Jbelief'];

            y_k_in_store = [y_k_in_store; y_k_in];

            %%
            i = j;

        theta1_sum_cl_in_store = [];
        
        if  i < tau_d 

            eta = eta*(tau_d)/i*0.1; 

            for j = 1:i
                theta1_sum_cl_in = (reward(theta, y_k_in_store(j)) - Reward_in_store(j)*ones(N,1)).*(-y_k_in_store(j).^2)./(1+y_k_in_store(j).^4.*eta); 
                
                theta1_sum_cl_in_store = [theta1_sum_cl_in_store; theta1_sum_cl_in'];
          
                para_conv_theta1 =  y_k_in_store(j).^2.*y_k_in_store(j).^2; 

                para_conv_theta1_in_taud = [para_conv_theta1_in_taud; para_conv_theta1];
            end
            
        else
            for j = (i-tau_d+1):i
                theta1_sum_cl_in = (reward(theta, y_k_in_store(j)) - Reward_in_store(j)*ones(N,1)).*(-y_k_in_store(j).^2)./(1+y_k_in_store(j).^4.*eta);
                
                theta1_sum_cl_in_store = [theta1_sum_cl_in_store; theta1_sum_cl_in'];
                
                para_conv_theta1 =  y_k_in_store(j).^2.*y_k_in_store(j).^2; 
            
                para_conv_theta1_in_taud = [para_conv_theta1_in_taud; para_conv_theta1];
            end
        end
    
        if  tau_d == 0
%      
            theta1_in = theta1_in - eta.*(Jbelief - RewardSim*ones(N,1)).*(-y_k_in_store(j).^2)./(1+y_k_in.^4.*eta);
            
        else

            theta1_sum_cl_in = sum(theta1_sum_cl_in_store', 2);
 
            
            theta1_in = theta1_in - eta.*theta1_sum_cl_in;  

            para_conv_theta1_in_taud_sum = abs(1 - eta.*sum(para_conv_theta1_in_taud));
        end

            %%
            theta_in = [theta1_in ];

            theta_in(theta_in < thetamin) = thetamin;
            theta_in(theta_in > thetamax) = thetamax;

            [future_mean_in, future_variance_in] = predicted_mean_variance(theta_in, y_k_in, eta);

            ref_in = 2./mean(theta1_in)./2;
 

            exploitation_gradient = (y_k_in - ref_in);

            exploration_gradient = (get_variance_gradient(theta_in, y_k_in, eta));

            dual_in = - delta*exploitation_gradient - delta*exploration_gradient;

            xi_k_in = xi_k_in + dual_in;

            x_k_in = ((A-B*K)*x_k_in' + B*(G+K*Psi)*xi_k_in)';

            x_k_in(x_k_in < xmin) = xmin;
            x_k_in(x_k_in > xmax) = xmax;

            y_k_in = C*x_k_in';

            dual_store = [dual_store xi_k_in];

        end

        elapsedTime_lksdcee_n8 = toc;


        v_kN = dual_store(n);

        xi_k = xi_k + alpha*(v_kN - xi_k);

        %%  decision making 
        ref = 2./mean(theta(:,1))./2; 

        u_k = -K*x_k' + B*(G+K*Psi)*xi_k; 
    
        x_k = ((A-B*K)*x_k' + B*(G+K*Psi)*xi_k)';
        x_k(x_k<xmin) = xmin;
        x_k(x_k>xmax) = xmax;
    
        y_k = C*x_k';

        %% transfer information into the inner loop
        xi_k_in = xi_k; 
        x_k_in = x_k; 

        %% store datasets
    
        xi_k_store_lksdcee_n8 = [xi_k_store_lksdcee_n8; xi_k];
        ref_store_lksdcee_n8 = [ref_store_lksdcee_n8; ref]; 
        eta_store = [eta_store; eta]; 
        v_kN_store = [v_kN_store; v_kN]; 

        uk_store_lksdcee_n8 = [uk_store_lksdcee_n8; u_k']; 

        %%    
        para_conv_theta1_store = [para_conv_theta1_store; para_conv_theta1_taud_sum];

    
    end
   
    %% Plot Function
    
    set(groot, 'defaulttextinterpreter','latex');
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');
    set(0, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [0 0 16 9]);
    
    
    %% 
    thetaMean=mean(theta_store_lksdcee_n8);
    thetaVariance=var(theta_store_lksdcee_n8);
    std_dev = sqrt(thetaVariance);
     
    Jbelief = -thetaMean(:,1).*y_k_store_lksdcee_n8.^2 + 2.*y_k_store_lksdcee_n8; 

    %% 
    figure(1)
    t=1:1:T;
    plot(t, thetaMean', 'r-', t, thetaStar(1)*ones(T,1), 'g--', 'LineWidth',1.5);
    ylabel('$\theta_1$','FontSize',12);
    my_legend = legend('$\hat{\theta}_1$','optimal reference');
    set(my_legend,'FontSize',12, 'Location', 'northeast');
    


    figure(11)
    patch([t fliplr(t)], [thetaMean-std_dev  fliplr(thetaMean+std_dev)], [0.3010, 0.7450, 0.9330])
    hold on
    plot(1:1:T, thetaMean,'LineWidth',1.5, 'Color', '[0.4940, 0.1840, 0.5560]');
    plot(1:1:T, ones(T,1), 'Color', 'r', 'LineWidth',1.5,'LineStyle', '--')
    box on
    xlim([0,T]);
    my_legend = legend('standard deviation of $\theta$','estimated mean of $\theta$', 'optimal $\theta^*$');
    set(my_legend,'FontSize',12);
    xlabel('$k$','FontSize',12);
    ylabel('$\theta$','FontSize',12);


    figure(2)
    plot(1:1:T, y_k_store_lksdcee_n8(:,1),'Color', 'k', 'LineWidth',1.5);
    hold on;
    plot(1:1:T, 1*ones(T,1),'Color', 'r', 'LineWidth',1.5,'LineStyle', '--')
    hold on;
    plot(1:1:T, ref_store_lksdcee_n8(:,1), 'g :', 'LineWidth',1.5);
    hold on;
    plot(1:1:T, xi_k_store_lksdcee_n8(:,1), 'c-.', 'LineWidth',1.5);
    hold on;
    plot(1:1:T, v_kN_store(:,1), 'b-.', 'LineWidth',1.5);
    xlim([0,T]);
    my_legend = legend('output $y$','optimal reference', 'generated reference', 'trajectory $\xi$', '$ v_{kn} $');
    set(my_legend,'FontSize',12,'Location', 'southeast');
    xlabel('$ n $','FontSize',12);
    ylabel('$y$','FontSize',12);
    

    figure(3)
    plot(1:1:T, x_k_store_lksdcee_n8(:,1:2), 'LineWidth', 1.5);
    hold on;
    % plot(1:1:T+1,ones(T+1,1),'Color', 'r', 'LineWidth',1.5,'LineStyle', '--')
    xlim([0,T]);
    my_legend = legend('state $x_1$', 'state $x_2$');
    set(my_legend,'FontSize',12);
    xlabel('$ n $','FontSize',12);
    ylabel('$ x $','FontSize',12);
    

    figure(4)
    plot(1:1:T, Reward_store_lksdcee_n8, 'Color', '[0, 0.4470, 0.7410]', 'LineWidth',1.5);
    hold on;
    plot(1:1:T, 1*ones(T,1),'Color', 'r', 'LineWidth',1.5,'LineStyle', '--')
    xlim([0,T]);
    % my_legend = legend('observed reward $J$', 'estimated reward $\hat{J}$', 'optimal reward','Location','southeast');
    my_legend = legend('observed reward $J$', 'optimal reward','Location','southeast');
    set(my_legend,'FontSize',12);
    xlabel('$ n $','FontSize',12);
    ylabel('$ J $','FontSize',12);


    figure(5)
    plot(1:1:T, uk_store_lksdcee_n8(:,1), '--g', 'LineWidth', 1.5);
    hold on;
    plot(1:1:T, uk_store_lksdcee_n8(:,2), ':r', 'LineWidth', 1.5);
    % plot(1:1:T+1,ones(T+1,1),'Color', 'r', 'LineWidth',1.5,'LineStyle', '--')
    xlim([0,T]);
    my_legend = legend('state $u_1$', 'state $u_2$');
    set(my_legend,'FontSize',12);
    xlabel('$ n $','FontSize',12);
    ylabel('$ u $','FontSize',12);


    figure(6)
    plot(1:1:T, y_k_store_lksdcee_n8(:,1),'Color', 'k', 'LineWidth',1.5);
    hold on;
    plot(1:1:T, 1*ones(T,1),'Color', 'r', 'LineWidth',1.5,'LineStyle', '--')
    hold on;
    plot(1:1:T, ref_store_lksdcee_n8(:,1), 'g :', 'LineWidth',1.5);
    hold on;
    plot(1:1:T, xi_k_store_lksdcee_n8(:,1), 'c-.', 'LineWidth',1.5);
    xlim([0,T]);
    my_legend = legend('output $y$','optimal reference', 'generated reference', 'trajectory $\xi$', '$ v_{kn} $');
    set(my_legend,'FontSize',12);
    xlabel('$ n $','FontSize',12);
    ylabel('$y$','FontSize',12);



    
