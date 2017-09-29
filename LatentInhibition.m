%% Latent Inhibition

% Note: I think the only way to make the preexposed CS have better timing
% (measured by a higher CR slope) is to change the width of x.

clear
close all
clc
%% 
% Two CSs are used: CS A (preexposed) and CS B (not preexposed). 
%% Trial Structure
% The cycle is composed of CS, US and ITI in this order as follows:
% 
%  
% 
% Choose the duration of each CS (in milliseconds) using the variables below:

CS_dur=5; % (ms)


%% Parameters

h=0.01;
tau_x=50;
alpha_t=0.1;
alpha_E=0.08;
mu=1;
sigma=0.35;
m=0.2;
H=4;

alpha_sigma=0.025;



%% 
% Enter the initial value of associability $\alpha$ in the variable _alpha_ini_ 
% below.
%% 
%% Phase 1: Preexposure
% CS A will be preexposed (A-). Nothing will happen to CS B.
% 
% Enter the number of preexposure trials in the variable _cycle_num_ below.

cycle_num=80;
%% 
% 

%---Create structure arrays to carry CS info
CS(1).name='preexposed CS';
CS(2).name='control CS';
CS(1).pre_rate=NaN(cycle_num,1);
CS(1).pre_rate(1)=1*10^(-3);
CS(1).pre_V=zeros(cycle_num,1);
CS(1).pre_alpha=zeros(cycle_num,1);
CS(1).pre_alpha(1)=0.4;
CS(1).pre_sigma=zeros(cycle_num,1);
CS(1).pre_sigma(1)=0.6;
%---

%---DDM constants
N=normrnd(0,1,ceil(CS_dur/h)*cycle_num,1); % noise for CS
%---

cycle_length=round(CS_dur/h);

% %---Associative strength and alpha
% V=zeros(1,cycle_num); % initialize vector to store V values for time steps
% h_IDBD=zeros(1,cycle_num); % initialize vector to store alpha values for trials
% beta=zeros(1,cycle_num); 
% beta(1)=log(0.05);
gamma=0.03;
% alpha=zeros(1,cycle_num);
% alpha(1)=0.1;
% %---

%---Initialize counter for DDM noise
counterDDM=0;
%---


for trial=1:cycle_num
    
    %--initialize values for timer and CS
    P=zeros(1,cycle_length);
    x=zeros(1,cycle_length);
    %--
    
    for t=1:cycle_length
        
        counterDDM=counterDDM+1; % update counter for random process in DDM
        
        A=CS(1).pre_rate(trial);
        V=CS(1).pre_V(trial);
        alpha=CS(1).pre_alpha(trial);
        sigma=CS(1).pre_sigma(trial);
        
        % min will take the minimum value: either DDM result or 3. This
        % caps the value of integrator at 3.
        P(t+1)=min(DDM( P(t), A, h, m, N(counterDDM) ), 3);
        
        % max ensures the minimum value the accumulator can reach is
        % 1*10^-6. This avoids division by zero later.
        P(t+1)=max(P(t+1), 1*10^(-3));
        
        %---Element (RBFs)
        x(t)=CStrace(P(t+1),mu,sigma,tau_x,1,x(t),h);
        %---
        
        
    end
    
    %---Slope Correction
    A=A+A*alpha_t*(1-P(t+1))/P(t+1); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    %---
    
    %---V update
    [V, alpha]=RW_PH(V,alpha,x(t),alpha_E,0*H,A,P(t+1),gamma);
    %[ V(trial+1),h_IDBD(trial+1),beta(trial+1) ]=RW_IDBD( V(trial),h_IDBD(trial),x(t),beta(trial),H,A,P(t),gamma);
    %---
    
    %---Sigma update
    sigma=sigma+alpha_sigma*(0.35-sigma);
    %---
    
    CS(1).pre_rate(trial+1)=A;
    CS(1).pre_V(trial+1)=V;
    CS(1).pre_alpha(trial+1)=alpha;
    CS(1).pre_sigma(trial+1)=sigma;
    
end
%% 
% The figure below shows the evolution of associability $\alpha$ during 
% the entire preexposure training (including US and ITI periods):


%% 
% The next figure shows activation of the CS A representation:

%% Phase 2: Conditioning
% Both CSs A and B will be separately associated with a US (A+,B+). The initial 
% value of $\alpha$ for CS A will be the final value from preexposure, whilst 
% for CS B it will be the initial value from A before preexposure.
% 
% The timing of the representation for CS A will be improved since it has 
% already been tuned during preexposure.
% 
% Enter the number of conditioning trials in the variable _cycle_num_ below.
% 
cycle_num=250;

%---structure arrays
CS(1).cond_rate=NaN(cycle_num,1);
CS(1).cond_rate(1)=CS(1).pre_rate(end);
CS(1).cond_V=zeros(cycle_num,1);
CS(1).cond_alpha=zeros(cycle_num,1);
CS(1).cond_alpha(1)=CS(1).pre_alpha(end);
CS(2).cond_rate=NaN(cycle_num,1);
CS(2).cond_rate(1)=1*10^(-3);
CS(2).cond_V=zeros(cycle_num,1);
CS(2).cond_alpha=zeros(cycle_num,1);
CS(2).cond_alpha(1)=CS(1).pre_alpha(1);
CS(1).CR=NaN(cycle_num,cycle_length);
CS(2).CR=NaN(cycle_num,cycle_length);
CS(1).cond_sigma=zeros(cycle_num,1);
CS(1).cond_sigma(1)=CS(1).pre_sigma(end);
CS(2).cond_sigma=zeros(cycle_num,1);
CS(2).cond_sigma(1)=CS(1).pre_sigma(1);
%---

%---DDM constants
N_pre=normrnd(0,1,ceil(CS_dur/h)*cycle_num,1); % noise for CS A
N_con=normrnd(0,1,ceil(CS_dur/h)*cycle_num,1); % noise for CS B
%---

counterDDM=0;

for trial=1:cycle_num
    
    %--initialize values for timer and CS
    P_pre=zeros(1,cycle_length);
    x_pre=zeros(1,cycle_length);
    P_con=zeros(1,cycle_length);
    x_con=zeros(1,cycle_length);
    %--
    
    for t=1:cycle_length
        
        counterDDM=counterDDM+1; % update counter for random process in DDM
        
        A_pre=CS(1).cond_rate(trial);
        V_pre=CS(1).cond_V(trial);
        alpha_pre=CS(1).cond_alpha(trial);
        A_con=CS(2).cond_rate(trial);
        V_con=CS(2).cond_V(trial);
        alpha_con=CS(2).cond_alpha(trial);
        sigma_con=CS(2).cond_sigma(trial);
        sigma_pre=CS(1).cond_sigma(trial);
        
        % min will take the minimum value: either DDM result or 3. This
        % caps the value of integrator at 3.
        P_pre(t+1)=min(DDM( P_pre(t), A_pre, h, m, N_pre(counterDDM) ), 3);
        P_con(t+1)=min(DDM( P_con(t), A_con, h, m, N_con(counterDDM) ), 3);

        % max ensures the minimum value the accumulator can reach is
        % 1*10^-6. This avoids division by zero later.
        P_pre(t+1)=max(P_pre(t+1), 1*10^(-3));
        P_con(t+1)=max(P_con(t+1), 1*10^(-3));
        
        %---Element (RBFs)
        x_pre(t)=CStrace(P_pre(t+1),mu,sigma,tau_x,1,x_pre(t),h);
        x_con(t)=CStrace(P_con(t+1),mu,sigma,tau_x,1,x_con(t),h);
        %---
        
        %---CR
        CS(1).CR(trial,t)=x_pre(t)*V_pre;
        CS(2).CR(trial,t)=x_con(t)*V_con;
        %---
       
        
    end
    
    %---Slope Correction
    A_pre=A_pre+A_pre*alpha_t*(1-P_pre(t+1))/P_pre(t+1); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    A_con=A_con+A_con*alpha_t*(1-P_con(t+1))/P_con(t+1); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    %---
    
    %---V update
    [V_pre, alpha_pre]=RW_PH(V_pre,alpha_pre,x_pre(t),alpha_E,H,A_pre,P_pre(t+1),gamma);
    [V_con, alpha_con]=RW_PH(V_con,alpha_con,x_con(t),alpha_E,H,A_con,P_con(t+1),gamma);
    %[ V_pre,h_IDBD(trial+1),beta(trial+1) ]=RW_IDBD( V(trial),h_IDBD(trial),x(t),beta(trial),H,A,P(t),gamma);
    %---
    
    %---Sigma update
    sigma_pre=sigma_pre+alpha_sigma*(0.35-sigma_pre);
    sigma_con=sigma_con+alpha_sigma*(0.35-sigma_con);
    %---
    
    CS(1).cond_rate(trial+1)=A_pre;
    CS(1).cond_V(trial+1)=V_pre;
    CS(1).cond_alpha(trial+1)=alpha_pre;
    CS(2).cond_rate(trial+1)=A_con;
    CS(2).cond_V(trial+1)=V_con;
    CS(2).cond_alpha(trial+1)=alpha_con;
    CS(1).cond_sigma(trial+1)=sigma_pre;
    CS(2).cond_sigma(trial+1)=sigma_con;
    
end

figure
plot(1:length(CS(1).cond_V),CS(1).cond_V,1:length(CS(2).cond_V),CS(2).cond_V,'LineWidth',5)
legend({'Preexposed CS','Control CS'},'Location','southeast','Box','off')
PlotProperties
xlabel('trial (phase 2)')
ylabel('associative strength')
%axis([0 250 0 1.5])

% alpha_preexposed=vertcat(CS(1).pre_alpha,CS(1).cond_alpha);
% 
% figure
% plot(alpha_preexposed,'LineWidth',3)
% xlabel('trial')
% ylabel('associability \alpha')
% PlotProperties

CR_pre_early_mean=mean(CS(1).CR(1:30,:),1);
CR_con_early_mean=mean(CS(2).CR(1:30,:),1);
figure
plot(h:h:length(CR_pre_early_mean)*h,CR_pre_early_mean,h:h:length(CR_con_early_mean)*h,CR_con_early_mean,'LineWidth',6)
PlotProperties
legend({'Preexposed CS','Control CS'},'Location','northwest','Box','off')
title('model')
yticks(0:0.1:0.5)
xlabel('time (s)')
ylabel('response strength')

CR_pre_late_mean=mean(CS(1).CR(170:200,:),1);
CR_con_late_mean=mean(CS(2).CR(170:200,:),1);
figure
plot(h:h:length(CR_pre_late_mean)*h,CR_pre_late_mean,h:h:length(CR_con_late_mean)*h,CR_con_late_mean,'LineWidth',6)
PlotProperties
legend({'Preexposed CS','Control CS'},'Location','southeast','Box','off')
title('model')
xlabel('time (s)')
ylabel('response strength')

save('LatentInhSim','CS','CR_pre_early_mean','CR_pre_late_mean','CR_con_early_mean','CR_con_late_mean','h')

clear
% % Here we plot the evolution of associative strength $V$ for both CSs.
% 
% plot(1:length(Vavg_trial_A),Vavg_trial_A,1:length(Vavg_trial_B),Vavg_trial_B,'LineWidth',3)
% PlotProperties
% xlabel('trials')
% ylabel('associative strength')
% legend('CS A','CS B')
% %% 
% % The next plot shows the evolution of associability $\alpha$ for both CSs 
% % during the entire Phase 2 (including US and ITI periods).
% 
% alphatrials_A=reshape(alphatrials_A',size(alphatrials_A,1)*size(alphatrials_A,2),[]);
% alphatrials_B=reshape(alphatrials_B',size(alphatrials_B,1)*size(alphatrials_B,2),[]);
% plot(h/1000:h/1000:length(alphatrials_A)*h/1000,alphatrials_A,h/1000:h/1000:length(alphatrials_B)*h/1000,alphatrials_B,'LineWidth',3)
% PlotProperties
% xlabel('time (sec)')
% ylabel('associability \alpha')
% legend('\alpha CS A','\alpha CS B')