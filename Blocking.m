% Blocking with different CS durations
clear
close all

cycle_num_ph1=120;
cycle_num_ph2=60;

CS1_dur=15;
CS2_dur=10;

%---parameters
h=0.01;
tau_x=50;
alpha_t=0.1;
alpha_E=0.08;
mu=1;
sigma=0.35;
m=0.2;
H=10;
%---

%---structure arrays
CS(1).name='blocking CS';
CS(2).name='blocked CS';
CS(3).name='control CS';
%-A
CS(1).A.ph1=zeros(cycle_num_ph1,1);
CS(1).A.ph2=zeros(cycle_num_ph2,1);
CS(1).A.ph1(1)=1*10^(-3);

CS(2).A.ph1=zeros(cycle_num_ph1,1);
CS(2).A.ph2=zeros(cycle_num_ph2,1);
CS(2).A.ph3=zeros(cycle_num_ph2,1);
CS(2).A.ph1(1)=1*10^(-3);

CS(3).A.ph3=zeros(cycle_num_ph2,1);
CS(3).A.ph3(1)=1*10^(-3);
%-V
CS(1).V.ph1=zeros(cycle_num_ph1,1);
CS(1).V.ph2=zeros(cycle_num_ph2,1);

CS(2).V.ph1=zeros(cycle_num_ph1,1);
CS(2).V.ph2=zeros(cycle_num_ph2,1);
CS(2).V.ph3=zeros(cycle_num_ph2,1);

CS(3).V.ph3=zeros(cycle_num_ph2,1);
%---


%% Phase 1


%---DDM constants
N=normrnd(0,1,ceil((CS1_dur/h)*cycle_num_ph1),1); % noise for CS A
%---

counterDDM=0;

cycle_length=round(CS1_dur/h);

for trial=1:cycle_num_ph1
    
    %--initialize values for timer and CS
    P=zeros(1,cycle_length);
    x=zeros(1,cycle_length);
    %--
    
    for t=1:cycle_length
        
        counterDDM=counterDDM+1; % update counter for random process in DDM
        
        A=CS(1).A.ph1(trial);
        V=CS(1).V.ph1(trial);
        
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
    V=RW(V,alpha_E,x(t),H,A,P(t+1));
    %---
    
    CS(1).A.ph1(trial+1)=A;
    CS(1).V.ph1(trial+1)=V;
    
end

figure('Name','Phase 1, V','NumberTitle','off')
plot(1:length(CS(1).V.ph1),CS(1).V.ph1,'LineWidth',6)
legend({'blocking CS'},'Box','off')
xlabel('trial')
ylabel('associative strength')
PlotProperties

% figure('Name','Phase 1, A','NumberTitle','off')
% plot(1:length(CS(1).A.ph1),CS(1).A.ph1,'LineWidth',6)
% legend('blocking CS')
% xlabel('trial')
% ylabel('slope A')
% PlotProperties

%% Phase 2

CS(1).A.ph2(1)=CS(1).A.ph1(end);
CS(1).V.ph2(1)=CS(1).V.ph1(end);
CS(2).A.ph2(1)=CS(1).A.ph1(1);

if CS1_dur>CS2_dur
    cycle_length=round(CS1_dur/h);
    CS1=ones(cycle_length,1);
    CS2=ones(cycle_length,1);
    CS2(1:round(CS1_dur/h)-round(CS2_dur/h))=zeros;
else
    cycle_length=round(CS2_dur/h);
    CS1=ones(cycle_length,1);
    CS2=ones(cycle_length,1);
    CS1(1:round(CS2_dur/h)-round(CS1_dur/h))=zeros;
end

%---DDM noise
N_1=normrnd(0,1,ceil(cycle_length*cycle_num_ph2),1); % noise for CS A
N_2=normrnd(0,1,ceil(cycle_length*cycle_num_ph2),1); % noise for CS A
%---

counterDDM=0;

for trial=1:cycle_num_ph2
    
    %--initialize values for timer and CS
    P_1=zeros(1,cycle_length);
    x_1=zeros(1,cycle_length);
    P_2=zeros(1,cycle_length);
    x_2=zeros(1,cycle_length);
    %--
    
    for t=1:cycle_length
        
        counterDDM=counterDDM+1; % update counter for random process in DDM
        
        A_1=CS(1).A.ph2(trial);
        V_1=CS(1).V.ph2(trial);
        
        A_2=CS(2).A.ph2(trial);
        V_2=CS(2).V.ph2(trial);

        
        % min will take the minimum value: either DDM result or 3. This
        % caps the value of integrator at 3.
        P_1(t+1)=CS1(t)*min(DDM( P_1(t), A_1, h, m, N_1(counterDDM) ), 3);
        P_2(t+1)=CS2(t)*min(DDM( P_2(t), A_2, h, m, N_2(counterDDM) ), 3);

        % max ensures the minimum value the accumulator can reach is
        % 1*10^-6. This avoids division by zero later.
        P_1(t+1)=max(P_1(t+1), 1*10^(-3));
        P_2(t+1)=max(P_2(t+1), 1*10^(-3));
        
        %---Element (RBFs)
        x_1(t)=CStrace(P_1(t+1),mu,sigma,tau_x,1,x_1(t),h);
        x_2(t)=CStrace(P_2(t+1),mu,sigma,tau_x,1,x_2(t),h);
        %---
        
        %---Timer selection
        P_att=min((P_1(t+1)-1)^2,(P_2(t+1)-1)^2);
        US_1_att=1*((P_1(t+1)-1)^2==P_att);
        US_2_att=1*((P_2(t+1)-1)^2==P_att);
        %---
        
        %---CR
        %CR(trial,t)=(x_1(t)*V_1(trial)*US_1_att+x_2(t)*V_2(trial)*US_2_att);
        %---
        
    end
    
    %---Slope Correction
    A_1=A_1+A_1*alpha_t*(1-P_1(t+1))/P_1(t+1); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    A_2=A_2+A_2*alpha_t*(1-P_2(t+1))/P_2(t+1); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    %---
    
    %---V update
    V_1=RW([V_1 V_2],alpha_E,US_1_att*[x_1(t) x_2(t)],H,A_1,P_1(t));
    V_2=RW([V_2 V_1],alpha_E,US_2_att*[x_2(t) x_1(t)],H,A_2,P_2(t));
    %---

    CS(1).A.ph2(trial+1)=A_1;
    CS(1).V.ph2(trial+1)=V_1;
    
    CS(2).A.ph2(trial+1)=A_2;
    CS(2).V.ph2(trial+1)=V_2;
end

%% Phase 3 - control
% CS2 is the same as in phase 2
% CS3 is a control (novel) CS that will take the place of CS1
% Otherwise phase 3 is exactly the same as phase 2

CS(2).A.ph3(1)=CS(1).A.ph1(1);

if CS1_dur>CS2_dur
    cycle_length=round(CS1_dur/h);
    CS3=ones(cycle_length,1);
    CS2=ones(cycle_length,1);
    CS2(1:round(CS1_dur/h)-round(CS2_dur/h))=zeros;
else
    cycle_length=round(CS2_dur/h);
    CS3=ones(cycle_length,1);
    CS2=ones(cycle_length,1);
    CS3(1:round(CS2_dur/h)-round(CS1_dur/h))=zeros;
end

%---DDM noise
N_3=normrnd(0,1,ceil(cycle_length*cycle_num_ph2),1); % noise for CS A
N_2=normrnd(0,1,ceil(cycle_length*cycle_num_ph2),1); % noise for CS A
%---

counterDDM=0;

for trial=1:cycle_num_ph2
    
    %--initialize values for timer and CS
    P_3=zeros(1,cycle_length);
    x_3=zeros(1,cycle_length);
    P_2=zeros(1,cycle_length);
    x_2=zeros(1,cycle_length);
    %--
    
    for t=1:cycle_length
        
        counterDDM=counterDDM+1; % update counter for random process in DDM
        
        A_3=CS(3).A.ph3(trial);
        V_3=CS(3).V.ph3(trial);
        
        A_2=CS(2).A.ph3(trial);
        V_2=CS(2).V.ph3(trial);

        
        % min will take the minimum value: either DDM result or 3. This
        % caps the value of integrator at 3.
        P_3(t+1)=CS3(t)*min(DDM( P_3(t), A_3, h, m, N_3(counterDDM) ), 3);
        P_2(t+1)=CS2(t)*min(DDM( P_2(t), A_2, h, m, N_2(counterDDM) ), 3);

        % max ensures the minimum value the accumulator can reach is
        % 1*10^-6. This avoids division by zero later.
        P_3(t+1)=max(P_3(t+1), 1*10^(-3));
        P_2(t+1)=max(P_2(t+1), 1*10^(-3));
        
        %---Element (RBFs)
        x_3(t)=CStrace(P_3(t+1),mu,sigma,tau_x,1,x_3(t),h);
        x_2(t)=CStrace(P_2(t+1),mu,sigma,tau_x,1,x_2(t),h);
        %---
        
        %---Timer selection
        P_att=min((P_3(t+1)-1)^2,(P_2(t+1)-1)^2);
        US_3_att=1*((P_3(t+1)-1)^2==P_att);
        US_2_att=1*((P_2(t+1)-1)^2==P_att);
        %---
        
        %---CR
        %CR(trial,t)=(x_1(t)*V_1(trial)*US_1_att+x_2(t)*V_2(trial)*US_2_att);
        %---
        
    end
    
    %---Slope Correction
    A_3=A_3+A_3*alpha_t*(1-P_3(t+1))/P_3(t+1); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    A_2=A_2+A_2*alpha_t*(1-P_2(t+1))/P_2(t+1); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    %---
    
    %---V update
    V_3=RW([V_3 V_2],alpha_E,US_3_att*[x_3(t) x_2(t)],H,A_3,P_3(t));
    V_2=RW([V_2 V_3],alpha_E,US_2_att*[x_2(t) x_3(t)],H,A_2,P_2(t));
    %---

    CS(3).A.ph3(trial+1)=A_3;
    CS(3).V.ph3(trial+1)=V_3;
    
    CS(2).A.ph3(trial+1)=A_2;
    CS(2).V.ph3(trial+1)=V_2;
end

figure('Name','Phase 2','NumberTitle','off')
plot(1:length(CS(2).V.ph2),CS(2).V.ph2,1:length(CS(2).V.ph3),CS(2).V.ph3,'LineWidth',6)
legend({'blocked CS','control CS'},'Location','southeast','Box','off')
title('model')
xlabel('trial')
ylabel('associative strength')
PlotProperties

save('BlockingSimLongBlocking','CS')