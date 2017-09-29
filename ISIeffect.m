%% ISI Effect

clear
close all

% Choose how many CSs are going to be used by entering their durations in 
% the vector _intervals._ For example, for two CSs, one that lasts 1000 ms and 
% the other 2000 ms enter: intervals=[1000 2000].

intervals=[5 10 20]; % duration of each CS (ms)

% Choose the other durations involved in this experiment:
cycle_num=150; % number of cycles per ISI


%---parameters
h=0.01;
tau_x=50;
alpha_t=0.2;
alpha_E=0.1;
mu=1;
sigma=0.3;
m=0.15;
H=30;
%---

%---Cell for saving CRs per ISI
CR_ISI=cell(1,length(intervals));
%---

%---Matrix for saving V per ISI
V_ISI=zeros(length(intervals),cycle_num);
%---

for ISI=1:length(intervals)
    
    %--Durations
    CS_dur=intervals(ISI); % (ms)
    cycle_length=round(CS_dur/h);
    CR=zeros(cycle_num,cycle_length);
    %--
    
    %---DDM constants
    A=1*10^(-3);
    N=normrnd(0,1,ceil(CS_dur/h)*cycle_num,1); % noise for CS
    %---
    
    %---Associative strength and alpha
    V=zeros(1,cycle_num); % initialize vector to store V values for trials
    %---
    
    %---Setting the presence or absence of CS and US
    CS=ones(1,cycle_length); % CS is 1 only during CS presentation
    %---
    
    %---Initialize counter for DDM noise
    counterDDM=0;
    %---
    
    for trial=1:cycle_num
        
        %--initialize values for timer and CS 
        P=zeros(1,cycle_length);
        x=zeros(1,cycle_length);
        z=0;
        %--
        
        for t=1:cycle_length
            
            counterDDM=counterDDM+1; % update counter for random process in DDM
            
            % min will take the minimum value: either DDM result or 3. This
            % caps the value of integrator at 3.
            P(t+1)=CS(t)*min(DDM( P(t), A, h, m, N(counterDDM) ), 3);
            
            % max ensures the minimum value the accumulator can reach is
            % 1*10^-6. This avoids division by zero later.
            P(t+1)=max(P(t+1), 1*10^(-3));
            
            %---Element (RBFs)
            x(t)=CStrace(P(t+1),mu,sigma,tau_x,CS(t),x(t),h);
            %---
            
            %---Response
            CR(trial,t)=CS(t)*x(t)*V(trial)*(x(t)*V(trial)>=0);
            %---
        end
        
        %---Slope Correction
        A=A+A*alpha_t*(1-P(t))/P(round(t)); % realistic correction rule, never fully converges. Only updates in rewarded trials.
        %---
         
        %---V update
        V(trial+1)=RW(V(trial),alpha_E,x(t),H,A,P(t));
        %V(trial+1)=RW(V(trial),alpha_E,alpha_I,x(t),H,A,P(t));
        %---
        
    end
    
    %---Save CRs from previous ISI
    CR_ISI{ISI}=CR;
    %---
    
    %---Save V from previous ISI
    V_ISI(ISI,:)=V(1:end-1);
    %---
    
end

%%


labels={'FI 5','FI 10','FI 20'};
figure('name','ISI effect')
hold on
%---Response frequency (average over trials)
Avg_CR=cell(1,length(intervals));
for i=1:length(intervals)
    Avg_CR{i}=mean(CR_ISI{1,i}(round(cycle_num/2):end,:));
    plot(h:h:length(Avg_CR{i})*h,Avg_CR{i},'LineWidth',6)
    PlotProperties
end
legend(labels,'Box','off')
xlabel('time (s)')
ylabel('response strength')

%---

figure('name','Superposition')
hold on
for i=1:length(intervals)
    plot((1:intervals(i)/h)/(intervals(i))*h,Avg_CR{i}/max(Avg_CR{i}),'LineWidth',6)
    PlotProperties
end
%legend(labels)
xlabel('proportion of interval')
ylabel('norm. resp. strength')

figure('name','Associative Strength')
hold on
for i=1:length(intervals)
    plot(1:cycle_num,V_ISI(i,:),'LineWidth',6)
    PlotProperties
end
%legend(labels)
xlabel('trials')
ylabel('associative strength')

save('ISIeffectSim','Avg_CR','V_ISI','intervals','CR_ISI','cycle_num','h')

clear

% figure('name','Non-Averaged Responses')
% hold on
% for i=1:20:cycle_num
%     plot(CR_ISI{1,1}(i,:),'LineWidth',6)
%     PlotProperties
% end