% Sobol_Indices_Calculation_Nested_Sampling_V4.m
% for model A1B1C1D1
clear;
rng('shuffle');

for id=1:10
    single_trial(id);
end


function single_trial(id)

%% parameter
Group_Number = 25000;
Group_Size = 4;

Input_Number = 11;
Output_Number = 9;

Sample_Total = 100000;

Method = 'Latin'; 

ODE_Plot = 0; % 1 yes; 0 no

output_dir = 'D:\CircadianRhythmicity\2019Fall\Sobol_Indices_Calculation_Version_for_Manuscript_V4\A1B1C1D1\Nested_Sampling\Table\';

%distribution of the half-life time of deadenylation in log10
mu_deA=log10(2);%median value is 2
% distribution of the time scale of polyadenylation in log10
mu_polyA=log10(2);%median value is 2
% distribution of the time scale of degradation in log10
mu_dgrd=log10(9);%median value is 9
sigma=0.2275;% standard deviation of the normal distribution



%% generate parameter space
[Xi,Xni,Xtot] = Nested_Sampling_Cube(Input_Number,Group_Number,Group_Size,Sample_Total,Method);
[Pi, Pni,Ptot] = Transform_Cube_to_Parameter(Xi, Xni,Xtot,mu_deA,mu_polyA,mu_dgrd,sigma);


%% solve ODE for Ptot

result_tot = zeros(Sample_Total,Input_Number+Output_Number);
parfor i=1:Sample_Total
    result_tot(i,:) = SolveODE(Ptot(i,:),Input_Number,Output_Number,ODE_Plot);
end

%% calculate Sobol indices
% calculate average variance
Avg_Variance = zeros(1,Output_Number);

for i = 1:Output_Number
    switch i
        case {1 2 4 5 7 8} % non-circular output
            Avg_Variance(i) = var(result_tot(:,Input_Number+i));            
        case {3 6 9} % circular output
            Avg_Variance(i) = 1 - circ_r(result_tot(:,Input_Number+i)*pi/12)^2; % 1- R^2
    end   
end

Avg_Variance=repmat(Avg_Variance,Input_Number,1);

%% Solve ODE for Single and Total Sobol
result_i = zeros(Group_Number*Group_Size,Input_Number+Output_Number,Input_Number);
result_ni = zeros(Group_Number*Group_Size,Input_Number+Output_Number,Input_Number);

for i=1:Input_Number
    tmp_i = Pi(:,:,i);
    tmp_ni = Pni(:,:,i);
    
    parfor j = 1: Group_Number*Group_Size
        result_i(j,:,i) = SolveODE(tmp_i(j,:),Input_Number,Output_Number,ODE_Plot)
        result_ni(j,:,i) = SolveODE(tmp_ni(j,:),Input_Number,Output_Number,ODE_Plot)
    end
end



%% calculate V1 and Vt
V1 = zeros(Input_Number,Output_Number);
Vt = zeros(Input_Number,Output_Number);

for i=1:Input_Number
    tmp_1 = result_i(:,:,i);
    tmp_t = result_ni(:,:,i);
    
    for j=1:Output_Number       
       
        group_mean_1 = zeros(1,Group_Number);
        group_mean_t = zeros(1,Group_Number);
        
        switch j               
            case {1 2 4 5 7 8} % non-circular output
                for k =1:Group_Number
                    % extract data for each group
                    data_1 = tmp_1((k-1)*Group_Size+1:k*Group_Size,Input_Number+j);
                    data_t = tmp_t((k-1)*Group_Size+1:k*Group_Size,Input_Number+j);                    

                    % take mean of each group
                    group_mean_1(k) = mean(data_1);
                    group_mean_t(k) = mean(data_t);
                end
                
                V1(i,j) = var(group_mean_1);
                Vt(i,j) = var(group_mean_t);
                             
            case {3 6 9} % circular output
               for k =1:Group_Number
                    % extract data for each group
                    data_1 = tmp_1((k-1)*Group_Size+1:k*Group_Size,Input_Number+j);
                    data_t = tmp_t((k-1)*Group_Size+1:k*Group_Size,Input_Number+j);                    

                    % take mean of each group
                    group_mean_1(k) = circ_r(data_1*pi/12)^2; 
                    group_mean_t(k) = circ_r(data_t*pi/12)^2;
                end
                
                V1(i,j) = mean(group_mean_1)-(1-Avg_Variance(i,j));
                Vt(i,j) = mean(group_mean_t)-(1-Avg_Variance(i,j)); 
                
                
        end
    end
 end

%% last step calculating S1 and St
S1 = V1./Avg_Variance;
St = 1-Vt./Avg_Variance;

% transform to table and add name tag
inputVariable = {'Atrsc','Ptrsc','KdeA','AdeA','PdeA','KpolyA','ApolyA','PpolyA','Kdgrd','Adgrd','Pdgrd'};
outputQuantity = {'meanls','rals','Pls','meanss','rass','Pss','meanL','raL','PL'};

S1 = array2table(S1);
St = array2table(St);

S1.Properties.VariableNames = outputQuantity;
S1.Properties.RowNames = inputVariable;

St.Properties.VariableNames = outputQuantity;
St.Properties.RowNames = inputVariable;

% save 2 table
writetable(S1,strcat(output_dir,'Trial','_',num2str(id),'_','SingleSobol_','GroupNumber',num2str(Group_Number),'_','GroupSize_',num2str(Group_Size),...
    '_','Method','_',Method,'.xlsx'),'WriteRowNames',true);

writetable(St,strcat(output_dir,'Trial','_',num2str(id),'_','TotalSobol_','GroupNumber',num2str(Group_Number),'_','GroupSize_',num2str(Group_Size),...
    '_','Method','_',Method,'.xlsx'),'WriteRowNames',true);

end


%% functions

function [Pi, Pni, Ptot] = Transform_Cube_to_Parameter(Xi, Xni, Xtot,mu_deA,mu_polyA,mu_dgrd,sigma)
hl_deA=makedist('normal','mu',mu_deA,'sigma',sigma);
hl_polyA=makedist('normal','mu',mu_polyA,'sigma',sigma);
hl_dgrd=makedist('normal','mu',mu_dgrd,'sigma',sigma);

dim_i = size(Xi);

Pi = Xi;
Pni = Xni;

for i=1:dim_i(3) % Number of input
    for j=1:dim_i(1) % Number Group * Number of sampling within each group
        % single sobol matrix
        Pi(j,1,i) = Pi(j,1,i);   % Atrsc
        Pi(j,2,i) = 24*Pi(j,2,i); % Ptrsc
        
        Pi(j,3,i) = log(2)/(10^icdf(hl_deA,Pi(j,3,i))); % KdeA
        Pi(j,4,i) = Pi(j,4,i); % AdeA
        Pi(j,5,i) = 24*Pi(j,5,i); % PdeA
        
        Pi(j,6,i) = log(2)/(10^icdf(hl_polyA,Pi(j,6,i))); % KpolyA       
        Pi(j,7,i) = Pi(j,7,i); % ApolyA
        Pi(j,8,i) = 24*Pi(j,8,i); % PpolyA
        
        Pi(j,9,i) = log(2)/(10^icdf(hl_dgrd,Pi(j,9,i))); % Kdgrd
        Pi(j,10,i) = Pi(j,10,i); % Adgrd
        Pi(j,11,i) = 24*Pi(j,11,i); % Pdgrd
        
        
        %% total sobol matrix
        Pni(j,1,i) = Pni(j,1,i);   % Atrsc
        Pni(j,2,i) = 24*Pni(j,2,i); % Ptrsc
        
        Pni(j,3,i) = log(2)/(10^icdf(hl_deA,Pni(j,3,i))); % KdeA
        Pni(j,4,i) = Pni(j,4,i); % AdeA
        Pni(j,5,i) = 24*Pni(j,5,i); % PdeA
        
                
        Pni(j,6,i) = log(2)/(10^icdf(hl_polyA,Pni(j,6,i))); % KpolyA
        Pni(j,7,i) = Pni(j,7,i); % ApolyA
        Pni(j,8,i) = 24*Pni(j,8,i); % PpolyA
        
        Pni(j,9,i) = log(2)/(10^icdf(hl_dgrd,Pni(j,9,i))); % Kdgrd       
        Pni(j,10,i) = Pni(j,10,i); % Adgrd
        Pni(j,11,i) = 24*Pni(j,11,i); % Pdgrd
        
    end
end


dim_tot = size(Xtot);
Ptot = Xtot;
for i=1:dim_tot(1)
    % single sobol matrix
    Ptot(i,1) = Ptot(i,1);   % Atrsc
    Ptot(i,2)= 24*Ptot(i,2); % Ptrsc

    Ptot(i,3) = log(2)/(10^icdf(hl_deA,Ptot(i,3))); % KdeA
    Ptot(i,4) = Ptot(i,4); % AdeA
    Ptot(i,5) = 24*Ptot(i,5); % PdeA

    Ptot(i,6) = log(2)/(10^icdf(hl_polyA,Ptot(i,6))); % KpolyA       
    Ptot(i,7) = Ptot(i,7); % ApolyA
    Ptot(i,8) = 24*Ptot(i,8); % PpolyA
    
    Ptot(i,9) = log(2)/(10^icdf(hl_dgrd,Ptot(i,9))); % Kdgrd       
    Ptot(i,10) = Ptot(i,10); % Adgrd
    Ptot(i,11) = 24*Ptot(i,11); % Pdgrd
        
end

end



function  [Xi, Xni, Xtot] = Nested_Sampling_Cube(Input_Number,Group_Number,Group_Size,Sample_Total,Method)

% Xi --- parameter set for single sobol indices
% Xni --- parameter set for total sobol indices
% Xtot --- parameter set for average variance calculation

% for each input variable in Pi, the ith column is fixed within each group
% while other columns are random
Xi = zeros(Group_Number*Group_Size,Input_Number,Input_Number);

switch Method
    case 'Naive'     
        for i = 1:Input_Number
            X_fix = reshape(repmat(rand(1,Group_Number),Group_Size,1),[],1);
            X_var = rand(Group_Number*Group_Size,Input_Number-1); % naive random
            Xi(:,:,i) = [X_var(:,1:i-1),X_fix, X_var(:,i:end)];
        end
        
    case 'Latin'          
        for i = 1:Input_Number
            X_fix = reshape(repmat(rand(1,Group_Number),Group_Size,1),[],1);
            X_var = lhsdesign(Group_Number*Group_Size,Input_Number-1); % call latin hyper cube
            Xi(:,:,i) = [X_var(:,1:i-1),X_fix, X_var(:,i:end)];
        end
end


% for each input variable in Pni, the ith column is random while other
% columns are fixed within each group
Xni = zeros(Group_Number*Group_Size,Input_Number,Input_Number);

switch Method
    case 'Naive'
        for i =1:Input_Number
            X_fix = zeros(Group_Number*Group_Size,Input_Number-1);
%             % Jing's code
%             for k = 1:Input_Number-1
%                 X_fix(:,k) = reshape(repmat(rand(1,Group_Number),Group_Size,1),[],1);
%             end

%           % XY's code
            tmp=rand(Group_Number,Input_Number-1);
            for k =1:Group_Number
                X_fix((k-1)*Group_Size+1:k*Group_Size,:) = repmat(tmp(k,:),Group_Size,1);
            end
        
            X_var = rand(Group_Number*Group_Size,1);
            Xni(:,:,i) = [X_fix(:,1:i-1),X_var,X_fix(:,i:end)];
        end
         
    case 'Latin'
         for i =1:Input_Number
             X_fix = zeros(Group_Number*Group_Size,Input_Number-1);
            % XY's code
            tmp=lhsdesign(Group_Number,Input_Number-1); % latin hyper cube
            for k =1:Group_Number
                X_fix((k-1)*Group_Size+1:k*Group_Size,:) = repmat(tmp(k,:),Group_Size,1);
            end
        
            X_var = rand(Group_Number*Group_Size,1);
            Xni(:,:,i) = [X_fix(:,1:i-1),X_var,X_fix(:,i:end)];           
             
             
         end
        
end


% for all inputs 

switch Method
    case 'Naive'
        Xtot = rand(Sample_Total,Input_Number);
    case 'Latin'
        Xtot = lhsdesign(Sample_Total,Input_Number);
end

end



%% given a parameter set, solve ODE, return the structure p containing
% parameter + result
function result = SolveODE(p,Input_Number,Output_Number,ODE_Plot)
result = zeros(1,Input_Number+Output_Number);
result(1:Input_Number) = p;

% simulation time interval
tspan = [0 700];
% intial value
x0 = [0;0];

Tol=1e-5;
opts = odeset('RelTol',Tol,'AbsTol',Tol);

sol = ode45(@(t,x) model_A1B1C1D1(t,x,p),tspan,x0,opts);
t = linspace(600,648,500); % extract a 48 h window
x=deval(sol,t);

L = x(1,:);% long-tailed
S = x(2,:);% short-tailed


%% interpolation to find the peak for LSR
LSR = L./S; % LS ratio
tls=linspace(600,648,2000);
lsq=interp1(t,LSR,tls,'pchip');
[maxls,maxlsindex]=max(lsq);
[minls,~]=min(lsq);

% calculate LSR rhythmicity
rals = (maxls-minls)/(maxls+minls);
meanls =  mean(lsq);
Pls=mod(tls(maxlsindex),24);

result(Input_Number+1) = meanls;
result(Input_Number+2) = rals;
result(Input_Number+3) = Pls;

%% interpolation to find the peak for total mRNA
ssmRNA = L+S; % total mRNA
tss=linspace(600,648,2000);
ssq=interp1(t,ssmRNA,tss,'pchip');
[maxss,maxssindex]=max(ssq);
[minss,~]=min(ssq);

% calculate total mRNA rhythmicity
result(Input_Number+4) =  mean(ssq);
result(Input_Number+5) = (maxss-minss)/(maxss+minss);
result(Input_Number+6)= mod(tss(maxssindex),24);


%% interpolation to find the peak for long-tailed mRNA
tL=linspace(600,648,2000);
Lq=interp1(t,L,tL,'pchip');
[maxL,maxLindex]=max(Lq);
[minL,~]=min(Lq);

result(Input_Number+7)=mean(Lq);
result(Input_Number+8)=(maxL-minL)/(maxL+minL);
result(Input_Number+9)=mod(tL(maxLindex),24);

if ODE_Plot == 1
    %test plot
    scatter(t,L./S,300,'.');
    hold on;
    scatter(t,L+S,300,'.');
    hold on;
    scatter(t,L,300,'.');
    legend('LSR','total mRNA','long-tailed');
end

end





function dxdt = model_A1B1C1D1(t,x,p)
% pass argument to parameter
omg = 2*pi/24;

Ktrsc = 1;
Atrsc = p(1);
Ptrsc= p(2);

KdeA= p(3);
AdeA = p(4);
PdeA = p(5);

KpolyA = p(6);
ApolyA = p(7);
PpolyA = p(8);

Kdgrd = p(9);
Adgrd = p(10);
Pdgrd = p(11);

dxdt = [0; 0];
L = x(1);
S = x(2);

dxdt(1) = Ktrsc*(1+Atrsc*cos(omg*(t-Ptrsc)))-KdeA*(1+AdeA*cos(omg*(t-PdeA)))*L+KpolyA*(1+ApolyA*cos(omg*(t-PpolyA)))*S;
dxdt(2) = KdeA*(1+AdeA*cos(omg*(t-PdeA)))*L -KpolyA*(1+ApolyA*cos(omg*(t-PpolyA)))*S -Kdgrd*(1+Adgrd*cos(omg*(t-Pdgrd)))*S;

end









