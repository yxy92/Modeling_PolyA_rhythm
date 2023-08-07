% Sobol_Indices_Calculation_Wiki_Latin_Sampling
% Based on V2 version, this version re-write the code such that 
% all parameters are stored using numberical array instead of struct array


clear;
rng('shuffle');

for id=1:10
    single_trial(id);
end

function single_trial(id)
%% parameter
output_dir = 'D:\CircadianRhythmicity\2019Fall\Sobol_Indices_Calculation_Version_for_Manuscript_V4\A1B1C1D1\Wiki_Sampling\Table\';

Input_Number = 11;
Output_Number = 15;

Sample_Size = 100000;

Method = 'Latin'; % 'Naive' or 'Latin'

ODE_Plot = 0; % 1 yes; 0 no

% 1) original parameter range
%distribution of the half-life time of deadenylation in log10
mu_deA=log10(2);%median value is 2
mu_polyA=log10(2);%median value is 2
% distribution of the time scale of degradation in log10
mu_dgrd=log10(9);%median value is 9
sigma=0.2275;% standard deviation of the normal distribution




%% generate parameter space
A= Sampling(Input_Number,Sample_Size,Method,mu_deA,mu_polyA,mu_dgrd,sigma);
B= Sampling(Input_Number,Sample_Size,Method,mu_deA,mu_polyA,mu_dgrd,sigma);
C= Sampling(Input_Number,Sample_Size,Method,mu_deA,mu_polyA,mu_dgrd,sigma);


%% mix parameter A and B to form AB
AB = zeros(Input_Number*Sample_Size,Input_Number);
for i =1:Input_Number
    AB((i-1)*Sample_Size+1:i*Sample_Size,:) = A;
    AB((i-1)*Sample_Size+1:i*Sample_Size,i) = B(:,i);
end


%% solve ODE for Ptot
parameter_tot = [A;B;AB;C];
result_tot = zeros(Sample_Size*(3+Input_Number),Input_Number+Output_Number);

parfor i=1:(3+Input_Number)*Sample_Size
    result_tot(i,:) = SolveODE(parameter_tot(i,:),Input_Number,Output_Number,ODE_Plot);
end


result_A = result_tot(1:Sample_Size,:);
result_B = result_tot(Sample_Size+1:2*Sample_Size,:);
result_AB = result_tot(2*Sample_Size+1:(2+Input_Number)*Sample_Size,:);
result_C = result_tot((2+Input_Number)*Sample_Size+1:(3+Input_Number)*Sample_Size,:);


%% calculate Sobol indices
% calculate average variance
Avg_Variance = zeros(1,Output_Number);

for i = 1:Output_Number
    Avg_Variance(i) = var(result_C(:,Input_Number+i));            
end

Avg_Variance=repmat(Avg_Variance,Input_Number,1);


%% calculate V1 and Vt
V1 = zeros(Input_Number,Output_Number);
Vt = zeros(Input_Number,Output_Number);

for i=1:Input_Number 
    for j=1:Output_Number
        V1(i,j) = dot(result_B(:,Input_Number+j),(result_AB((i-1)*Sample_Size+1:i*Sample_Size,Input_Number+j)-result_A(:,Input_Number+j)))/Sample_Size;
        Vt(i,j) = norm(result_A(:,Input_Number+j)-result_AB((i-1)*Sample_Size+1:i*Sample_Size,Input_Number+j))^2/(2*Sample_Size);       
    end
 end

%% last step calculating S1 and St
S1 = V1./Avg_Variance;
St = Vt./Avg_Variance;

% transform to table and add name tag
inputVariable = {'Atrsc','Ptrsc','KdeA','AdeA','PdeA','KpolyA','ApolyA','PpolyA','Kdgrd','Adgrd','Pdgrd'};
outputQuantity = {'meanls','rals','Pls','meanss','rass','Pss','meanL','raL','PL','meanS','raS','PS','meanSL','raSL','PSL'};

S1 = array2table(S1);
St = array2table(St);

S1.Properties.VariableNames = outputQuantity;
S1.Properties.RowNames = inputVariable;

St.Properties.VariableNames = outputQuantity;
St.Properties.RowNames = inputVariable;

%% save 2 table
writetable(S1,strcat(output_dir,'Trial','_',num2str(id),'SL_SingleSobol_','SampleSize','_',num2str(Sample_Size),'Method','_',Method,'.xlsx'),'WriteRowNames',true);
writetable(St,strcat(output_dir,'Trial','_',num2str(id),'SL_TotalSobol_','SampleSize','_',num2str(Sample_Size),'Method','_',Method,'.xlsx'),'WriteRowNames',true);




end



function  X = Sampling(Input_Number,Sample_Size,Method,mu_deA,mu_polyA,mu_dgrd,sigma)
% X --- a parameter set 

hl_deA=makedist('normal','mu',mu_deA,'sigma',sigma);
hl_polyA=makedist('normal','mu',mu_polyA,'sigma',sigma);
hl_dgrd=makedist('normal','mu',mu_dgrd,'sigma',sigma);

switch Method
    case 'Naive'
        X = rand(Sample_Size,Input_Number);        
    case 'Latin'
        X = lhsdesign(Sample_Size,Input_Number);
end


% transform sample cube to parameters
for i=1:Sample_Size
    X(i,1) = X(i,1);   % Atrsc
    X(i,2) = 24*X(i,2); % Ptrsc

    X(i,3) = log(2)/(10^icdf(hl_deA,X(i,3))); % KdeA
    X(i,4) = X(i,4); % AdeA
    X(i,5) = 24*X(i,5); % PdeA

    X(i,6) = log(2)/(10^icdf(hl_polyA,X(i,6))); % KpolyA       
    X(i,7) = X(i,7); % ApolyA
    X(i,8) = 24*X(i,8); % PpolyA
    
    X(i,9) = log(2)/(10^icdf(hl_dgrd,X(i,9))); % Kdgrd       
    X(i,10) = X(i,10); % Adgrd
    X(i,11) = 24*X(i,11); % Pdgrd
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


%% interpolation to find the peak for short-tailed mRNA
tS=linspace(600,648,2000);
Sq=interp1(t,S,tS,'pchip');
[maxS,maxSindex]=max(Sq);
[minS,~]=min(Sq);

result(Input_Number+10)=mean(Sq);
result(Input_Number+11)=(maxS-minS)/(maxS+minS);
result(Input_Number+12)=mod(tS(maxSindex),24);

%% interpolation to find the peak for SLR 
tSL=linspace(600,648,2000);
SLq=interp1(t,S./L,tSL,'pchip');
[maxSL,maxSLindex]=max(SLq);
[minSL,~]=min(SLq);

result(Input_Number+13)=mean(SLq);
result(Input_Number+14)=(maxSL-minSL)/(maxSL+minSL);
result(Input_Number+15)=mod(tSL(maxSLindex),24);


if ODE_Plot == 1
    %test plot
    scatter(t,S,300,'.');
    hold on;
    scatter(t,S./L,300,'.');
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









