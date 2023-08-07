% errorbar_sobol_indices_big_single_together_total.m
% --------------------------------------------
% For sample variance, the degree of freedom is n-1 instead of n
%
%
% Input dir/
% D:\CircadianRhythmicity\2019Fall\Sobol_Indices_Calculation_Version_for_Manuscript_V4\A1B1C1D1\Wiki_Sampling\Table

clear;

%% read source table from directories

input_dir = 'D:\CircadianRhythmicity\2019Fall\Sobol_Indices_Calculation_Version_for_Manuscript_V4\A1B1C1D1\Wiki_Sampling\Table\';
output_dir = 'D:\CircadianRhythmicity\2019Fall\Sobol_Indices_Calculation_Version_for_Manuscript_V4\A1B1C1D1\Wiki_Sampling\Supplemental_Sobol_Figure\';
%% box plot of sobol indices results from 3 trials
trial_num = 10;
q = 4;    

Input_Number = 11;
Output_Number = 9;
    
figure_name = 'ErrorBar_Single_Total_Sobol.svg';


Single_Sobol = zeros(Input_Number,Output_Number,trial_num);
Total_Sobol = zeros(Input_Number,Output_Number,trial_num);

for trial=1:trial_num 

    % single sobol indices
    single_filename = strcat(input_dir,'Trial_',num2str(trial),'SingleSobol_SampleSize_100000Method_Latin.xlsx');

    single_T = readtable(single_filename,'ReadRowNames',true);
    single_T = table2array(single_T);
    
    Single_Sobol(:,:,trial) = single_T;

    % total sobol indices
    total_filename =  strcat(input_dir,'Trial_',num2str(trial),'TotalSobol_SampleSize_100000Method_Latin.xlsx');

    total_T = readtable(total_filename,'ReadRowNames',true);
    total_T = table2array(total_T);
    
    Total_Sobol(:,:,trial) = total_T;
end

%% plotting stuff   
inputVariable = {'Atrsc','Ptrsc','KdeA','AdeA','PdeA','KpolyA','ApolyA','PpolyA','Kdgrd','Adgrd','Pdgrd'};
outputQuantity = {'mean L/S','amp. of L/S','phase of L/S','mean L+S','amp. of L+S','phase of L+S','mean L','amp. of L','phase of L'};
outputQuantityFile = {'meanls','rals','Pls','meanss','rass','Pss','meanL','raL','PL'};

% box plot for single sobol indices
for i=1:length(outputQuantity)
    
    single_tmp = Single_Sobol(:,i,:);
    single_tmp = reshape(single_tmp,[Input_Number,trial_num]);
    single_tmp = single_tmp';
    
    
    total_tmp = Total_Sobol(:,i,:);
    total_tmp = reshape(total_tmp,[Input_Number,trial_num]);
    total_tmp = total_tmp';
    
    % change order of dgrd and deA
    single_tmp = single_tmp(:,[1 2 9 10 11 3 4 5 6 7 8]);
    total_tmp = total_tmp(:,[1 2 9 10 11 3 4 5 6 7 8]);
    
     
    Sobol = zeros(10,22);
    for k=1:11
        Sobol(:,2*k-1) = single_tmp(:,k);
        Sobol(:,2*k) = total_tmp(:,k);
    end    

    
    figure;
    x=1:1:Input_Number*2;
    y = mean(Sobol);
    err = std(Sobol);
    b1=errorbar(x,y,err,'LineWidth',2,'LineStyle','None','Marker','.','MarkerSize',20,'CapSize',14,...
        'MarkerFaceColor',[0,0,0],'MarkerEdgeColor',[0,0,0]);
    
    % cumstomize text on top of each bar
    for j=1:Input_Number*2
        if mod(j,2)==1
            text(j,y(j)+err(j),'S','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',20); % Single
        else
            text(j,y(j)+err(j),'T','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',20); % Total
        end
    end
    
    
    b1.Color = [0,0,0];
    
    ax = gca;
    ax.TickLength = [0.01,0.08];
    ax.YLim = [0 1];
    ax.LineWidth = 2;
    ax.FontSize = 36;
    ax.XLim  = [0 23];
    xticks(1:1:23)
    yticks(0:0.2:1)
    set(gca,'xticklabel',{[]});

    title(outputQuantity{i});
    saveas(gcf,strcat(output_dir,outputQuantityFile{i},figure_name));
end





