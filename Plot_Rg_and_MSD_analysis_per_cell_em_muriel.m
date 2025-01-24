%% Analyze and Plot Rg and MSD per cell
%% Abdullah R. Chaudhary

close all; clear all; clc; 
set(0,'DefaultFigureWindowStyle','docked')
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Statistical_testing/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/AC_codes_epmodified_20200603/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/plotSpread/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/violin/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/+NBP/');


plus_processive_runs=[];

colour_CTRL=[0.902 0.902 0.988];
colour_DCXKO=[0.008	0.737 0.745];

for k_choose = 1:2

if k_choose == 1    % CTRL
    col1=colour_CTRL;
    pth1='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Muriel xml trackmate files/Compiled Lysotracker/CTRL0/mats_2/';
    addpath(pth1);
    fl1='Rg_per_cell.mat';
    fl2='MSD_per_cell.mat';
    fl3='frac_proc_out_per_cell.mat';
    fl4='frac_proc_in_per_cell.mat';
    fl5='frac_diff_per_cell.mat';
    fl6='frac_stat_per_cell.mat';
    fl7='frac_diff_out_per_cell.mat';
    fl8='frac_diff_in_per_cell.mat';
    fl9='frac_proc_per_cell.mat';
    an=0;
    DT=0.5;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;
    msd_all=[];
    mean_rg_all=[];
    slp_all=[];
elseif k_choose == 2    % DCXKO
    col1=colour_DCXKO;
    pth1='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Muriel xml trackmate files/Compiled Lysotracker/DCXKO1/mats_2/';
    addpath(pth1);
    fl1='Rg_per_cell.mat';
    fl2='MSD_per_cell.mat';
    fl3='frac_proc_out_per_cell.mat';
    fl4='frac_proc_in_per_cell.mat';
    fl5='frac_diff_per_cell.mat';
    fl6='frac_stat_per_cell.mat';
    fl7='frac_diff_out_per_cell.mat';
    fl8='frac_diff_in_per_cell.mat';
    fl9='frac_proc_per_cell.mat';
    an=0;
    DT=0.5;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;
    msd_all=[];
    mean_rg_all=[];
    slp_all=[];
    
end

load(fullfile(pth1,fl1));
load(fullfile(pth1,fl2));
load(fullfile(pth1,fl3)); %Emily added
load(fullfile(pth1,fl4)); %Emily added 20220927
load(fullfile(pth1,fl5)); %Emily added 20220927
load(fullfile(pth1,fl6)); %Emily added 20220927
load(fullfile(pth1,fl7)); %Emily added 20221020
load(fullfile(pth1,fl8)); %Emily added 20221020
load(fullfile(pth1,fl9));

%Emily adding in a longer minimum run length version 20221024
proc_frac_t{k_choose}=t_proc_per_cell;
plus_proc_frac_t{k_choose}=t_proc_out_per_cell;
minus_proc_frac_t{k_choose}=t_proc_in_per_cell;
diff_frac_t{k_choose}=t_diff_per_cell;
plus_diff_frac_t{k_choose}=t_diff_out_per_cell; %Emily added 20221020
minus_diff_frac_t{k_choose}=t_diff_in_per_cell; %Emily added 20221020
stat_frac_t{k_choose}=t_stat_per_cell;
[timek2,msd2] = average_MSD_per_cell_emmodified4kchoose(res1,col1,k_choose);

% MSD Analysis:
 for ik=1:numel(res1)
    timek=res1{ik}.timek;
    msd=res1{ik}.msd;
    logt=res1{ik}.logtk;
    logmsd=res1{ik}.log_msd;
    alph=res1{ik}.slp;
    if isnan(msd) %Emily added because the msd had some cells without any data
        k_choose %1 for CTRL 2 for DCXKO, helps identify the empty cell
        ik
        warning='Empty cell in msd analysis'
    else
    % Figure(1)
    figure(k_choose*1e1), hold on, 
    p1=plot(timek,msd,'Color',col1);
    p1.Color(4)=0.25;
    xlabel('Time interval (sec)'); ylabel('MSD (\mum^2)'); 
    publication_fig(0,0,1);
    
    figure(3e1), hold on, 
    p1=plot(timek,msd,'Color',col1);
    p1.Color(4)=0.25;
    xlabel('Time interval (sec)'); ylabel('MSD (\mum^2)'); 
    publication_fig(0,0,1);
    
    % Figure(1)
    figure(k_choose*1e3), hold on, 
    p1=plot(logt,logmsd,'Color',col1);
    p1.Color(4)=0.25;
    xlabel('Log[Time interval (sec)]'); ylabel('Log[MSD (\mum^2)]'); 
    publication_fig(0,0,1);
    
    figure(3e3), hold on, 
    p1=plot(logt,logmsd,'Color',col1);
    p1.Color(4)=0.25;
    xlabel('Log[Time interval (sec)]'); ylabel('Log[MSD (\mum^2)]'); 
    publication_fig(0,0,1);
    
    msd_all=[msd_all,msd];
    slp_all=[slp_all,alph];
    
    rg_per_cell=rg_all{ik};
    if isempty(rg_per_cell) %Emily added because the Rg had some cells without any data
        k_choose %1 for CTRL 2 for DCXKO, helps identify the empty cell
        ik
        warning='Empty cell in Rg analysis'
        %clear rg_per_cell
    else
    [mean_rg,pci]=mle(rg_per_cell,'distribution','exp'); 
    mean_rg_all=[mean_rg_all,mean_rg];
    
    
    clear rg_per_cell
    end
    end
end

msd_data{k_choose}=msd_all';
rg_data{k_choose}=mean_rg_all';
timek_f{k_choose}=timek2;
MSD_f{k_choose}=msd2;
alph_f{k_choose}=slp_all; %slopes from the fits to the MSD curve

clear msd_all mean_rg_all rg_per_cell timek2 msd2 slp_all outward_runs_per_cell;
%rmpath(pth1)
%rmpath(pth2)
end

label_ctrl=repmat({'CTRL Plus'}, numel(alph_f{1}),1);
label_dcxko=repmat({'DCXKO Plus'}, numel(alph_f{2}),1);
label_both=[label_dcxko;label_ctrl];
alph_f_flipped=flip(alph_f);

figure('Name','Alpha_per_cell','NumberTitle','off'), hold on,
violin(alph_f_flipped, 'xlabel',{'DCXKO','CTRL'},'facecolor',[colour_DCXKO;colour_CTRL],'facealpha',0.3,'edgecolor','k','mc',[],'medc',[]);
boxplot([alph_f{2}; alph_f{1}], label_both, 'PlotStyle', 'compact', 'Colors', 'k', 'Symbol',' ')
% h1=notBoxPlot_nodotorline(alph_f{1},1,'style','line');
% set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_CTRL,'MarkerSize', 5);
% h2=notBoxPlot_nodotorline(alph_f{2},2,'style','line');
% set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_DCXKO,'MarkerSize', 5);
publication_fig(0,0,1)
axisHandle = gca; 
set(gca,'FontSize',24);
    set(gca, 'FontName', 'Arial');
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
    set(gca,'Box','on','LineWidth',2);
    set(gca, 'XTickLabel', {'DCXKO', 'CTRL'});
xticks([1 2]);
ylim([0 2.5]);
ylabel('\alpha');
% ylim([0.5 2.25]);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');
view([90 -90]);

rg_data_flipped=flip(rg_data);

figure('Name','Rg_per_cell','NumberTitle','off'), hold on,
violin(rg_data_flipped, 'xlabel',{'DCXKO','CTRL'},'facecolor',[colour_DCXKO;colour_CTRL],'facealpha',0.3,'edgecolor','k','mc',[],'medc',[]);
boxplot([rg_data{2}; rg_data{1}], label_both, 'PlotStyle', 'compact', 'Colors', 'k', 'Symbol',' ')
% h1=notBoxPlot_nodotorline(rg_data{1},1,'style','line');
% set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_CTRL, 'MarkerSize', 5);
% h2=notBoxPlot_nodotorline(rg_data{2},2,'style','line');
% set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_DCXKO,'MarkerSize', 5);
% xlim([0.5 4.5]);
% ylim([-0.25 2]);
publication_fig(0,0,1)
axisHandle = gca; 
set(gca,'FontSize',24);
    set(gca, 'FontName', 'Arial');
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
    set(gca,'Box','on','LineWidth',2);
    set(gca, 'XTickLabel', {'DCXKO', 'CTRL'});
xticks([1 2]);
ylim([0 3]);
ylabel('Radius of Gyration (\mum)');
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off')
view([90 -90]);

avg_plus_proc_frac_ctrl=mean(plus_proc_frac_t{1});
avg_plus_proc_frac_DCXKO=mean(plus_proc_frac_t{2});

plus_proc_frac_t_flipped=flip(plus_proc_frac_t);

figure('Name','Directional_bias','NumberTitle','off'), hold on,
violin(plus_proc_frac_t_flipped, 'xlabel',{'DCXKO','CTRL'},'facecolor',[colour_DCXKO;colour_CTRL],'facealpha',0.3,'edgecolor','k','mc',[],'medc',[]);
boxplot([plus_proc_frac_t{2}; plus_proc_frac_t{1}], label_both, 'PlotStyle', 'compact', 'Colors', 'k', 'Symbol',' ')
% h1=notBoxPlot_nodotorline(plus_proc_frac_t{1},1,'style','line');
% set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_CTRL, 'MarkerSize', 5);
% h2=notBoxPlot_nodotorline(plus_proc_frac_t{2},2,'style','line');
% set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_DCXKO,'MarkerSize', 5);
% ylabel('Fraction of time of processive runs moving outwards');
axisHandle = gca; 
set(gca,'FontSize',24);
    set(gca, 'FontName', 'Arial');
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
    set(gca,'Box','on','LineWidth',2);
    set(gca, 'XTickLabel', {'DCXKO', 'CTRL'});
xticks([1 2]);
ylim([0 1]);
ylabel('Directional Bias');
publication_fig(0,0,1)
view([90 -90]);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');

plus_diff_frac_t_flipped=flip(plus_diff_frac_t);

%Emily added for testing if min_lr was set correctly 20221020
figure('Name','Diff_runs_out','NumberTitle','off'), hold on,
violin(plus_diff_frac_t_flipped, 'xlabel',{'DCXKO','CTRL'},'facecolor',[colour_DCXKO;colour_CTRL],'facealpha',0.3,'edgecolor','k','mc',[],'medc',[]);
boxplot([plus_diff_frac_t{2}; plus_diff_frac_t{1}], label_both, 'PlotStyle', 'compact', 'Colors', 'k', 'Symbol',' ')
% h5=notBoxPlot_nodotorline(plus_diff_frac_t{1},1,'style','line');
% set(h5.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_CTRL, 'MarkerSize', 5);
% h6=notBoxPlot_nodotorline(plus_diff_frac_t{2},2,'style','line');
% set(h6.data,'MarkerFaceColor','none','MarkerEdgeColor',colour_DCXKO,'MarkerSize', 5);
axisHandle = gca; 
set(gca,'FontSize',24);
    set(gca, 'FontName', 'Arial');
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
    set(gca,'Box','on','LineWidth',2);
    set(gca, 'XTickLabel', {'DCXKO', 'CTRL'});
xticks([1 2]);
ylim([0 1]);
ylabel('Fraction of time of diffusive runs moving outwards');
publication_fig(0,0,1);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off');
view([90 -90]);

%Calculating median to plot bar plots at longer minimum run lengths
plus_runs_CTRL=mean(plus_proc_frac_t{1});
plus_runs_DCXKO=mean(plus_proc_frac_t{2});
minus_runs_CTRL=mean(minus_proc_frac_t{1});
minus_runs_DCXKO=mean(minus_proc_frac_t{2});
proc_runs_CTRL=mean(proc_frac_t{1});
proc_runs_DCXKO=mean(proc_frac_t{2});
diff_frac_t_CTRL=mean(diff_frac_t{1});
diff_frac_t_DCXKO=mean(diff_frac_t{2});
stat_frac_t_CTRL=mean(stat_frac_t{1});
stat_frac_t_DCXKO=mean(stat_frac_t{2});

%SEM for the different types of runs
error_proc_CTRL=std(proc_frac_t{1})/(sqrt(length(proc_frac_t{1})));
error_diff_CTRL=std(diff_frac_t{1})/(sqrt(length(diff_frac_t{1})));
error_stat_CTRL=std(stat_frac_t{1})/(sqrt(length(stat_frac_t{1})));
error_proc_DCXKO=std(proc_frac_t{2})/(sqrt(length(proc_frac_t{2})));
error_diff_DCXKO=std(diff_frac_t{2})/(sqrt(length(diff_frac_t{2})));
error_stat_DCXKO=std(stat_frac_t{2})/(sqrt(length(stat_frac_t{2})));

%creating a variable with all the errors for all the conditions
all_err_bars=[error_proc_CTRL error_proc_DCXKO ; error_diff_CTRL error_diff_DCXKO ; error_stat_CTRL error_stat_DCXKO ];

values=[proc_runs_CTRL proc_runs_DCXKO ; diff_frac_t_CTRL diff_frac_t_DCXKO; stat_frac_t_CTRL stat_frac_t_DCXKO];
figure('Name','Bar plot proc stat diff ','NumberTitle','off'), hold on,
b2=bar(values ,'grouped');
b2(1).FaceColor=colour_CTRL;
b2(2).FaceColor=colour_DCXKO;
b2(1).FaceAlpha=0.5;
b2(2).FaceAlpha=0.5;
hold on
[ngroups,nbars]= size(values );
groupwidth= min(0.8, nbars/(nbars + 1.5));
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
   % Calculate center of each bar
   x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
   errorbar(x, values (:,i), all_err_bars (:,i), 'k', 'linestyle', 'none','LineWidth',2);
end
scatter(0.8571-0.05+0.1*rand(numel(proc_frac_t{1}),1),proc_frac_t{1},35,'MarkerFaceColor','none','MarkerEdgeColor',colour_CTRL*0.9,'LineWidth',2);
scatter(1.1429-0.05+0.1*rand(numel(proc_frac_t{2}),1),proc_frac_t{2},35,'MarkerFaceColor','none','MarkerEdgeColor',colour_DCXKO*0.9,'LineWidth',2);
scatter(1.8571-0.05+0.1*rand(numel(diff_frac_t{1}),1),diff_frac_t{1},35,'MarkerFaceColor','none','MarkerEdgeColor',colour_CTRL*0.9,'LineWidth',2);
scatter(2.1429-0.05+0.1*rand(numel(diff_frac_t{2}),1),diff_frac_t{2},35,'MarkerFaceColor','none','MarkerEdgeColor',colour_DCXKO*0.9,'LineWidth',2);
scatter(2.8571-0.05+0.1*rand(numel(stat_frac_t{1}),1),stat_frac_t{1},35,'MarkerFaceColor','none','MarkerEdgeColor',colour_CTRL*0.9,'LineWidth',2);
scatter(3.1429-0.05+0.1*rand(numel(stat_frac_t{2}),1),stat_frac_t{2},35,'MarkerFaceColor','none','MarkerEdgeColor',colour_DCXKO*0.9,'LineWidth',2);
xticks([1 2 3]);
xticklabels({'Processive', 'Diffusive', 'Stationary'});
ylabel('Fraction of runs');
publication_fig(0,0,1)
box('off');

%Bootstrapping Rg and Alpha values

% bootstrapped_all=[];
% 
% nsamp=min([numel(rg_all{1}),numel(rg_all{2})])
% 
% for dataset=1:2
% bstrp_CTRL_lyso=Loic_bootstrap_code_04092019_em(rg_all{dataset},nsamp,1000,5);
% bootstrapped_dataset{dataset}=bstrp_CTRL_lyso';
% end
% 
% 
% difference_CTRL_DCXKO = bootstrapped_dataset{1}.bstrap_means - bootstrapped_dataset{2}.bstrap_means;
% difference_fake = bootstrapped_dataset{1}.bstrap_means - bootstrapped_dataset{1}.bstrap_means; %testing the bootstrapping algorithm
% 
% figure('Name','Rg_bootstrapping_histogram','NumberTitle','off'), hold on,
% histogram(difference_CTRL_DCXKO,'facecolor',colour_CTRL, 'edgecolor','none','facealpha',0.5), hold on,
% xlabel('Radius of Gyration (\mum)');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);
% 
% % figure('Name','Rg_bootstrapping_histogram_fake','NumberTitle','off'), hold on,
% % histogram(difference_fake,'facecolor',[0 0 1], 'edgecolor','none','facealpha',0.5), hold on,
% % xlabel('Radius of Gyration (\mum)');
% % ylabel('Number of Trajectories');
% % publication_fig(0,0,1);
% 
% figure('Name','Rg_bootstrapping_histogram_CTRL_DCXKO','NumberTitle','off'), hold on,
% histogram(difference_CTRL_DCXKO,'facecolor',colour_CTRL, 'edgecolor','none','facealpha',0.5), hold on,
% xlabel('Radius of Gyration (\mum)');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);

% nsamp=min([numel(alph_f{1}),numel(alph_f{2})])
% 
% for dataset=1:2
% bstrp_CTRL_lyso_alpha=Loic_bootstrap_code_04092019_em(alph_f{dataset},nsamp,1000,5);
% bootstrapped_dataset_alpha{dataset}=bstrp_CTRL_lyso_alpha';
% end
% 
% 
% difference_CTRL_DCXKO_alpha = bootstrapped_dataset_alpha{1}.bstrap_means - bootstrapped_dataset_alpha{2}.bstrap_means;
% difference_fake_alpha = bootstrapped_dataset_alpha{1}.bstrap_means - bootstrapped_dataset_alpha{1}.bstrap_means; %testing the bootstrapping
% 
% figure('Name','Bootstrapping_histogram_alpha','NumberTitle','off'), hold on,
% histogram(difference_CTRL_DCXKO_alpha,'facecolor',colour_CTRL, 'edgecolor','none','facealpha',0.5), hold on,
% xlabel('\alpha');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);
% 
% figure('Name','Bootstrapping_fake_alpha','NumberTitle','off'), hold on,
% histogram(difference_fake_alpha,'facecolor',[0 0 1], 'edgecolor','none','facealpha',0.5), hold on,
% xlabel('\alpha');
% ylabel('Number of Trajectories');
% xlim([-1 1]);
% ylim([0 2000]);
% publication_fig(0,0,1);
% 
% 
% figure('Name','Bootstrapping_histogram_CTRL_DCXKO_alpha','NumberTitle','off'), hold on,
% histogram(difference_CTRL_DCXKO_alpha,'facecolor',colour_CTRL, 'edgecolor','none','facealpha',0.5), hold on,
% xlabel('\alpha');
% ylabel('Number of Trajectories');
% publication_fig(0,0,1);

%T testing
[CTRL_DCXKO_rg,p_CTRL_DCXKO_rg]=ttest2(rg_data{1}, rg_data{2})
[CTRL_DCXKO_alpha,p_CTRL_DCXKO_alpha]=ttest2(alph_f{1}, alph_f{2})

%testing for minimum run length 1um
[CTRL_DCXKO_dirbias ,p_CTRL_DCXKO_dirbias]=ranksum(plus_proc_frac_t{1}, plus_proc_frac_t{2}) % changed based on distribution shape
[CTRL_DCXKO_proc,p_CTRL_DCXKO_proc]=ttest2(proc_frac_t{1}, proc_frac_t{2})
[CTRL_DCXKO_diff,p_CTRL_DCXKO_diff]=ttest2(diff_frac_t{1}, diff_frac_t{2})
[CTRL_DCXKO_stat,p_CTRL_DCXKO_stat]=ttest2(stat_frac_t{1}, stat_frac_t{2})

titles={'CTRL vs DCXKO'};
rg_stats=p_CTRL_DCXKO_rg;
alpha_stats=p_CTRL_DCXKO_alpha;
dirbias_stats=p_CTRL_DCXKO_dirbias;
proc_stats=p_CTRL_DCXKO_proc;
diff_stats=p_CTRL_DCXKO_diff;
stat_stats=p_CTRL_DCXKO_stat;

%Auto-saving stats and figures
% 
fileprefix='20250123_lyso_DCXKO_per_cell_';

tempdir_1 = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Muriel xml trackmate files/';   % Your destination folder
FolderName_1 = tempdir_1;   % Your destination folder
%Saving stats to a table
stats_table=table(titles, rg_stats, alpha_stats,dirbias_stats, proc_stats, diff_stats, stat_stats);
writetable(stats_table,fullfile(FolderName_1, [fileprefix,'.csv']),'WriteRowNames',true);

tempdir = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Muriel xml trackmate files/';  % Your destination folder
FolderName = tempdir;  % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
 FigHandle = FigList(iFig);
 FigName  = get(FigHandle, 'Name');
 savefig(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.fig'])); 
%  saveas(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.png']));
 saveas(FigHandle, fullfile(FolderName, [append(fileprefix,FigName), '.svg']));
end