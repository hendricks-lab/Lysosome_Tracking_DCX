%% Analyze velocties and run length of processive events:
% with each new dataset you run in this code, change the directories. If you changed the exposure time between
% experiments, make sure you also change DT (after each directory is defined). 
%To have specific filenames for the figures generated and automatically saved, change the variable q at the end of the script. 

%Emily modifying Sept 17/22 to fix some aesthetics

clear all; close all; clc; 
cd(strtok(userpath, pathsep));
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/wavelet_deriv/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/General_codes/General-codes/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/raacampbell-notBoxPlot-7d90c27/code/+NBP/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/violin/'); %Emily added 20220808
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/Statistical_testing/'); %Emily added 20220809

set(0,'defaulttextinterpreter','tex');
set(0,'DefaultFigureWindowStyle','docked');
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

for k_choose = 1:2
    
if k_choose == 1 % Control
    cd('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Muriel xml trackmate files/Compiled Lysotracker/CTRL0/mats_2/');
    DT=0.5;
    rev_rate=[];
    proc_run=[];
    proc_time=[];
    diff_run=[];
    diff_time=[];
    event_dur_v2=[];
    avg_lr_per_traj1p=[];
    avg_lr_per_traj1m=[];
    avg_vel_per_traj1p=[];
    avg_vel_per_traj1m=[];
    pe=[];
    dc=[];
    all_pos=[];
    delt=[];
    msd_pts=[];
    kcol=[0.5 0.5 0.5];
    alpha=[];
    
elseif k_choose == 2 %DCX KO
    cd('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Muriel xml trackmate files/Compiled Lysotracker/DCXKO1/mats_2/');
    DT=0.5;
    rev_rate=[];
    proc_run=[];
    proc_time=[];
    diff_run=[];
    diff_time=[];
    event_dur_v2=[];
    avg_lr_per_traj1p=[];
    avg_lr_per_traj1m=[];
    avg_vel_per_traj1p=[];
    avg_vel_per_traj1m=[];
    pe=[];
    dc=[];
    all_pos=[];
    delt=[];
    msd_pts=[];
    kcol=[0 0 1];
    alpha=[];
end
min_lr=1.0;

load('position');

for i=1:numel(position)
    pos1=position{i}; 
    time=[1:numel(pos1)]'.*DT;
    diff_coeff=Diffusion_coefficient(i);
    event_dur=time(end)-time(1);
    pos=pos1;%smooth(pos1,10,'sgolay',1);
    res=analyze_run_length_reversals_v4(time,pos,min_lr,0,k_choose);
    
    rev_rate=[rev_rate,(numel(res.proc_run)./(event_dur))];%[rev_rate,res.reversal_rate];
    proc_run=[proc_run;res.proc_run];
    proc_time=[proc_time;res.proc_time];
    diff_run=[diff_run;res.diff_run];
    diff_time=[diff_time;res.diff_time];
    event_dur_v2=[event_dur_v2,event_dur];
    
    avg_lr_per_traj1p=[avg_lr_per_traj1p;res.avg_lrplus];
    avg_lr_per_traj1m=[avg_lr_per_traj1m;res.avg_lrminus];
    avg_vel_per_traj1p=[avg_vel_per_traj1p;res.avg_velplus];
    avg_vel_per_traj1m=[avg_vel_per_traj1m;res.avg_velminus];
    
    %pe1=func_peclet_number(res.run_bw_rev,res.time_bw_rev,Diffusion_coefficient(i));
    %pe=[pe,pe1];
%     dc=[dc,diff_coeff];
    all_pos=[all_pos;pos];
    %Emily recommented the section below because she got weird errors and I
    %think we ended up using the per cell alpha values anyway. 
%     [deltat, msdpts, sem, log_deltat, log_msdpts, alpha_1, DiffCoef] = MSD_2d(xk1{i}, yk1{i}, DT, k_choose); %Emily uncommented 20220809, fits first 75 frames with time averaged MSD
%     delt=[delt,deltat]; %Emily uncommented 20220809 
%     msd_pts=[msd_pts,msdpts]; %Emily uncommented 20220809
%     alpha=[alpha,alpha_1]; %Emily added 20220815

    clear res pe1 deltat msdpts sem log_deltat log_msdpts alpha_1 Diff_Coef
end

display(i)

proc_runs{k_choose}=proc_run;
proc_times{k_choose}=proc_time;
rev_rates{k_choose}=rev_rate;
diff_runs{k_choose}=diff_run;
diff_times{k_choose}=diff_time;
avg_lr_per_trajp{k_choose}=avg_lr_per_traj1p;
avg_vel_per_trajp{k_choose}=avg_vel_per_traj1p;
avg_lr_per_trajm{k_choose}=-1.*avg_lr_per_traj1m;
avg_vel_per_trajm{k_choose}=-1.*avg_vel_per_traj1m;
%peclet{k_choose}=pe;
diff_coeff_fn{k_choose}=dc;
pos_all{k_choose}=all_pos;
all_delt{k_choose}=delt;
all_msd{k_choose}=msd_pts;
all_alpha{k_choose}=alpha;
% xk_all{k_choose}=xk1;
% yk_all{k_choose}=yk1;
posit{k_choose}=position;
duration_traj{k_choose}=event_dur_v2;

clear position xk_traj yk_traj
clear pos time 
proc_run=[];
proc_time=[];
rev_rate=[];
diff_run=[];
diff_times=[];
pe=[];
dc=[];
all_pos=[];
event_dur_v2=[];
avg_lr_per_traj1p=[];
avg_lr_per_traj1m=[];
avg_vel_per_traj1p=[];
avg_vel_per_traj1m=[];
end

%% Analyze Reversal Rate:
% xk_rev=-1:0.01:1;
% rev_rate_ctrl=rev_rates{1};
% rev_rate_s421d=rev_rates{2};
% hist_xkrev_ctrl=hist(rev_rate_ctrl,xk_rev);
% hist_xkrev_s421d=hist(rev_rate_s421d,xk_rev);
% figure(1), hold on, 
% bar(xk_rev,hist_xkrev_s421d./sum(hist_xkrev_s421d),'BarWidth',1,'FaceColor','g','EdgeColor','g');
% hold on, 
% stairs(xk_rev,hist_xkrev_ctrl./sum(hist_xkrev_ctrl),'Color','k','LineWidth',2);
% xlabel('Reversal Rate (# of reversals/s)'); ylabel('Fraction of Events'); 
% publication_fig(0,0,1);
% % ./sum(hist_xkrev_ctrl+hist_xkrev_s421d)

%% Analyze Processive Runs:
proc_runs_ctrl=proc_runs{1};
proc_times_ctrl=proc_times{1};
proc_runs_DCXKO=proc_runs{2};
proc_times_DCXKO=proc_times{2};


proc_plus_run_ctrl=proc_runs_ctrl(find(proc_runs_ctrl>0));
proc_plus_time_ctrl=proc_times_ctrl(find(proc_runs_ctrl>0));
proc_minus_run_ctrl=proc_runs_ctrl(find(proc_runs_ctrl<0));
proc_minus_time_ctrl=proc_times_ctrl(find(proc_runs_ctrl<0));

proc_plus_run_DCXKO=proc_runs_DCXKO(find(proc_runs_DCXKO>0));
proc_plus_time_DCXKO=proc_times_DCXKO(find(proc_runs_DCXKO>0));
proc_minus_run_DCXKO=proc_runs_DCXKO(find(proc_runs_DCXKO<0));
proc_minus_time_DCXKO=proc_times_DCXKO(find(proc_runs_DCXKO<0));

total_time_ctrl=sum(proc_plus_time_ctrl)+sum(proc_minus_time_ctrl);
total_time_DCXKO=sum(proc_plus_time_DCXKO)+sum(proc_minus_time_DCXKO);

fraction_plus_ctrl=sum(proc_plus_time_ctrl)./total_time_ctrl;
fraction_minus_ctrl=sum(proc_minus_time_ctrl)./total_time_ctrl;
fraction_plus_DCXKO=sum(proc_plus_time_DCXKO)./total_time_DCXKO;
fraction_minus_DCXKO=sum(proc_minus_time_DCXKO)./total_time_DCXKO;

sem_total_time_ctrl=std(total_time_ctrl)/sqrt(length(total_time_ctrl));
sem_total_time_DCXKO=std(total_time_DCXKO)/sqrt(length(total_time_DCXKO));
%Calculating the standard error of the mean (SEM) of each processive
%motility event time
sem_proc_plus_time_ctrl=std(proc_plus_time_ctrl)/sqrt(length(proc_plus_time_ctrl));
sem_proc_minus_time_ctrl=std(proc_minus_time_ctrl)/sqrt(length(proc_minus_time_ctrl));
sem_proc_plus_time_DCXKO=std(proc_plus_time_DCXKO)/sqrt(length(proc_plus_time_DCXKO));
sem_proc_minus_time_DCXKO=std(proc_minus_time_DCXKO)/sqrt(length(proc_minus_time_DCXKO));
%Error propagation for dividing two variables to get the SEM for the
%fraction
sem_fraction_plus_ctrl=sem_total_time_ctrl/mean(total_time_ctrl)+sem_proc_plus_time_ctrl/mean(proc_plus_time_ctrl); %Error propagation for dividing two variables
sem_fraction_minus_ctrl=sem_total_time_ctrl/mean(total_time_ctrl)+sem_proc_minus_time_ctrl/mean(proc_minus_time_ctrl); 
sem_fraction_plus_DCXKO=sem_total_time_DCXKO/mean(total_time_DCXKO)+sem_proc_plus_time_DCXKO/mean(proc_plus_time_DCXKO); 
sem_fraction_minus_DCXKO=sem_total_time_DCXKO/mean(total_time_DCXKO)+sem_proc_minus_time_DCXKO/mean(proc_minus_time_DCXKO); 


figure('Name','Fraction of Processive Runs','NumberTitle','off'), hold on,
bar(1.25,fraction_plus_ctrl,'BarWidth',0.2,'FaceColor',[0.90 0.9 0.98],'EdgeColor','k','Facealpha',1);
hold on, 
er1 = errorbar(1.25,fraction_plus_ctrl,sem_fraction_plus_ctrl,sem_fraction_plus_ctrl);    
er1.Color = [0 0 0];                            
er1.LineStyle = 'none';  
bar(1.25,-1.*fraction_minus_ctrl,'BarWidth',0.2,'FaceColor',[0.9 0.9 0.98],'EdgeColor','k','Facealpha',0.3);
er2 = errorbar(1.25,-1*fraction_minus_ctrl,-1*sem_fraction_minus_ctrl,-1*sem_fraction_minus_ctrl);    
er2.Color = [0 0 0];                            
er2.LineStyle = 'none'; 
hold on, 
bar(1.5,fraction_plus_DCXKO,'BarWidth',0.2,'FaceColor',[0.008 0.737 0.745],'EdgeColor',[0 0 0],'Facealpha',1);
er3 = errorbar(1.5,fraction_plus_DCXKO,sem_fraction_plus_DCXKO,sem_fraction_plus_DCXKO);    
er3.Color = [0 0 0];                            
er3.LineStyle = 'none'; 
hold on, 
bar(1.5,-1.*fraction_minus_DCXKO,'BarWidth',0.2,'FaceColor', [0.008	0.737	0.745],'EdgeColor',[0 0 0],'Facealpha',0.3);
er4 = errorbar(1.5,-1*fraction_minus_DCXKO,-1*sem_fraction_minus_DCXKO,-1*sem_fraction_minus_DCXKO);    
er4.Color = [0 0 0];                            
er4.LineStyle = 'none'; 
hold on, 
xlim([1 2]);
set(gca,'ygrid','on');
publication_fig(0,0,1);

%Emily making a violin plot of the plus/minus fractions
all_fractions={fraction_plus_ctrl,fraction_minus_ctrl,fraction_plus_DCXKO,fraction_minus_DCXKO};
all_fractions_flipped=flip(all_fractions);
plus_fractions={fraction_plus_ctrl,fraction_plus_DCXKO};
plus_fractions_flipped=flip(plus_fractions);
minus_fractions={fraction_minus_ctrl,fraction_minus_DCXKO};
minus_fractions_flipped=flip(minus_fractions);

% label_p_frac_ctrl=repmat({'CTRL+'}, numel(fraction_plus_ctrl),1);
% label_p_frac_dcxko=repmat({'DCXKO+'}, numel(fraction_minus_ctrl),1);
% label_m_frac_ctrl=repmat({'CTRL-'}, numel(fraction_plus_DCXKO),1);
% label_m_frac_dcxko=repmat({'DCXKO-'}, numel(fraction_minus_DCXKO),1);
% label_all_fracs=[label_m_frac_dcxko; label_p_frac_dcxko; label_m_frac_ctrl; label_p_frac_ctrl];

figure('Name','Fraction of Processive Runs','NumberTitle','off'), hold on,
violin(all_fractions_flipped,'xlabel',{'DCXKO -','DCXKO +','CTRL -','CTRL +'},'facecolor',[0.008 0.737 0.745; 0.008 0.737 0.745; 0.9	0.9	0.98; 0.90	0.9	0.98],'facealpha',1,'edgecolor','k','medc','k--');
% boxplot([all_fractions_flipped{1}; all_fractions_flipped{2}; all_fractions_flipped{3}; all_fractions_flipped{4}], 'PlotStyle', 'compact', 'Colors', 'k', 'Symbol',' ')
ylabel('Fraction of Processive Runs');
view([90 -90]);
publication_fig(0,0,1);

% label_p_fracs=[label_p_frac_dcxko; label_p_frac_ctrl];

figure('Name','Fraction of Plus-Ended Processive Runs','NumberTitle','off'), hold on,
violin(plus_fractions_flipped,'xlabel',{'DCXKO ','CTRL'},'facecolor',[ 0.008 0.737 0.745; 0.90	0.9	0.98],'facealpha',1,'edgecolor','k','medc','k--');
ylabel('Fraction of Processive Runs');
view([90 -90]);
publication_fig(0,0,1);

% label_m_fracs=[label_m_frac_dcxko; label_m_frac_dcxko];

figure('Name','Fraction of Minus-Ended Processive Runs','NumberTitle','off'), hold on,
violin(minus_fractions_flipped,'xlabel',{'DCXKO ','CTRL'},'facecolor',[ 0.008 0.737 0.745; 0.9	0.9	0.98],'facealpha',1,'edgecolor','k','medc','k--');
ylabel('Fraction of Processive Runs');
view([90 -90]);
publication_fig(0,0,1);

%% Plot CDF Plots:

figure('Name','CDF plus','NumberTitle','off'), hold on,
h1=cdfplot(proc_plus_run_ctrl);
h1.Color=[0.902 0.902 0.988];
hold on, 
h2=cdfplot(proc_plus_run_DCXKO);
h2.Color=[0.008 0.737 0.745];
xlim([0 10]);
xlabel('LOG[L_{REV,+}]'); ylabel('CDF'); 
title('Plus-ended runs');
publication_fig(0,0,1);

%figure(4), hold on, 
figure('Name','CDF minus','NumberTitle','off'), hold on,
h1=cdfplot(-1.*proc_minus_run_ctrl);
h1.Color=[0.902 0.902 0.988];
hold on, 
h2=cdfplot(-1.*proc_minus_run_DCXKO);
h2.Color=[0.008 0.737 0.745];
hold on, 
xlim([0 10]);
xlabel('LOG[L_{REV,-}]'); ylabel('CDF'); 
title('Minus-ended runs');
publication_fig(0,0,1);

%% Determine Dispersion by measuring variance:
pos_ctrl=pos_all{1};
pos_DCXKO=pos_all{2};

plus_pos_ctrl=pos_ctrl(find(pos_ctrl>0));
minus_pos_ctrl=pos_ctrl(find(pos_ctrl<0));

plus_pos_DCXKO=pos_DCXKO(find(pos_DCXKO>0));
minus_pos_DCXKO=pos_DCXKO(find(pos_DCXKO<0));


plus_pos_dist=0:0.05:11;
minus_pos_dist=-11:0.05:0;

hist_plus_ctrl=hist(plus_pos_ctrl,plus_pos_dist);
hist_plus_DCXKO=hist(plus_pos_DCXKO,plus_pos_dist);
hist_minus_ctrl=hist(minus_pos_ctrl,minus_pos_dist);
hist_minus_DCXKO=hist(minus_pos_DCXKO,minus_pos_dist);

%figure(5), hold on, 
figure('Name','Position frequency plus','NumberTitle','off'), hold on,
stairs(plus_pos_dist,hist_plus_ctrl./sum(hist_plus_ctrl+hist_plus_DCXKO),'Color',[0.902 0.902 0.988],'LineWidth',2);
hold on, 
stairs(plus_pos_dist,hist_plus_DCXKO./sum(hist_plus_ctrl+hist_plus_DCXKO),'Color',[0.008 0.737 0.745],'LineWidth',2);
hold on, 
xlabel('Position (Plus)'); ylabel('Frequency'); 
publication_fig(0,0,1);

%figure(6), hold on, 
figure('Name','Position frequency minus','NumberTitle','off'), hold on,
stairs(minus_pos_dist,hist_minus_ctrl./sum(hist_minus_ctrl+hist_minus_DCXKO),'Color',[0.902 0.902 0.988],'LineWidth',2);
hold on, 
stairs(minus_pos_dist,hist_minus_DCXKO./sum(hist_minus_ctrl+hist_minus_DCXKO),'Color',[0.008 0.737 0.745],'LineWidth',2);
set(gca, 'YScale', 'log');
xlabel('Position (Minus)'); ylabel('Frequency'); 
publication_fig(0,0,1);

%figure(1000), hold on, 
figure('Name','Dispersion','NumberTitle','off'), hold on,
bar(1,std(plus_pos_ctrl).*1000,'BarWidth',0.5,'FaceColor',[0.902 0.902 0.988],'Facealpha',1);
hold on, 
errorbar(1,std(plus_pos_ctrl).*1000,0.032.*1000,'k-','LineWidth',2)
hold on,
bar(2,std(plus_pos_DCXKO).*1000,'BarWidth',0.5,'FaceColor',[0.008 0.737 0.745],'Facealpha',1);
hold on, 
errorbar(2,std(plus_pos_DCXKO).*1000,0.032.*1000,'k-','LineWidth',2)
hold on,
bar(1,-1*std(minus_pos_ctrl).*1000,'BarWidth',0.5,'FaceColor',[0.902 0.902 0.988],'Facealpha',0.3);
hold on, 
errorbar(1,-1*std(minus_pos_ctrl).*1000,0.057.*1000,'k-','LineWidth',2)
hold on,
bar(2,-1*std(minus_pos_DCXKO).*1000,'BarWidth',0.5,'FaceColor',[0.008 0.737 0.745],'Facealpha',0.3);
hold on, 
errorbar(2,-1*std(minus_pos_DCXKO).*1000,0.057.*1000,'k-','LineWidth',2)
ylabel('Average Dispersion (nm)');
publication_fig(0,0,1);
xlim([0 5])
ylim([-1500 1200])
set(gca,'Ygrid','on');

%% Velocities:
plus_vel_ctrl=proc_plus_run_ctrl./proc_plus_time_ctrl;
plus_vel_DCXKO=proc_plus_run_DCXKO./proc_plus_time_DCXKO;
minus_vel_ctrl=-1.*proc_minus_run_ctrl./proc_minus_time_ctrl;
minus_vel_DCXKO=-1.*proc_minus_run_DCXKO./proc_minus_time_DCXKO;

xplus_vel=0:0.05:2;
xminus_vel=-2:0.05:0;

xhist_plus_vel_ctrl=hist(plus_vel_ctrl,xplus_vel); 
xhist_minus_vel_ctrl=hist(minus_vel_ctrl,xminus_vel);
xhist_plus_vel_DCXKO=hist(plus_vel_DCXKO,xplus_vel);
xhist_minus_vel_DCXKO=hist(minus_vel_DCXKO,xminus_vel);


%figure(7), hold on, 
figure('Name','Processive velocity distribution','NumberTitle','off'), hold on,
plot(xplus_vel,xhist_plus_vel_ctrl./sum(xhist_plus_vel_ctrl+xhist_plus_vel_DCXKO),'Color',[0.902 0.902 0.988]);
hold on, 
plot(xplus_vel,xhist_plus_vel_DCXKO./sum(xhist_plus_vel_ctrl+xhist_plus_vel_DCXKO),'Color',[0.008 0.737 0.745]);
hold on, 
plot(xminus_vel,xhist_minus_vel_ctrl./sum(xhist_plus_vel_ctrl+xhist_plus_vel_DCXKO),'Color',[0.902 0.902 0.988]);
hold on, 
plot(xminus_vel,xhist_minus_vel_DCXKO./sum(xhist_plus_vel_ctrl+xhist_plus_vel_DCXKO), 'Color',[0.008 0.737 0.745]);
xlabel('Velocity (\mum/s)'); ylabel('Number of Processive Runs'); 
xlim([-2.2 2.2]);
publication_fig(0,0,1);

%% Travel Distance vs Velocity (um/s)

%figure(8), hold on, subplot(1,2,1)
figure('Name','Run length velocity plus','NumberTitle','off'), hold on,% subplot(1,2,1)
scatter(proc_plus_run_ctrl,plus_vel_ctrl,'o','MarkerFaceColor',[0.902 0.902 0.988],'MarkerEdgeColor',[0.902 0.902 0.988], 'MarkerFaceAlpha',0.3);
hold on,
scatter(proc_plus_run_DCXKO,plus_vel_DCXKO,'o','MarkerFaceColor',[0.008 0.737 0.745],'MarkerEdgeColor',[0.008 0.737 0.745],'MarkerFaceAlpha',0.3);
xlabel('Travel Distance, L_R (\mum)'); ylabel('Velocity (\mum/s)');
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca,'LineWidth',2);
set(gca,'FontSize',24);
set(gca, 'FontName', 'Arial');
set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
set(gca,'Box','on');
xlim([0 8]); ylim([0 2]);
    

%figure(9), hold on, subplot(1,2,1)
figure('Name','Run length velocity minus','NumberTitle','off'), hold on,% subplot(1,2,1)
scatter(-1.*proc_minus_run_ctrl,minus_vel_ctrl,'o','MarkerFaceColor',[0.902 0.902 0.988],'MarkerEdgeColor',[0.902 0.902 0.988],'MarkerFaceAlpha',0.3);
hold on,
scatter(-1.*proc_minus_run_DCXKO,minus_vel_DCXKO,'o','MarkerFaceColor',[0.008 0.737 0.745],'MarkerEdgeColor',[0.008 0.737 0.745],'MarkerFaceAlpha',0.3);
xlabel('Travel Distance, L_R (\mum)'); ylabel('Velocity (\mum/s)');
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca,'LineWidth',2);
set(gca,'FontSize',24);
set(gca, 'FontName', 'Arial');
set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
set(gca,'Box','on');
xlim([0 8]); ylim([0 2]);
    

proc_rev_rate_ct=rev_rates{1}(find(rev_rates{1}>0));
proc_rev_rate_tr=rev_rates{2}(find(rev_rates{2}>0));

%Emily changes start here, Aug 5/22
all_run_lengths={proc_plus_run_ctrl,-1.*proc_minus_run_ctrl,proc_plus_run_DCXKO,-1.*proc_minus_run_DCXKO};
all_run_lengths_flipped={-1.*proc_minus_run_DCXKO,proc_plus_run_DCXKO,-1.*proc_minus_run_ctrl,proc_plus_run_ctrl};
plus_run_lengths={proc_plus_run_ctrl,proc_plus_run_DCXKO};
minus_run_lengths={-1.*proc_minus_run_ctrl,-1.*proc_minus_run_DCXKO};
plus_run_lengths_flipped=flip(plus_run_lengths);
minus_run_lengths_flipped=flip(minus_run_lengths);

label_p_rl_ctrl=repmat({'CTRL Plus'}, numel(proc_plus_run_ctrl),1);
label_p_rl_dcxko=repmat({'DCXKO Plus'}, numel(proc_plus_run_DCXKO),1);
label_m_rl_ctrl=repmat({'CTRL Minus'}, numel(proc_minus_run_ctrl),1);
label_m_rl_dcxko=repmat({'DCXKO Minus'}, numel(proc_minus_run_DCXKO),1);
label_all_rls=[label_m_rl_dcxko; label_p_rl_dcxko; label_m_rl_ctrl; label_p_rl_ctrl];

mean_p_rl_ctrl=mean(proc_plus_run_ctrl);
mean_p_rl_dcxko=mean(proc_plus_run_DCXKO);
mean_m_rl_ctrl=mean(proc_minus_run_ctrl);
mean_m_rl_dcxko=mean(proc_minus_run_DCXKO);

mean_p_vel_ctrl=mean(plus_vel_ctrl);
mean_p_vel_dcxko=mean(plus_vel_DCXKO);
mean_m_vel_ctrl=mean(minus_vel_ctrl);
mean_m_vel_dcxko=mean(minus_vel_DCXKO);

figure('Name','Run Length','NumberTitle','off'), hold on,
violin(all_run_lengths_flipped,'xlabel',{'DCXKO -','DCXKO +','CTRL -','CTRL +'},'facecolor',[0.008 0.737 0.745; 0.008 0.737 0.745; 0.9	0.9	0.98; 0.90	0.9	0.98],'facealpha',1,'edgecolor','k','mc',[],'medc',[]);
boxplot([all_run_lengths_flipped{1};all_run_lengths_flipped{2};all_run_lengths_flipped{3};all_run_lengths_flipped{4}], label_all_rls, 'PlotStyle', 'compact', 'Colors', 'k', 'Symbol',' ')
set(gca,'FontSize',24);
    set(gca, 'FontName', 'Arial');
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
    set(gca,'Box','on','LineWidth',2);
    set(gca, 'XTickLabel', {'DCXKO-','DCXKO+', 'CTRL-', 'CTRL+'});
xticks([1 2 3 4]);
ylabel('Run Length (\mum)');
view([90 -90]);
publication_fig(0,0,1);

label_p_rls=[label_p_rl_dcxko; label_p_rl_ctrl];

figure('Name','Plus-Ended Run Lengths','NumberTitle','off'), hold on,
violin(plus_run_lengths_flipped,'xlabel',{'DCXKO ','CTRL'},'facecolor',[0.008 0.737 0.745; 0.90 0.9 0.98],'facealpha',1,'edgecolor','k','mc',[],'medc',[]);
boxplot([plus_run_lengths_flipped{1};plus_run_lengths_flipped{2}], label_p_rls, 'PlotStyle', 'compact', 'Colors', 'k', 'Symbol',' ');
set(gca,'FontSize',24);
    set(gca, 'FontName', 'Arial');
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
    set(gca,'Box','on','LineWidth',2);
    set(gca, 'XTickLabel', {'DCXKO','CTRL'});
xticks([1 2]);
ylabel('Plus-Ended Run Lengths (\mum)');
view([90 -90]);
publication_fig(0,0,1);

label_m_rls=[label_m_rl_dcxko; label_m_rl_ctrl];

figure('Name','Minus-Ended Run Lengths','NumberTitle','off'), hold on,
violin(minus_run_lengths_flipped,'xlabel',{'DCXKO ','CTRL'},'facecolor',[0.008 0.737 0.745; 0.90	0.9	0.98],'facealpha',1,'edgecolor','k','mc',[],'medc',[]);
boxplot([minus_run_lengths_flipped{1};minus_run_lengths_flipped{2}], label_m_rls,'PlotStyle', 'compact', 'Colors', 'k', 'Symbol',' ');
set(gca,'FontSize',24);
    set(gca, 'FontName', 'Arial');
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
    set(gca,'Box','on','LineWidth',2);
    set(gca, 'XTickLabel', {'DCXKO','CTRL'});
xticks([1 2]);
ylabel('Minus-Ended Run Lengths (\mum)');
view([90 -90]);
publication_fig(0,0,1);

plus_run= min_lr:0.1:3;
minus_run=min_lr:0.1:3;
hist_plus_run_ctrl=hist(proc_plus_run_ctrl,numel(plus_run));
hist_minus_run_ctrl=hist(-1.*proc_minus_run_ctrl,numel(minus_run));
hist_plus_run_DCXKO=hist(proc_plus_run_DCXKO,numel(plus_run));
hist_minus_run_DCXKO=hist(-1.*proc_minus_run_DCXKO,numel(minus_run));

fraction_anterograde_runs_vs_proc_ctrl=hist_plus_run_ctrl./(numel(proc_plus_run_ctrl));
fraction_anterograde_runs_vs_proc_DCXKO=hist_plus_run_DCXKO./(numel(proc_plus_run_DCXKO));

fraction_retrograde_runs_vs_proc_ctrl=hist_minus_run_ctrl./(numel(proc_minus_run_ctrl));
fraction_retrograde_runs_vs_proc_DCXKO=hist_minus_run_DCXKO./(numel(proc_minus_run_DCXKO));

%Exponential fit
exp_fit_plus_ctrl=fit(plus_run.', (hist_plus_run_ctrl./(numel(proc_plus_run_ctrl))).', 'exp1');
exp_fit_minus_ctrl=fit(minus_run.', (hist_minus_run_ctrl./(numel(proc_minus_run_ctrl))).', 'exp1');
exp_fit_plus_DCXKO=fit(plus_run.', (hist_plus_run_DCXKO./(numel(proc_plus_run_DCXKO))).', 'exp1');
exp_fit_minus_DCXKO=fit(minus_run.', (hist_minus_run_DCXKO./(numel(proc_minus_run_DCXKO))).', 'exp1');

%Plotting the fraction of processive runs vs run length with the
%exponential fits (Emily)
figure('Name','Fraction of Processive Runs vs Run Length Plus', 'Numbertitle', 'off'), hold on
data1_p=plot(plus_run, hist_plus_run_ctrl./(numel(proc_plus_run_ctrl)),'.','MarkerFaceColor',[0.90	0.9	0.98],'MarkerEdgeColor',[0.90	0.9	0.98],'MarkerSize',10);
data2_p=plot(plus_run, hist_plus_run_DCXKO./(numel(proc_plus_run_DCXKO)),'.','MarkerFaceColor',[0.008 0.737 0.745],'MarkerEdgeColor',[0.008 0.737 0.745],'MarkerSize',10);
fit1_p=plot(exp_fit_plus_ctrl,plus_run, hist_plus_run_ctrl./(numel(proc_plus_run_ctrl)));
fit2_p=plot(exp_fit_plus_DCXKO,plus_run, hist_plus_run_DCXKO./(numel(proc_plus_run_DCXKO)));
set(fit1_p,'color',[0.9 0.9 0.98]);
set(fit2_p,'color',[0.008 0.737 0.745]);
xlabel('Run Length (\mum)');
ylabel('Fraction of Processive Runs');
set(gca,'LineWidth',2);
set(gca,'FontSize',24);
set(gca, 'FontName', 'Arial');
set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off')

figure('Name','Fraction of Processive Runs vs Run Length Minus', 'Numbertitle', 'off'), hold on
data1_m=plot(minus_run, hist_minus_run_ctrl./(numel(proc_minus_run_ctrl)),'.','MarkerFaceColor',[0.90	0.9	0.98],'MarkerEdgeColor',[0.90	0.9	0.98],'MarkerSize',10);
data2_m=plot(minus_run, hist_minus_run_DCXKO./(numel(proc_minus_run_DCXKO)),'.','MarkerFaceColor',[0.008 0.737 0.745],'MarkerEdgeColor',[0.008 0.737 0.745],'MarkerSize',10);
fit1_m=plot(exp_fit_minus_ctrl,minus_run, hist_minus_run_ctrl./(numel(proc_minus_run_ctrl)));
fit2_m=plot(exp_fit_minus_DCXKO,minus_run, hist_minus_run_DCXKO./(numel(proc_minus_run_DCXKO)));
set(fit1_m,'color',[0.9 0.9 0.98]);
set(fit2_m,'color',[0.008 0.737 0.745]);
xlabel('Run Length (\mum)');
ylabel('Number of Processive Runs');
set(gca,'LineWidth',2);
set(gca,'FontSize',24);
set(gca, 'FontName', 'Arial');
set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
hFig=findall(0,'type','figure');
hLeg=findobj(hFig(1,1),'type','legend');
set(hLeg,'visible','off')

%Setting up variables for boostrapping run length fraction
run_length_fraction_ctrl_p=[hist_plus_run_ctrl;plus_run].';
run_length_fraction_DCXKO_p=[hist_plus_run_DCXKO;plus_run].';
run_length_fraction_ctrl_m=[hist_minus_run_ctrl;abs(minus_run)].';
run_length_fraction_DCXKO_m=[hist_minus_run_DCXKO;abs(minus_run)].';

% Emily plotting average velocity 20220808

all_velocities={avg_vel_per_trajp{1},avg_vel_per_trajm{1}, avg_vel_per_trajp{2}, avg_vel_per_trajm{2}};
all_velocities_flipped=flip(all_velocities);
plus_velocities_flipped=flip(avg_vel_per_trajp);
minus_velocities_flipped=flip(avg_vel_per_trajm);

label_p_vel_ctrl=repmat({'CTRL Plus'}, numel(avg_vel_per_trajp{1}),1);
label_p_vel_dcxko=repmat({'DCXKO Plus'}, numel(avg_vel_per_trajp{2}),1);
label_m_vel_ctrl=repmat({'CTRL Minus'}, numel(avg_vel_per_trajm{1}),1);
label_m_vel_dcxko=repmat({'DCXKO Minus'}, numel(avg_vel_per_trajm{2}),1);
label_all_vels=[label_m_vel_dcxko; label_p_vel_dcxko; label_m_vel_ctrl; label_p_vel_ctrl];

figure('Name','Average Velocity','NumberTitle','off'), hold on,
violin(all_velocities_flipped,'xlabel',{'DCXKO -','DCXKO +','CTRL -','CTRL +'},'facecolor',[0.008 0.737 0.745; 0.008 0.737 0.745; 0.9	0.9	0.98; 0.90	0.9	0.98],'facealpha',1,'edgecolor','k','mc',[],'medc',[]);
boxplot([all_velocities_flipped{1}; all_velocities_flipped{2}; all_velocities_flipped{3}; all_velocities_flipped{4}], label_all_vels, 'PlotStyle', 'compact', 'Colors', 'k', 'Symbol',' ')
% h1=notBoxPlot_nodotorline(all_velocities{1},4,'style','line'); %need to put ctrl in position 2 so that it's on top when the view is flipped
% set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',[0.6	0.6	0.65], 'MarkerSize', 5);
% h2=notBoxPlot_nodotorline(all_velocities{2},3,'style','line'); %need to put ctrl in position 2 so that it's on top when the view is flipped
% set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',[0.6	0.6	0.65], 'MarkerSize', 5);
% h3=notBoxPlot_nodotorline(all_velocities{3},2,'style','line'); %need to put ctrl in position 2 so that it's on top when the view is flipped
% set(h3.data,'MarkerFaceColor','none','MarkerEdgeColor',[0.10	0.55	0.55], 'MarkerSize', 5);
% h4=notBoxPlot_nodotorline(all_velocities{4},1,'style','line'); 
% set(h4.data,'MarkerFaceColor','none','MarkerEdgeColor',[0.10	0.55	0.55],'MarkerSize', 5);
ylim([0 7])
axisHandle = gca; 
set(gca,'FontSize',24);
    set(gca, 'FontName', 'Arial');
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
    set(gca,'Box','on','LineWidth',2);
    set(gca, 'XTickLabel', {'DCXKO-','DCXKO+', 'CTRL-', 'CTRL+'});
xticks([1 2 3 4]);
ylabel('Average Velocity (\mum/s)');
view([90 -90]);
publication_fig(0,0,1);

plus_proc_vels_ctrl=avg_vel_per_trajp{1};
plus_proc_vels_dcxko=avg_vel_per_trajp{2};
plus_proc_vels_all=[plus_proc_vels_dcxko;plus_proc_vels_ctrl];

label_p_vel_all=[label_p_vel_dcxko; label_p_vel_ctrl];

figure('Name','Average Velocity Plus','NumberTitle','off'), hold on,
violin(plus_velocities_flipped,'xlabel',{'DCXKO ','CTRL'},'facecolor',[ 0.008 0.737 0.745; 0.90	0.9	0.98],'facealpha',1,'edgecolor','k','mc',[],'medc',[]);
boxplot(plus_proc_vels_all, label_p_vel_all,'PlotStyle', 'compact', 'Colors', 'k', 'Symbol',' ')
% h1=notBoxPlot_nodotorline(avg_vel_per_trajp{1},2,'style','line'); %need to put ctrl in position 2 so that it's on top when the view is flipped
% set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',[0.6	0.6	0.65], 'MarkerSize', 5);
% h2=notBoxPlot_nodotorline(avg_vel_per_trajp{2},1,'style','line'); %need to put ctrl in position 2 so that it's on top when the view is flipped
% set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',[0.10	0.55	0.55], 'MarkerSize', 5);
axisHandle = gca; 
set(gca,'FontSize',24);
    set(gca, 'FontName', 'Arial');
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
    set(gca,'Box','on','LineWidth',2);
    set(gca, 'XTickLabel', {'DCXKO', 'CTRL'});
xticks([1 2]);
ylim([0 5]);
ylabel('Average Velocity Plus (\mum/s)');
view([90 -90]);

minus_proc_vels_ctrl=avg_vel_per_trajm{1};
minus_proc_vels_dcxko=avg_vel_per_trajm{2};
minus_proc_vels_all=[minus_proc_vels_dcxko;minus_proc_vels_ctrl];

label_m_vel_all=[label_m_vel_dcxko; label_m_vel_ctrl];

figure('Name','Average Velocity Minus','NumberTitle','off'), hold on,
violin(minus_velocities_flipped,'xlabel',{'DCXKO ','CTRL'},'facecolor',[ 0.008 0.737 0.745; 0.90	0.9	0.98],'facealpha',1,'edgecolor','k','mc',[],'medc',[]);
boxplot(minus_proc_vels_all,label_m_vel_all, 'PlotStyle', 'compact', 'Colors', 'k', 'Symbol',' ')
% h1=notBoxPlot_nodotorline(avg_vel_per_trajm{1},2,'style','line'); %need to put ctrl in position 2 so that it's on top when the view is flipped
% set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',[0.6	0.6	0.65], 'MarkerSize', 5);
% h2=notBoxPlot_nodotorline(avg_vel_per_trajm{2},1,'style','line'); %need to put ctrl in position 2 so that it's on top when the view is flipped
% set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',[0.10	0.55	0.55], 'MarkerSize', 5);
% ylabel('Average Velocity Minus (\mum/s)');
axisHandle = gca; 
set(gca,'FontSize',24);
    set(gca, 'FontName', 'Arial');
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
    set(gca,'Box','on','LineWidth',2);
    set(gca, 'XTickLabel', {'DCXKO', 'CTRL'});
ylabel('Average Velocity Minus (\mum/s)');
ylim([0 7]);
xticks([1 2]);
view([90 -90]);

%Emily plotting Alpha values
% ctrl_alpha=all_alpha{1}(~isnan(all_alpha{1}));
% DCXKO_alpha=all_alpha{2}(~isnan(all_alpha{2}));
% all_alphas_flipped={DCXKO_alpha,ctrl_alpha};
% 
% figure('Name','Alpha','NumberTitle','off'), hold on,
% violin(all_alphas_flipped,'xlabel',{'DCXKO','CTRL'},'facecolor',[ 0.008 0.737 0.745; 0.90	0.9	0.98],'facealpha',1,'edgecolor','k','medc','k--');
% h1=notBoxPlot_nodotorline(DCXKO_alpha,1,'style','line'); %need to put ctrl in position 2 so that it's on top when the view is flipped
% set(h1.data,'MarkerFaceColor','none','MarkerEdgeColor',[0.10	0.55	0.55], 'MarkerSize', 5);
% h2=notBoxPlot_nodotorline(ctrl_alpha,2,'style','line'); %need to put ctrl in position 2 so that it's on top when the view is flipped
% set(h2.data,'MarkerFaceColor','none','MarkerEdgeColor',[0.6	0.6	0.65], 'MarkerSize', 5);
% ylabel('\alpha');
% view([90 -90]);
% publication_fig(0,0,1);

%% Bootstrapping for statistical analysis
% Number of samples calculation for all bootstrapped samples
nsamp_f_p=min([numel(fraction_plus_ctrl);numel(fraction_plus_DCXKO)]);
nsamp_f_m=min([numel(fraction_minus_ctrl);numel(fraction_minus_DCXKO)]);
nsamp_rl=min([numel(proc_plus_run_ctrl);numel(proc_plus_run_DCXKO);numel(proc_minus_run_ctrl);numel(proc_minus_run_DCXKO)]);
nsamp_vel=min([numel(avg_vel_per_trajp{1});numel(avg_vel_per_trajp{2});numel(avg_vel_per_trajm{1});numel(avg_vel_per_trajm{2})]);
% nsamp_alpha=min([numel(ctrl_alpha);numel(DCXKO_alpha)]);
% Bootstrap calculations for all bootstrapped samples
%Fraction
ctrl_bstrp_f_p=Loic_bootstrap_code_04092019_em_single_variable(fraction_plus_ctrl.', nsamp_f_p, 1000,0.05);
DCXKO_bstrp_f_p=Loic_bootstrap_code_04092019_em_single_variable(fraction_plus_DCXKO.', nsamp_f_p,1000,0.05);
ctrl_bstrp_f_m=Loic_bootstrap_code_04092019_em_single_variable(fraction_minus_ctrl.', nsamp_f_m,1000,0.05);
DCXKO_bstrp_f_m=Loic_bootstrap_code_04092019_em_single_variable(fraction_minus_DCXKO.', nsamp_f_m,1000,0.05);
%Run Length Fraction
ctrl_bstrp_rl_frac_p=Loic_bootstrap_code_04092019_OT_em(run_length_fraction_ctrl_p, 1000,0.05,3,5);
DCXKO_bstrp_rl_frac_p=Loic_bootstrap_code_04092019_OT_em(run_length_fraction_DCXKO_p, 1000,0.05,3,5);
ctrl_bstrp_rl_frac_m=Loic_bootstrap_code_04092019_OT_em(run_length_fraction_ctrl_m, 1000,0.05,3,5);
DCXKO_bstrp_rl_frac_m=Loic_bootstrap_code_04092019_OT_em(run_length_fraction_DCXKO_m, 1000,0.05,3,5);
%Run Length
ctrl_bstrp_rl_p=Loic_bootstrap_code_04092019_em_single_variable(proc_plus_run_ctrl.',nsamp_rl,1000,5);
DCXKO_bstrp_rl_p=Loic_bootstrap_code_04092019_em_single_variable(proc_plus_run_DCXKO.',nsamp_rl,1000,5);
ctrl_bstrp_rl_m=Loic_bootstrap_code_04092019_em_single_variable(proc_minus_run_ctrl.',nsamp_rl,1000,5);
DCXKO_bstrp_rl_m=Loic_bootstrap_code_04092019_em_single_variable(proc_minus_run_DCXKO.',nsamp_rl,1000,5);
%Velocity
ctrl_bstrp_vel_p=Loic_bootstrap_code_04092019_em_single_variable(avg_vel_per_trajp{1}(~isnan(avg_vel_per_trajp{1})).',nsamp_vel,1000,5);
DCXKO_bstrp_vel_p=Loic_bootstrap_code_04092019_em_single_variable(avg_vel_per_trajp{2}(~isnan(avg_vel_per_trajp{2})).',nsamp_vel,1000,5);
ctrl_bstrp_vel_m=Loic_bootstrap_code_04092019_em_single_variable(avg_vel_per_trajm{1}(~isnan(avg_vel_per_trajm{1})).',nsamp_vel,1000,5);
DCXKO_bstrp_vel_m=Loic_bootstrap_code_04092019_em_single_variable(avg_vel_per_trajm{2}(~isnan(avg_vel_per_trajm{2})).',nsamp_vel,1000,5);
%Alpha
% ctrl_bstrp_alpha=Loic_bootstrap_code_04092019_em_single_variable(ctrl_alpha.',nsamp_alpha,1000,5);
% DCXKO_bstrp_alpha=Loic_bootstrap_code_04092019_em_single_variable(DCXKO_alpha.',nsamp_alpha,1000,5);

% Calculating bootstrap means difference for all bootstrapped samples
%Fraction Plus or Minus
diff_f_DCXKO_CTRL_p=DCXKO_bstrp_f_p.bstrap_means-ctrl_bstrp_f_p.bstrap_means;
diff_f_CTRL_p_m=ctrl_bstrp_f_p.bstrap_means-abs(ctrl_bstrp_f_m.bstrap_means);
diff_f_DCXKO_CTRL_m=abs(DCXKO_bstrp_f_m.bstrap_means)-abs(ctrl_bstrp_f_m.bstrap_means);
diff_f_DCXKO_p_m=DCXKO_bstrp_f_p.bstrap_means-abs(DCXKO_bstrp_f_m.bstrap_means);
%Run Length Fraction
diff_rlf_DCXKO_CTRL_p=cell2mat(DCXKO_bstrp_rl_frac_p.bstrap_means)-cell2mat(ctrl_bstrp_rl_frac_p.bstrap_means);
diff_rlf_CTRL_p_m=cell2mat(ctrl_bstrp_rl_frac_p.bstrap_means)-abs(cell2mat(ctrl_bstrp_rl_frac_m.bstrap_means));
diff_rlf_DCXKO_CTRL_m=abs(cell2mat(DCXKO_bstrp_rl_frac_m.bstrap_means))-abs(cell2mat(ctrl_bstrp_rl_frac_m.bstrap_means));
diff_rlf_DCXKO_p_m=cell2mat(DCXKO_bstrp_rl_frac_p.bstrap_means)-abs(cell2mat(DCXKO_bstrp_rl_frac_m.bstrap_means));
%Run Length
diff_rl_DCXKO_CTRL_p=DCXKO_bstrp_rl_p.bstrap_means-ctrl_bstrp_rl_p.bstrap_means;
diff_rl_CTRL_p_m=ctrl_bstrp_rl_p.bstrap_means-abs(ctrl_bstrp_rl_m.bstrap_means);
diff_rl_DCXKO_CTRL_m=abs(DCXKO_bstrp_rl_m.bstrap_means)-abs(ctrl_bstrp_rl_m.bstrap_means);
diff_rl_DCXKO_p_m=DCXKO_bstrp_rl_p.bstrap_means-abs(DCXKO_bstrp_rl_m.bstrap_means);
%Velocity
diff_vel_DCXKO_CTRL_p=DCXKO_bstrp_vel_p.bstrap_means-ctrl_bstrp_vel_p.bstrap_means;
diff_vel_CTRL_p_m=ctrl_bstrp_vel_p.bstrap_means-abs(ctrl_bstrp_vel_m.bstrap_means);
diff_vel_DCXKO_CTRL_m=abs(DCXKO_bstrp_vel_m.bstrap_means)-abs(ctrl_bstrp_vel_m.bstrap_means);
diff_vel_DCXKO_p_m=DCXKO_bstrp_vel_p.bstrap_means-abs(DCXKO_bstrp_vel_m.bstrap_means);
%Alpha
% diff_alpha_DCXKO_CTRL=DCXKO_bstrp_alpha.bstrap_means-ctrl_bstrp_alpha.bstrap_means;

% Plotting Bootstrap Difference histograms for all bootstrapped samples
%Fraction Plus/Minus
% figure, histogram(diff_f_DCXKO_CTRL_p);
% line([0 0],get(gca, 'Ylim'),'color','r', 'LineWidth', 2);
% xlabel('DCXKO-CTRL Plus Fraction'), ylabel('Frequency');
% publication_fig(0,0,1);
% figure, histogram(diff_f_CTRL_p_m);
% line([0 0],get(gca, 'Ylim'),'color','r', 'LineWidth', 2);
% xlabel('CTRL Plus-CTRL Minus Fraction'), ylabel('Frequency');
% publication_fig(0,0,1);
% figure, histogram(diff_f_DCXKO_CTRL_m);
% line([0 0],get(gca, 'Ylim'),'color','r', 'LineWidth', 2);
% xlabel('DCXKO-CTRL Minus Fraction'), ylabel('Frequency');
% publication_fig(0,0,1);
% figure, histogram(diff_f_DCXKO_p_m);
% line([0 0],get(gca, 'Ylim'),'color','r', 'LineWidth', 2);
% xlabel('DCXKO Plus-DCXKO Minus Fraction'), ylabel('Frequency');
% publication_fig(0,0,1);
%Run Length Fraction
figure, histogram(diff_rlf_DCXKO_CTRL_p);
line([0 0],get(gca, 'Ylim'),'color','r', 'LineWidth', 2);
xlabel('DCXKO-CTRL Run Length Fraction Plus'), ylabel('Frequency');
publication_fig(0,0,1);
figure, histogram(diff_rlf_CTRL_p_m);
line([0 0],get(gca, 'Ylim'),'color','r', 'LineWidth', 2);
xlabel('CTRL Plus-CTRL Minus Run Length Fraction'), ylabel('Frequency');
publication_fig(0,0,1);
figure, histogram(diff_rlf_DCXKO_CTRL_m);
line([0 0],get(gca, 'Ylim'),'color','r', 'LineWidth', 2);
xlabel('DCXKO-CTRL Run Length Fraction Minus'), ylabel('Frequency');
publication_fig(0,0,1);
figure, histogram(diff_rlf_DCXKO_p_m);
line([0 0],get(gca, 'Ylim'),'color','r', 'LineWidth', 2);
xlabel('DCXKO Plus-DCXKO Minus Run Length Fraction'), ylabel('Frequency');
publication_fig(0,0,1);
%Run Length
figure, histogram(diff_rl_DCXKO_CTRL_p);
line([0 0],get(gca, 'Ylim'),'color','r', 'LineWidth', 2);
xlabel('DCXKO-CTRL Plus Run Length'), ylabel('Frequency');
publication_fig(0,0,1);
figure, histogram(diff_rl_CTRL_p_m);
line([0 0],get(gca, 'Ylim'),'color','r', 'LineWidth', 2);
xlabel('CTRL Plus-CTRL Minus Run Length'), ylabel('Frequency');
publication_fig(0,0,1);
figure, histogram(diff_rl_DCXKO_CTRL_m);
line([0 0],get(gca, 'Ylim'),'color','r', 'LineWidth', 2);
xlabel('DCXKO-CTRL Minus Run Length'), ylabel('Frequency');
publication_fig(0,0,1);
figure, histogram(diff_rl_DCXKO_p_m);
line([0 0],get(gca, 'Ylim'),'color','r', 'LineWidth', 2);
xlabel('DCXKO Plus-DCXKO Minus Run Length'), ylabel('Frequency');
publication_fig(0,0,1);
%Velocity
figure, histogram(diff_vel_DCXKO_CTRL_p);
line([0 0],get(gca, 'Ylim'),'color','r', 'LineWidth', 2);
xlabel('DCXKO-CTRL Plus'), ylabel('Frequency');
publication_fig(0,0,1);
figure, histogram(diff_vel_CTRL_p_m);
line([0 0],get(gca, 'Ylim'),'color','r', 'LineWidth', 2);
xlabel('CTRL Plus-CTRL Minus'), ylabel('Frequency');
publication_fig(0,0,1);
figure, histogram(diff_vel_DCXKO_CTRL_m);
line([0 0],get(gca, 'Ylim'),'color','r', 'LineWidth', 2);
xlabel('DCXKO-CTRL Minus'), ylabel('Frequency');
publication_fig(0,0,1);
figure, histogram(diff_vel_DCXKO_p_m);
line([0 0],get(gca, 'Ylim'),'color','r', 'LineWidth', 2);
xlabel('DCXKO Plus-DCXKO Minus'), ylabel('Frequency');
publication_fig(0,0,1);
% Alpha
% figure, histogram(diff_alpha_DCXKO_CTRL);
% line([0 0],get(gca, 'Ylim'),'color','r', 'LineWidth', 2);
% xlabel('DCXKO-CTRL Alpha'), ylabel('Frequency');
% publication_fig(0,0,1);
% Confidence interval evaluation for all bootstrapped samples
%Fraction
ci_bstrp_f_ctrl_p=quantile(ctrl_bstrp_f_p.bstrap_means, [0.05 0.95]);
ci_bstrp_f_DCXKO_p=quantile(DCXKO_bstrp_f_p.bstrap_means, [0.05 0.95]);
ci_bstrp_f_ctrl_m=quantile(ctrl_bstrp_f_m.bstrap_means, [0.05 0.95]);
ci_bstrp_f_DCXKO_m=quantile(DCXKO_bstrp_f_m.bstrap_means, [0.05 0.95]);
%Run Length Fraction
ci_bstrp_ctrl_rlf_p=quantile(cell2mat(ctrl_bstrp_rl_frac_p.bstrap_means), [0.05 0.95]);
ci_bstrp_DCXKO_rlf_p=quantile(cell2mat(DCXKO_bstrp_rl_frac_p.bstrap_means), [0.05 0.95]);
ci_bstrp_ctrl_rlf_m=quantile(cell2mat(ctrl_bstrp_rl_frac_m.bstrap_means), [0.05 0.95]);
ci_bstrp_DCXKO_rlf_m=quantile(cell2mat(DCXKO_bstrp_rl_frac_p.bstrap_means), [0.05 0.95]);
%Run Length
ci_bstrp_rl_ctrl_p=quantile(ctrl_bstrp_rl_p.bstrap_means, [0.05 0.95]);
ci_bstrp_rl_DCXKO_p=quantile(DCXKO_bstrp_rl_p.bstrap_means, [0.05 0.95]);
ci_bstrp_rl_ctrl_m=quantile(ctrl_bstrp_rl_m.bstrap_means, [0.05 0.95]);
ci_bstrp_rl_DCXKO_m=quantile(DCXKO_bstrp_rl_m.bstrap_means, [0.05 0.95]);
%Velocity
ci_bstrp_vel_ctrl_p=quantile(ctrl_bstrp_vel_p.bstrap_means, [0.05 0.95]);
ci_bstrp_vel_DCXKO_p=quantile(DCXKO_bstrp_vel_p.bstrap_means, [0.05 0.95]);
ci_bstrp_vel_ctrl_m=quantile(ctrl_bstrp_vel_m.bstrap_means, [0.05 0.95]);
ci_bstrp_vel_DCXKO_m=quantile(DCXKO_bstrp_vel_m.bstrap_means, [0.05 0.95]);
%Alpha
% ci_alpha_bstrp_ctrl_p=quantile(ctrl_bstrp_alpha.bstrap_means, [0.05 0.95]);
% ci_alpha_bstrp_DCXKO_p=quantile(DCXKO_bstrp_alpha.bstrap_means, [0.05 0.95]);

% P value evaluation for all bootstrapped samples
%Fraction (right now this doesn't make sense because it's a single
%value)-working on it!
p_f_ctrl_DCXKO_p=numel(find(diff_f_DCXKO_CTRL_p>0))/numel(diff_f_DCXKO_CTRL_p)
p_f_ctrl_p_m=numel(find(diff_f_CTRL_p_m>0))/numel(diff_f_CTRL_p_m)
p_f_ctrl_DCXKO_m=numel(find(diff_f_DCXKO_CTRL_m>0))/numel(diff_f_DCXKO_CTRL_m)
p_f_DCXKO_p_m=numel(find(diff_f_DCXKO_p_m>0))/numel(diff_f_DCXKO_p_m)
%Run Length Fraction
p_ctrl_DCXKO_rlf_p=numel(find(diff_rlf_DCXKO_CTRL_p>0))/numel(diff_rlf_DCXKO_CTRL_p)
p_ctrl_rlf_p_m=numel(find(diff_rlf_CTRL_p_m>0))/numel(diff_rlf_CTRL_p_m)
p_ctrl_DCXKO_rlf_m=numel(find(diff_rlf_DCXKO_CTRL_m>0))/numel(diff_rlf_DCXKO_CTRL_m)
p_DCXKO_rlf_p_m=numel(find(diff_rlf_DCXKO_p_m>0))/numel(diff_rlf_DCXKO_p_m)
%Run Length
p_rl_ctrl_DCXKO_p=numel(find(diff_rl_DCXKO_CTRL_p>0))/numel(diff_rl_DCXKO_CTRL_p)
p_rl_ctrl_p_m=numel(find(diff_rl_CTRL_p_m>0))/numel(diff_rl_CTRL_p_m)
p_rl_ctrl_DCXKO_m=numel(find(diff_rl_DCXKO_CTRL_m>0))/numel(diff_rl_DCXKO_CTRL_m)
p_rl_DCXKO_p_m=numel(find(diff_rl_DCXKO_p_m>0))/numel(diff_rl_DCXKO_p_m)
%Velocity
p_vel_ctrl_DCXKO_p=numel(find(diff_vel_DCXKO_CTRL_p>0))/numel(diff_vel_DCXKO_CTRL_p)
p_vel_ctrl_p_m=numel(find(diff_vel_CTRL_p_m>0))/numel(diff_vel_CTRL_p_m)
p_vel_ctrl_DCXKO_m=numel(find(diff_vel_DCXKO_CTRL_m>0))/numel(diff_vel_DCXKO_CTRL_m)
p_vel_DCXKO_p_m=numel(find(diff_vel_DCXKO_p_m>0))/numel(diff_vel_DCXKO_p_m)
%Alpha
% p_alpha_ctrl_DCXKO=numel(find(diff_alpha_DCXKO_CTRL>0))/numel(diff_alpha_DCXKO_CTRL)

%% Saving all the Figures

filesavename='20250123_Muriel';

tempdir = '/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Muriel xml trackmate files/Compiled Lysotracker/';   % Your destination folder
FolderName = tempdir;   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  savefig(FigHandle, fullfile(FolderName, [append(filesavename,FigName), '.fig'])); 
  saveas(FigHandle, fullfile(FolderName, [append(filesavename,FigName), '.svg']));
end


