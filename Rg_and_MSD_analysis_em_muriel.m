% Final analysis:

close all; clear all; clc; 
set(0,'DefaultFigureWindowStyle','docked')
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/AC_codes_epmodified_20200603/');

min_lr=0.4;
min_disp=0;

%% Savitsky Golay smoothing
span=10;    % 10 is good, 25 is good as well
pwr=1;      % 1 is good, 2 is good

for k_choose = 1:2

if k_choose == 1    % Control
    kcol='k';%[0.5 0.5 0.5];
    cd('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Muriel xml trackmate files/Compiled Lysotracker/CTRL0/mats_2/');
    fl='*Copy*.mat';
    an=0;
    DT=0.5; %delta t, the exposure time-modify this based on your experiment
    save_dir='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Muriel xml trackmate files/Compiled Lysotracker/CTRL0/mats_2/';
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    rg=[];
    alpha=[];
    Diff_coeff=[];
    fa=1;
    
elseif k_choose == 2    % DCX knockout
    kcol = 'b';
    cd('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Muriel xml trackmate files/Compiled Lysotracker/DCXKO1/mats_2/');
    fl='*Copy*.mat';
    an=0.030;
    DT=0.5; %delta t, the exposure time-modify this based on your experiment
    save_dir='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Muriel xml trackmate files/Compiled Lysotracker/DCXKO1/mats_2/';
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    rg=[];
    alpha=[];
    Diff_coeff=[];
    fa=0.2;
end

fls=dir(fl);
kpr=0;
    
for k=1:numel(fls)
    load(fls(k).name);
    display(fls(k).name);
    for j=1:numel(ab)
        %display(j)
    if (numel(ab(j).position)>15) %&& (mean(ab(j).medint>3000))
    kf=kf+1;
    pos{kf}=ab(j).position;
    tim{kf}=[1:numel(ab(j).position)]'.*DT;%ab(j).tk;
    xk_smooth=smooth(ab(j).xk,span,'sgolay',pwr);
    yk_smooth=smooth(ab(j).yk,span,'sgolay',pwr);
    %res{j}=analyze_run_length_reversals_v2(tim{kf},pos{kf},min_lr);
    %proc_runs=[proc_runs;res{j}.proc_run];
    position_smooth=ab(j).position;%smooth(ab(j).position,span,'sgolay',pwr);
    
    %% Radius of Gyration:
    rg_1=func_Rg_Linda_v2(ab(j).xk,ab(j).yk);
    rg=[rg,rg_1];
    
    %% Mean Squared Displacement (2D)
    [deltat, msdpts, sem, log_deltat, log_msdpts, alpha_1, DiffCoef] = MSD_2d (xk_smooth, yk_smooth, DT, k_choose);
    
%     if iscolumn(res{j}.proc_time)==0
%     res{j}.proc_time=res{j}.proc_time';
%     proc_times=[proc_times;res{j}.proc_time];
%     else
%     proc_times=[proc_times;res{j}.proc_time];
%     end
%     
    if alpha_1>0
    alpha=[alpha,alpha_1];
    delT{kf}=deltat;
    MSD{kf}=msdpts;
    xk1{kf}=xk_smooth;
    yk1{kf}=yk_smooth;
    dat_pos_tim{kf}=position_smooth;
    Diff_coeff=[Diff_coeff,DiffCoef];
    else
        display(1)
    end
    
    else
    end
    clear res
    end
end
    
    dat_pos_tim=dat_pos_tim(~cellfun('isempty',dat_pos_tim));
    delT=delT(~cellfun('isempty',delT));
    MSD=MSD(~cellfun('isempty',MSD));
    xk1=xk1(~cellfun('isempty',xk1));
    yk1=yk1(~cellfun('isempty',yk1));
    
%     f = fit(rg_x',rg_hist','exp1');
%     avg=abs(1./f.b);
%     display(avg)
    
    rg_kch{k_choose}=rg;
    alpha_kch{k_choose}=alpha;
    Diff_coeff_kch{k_choose}=Diff_coeff;
    position_dat{k_choose}=dat_pos_tim;
    xk_dat{k_choose}=xk1;
    yk_dat{k_choose}=yk1;
    
%     figure(k_choose.*1000), hold on, 
%     cellfun(@(x) plot(x, 'Color',kcol), MSD);
%     publication_fig(0,0,1);
%     xlabel('\tau (Frame)'); ylabel('Mean Squared Displacement (\mum^2)');

%     figure(k_choose.*10000), hold on, 
%     cellfun(@(x) plot(x, 'Color',kcol), dat_pos_tim);
%     publication_fig(0,0,1);
%     xlabel('Position (\mum)'); ylabel('Time (Frames)');

    display(mean(alpha))
    
    clear xk1 yk1 delT MSD rg_hist pos tim dat_pos_tim xk_yk frac_plus_time frac_minus_time x_hist_plus x_hist_minus rg_hist rg_x;
    
    position=position_dat{k_choose};
    Diffusion_coefficient=Diff_coeff_kch{k_choose};
    xk1=xk_dat{k_choose};
    yk1=yk_dat{k_choose};
    
    save([save_dir, 'position'],'position','Diffusion_coefficient','xk1','yk1');
    
    clear position Diffusion_coefficient xk1 yk1
    
    proc_runs=[];
    proc_times=[];
    nvproc=[];
    rg=[];
    alpha=[];
    Diff_coeff=[];
end

%% Plot alpha distribution:
alpha_ctrl=alpha_kch{1};
alpha_DCXKO=alpha_kch{2};

alpha_x=0.2:0.04:2;
alpha_hist_ctrl=hist(alpha_ctrl,alpha_x);
alpha_hist_DCXKO=hist(alpha_DCXKO,alpha_x);
figure(2), hold on, 
bar(alpha_x, alpha_hist_DCXKO./sum(alpha_hist_DCXKO), 'BarWidth',1,'FaceColor',[0, 0.75, 0.75],'edgecolor',[0, 0.75, 0.75],'Facealpha',1);
hold on,
stairs(alpha_x, alpha_hist_ctrl./sum(alpha_hist_ctrl), 'LineWidth',2,'Color','k');
hold on, 
xlabel('\alpha'); ylabel('Number of Trajectories');
xlim([0 2.1]);
publication_fig(0,0,1);

%% Plot Rg distribution:
rg_ctrl=rg_kch{1};
rg_DCXKO=rg_kch{2};

%rg_x=0.01:0.03:2.5;
rg_x=0:0.03:2.5;
rg_hist_ctrl=hist(rg_ctrl,rg_x);
rg_hist_DCXKO=hist(rg_DCXKO,rg_x);
figure(1), hold on, 
bar(rg_x, rg_hist_DCXKO./sum(rg_hist_DCXKO), 'BarWidth',1,'FaceColor',[0, 0.75, 0.75],'edgecolor',[0, 0.75, 0.75],'Facealpha',1);
hold on,
stairs(rg_x, rg_hist_ctrl./sum(rg_hist_ctrl), 'LineWidth',2,'Color','k');
hold on, 
xlabel('Radius of Gyration (\mum)'); ylabel('Number of Trajectories');
xlim([0 2.1]);
publication_fig(0,0,1);

