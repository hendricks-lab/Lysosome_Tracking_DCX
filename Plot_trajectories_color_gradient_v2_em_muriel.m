%% Analyze Trajectories by color gradient

%Emily modifying to suit Muriel's data on Sept 17/22

close all; clear all; clc; 
set(0,'DefaultFigureWindowStyle','docked')
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/AC_codes_epmodified_20200603/');
addpath('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Analysis/New_General_codes/plotSpread/');

min_lr=1.0;
 
n_int=1;

% Create a color gradient matrix:
min_disp=0:n_int:13;
cmapp1=gray.*[0.902 0.902 0.988]; %for CTRL
% cmapp1=gray.*[0.008 0.737 0.745]; % for DCXKO

cmapp2=flipud(cmapp1);  % Color map goes from dark to light
gradient_matsz=size(cmapp2);
resiz_mat_1=round(gradient_matsz(1)./numel(min_disp));
resize_mat=1:resiz_mat_1:gradient_matsz(1);
cmapp5=[cmapp2(resize_mat,1),cmapp2(resize_mat,2),cmapp2(resize_mat,3)];
cmapp6=flipud(cmapp5);
cmapp6=[cmapp6];

tic
for k_choose = 1

% if k_choose == 1    % WT - Lysotracker
    col1=[0.902 0.902 0.988]; % for CTRL
%     col1=[0.008 0.737 0.745]; % for DCXKO
    cd('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Muriel xml trackmate files/Compiled Lysotracker/CTRL0/mats_2/');
%     cd('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Muriel xml trackmate files/Compiled Lysotracker/DCXKO1/mats/');
    fl='position.mat';
    an=0;
    DT=0.5;
    kf=0;
    proc_runs=[];
    proc_times=[];
    rng=[];
    tlp=-0.5;
    all_pos=[];
    
% end
save_fig_final='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Muriel xml trackmate files/Compiled Lysotracker/';

fls=dir(fl);
    
for k=1:numel(fls)
    display(fls(k).name)
    load(fls(k).name);
    %ab = ab(~cellfun(@isempty,{ab.position}));
    position=position(~cellfun(@isempty,(position)));
    for j=1:numel(position)
        kf=kf+1;
        range_disp(kf)=range(position{j});
        dat_pos_tim{kf}=position{j};
        %all_pos=[all_pos;ab(j).position];
    end
    
    for id=1:(numel(min_disp)+1)
        if id<min_disp(end)
        jk = find(range_disp > (min_disp(id)-(n_int./2)) & range_disp < (min_disp(id)+(n_int./2)));
        display(min_disp(id))
        display(min_disp(id)+(n_int./2))
        jk_ind{id}=jk;
        range_find{id}=range_disp(jk);
        data_position{id}=dat_pos_tim(jk);
        elseif id>min_disp(end)
        jk=find(range_disp > (min_disp(end)));
        jk_ind{id}=jk;
        range_find{id}=range_disp(jk);
        data_position{id}=dat_pos_tim(jk);
        end
    end
     
end
    
    dat_pos_tim(~cellfun('isempty',dat_pos_tim));
    data_position(cellfun('isempty', data_position)) = {{[0,0]}};
    all_pos=cell2mat(dat_pos_tim');
    
    %% Bin all the position points and plot them:
    pos_x=-26:1:26;
    pos_bin=hist(all_pos,pos_x);
    
    figure(2), hold on,
    bar(pos_x,pos_bin,'FaceColor',col1,'Facealpha',1,'BarWidth',0.5);
    set(gca, 'YScale', 'log');
    xlabel('Position (\mum)'); ylabel('Frequency'); 
    publication_fig(0,0,1);
    xlim([-26 26]);
    
    display(numel(dat_pos_tim))
    
    data_position=flipud(data_position');
%     
    for k2=1:numel(data_position)
        hold on,
        figure(1); %hold on,
        %ax = axes;
        %set(ax, 'ColorOrder', [get(ax,'ColorOrder'); 0,0,0]);
        %hold(ax, 'on')
        cellfun(@(x) plot(x, 'Color',cmapp6(k2,:)), data_position{k2});
        axisHandle = gca; 
        set(findall(gca, 'Type', 'Line'),'LineWidth',2);
        set(gca,'LineWidth',2);
        set(gca,'FontSize',24);
        set(gca, 'FontName', 'Arial');
        set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
        set(gca,'Box','on');
        xlabel('Time (Frame)'); ylabel('Position (\mum)');
        ylim([-26 26]); xlim([0 250]);
        pbaspect([2 1 1]);
        hold on,
       
        
        hold on,
        figure(3); %hold on,
        
        subplot(1,7,7)
        plot(pos_x,pos_bin,'Color', col1);
        pbaspect([5.4 1 1]);
%         histogram(all_pos,1000, 'EdgeColor',col1,'FaceColor','none','Normalization', 'Probability');
%         histogram(pos_bin,pos_x,'EdgeColor',col1,'FaceColor','none','Normalization','Probability');
%         stairs(pos_x,pos_bin,'Color',col1, 'LineWidth', 2);
%         bar(pos_x,pos_bin,'FaceColor',col1,'BarWidth',0.5);
        set(gca, 'YScale', 'log');
        ylabel('Frequency'); 
        axisHandle = gca; 
        set(findall(gca, 'Type', 'Line'),'LineWidth',2);
        set(gca,'Box','off');
        set(gca,'FontSize',17);
        set(gca, 'FontName', 'Arial');
        set(gca, 'Color', 'none');
        set(gca,'XColor','none');
        set(gca,'YColor','none');
        xlim([-26 26]);
        view([90 -90]);

        subplot(1,7,[1,2,3,4,5,6])
        cellfun(@(x) plot(x, 'Color',cmapp6(k2,:)), data_position{k2});
        axisHandle = gca; 
        set(findall(gca, 'Type', 'Line'),'LineWidth',2);
        set(gca,'LineWidth',2);
        set(gca,'FontSize',17);
        set(gca, 'FontName', 'Arial');
        set(gca,'XColor',[0 0 0],'YColor',[0 0 0]);
        set(gca,'Box','on');
      
        xlabel('Time (Frame)'); ylabel('Position (\mum)');
        ylim([-26 26]); xlim([0 250]);
        pbaspect([2 1 1]);
        hold on,
        
        
    end
    fl_nm=['20230113_CTRL_', 'all_Trajectories_color_coded'];
    saveas(gca, fullfile(save_fig_final, fl_nm), 'svg');
    saveas(gca, fullfile(save_fig_final, fl_nm), 'fig');
%     fl_nm=['2023_DCXKO', fullfile(save_fig_final, fl_nm), 'svg');
%     saveas(gca,fl_nm,'png')

    clear res all_pos pos_x pos_bin;
    clear pos tim dat_pos_tim xk_yk lr_p_2 lr_m_2 t_p_2 t_m_2 data_position range_find jk_ind dat_pos_tim;
    clear range_disp
    
    proc_runs=[];
    proc_times=[];
    nvproc=[];
    rng=[];
    all_pos=[];
    clear fl_nm
    %close all
end

toc
