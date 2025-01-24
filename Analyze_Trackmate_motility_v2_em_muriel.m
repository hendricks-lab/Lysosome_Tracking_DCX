%% Abdullah R. Chaudhary
%% This code will only run .csv files
clear all; close all; clc; 
%this code doesn't work for me :P
addpath('/Volumes/Emily_htt/New_General_codes/');
set(0,'DefaultFigureWindowStyle','docked');
addpath('/Volumes/Emily_htt/New_General_codes/');

tic
for k_choose = 1:2
%k_choose = 1;

if k_choose == 1    % Control
    cd('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Muriel xml trackmate files/Compiled Lysotracker/CTRL0/csvs/');
    save_dir='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Muriel xml trackmate files/Compiled Lysotracker/CTRL0/mats_2/';
elseif k_choose == 2    % DCX knockout
    cd('/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Muriel xml trackmate files/Compiled Lysotracker/DCXKO1/csvs/');
    save_dir='/Users/emilyprowse/Documents/McGill/Hendricks_Lab/Muriel xml trackmate files/Compiled Lysotracker/DCXKO1/mats_2/';
end

% Asterik indicates all the files with names starting with Spots and with
% extension .csv
Mot_file='Copy*.csv';
DT=0.5;    % single channel exposure time 500ms for Muriel
%% Variables defined
tm_bw_rev=[];
rn_bw_rev=[];
dat=[];

dmot=dir(Mot_file); % Pick the trajectory files in .mat format

for k=1:length(dmot) % k is a vector of length dmot (no k=4)
    display(k)
dat=readmatrix(dmot(k).name); %load everything in dmot
[filepath,name,ext] = fileparts(dmot(k).name);
display(dmot(k).name);  % display
cmap = colormap(lines(100));
Ndat=dat(end,10);
kp=0;

for kd=1:Ndat
    kp=kp+1;
    jkd=find(dat(:,10)==kd); %indices of track_k

    xk=dat(jkd,12); %x position (um)
    yk=dat(jkd,13); %y position (um)
    timek=[1:numel(xk)].*DT;    % Time (s)
    origin_x=dat(3,15); % Point of origin in x-axis
    origin_y=dat(3,16); % Point of origin in y-axis
    
    %if range(xk)>=0.5 || range(yk)>=0.5
    % Determine median position from kcluster analysis:
    [idx,kclust]=kmeans([xk,yk],1);
    dx_med=(origin_x-kclust(1));
    dy_med=(origin_y-kclust(2));
    
    % Vectors of x and y components:
    dx=[0;diff(xk)]; dy=[0;diff(yk)];
    dx_origin=(origin_x-xk); dy_origin=(origin_y-yk);
    
    for j=1:numel(dx)
    dt_prd(j)=dot([dx(j),dy(j)],[dx_origin(j),dy_origin(j)],2);
    end
    
    mag_vect=sqrt(dx_origin.^2 + dy_origin.^2);
    dt_prd_norm=dt_prd'./mag_vect;
    position=cumsum(dt_prd_norm);
    
    ab(kp).position=position.*-1;
    ab(kp).xk=xk;
    ab(kp).yk=yk;
    ab(kp).tk=timek;
    ab(kp).med_clust_x=dx_med;
    ab(kp).med_clust_y=dy_med;
    ab(kp).origin_x=origin_x;
    ab(kp).origin_y=origin_y;
    clear dt_prd;
    %else
%     ab(kp).position=[];
%     ab(kp).xk=[];
%     ab(kp).yk=[];
%     ab(kp).tk=[];
%     ab(kp).med_clust_x=[];
%     ab(kp).med_clust_y=[];
%     ab(kp).origin_x=[];
%     ab(kp).origin_y=[];
    %end
  
end
    save([save_dir, num2str(k_choose) name],'ab');
    
    clear ab;
end


end
toc


