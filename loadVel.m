%% Load and manipulate all the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loads however many timesteps specified by variable 'times' and puts data
% into a cell array for each velocity component.  Timesteps are cell array
% entries in ascending order.

clear all
close all

times=1; % number of timesteps.
firstTime=460;
velLoc='/Users/Kshitij/Desktop/ALES/DNS_Data/'; % change to your directory

nx=[256,256,256];                          %computational grid dimensions
lx=[0.259,0.259,0.259]/100;                %physical domain size
dx=lx./nx;                                 %physical grid cell width 
x=linspace(0+dx(1)/2,lx(1)-dx(1)/2,nx(1)); %1D x physical grid
y=linspace(0+dx(2)/2,lx(2)-dx(2)/2,nx(2)); %1D y physical grid
z=linspace(0+dx(3)/2,lx(3)-dx(3)/2,nx(3)); %1D z physical grid

for tt = 1:times
    time=num2str(firstTime + (tt-1));
    disp(['Loading velocities from time ' time])
    %Read in HST velocity fields
    [fid,errmsg]=fopen([velLoc,'Velocity1_0',time,'.bin']);
    tmp=single(fread(fid,nx(1)*nx(2)*nx(3),'single','b'));
    fclose(fid);
    Uhit{tt}=reshape(tmp,nx(1),nx(2),nx(3))/100;
    clear tmp ;

    fid=fopen([velLoc,'Velocity2_0',time,'.bin']);
    tmp=single(fread(fid,nx(1)*nx(2)*nx(3),'single','b'));
    fclose(fid);
    Vhit{tt}=reshape(tmp,nx(1),nx(2),nx(3))/100;
    clear tmp ;
    
    fid=fopen([velLoc,'Velocity3_0',time,'.bin']);
    tmp=single(fread(fid,nx(1)*nx(2)*nx(3),'single','b'));
    fclose(fid);
    Whit{tt}=reshape(tmp,nx(1),nx(2),nx(3))/100;
    clear tmp ;
end

%% Visualize a data slice %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vizTime=1;
vizSlice=256;
clims=[-15,15];
figure(2)
subplot(131)
imagesc(squeeze(Uhit{vizTime}(:,:,vizSlice)),clims);
title('U Velocity')
axis square
subplot(132)
imagesc(squeeze(Vhit{vizTime}(:,:,vizSlice)),clims);
title('V Velocity')
axis square
subplot(133)
imagesc(squeeze(Whit{vizTime}(:,:,vizSlice)),clims);
title('W Velocity')
axis square
