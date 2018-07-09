
close all
clear all
% this script tries to do a budget for the estuary test case model
% post-processing
% read the netcdf file to get the total biomass of all tracers 

% the /1000/1000/1000 is to get from mg to tonnes


%%%% OUTLINE OF THE ROUTINE


%%  1:) READ THE model output (netcdf) and sum the tracer (for water colum sediment and benthic)
%   over a predefine region
%   sum the water column the sediment and the benthic to get the total mass of tracer

%%  2:) compute the fluxes of tracer trought the boundary 

%%  3:) do some plotting


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% get the sum of all the tracer in the region (whole domain or region

%regio ij list of the reg1
31 10
32 10
33 10
34 10
31 11
32 11
33 11
34 11
31 12
32 12
33 12
34 12
31 13
32 13
33 13
34 13
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% First initialise the eco struct
s = eco_init;
% Optional name
s.name = 'estuary_test_case ';
% Specify files
s.files = {'/home/mon126/ems_test_case/bgc_test_case/3d/esturary/out/out1.nc'};
% Specify variables names 
s.vars = {...
     'NO3'  ...          % 1
    ,'NH4'  ...          % 2
    ,'DIP'  ...          % 3
    ,'DOR_N'...          % 4  
    ,'DOR_P'...          % 5
    ,'DetR_N'...         % 6 
    ,'DetPL_N'...        % 7
    ,'DetBL_N'...        % 8 
    ,'DetR_C' ...        % 9  
    ,'DetR_P' ...        % 10
    ,'PhyS_N'...         % 11
    ,'PhyL_N'...         % 12
    ,'ZooS_N'...         % 13
    ,'ZooL_N'...         % 14
    ,'PhyL_N_pr'...      % 15
    ,'PhyS_N_pr'...      % 16
    ,'ZooL_N_rm'...      % 17
    ,'ZooS_N_rm'...      % 18
    ,'Den_fl'...         % 19
    ,'NH4_pr'...         % 20
    ,'MPB_N_pr' ...      % 21
    ,'PhyD_N_pr'...      % 22
    ,'Oxy_pr'...         % 23
    ,'Nfix'...           % 24
    ,'TN'        ...     % 25
    ,'TP'  ...           % 26
    ,'MPB_N'...          % 27
    ,'ptras'...          %28
    ,'PhyD_N'...         % 29
    }

% Specify what function to perform
s.function = 'sum'; % mean is the default
% Call the post processor to calculate means
s = eco_post_process(s); % This updates the struct s
totals=s.region(1).data{1};
clear s


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do the same for sediments
s1 = eco_init;
% Optional name
s1.name = 'estuary_test_case ';
s1.files = {'/home/mon126/ems_test_case/bgc_test_case/3d/esturary/out/out1.nc'};
% Specify variables names 
s1.vars = {...
     'NO3_sed'  ...          % 1
    ,'NH4_sed'  ...          % 2
    ,'DIP_sed'  ...          % 3
    ,'DOR_N_sed'...          % 4  
    ,'DOR_P_sed'...          % 5
    ,'DetR_N_sed'...         % 6 
    ,'DetPL_N_sed'...        % 7
    ,'DetBL_N_sed'...        % 8 
    ,'DetR_C_sed' ...        % 9  
    ,'DetR_P_sed' ...        % 10
    ,'PhyS_N_sed'...         % 11
    ,'PhyL_N_sed'...         % 12
    ,'ZooS_N_sed'...         % 13
    ,'ZooL_N_sed'...         % 14
    ,'Den_fl_sed'...         % 15
    ,'NH4_pr_sed'...         % 16
    ,'TN_sed'...             % 17
    ,'TP_sed'...             % 18
    ,'MPB_N_sed'...          % 19
    ,'PhyL_N_pr_sed'...      % 20
    ,'PhyS_N_pr_sed'...      % 21
    ,'ZooL_N_rm_sed'...      % 22
    ,'ZooS_N_rm_sed'...      % 23
    ,'PhyD_N_sed'...         % 24
    }
% Specify what function to perform
s1.function = 'sum'; % mean is the default
% Call the post processor to calculate means
s1= eco_post_process(s1); % This updates the struct s
totals_sed=s1.region(1).data{1};

% get the time in matlab format
time_t=s1.t./86400+datenum('01/01/2000');

clear s1
clear s
%% do the same for benthic tracer
s2 = eco_init;
% Optional name
s2.name = 'estuary_test_case ';
s2.files = {'/home/mon126/ems_test_case/bgc_test_case/3d/esturary/out/out1.nc'};
% Specify variables names 
s2.vars = {...
     'MA_N'  ...          % 1
    ,'SG_N'  ...          % 2
    ,'EpiTN' ...          % 3
    ,'EpiTP' ...          % 4
    }
% Specify what function to perform
s2.function = 'sum'; % mean is the default
% Call the post processor to calculate means
s2= eco_post_process(s2); % This updates the struct s
totals_epi=s2.region(1).data{1};



% extract wc tracers
%TNm_wc     =totals(:,25)./(1000*1000*1000); % 17 is the column  number corresponding to the data
no3_wc     =totals(:,1)./(1000*1000*1000); % 17 is the column  number corresponding to the data
nh4_wc     =totals(:,2)./(1000*1000*1000); % 17 is the column  number corresponding to the data
dor_n_wc   =totals(:,4)./(1000*1000*1000); % 17 is the column  number corresponding to the data
detpl_n_wc =totals(:,7)./(1000*1000*1000); % 17 is the column  number corresponding to the data
detbl_n_wc =totals(:,8)./(1000*1000*1000); % 17 is the column  number corresponding to the data
detr_n_wc  =totals(:,6)./(1000*1000*1000); % 17 is the column  number corresponding to the data
phys_n_wc  =totals(:,11)./(1000*1000*1000); % 17 is the column  number corresponding to the data
phyl_n_wc  =totals(:,12)./(1000*1000*1000); % 17 is the column  number corresponding to the data
zool_n_wc  =totals(:,14)./(1000*1000*1000); % 17 is the column  number corresponding to the data
zoos_n_wc  =totals(:,13)./(1000*1000*1000); % 17 is the column  number corresponding to the data
mpb_n_wc   =totals(:,27)./(1000*1000*1000); % 17 is the column  number corresponding to the data
ptras_wc   =totals(:,28)./(1000*1000*1000); % 17 is the column  number corresponding to the data
dip_wc     =totals(:,3)./(1000*1000*1000); % 17 is the column  number corresponding to the data
phyd_n_wc  =totals(:,29)./(1000*1000*1000); % 17 is the column  number corresponding to the data

% exctract wc fluxes tracer
phys_prod_n_wc  =totals(:,16)./(1000*1000*1000);
phyl_prod_n_wc  =totals(:,15)./(1000*1000*1000);
zoos_rem_n_wc   =totals(:,18)./(1000*1000*1000);
zool_rem_n_wc   =totals(:,17)./(1000*1000*1000);

Den_fl_wc       =totals(:,19)./(1000*1000*1000);
NH4_pr_wc       =totals(:,20)./(1000*1000*1000);


% extract sed tracers
dip_sed             =totals_sed(:,3)./(1000*1000*1000); % 17 is the column  number corresponding to the data
no3_sed             =totals_sed(:,1)./(1000*1000*1000); % 17 is the column  number corresponding to the data
nh4_sed             =totals_sed(:,2)./(1000*1000*1000); % 17 is the column  number corresponding to the data
dor_n_sed           =totals_sed(:,4)./(1000*1000*1000); % 17 is the column  number corresponding to the data
detpl_n_sed         =totals_sed(:,7)./(1000*1000*1000); % 17 is the column  number corresponding to the data
detbl_n_sed         =totals_sed(:,8)./(1000*1000*1000); % 17 is the column  number corresponding to the data
detr_n_sed          =totals_sed(:,6)./(1000*1000*1000); % 17 is the column  number corresponding to the data
%TNm_sed             =totals_sed(:,17)./(1000*1000*1000); % 17 is the column  number correspond
mpb_n_sed           =totals_sed(:,19)./(1000*1000*1000); % 17 is
phyd_n_sed          =totals_sed(:,24)./(1000*1000*1000); % 17 is the column  number corresponding to the data

phys_n_sed  =totals_sed(:,11)./(1000*1000*1000); % 17 is the column  number corresponding to the data
phyl_n_sed  =totals_sed(:,12)./(1000*1000*1000); % 17 is the column  number corresponding to the data
zool_n_sed  =totals_sed(:,14)./(1000*1000*1000); % 17 is the column  number corresponding to the data
zoos_n_sed  =totals_sed(:,13)./(1000*1000*1000); % 17 is the column  number corresponding to the data

phys_prod_n_sed   =totals_sed(:,21)./(1000*1000*1000);
phyl_prod_n_sed   =totals_sed(:,20)./(1000*1000*1000);
zoos_rem_n_sed    =totals_sed(:,23)./(1000*1000*1000);
zool_rem_n_sed    =totals_sed(:,22)./(1000*1000*1000);


Den_fl_sed    =totals_sed(:,15)./(1000*1000*1000);
NH4_pr_sed    =totals_sed(:,16)./(1000*1000*1000);


% extract benthic tracers
MA_epi               =totals_epi(:,1)./(1000*1000*1000); % 17 is the column  number corresponding to the data
SG_epi               =totals_epi(:,2)./(1000*1000*1000); % 17 is the column  number corresponding to the data
%TNm_epi              =totals_epi(:,3)./(1000*1000*1000); % 17 is
%TPm_epi              =totals_epi(:,4)./(1000*1000*1000); % 17 is



%% sum the sedi and wc tracer
no3_total=no3_sed+no3_wc;
nh4_total=nh4_sed+nh4_wc;
dip_total=dip_sed+dip_wc;

detpl_n_total  =detpl_n_sed+detpl_n_wc;
detr_n_total   =detr_n_sed+detr_n_wc;
detbl_n_total  =detbl_n_sed+detbl_n_wc;

dor_n_total   =dor_n_sed+dor_n_wc;

phys_n_total   =  phys_n_sed+phys_n_wc;
phyl_n_total   =  phyl_n_sed+phyl_n_wc;
zool_n_total   =  zool_n_sed+zool_n_wc;
zoos_n_total   =  zoos_n_sed+zoos_n_wc;

phys_prod_n_total    = phys_prod_n_sed+phys_prod_n_wc;
phyl_prod_n_total    = phyl_prod_n_sed+phyl_prod_n_wc;

zoos_rem_n_total   =  zoos_rem_n_sed+zoos_rem_n_wc;
zool_rem_n_total   =  zool_rem_n_sed+zool_rem_n_wc;



%TNm_total=TNm_sed+TNm_wc+TNm_epi;

%% do my own TN calculation

TN_wc=zool_n_wc+...
    zoos_n_wc+...
    phyl_n_wc+...
    phys_n_wc+....
    detbl_n_wc+...
    detpl_n_wc+...
    detr_n_wc+...
    no3_wc+...
    nh4_wc +...
    dor_n_wc+...
    mpb_n_wc;


TN_sed=zool_n_sed+...
    zoos_n_sed+...
    phyl_n_sed+...
    phys_n_sed+....
    detbl_n_sed+...
    detpl_n_sed+...
    detr_n_sed+...
    no3_sed+...
    nh4_sed +...
    dor_n_sed+...
    mpb_n_sed;

TN_epi=MA_epi+SG_epi;

TN_total=TN_wc+TN_sed+TN_epi;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% get the fluxes trought the boundaries (using the NRTSTATS  section fluxes diagnostic in SHOC)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% header number of line in the ts files has to be done manualy 
header_totals=129; % total ts files
header_bdr1=174;    % open boundaries 
header_bdr2=150;    % open boundaries 

ocean_B=dlmread('/home/mon126/ems_test_case/bgc_test_case/3d/esturary/out/bdr1.ts','',header_bdr1); % this is the open ocean boundary (need to adapt the 78 to the header size on the ts files)
river_B=dlmread('/home/mon126/ems_test_case/bgc_test_case/3d/esturary/out/bdr2.ts','',header_bdr2); % this is the open ocean boundary  .....

% get the time
time_b=ocean_B(:,1)+datenum('01/01/1990');

%get the flux odf all the variable trought the boundaries
no3_ocean_B=ocean_B(:,2)./(1000*1000*1000);  % 3 is the column number corresponding to the data
no3_river_B=river_B(:,2)./(1000*1000*1000);

dip_ocean_B=ocean_B(:,13)./(1000*1000*1000);  % 3 is the column number corresponding to the data
dip_river_B=river_B(:,13)./(1000*1000*1000);

nh4_ocean_B=ocean_B(:,3)./(1000*1000*1000);  % 3 is the column number corresponding to the data
nh4_river_B=river_B(:,3)./(1000*1000*1000);

DetR_N_ocean_B=ocean_B(:,5)./(1000*1000*1000);  % 6 is the column number corresponding to the data
DetR_N_river_B=river_B(:,5)./(1000*1000*1000);

DetPL_N_ocean_B=ocean_B(:,8)./(1000*1000*1000);  % 6 is the column number corresponding to the data
DetPL_N_river_B=river_B(:,8)./(1000*1000*1000);

DetBL_N_ocean_B=ocean_B(:,7)./(1000*1000*1000);  % 6 is the column number corresponding to the data
DetBL_N_river_B=river_B(:,7)./(1000*1000*1000);

DOR_N_ocean_B=ocean_B(:,6)./(1000*1000*1000);  % 6 is the column number corresponding to the data
DOR_N_river_B=river_B(:,6)./(1000*1000*1000);

phys_N_ocean_B=ocean_B(:,9)./(1000*1000*1000);  % 6 is the column number corresponding to the data
phys_N_river_B=river_B(:,9)./(1000*1000*1000);

phyl_N_ocean_B=ocean_B(:,10)./(1000*1000*1000);  % 6 is the column number corresponding to the data
phyl_N_river_B=river_B(:,10)./(1000*1000*1000);

zoos_N_ocean_B=ocean_B(:,11)./(1000*1000*1000);  % 6 is the column number corresponding to the data
zoos_N_river_B=river_B(:,11)./(1000*1000*1000);

zool_N_ocean_B=ocean_B(:,12)./(1000*1000*1000);  % 6 is the column number corresponding to the data
zool_N_river_B=river_B(:,12)./(1000*1000*1000);

ptras_ocean_B=ocean_B(:,18)./(1000*1000*1000);  % 6 is the column number corresponding to the data
ptras_river_B=river_B(:,18)./(1000*1000*1000);

mpb_ocean_B=ocean_B(:,14)./(1000*1000*1000);  
mpb_river_B=river_B(:,14)./(1000*1000*1000);

TN_ocean_B=zool_N_ocean_B+zoos_N_ocean_B+phyl_N_ocean_B+phys_N_ocean_B+DetBL_N_ocean_B...
  +  DetPL_N_ocean_B+DetR_N_ocean_B+ no3_ocean_B+nh4_ocean_B +DOR_N_ocean_B+mpb_ocean_B;

TN_river_B=zool_N_river_B+zoos_N_river_B+phyl_N_river_B+phys_N_river_B+DetBL_N_river_B...
  +  DetPL_N_river_B+DetR_N_river_B+no3_river_B+nh4_river_B +DOR_N_river_B+mpb_river_B;



%%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% plot the total biomass i the model domain
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
figure
% plot the total N in the whole model domain
subplot(2,1,1),plot(time_t,TN_sed,...
                    time_t,phys_n_sed,...
                    time_t,phyl_n_sed,...
                    time_t,zool_n_sed,...
                    time_t,zoos_n_sed,...
                    time_t,detr_n_sed,...
                    time_t,detpl_n_sed,...
                    time_t,detbl_n_sed)
legend('TN','phys','phyl','zool','zoos','detr','detrpl','detbl')    
title('Total biomass in the sed region')
  datetick  
  
 subplot(2,1,2),plot(time_t,no3_sed,...
                      time_t,nh4_sed,...
                      time_t,dor_n_sed)
legend('no3','nh4','dor')          
datetick
print -dpsc total_mass_sed_domain.ps

figure
% plot the total N in the whole model domain
subplot(2,1,1),plot(time_t,TN_wc,...
                    time_t,phys_n_wc,...
                    time_t,phyl_n_wc,...
                    time_t,zool_n_wc,...
                    time_t,zoos_n_wc,...
                    time_t,detr_n_wc,...
                    time_t,detpl_n_wc,...
                    time_t,detbl_n_wc)
legend('TN','phys','phyl','zool','zoos','detr','detrpl','detbl')    
title('Total biomass in the wc region')
  datetick  
  
 subplot(2,1,2),plot(time_t,no3_wc,...
                      time_t,nh4_wc,...
                      time_t,dor_n_wc)
legend('no3','nh4','dor')          
datetick
print -dpsc total_mass_wc_domain.ps


% plot the total fluxes in the model domain

figure
%------------------------------phy prod and zoo removal
plot(time_t,phys_prod_n_total,...
     time_t,zoos_rem_n_total,...
     time_t,phyl_prod_n_total,... 
     time_t,zool_rem_n_total);

legend(['phys prod     ',num2str(sum(phys_prod_n_total))],...
       ['zoos removal  ',num2str(sum(zoos_rem_n_total))],...
       ['phyl prod     ',num2str(sum(phyl_prod_n_total))],... 
       ['zool removal  ',num2str(sum(zool_rem_n_total))])
datetick
print -dpsc produc_fluxes.ps
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% plot the input/ output of tracers trought the boundaries
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%-------------------------------------------------------------------------
figure
%------------------------------DIN
subplot(4,1,1),plot(time_b,no3_ocean_B*24); 
title('Flux of No3 trought the ocean boundary in Tonne per day');
datetick
subplot(4,1,2),plot(time_b,no3_river_B*24); 
title('Flux of no3 trought the river  in Tonne per day');
datetick
subplot(4,1,3),plot(time_b,nh4_ocean_B*24); 
title('Flux of nh4 trought the ocean boundary in Tonne per day');
datetick
subplot(4,1,4),plot(time_b,nh4_river_B*24); 
title('Flux of nh4 trought the river  in Tonne per day');
datetick
print -dpsc input_output_din.ps


figure
%-----------------------------Det
subplot(4,1,1),plot(time_b,DetR_N_ocean_B*24); 
title('Flux of DetR_N trought the ocean boundary in Tonne per day');
datetick
subplot(4,1,2),plot(time_b,DetR_N_river_B*24); 
title('Flux of DetR_N trought the river  in Tonne per day');
datetick
subplot(4,1,3),plot(time_b,DetPL_N_ocean_B*24,time_b,DetBL_N_ocean_B*24); 
title('Flux of DetPL and detBL trought the ocean boundary in Tonne per day');
datetick
subplot(4,1,4),plot(time_b,DetPL_N_river_B*24,time_b,DetBL_N_river_B*24); 
title('Flux of DetPL and detBL  trought the river  in Tonne per day');
datetick
print -dpsc input_output_det.ps


figure
%-----------------------------phy zoo
subplot(4,1,1),plot(time_b,phys_N_ocean_B*24,time_b,phyl_N_ocean_B*24); 
title('Flux of phy trought the ocean boundary in Tonne per day');
datetick
subplot(4,1,2),plot(time_b,phyl_N_river_B*24,time_b,phys_N_river_B*24); 
title('Flux of phy trought the river  in Tonne per day');
datetick
subplot(4,1,3),plot(time_b,zoos_N_ocean_B*24,time_b,zool_N_ocean_B*24); 
title('Flux of DetPL and zoo trought the ocean boundary in Tonne per day');
datetick
subplot(4,1,4),plot(time_b,zool_N_river_B*24,time_b,zoos_N_river_B*24); 
title('Flux of zoo and zoo trought the river  in Tonne per day');
datetick
print -dpsc input_output_phyzoo.ps



figure
%-----------------------------TN
sum1=num2str(sum(TN_ocean_B));

subplot(3,1,1),plot(time_b,TN_ocean_B*24);
title('Flux of TN  trought the ocean  in Tonne per day');
txtar =  annotation('textbox',[0 0  1 1],'string',['total of  ',sum1],'FontSize',14,'FitBoxToText','on');
datetick

sum1=num2str(sum(TN_river_B));

subplot(3,1,2),plot(time_b,TN_river_B)
title('Flux of TN  trought the river  in Tonne per day');
txtar =  annotation('textbox',[0.0 0.4 0.028 0.086],'string',['total of  ',sum1],'FontSize',14,'FitBoxToText','on');
datetick


subplot(3,1,3),plot(time_t,TN_total,time_t,TN_epi,time_t,TN_wc,time_t,TN_sed);
title(' TN in Tonnes');
legend('TM total','TN epi','TN wc','TN sed')
txtar =  annotation('textbox',[0.0 0.2 0.028 0.086],'string',['Total TN diff  ',num2str(sum(TN_total(end)-TN_total(1)))],'FontSize',14,'FitBoxToText','on');
datetick

print -dpsc input_output_tn.ps



figure  % denitr remi pp
%-----------------------------n budget
sum1=num2str(sum(Den_fl_sed+Den_fl_wc)); % denitrification nitrification 
subplot(3,1,1),plot(time_t,Den_fl_sed+Den_fl_wc);
title('Flux of Den_fl in Tonne per day');
txtar =  annotation('textbox',[0 0  1 1],'string',['total of  ',sum1],'FontSize',14,'FitBoxToText','on');
datetick

sum1=num2str(sum(NH4_pr_sed+NH4_pr_wc));  % remineralization detritus to DIN
subplot(3,1,2),plot(time_b,TN_river_B)
title('remineralization detritus to DIN Tonne per day');
txtar =  annotation('textbox',[0.0 0.4 0.028 0.086],'string',['total of  ',sum1],'FontSize',14,'FitBoxToText','on');
datetick

sum1=num2str(sum(phyl_prod_n_sed+phys_prod_n_sed+phyl_prod_n_wc+phys_prod_n_wc));  %  
subplot(3,1,3),plot(time_t,TN_total);
title(' DIN uptkake by phyto Tonnes per day');
datetick
txtar =  annotation('textbox',[0.0 0.2 0.028 0.086],'string',['total of  ',sum1],'FontSize',14,'FitBoxToText','on');

print -dpsc N_budget_1.ps


figure
%-----------------------------n budget phyto production and zoo grazing
sum1=num2str(sum(zoos_rem_n_total)); % 
sum2=num2str(sum(zool_rem_n_total)); % 

subplot(2,1,1),plot(time_t,zoos_rem_n_total,time_t,zool_rem_n_total);
legend('Zoos','Zool')
title('Flux of   zoos and zool removal of phyl and phys in tonnes');
txtar =  annotation('textbox',[0 0  1 1],'string',['total zoos of  ',sum1,'  and  ',sum2],'FontSize',14,'FitBoxToText','on');
datetick

clear sum1 sum2
sum1=num2str(sum(phys_prod_n_total));  %  
sum2=num2str(sum(phyl_prod_n_total));  %  

subplot(2,1,2),plot(time_t,phys_prod_n_total,time_t,phyl_prod_n_total)
legend('phys','phyl')
title(' Flux of   Phys and Phyl production Tonne per day');
txtar =  annotation('textbox',[0.0 0.4 0.028 0.086],'string',['total of  ',sum1,'  and  ',sum2],'FontSize',14,'FitBoxToText','on');
datetick
print -dpsc N_budget_2.ps


