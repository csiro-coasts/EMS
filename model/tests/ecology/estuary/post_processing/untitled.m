% this scipt attmep to do a budget for the estuary test case model
% post-processing


%--------------------------------------------------------------------------
%% read the TOTAL mass of the traces (using the TOTALS diagnostic in SHOC)
%--------------------------------------------------------------------------

totals=dlmread('totals.ts',' ',99);%(need to adapt the 99 to the header size on the ts files)

DIN_total=totals(:,17);
time_t=totals(:,1)+datenum('01/01/1990');
% plot the total DIN in the model domain

plot(time_t,DIN/(1000*1000*1000))
datetick

%--------------------------------------------------------------------------
%% get the fluxes trought the boundaries (using the NRTSTATS  section fluxes diagnostic in SHOC)
%--------------------------------------------------------------------------
ocean_B=dlmread('out/bdr1.ts',' ',78); % this is the open ocean boundary (need to adapt the 78 to the header size on the ts files)
river_B=dlmread('out/bdr2.ts',' ',54); % this is the open ocean boundary  .....

time_b=ocean_B(:,1)+datenum('01/01/1990');
DIN_ocean_B=ocean_B(:,7);
DIN_river_B=river_B(:,7);

figure
plot(time_b,DIN_ocean_B./(1000*1000*1000),time_b,DIN_river_B./(1000*1000*1000))
datetick

hold on
plot(tine_t,diff(DIN);
pl

