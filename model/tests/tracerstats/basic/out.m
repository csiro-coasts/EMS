% Matlab script to graph basic test output

ncf = netcdf('out/all.nc');
% WC
salt   = ncf{'salt'}; s = salt(:,1,1,1);
m_salt = ncf{'HM_salt'};  ms = m_salt(:,1,1,1);
v_salt = ncf{'VAR_salt'}; vs = v_salt(:,1,1,1);
std_salt = ncf{'STDEV_salt'}; ss = std_salt(:,1,1,1);
max_salt = ncf{'MAX_salt'}; max_s = max_salt(:,1,1,1);
min_salt = ncf{'MIN_salt'}; min_s = min_salt(:,1,1,1);
tands  = ncf{'SUM_TS'}; ts = tands(:,1,1,1);
temper = ncf{'temp'}; tt = temper(:,1,1,1);
cv_ts = ncf{'COV_TS'}; cv_ts = cv_ts(:,1,1,1);

t = ncf{'t'}(:);

% This bit of code verifies that the mean is in fact numerically
% correct.  We know that dt is 100 secs and that the averaging
% interval is 1 hour which makes 100.041667 (i.e. 1/24) the
% timestep of the mean. This in turn is index 3600/100 = 36 + 1
index = 37;
exp_mean_val = mean(s(2:index)); % 2 as tracerstats starts at t1
                                 % initial value
act_mean_val = ms(index);

% Do the test
if (exp_mean_val ~= act_mean_val)
  disp('---------------------------------');
  disp('!!!!!ERROR: MEAN TEST FAILED!!!!!');
  disp('---------------------------------');
else
  disp('MEAN TEST PASSED');
end

figure(1);
clf

% Means
subplot(321); 
plot(t,ms,t,s);
grid on;
xlabel('time (days)');
ylabel('Salinity');
legend('mean(salt)','salt');
title('Mean');

% Variance/Stdev
subplot(322); 
plot(t,vs,t,s,t,ss.^2+1);
grid on;
xlabel('time (days)');
ylabel('Salinity');
legend('variance(salt)','salt','stdev^2(salt)+1'); % plus one so
                                                   % that you can
                                                   % distiguish
                                                   % between the
                                                   % two graphs
title('Variance/Standard deviation');

% Max
subplot(323); 
plot(t,max_s,'o',t,s);
grid on;
xlabel('time (days)');
ylabel('Salinity');
legend('max(salt)','salt');
title('Max');

% Min
subplot(324);
plot(t,min_s,'x',t,s);
grid on;
xlabel('time (days)');
ylabel('Salinity');
legend('min(salt)','salt');
title('Min');

% Sum
subplot(325);
plot(t,ts,'x',t,s,t,tt,'+');
grid on;
xlabel('time (days)');
ylabel('Salinity + temp');
legend('sum(salt,temp)','salt','temp');
title('Sum');

% Covariance
subplot(326);
plot(t,cv_ts,'x',t,s, t,tt,'+');
grid on;
xlabel('time (days)');
ylabel('');
legend('cov(salt,temp)','salt','temp');
title('Covariance (basicaly uncorrelated)');


%
% SED
%
salt   = ncf{'salt_sed'}; s = salt(:,1,1,1);
m_salt = ncf{'HM_salt_sed'};  ms = m_salt(:,1,1,1);
v_salt = ncf{'VAR_salt_sed'}; vs = v_salt(:,1,1,1);
std_salt = ncf{'STDEV_salt_sed'}; ss = std_salt(:,1,1,1);
max_salt = ncf{'MAX_salt_sed'}; max_s = max_salt(:,1,1,1);
min_salt = ncf{'MIN_salt_sed'}; min_s = min_salt(:,1,1,1);
tands  = ncf{'SUM_TS_sed'}; ts = tands(:,1,1,1);

t = ncf{'t'}(:);

figure(2);
clf

% Means
subplot(321); 
plot(t,ms,t,s);
grid on;
xlabel('time (days)');
ylabel('Salinity');
legend('mean(salt)','salt');
title('Mean');

% Variance/Stdev
subplot(322); 
plot(t,vs,t,s,t,ss.^2+1);
grid on;
xlabel('time (days)');
ylabel('Salinity');
legend('variance(salt)','salt','stdev^2(salt)+1'); % plus one so
                                                   % that you can
                                                   % distiguish
                                                   % between the
                                                   % two graphs
title('Variance/Standard deviation');

% Max
subplot(323); 
plot(t,max_s,'o',t,s);
grid on;
xlabel('time (days)');
ylabel('Salinity');
legend('max(salt)','salt');
title('Max');

% Min
subplot(324);
plot(t,min_s,'x',t,s);
grid on;
xlabel('time (days)');
ylabel('Salinity');
legend('min(salt)','salt');
title('Min');

% Sum
subplot(325);
plot(t,ts,'x',t,s,t,tt,'+');
grid on;
xlabel('time (days)');
ylabel('Salinity + temp');
legend('sum(salt,temp)','salt','temp');
title('Sum');

% EOF

