%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Firm price adjustment problem
% Implementing baseline simulation
% Firm_Pricing_Simulate.m
%
% S. Terry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('%%% Implement baseline simulation')
disp(' ')
tic

%%%%%%%%%% simulate the exogenous process m using inverse CDF sampling
disp('Simulating exogenous process')

%draw the U(0,1) random shocks for simulation of exogenous m
%recall that 'Ttot' is the total number of periods in simulation
firmushocks = rand(firmnum,Ttot);

%form cumulative probability bins for comparison
%for a matrix A, cumsum(A,2) returns the cumulative sum along the rows of A
%each row of 'pr_thresh_m' is the conditional CDF of m' given m
pr_thresh_m = cumsum(pr_mat_m,2);

%initialize each firm's m in the first period to specified value
msimpos = zeros(firmnum,Ttot);
t=1;
msimpos(:,t) = minit;

%loop over firms and periods
for firmct=1:firmnum
    for t=2:Ttot

      %what was the last value of m?
      mctmin1 = msimpos(firmct,t-1);

      %what is today's uniform shock?
      tshockval = firmushocks(firmct,t);

      %compare the shock to the appropriate thresholds
      lessvec = find(tshockval<=pr_thresh_m(mctmin1,:));
      mct = lessvec(1);

      msimpos(firmct,t) = mct;

    end
end

msim = mgrid(msimpos);

%%%%%%%%%% simulate the endogenous pricing process using the policy function
disp('Simulating endogenous process')

%initialize each firm's price (p_{-1}) in the first period to specified value
psimpos = zeros(firmnum,Ttot);
t=1;
psimpos(:,t) = pinit;

%loop over firms and periods
for firmct=1:firmnum
   for t=2:Ttot

      %extract states
      pmin1ct = psimpos(firmct,t-1);
      mct = msimpos(firmct,t-1);

      %store policy
      psimpos(firmct,t) = polind(pmin1ct,mct);

   end
end

psim = pgrid(psimpos);

%%%%%%%%%% throw away burn-in period
%ergsamp = (Terg+2):(Ttot);
ergsamp = (Terg+1):(Ttot-1);
msimpos = msimpos(:,ergsamp);
msim = msim(:,ergsamp);
psimpos = psimpos(:,ergsamp);
psim = psim(:,ergsamp);

%%%%%%%%%% some simple diagnostics
disp(' ')
disp('Quick simulation check: p lb, min(psim), max(psim), p ub')
disp(num2str([pgrid(1) min(psim(:)) max(psim(:)) pgrid(pnum)]))
disp(' ')


%%%%%%%%%% plot some sample paths and distributions
figure;
plot((1:Tsim),psim(1,:),'b',...
    (1:Tsim),psim(2,:),'g',...
    (1:Tsim),psim(3,:),'r',...
    'LineWidth',lwidnum)
xlabel('Time t')
title('Firm Price')
ylabel('p_t')
axis([1 Tsim min(psim(:))*0.95 max(psim(:))*1.05])
set(gca,'FontSize',fsizenum)

figure;
plot((1:Tsim),msim(1,:),'b',...
    (1:Tsim),msim(2,:),'g',...
    (1:Tsim),msim(3,:),'r',...
    'LineWidth',lwidnum)
xlabel('Time t')
title('Firm Marginal Cost')
ylabel('m_t')
axis([1 Tsim min(msim(:))*0.95 max(msim(:))*1.05])
set(gca,'FontSize',fsizenum)

figure;
hist(msim(:))
title('Simulated Distribution: Firm Marginal Cost')
ylabel('Simulated Frequency')
xlabel('m')
set(gca,'FontSize',fsizenum)

figure;
hist(psim(:))
title('Simulated Distribution: Firm Price')
ylabel('Simulated Frequency')
xlabel('p')
set(gca,'FontSize',fsizenum)

figure;
scatter(msim(:),psim(:),'b')
title('Simulated Scatter Plot')
ylabel('Firm Price p')
xlabel('Firm Marginal Cost m')
set(gca,'FontSize',fsizenum)

plag = psim(:,1:(end-1));
pnow = psim(:,2:end);
figure;
scatter(plag(:),pnow(:),'b'); hold on;
plot(pgrid,pgrid,'k','LineWidth',lwidnum);
axis([min(psim(:))*0.95 max(psim(:))*1.05 min(psim(:))*0.95 max(psim(:))*1.05])
title('Simulated Scatter Plot')
ylabel('Firm Price p')
xlabel('Firm Lagged Price p_{-1}')
set(gca,'FontSize',fsizenum)

toc
disp(' ')
