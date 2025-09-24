%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read_firm_price.m
%
% Read and process output from Fortran
%
% Stephen Terry
%
% This Version: August 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc

%some formatting stuff
lwidnum = 2;
fsizenum = 12;


constants = importdata('constantvec.txt');
mnum = constants(1);
pnum = constants(2);
Ttot = constants(3);
firmnum = constants(4);
statenum = constants(5);
Terg = constants(6);
Tsim = constants(7);

pr_mat_m = importdata('pr_mat_m.txt');
mgrid = importdata('mgrid.txt');
pgrid = importdata('pgrid.txt');

distp = importdata('distp.txt');
distm = importdata('distm.txt');


distvec = importdata('dist.txt');
Vvec = importdata('V.txt');
polvec = importdata('pol.txt');

dist = zeros(pnum,mnum);
V = dist;
pol = dist;

ct=0;
for pmin1ct=1:pnum
   for mct=1:mnum
       ct = ct + 1;
       dist(pmin1ct,mct) = distvec(ct);
       V(pmin1ct,mct) = Vvec(ct);
       pol(pmin1ct,mct) = polvec(ct);
   end
end

msimvec = importdata('msim.txt');
psimvec = importdata('psim.txt');

msim = zeros(firmnum,Ttot);
psim = msim;

ct=0;
for firmct=1:firmnum
   for t=1:Ttot
      ct=ct+1;
      msim(firmct,t) = msimvec(ct);
      psim(firmct,t) = psimvec(ct);
   end
end

%%%%%%%%%% throw away burn-in period
ergsamp = (Terg+1):(Ttot-1);
msim = msim(:,ergsamp);
psim = psim(:,ergsamp);

figure;
plot(pgrid,V(:,1),'b',...
    pgrid,V(:,floor(mnum/2)),'g',...
    pgrid,V(:,mnum),'r',...
    'LineWidth',lwidnum)
xlabel('Old Price p_{-1}')
ylabel('Firm Value V(p_{-1},m)')
axis([pgrid(1) pgrid(pnum) min(V(:))*0.95 max(V(:))*1.05])
set(gca,'FontSize',fsizenum)
legend('Low m','Medium m','High m','FontSize',fsizenum,'Location','northwest')
legend boxoff

figure
plot(pgrid,pol(:,2),'b',...
    pgrid,pol(:,floor(mnum/2)),'g',...
    pgrid,pol(:,mnum),'r',...
    pgrid,pgrid,'k',...
    'LineWidth',lwidnum)
xlabel('Old Price p_{-1}')
ylabel('New Price p(p_{-1},m)')
axis([pgrid(1) pgrid(pnum) pgrid(1) pgrid(pnum)])
set(gca,'FontSize',fsizenum)
legend('Low m','Medium m','High m','Identity','FontSize',fsizenum,'Location','northwest')
legend boxoff


figure;
plot(pgrid,distp,'b',...
    'LineWidth',lwidnum)
xlabel('p')
ylabel('Weight')
title('Ergodic Distribution for Firm Price')
set(gca,'FontSize',fsizenum)

figure;
plot(mgrid,distm,'b',...
    'LineWidth',lwidnum)
xlabel('m')
ylabel('Weight')
title('Ergodic Distribution for Firm Marginal Cost')
set(gca,'FontSize',fsizenum)

figure; 
mesh(mgrid,pgrid,dist);
xlabel('Firm Marginal Cost m')
ylabel('Firm Lagged Price p_{-1}')
title('Joint Ergodic Distribution')
zlabel('Weight')
set(gca,'FontSize',fsizenum)



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




