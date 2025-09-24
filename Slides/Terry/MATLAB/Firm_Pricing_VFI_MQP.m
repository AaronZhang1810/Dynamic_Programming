%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Firm price adjustment problem
% Implementing VFI solution with MQP bounds
% Firm_Pricing_VFI_MQP.m
%
% S. Terry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('%%% Implement VFI solution with MQP bounds')
disp(' ')
tic

%%%%%%%%%% do the VFI loop

%initialize guess for value function
Vold = zeros(statenum,1);
polindold = zeros(statenum,1);
polold = polindold;

for vfct=1:maxit

    %based on guess for value function, form continuation value EVmat
    %with (i,j) entry equal to the EV(') at state i given policy j
    EVmat = kron((reshape(Vold,pnum,mnum)*(pr_mat_m'))',ones(pnum,1));

    %form RHS matrix with static and continutation payoff for state i and
    %policy j
    RHSmat = Rmat + beta*EVmat;

    %actually do the optimization
    [V,polind] = max(RHSmat,[],2);
    pol = pgrid(polind);

    %now, compute the error in the value function from one iteration to the next
    solerr = max(abs(V(:)-Vold(:)));
    polerr = max(abs(pol(:)-polold(:)));

    %display some diagnostics
    if (mod(vfct,5)==1)
        disp(['For iter ' num2str(vfct) ', VF error = ' num2str(solerr)])
    end

    %if VF has converged, break loop
    if (solerr<soltol)
        break
    end

    %otherwise, compute MQP bounds ot use in update of VF and continue

    minshift = (beta/(1-beta))*min(V(:)-Vold(:));
    maxshift = (beta/(1-beta))*max(V(:)-Vold(:));
    Vold = V + (maxshift+minshift)/2;
    polindold=polind;
    polold = pol;

end

%%%%%%%%%% some simple diagnostics
disp(' ')
disp('Quick policy check: p lb, min(pol), max(pol), p ub')
disp(num2str([pgrid(1) min(pol) max(pol) pgrid(pnum)]))
disp(' ')

%%%%%%%%%% some simple plots
V = reshape(V,pnum,mnum);
pol = reshape(pol,pnum,mnum);
polind = reshape(polind,pnum,mnum);

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


toc
disp(' ')
