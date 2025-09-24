%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Firm price adjustment problem
% Implementing VFI solution with Howard improvement
% Firm_Pricing_VFI_Howard.m
%
% S. Terry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('%%% Implement VFI with Howard improvement')
disp(' ')
tic

%%%%%%%%%% do the VFI loop

%initialize guess for value function & policy
Vold = zeros(statenum,1);
polindold = zeros(statenum,1); % index for the policy function
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
    disp(['For iter ' num2str(vfct) ', VF error = ' num2str(solerr) ' & policy error = ' num2str(polerr)])

    %if VF has converged, break loop
    if (polerr<soltol)
        break
    end

    %otherwise, update VF and continue
    Vold = V;
        
    polindold=polind;
    polold = pol;
    
    %apply Howard improvement to skip some iterations
    %for each vfct step, do howardnum pushforwards of the current policy
    for howct=1:howardnum

        %based on guess for value function, form continuation value EVmat
        %with (i,j) entry equal to the EV(') at state i given policy j
        EVmat = kron((reshape(Vold,pnum,mnum)*(pr_mat_m'))',ones(pnum,1));

        %form RHS matrix with static and continutation payoff for state i and
        %policy j
        RHSmat = Rmat + beta*EVmat;

        %extract the appropriate elements of RHSmat implied by polindold
        polind_indexes = sub2ind(size(RHSmat),(1:statenum)',polindold);

        %update value function
        V = RHSmat(polind_indexes);

        %now, go to next step assuming this is the starting point for Vold
        Vold = V;

    end
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
