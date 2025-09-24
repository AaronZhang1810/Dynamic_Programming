%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Firm price adjustment problem
% Implementing nonstochastic simulation
% Firm_Pricing_Simulate_Dist.m
%
% S. Terry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('%%% Implement nonstochastic simulation')
disp(' ')
tic

%%%%%%%%%initialize distribution
distold = zeros(pnum,mnum);
distold(:)=1;
distold = distold/sum(distold(:));

dist = 0*distold;

for distct=1:maxit

    for pmin1ct=1:pnum
        for mct=1:mnum
            if (distold(pmin1ct,mct)>0)
                %extract policy
                pct = polind(pmin1ct,mct);

                %push forward the distribution
                dist(pct,:) = dist(pct,:) + pr_mat_m(mct,:)*distold(pmin1ct,mct);
            end
        end
    end

    %compute error
    disterr = max(abs(dist(:)-distold(:)));

    %display some diagnostics
    if (mod(distct,5)==1)
        disp(['For iter ' num2str(distct) ', dist error = ' num2str(disterr)])
    end

    %if dist has converged, break loop
    if (disterr<soltol)
        break
    end


    distold = dist/sum(dist(:));
    dist = 0*distold;

end


%%%%%%%%% now compute marginal distributions
pdist = sum(dist,2);
mdist = sum(dist,1)';

%%%%%%%%%% some simple diagnostics
disp(' ')
disp('Quick distribution bounds check: wgt on lb, ub')
disp(['p: ' num2str([pdist(1) pdist(pnum)])])
disp(['m: ' num2str([mdist(1) mdist(mnum)])])
disp(' ')


%%%%%%%%%% plot some distributions

figure;
plot(pgrid,pdist,'b',...
    'LineWidth',lwidnum)
xlabel('p')
ylabel('Weight')
title('Ergodic Distribution for Firm Price')
set(gca,'FontSize',fsizenum)

figure;
plot(mgrid,mdist,'b',...
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


toc
disp(' ')
