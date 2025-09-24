// Stochastic case model (linearized version)

var c k a;

varexo eps;

parameters alpha delta beta gamma rho;
alpha = 0.33;
beta = 0.95;
delta = 0.1;
gamma = 5;
rho = 0.95;


model(linear);
# kss = ((1 - beta*(1 - delta))/(alpha*beta))^(1/(alpha-1));
# css = kss^alpha - delta*kss;
# Rss = alpha*kss^(alpha - 1);

% eq(1)
c(+1) - beta*Rss/gamma*a(+1) -beta*(alpha-1)*Rss/gamma*k = c;

% eq(2)
k = kss^(alpha-1)*a - css/kss*c + 1/beta*k(-1);

% eq(3)
a = rho*a(-1) + eps;

end;

steady;
check;

shocks;
var eps;
stderr 0.01;
end;

stoch_simul(irf=10, order=1);

