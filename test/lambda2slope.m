% %%%%%% u(lambda)
% 
K = 0.5;

Nx = 100000;

Ns = 1000;

R = 100;

lambda = linspace(0,R,Nx);

d = potential(0)*ones(1,Nx);
smin = zeros(1,Nx);
for s = 1:Ns;
    dn = potential(s)-lambda*s;
    J = find(dn<d);
    d(J) = dn(J);
    smin(J) = s*ones(1,length(J));
end

z = exp(-K*(potential(0)-d));
u = zeros(1,Nx);
for s = 1:Ns;
    u = u + s*(exp(-K*(potential(s)-lambda*s-d))...
        -exp(-K*(potential(s)+lambda*s-d)));
    z = z + exp(-K*(potential(s)-lambda*s-d))...
        +exp(-K*(potential(s)+lambda*s-d));
end
u = u./z;

u = u';
lambda = lambda';

save lambdaofu_1 u lambda 
C = 4*pi*exp(-pi^2/K)/K;
plot(u,lambda-2*u - C*sin(2*pi*u));
% % figure(1)
% % subplot(2,1,1)
% % plot(lambda,-d+log(z)/K)
% % subplot(2,1,2)
% % plot(lambda,u);
% % hold on
% % plot(lambda(1:Nx-1),(log(z(2:Nx))-log(z(1:Nx-1)))*Nx/K)
% % hold off
% 
% % subplot(2,1,2)
% % plot(lambda,sinh(lambda)./(cosh(K) - cosh(lambda)));
% 
% U = linspace(0,u(Nx),Nx);
% 
% Lambda = interp1(u',lambda', U');
% 
% figure(1)
% plot(U,Lambda);
