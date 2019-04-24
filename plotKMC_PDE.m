A = load('parameters.txt');
isdH = A(1); L = A(2); K = A(3); nsmpls = A(4); numTimes = A(5);
times = zeros(1, numTimes);
for i = 1:numTimes
    times(i) = A(5+i);
end

hKMC = load('h.txt');
dx = 1/L;
x = [0:dx:1-dx]';
figure; hold;
for i = 1:numTimes+1
    plot(x, hKMC(i,:));
end


%%

e = ones(L,1);

load lambdaofu_2

Lap = spdiags([e e -2*e e e], [-L+1, -1, 0, 1, L-1], L, L)/(dx*dx);

Dp = spdiags([e -e e], [-L+1, 0, 1], L, L)/(dx);
Dm = spdiags([-e e -e], [-1, 0, L-1], L, L)/(dx);

h0 = sin(2*pi*x);
p=2;
Tf = times(numTimes);
%pdehandle = @(t,h) smoothpde_anya(h,Lap,Dp,Dm,u,lambda,K);
%pdehandle = @(t,h)anyapde(t,h,Lap,Dp,Dm,K);
%pdehandle = @(t,h)smoothpde(t,h,Lap,Dp,Dm,u,K*lambda);
pdehandle2 = @(t,h)smoothpde_sigmad(t,h,Lap,Dp,Dm,K,p);
%pdehandle2 = @(t,h)smoothpde_Langevin(t,h,Lap,Dp,Dm,u,K*lambda);
%pdehandle = @(t,h)smoothpde_fsp(t,h,Lap,Dp,Dm,u,K*lambda);
%pdehandle2 = @(t,h)smoothpde_sigmadp1(t,h,Lap,Dp,Dm,K,p);
%pdehandle = @(t,h)roughpde(t,h,Dp,Dm,K,p);
%jachandle = @(t,h)roughpdejac(t,h,Dp,Dm,K,p);



options = odeset('RelTol', 1e-6,'AbsTol',1e-6); %,'Jacobian',jachandle);

tic
[T,Y] = ode15s(pdehandle2,[0 Tf],h0,options);
toc
h = Y(length(T),:);


figure;
plot(x, h, x, hKMC(numTimes+1,:));



