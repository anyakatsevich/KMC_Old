%clear all

Ti = 1E-4;

N = 600;
dx = 1/N;

p = 2;

K = 1.5; %0.5;

x1 = [0:dx:1-dx]';
N1 = length(x1);

e = ones(N1,1);

load lambdaofu_2

Lap = spdiags([e e -2*e e e], [-N1+1, -1, 0, 1, N1-1], N1, N1)/(dx*dx);

Dp = spdiags([e -e e], [-N1+1, 0, 1], N1, N1)/(dx);
Dm = spdiags([-e e -e], [-1, 0, N1-1], N1, N1)/(dx);

h0 = sin(2*pi*x1);

%h0 = max(sin(2*pi*x),zeros(N1,1));
%h0 = sin(pi*x1);

% h0 = 0*sin(2*pi*x1);
%  for j = 1:N1
%      if j>1 &&  j < round(N1/2)
%       h0(j) = exp(8)*exp(-(1./x1(j)) - 1./(1/2-x1(j)));
%      end
%  end
% h0 = .1*h0;
%h0 = exp(-(x1-.5).^2/.1);

%pdehandle = @(t,h) smoothpde_anya(h,Lap,Dp,Dm,u,lambda,K);
%pdehandle = @(t,h)anyapde(t,h,Lap,Dp,Dm,K);
%pdehandle = @(t,h)smoothpde(t,h,Lap,Dp,Dm,u,K*lambda);
pdehandle2 = @(t,h)smoothpde_sigmad(t,h,Lap,Dp,Dm,K,p);
%pdehandle2 = @(t,h)smoothpde_Langevin(t,h,Lap,Dp,Dm,u,K*lambda);
%pdehandle = @(t,h)smoothpde_fsp(t,h,Lap,Dp,Dm,u,K*lambda);
%pdehandle2 = @(t,h)smoothpde_sigmadp1(t,h,Lap,Dp,Dm,K,p);
%pdehandle = @(t,h)roughpde(t,h,Dp,Dm,K,p);
%jachandle = @(t,h)roughpdejac(t,h,Dp,Dm,K,p);


% Dh = Dp*h0;
% 
% v = sign(Dh).*interp1(u,lambda,abs(Dh));
% 
% dh = -0.5*Lap*(Dm*v);
% 
% plot(x,Dm*v,x,2*Lap*h0);

% eps = 1e-12;
% 
% pert = randn(N,1);
% 
% dfdh1 = (pdehandle(0,h0+eps*pert)-pdehandle(0,h0))/eps
% 
% dfdh2 = jachandle(0,h0);
% 
% dfdh2*pert
% 
% stop;

options = odeset('RelTol', 1e-6,'AbsTol',1e-6); %,'Jacobian',jachandle);

% h = h0;
% dt = 1e-12;
% for i = 1:1e4
%     f = pdehandle(0,h);
%     Df = jachandle(0,h);
%     z = (speye(N)-0.5*dt*Df)\f;
%     h = h + dt*z;
% end


tic
[T,Y] = ode15s(pdehandle2,[0 Ti],h0,options);
toc
h = Y(length(T),:);
figure; plot(x1, h - h0', 'r', x1, hFinal/L - h0');
% N=200;
% x1 = 0:1/N:1-1/N;
%hKMC_new = load('hFinal_1.txt');
%hKMC_new2 = load('hFinal.txt');
%hKMC_old = load('h.txt');
%hKMC_old = hKMC_old(2,:);
%hKMC_new = hKMC_new*N;
%hKMC_new2 = hKMC_new2*N;

%figure;
%plot(x1, h-h0', 'r',x1, hKMC_new - h0','b', x1, hKMC_new2 - h0', 'g',x1, hKMC_old - h0', 'k');
%plot(x1, hKMC_new - h0','b', x1, hKMC_new2 - h0', 'g');

  % Plot regular time slices
%  for j = 1:length(T)
%      h = Y(j,:)';
%      plot(x1,h,'-g','LineWidth',2)
%      xlim([0 1])
%      %ylim([-1 1])
%      xlabel('x');
%      ylabel('Crystal Height');
%      drawnow;
%  end
% title('p = 1,T = 1E-4 vs dH rates L = 200, 1 samples');
% z = load('C:\Users\Anya\C Projects\Crystals\p=2_L=200_T=5x1E-4_nsmpls =1.dat');
% M = 200;
% x = [0:1/M:1-1/M]';
% hold on
% plot(x, z, 'k-s', 'MarkerSize', 2);
%title('dH Rates, p = 1.5, T = 2x1e-4, L = 100, 4 samples');
% tic
% [T1,Y1] = ode15s(pdehandle,[0 Ti],h/(max(h)),options);
% toc
% h1 = Y1(length(T1),:);
% 
% tic
% [T2,Y2] = ode15s(pdehandle,[0 Ti],h1/(max(h1)),options);
% toc
% h2 = Y2(length(T2),:);
% 
% tic
% [T3,Y3] = ode15s(pdehandle,[0 Ti],h2/(max(h2)),options);
% toc
% h3 = Y3(length(T3),:);
% 
% tic
% [T4,Y4] = ode15s(pdehandle,[0 Ti],h3/(max(h3)),options);
% toc
% h4 = Y4(length(T4),:);
% 
% tic
% [T5,Y5] = ode15s(pdehandle,[0 Ti],h4/(max(h4)),options);
% toc
% h5 = Y5(length(T5),:);
% 
% plot(x(1:end),h3(1:end),'-sk',x(1:end),h4(1:end),'-*b',x(1:end),h5(1:end),'-^c','LineWidth',2);
% xlabel('$x$','interpreter','latex','FontSize',18);
% ylabel('$h$','interpreter','latex','FontSize',18);
% %legend('N = \infty', 'N = 400');
% title('(c)','FontSize',18);
% set(gca,'FontSize',18,'TickLength',[.02 0])

% kmc50 = load('kmcrough_p2_K1pt5_L50_T1em25_heights.txt');
% 
% figure; plot(x(1:end),h(1:end),'k',x,kmc50,'-.ok')
% axis([0,1,-1,1]);
% legend('N=\infty','N=50')


%,x(1:end),h3(1:end),'-sk',x(1:end),h4(1:end),'-*b',x(1:end),h5(1:end),'-^c','LineWidth',2);
%legend('h(T)','Evolution of h1 = h(T)/max(h) by T','Evolution of h2 = h2(T)/max(h2) by T' );

%title('p=1');
%exportfig(gcf,'self-similar_rough_p2_K1pt5_T1e-20it.eps','color')
% kmc50 = load('/Users/JLM/Dropbox/kmckrug/data/kmcrough_p1pt5_T1e-12_L50.dat');
% x50 = [0:1/50:1-1/50];
% % kmc100 = load('/Users/JLM/Dropbox/kmckrug/data/kmcrough_p1.2_T1e-1_L100.dat');
% % x100 = [0:1/100:1-1/100];
% % kmc200 = load('/Users/JLM/Dropbox/kmckrug/data/kmcrough_p1.2_T1e-1_L200.dat');
% % x200 = [0:1/200:1-1/200];
% % kmc400 = load('/Users/JLM/Dropbox/kmckrug/data/kmcrough_p1.2_T1e-1_L400.dat');
% % x400 = [0:1/400:1-1/400];
% 
% 
% plot(x,h,'k',x50,kmc50,'ko-.',x100,kmc100,'ks-.',x200,kmc200,'kd-.',x400,kmc400,'k^-.');
% legend('N = \infty','N = 50','N = 100','N = 200','N = 400');
% % plot(x,h,'k',x400,kmc400,'k-.');
% % legend('N = \infty','N = 400');
% axis([0,1,-1,1]);
% xlabel('x');
% ylabel('height');
%kmcrough_p2_K1pt5_L50_T1em25_heights.txt

