function dh = smoothpde2d(t,h,Lap2d,Dp2d1,Dp2d2,Dm2d1,Dm2d2,u,lambda)

Dh1 = Dp2d1*h;
Dh2 = Dp2d2*h;

v1 = sign(Dh1).*interp1q(u,lambda,abs(Dh1));
v2 = sign(Dh2).*interp1q(u,lambda,abs(Dh2));
x
dh = -0.25*Lap2d*(Dm2d1*v1+Dm2d2*v2);

% dh = -Lap*(Lap*h);