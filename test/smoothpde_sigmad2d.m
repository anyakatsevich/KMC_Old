function dh = smoothpde_sigmad2d(t,h,Lap2d,Dp2d1,Dp2d2,Dm2d1,Dm2d2,K,p)

Dh1 = Dp2d1*h;
Dh2 = Dp2d2*h;



v1 = p*K*sign(Dh1).*(abs(Dh1).^(p-1));
v2 = p*K*sign(Dh2).*(abs(Dh2).^(p-1));

dh = -0.25*Lap2d*(Dm2d1*v1+Dm2d2*v2);

