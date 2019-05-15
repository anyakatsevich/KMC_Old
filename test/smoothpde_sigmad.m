function dh = smoothpde_sigmad(t,h,Lap,Dp,Dm,K,p)


Dh = Dp*h;

v = p*K*sign(Dh).*(abs(Dh).^(p-1));

dh = -0.5*Lap*Dm*v;

