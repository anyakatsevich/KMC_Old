L_loc = 200;
np = 2;
L = L_loc*np;
hInit = zeros(1, L);
hFinal = zeros(1, L);

for i = 0:np-1
    htemp = load(strcat('hInit0',num2str(i),'.txt'));
    hInit((L_loc*i+1):(L_loc*i + L_loc)) = htemp;
    htemp = load(strcat('hFinal0',num2str(i),'.txt'));
    hFinal((L_loc*i+1):(L_loc*i + L_loc)) = htemp;
end

figure(1); 
plot(hFinal/400); 
hold on; 
plot(hInit/400); 
legend('Final', 'Init');
hold off;


figure(2);
plot(hInit - hFinal);