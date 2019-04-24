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





