clear all

gauss = @(x,a,b,c) a*exp(-(((x-b).^2)/(2*c.^2)));

% time line
x = 0:.001:20;
lengthX = length(x);
% gaussian params
amp = 1; 
center = 0.5; 
sig = 0.05;
r = 1 + 20*rand(30,1);
% train of gaussian pulses
y = gauss(x, amp, r, sig);

% overlapping signals
ys = sum(y);


% Samplings
nm = 10;

% m0: It is the signal that we expected to get with n samples
step = fix(1000/nm);
m = movsum(ys, [0,step-1]);
m0 = zeros(size(x));
ini = 1;
while ini < lengthX-step
    m0(ini:ini+step) = m(ini);
    ini = ini + step;
end

% Generamos nm muestreos
m = movsum(ys, [0,999]);
M = zeros(nm, lengthX+1000);

for i=0:19
    for j=1:nm
        ini = 1+i*1000+step*(j-1);
        M(j,ini:ini+1000) = m(ini);
    end
end
M = M(:,1:lengthX);

%%
% Original signal plot and nm samples 
figure
% train of gaussian pulses
subplot(nm+1,1,1)
plot(x,y)
% xlabel('time(ms)')
title('Original signals and samplings')

for j=1:nm
    subplot(nm+1,1,j+1)
    plot(x, M(j,:))
end


%% 
figure;
subplot(3,1,1)
plot(x,y)
title('Original signals - train of gaussian pulses')

subplot(3,1,2)
plot(x,ys)
title('Overlapping original signals')

subplot(3,1,3)
plot(x,m0)
title(['Sampling with a wide window ' num2str(step/1000) ' (Expected to get)'])
