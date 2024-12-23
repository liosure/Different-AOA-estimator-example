clear;close;
a = @(theta,K) exp(1j*2*pi*(0:K-1)'*theta);
L = 100;  % symbols
theta = [-5,10];
f = sin(theta/180*pi)/2;  % f = cos(theta)
scale = 1024;
SNR = 1e-1;
p_dft = zeros(4,scale);
i=1;
for K = [4,8,16,32]*2 % number of antennas
A = a(f,K);  % A = [a(f1), ..., a(fK)]
Data = randi([0,3],2,L);
s = 1/sqrt(2)*qammod(Data,4,'gray');  % 3 sources
n = sqrt(SNR/2)*(randn(K,L)+1j*randn(K,L));
y = A*s+n;
p_dft(i,:) = sum(abs(fft(y,scale)),2);
p_dft(i,:) = p_dft(i,:)./max(p_dft(i,:));
eval(['label',num2str(i),' = [',char("'K=',"),'num2str(K)];']);
i = i+1;
end

f = figure(1);
axes1 = axes('Parent',f);

hold(axes1,'on')
for i = 1:4
plot(axes1,asin(2*(-0.5:1/1024:0.5-1/1024))/pi*180,...
    [p_dft(i,scale/2:end),p_dft(i,1:scale/2-1)]);
end
stem(theta,[1,1],'ro')
axis([-90,90,0,1])
legend(label1,label2,label3,label4,'真实值')
xlabel('入射角/°')
ylabel('归一化功率')
grid on

axes2 =axes('Parent',f,...
    'Position',[0.17 0.59 0.24 0.28]);
hold(axes2,'on')
for i = 1:4
plot(axes2,asin(2*(-0.5:1/1024:0.5-1/1024))/pi*180,...
    [p_dft(i,scale/2:end),p_dft(i,1:scale/2-1)]);
end
stem(theta,[1,1],'ro')
axis([8,12,0.9,1.05])
xlabel('入射角/°')
ylabel('归一化功率')
grid on