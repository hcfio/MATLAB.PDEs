% lkdv.m
%%%%%%%%%%%%%%%%%%%%%%
% initial data
N=256;
M=1000;
u=zeros(N,N);
for p=33:96
    u(p,1)=1/4;
    u(p+128,1)=1/8;
end

I=eye(N);
t=1/100;

% sqrt(M) times Fourier matrix
omega=exp(i*2*pi/N);
F=zeros(N);
for k=1:N
    for l=1:N
    F(k,l)=power(omega,-(k-1)*(l-1));
    end
end

D1=zeros(N,N);
D3=zeros(N,N);
for k=1:N/2
    D1(k,k)=k-1;
    D1(k+N/2,k+N/2)=k-N/2;
    D3(k,k)=power(k-1,3);
    D3(k+N/2,k+N/2)=power(k-N/2,3);
end

G=F'*(i*8*power(pi,3)*t*power(N,-3)*D3)*F/N;
A=I+G+G^2/2+G^3/6+G^4/24;


for l=2:M
    u(:,l)=A*u(:,l-1); 
end

v=zeros(M,N);
for l=1:M
    v(l,:)=u(:,M-l+1);
end

% plotting
subplot(1,2,1)
imagesc(real(v));
colorbar;
colormap(jet(50))
yticks([1 1000])
yticklabels({'10','0'})
ylabel('time t')
xticks([1 256])
xticklabels({'0','256'})
xlabel('position x')
title('IVP for the linearlized KdV equation')

subplot(1,2,2)
mesh(real(v'));
view(100,30)
xticks([1 1000])
xticklabels({'10','0'})
xlabel('time t')
yticks([1 256])
yticklabels({'0','256'})
ylabel('position x')
title('IVP for the linearlized KdV equation')
