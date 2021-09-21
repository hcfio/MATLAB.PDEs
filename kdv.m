% kdv.m
%%%%%%%%%%%%%%%%%%%%%%
% initial data
N=256;
M=1000;
u=zeros(N,M);
v=zeros(N,M);
z=zeros(2*N,M);
w=zeros(N,M);
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
G=i*8*power(pi,3)*t*power(N,-3)*D3;
A=I+G+G^2/2+G^3/6+G^4/24;
B=-3*i*2*pi*t*power(N,-1)*(I+G/2+G^2/6+G^3/24+G^4/120)*D1;

v(:,1)=F*u(:,1);
z(:,1)=[v(:,1);v(:,1)];
for k=1:N
    for q=1:N
        w(k,1)=w(k,1)+z(N+k-q+1,1)*z(q,1)/N;
    end
end



for l=2:M
    v(:,l)=A*v(:,l-1)+B*w(:,l-1);
    z(:,l)=[v(:,l);v(:,l)];
    for k=1:N
        for q=1:N
            w(k,l)=w(k,l)+z(N+k-q+1,l)*z(q,l)/N;
        end
    end
end

U=F'*v/N;
V=zeros(M,N);
for l=1:M
    V(l,:)=U(:,M-l+1);
end

% plotting
subplot(1,2,1)
imagesc(real(V));
colorbar;
colormap(jet(50))
yticks([1 1000])
yticklabels({'10','0'})
ylabel('time t')
xticks([1 256])
xticklabels({'0','256'})
xlabel('position x')
title('IVP for the KdV equation')

subplot(1,2,2)
mesh(real(V'));
view(100,30)
xticks([1 1000])
xticklabels({'10','0'})
xlabel('time t')
yticks([1 256])
yticklabels({'0','256'})
ylabel('position x')
title('IVP for the KdV equation')