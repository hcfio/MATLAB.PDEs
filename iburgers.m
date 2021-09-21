% iburgers.m
%%%%%%%%%%%%%%%%%%%%%%
% initial data
N=256;
M=1000;
u=zeros(N,M);
v=zeros(N,M);
dv=zeros(N,M);
z=zeros(2*N,M);
w=zeros(N,M);
for p=1:32
    u(p,1)=5*sin(pi*p/64);
end
for p=33:192
    u(p,1)=5*cos(pi*(p-32)/320);
end

I=eye(N);
t=1/200;

% sqrt(M) times Fourier matrix
omega=exp(i*2*pi/N);
F=zeros(N);
for k=1:N
    for l=1:N
    F(k,l)=power(omega,-(k-1)*(l-1));
    end
end


D1=zeros(N,N);
for k=1:N/2
    D1(k,k)=k-1;
    D1(k+N/2,k+N/2)=k-N/2;
end

v(:,1)=F*u(:,1);
dv(:,1)=2*pi*i*D1*v(:,1)/N;
z(:,1)=[v(:,1);v(:,1)];
for k=1:N
    for q=1:N
        w(k,1)=w(k,1)+z(N+k-q+1,1)*dv(q,1)/N;
    end
end



for l=2:M
    v(:,l)=v(:,l-1)-t*w(:,l-1);
    dv(:,l)=2*pi*i*D1*v(:,l)/N;
    z(:,l)=[v(:,l);v(:,l)];
    for k=1:N
        for q=1:N
            w(k,l)=w(k,l)+z(N+k-q+1,l)*dv(q,l)/N;
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
yticklabels({'5','0'})
ylabel('time t')
xticks([1 256])
xticklabels({'0','256'})
xlabel('position x')
title('IVP for the inviscid Burgers equation')

subplot(1,2,2)
mesh(real(V'));
view(100,30)
xticks([1 1000])
xticklabels({'5','0'})
xlabel('time t')
yticks([1 256])
yticklabels({'0','256'})
ylabel('position x')
title('IVP for the inviscid Burgers equation')
