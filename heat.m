% heat.m
%%%%%%%%%%%%%%%%%%%%%%
% initial data
N=256;
M=1000;
u=zeros(N,M);
for p=65:192
    u(p,1)=1;
end

I=eye(N);
t=1/10;

% sqrt(M) times Fourier matrix
omega=exp(i*2*pi/N);
F=zeros(N);
for k=1:N
    for l=1:N
    F(k,l)=power(omega,-(k-1).*(l-1));
    end
end

w=zeros(N,1);
for k=1:N/2
    w(k,1)=power(k-1,2);
    w(k+N/2,1)=power(N/2-k,2);
end
D=diag(w);

G=-(4*power(pi,2)*t*power(N,-3))*F'*D*F;
A=I+G+G^2/2;

for n=2:M
    u(:,n)=A*u(:,n-1);
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
yticklabels({'100','0'})
ylabel('time')
xticks([1 256])
xticklabels({'0','256'})
xlabel('x')
title('heat flow')

subplot(1,2,2)
mesh(real(v'));
view(120,30)
xticks([1 1000])
xticklabels({'100','0'})
xlabel('time')
yticks([1 256])
yticklabels({'0','256'})
ylabel('x')
title('heat flow')