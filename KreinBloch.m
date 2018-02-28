% Krein matrix for 5th-order KdV (water wave problem): Bloch-wave spectra 

clear all
% pack

global L N x b c ks
global D1 D2 D4

TIME = clock;

b = -8/15; c = b + 2/15;
% mu = 0.3585;                   % Bloch wave decomposition parameter- -1/2<mu<1/2
% mu = 1 - .2*sqrt(10);
mu = 0.20;

% parameters for Krein matrix calculation

zm = -0.21; zp = 0.21;
Nz = floor( (zp-zm)*10000 );
z = linspace(zm,zp,Nz+1);      % interval over which to evaluation matrix
zsize = length(z);

% the differentiation matrices - Fourier spectral

N = 96;                     % number of grid points
L = pi;                      % spatial interval is [-L,+L]
dx = 2*L/N;                  % spatial step 
x = linspace(-L+dx,L,N)';    % discretization of interval 
[~,D1temp] = fourdif(N,1);
D1 = (pi/L)*D1temp;
[~,D2temp] = fourdif(N,2);
D2 = (pi/L)^2*D2temp;
[~,D3temp] = fourdif(N,3);
D3 = (pi/L)^3*D3temp;
[~,D4temp] = fourdif(N,4);
D4 = (pi/L)^4*D4temp;

% compute the steady-state solution

kstart = 1 - 0.01;     % period increases as magnitude increases (p=3,+/-) - stable
kstop = 1 - 0.005;      % small solution (p=3, +/-) 
kstep = (kstop - kstart)/12;
q0 = 0.4*cos(x);  q1 = q0;

options = optimset('Display','off','MaxFunEvals',1e6,'MaxIter',500,'TolFun',1e-5);
kint=kstart:kstep:kstop;

for ks=kint;
    [q,fval,exitflag,output] = fsolve(@KDVsteady,q1,options); 
    q1 = q;   
end
k = ks;

figure(3)
plot(x,q0,x,q1)

% the linear operators


q = 0*x; k = 1;
A1t = D1 + 1i*mu*eye(N); % A1 = (A1-A1')/2;
A1 = inv(A1t);  A1 = (A1-A1')/2;
L1 = 2/15*k^2*(D4 + 1i*4*mu*D3 - 6*mu^2*D2 - 1i*4*mu^3*D1 + mu^4*eye(N));
L2 = k*(-b*eye(N) + diag(q))*(D2 + 1i*2*mu*D1 - mu^2*eye(N));
L3 = -c*eye(N) + 3*eye(N)*diag(q) + k*diag(D2*q);
L4 = k*diag(D1*q)*(D1 + 1i*mu*eye(N));
A0 = L1 + L2 + L3 + L4; A0 = (A0+A0')/2;

% construct negative subspace and projections

[VA0,evalsA0] = eig(A0);
evalsA0 = diag(evalsA0);
[evalsA0,i] = sort(evalsA0);
numA0 = 0;
for i=1:N
    if (abs(evalsA0(i))<0.001)
        numA0 = numA0;
    elseif evalsA0(i)<0
        numA0 = numA0+1;
    end
end   

% block matrices

evalsA0abs = abs(evalsA0);
D0m12 = diag(sqrt(1./evalsA0abs));
A1new = D0m12*VA0'*A1*VA0*D0m12;
A0m = -eye(numA0);
A0p = eye(N-numA0);
A1m = A1new(1:numA0,1:numA0);
A1p = A1new(numA0+1:N,numA0+1:N);
A1L = A1new(numA0+1:N,1:numA0);

% Krein eigenvalues

Kreineval = zeros(numA0,zsize);
for kk=1:zsize
    Ks = zeros(numA0);  
    Ks = (A0p + z(kk)*1i*A1p)\A1L;
    Ks = A1L'*Ks;
    Ks = A0m + z(kk)*1i*A1m  - z(kk)^2*Ks;
%     Ks = -z(kk)*(Ks' + Ks)/2;
    Kreineval(:,kk) = eig(Ks);
end  

poles = polyeig(eye(N-numA0),1i*A1p); poles = real(poles);
stareig = polyeig(A0,1i*A1);
falsezero = 0;

% plot the Krein eigenvalues

ymin = -0.35; ymax = 0.35;

figure(1)
plot(z, Kreineval ,'b.')
hold on
for nn=1:N-numA0
    line([poles(nn) poles(nn)],[ymin ymax],'Color','g','LineWidth',1)
end    
hold on
plot(z,0*z,'k-','Linewidth',1)
%plot(real(poles),imag(poles),'rx','Linewidth',1,'MarkerSize',8)
hold on
plot(real(stareig),imag(stareig),'ro','MarkerSize',6,'Linewidth',1)
hold on
plot(real(falsezero),imag(falsezero),'rx','MarkerSize',8,'Linewidth',1)
hold off
%grid on
%xlabel z, ylabel r(z)
% axis([zmin zmax ymin ymax]), set(gca,'xtick',zmin:3.0:zmax,'ytick',ymin:10:ymax)
% axis([zm zp ymin ymax]), set(gca,'xtick',zm:0.1:zp,'ytick',-100:0.1:100)
axis([zm zp -0.18 0.18]), set(gca,'xtick',-3:0.1:3,'ytick',-100:0.1:100) % mu=0.2
% hold on
% %plot(lambda, klambda2,'r.')
% plot([xmin xmax],[0 0])
% plot([0 0],[ymin ymax])
% title('K(\lambda)','fontweight','bold','fontsize',12,'fontname','Times New Roman')
% axis([xmin xmax ymin ymax]); 
xlabel('z','fontsize',10,'fontname','Times New Roman') 
% ylabel('r(z)','fontsize',10,'fontname','Times New Roman')

figure(2)
subplot(2,1,1)
plot(z, Kreineval ,'b.')
hold on
for nn=1:N-numA0
    line([poles(nn) poles(nn)],[ymin ymax],'Color','g','LineWidth',1)
end    
hold on
plot(z,0*z,'k-','Linewidth',1)
%plot(real(poles),imag(poles),'rx','Linewidth',1,'MarkerSize',8)
hold on
plot(real(stareig),imag(stareig),'ro','MarkerSize',6,'Linewidth',1)
hold on
plot(real(falsezero),imag(falsezero),'rx','MarkerSize',8,'Linewidth',1)
hold off
%grid on
%xlabel z, ylabel r(z)
% axis([zmin zmax ymin ymax]), set(gca,'xtick',zmin:3.0:zmax,'ytick',ymin:10:ymax)
% axis([zm zp ymin ymax]), set(gca,'xtick',zm:0.1:zp,'ytick',-100:0.1:100)
axis([zm zp -0.18 0.18]), set(gca,'xtick',-3:0.1:3,'ytick',-100:0.1:100) % mu=0.2
% hold on
% %plot(lambda, klambda2,'r.')
% plot([xmin xmax],[0 0])
% plot([0 0],[ymin ymax])
% title('K(\lambda)','fontweight','bold','fontsize',12,'fontname','Times New Roman')
% axis([xmin xmax ymin ymax]); 
xlabel('z','fontsize',10,'fontname','Times New Roman') 
% ylabel('r(z)','fontsize',10,'fontname','Times New Roman')
subplot(2,1,2)
plot(z, Kreineval ,'b.')
hold on
for nn=1:N-numA0
    line([poles(nn) poles(nn)],[ymin ymax],'Color','g','LineWidth',1)
end    
hold on
plot(z,0*z,'k-','Linewidth',1)
%plot(real(poles),imag(poles),'rx','Linewidth',1,'MarkerSize',8)
hold on
plot(real(stareig),imag(stareig),'ro','MarkerSize',6,'Linewidth',1)
hold off
%grid on
%xlabel z, ylabel r(z)
% axis([zmin zmax ymin ymax]), set(gca,'xtick',zmin:3.0:zmax,'ytick',ymin:10:ymax)
% axis([zm zp ymin ymax]), set(gca,'xtick',zm:0.1:zp,'ytick',-100:0.1:100)
% axis([-0.117 -0.113 -0.006 0.006]), set(gca,'xtick',-3:0.01:3,'ytick',-100:0.005:100) % mu=0.36
axis([0.085 0.135 -0.026 0.026]), set(gca,'xtick',-3:0.01:3,'ytick',-100:0.01:100) % mu=0.36
% hold on
% %plot(lambda, klambda2,'r.')
% plot([xmin xmax],[0 0])
% plot([0 0],[ymin ymax])
% title('K(\lambda)','fontweight','bold','fontsize',12,'fontname','Times New Roman')
% axis([xmin xmax ymin ymax]); 
xlabel('z','fontsize',10,'fontname','Times New Roman')

% figure(2)
% plot(real(stareig),imag(stareig),'bo','MarkerSize',8)
% hold on
% plot(z,0*z,'k-','Linewidth',1)
% hold off
% %grid on
% axis([zm zp ymin ymax]), set(gca,'xtick',zm:0.1:zp,'ytick',-100:0.1:100)
% xlabel('Re(z)','fontsize',12,'fontname','Times New Roman')
% ylabel('Im(z)','fontsize',12,'fontname','Times New Roman')

TOTtime = etime(clock,TIME);
['total time = ', num2str(TOTtime/60,'%8.3f'), ' minutes']