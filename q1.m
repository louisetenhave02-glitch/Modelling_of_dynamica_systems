clear;
close all;
clc;
p=params();
C2K=273.15;
if ~isfolder(fullfile(pwd,'figures'))
    mkdir(fullfile(pwd,'figures'));
end
set(groot,'defaultAxesFontSize',11,'defaultLineLineWidth',1.3,'defaultAxesXGrid','on','defaultAxesYGrid','on');
set(groot,'defaultTextInterpreter','latex','defaultAxesTickLabelInterpreter','latex','defaultLegendInterpreter','latex');

ellipse=struct('xc',0.11,'yc',0.11,'ax',0.06,'ay',0.04,'T',8);
heart=struct('xc',0.11,'yc',0.12,'scale',0.10,'T',8);

th=[0.3;-1.2];
dth=[2;5];
J=p.motor.r_m*p.K;
th2=J\kin(th,p);
dth2=J\kin(dth,p);
show('Q1.1 roundtrip error theta',max(abs([th-th2;dth-dth2])),'rad (rad/s)','%.2e');
dl=kin([1;0],p);
dr=kin([0;1],p);
show('Q1.1 [dx dy] per rad theta_l',dl(1)*1e3,sprintf('mm  (dy = %.3f mm)',dl(2)*1e3));
show('Q1.1 [dx dy] per rad theta_r',dr(1)*1e3,sprintf('mm  (dy = %.3f mm)',dr(2)*1e3));

A=zeros(6);
B=zeros(6,2);
I6=eye(6);
I2=eye(2);
for k=1:6
    A(:,k)=f_full(I6(:,k),[0;0],p);
end
for k=1:2
    B(:,k)=f_full(zeros(6,1),I2(:,k),p);
end
lam=eig(A);
disp('Q1.2 A =');
disp(A);
disp('Q1.2 B =');
disp(B);
disp('Q1.2 eigenvalues A =');
disp(lam);
show('Q1.2 tau_el = L_m/R_m',p.motor.L_m/p.motor.R_m*1e3,'ms');
lam=lam(abs(lam)>1e-9);
show('Q1.2 slowest pole (nonzero)',-max(real(lam)),'1/s');
show('Q1.2 fastest pole',-min(real(lam)),'1/s');

opts=odeset(p.opts,'MaxStep',1e-3);
[t,x]=ode45(@(t,x)f_full(x,[1;0],p),[0 0.2],zeros(6,1),opts);
pos=kin(x(:,1:2)',p);
vel=kin(x(:,3:4)',p);
fig=figure('Position',[100 100 800 650]);
subplot(2,2,1);
plot(t,x(:,3:4));
ylabel('$\omega$ [rad/s]');
xlabel('$t$ [s]');
legend('$\omega_l$','$\omega_r$','Location','best');
subplot(2,2,2);
plot(t,x(:,5:6));
ylabel('$i_m$ [A]');
xlabel('$t$ [s]');
legend('$i_l$','$i_r$','Location','best');
subplot(2,2,3);
plot(t,vel'*1e3);
ylabel('velocity [mm/s]');
xlabel('$t$ [s]');
legend('$\dot x$','$\dot y$','Location','best');
subplot(2,2,4);
plot(pos(1,:)*1e3,pos(2,:)*1e3);
axis equal;
xlabel('$x$ [mm]');
ylabel('$y$ [mm]');
title('step $u_l = 1$ V, $u_r = 0$');
save_plot(fig,'Q1_2_step_response');
show('Q1.2 final velocity xdot (u_l = 1 V)',vel(1,end)*1e3,'mm/s');
show('Q1.2 final velocity ydot (u_l = 1 V)',vel(2,end)*1e3,'mm/s');

u=@(t)[2*sin(15*t);1.5*cos(9*t)-0.5];
us=@(t)p.K*u(t);
tg=linspace(0,1,1001);
[~,xf]=ode45(@(t,x)f_full(x,u(t),p),tg,zeros(6,1),p.opts);
[~,xx]=ode45(@(t,z)f_axis(z,[1 0]*us(t),1,p),tg,zeros(3,1),p.opts);
[~,xy]=ode45(@(t,z)f_axis(z,[0 1]*us(t),2,p),tg,zeros(3,1),p.opts);

pos=kin(xf(:,1:2)',p);
is=p.K*xf(:,5:6)';
err=[pos(1,:)'-xx(:,1),is(1,:)'-xx(:,3),pos(2,:)'-xy(:,1)];
show('Q1.3 max |x_full - x_dec|',max(abs(err(:,1)))*1e3,'mm','%.2e');
show('Q1.3 max |y_full - y_dec|',max(abs(err(:,3)))*1e3,'mm','%.2e');
show('Q1.3 max |i_x,full - i_x,dec|',max(abs(err(:,2))),'A','%.2e');
show('Q1.3 max |x| range',max(abs(xx(:,1)))*1e3,'mm');

fig=figure('Position',[100 100 800 500]);
subplot(2,1,1);
plot(tg,pos(1,:)*1e3,tg,xx(:,1)*1e3,'--',tg,pos(2,:)*1e3,tg,xy(:,1)*1e3,'--');
ylabel('position [mm]');
legend('$x$ full','$x$ decoupled','$y$ full','$y$ decoupled','Location','best');
subplot(2,1,2);
plot(tg,err(:,1)*1e3,tg,err(:,3)*1e3);
ylabel('difference [mm]');
xlabel('$t$ [s]');
legend('$x$','$y$','Location','best');
save_plot(fig,'Q1_3_decoupling_check');

shapes={'ellipse','heart'};
paths={ellipse,heart};
for k=1:2
    sh=shapes{k};
    sp=paths{k};
    t=linspace(0,sp.T,4001);
    r=ref_path(sh,t,sp);
    [ux,ix]=feedforward(r,p);
    u=p.K\ux;
    x0=[J\r.p(:,1);J\r.v(:,1);p.K\ix(:,1)];
    u_fun=@(t)p.K\feedforward(ref_path(sh,t,sp),p);
    opts=odeset('RelTol',1e-6,'AbsTol',1e-9,'MaxStep',sp.T/1000);
    [t,x]=ode45(@(t,x)f_full(x,u_fun(t),p),t,x0,opts);
    pos=kin(x(:,1:2)',p);
    if k==1
        ell_t=t;
        ell_pos=pos;
    end

    fig=figure('Position',[100 100 900 380]);
    subplot(1,2,1);
    plot(r.p(1,:)*1e3,r.p(2,:)*1e3,'k--',pos(1,:)*1e3,pos(2,:)*1e3);
    axis equal;
    xlabel('$x$ [mm]');
    ylabel('$y$ [mm]');
    legend('reference','simulation','Location','best');
    subplot(1,2,2);
    plot(t,u');
    xlabel('$t$ [s]');
    ylabel('$u$ [V]');
    legend('$u_l$','$u_r$','Location','best');
    save_plot(fig,['Q1_4_' sh]);
    show(sprintf('Q1.4 %s max tracking error',sh),max(vecnorm(pos-r.p))*1e6,'um','%.3g');
    show(sprintf('Q1.4 %s max |u|',sh),max(abs(u(:))),'V','%.3g');
    show(sprintf('Q1.4 %s max speed',sh),max(vecnorm(r.v))*1e3,'mm/s','%.3g');
end

n=p.noz;
xn0=[n.T_e;0;0];
[t1,x1]=ode15s(@(t,x)f_nozzle(x,n.u_print,p),[0 1000],xn0,p.opts_stiff);
[t2,x2]=ode45(@(t,x)f_nozzle(x,n.u_print,p),[0 1e-3],xn0,p.opts);
T_ss=n.T_e+n.R_tn*n.u_print^2/n.R_n;
show('Q1.5 u_print',n.u_print,'V');
show('Q1.5 T_ss analytical',T_ss-C2K,'degC');
show('Q1.5 T_n(1000 s) simulation',x1(end,1)-C2K,'degC');
show('Q1.5 thermal time constant R_tn*C_tn',n.R_tn*n.C_tn,'s');
show('Q1.5 LC damping ratio: (1/2R)sqrt(L/C)',sqrt(n.L_n/n.C_n)/(2*n.R_n),'-');

fig=figure('Position',[100 100 900 380]);
subplot(1,2,1);
plot(t1,x1(:,1)-C2K);
hold on;
plot([0 t1(end)],[T_ss T_ss]-C2K,'k--');
xlabel('$t$ [s]');
ylabel('$T_n$ [$^\circ$C]');
legend('simulation','$T_{ss}$ analytical','Location','southeast');
subplot(1,2,2);
plot(t2*1e3,x2(:,2),t2*1e3,x2(:,3)/n.R_n,'--');
xlabel('$t$ [ms]');
ylabel('current [A]');
legend('$i_n$','$u_C/R_n$','Location','southeast');
save_plot(fig,'Q1_5_nozzle_step');

s_cum=[0 cumsum(vecnorm(diff(ell_pos,1,2)))];
L_path=s_cum(end);
t_print=ell_t(end);
s_print=@(t)interp1(ell_t,s_cum,min(t,t_print));
res=filament(s_print,t_print,L_path,p.N_cells,p);
res2=filament(s_print,t_print,L_path,2*p.N_cells,p);
c=cell_coeffs(res.dx,p);
T_inf=(c.Gb*p.fil.T_b+c.Ge*n.T_e)/(c.Gb+c.Ge);

show('Q1.6 path length',L_path*1e3,'mm');
show('Q1.6 print time',t_print,'s');
show('Q1.6 dx',res.dx*1e3,'mm');
show('Q1.6 dt (used)',res.dt,'s');
show('Q1.6 dt_max stability',res.dt_max,'s');
show('Q1.6 dt_dep (1 cell per step)',res.dt_dep,'s');
show('Q1.6 local time constant C/(G_b+G_e)',c.C/(c.Gb+c.Ge),'s');
show('Q1.6 equilibrium temperature (bed+air)',T_inf-C2K,'degC');
show('Q1.6 cooling time after print start',res.t_cool,'s');
show('Q1.6 cooling time after print end',res.t_cool-t_print,'s');
show('Q1.6 cooling time after print end (2N)',res2.t_cool-t_print,'s');
show('Q1.6 mean final temperature',mean(res.T(:,end))-C2K,'degC');

Tmax=NaN(size(res.t));
Tmin=Tmax;
Tmean=Tmax;
for k=1:numel(res.t)
    Tk=res.T(~isnan(res.T(:,k)),k);
    if isempty(Tk)
        continue
    end
    Tmax(k)=max(Tk);
    Tmin(k)=min(Tk);
    Tmean(k)=mean(Tk);
end

fig=figure('Position',[100 100 800 420]);
imagesc(res.t,res.s*1e3,res.T-C2K);
axis xy;
cb=colorbar;
ylabel(cb,'$T$ [$^\circ$C]','Interpreter','latex');
set(gca,'Color',[1 1 1]);
xlabel('$t$ [s]');
ylabel('$s$ along path [mm]');
title('$T(s,t)$; white = not printed yet');
colormap(hot);
save_plot(fig,'Q1_6_heatmap');

fig=figure('Position',[100 100 800 400]);
plot(res.t,Tmax-C2K,res.t,Tmin-C2K,res.t,Tmean-C2K);
hold on;
plot([res.t(1) res.t(end)],[T_inf T_inf]-C2K,'k--');
yl=ylim;
plot([1;1]*[t_print res.t_cool],yl',':','Color',[0.4 0.4 0.4],'HandleVisibility','off');
ylim(yl);
xlabel('$t$ [s]');
ylabel('$T$ [$^\circ$C]');
legend('max','min','mean','$T_\infty$ (bed + air)','Location','northeast');
save_plot(fig,'Q1_6_min_max_mean');


function p=params()
p.K=0.5*[1 1;1 -1];

p.motor.r_m=6e-3;
p.motor.L_m=0.238e-3;
p.motor.R_m=2.32;
p.motor.J_m=10.8e-7;
p.motor.K_m=23.4e-3;

p.head.m=0.3;
p.head.m_g=0;
p.head.mu_r=0;
p.head.mu_lx=0;
p.head.mu_ly=0;

p.noz.R_n=24^2/40;
p.noz.L_n=1e-3;
p.noz.C_n=1e-6;
p.noz.C_tn=1.5*(0.02*0.02*0.01)*2.7e3*900;
p.noz.R_tn=7.5;
p.noz.T_e=25+273.15;
p.noz.T_print=210+273.15;
p.noz.u_print=sqrt((p.noz.T_print-p.noz.T_e)*p.noz.R_n/p.noz.R_tn);

p.fil.a=0.4e-3;
p.fil.rho=1300;
p.fil.c=1800;
p.fil.gamma=0.13;
p.fil.h_b=250;
p.fil.h_e=62;
p.fil.T_b=60+273.15;
p.fil.eps=1;

p.N_cells=300;
p.opts=odeset('RelTol',1e-8,'AbsTol',1e-11);
p.opts_stiff=odeset('RelTol',1e-6,'AbsTol',1e-6,'InitialStep',1e-7);
end

function show(name,val,unit,fmt)
if nargin<4
    fmt='%.4g';
end
fprintf(['%-40s ' fmt ' %s\n'],name,val,unit);
end

function save_plot(fig,name)
set(fig,'PaperPositionMode','auto');
print(fig,fullfile('figures',[name '.pdf']),'-dpdf');
print(fig,fullfile('figures',[name '.png']),'-dpng','-r200');
end

function xy=kin(th,p)
xy=p.motor.r_m*p.K*th;
end

function dx=f_full(x,u,p)
m=p.motor;
J=m.r_m*p.K;
M=m.J_m*eye(2)+J'*diag([p.head.m,p.head.m+p.head.m_g])*J;
D=p.head.mu_r*eye(2)+J'*diag([p.head.mu_lx,p.head.mu_ly])*J;
om=x(3:4);
i=x(5:6);
dx=[om;
    M\(m.K_m*i-D*om);
    (u(:)-m.R_m*i-m.K_m*om)/m.L_m];
end

function a=axis_coeffs(p)
r=p.motor.r_m;
a.G=2*p.motor.K_m/r;
a.kE=p.motor.K_m/r;
a.M=p.head.m+2*p.motor.J_m/r^2+[0;p.head.m_g];
a.c=2*p.head.mu_r/r^2+[p.head.mu_lx;p.head.mu_ly];
end

function dz=f_axis(z,u,s,p)
a=axis_coeffs(p);
dz=[z(2);
    (a.G*z(3)-a.c(s)*z(2))/a.M(s);
    (u-p.motor.R_m*z(3)-a.kE*z(2))/p.motor.L_m];
end

function [us,is]=feedforward(r,p)
a=axis_coeffs(p);
is=(a.M.*r.a+a.c.*r.v)/a.G;
dis=(a.M.*r.j+a.c.*r.a)/a.G;
us=p.motor.L_m*dis+p.motor.R_m*is+a.kE*r.v;
end

function r=ref_path(shape,t,sp)
t=t(:)';
w=2*pi/sp.T;
if strcmp(shape,'ellipse')
    c=cos(w*t);
    s=sin(w*t);
    r.p=[sp.xc+sp.ax*c;sp.yc+sp.ay*s];
    r.v=[-sp.ax*w*s;sp.ay*w*c];
    r.a=[-sp.ax*w^2*c;-sp.ay*w^2*s];
    r.j=[sp.ax*w^3*s;-sp.ay*w^3*c];
else
    s=w*t;
    S=sp.scale/32;
    A=16*S;
    x0=A*sin(s).^3;
    x1=A*w*3*sin(s).^2.*cos(s);
    x2=A*w^2*(6*sin(s).*cos(s).^2-3*sin(s).^3);
    x3=A*w^3*(6*cos(s).^3-21*sin(s).^2.*cos(s));
    ck=[13 -5 -2 -1];
    y0=0;
    y1=0;
    y2=0;
    y3=0;
    for k=1:4
        y0=y0+ck(k)*cos(k*s);
        y1=y1-ck(k)*k*sin(k*s);
        y2=y2-ck(k)*k^2*cos(k*s);
        y3=y3+ck(k)*k^3*sin(k*s);
    end
    r.p=[sp.xc+x0;sp.yc+S*y0];
    r.v=[x1;S*w*y1];
    r.a=[x2;S*w^2*y2];
    r.j=[x3;S*w^3*y3];
end
end

function dx=f_nozzle(x,u,p)
n=p.noz;
dx=[(x(3)^2/n.R_n-(x(1)-n.T_e)/n.R_tn)/n.C_tn;
    (u-x(3))/n.L_n;
    (x(2)-x(3)/n.R_n)/n.C_n];
end

function c=cell_coeffs(dx,p)
a=p.fil.a;
c.C=p.fil.rho*p.fil.c*a^2*dx;
c.Gax=p.fil.gamma*a^2/dx;
c.Gb=p.fil.h_b*a*dx;
c.Ge=3*p.fil.h_e*a*dx;
c.Gend=p.fil.h_e*a^2;
end

function res=filament(s_print,t_print,L_path,N,p)
dx=L_path/N;
c=cell_coeffs(dx,p);
dt_max=c.C/(2*max(c.Gax,c.Gend)+c.Gb+c.Ge);
dt_dep=t_print/N;
dt=min(0.9*dt_max,dt_dep);

T=NaN(N,1);
active=false(N,1);
t=0;
t_store=0;
t_hist=[];
T_hist=[];
idx=(1:N)';
while t<t_print+1e4
    if t<=t_print
        n_act=min(N,floor(s_print(t)/dx+1e-9));
    else
        n_act=N;
    end
    new=find(~active(1:n_act));
    T(new)=p.noz.T_print;
    active(new)=true;

    if all(active)
        left=mod(idx-2,N)+1;
        right=mod(idx,N)+1;
        hasL=true(N,1);
        hasR=true(N,1);
    else
        left=max(idx-1,1);
        right=min(idx+1,N);
        hasL=idx>1&active(left);
        hasR=idx<N&active(right);
    end
    dTl=T(left)-T;
    dTl(isnan(dTl))=0;
    dTr=T(right)-T;
    dTr(isnan(dTr))=0;
    q=c.Gax*(hasL.*dTl+hasR.*dTr)+c.Gb*(p.fil.T_b-T)+c.Ge*(p.noz.T_e-T) ...
        +c.Gend*(~hasL+~hasR).*(p.noz.T_e-T);
    T(active)=T(active)+dt*q(active)/c.C;
    t=t+dt;

    if t>=t_store
        t_hist(end+1)=t;
        T_hist(:,end+1)=T;
        t_store=t_store+t_print/200*(1+9*(t>t_print));
    end
    if t>t_print&&all(active)&&max(T)-min(T)<p.fil.eps
        t_hist(end+1)=t;
        T_hist(:,end+1)=T;
        break
    end
end
res=struct('t',t_hist,'T',T_hist,'s',(idx-0.5)*dx,'t_cool',t, ...
    'dt',dt,'dt_max',dt_max,'dt_dep',dt_dep,'dx',dx);
end
