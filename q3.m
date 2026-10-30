clear;
close all;
clc;
p=params();
L=layout();
C2K=273.15;
if ~isfolder(fullfile(pwd,'figures'))
    mkdir(fullfile(pwd,'figures'));
end
set(groot,'defaultAxesFontSize',11,'defaultLineLineWidth',1.3,'defaultAxesXGrid','on','defaultAxesYGrid','on');
set(groot,'defaultTextInterpreter','latex','defaultAxesTickLabelInterpreter','latex','defaultLegendInterpreter','latex');
ellipse=struct('xc',0.11,'yc',0.11,'ax',0.06,'ay',0.04,'T',8);

H=filament_automaton(ellipse,3,p);
sim=hybrid_sim(H,115,p);

for k=1:numel(sim.events)
    ev=sim.events(k);
    show(sprintf('Q3.1 event %02d %s (%s->%s)',k,ev.name,ev.from,ev.to),ev.t,'s','%.3f');
end

modes={'printing','stop_save','move_wipe','wipe','heat_off','load','heat_on','resume'};
[~,qi]=ismember(sim.q,modes);
pos=kin(sim.x(:,L.th)',p);
t_ev=[sim.events.t];

fig=figure('Position',[100 100 900 800]);
subplot(4,1,1);
stairs(sim.t,qi);
set(gca,'YTick',1:8,'YTickLabel',strrep(modes,'_','\_'));
ylim([0.5 8.5]);
ylabel('mode');
subplot(4,1,2);
plot(sim.t,sim.x(:,L.Tn)-C2K);
ylabel('$T_n$ [$^\circ$C]');
subplot(4,1,3);
plot(sim.t,pos'*1e3);
ylabel('position [mm]');
legend('$x$','$y$','Location','best');
subplot(4,1,4);
plot(sim.t,sim.x(:,L.dist)*1e3,sim.t,sim.x(:,L.spool)*1e3,'--');
ylabel('[mm]');
xlabel('$t$ [s]');
legend('printed distance','spool empty at','Location','best');
for s=1:4
    subplot(4,1,s);
    vlines(t_ev);
end
save_plot(fig,'Q3_1_hybrid_filament');

pr=strcmp(sim.q,'printing');
fig=figure('Position',[100 100 500 450]);
plot(pos(1,:)*1e3,pos(2,:)*1e3,'Color',[0.7 0.7 0.7]);
hold on;
plot(pos(1,pr)*1e3,pos(2,pr)*1e3,'.','MarkerSize',4);
plot(0,0,'ks','MarkerSize',8);
axis equal;
xlabel('$x$ [mm]');
ylabel('$y$ [mm]');
legend('all motion','printing','wipe tower','Location','best');
save_plot(fig,'Q3_1_path');

r=ref_path('ellipse',sim.x(pr,L.phi)',ellipse);
show('Q3.1 max tracking error while printing',max(vecnorm(pos(:,pr)-r.p))*1e3,'mm','%.3g');
show('Q3.1 max motor current',max(max(abs(sim.x(:,L.im)))),'A','%.3g');

H=square_automaton([0.06;0.06],0.10,p);
sim=hybrid_sim(H,10,p);
pos=kin(sim.x(:,L.th)',p);

for k=1:4
    in=strcmp(sim.q,sprintf('side%d',k));
    d=1+mod(k,2);
    show(sprintf('Q3.2 max deviation side %d',k),max(abs(pos(d,in)-H.corners(d,k)))*1e6,'um','%.3g');
end
for k=1:numel(sim.events)
    pc=pos(:,find(sim.t>=sim.events(k).t,1));
    show(sprintf('Q3.2 corner error at event %d',k),norm(pc-H.corners(:,k+1))*1e6,'um','%.3g');
end

fig=figure('Position',[100 100 900 400]);
subplot(1,2,1);
plot(H.corners(1,:)*1e3,H.corners(2,:)*1e3,'k--',pos(1,:)*1e3,pos(2,:)*1e3);
axis equal;
xlabel('$x$ [mm]');
ylabel('$y$ [mm]');
legend('ideal','simulation','Location','best');
subplot(1,2,2);
plot(sim.t,sim.x(:,L.om));
ylabel('$\omega$ [rad/s]');
xlabel('$t$ [s]');
legend('$\omega_l$','$\omega_r$','Location','best');
vlines([sim.events.t]);
save_plot(fig,'Q3_2_square');


function p=params()
p.K=0.5*[1 1;1 -1];

p.motor.r_m=6e-3;
p.motor.L_m=0.238e-3;
p.motor.R_m=2.32;
p.motor.J_m=10.8e-7;
p.motor.K_m=23.4e-3;
p.motor.u_max=24;

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
p.noz.T_high=205+273.15;
p.noz.T_low=180+273.15;
p.noz.u_max=24;
p.noz.u_print=sqrt((p.noz.T_print-p.noz.T_e)*p.noz.R_n/p.noz.R_tn);

M=p.head.m+2*p.motor.J_m/p.motor.r_m^2;
g=2*p.motor.K_m/(p.motor.r_m*p.motor.R_m);
be=2*p.motor.K_m^2/(p.motor.r_m^2*p.motor.R_m);
w=2*pi*10;
p.Kp=M*w^2/g;
p.Kd=max(0,(2*M*w-be)/g);

p.T_move=1.5;
p.A_wipe=5e-3;
p.T_wipe=1;
p.L_wipe=4*p.A_wipe;
p.T_load=5;
p.L_spool0=0.25;
p.L_spool=0.15;
p.om_eps=0.05;
p.e_tol=5e-6;

p.opts=odeset('RelTol',1e-6,'AbsTol',1e-6,'InitialStep',1e-7,'MaxStep',0.02);
end

function L=layout()
L.mech=1:6;
L.th=1:2;
L.om=3:4;
L.im=5:6;
L.noz=7:9;
L.Tn=7;
L.clk=10;
L.dist=11;
L.phi=12;
L.sv=13:14;
L.a=15:16;
L.dmark=17;
L.spool=18;
L.n=18;
end

function show(name,val,unit,fmt)
fprintf(['%-40s ' fmt ' %s\n'],name,val,unit);
end

function save_plot(fig,name)
set(fig,'PaperPositionMode','auto');
print(fig,fullfile('figures',[name '.pdf']),'-dpdf');
print(fig,fullfile('figures',[name '.png']),'-dpng','-r200');
end

function vlines(tv)
yl=ylim;
hold on;
plot([1;1]*tv,yl',':','Color',[0.4 0.4 0.4],'HandleVisibility','off');
ylim(yl);
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

function dx=f_nozzle(x,u,p)
n=p.noz;
dx=[(x(3)^2/n.R_n-(x(1)-n.T_e)/n.R_tn)/n.C_tn;
    (u-x(3))/n.L_n;
    (x(2)-x(3)/n.R_n)/n.C_n];
end

function dx=sys(x,u,p)
L=layout();
dx=zeros(L.n,1);
dx(L.mech)=f_full(x(L.mech),u(1:2),p);
dx(L.noz)=f_nozzle(x(L.noz),u(3),p);
dx(L.clk)=1;
dx(L.dist)=u(4)*norm(kin(x(L.om),p));
dx(L.phi)=u(5);
end

function [us,is]=feedforward(r,p)
r_m=p.motor.r_m;
G=2*p.motor.K_m/r_m;
kE=p.motor.K_m/r_m;
M=p.head.m+2*p.motor.J_m/r_m^2+[0;p.head.m_g];
c=2*p.head.mu_r/r_m^2+[p.head.mu_lx;p.head.mu_ly];
is=(M.*r.a+c.*r.v)/G;
dis=(M.*r.j+c.*r.a)/G;
us=p.motor.L_m*dis+p.motor.R_m*is+kE*r.v;
end

function u=motion_input(x,r,p)
L=layout();
us=feedforward(r,p)+p.Kp*(r.p-kin(x(L.th),p))+p.Kd*(r.v-kin(x(L.om),p));
u=p.K\us;
u=min(max(u,-p.motor.u_max),p.motor.u_max);
end

function r=ref_path(shape,t,sp)
t=t(:)';
if strcmp(shape,'ellipse')
    w=2*pi/sp.T;
    c=cos(w*t);
    s=sin(w*t);
    r.p=[sp.xc+sp.ax*c;sp.yc+sp.ay*s];
    r.v=[-sp.ax*w*s;sp.ay*w*c];
    r.a=[-sp.ax*w^2*c;-sp.ay*w^2*s];
    r.j=[sp.ax*w^3*s;-sp.ay*w^3*c];
else
    tau=min(max(t/sp.T,0),1);
    in=(t>=0)&(t<=sp.T);
    d=sp.p1-sp.p0;
    r.p=sp.p0+d*(10*tau.^3-15*tau.^4+6*tau.^5);
    r.v=d*((30*tau.^2-60*tau.^3+30*tau.^4)/sp.T.*in);
    r.a=d*((60*tau-180*tau.^2+120*tau.^3)/sp.T^2 .*in);
    r.j=d*((60-360*tau+360*tau.^2)/sp.T^3 .*in);
end
end

function r=hold_at(pt)
r=struct('p',pt,'v',[0;0],'a',[0;0],'j',[0;0]);
end

function r=move(x,target,p)
L=layout();
r=ref_path('minjerk',x(L.clk),struct('p0',x(L.a),'p1',target,'T',p.T_move));
end

function d=still_at(x,target,p)
L=layout();
d=max(max(abs(x(L.om)))-p.om_eps,norm(kin(x(L.th),p)-target)-p.e_tol);
end

function q_next=fsm_step(q,e)
next={'printing',{'C','F'},'stop_save'
    'stop_save',{'stopped'},'move_wipe'
    'move_wipe',{'arrived'},'wipe'
    'wipe',{'wiped'},'heat_off'
    'heat_off',{'T_low'},'load'
    'load',{'loaded'},'heat_on'
    'heat_on',{'T_high'},'resume'
    'resume',{'returned'},'printing'};
k=find(strcmp(next(:,1),q));
if any(strcmp(e,next{k,2}))
    q_next=next{k,3};
else
    q_next=q;
end
end

function sim=hybrid_sim(H,t_end,p)
t=0;
x=H.x0;
q=H.q0;
sim.t=[];
sim.x=[];
sim.q={};
sim.events=struct('t',{},'from',{},'name',{},'to',{});
while t<t_end
    k_ext=find(H.t_ext>t,1);
    t_stop=min([t_end,H.t_ext(k_ext)]);
    [~,name,dir]=H.guard(q,x);
    opts=odeset(p.opts,'Events',@(tt,xx)ev_fun(H.guard(q,xx),dir));
    [ts,xs,te,~,ie]=ode15s(@(tt,xx)sys(xx,H.input(q,xx),p),[t t_stop],x,opts);

    sim.t=[sim.t;ts];
    sim.x=[sim.x;xs];
    sim.q=[sim.q;repmat({q},numel(ts),1)];
    t=ts(end);
    x=xs(end,:)';
    if ~isempty(ie)
        t=te(end);
    elseif ~isempty(k_ext)&&abs(t-H.t_ext(k_ext))<1e-12
        name='C';
    else
        break
    end
    [q_next,x]=H.jump(q,name,x);
    sim.events(end+1)=struct('t',t,'from',q,'name',name,'to',q_next);
    q=q_next;
end
end

function [value,isterminal,direction]=ev_fun(g,dir)
if isempty(g)
    value=1;
    isterminal=0;
    direction=0;
else
    value=g;
    isterminal=1;
    direction=dir;
end
end

function H=filament_automaton(sp,t_C,p)
L=layout();
r0=ref_path('ellipse',0,sp);
[~,is0]=feedforward(r0,p);
x0=zeros(L.n,1);
x0(L.mech)=[(p.motor.r_m*p.K)\r0.p;(p.motor.r_m*p.K)\r0.v;p.K\is0];
x0(L.noz)=[p.noz.T_print;p.noz.u_print/p.noz.R_n;p.noz.u_print];
x0(L.sv)=r0.p;
x0(L.a)=r0.p;
x0(L.spool)=p.L_spool0;

H.q0='printing';
H.x0=x0;
H.t_ext=t_C;
H.input=@(q,x)fil_input(q,x,sp,p);
H.guard=@(q,x)fil_guard(q,x,p);
H.jump=@(q,e,x)fil_jump(q,e,x,p);
end

function u=fil_input(q,x,sp,p)
L=layout();
tower=[0;0];
u_n=p.noz.u_print;
extrude=0;
phase=0;
switch q
    case 'printing'
        r=ref_path('ellipse',x(L.phi),sp);
        extrude=1;
        phase=1;
    case 'stop_save'
        r=hold_at(x(L.sv));
    case 'move_wipe'
        r=move(x,tower,p);
    case 'wipe'
        w=2*pi/p.T_wipe;
        A=p.A_wipe;
        c=x(L.clk);
        r=struct('p',tower+[A*(1-cos(w*c));0],'v',[A*w*sin(w*c);0], ...
            'a',[A*w^2*cos(w*c);0],'j',[-A*w^3*sin(w*c);0]);
        extrude=1;
    case {'heat_off','load'}
        r=hold_at(x(L.a));
        u_n=0;
    case 'heat_on'
        r=hold_at(x(L.a));
        u_n=p.noz.u_max;
    case 'resume'
        r=move(x,x(L.sv),p);
end
u=[motion_input(x,r,p);u_n;extrude;phase];
end

function [g,name,dir]=fil_guard(q,x,p)
L=layout();
switch q
    case 'printing'
        g=x(L.dist)-x(L.spool);
        name='F';
        dir=1;
    case 'stop_save'
        g=still_at(x,x(L.sv),p);
        name='stopped';
        dir=-1;
    case 'move_wipe'
        g=still_at(x,[0;0],p);
        name='arrived';
        dir=-1;
    case 'wipe'
        g=x(L.dist)-x(L.dmark)-p.L_wipe;
        name='wiped';
        dir=1;
    case 'heat_off'
        g=x(L.Tn)-p.noz.T_low;
        name='T_low';
        dir=-1;
    case 'load'
        g=x(L.clk)-p.T_load;
        name='loaded';
        dir=1;
    case 'heat_on'
        g=x(L.Tn)-p.noz.T_high;
        name='T_high';
        dir=1;
    case 'resume'
        g=still_at(x,x(L.sv),p);
        name='returned';
        dir=-1;
end
end

function [q_next,x]=fil_jump(q,e,x,p)
L=layout();
q_next=fsm_step(q,e);
if strcmp(q_next,q)
    return
end
pm=kin(x(L.th),p);
x(L.clk)=0;
x(L.a)=pm;
if strcmp(q_next,'stop_save')
    x(L.sv)=pm;
end
if strcmp(q_next,'wipe')
    x(L.dmark)=x(L.dist);
end
if strcmp(q_next,'heat_on')
    x(L.spool)=x(L.dist)+p.L_spool;
end
end

function H=square_automaton(c0,side,p)
L=layout();
corners=c0+side*[0 1 1 0 0;0 0 1 1 0];
x0=zeros(L.n,1);
x0(L.th)=(p.motor.r_m*p.K)\c0;
x0(L.noz)=[p.noz.T_print;p.noz.u_print/p.noz.R_n;p.noz.u_print];
x0(L.a)=corners(:,1);
x0(L.sv)=corners(:,2);

H.q0='side1';
H.x0=x0;
H.corners=corners;
H.t_ext=[];
H.input=@(q,x)sq_input(q,x,p);
H.guard=@(q,x)sq_guard(q,x,p);
H.jump=@(q,e,x)sq_jump(q,x,corners);
end

function u=sq_input(q,x,p)
L=layout();
if strcmp(q,'done')
    r=hold_at(x(L.sv));
else
    r=move(x,x(L.sv),p);
end
u=[motion_input(x,r,p);p.noz.u_print;1;0];
end

function [g,name,dir]=sq_guard(q,x,p)
L=layout();
if strcmp(q,'done')
    g=[];
    name='';
    dir=[];
else
    g=still_at(x,x(L.sv),p);
    name='corner';
    dir=-1;
end
end

function [q_next,x]=sq_jump(q,x,corners)
L=layout();
k=sscanf(q,'side%d');
if k<4
    q_next=sprintf('side%d',k+1);
else
    q_next='done';
end
x(L.clk)=0;
x(L.a)=corners(:,k+1);
x(L.sv)=corners(:,min(k+2,5));
end
