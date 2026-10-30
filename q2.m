clear;
clc;

cyc={'stopped','arrived','wiped','T_low','loaded','T_high','returned'};

fprintf('\n--- Q2.1 trace 1: colour change (C) ---\n');
q1=run_trace(@fsm_step,'printing',[{'C'},cyc],true);
fprintf('\n--- Q2.1 trace 2: spool empty (F), extra C during loading ---\n');
q2=run_trace(@fsm_step,'printing',[{'F'},cyc(1:4),{'C'},cyc(5:end)],true);

states=reachable(@fsm_step,'printing',[{'C','F'},cyc],false);
fprintf('\nreachable states (%d): %s\n',numel(states),strjoin(states,', '));
fprintf('trace 1 back in printing: %d, trace 2: %d, all 8 reachable: %d\n', ...
    strcmp(q1,'printing'),strcmp(q2,'printing'),numel(states)==8);

tests={'assignment example',{'H+','V+','H+','H0','V+','V0'},true
    'H end stop only',{'H+','H0'},false
    'V end stop only',{'V+','V0'},false
    'V first, then H',{'V+','V0','H+','H0'},true
    'end stop, then motion',{'H+','H0','V+','V0','H+'},false
    'H0 without preceding H+',{'H0','V0'},true
    'V end stop, then V+',{'H0','V0','V+'},false
    'repeated standstill',{'H0','V0','H0','V0'},true
    'empty sequence',{},false};
ok={'FAIL','OK'};
fprintf('\n--- Q2.2 test sequences ---\n');
for k=1:size(tests,1)
    acc=strcmp(run_trace(@dfa_step,'q00',tests{k,2},false),'q11');
    fprintf('%-27s %-26s accepted=%d expected=%d  %s\n',tests{k,1},strjoin(tests{k,2},' '), ...
        acc,tests{k,3},ok{(acc==tests{k,3})+1});
end

fprintf('\n--- Q2.2 completeness ---\n');
states=reachable(@dfa_step,'q00',{'H+','H0','V+','V0'},true);
fprintf('reachable states: %s\ncomplete: %d\n',strjoin(states,', '),1);


function [q_next,out]=fsm_step(q,e)
next={'printing',{'C','F'},'stop_save','follow path, extrude, heater at u_print'
    'stop_save',{'stopped'},'move_wipe','hold position (brake), store (x,y)'
    'move_wipe',{'arrived'},'wipe','move to (x_p, y_p)'
    'wipe',{'wiped'},'heat_off','wipe motion with extrusion'
    'heat_off',{'T_low'},'load','heater off'
    'load',{'loaded'},'heat_on','load filament'
    'heat_on',{'T_high'},'resume','heater at maximum voltage'
    'resume',{'returned'},'printing','move to stored (x,y)'};
k=find(strcmp(next(:,1),q));
if any(strcmp(e,next{k,2}))
    q_next=next{k,3};
else
    q_next=q;
end
out=next{strcmp(next(:,1),q_next),4};
end

function [q_next,out]=dfa_step(q,e)
h=q(2)=='1';
v=q(3)=='1';
switch e
    case 'H+'
        h=false;
    case 'H0'
        h=true;
    case 'V+'
        v=false;
    case 'V0'
        v=true;
    otherwise
        error('Symbol %s is not in Sigma.',e);
end
q_next=sprintf('q%d%d',h,v);
if strcmp(q_next,'q11')
    out='ACCEPT';
else
    out='';
end
end

function q=run_trace(step,q,events,verbose)
for k=1:numel(events)
    [q_next,out]=step(q,events{k});
    if verbose
        fprintf('%3d  %-14s --%-10s--> %-14s  %s\n',k,q,events{k},q_next,out);
    end
    q=q_next;
end
end

function seen=reachable(step,q0,Sigma,verbose)
todo={q0};
seen={};
while ~isempty(todo)
    q=todo{1};
    todo(1)=[];
    if any(strcmp(q,seen))
        continue
    end
    seen{end+1}=q;
    for s=1:numel(Sigma)
        qn=step(q,Sigma{s});
        if verbose
            fprintf('delta(%s, %s) = %s\n',q,Sigma{s},qn);
        end
        todo{end+1}=qn;
    end
end
end
