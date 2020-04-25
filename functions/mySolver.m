function [t,y] = mySolver(ddtfn,t0tf,y0,opts,param)
% my ode solver
% y0 - IC
% t0tf - [starttime stoptime]
% param - any params for the ddt function
% ddtfn - name of ddt function

ddt = str2func(ddtfn);

dt = param.dt;

NT = ceil(diff(t0tf)/dt)+1;

t = linspace(t0tf(1),t0tf(end),NT);

y = zeros(length(y0),length(t));
y(:,1) = y0;

checkpoints = floor(linspace(1,NT,11));
ci = 1;
for n=2:NT
    y(:,n) = y(:,n-1) + dt*ddt(t(n-1),y(:,n-1),param);
    if n == checkpoints(ci)
        display([num2str(10*(n-1)) '% done'])
    end
end