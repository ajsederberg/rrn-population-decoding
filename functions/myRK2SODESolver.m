function [t,y] = myRK2SODESolver(ddtfn,t0tf,y0,dt,sig2)
% myRK2SODESolver - updated 9/28/19. Not actually RK2. Just stochastic
% Euler for now. 
% ddtfn - handle to the ode function with parameters already accounted for.
% 

% if sig2 is a vector of poisitive numbbers:

if size(sig2, 2) == 1 && all(sig2 > 0)
    % y0 - IC
    % t0tf - [starttime stoptime]
    % param - any params for the ddt function
    % ddtfn - name of ddt function

    % ddt = str2func(ddtfn);

    % dt = param.dt;
    % sig2: variance of noise

    NT = ceil(diff(t0tf)/dt)+1;

    t = linspace(t0tf(1),t0tf(end),NT);

    y = zeros(length(y0),length(t));
    y(:,1) = y0;

    checkpoints = floor(linspace(1,NT,11));
    ci = 1;
    for n=2:NT
        y(:,n) = y(:,n-1) + dt*ddtfn(t(n-1),y(:,n-1)) + sqrt(dt*sig2).*randn(length(y0), 1);
        if n == checkpoints(ci)
            disp([num2str(10*(n-1)) '% done'])
        end
    end

    % return withh time along the first dimensions. 
    y = y';
    t = t';

else
    % y0 - IC
    % t0tf - [starttime stoptime]
    % param - any params for the ddt function
    % ddtfn - name of ddt function

    % ddt = str2func(ddtfn);

    % dt = param.dt;
    % sig2: variance of noise

    NT = ceil(diff(t0tf)/dt)+1;

    t = linspace(t0tf(1),t0tf(end),NT);

    y = zeros(length(y0),length(t));
    y(:,1) = y0;
    % number of noise components
    n_x_noise = size(sig2, 2);
    x_n = zeros(n_x_noise, length(t));
    tau_noise = 3;  

    checkpoints = floor(linspace(1,NT,11));
    ci = 1;
    for n=2:NT
        % integrate o-u process 
        x_n(:, n) = x_n(:, n - 1)*(1 - dt/tau_noise)+sqrt(dt)*randn(n_x_noise, 1);
        
        % each unit couples to o-u process as determined by projection
        % matrix
        y(:, n) = y(:,n-1) + dt*(ddtfn(t(n-1),y(:, n-1)) + sig2*x_n(:, n));
%         y(:,n) = y(:,n-1) + dt*ddtfn(t(n-1),y(:,n-1)) + sqrt(dt*sig2).*randn(length(y0), 1);
        if n == checkpoints(ci)
            disp([num2str(10*(n-1)) '% done'])
        end
    end

    % return withh time along the first dimensions. 
    y = y';
    t = t';
    
end