%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIR_Model.m
%Christina Edholm 
%
% This code is the ODE model we are using, first we establish parameters
% then write out the system of ODEs. 
%
% 
%
% November 2, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function dydt = f(t,y,z)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% COVID - 19   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Christina Edholm 

%Parameters
%We are approximating:
% z = [beta] 

beta = z(1);        % susc to inf
gamma = 1/6;       % inf to rec

% Initiate DE variables
dydt = zeros(4,1);

% Differential Equations -----------------------

S = y(1);
I = y(2);
R = y(3);

% Start DEs------
% dS/dt
dydt(1) = -beta*S*I;
% dE/dt
dydt(2) = beta*S*I-gamma*I;
% dI/dt
dydt(3) = gamma*I;
% %Keep track of cumulative infections - what we have data for to match
% with
dydt(4) = beta*S*I;

end



