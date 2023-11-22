%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SIR_MultiStart.m
%Christina Edholm
%
%
%
% This code fits a SIR model for LA county data on COVID-19 - https://github.com/datadesk/california-coronavirus-data and http://publichealth.lacounty.gov/media/coronavirus/
%
% To run this code you need to give how many Multistart runs you want
% NoStartPoints. Also pick which time frame you are fitting during the
% COVID-19 pandemic.
%
% This code calls the SIR_Model.m -- ODE equations
% 
% This code will plot results.
% 
% You can change the LowerBounds and UpperBounds below for the parameters
% we are estimating, beta. This can be changed!
% 
%
% November 2, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all %clear the workspace

%% How many runs for parameter estimation do you want:

NoStartPoints=10;

%% Given Data on La County - repeated below in COVID_RUN_ODE45 -- http://publichealth.lacounty.gov/media/coronavirus/

% LA County Data - April 20, 2020 (START DATE)  through May 16, 2020
    ConfCase=[13823	15165	16449	17567	18545	19159	19567	20460	21017	22522	23233	24262	24936	25699	26238	27866	28665	29526	30334	31241	31703	32269	33247	34552	35447	36324	37374];
    Deaths=[619	666	732	798	850	896	916	948	1004	1065	1119	1174	1212	1231	1260	1315	1369	1420	1470	1515	1531	1570	1617	1663	1711	1752	1793];
    Hpatient=[2393	2471	2448	2458	2406	2417	2462	2549	2543	2550	2364	2395	2380	2299	2371	2370	2292	2323	2276	2221	2270	2310	2248	2351	2207	2201	2086];
    ten_dayCumulative = 5370;
    CumulativeTo10DayStart = 8453;

    %Total Population considered for outbreak
    TotalPopulation = 9651332; 
    
% % LA County Data - June 13, 2020 through June 30, 2020
%     ConfCase = [74123,74717,76679,78486,80574,82623,84549,86072,86934,89210,91562,93959,96641,99008,100526,101070,101768,101945];
%     Deaths = [2833,2857,2886,2909,2946,2968,3001,3030,3053,3075,3090,3119,3138,3154,3170,3191,3210,3224];
%     Hpatient_rev = [2538, 2511, 2465, 2316, 2296, 2237, 2283, 2259, 2090, 2001, 1969, 1920, 1900, 1958, 1947, 1956, 1797, 1768];
%     Hpatient = fliplr(Hpatient_rev);
%     ten_dayCumulative = 13615;
%     CumulativeTo10DayStart = 60508;

day=length(ConfCase);           %%how many days of COVID-19 data

tspan = 1:1:day;                %%vector how many days you are simulating, will be inputed into ode45

I0 = ten_dayCumulative;                %%Add together all the cases 10 days prior to the start date
R0 = .95*20*CumulativeTo10DayStart;    %% Changed to a percentage of cumulative infected up to 10 days before start date
CI0 = ConfCase(1);              %%total infections

S0 = TotalPopulation-I0-R0;


initialvalues = [S0,I0,R0,CI0]; 

%% Upper and Lower bounds for parameters you are fitting.

%Parameter vector we are approximating:
% z = [beta] 

LowerBounds=[1.00e-9];      %Lowerbounds for the parameters you are estimating
UpperBounds=[1.00e-07];       %Upperbounds for the parameters you are estimating

xstart=.5*(LowerBounds+UpperBounds);                            %What initial parameter values you want to start with

%% MultiStart and fmincon - Fitting Part - Parallelization - not many comments ask Prof. Edholm for clarification if you want.

% Here we set-up the optimization problem, specificying we will use fmincon
% as the local solver, and the what model we want to minimize along with
% the specific measure down below in the SIR_RUN_ODE45 function. We give
% initial conditions and the bounds.
problem = createOptimProblem('fmincon','objective',@SIR_RUN_ODE45...
,'x0',xstart,'lb',LowerBounds,'ub',UpperBounds);%,'Aineq',A,'bineq',b)%,'Aeq',Aeq,'beq',beq);

problem.options = optimoptions(problem.options,'MaxFunEvals',9999,'MaxIter',9999);%,'TolFun',0,'TolCon',0)
%problem.options = optimoptions(problem.options,'MaxFunEvals',inf,'MaxIter',inf,'TolFun',1e-10,'TolCon',0,'TolX',0,'MaxFunEvals',999999)

numstartpoints=NoStartPoints;                               % How many runs do you want?

% %  ms=MultiStart('Display','iter');                       %defines a multistart problem without parallel

ms=MultiStart('UseParallel',true,'Display','iter');         %defines a parallel multistart problem

%parpool %accesses the cores for parallel on your computer (laptop goes for 2-8, can be more specific)

[b,fval,exitflag,output,manymins]=run(ms,problem,numstartpoints);  %runs the multistart 

% the following takes solutions from manymins and makes a matrix out of them


for i=1:length(manymins)
    SIRParameters(i,:)=manymins(i).X;       %what are the parameter values
end

for i=1:length(manymins)
    fvalues(i)=manymins(i).Fval;            %the minimization error
end

for i=1:length(manymins)
    ExitFlags(i)=manymins(i).Exitflag;      %how "good" is the solution, we want 1 or 2.
end


%delete(gcp('nocreate'))  %turns off the parallel feature




%% Plot the "best" solution



%%Outputs state variables for "best" fit
[t,y] = ode45(@(t,y) SIR_Model(t,y,SIRParameters(1,:)),tspan,initialvalues);

       S = y(:,1);
        
       I = y(:,2);
        
       R = y(:,3);
        
       CI = y(:,4);


    
   figure(1)
    tiledlayout(1,4)
    nexttile
     hold all
    plot(tspan,CI)
    scatter(tspan,ConfCase, 'filled')
    title('Cumulative Cases')
    xlabel('Days')
    nexttile
    plot(tspan,S)
    title('Susceptible')
    xlabel('Days')
    nexttile
    plot(tspan,I)
    title('Infectious')
    xlabel('Days')
    nexttile
    plot(tspan,R)
    title('Recovered')
    xlabel('Days')


       
        
function value=SIR_RUN_ODE45(z)     %This is the function which sets up the ODE problem and also what we are minimizing
%% LA County Data April 20-May 16, 2020 
    
% LA County Data - April 20, 2020 (START DATE)  through May 16, 2020
    ConfCase=[13823	15165	16449	17567	18545	19159	19567	20460	21017	22522	23233	24262	24936	25699	26238	27866	28665	29526	30334	31241	31703	32269	33247	34552	35447	36324	37374];
    Deaths=[619	666	732	798	850	896	916	948	1004	1065	1119	1174	1212	1231	1260	1315	1369	1420	1470	1515	1531	1570	1617	1663	1711	1752	1793];
    Hpatient=[2393	2471	2448	2458	2406	2417	2462	2549	2543	2550	2364	2395	2380	2299	2371	2370	2292	2323	2276	2221	2270	2310	2248	2351	2207	2201	2086];
    ten_dayCumulative = 5370;
    CumulativeTo10DayStart = 8453;

    %Total Population considered for outbreak
    TotalPopulation = 9651332; 
    
% % LA County Data - June 13, 2020 through June 30, 2020
%     ConfCase = [74123,74717,76679,78486,80574,82623,84549,86072,86934,89210,91562,93959,96641,99008,100526,101070,101768,101945];
%     Deaths = [2833,2857,2886,2909,2946,2968,3001,3030,3053,3075,3090,3119,3138,3154,3170,3191,3210,3224];
%     Hpatient_rev = [2538, 2511, 2465, 2316, 2296, 2237, 2283, 2259, 2090, 2001, 1969, 1920, 1900, 1958, 1947, 1956, 1797, 1768];
%     Hpatient = fliplr(Hpatient_rev);
%     ten_dayCumulative = 13615;
%     CumulativeTo10DayStart = 60508;

day=length(ConfCase);           %%how many days of COVID-19 data

tspan = 1:1:day;                %%vector how many days you are simulating, will be inputed into ode45

I0 = ten_dayCumulative;                %%Add together all the cases 10 days prior to the start date
R0 = .95*20*CumulativeTo10DayStart;    %% Changed to a percentage of cumulative infected up to 10 days before start date
CI0 = ConfCase(1);              %%total infections

S0 = TotalPopulation-I0-R0;


initialvalues = [S0,I0,R0,CI0]; 

%% Run the ode45 solver

[t,y] = ode45(@(t,y) SIR_Model(t,y,z),tspan,initialvalues);
        
       CI = y(:,4);
       

%Difference Data and Model
% Make sure everything is the same size
%
diff = CI - reshape(ConfCase,size(CI));

value=norm(diff,2);

end