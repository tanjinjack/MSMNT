% MSMNT method as proposed by Chang LIU in in her thesis.
% Thesis link: https://research.tue.nl/nl/publications/in-situ-characterization-of-the-acoustic-impedance-of-vegetated-r
% The code can be understood to have 2 parts.
% First, some parameters are defined that are to be used by the function
% processImpulseResponse(), the key output being the level difference
% between the two microphones.
% Then, few more parameters are defined (paramMNT) before the minimization
% routine starts (via patternsearch).

clear all
close all

% see processImpulseResponse() for detail of each variable
imp_param.datafolder      = 'Data/Soft_foam/';%'Data/Hard_foam/'%'Data/Concrete/'
imp_param.n_layer         = 1; 
imp_param.guess_thickness = 0.1; 
imp_param.n_meas          = 5;
imp_param.fs              = 48000;
imp_param.n_config        = 3;
imp_param.flow            = 200;
imp_param.fhigh           = 2500;

% here we define our initial guess parameters, lower and upper boundaries
% of the search space.
% for material where n_layer = 0, we will have two parameters, in the
% following order: flow resistivity, porosity.
% for material where n_layer = 1, we will have three parameters, in the
% following order: flow resistivity, porosity, thickness.
% for material where n_layer = 2, we could have seven parameters, in the
% following order: NOT IMPLEMENTED YET.
guess_param = [3.5e3 0.97 0.13];
lowlim = [1e3 0.8 0.1];
highlim =[10e3 0.99 0.2];

% if imp_param.n_layer == 0
%     guess_param = [2e4 0.01];
%     lowlim = [2e4 0.01];
%     highlim =[1e5 0.3];
% elseif imp_param.n_layer == 1
%     guess_param = [3.5e3 0.97 0.13];
%     lowlim = [1e3 0.8 0.1];
%     highlim =[10e3 0.99 0.2];
% %     guess_param = [4e3 0.8 0.11];
% %     lowlim = [1e3 0.8 0.1];
% %     highlim =[20e3 0.99 0.2];
% else
%     disp('not sure how to start guessing');
%     return
% end

% option to let patternsearch to adjust the physical configurations 
% slightly for potential better match. unique of result might be affected.
param.var_dim = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% you should not be changing anything below this line unless you know
% what exactly you are doing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUT_IR = processImpulseResponse(imp_param);

% see solveMNT for details about variables in paramMNT
paramMNT.n_layer      = imp_param.n_layer;
paramMNT.hs           = OUT_IR.config_dim(1,:);
paramMNT.hr1          = OUT_IR.config_dim(2,:);
paramMNT.hr2          = OUT_IR.config_dim(3,:);
paramMNT.dsr          = OUT_IR.config_dim(4,:);
paramMNT.n_config     = length(paramMNT.hs);
paramMNT.flow         = imp_param.flow;
paramMNT.fhigh        = imp_param.fhigh;
paramMNT.DeltaL_meas  = OUT_IR.DeltaL_meas;
paramMNT.f_resolution = 2;

% if we allow slight flexibility in the physical configurations, we will
% have 4 more search parameters. guess_param, lowlim and highlim are then
% adjusted accordingly.
if param.var_dim == 1
    guess_param = [guess_param 0.01 0.01 0.01 0.01];
    lowlim      = [lowlim -0.05 -0.05 -0.05 -0.1];
    highlim     = [highlim 0.05 0.05 0.05 0.1];
end

% we set some optimization parameters for patternsearch
options = psoptimset('MaxIter',10000,'MaxFunEvals',1000000, ...
    'CompletePoll','on','Display','final','Vectorized','off', ...
    'TolX',1e-6,'TolFun',1e-6);

%% this is where the minimization routine takes place.

% our minimization function, func is basically the accumulated error when
% one compares the predicted and measured level differences, as defined
% below.
func = @(X) sum(sum(abs(OUT_IR.DeltaL_meas-solveMNT(X,paramMNT)))) ...
            /paramMNT.n_config;

% then we do the search
[guessed_param,fval,exitflag,output]= patternsearch(func,guess_param, ...
    [],[],[],[],lowlim, ...
    highlim,[],options);

% and output the result when completed
disp(['total error: ', num2str(fval), ...
    ' air flow resistivity: ',num2str(guessed_param(1)/1000), ...
    ' kPas/m2 and porosity: ', num2str(guessed_param(2)), ...
    ' and thickness: ', num2str(guessed_param(3)*1000),'mm'])

return
%%
% make some comparison plot if necessary
[DeltaL_pred,extraMNTout] = solveMNT(sigma,paramMNT,opt);

figure
for i = 1:3
    subplot(3,1,i)
    semilogx(OUT_IR.foct,DeltaL_pred(:,i),'r--d')
    hold on
    semilogx(OUT_IR.foct,OUT_IR.DeltaL_meas(:,i),'k-*')
     legend('pred','meas')
     grid on
    xlabel('Frequency (Hz)')
    ylabel('Level difference (dB)')
end


