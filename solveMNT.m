function [DeltaL_pred,extraOutput] = solveMNT(inp,param)
% OUTPUT:
% DeltaL_pred, the theoretical level difference based on the geometry 
%              and surface impedance of the ground. This is a key
%              output for the minimization routine.
%
% extraOutput contains the following:
% extraOutput.Impedance, gives the surface impedance based on the input 
%                        material parameters.
% extraOutput.alpha    , gives the surface absorption coefficient based
%                        on the input material parameters.
% extraOutput.freq     , gives the frequency values used to represent
%                        extraOutput.Impedance and extraOutput.alpha.
%
%
% INPUT:
% inp  , is a vector of variable length, depending on the case. 
%
%          For n_layer = 0
%              inp(1) represents flow resistivity,
%              inp(2) represents porosity,
%          if we allow slight adjustment of the configuration parameters
%          i.e. param.var_dim == 1
%          input(3:6) represents the slight variance for height of the
%          source, height of mic 1 and 2, and distance between the 
%          source and the mic.
%
%          For n_layer = 1
%              inp(1) represents flow resistivity,
%              inp(2) represents porosity,
%              inp(3) represents the material thickness 
%          if we allow slight adjustment of the configuration parameters
%          i.e. param.var_dim == 1
%          input(4:7) represents the slight variance for height of the
%          source, height of mic 1 and 2, and distance between the 
%          source and the mic.
%
%
% param contains several parameters"
%
% param.n_layer      , same as imp_param.n_layer, defines the number of
%                      layer of the material in concern.
% param.hs           , as extracted in processImpulseResponse,() height 
%                      of source.
% param.hr1          , as extracted in processImpulseResponse, height 
%                      of microphone 1.
% param.hr2          , as extracted in processImpulseResponse, height
%                      of microphone 2.
% param.dsr          , as extracted in processImpulseResponse, distance 
%                      between the source and the microphones.
% param.n_config     , number of the configurations available, ideally 
%                      to be equaled to number of material properties 
%                      to be solved.
% param.flow         , = imp_param.flow, see processImpulseResponse().
% param.fhigh        , = imp_param.fhigh, see processImpulseResponse().
% param.DeltaL_meas  , = OUT_IR.DeltaL_meas, see processImpulseResponse().
% param.f_resolution , resolution of the frequency vector for the 
%                      computation. Smaller number results in better 
%                      accuracy but longer computation time.
% param.var_dim      , allows adjustment of configuration parameters via
%                      patternsearch(). Unique of result cannot be
%                      guaranteed.

f_resolution = param.f_resolution;
flo_resist   = inp(1);
porosity     = inp(2);

if param.n_layer == 1
     thick1 = inp(3);   
     if param.var_dim == 1
        hs_diff  = inp(4); %variability of source height.
        hr1_diff = inp(5); %variability of mic1 height.
        hr2_diff = inp(6); %variability of mic2 height
        dsr_diff = inp(7); %variability of distance between source and mic.
    else
        hs_diff  = 0;
        hr1_diff = 0;
        hr2_diff = 0;
        dsr_diff = 0;
    end       
     
elseif param.n_layer == 0
    thick1 = 0;
    if param.var_dim == 1
        hs_diff  = inp(3);
        hr1_diff = inp(4);
        hr2_diff = inp(5);
        dsr_diff = inp(6);
    else
        hs_diff  = 0;
        hr1_diff = 0;
        hr2_diff = 0;
        dsr_diff = 0;
    end
end

n_config     = param.n_config;
P0           = 1.42e5; % adiabatic bulk modulus
Pr           = 0.713; % Prandtl number for air
gamma        = 1.4; %specific heat ratio
temperature  = 19.6; %temperature [C]
c1           = 331.3+0.606*temperature;  % speed of sound [m/s]
rho_air      = 1.2;  % density of air [kg/m^3]
omega        = 2*pi*(param.flow:f_resolution:param.fhigh); %angular frequency
k1           = omega/c1; % wave number of air

%% the slit-pore model 
% see Chang's thesis pg 13 (note error in Eq 1.46
% or Attenborough, Bashir and Taherzadeh "Outdoor ground impedance models"
% JASA 129(5), 2011.
Tort     = 1/sqrt(porosity); %Tortuosity
lambda   = sqrt((3*rho_air *omega *Tort)/(porosity*flo_resist));
lnpr     = lambda.*sqrt(Pr);
Gs_lamb  = 1 - tanh(lambda*sqrt(-1i))./(lambda*sqrt(-1i));
Gs_lnpr  = 1 - tanh(lnpr*sqrt(-1i))./(lnpr*sqrt(-1i));
rho_lamb = rho_air ./ Gs_lamb;
C_lamb   = (gamma-(gamma-1)*Gs_lnpr)/(gamma*P0);
K        = omega.*sqrt(Tort * rho_lamb .* C_lamb );
Zc       = sqrt(Tort/porosity^2 * rho_lamb ./ C_lamb)/(rho_air * c1);
N        = K./k1;
M        = rho_air./rho_lamb;
% M        = 1./(N.*Zc); %alternative formulation

hs  = param.hs + hs_diff - thick1;  %corrected height of source
hr1 = param.hr1 + hr1_diff - thick1; %corrected height of mic 1
hr2 = param.hr2 + hr2_diff - thick1; %corrected height of mic 2
dsr = param.dsr + dsr_diff; %distance between source and mic


%corrected R1 and R2
R1{1}    = sqrt( ( hs -  hr1).^2 +  dsr.^2 ); %direct path
R2{1}    = sqrt( ( hs +  hr1).^2 +  dsr.^2 ); %reflected path
Theta{1} = atan(  dsr./( hs +  hr1) ); %reflection angle
R1{2}    = sqrt( ( hs -  hr2).^2 +  dsr.^2 ); %direct path
R2{2}    = sqrt( ( hs+  hr2).^2 +  dsr.^2 ); %reflected path
Theta{2} = atan(  dsr./( hs +  hr2) ); %reflection angle
 
for j = 1:2 %total number of mic
    
    % calculate Beta_e (effective acoustic admittance)
    coef = sqrt( N.'.^2 - sin(Theta{j}).^2 );
    if param.n_layer == 0
        Be{j} = M.'.*coef; %for semi-infinite ground
    elseif param.n_layer == 1
        Be{j} = (-1i * M).'.* coef .* tan( k1.'.* thick1 .* coef ); %for hard-backed layer
    else
        Be{j} = 0;
        disp('select correct n_layer')
        % NOTE: 2 layers not implemented yet.
        return
    end
    
    %calculate Q - spherical wave reflection coefficient
    % Salomons, Computational atmospheric acoustics, Springer, 2001, pg 134
    % or Li, Waters-Fuller, Attenborough "Sound propagation from a point
    % source over extended-rection ground", JASA 104(2), 1998
    % but first we need Rp and d.
    Rp = (cos(Theta{j}) - Be{j})./(cos(Theta{j}) + Be{j});
    d  = sqrt( 1i*k1.'.*R2{j}/2 ).* (Be{j} + cos(Theta{j}));
     
    % NOTE: There are a few ways to calculate the F(d) function. 
    % Three ways are provided here but by experience, the Attenborough
    % implementation appears to be the fastest and most reliable. Chang's
    % implementation is a variant of Salomons' implementation.
    Fd = compute_Fd_Attenborough(d);
%     Fd = compute_Fd_Salomons(d);
%     Fd = compute_Fd_Chang(d);
    
    Q{j} = Rp + (1-Rp).*Fd;
    
    % and we can calculate the theoretical pressure level of each mic.
    % note that there isn't a multiplier of 1/(4pi) in this case since we
    % are interested only in the level difference.
    L_pred_freq{j} = abs(   exp( 1i.*k1.'.*R1{j} )./R1{j} ...
                    + exp( 1i.*k1.'.*R2{j} ).* Q{j}./R2{j}   ) ;
    
    % we convert the difference into octave band
    for i = 1:n_config                
    [foct,L_pred{j}(:,i)] = compute_octave_band_SPL(omega/2/pi,L_pred_freq{j}(:,i), ...
                                    param.flow,param.fhigh,'1/3');
    end
end

% and finally we substract mic 1 from mic 2.
DeltaL_pred = L_pred{1} - L_pred{2};
        
% we also prepare some extra outputs
Impedance             = Zc.*coth(-1i*K*thick1);   
extraOutput.Impedance = Impedance;
extraOutput.alpha     = 1-abs(((Impedance-1)./(Impedance+1))).^2;
extraOutput.freq      = omega/2/pi;

end

function Fd = compute_Fd_Salomons(d)
    FdA = zeros(size(d));
    FdB = zeros(size(d));
    FdC = zeros(size(d));
%for first criterion
    critA = (abs(d.^2) < 8);
    if (any(any(critA)))
        summand = 0;
        i_iter = 1;
        temp = 1;
        while max(max(abs(temp))) > 1e-6
            temp =(d.^2).^i_iter./factorial(i_iter)./(2*i_iter+1);
            i_iter = 1 + i_iter;
            summand = summand + temp;
            if i_iter > 500
                disp('500 iteration reached when trying to compute erfc')
                break
            end
        end
        FdA = 1 + 1i*d*sqrt(pi).*exp(-d.^2) .* (1+ 2*1i*d/sqrt(pi).* summand);
    end

    %for second criterion
    critB = (d>=0 & abs(d.^2) >= 8);
    if (any(any(critB)))
        summand = 0;
        for i_iter = 1:8
            temp = (-1)^i_iter*(2*i_iter-1)./(2*(-1i*d).^2).^i_iter;
            summand = summand + temp;
        end
        FdB = - summand;

    end
    
    %for third criterion
    critC = (d<0 & abs(d.^2) >= 8);
    if (any(any(critC)))
        summand = 0;
        for i_iter = 1:8
            temp = (2*i_iter-1)./(2*(d).^2).^i_iter;
            summand = summand + temp;
        end            
        Heav = imag(d) > 0;
        FdC = 2*1i*d*sqrt(pi).*exp(-d.^2).*Heav - summand;

    end
    
    FdA(isnan(FdA)) = 0;
    FdB(isnan(FdB)) = 0;
    FdC(isnan(FdC)) = 0;
    
    Fd = FdA .* critA + FdB .* critB + FdC .* critC;
end


function Fd = compute_Fd_Attenborough(w)
    
    FwA = zeros(size(w));
    FwB = zeros(size(w));

    critA = (abs(w)<1);
    if any(any(critA))
        FwA = 1 + 1i * sqrt(pi) * w .* exp(-w.^2);
    end

    critB = (abs(w)>=1);
    if any(any(critB))
        HImw = -imag(w) > 0;
        FwB = 2*1i * sqrt(pi) * w .* exp(-w.^2) .* HImw - 1./(2*w.^2);
    end

    Fd = FwA .* critA + FwB .* critB ; 

end


function Fd = compute_Fd_Chang(d)
    FdA = zeros(size(d));
    FdB = zeros(size(d));
    
    critA = abs(d.^2) > 8 ; 
    if any(any(critA))
        summand = 0;
        for i_iter = 1:8
            temp = (2*i_iter-1)./(2*(d).^2).^i_iter;
            summand = summand + temp;
        end            
        FdA = 2*1i*d*sqrt(pi).*exp(-d.^2).*summand;   
    end
    
    critB = abs(d.^2) <= 8 ;
    if any(any(critB))
        summand = 0;
        for i_iter = 1:20
            temp = (d.^2).^i_iter/factorial(i_iter)/(2*i_iter+1);
            summand = summand + temp;
        end
        FdB = 1 + 1i*d*sqrt(pi).*exp(-d.^2) .* (1 + 2*1i/sqrt(pi) * d .* summand);
    end
    FdA(isnan(FdA)) = 0;
    FdB(isnan(FdB)) = 0;
    
    Fd = FdA .* critA + FdB .* critB;
    
end
