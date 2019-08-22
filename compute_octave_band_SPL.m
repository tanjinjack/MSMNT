function [fout,SPL] = compute_octave_band_SPL(inpfreq,inpSPL,fl,fh,bandtype)
%generate octave values between lower bound "fl" and upper bound "fh" for a 
%vector of "inpSPL" with corresponding frequency of "inpfreq".
%input '1/3', '1/2' or '1' to choose octave band type. As of 8/8/19, only
%'1/3' is implemented.
%the octave band is calculated based on essentially rectangular window
%type, instead of any sort of bandpass filter.
%note default ref level is used as 2e-5 here.
cf_use = [];
reflevel = 2e-5;

if (bandtype == '1/3' | bandtype == 1/3)
    cf_all = [20 25 31.5 40 50 63 80 100 125 160 200 250 315 400 500 ...
        630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000 ...
        12500 16000 20000];
    cf_use = cf_all(cf_all >= fl & cf_all <= fh);
    cf_low = 2^(-1/6)*cf_use;
    cf_high= 2^(+1/6)*cf_use;
 
elseif (bandtype == '1' | bandtype == 1)
    disp('One octave band is not implemented yet')
    return
elseif (bandtype == '1/2' | bandtype == 1/2)
    disp('Half octave band is not implemented yet') 
    return
else
    disp('Please input only ''1/3'', ''1/2'' or ''1''')
    return
end

if isempty(cf_use)
    disp('No center frequency is being used. Check your freq bounds.')
end

for i = 1:length(cf_use)
    i_fl = find(inpfreq >= cf_low(i),1);
    i_fh = find(inpfreq < cf_high(i),1,'last');
    %formula based on "octave" on MATLAB file exchange
    %https://nl.mathworks.com/matlabcentral/fileexchange/69-octave
    RMS_SPL = sqrt(sum(((inpSPL(i_fl:i_fh)).^2)/(i_fh-i_fl+1)));
    SPL(i) = 20*log10(RMS_SPL/reflevel);
end

fout = cf_use;

end