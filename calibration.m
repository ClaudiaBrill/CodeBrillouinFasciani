function [calibration_fit, par]=calibration(y_calib, threshold, distance_btw_peaks)
%this program calculates the parameters for px-to-GHz conversion.
cm2Ghz = 29.9702547; %%VIPA FSR

xlim1 = 1;
xlim2 = length(y_calib);
[y_peaks, x_peaks] = find_peaks(y_calib(xlim1:xlim2),xlim1:xlim2, distance_btw_peaks, threshold);



%%
nPeaks=length(x_peaks);
sigma=5;
[par_spectrum, ~] = fitting_part(xlim1:xlim2, y_calib(xlim1:xlim2), y_peaks, repmat(sigma, 1,nPeaks), x_peaks);

a_el = par_spectrum(1:nPeaks);                    %only heights of the peaks
c_el = par_spectrum(2*nPeaks+1:3*nPeaks);         %only x of the peaks

%%this automatically detects the reference Rayleigh peak
[~,index_first_max] = max(a_el);
index_second_max = index_first_max - 3;  %%%3 is needed to shift in on the left (but there are 2 brillouins between)
if index_second_max > 0 && isempty(find(a_el > a_el(index_second_max) & a_el < a_el(index_first_max), 1 ))
    peak_reference = min([index_first_max, index_second_max]);
else
    peak_reference = index_first_max;
end

deltapixel_calib = c_el - c_el(peak_reference);

shift_b = 7.50; %%Brillouin shift of water in backscattering, lambda = 532 nm
deltaGHz_calib = @(par) calculate_GHz_vector(deltapixel_calib, peak_reference, shift_b, par(1));
starting_points = [cm2Ghz, 10^-5, 10^-1, 10^-10];

sy = cm2Ghz/deltapixel_calib(peak_reference+3); %%errors

calibration_fit = @(par, x) par(2)*x.^2 + par(3)*x + par(4);
weight = 1./repmat(sy,1, length(deltapixel_calib));

fun = @(par) weight.*(deltaGHz_calib(par) - (calibration_fit(par, deltapixel_calib)));

options = optimoptions(@lsqnonlin,'Display','iter','Algorithm','levenberg-marquardt',...
    'FiniteDifferenceType','central','FunctionTolerance',1.0e-16,'MaxIterations',10,'MaxFunctionEvaluations',1.0e8,'StepTolerance',1.0e-10);

[par,~,residual,~,~,~,jacobian] = lsqnonlin(fun, starting_points, [], [], options);
CI = nlparci(par,residual,'jacobian',jacobian,'alpha',0.32);                %CI = 68% confidence interval = 1 sigma



end





function deltaGHz = calculate_GHz_vector(deltapixel, peak_reference, shift_b, FSR)
%%This function creates automatically the vector in GHz from the reference
%%that is in pixels.

deltapix_brill = deltapixel(peak_reference+1);
deltapix_FSR = deltapixel(peak_reference+3);

Npeaks = length(deltapixel);
deltaGHz = NaN(1, Npeaks);

deltaGHz(peak_reference) = 0;
deltaGHz(peak_reference+1) = shift_b;
deltaGHz(peak_reference+2) = FSR - shift_b;
deltaGHz(peak_reference+3) = FSR;

tolerance = 10;
if peak_reference ~= 1
    if (abs(deltapixel(peak_reference-1)) > deltapix_brill-tolerance) && (abs(deltapixel(peak_reference-1)) < deltapix_brill+tolerance)
        deltaGHz(peak_reference-1) = -shift_b;
    else
        deltaGHz(peak_reference-1) = -FSR;
    end
end
if peak_reference+4 <= Npeaks
    if (deltapixel(peak_reference+4) > deltapix_brill+deltapix_FSR-tolerance) && (deltapixel(peak_reference+4) < deltapix_brill+deltapix_FSR+tolerance)
        deltaGHz(peak_reference+4) = FSR + shift_b;
    else
        deltaGHz(peak_reference+4) = 2*FSR;
    end
end

for k=1:Npeaks
    if isnan(deltaGHz(k))
        deltaGHz(k) = FSR * ceil((deltapixel(k))/deltapix_FSR);
    end
end

end
