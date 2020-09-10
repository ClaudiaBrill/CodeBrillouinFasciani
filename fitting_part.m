function [par, CI] = fitting_part(xaxis, yaxis, amplitude, sigma, central_position)
[xData, yData] = prepareCurveData( xaxis, yaxis );
nPeaks=length(central_position);
b = sigma/2;
c = reshape(central_position, [1,nPeaks]);
off = 180; %%camera offset
a = reshape(amplitude, [1,nPeaks])- off;

starting_points = [a, b, c, off];
sy = sqrt(abs(yData'-off)) ;  %%Poissonian error
f_fit = @(par) SumOfLorentzian(par, xData) + par(end);

lb = [];
ub = [];
options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt', 'FiniteDifferenceType','central','FunctionTolerance',1.0e-16,'MaxIterations',1e2,'MaxFunctionEvaluations',1.0e8,'StepTolerance',1.0e-16);
weight =  1./sy;
fun   = @(par) weight.*(yData' - reshape(f_fit(par), [1, length(yData)]));

[par,~,residual,~,~,~,jacobian] = lsqnonlin(fun,starting_points,lb,ub,options);
CI = nlparci(par,residual,'jacobian',jacobian,'alpha',0.32);                %CI = 68% confidence interval
chi_square= sum(fun(par).^2)/(length(yData)-length(par));

    function res = SumOfLorentzian(par, x)
        summm = 0.0 ;
        nPeaks = (length(par)-1)/3; %%%sum of N independent lorentzians
        a = par(1:nPeaks);
        b = par((nPeaks+1):2*nPeaks);
        c = par((2*nPeaks+1):3*nPeaks);
        
        L = @(a,b,c) a*b^2./((x-c).^2+b^2);
        
        for i=1:nPeaks
            summm = summm + L(a(i),b(i),c(i));
        end
        res = summm  ;
    end

%%this should be commented to avoid hundreds of graphs while fitting all the data :)
figure
errorbar(xData, yData, sy, 'o');
hold on
plot(xData, f_fit(par), 'LineWidth', 1.5);
plot(xData, f_fit(starting_points), '--', 'LineWidth', 1.5);
end
