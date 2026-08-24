function plot_iter2D(roacal,iter,rhoa)
%Plotting iteration Summary of this function goes here
%   Detailed explanation goes here


rms_error = norm(roacal-rhoa)/sqrt(length(rhoa));
%title(['\bf \fontsize{12}\fontname{Times}RESPONS - LM || RMS : ', num2str(rms_error),' || iter : ', num2str(iter)]);

plot(iter,rms_error)


end

