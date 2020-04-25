
function [area_roc, false_alarm, hit, dprime, u,s] = auroc (a_dist, b_dist, binsize, titlestr)

minvalue = min(min(a_dist),min(b_dist));
maxvalue = max(max(a_dist),max(b_dist));
if binsize == 0
    num_bins = 100;
    binsize = abs(((maxvalue+abs(maxvalue))-(minvalue-abs(minvalue))))/num_bins;
end;
xvalue = (minvalue-abs(minvalue)):binsize:(maxvalue+abs(maxvalue));

%ROC on the distributions
zi = 1;
%if mean(a_dist)<=mean(b_dist)
    noise_dist = a_dist;
    signal_dist = b_dist;
%else
%    noise_dist = b_dist;
%    signal_dist = a_dist;
%end;

num_bins1 = round((max(noise_dist)-min(noise_dist))/binsize);
num_bins2 = round((max(signal_dist)-min(signal_dist))/binsize);
 

%   figure, histfit(signal_dist, num_bins2), hold on,
%   h = get(gca,'Children');
%   set(h(2),'FaceColor','r'); set(h(1), 'Color','r'), hold on, %alpha(0.5), hold on, 
%   h = findobj(gca,'Type','patch'), hold on,
%   set(h,'FaceColor','r'), hold on, %'EdgeColor','w', 'AlphaData', 0.7 * get(h,'AlphaData'))
%  % plot(xvalue,normpdf(xvalue, mean(signal_dist), std(signal_dist)), 'r', 'LineWidth', 3), hold on,
%   
%   histfit(noise_dist, num_bins1), hold on,
%   h = get(gca,'Children');
%   set(h(2),'FaceColor','b'); set(h(1), 'Color','b'), hold on, %alpha(0.5), hold on, 
%   
%   title(titlestr);

for z = (minvalue-abs(minvalue)):binsize:(maxvalue+abs(maxvalue))
    false_alarm(zi) = length(nonzeros(noise_dist>=z))/length(noise_dist);
    hit(zi) = length(nonzeros(signal_dist>=z))/length(signal_dist);
    zi=zi+1;
end;

%figure, plot(false_alarm, hit, 'LineWidth', 3), hold on,
if maxvalue == minvalue
    area_roc = -0.1;
else
area_roc = trapz(fliplr(false_alarm), fliplr(hit));
end;
u = abs(mean(signal_dist)-mean(noise_dist));
s = sqrt((std(signal_dist)^2+std(noise_dist)^2)/2);
dprime = u/s;
%dprime = abs(mean(signal_dist)-mean(noise_dist))/std(signal_dist)/std(noise_dist);
