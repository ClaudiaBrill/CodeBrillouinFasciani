function [amplitude, position] = find_peaks(yy, xx, distance_btw_peaks, threshold)
funct = (yy-min(yy))/max(yy-min(yy));

[amplitude, position] = findpeaks(funct, xx, 'MinPeakHeight', threshold, 'MinPeakDistance', distance_btw_peaks);
amplitude = amplitude.*max(yy-min(yy)) + min(yy);
end
