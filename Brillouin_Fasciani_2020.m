%%%%CODE FOR BRILLOUIN SPECTRA MEASUREMENTS%%
%%Used for Brillouin Microscopy measurements in the paper Fasciani et al., to appear in Nature Genetics
%%%Claudia Testi, PhD
%%Post-doc @ Istituto Italiano di Tecnologia
%%claudia.testi@iit.it
%%Please cite if used!

%%Done in Matlab 2017a
%%The input is y, the matrix of the sample containing per each point a vector of Brillouin spectra 
%%(Rayleigh, Brillouin_Stokes, Brillouin_Anti-Stokes, Rayleigh).
%%y_calib is the calibration curve for px-to-GHz conversion on a liquid whose Brillouin shift is
%%known (in our case, water = 7.50 GHz in backscattering, lambda = 532 nm)

close all
clear
clc
cm2Ghz = 29.9702547; %%VIPA FSR


load('data.mat'); %%y cell matrix
load('y_calib.mat'); %%calibration curve vector

distance_btw_peaks = 30;
threshold_for_peaks = 0.01;
[calibration_fit, par_calibration] = calibration(y_calib, threshold_for_peaks, distance_btw_peaks); %this program calculates the parameters for px-to-GHz conversion


last_x=size(y,1); %%x size of the sample
last_y=size(y,2); %%y size of the sample

xlim1 = 1;
xlim2 = size(y{1,1},1);
n_points = 60;
xbrill1 = round(xlim2/2 - n_points);
xbrill2 = round(xlim2/2 + n_points);
elastic_range_px=[(xlim1:xbrill1),(xbrill2:xlim2)];
brill_range_px = xbrill1:xbrill2;

distance_btw_peaks_rayleigh = 150;
distance_btw_peaks_brill = 70;
threshold_findpeaks_el = 0.5;
threshold_findpeaks_brill = 0.2;


xxx = xlim1:xlim2;
c_brill_left = 12;
c_brill_right = 14;
shift_sample = 7.50;

db_left = zeros(last_x,last_y); db_right = zeros(last_x,last_y);
s_db_left = zeros(last_x,last_y); s_db_right = zeros(last_x,last_y);
db_wm = zeros(last_x,last_y); s_db_wm = zeros(last_x,last_y);

max_value = 1*10^4; %%%if Rayleigh peak > max_value, skip it.

for iii=1:last_x
    for jjj=1:last_y
        if max(y{iii,jjj}) > max_value   %%%if Rayleigh peak > max_value, skip it.
            continue
        end
        [elastic_peaks_ampl, elastic_peaks_position] = find_peaks(y{iii,jjj}(elastic_range_px), elastic_range_px, distance_btw_peaks_rayleigh, threshold_findpeaks_el);
        [brill_peaks_ampl, ~] = find_peaks(y{iii,jjj}(brill_range_px), brill_range_px, distance_btw_peaks_brill, threshold_findpeaks_brill);
        
        %%Fit for pixel->GHz conversion
        amplitude_0 = elastic_peaks_ampl;
        sigma_0 = [6, 6];
        xposition_0 = elastic_peaks_position;
        par_pixel = fitting_part(xxx, y{iii,jjj}(xxx), amplitude_0, sigma_0, xposition_0);
        nPeaks = length(amplitude_0);
        c_pixel = par_pixel(2*nPeaks+1);
        x_data_GHz_all = calibration_fit(par_calibration, ((1:length(y{1,1}))- c_pixel));
        
        %%Final fit in GHz
        amplitude_0 = [elastic_peaks_ampl(1), brill_peaks_ampl(1), 190, brill_peaks_ampl(2), elastic_peaks_ampl(2)];
        xposition_0 = [0, shift_sample, 15, cm2Ghz-shift_sample, cm2Ghz];
        x_data_GHz= x_data_GHz_all(xxx);
        sigma_0 = [1, 1.5, 10, 1.5, 1];
        [par_all, CI] = fitting_part(x_data_GHz, y{iii,jjj}(xxx), amplitude_0, sigma_0, xposition_0);
        
        db_left(iii,jjj) = par_all(c_brill_left)-par_all(c_brill_left-1);
        db_right(iii,jjj) = par_all(c_brill_right+1)-par_all(c_brill_right);
        s_db_left(iii,jjj) = (CI(c_brill_left,2)-CI(c_brill_left,1))/2;
        s_db_right(iii,jjj) = (CI(c_brill_right,2)-CI(c_brill_right,1))/2;
        [db_wm(iii,jjj), s_db_wm(iii,jjj)] = weighted_mean(par_all(c_brill_left)-par_all(c_brill_left-1), par_all(c_brill_right+1)-par_all(c_brill_right), s_db_left(iii, jjj), s_db_right(iii, jjj));
    end
end


c_axis_limits = [7.45, 8]; %Brillouin map range colorbar, in GHz

brillouinmap = db_wm;
figure
imagesc(db_wm); colorbar; colormap(jet(100)); pbaspect([last_x last_y 1])
caxis([c_axis_limits(1) c_axis_limits(2)])

%%%this fills eventual points that have saturated (their db_shift is 0)
threshold = 0;
imFilt = smoothWithNeighbours2(brillouinmap,[1,last_x],[1,last_y],threshold);
figure
imagesc(imFilt); colorbar; colormap(jet(100)); pbaspect([last_x last_y 1])
caxis([c_axis_limits(1) c_axis_limits(2)])


%%2x resizing, 3x3 mean filter
finalmap = calculate_final_map(imFilt, [1,last_x],[1,last_y], c_axis_limits);
figure
imagesc(finalmap); colorbar; colormap(jet(100)); pbaspect([last_x last_y 1])
caxis([c_axis_limits(1) c_axis_limits(2)])


%%This segments different regions of the maps. 
%%We consider a “nucleus” all the pixels above 7.70 GHz; the program then connects different parts of the images 
%%having the same shift
threshold_nucleus = 7.70;
nucleus_cropped_shifts = nucleus_contour2(finalmap, threshold_nucleus, c_axis_limits);


%Final results of Brillouin shift +- SD
mean_db_nucleus = mean(nucleus_cropped_shifts);
std_db_nucleus = std(nucleus_cropped_shifts);







%%

function [wm, s_wm] = weighted_mean(x1, x2, s1, s2)
%%weighted mean between 2 values, using as weights their SD.
weight = 1./[s1, s2].^2;
wm = (x1*weight(1) + x2*weight(2))/sum(weight);
s_wm = sqrt(1/sum(weight));
end


function finalmap = calculate_final_map(im, nii, njj, c_axis_limits)
ratio = 2 ;
h = 1/3*ones(3,1);
H = h*h';

im_resize2x = imresize(im,ratio,'bilinear');
im2 = filter2(H,im_resize2x);

finalmap = im2(nii(1)+1:ratio*nii(2)-1,njj(1)+1:ratio*njj(2)-1);

figure
imagesc(finalmap)
colorbar; colormap(jet(100)); pbaspect([size(finalmap,1) size(finalmap,2) 1])
caxis([c_axis_limits(1), c_axis_limits(2)])
end


function nucleus_cropped_shifts = nucleus_contour2(to_analyze, threshold_nucleus, axis_values)
mask = to_analyze >= threshold_nucleus;

disk_to_close = 2;
mask_closed = imclose(mask, strel('diamond', disk_to_close));

area_to_remove = 100;
mask_cleaned = bwareaopen(mask_closed, area_to_remove);

figure
subplot(1,3,1); imshow(mask)
subplot(1,3,2); imshow(mask_closed)
subplot(1,3,3); imshow(mask_cleaned)

nucleus_cropped_shifts = to_analyze(mask_cleaned);

figure
imagesc(to_analyze)
colorbar; colormap(jet(100)); pbaspect([size(to_analyze,1) size(to_analyze,2) 1])
caxis([axis_values(1) axis_values(2)])
hold on
contour(mask_cleaned, 1, 'LineWidth', 1);

end


