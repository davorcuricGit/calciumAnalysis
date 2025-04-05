function [validPixels, dmask] = get_universal_mask(setpathsloc);
%this mask is going to be same for each subject
dmask = getfield(load([setpathsloc '/maskDavor.mat'], 'Mask_Davor'), 'Mask_Davor');
dmask = single(dmask/255);
dmask(dmask ~= 1) = NaN;
validPixels = single(find(dmask == 1));
end