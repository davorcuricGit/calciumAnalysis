function mask = get_universal_mask(recLoc, dorsalMaps);
load([ recLoc '/transform.mat'])
mask = output_warped(mask,tform,dorsalMaps);
end