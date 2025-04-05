setpathsloc = '/home/dcuric/Documents/WFOIanalysisfunctions/';
run([setpathsloc 'setpaths.m'])
dorsalMapLoc = [setpathsloc 'Auxillary/allenDorsalMap_donovan.mat'];
load([setpathsloc 'Auxillary/allenDorsalMap_donovan.mat'],'dorsalMaps');
clear dorsalMapLoc