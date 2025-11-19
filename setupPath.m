%Add Paths for simulations

%Get absolute path to the folder containing this file
basePath = [fileparts(mfilename('fullpath')) filesep];

%Add paths
addpath([basePath '/BiGAMP']) %BiGAMP code
addpath([basePath '/supportFunction']) %support functions
addpath([basePath '/supportClass']) %support classes 
