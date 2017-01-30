function err = update_aacgmv2coefs( fdir )
% 
% Update the spherical harmonic coefficients file
% aacgmv2coefs.mat. At time of writing this function, the latest
% coefficients are availble in
% https://engineering.dartmouth.edu/superdarn/aacgm.html 
% Download the package aacgm_coeffs-12.tar and untar its
% contents.
%
% The updated coefficient file aacgmv2coefs.mat is written in the
% current working directory. 
% 
% INPUT:
%  fdir   path to the directory with aacgm coefficient files. The
%         file names MUST end with yyyy.asc, where yyyy is year,
%         and the directory fdir MUST NOT contain any other files
%         than valid coefficient files. 
% 
% OUTPUT:
%  err    0 if the file was successfully created
% 
% 
% IV 2016
% 

err = 1;

NCOORD =  3; %/* xyz */
POLYORD =  5; %/* quartic polynomial fit in altitude */
NFLAG = 2; %/* 0: geo to AACGM, 1: AACGM to geo */
SHORDER = 10; %/* order of Spherical Harmonic expansion */
AACGM_KMAX = ((SHORDER+1)*(SHORDER+1));   %/* number of SH coefficients */
MAXALT = 2000; % maximum altitude 
RE = 6371.2; % Earth radius in km


% list coefficient files
coeffiles = dir(fullfile(fdir,'*.asc'));

% number of files
nf = length(coeffiles);

if nf==0
    error('no index files found')
end

% allocate the table
aacgmv2coefs = zeros(nf,AACGM_KMAX,NCOORD,POLYORD,NFLAG);

% a vector for the years
aacgmv2years = zeros(nf,1);

% read the values from files
for y=1:nf
    fid = fopen(fullfile(fdir,coeffiles(y).name),'r');
    dd = fscanf(fid,'%lf');
    aacgmv2coefs(y,:,:,:,:) = reshape(dd,[AACGM_KMAX,NCOORD, ...
                        POLYORD,NFLAG]);
    aacgmv2years(y) = str2num(coeffiles(y).name(end-7:end-4));
end

save('aacgmv2coefs.mat','aacgmv2coefs','aacgmv2years','AACGM_KMAX','SHORDER','NFLAG','POLYORD','NCOORD','MAXALT','RE');

err = 0;

end
