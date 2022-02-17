function wget_dload(wdir,fname)

%--------------------------------------------------------------------------
% Download files using wget
% Files will be saved in the directory wdir
% 
% USAGE: wget_dload(wdir,fnames);
% 
% INPUT:
% wdir = directory where to save files (string)
% fname = name of text file
% 
% OUTPUT:
% None
% 
% Cara Manning
% modified Aug 29, 2021
%--------------------------------------------------------------------------

cd(wdir);

%download the data through command line using wget
    %first, change system directory
       system(['cd ', wdir, '\']); 
    %now download
       disp(['downloading ', fname]);
       system(['wget --no-check-certificate -i', fname]); 
