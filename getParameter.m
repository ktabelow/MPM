% given a nii file, search for the parameter (TE or TR or FA) in the header
% 
% function [res] = getParameter(FilesCellArray, StringToLookFor)
%
% 2015/12/08
%
% written by C. D'Alonzo
% 
% search for the values of the parameter indicated in the StringToLookFor
% in the description of the nii file given as single filename or cell array
% of filenames
%
% Input:
% FilesCellArray    - a cell array containing one or more filenames of nii
%                     images
% StringToLookFor   - a string containing the parameter that has to be
%                     found ('TE' or 'TR' or 'FA')
% Output:
% res               - vector containing all the values of the parameter of
%                     interest (one for each file)
% 
% ========================================================================


function [res] = getParameter(FilesCellArray, StringToLookFor)

 V = spm_vol(FilesCellArray);
 if ~isstruct(V)
     res = zeros(1,length(V));
 
 for i=1:length(V)
     tokens = regexpi(V{i}.descrip, '([A-Z]{2})=([.0-9]+)([a-z]{2,})','tokens');
     for k=1:length(tokens)
         if strcmp(tokens{k}{1},StringToLookFor)
             res(i) = str2double(tokens{k}{2});
         end
     end
 end
 else
     tokens = regexpi(V.descrip, '([A-Z]{2})=([.0-9]+)([a-z]{2,})','tokens');
     for k=1:length(tokens)
         if strcmp(tokens{k}{1},StringToLookFor)
             res = str2double(tokens{k}{2});
         end
     end
 end
 
 
 
 

end