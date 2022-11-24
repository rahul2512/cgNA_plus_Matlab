% ----------------
% Copyright 2015 Jaroslaw Glowacki,  Lennart de Bruin, Thomas Zwahlen
% (current)
% jarek (dot) glowacki (at) gmail (dot) com
%
% This file is part of cgDNArecon.
%
% cgDNArecon is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% cgDNArecon is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with cgDNArecon.  If not, see <http://www.gnu.org/licenses/>.
% ----------------
%
% The function saves a cgDNA parameter set given as argument in a text file
% formated for use with the cgDNAmc C++ Monte Carlo code.
%
% Arguments:
%     paramset    The cgDNA parameter set structure as defined by the cgDNA
%                 package
%     filename    Name (path) for the file to save the parameter set in

function Mat2Text_ps( paramset, filename )

number_format = '%16.12f ';

file_id = fopen( [ filename '.txt' ] , 'w');

stiff_int = paramset.stiff_int;
stiff_end5 = paramset.stiff_end5;
stiff_end3 = paramset.stiff_end3;

sigma_int = paramset.sigma_int;
sigma_end5 = paramset.sigma_end5;
sigma_end3 = paramset.sigma_end3;

Dimers = fieldnames(stiff_int);
% Write inner dimer data
for k = 1:length(Dimers)
    dimer = Dimers{k};
    fprintf(file_id, '%s\n', dimer);
    K = stiff_int.(dimer) ; 
    sigma = sigma_int.(dimer) ;
    for r = 1:42
        fprintf(file_id, number_format,K(r,:));
        fprintf(file_id, '\n');
    end
    fprintf(file_id, '\n');
    fprintf(file_id, number_format,sigma);
    fprintf(file_id, '\n\n');
    
end

SpecialDimers = {'MN','NM','AM','TM','GM','CM','NT','NA','NC','NG','HI','IH',...
    'AH','TH','GH','CH','IT','IA','IC','IG'} ;

for k = 1:length(SpecialDimers)
    dimer = SpecialDimers{k};
    fprintf(file_id, '%s\n', dimer);
    
    if isfield(stiff_int,dimer)
        K = stiff_int.(dimer) ;
        sigma = sigma_int.(dimer) ;
    else
        K = zeros(42,42) ;
        sigma = zeros(42,1) ;
    end
    for r = 1:42
        fprintf(file_id, number_format, K(r,:));
        fprintf(file_id, '\n');
    end
    fprintf(file_id, '\n');
    fprintf(file_id, number_format, sigma);
    fprintf(file_id, '\n\n');
end

% Write outer/front dimer data
for k = 1:length(Dimers) % only nonmethylated dimers appear at the ends
    dimer = Dimers{k};
    K = zeros(42,42) ;
    sigma = zeros(42,1) ;
    if isfield(stiff_end5,dimer)
        K(1:36,1:36) = stiff_end5.(dimer) ;
        sigma(1:36) = sigma_end5.(dimer) ;
    end
    
    fprintf(file_id, '%s\n', strcat(dimer,'_end5'));
    for r = 1:42
        fprintf(file_id, number_format, K(r, :));
        fprintf(file_id, '\n');
    end
    fprintf(file_id, '\n');
    fprintf(file_id, number_format, sigma);
    fprintf(file_id, '\n\n');
    
end

% Write outer/back dimer data
for k = 1:length(Dimers) % only nonmethylated dimers appear at the ends
    dimer = Dimers{k};
    K = zeros(42,42) ;
    sigma = zeros(42,1) ;
    if isfield(stiff_end3,dimer)
        K(end-35:end,end-35:end) = stiff_end3.(dimer) ;
        sigma(end-35:end) = sigma_end3.(dimer) ;
    end
    
    fprintf(file_id, '%s\n', strcat(dimer,'_end3'));
    for r = 1:42
        fprintf(file_id, number_format, K(r, :));
        fprintf(file_id, '\n');
    end
    fprintf(file_id, '\n');
    fprintf(file_id, number_format, sigma);
    fprintf(file_id, '\n\n');
    
end

fclose(file_id);

end

