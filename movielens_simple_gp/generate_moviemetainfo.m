% Download Movielens 100K from 
% http://grouplens.org/datasets/movielens/
% Extract the data in a folder.
% This function extracts the following
% (1) Item metadata, title, release date.
% (2) Ratings and timestamps
% There are 943 users and 1682 movies
% See README for details
%
% Written by Emtiyaz, EPFL, pl8787
% Modified on March 11, 2014
function generate_moviemetainfo(item_file, save_file)
    % movie information
    fileName = item_file;
    out = importdata(fileName,'|');
    data = out.data;

    movieInformation = out.textdata; % title, release date etc
    movieMetaData = sparse(data); % genres information
    
    save(save_file, 'movieMetaData');
end