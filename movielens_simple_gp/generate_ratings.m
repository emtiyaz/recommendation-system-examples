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
function generate_ratings(train_file, test_file, save_file)
    nUsers = 943;
    nMovies = 1682;
    % read user ratings in a sparse matrix (nUser X nMovies)
    fileName = train_file;
    data = importdata(fileName);
    % user ratings
    ratings_train = sparse(data(:,2), data(:,1), data(:,3), nMovies, nUsers);

    fileName = test_file;
    data = importdata(fileName);
    % user ratings
    ratings_test = sparse(data(:,2), data(:,1), data(:,3), nMovies, nUsers);
    
    save(save_file, 'ratings_train', 'ratings_test');
end