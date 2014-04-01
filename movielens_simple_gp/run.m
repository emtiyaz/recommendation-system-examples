% A simple recommendation system
% Written by Emtiyaz, EPFL, pl8787
% Modified on March 11, 2014

movielensdata;

fprintf('Dataset U1\n');
recommend0('movielens100k1','moviemetainfo');

fprintf('Dataset U2\n');
recommend0('movielens100k2','moviemetainfo');

fprintf('Dataset U3\n');
recommend0('movielens100k3','moviemetainfo');

fprintf('Dataset U4\n');
recommend0('movielens100k4','moviemetainfo');

fprintf('Dataset U5\n');
recommend0('movielens100k5','moviemetainfo');
