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

clear all

generate_ratings('..\\ml-100k\\u1.base','..\\ml-100k\\u1.test','movielens100k1');
generate_ratings('..\\ml-100k\\u2.base','..\\ml-100k\\u2.test','movielens100k2');
generate_ratings('..\\ml-100k\\u3.base','..\\ml-100k\\u3.test','movielens100k3');
generate_ratings('..\\ml-100k\\u4.base','..\\ml-100k\\u4.test','movielens100k4');
generate_ratings('..\\ml-100k\\u5.base','..\\ml-100k\\u5.test','movielens100k5');

generate_moviemetainfo('..\\ml-100k\\u.item','moviemetainfo');
