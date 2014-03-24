function array2table(fileName, mean_, std_, rowLabels, colLabels)
% prints mean and standarad deviation in a latex table form

  [nRows, nCols] = size(mean_);

  filename = [fileName '.tex']; 
  filename
  fid = fopen(filename,'w');

  % get the {c|c|...|c} 
  cols = '{c';
  for i = 1:nCols
    cols = sprintf('%s|c', cols);
  end
  cols = [cols '}'];

  % first line of table
  fprintf(fid, '\\begin{tabular}[!h]%s\n',cols);
  % print column names
  for j = 1:nCols
    fprintf(fid, ' \t& %s ', colLabels{j});
  end
  fprintf(fid, '\\\\\n');
  fprintf(fid, '\\hline\\\\\n');

  % print each row
  for i = 1:nRows
    fprintf(fid, '%s ', rowLabels{i});
    for j = 1:nCols
      fprintf(fid, '\t & %2.2f (%2.2f) ', mean_(i,j), std_(i,j));
    end
    if i~=nRows
      fprintf(fid, '\\\\\n');
    end
  end
  fprintf(fid, '\n\\end{tabular}');
  fclose(fid);

