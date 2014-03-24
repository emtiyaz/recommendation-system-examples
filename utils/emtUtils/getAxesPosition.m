function rect = getAxesPosition(options)
%options [nRows, nCols, xBottom, yBottom, xTop, yTop, xGap, yGap] =  myProcessOptions(options, 'nRows',2,'nCols',2,'xBottom',0.1,'yBottom',0.1, 'xTop', 0.9, 'yTop',0.9,'xGap',0.05,'yGap',0.05); 

[nRows, nCols, xBottom, yBottom, xTop, yTop, xGap, yGap] =  myProcessOptions(options, 'nRows',2,'nCols',2,'xBottom',0.1,'yBottom',0.1, 'xTop', 0.9, 'yTop',0.9,'xGap',0.05,'yGap',0.05); 

W = xTop - xBottom;
H = yTop - yBottom; 
height = (H - (nRows -1)*yGap)/nRows;
width = (W - (nCols -1)*xGap)/nCols;

count = 1;
for i = nRows:-1:1
for j = 1:nCols
  x = xBottom + (j-1)*(width + xGap);
  y = yBottom + (i-1)*(height+ yGap);
  rect(count,:) = [x y width height];
  count = count +1;
end
end

