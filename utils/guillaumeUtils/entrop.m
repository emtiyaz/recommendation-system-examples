function h=entrop(x)
I = (x==0 | x==1);
h = -x.*log(x) - (1-x).*log(1-x);
h(I) = 0;
