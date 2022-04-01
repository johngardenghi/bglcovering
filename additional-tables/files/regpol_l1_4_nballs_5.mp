prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(5);
u := 20cm;
path pols[];
pair p[];
R := .3261408863u;
D := 2*R;

pols[1] = 
(        0.0000000000u,        0.0000000000u)--
(        1.0000000000u,        0.0000000000u)--
(        1.0000000000u,        1.0000000000u)--
(        0.0000000000u,        1.0000000000u)--
(        0.0000000000u,        0.0000000000u)--
cycle;
fill pols[1] withcolor 0.9Dandelion;
draw (        0.0000000000u,        0.0000000000u)--(        1.0000000000u,        0.0000000000u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,        0.0000000000u)--(        1.0000000000u,        1.0000000000u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,        1.0000000000u)--(        0.0000000000u,        1.0000000000u) withpen pencircle scaled 1.0;
draw (        0.0000000000u,        1.0000000000u)--(        0.0000000000u,        0.0000000000u) withpen pencircle scaled 1.0;

x1 = .7905173615u;
y1 = .2499933606u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .2905000396u;
y2 = .1482484619u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .7905173613u;
y3 = .7500066393u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .2549096879u;
y4 = .5000000000u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .2905000395u;
y5 = .8517515382u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.7905173615u,        0.2499933606u)--(        1.0000000000u,        0.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7905173615u,        0.2499933606u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7905173615u,        0.2499933606u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.7905173615u,        0.2499933606u)--(        0.5810188742u,        0.0000368413u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2905000396u,        0.1482484619u)--(        0.5810188742u,        0.0000368413u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7905173615u,        0.2499933606u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7905173615u,        0.2499933606u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2905000396u,        0.1482484619u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2905000396u,        0.1482484619u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5810188742u,        0.0000368413u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5810188742u,        0.0000368413u) withcolor Green;
draw (        0.7905173615u,        0.2499933606u)--(        0.5810786877u,        0.4999999999u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.7905173613u,        0.7500066393u)--(        0.5810786877u,        0.4999999999u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7905173615u,        0.2499933606u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7905173615u,        0.2499933606u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7905173613u,        0.7500066393u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7905173613u,        0.7500066393u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5810786877u,        0.4999999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5810786877u,        0.4999999999u) withcolor MidnightBlue;
draw (        0.7905173615u,        0.2499933606u)--(        0.9999560351u,        0.5000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.7905173613u,        0.7500066393u)--(        0.9999560351u,        0.5000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7905173615u,        0.2499933606u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7905173615u,        0.2499933606u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7905173613u,        0.7500066393u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7905173613u,        0.7500066393u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9999560351u,        0.5000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9999560351u,        0.5000000000u) withcolor Green;
draw (        0.7905173615u,        0.2499933606u)--(        0.5810505734u,        0.4999764450u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2549096879u,        0.5000000000u)--(        0.5810505734u,        0.4999764450u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7905173615u,        0.2499933606u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7905173615u,        0.2499933606u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2549096879u,        0.5000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2549096879u,        0.5000000000u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5810505734u,        0.4999764450u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5810505734u,        0.4999764450u) withcolor MidnightBlue;
draw (        0.2905000396u,        0.1482484619u)--(        0.0000000000u,        0.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2905000396u,        0.1482484619u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2905000396u,        0.1482484619u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.2905000396u,        0.1482484619u)--(        0.0000188009u,        0.2965337535u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2549096879u,        0.5000000000u)--(        0.0000188009u,        0.2965337535u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2905000396u,        0.1482484619u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2905000396u,        0.1482484619u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2549096879u,        0.5000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2549096879u,        0.5000000000u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000188009u,        0.2965337535u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000188009u,        0.2965337535u) withcolor Green;
draw (        0.7905173613u,        0.7500066393u)--(        1.0000000000u,        1.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7905173613u,        0.7500066393u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7905173613u,        0.7500066393u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        1.0000000000u) withcolor Yellow;
draw (        0.7905173613u,        0.7500066393u)--(        0.5810505734u,        0.5000235548u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2549096879u,        0.5000000000u)--(        0.5810505734u,        0.5000235548u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7905173613u,        0.7500066393u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7905173613u,        0.7500066393u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2549096879u,        0.5000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2549096879u,        0.5000000000u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5810505734u,        0.5000235548u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5810505734u,        0.5000235548u) withcolor MidnightBlue;
draw (        0.7905173613u,        0.7500066393u)--(        0.5810188741u,        0.9999631587u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2905000395u,        0.8517515382u)--(        0.5810188741u,        0.9999631587u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7905173613u,        0.7500066393u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7905173613u,        0.7500066393u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2905000395u,        0.8517515382u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2905000395u,        0.8517515382u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5810188741u,        0.9999631587u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5810188741u,        0.9999631587u) withcolor Green;
draw (        0.2549096879u,        0.5000000000u)--(        0.0000188009u,        0.7034662464u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2905000395u,        0.8517515382u)--(        0.0000188009u,        0.7034662464u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2549096879u,        0.5000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2549096879u,        0.5000000000u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2905000395u,        0.8517515382u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2905000395u,        0.8517515382u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000188009u,        0.7034662464u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000188009u,        0.7034662464u) withcolor Green;
draw (        0.2905000395u,        0.8517515382u)--(        0.0000000000u,        1.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2905000395u,        0.8517515382u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2905000395u,        0.8517515382u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        1.0000000000u) withcolor Yellow;
% r =         0.3261408863
endfig;
end;
