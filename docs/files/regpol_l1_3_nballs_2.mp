prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(2);
u := 20cm;
path pols[];
pair p[];
R := .4999298764u;
D := 2*R;

pols[1] = 
(        0.0000000000u,        0.0000000000u)--
(        1.0000000000u,        0.0000000000u)--
(        0.5000000000u,        0.8660254038u)--
(        0.0000000000u,        0.0000000000u)--
cycle;
fill pols[1] withcolor 0.9Dandelion;
draw (        0.0000000000u,        0.0000000000u)--(        1.0000000000u,        0.0000000000u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,        0.0000000000u)--(        0.5000000000u,        0.8660254038u) withpen pencircle scaled 1.0;
draw (        0.5000000000u,        0.8660254038u)--(        0.0000000000u,        0.0000000000u) withpen pencircle scaled 1.0;

x1 = .5000350642u;
y1 = .0000607234u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .2500701197u;
y2 = .4330127058u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.5000350642u,        0.0000607234u)--(        1.0000000000u,        0.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000350642u,        0.0000607234u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000350642u,        0.0000607234u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.5000350642u,        0.0000607234u)--(        0.0001051878u,        0.0000607291u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2500701197u,        0.4330127058u)--(        0.0001051878u,        0.0000607291u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000350642u,        0.0000607234u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000350642u,        0.0000607234u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2500701197u,        0.4330127058u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500701197u,        0.4330127058u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0001051878u,        0.0000607291u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0001051878u,        0.0000607291u) withcolor Green;
draw (        0.5000350642u,        0.0000607234u)--(        0.7499999961u,        0.4330127001u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2500701197u,        0.4330127058u)--(        0.7499999961u,        0.4330127001u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000350642u,        0.0000607234u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000350642u,        0.0000607234u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2500701197u,        0.4330127058u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500701197u,        0.4330127058u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7499999961u,        0.4330127001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7499999961u,        0.4330127001u) withcolor Green;
draw (        0.2500701197u,        0.4330127058u)--(        0.5000000000u,        0.8660254038u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2500701197u,        0.4330127058u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500701197u,        0.4330127058u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.8660254038u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.8660254038u) withcolor Yellow;
% r =         0.4999298764
endfig;
end;
