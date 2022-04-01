prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(3);
u := 20cm;
path pols[];
pair p[];
R := .5038736962u;
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

x1 = .7500093585u;
y1 = .5625147072u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .5000000000u;
y2 = .0625093597u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .2499906416u;
y3 = .5625147073u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.7500093585u,        0.5625147072u)--(        1.0000000000u,        1.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7500093585u,        0.5625147072u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7500093585u,        0.5625147072u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        1.0000000000u) withcolor Yellow;
draw (        0.7500093585u,        0.5625147072u)--(        0.9999812814u,        0.1250187192u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5000000000u,        0.0625093597u)--(        0.9999812814u,        0.1250187192u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7500093585u,        0.5625147072u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7500093585u,        0.5625147072u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.0625093597u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.0625093597u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9999812814u,        0.1250187192u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9999812814u,        0.1250187192u) withcolor Green;
draw (        0.7500093585u,        0.5625147072u)--(        0.5000000001u,        0.9999893036u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2499906416u,        0.5625147073u)--(        0.5000000001u,        0.9999893036u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7500093585u,        0.5625147072u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7500093585u,        0.5625147072u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2499906416u,        0.5625147073u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2499906416u,        0.5625147073u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.9999893036u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.9999893036u) withcolor Green;
draw (        0.5000000000u,        0.0625093597u)--(        0.0000000000u,        0.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.0625093597u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.0625093597u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.5000000000u,        0.0625093597u)--(        1.0000000000u,        0.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.0625093597u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.0625093597u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.5000000000u,        0.0625093597u)--(        0.0000187186u,        0.1250187193u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2499906416u,        0.5625147073u)--(        0.0000187186u,        0.1250187193u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.0625093597u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.0625093597u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2499906416u,        0.5625147073u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2499906416u,        0.5625147073u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000187186u,        0.1250187193u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000187186u,        0.1250187193u) withcolor Green;
draw (        0.2499906416u,        0.5625147073u)--(        0.0000000000u,        1.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2499906416u,        0.5625147073u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2499906416u,        0.5625147073u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        1.0000000000u) withcolor Yellow;
% r =         0.5038736962
endfig;
end;
