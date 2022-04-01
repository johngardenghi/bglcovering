prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(2);
u := 20cm;
path pols[];
pair p[];
R := .5589895213u;
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

x1 = .5000000001u;
y1 = .2499999999u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .5000000000u;
y2 = .7500000000u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.5000000001u,        0.2499999999u)--(        0.0000000000u,        0.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2499999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2499999999u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.5000000001u,        0.2499999999u)--(        1.0000000000u,        0.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2499999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2499999999u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.5000000001u,        0.2499999999u)--(        0.0000307161u,        0.4999999999u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5000000000u,        0.7500000000u)--(        0.0000307161u,        0.4999999999u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2499999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2499999999u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.7500000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.7500000000u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000307161u,        0.4999999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000307161u,        0.4999999999u) withcolor Green;
draw (        0.5000000001u,        0.2499999999u)--(        0.9999692840u,        0.5000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5000000000u,        0.7500000000u)--(        0.9999692840u,        0.5000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2499999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2499999999u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.7500000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.7500000000u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9999692840u,        0.5000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9999692840u,        0.5000000000u) withcolor Green;
draw (        0.5000000000u,        0.7500000000u)--(        1.0000000000u,        1.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.7500000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.7500000000u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        1.0000000000u) withcolor Yellow;
draw (        0.5000000000u,        0.7500000000u)--(        0.0000000000u,        1.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.7500000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.7500000000u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        1.0000000000u) withcolor Yellow;
% r =         0.5589895213
endfig;
end;
