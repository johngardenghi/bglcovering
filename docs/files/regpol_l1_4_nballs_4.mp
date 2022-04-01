prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(4);
u := 20cm;
path pols[];
pair p[];
R := .3535332969u;
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

x1 = .7500000001u;
y1 = .7499999999u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .7500000000u;
y2 = .2500000000u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .2499999999u;
y3 = .7500000000u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .2500000000u;
y4 = .2500000000u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.7500000001u,        0.7499999999u)--(        1.0000000000u,        1.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7500000001u,        0.7499999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7500000001u,        0.7499999999u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        1.0000000000u) withcolor Yellow;
draw (        0.7500000001u,        0.7499999999u)--(        0.9999715825u,        0.4999999999u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.7500000000u,        0.2500000000u)--(        0.9999715825u,        0.4999999999u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7500000001u,        0.7499999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7500000001u,        0.7499999999u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7500000000u,        0.2500000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7500000000u,        0.2500000000u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9999715825u,        0.4999999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9999715825u,        0.4999999999u) withcolor Green;
draw (        0.7500000001u,        0.7499999999u)--(        0.5000284176u,        0.5000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.7500000000u,        0.2500000000u)--(        0.5000284176u,        0.5000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7500000001u,        0.7499999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7500000001u,        0.7499999999u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7500000000u,        0.2500000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7500000000u,        0.2500000000u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000284176u,        0.5000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000284176u,        0.5000000000u) withcolor MidnightBlue;
draw (        0.7500000001u,        0.7499999999u)--(        0.5000000000u,        0.5000284176u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2499999999u,        0.7500000000u)--(        0.5000000000u,        0.5000284176u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7500000001u,        0.7499999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7500000001u,        0.7499999999u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2499999999u,        0.7500000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2499999999u,        0.7500000000u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.5000284176u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.5000284176u) withcolor MidnightBlue;
draw (        0.7500000001u,        0.7499999999u)--(        0.5000000000u,        0.9999715822u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2499999999u,        0.7500000000u)--(        0.5000000000u,        0.9999715822u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7500000001u,        0.7499999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7500000001u,        0.7499999999u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2499999999u,        0.7500000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2499999999u,        0.7500000000u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.9999715822u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.9999715822u) withcolor Green;
draw (        0.7500000000u,        0.2500000000u)--(        1.0000000000u,        0.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7500000000u,        0.2500000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7500000000u,        0.2500000000u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.7500000000u,        0.2500000000u)--(        0.5000000000u,        0.0000284176u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2500000000u,        0.2500000000u)--(        0.5000000000u,        0.0000284176u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7500000000u,        0.2500000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7500000000u,        0.2500000000u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2500000000u,        0.2500000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500000000u,        0.2500000000u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.0000284176u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.0000284176u) withcolor Green;
draw (        0.7500000000u,        0.2500000000u)--(        0.5000000000u,        0.4999715824u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2500000000u,        0.2500000000u)--(        0.5000000000u,        0.4999715824u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7500000000u,        0.2500000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7500000000u,        0.2500000000u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2500000000u,        0.2500000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500000000u,        0.2500000000u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4999715824u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4999715824u) withcolor MidnightBlue;
draw (        0.2499999999u,        0.7500000000u)--(        0.0000000000u,        1.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2499999999u,        0.7500000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2499999999u,        0.7500000000u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        1.0000000000u) withcolor Yellow;
draw (        0.2499999999u,        0.7500000000u)--(        0.4999715824u,        0.5000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2500000000u,        0.2500000000u)--(        0.4999715824u,        0.5000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2499999999u,        0.7500000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2499999999u,        0.7500000000u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2500000000u,        0.2500000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500000000u,        0.2500000000u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4999715824u,        0.5000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999715824u,        0.5000000000u) withcolor MidnightBlue;
draw (        0.2499999999u,        0.7500000000u)--(        0.0000284176u,        0.4999999999u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2500000000u,        0.2500000000u)--(        0.0000284176u,        0.4999999999u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2499999999u,        0.7500000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2499999999u,        0.7500000000u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2500000000u,        0.2500000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500000000u,        0.2500000000u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000284176u,        0.4999999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000284176u,        0.4999999999u) withcolor Green;
draw (        0.2500000000u,        0.2500000000u)--(        0.0000000000u,        0.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2500000000u,        0.2500000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500000000u,        0.2500000000u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor Yellow;
% r =         0.3535332969
endfig;
end;
