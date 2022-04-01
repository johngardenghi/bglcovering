prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(1);
u := 20cm;
path pols[];
pair p[];
R := .5772801448u;
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

x1 = .5000000001u;
y1 = .2886751345u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.5000000001u,        0.2886751345u)--(        0.0000000000u,        0.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751345u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751345u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.5000000001u,        0.2886751345u)--(        1.0000000000u,        0.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751345u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751345u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.5000000001u,        0.2886751345u)--(        0.5000000000u,        0.8660254038u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751345u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751345u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.8660254038u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.8660254038u) withcolor Yellow;
% r =         0.5772801448
endfig;
end;
