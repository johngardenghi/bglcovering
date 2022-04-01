prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(40);
u := 20cm;
path pols[];
pair p[];
R := .0932071842u;
D := 2*R;

pols[1] = 
(        0.0000000000u,        0.0000000000u)--
(        0.1111111111u,        0.0000000000u)--
(        0.1666666667u,        0.0962250449u)--
(        0.0962250449u,        0.1666666667u)--
(        0.0000000000u,        0.1111111111u)--
(        0.0000000000u,        0.0000000000u)--
cycle;
pols[2] = 
(        0.1666666667u,        0.0962250449u)--
(        0.2222222222u,        0.0000000000u)--
(        0.3333333333u,        0.0000000000u)--
(        0.3888888889u,        0.0962250449u)--
(        0.3333333333u,        0.1924500897u)--
(        0.1666666667u,        0.0962250449u)--
cycle;
pols[3] = 
(        0.1924500897u,        0.3333333333u)--
(        0.0962250449u,        0.3888888889u)--
(        0.0000000000u,        0.3333333333u)--
(        0.0000000000u,        0.2222222222u)--
(        0.0962250449u,        0.1666666667u)--
(        0.1924500897u,        0.3333333333u)--
cycle;
pols[4] = 
(        0.0962250449u,        0.1666666667u)--
(        0.1666666667u,        0.0962250449u)--
(        0.3333333333u,        0.1924500897u)--
(        0.1924500897u,        0.3333333333u)--
(        0.0962250449u,        0.1666666667u)--
cycle;
pols[5] = 
(        0.1924500897u,        0.3333333333u)--
(        0.3333333333u,        0.1924500897u)--
(        0.4444444444u,        0.1924500897u)--
(        0.5000000000u,        0.2886751346u)--
(        0.2886751346u,        0.5000000000u)--
(        0.1924500897u,        0.4444444444u)--
(        0.1924500897u,        0.3333333333u)--
cycle;
pols[6] = 
(        0.6666666667u,        0.1924500897u)--
(        0.6111111111u,        0.0962250449u)--
(        0.6666666667u,        0.0000000000u)--
(        0.7777777778u,        0.0000000000u)--
(        0.8333333333u,        0.0962250449u)--
(        0.6666666667u,        0.1924500897u)--
cycle;
pols[7] = 
(        0.8333333333u,        0.0962250449u)--
(        0.8888888889u,        0.0000000000u)--
(        1.0000000000u,        0.0000000000u)--
(        1.0000000000u,        0.1111111111u)--
(        0.9037749551u,        0.1666666667u)--
(        0.8333333333u,        0.0962250449u)--
cycle;
pols[8] = 
(        0.9037749551u,        0.1666666667u)--
(        1.0000000000u,        0.2222222222u)--
(        1.0000000000u,        0.3333333333u)--
(        0.9037749551u,        0.3888888889u)--
(        0.8075499103u,        0.3333333333u)--
(        0.9037749551u,        0.1666666667u)--
cycle;
pols[9] = 
(        0.9037749551u,        0.1666666667u)--
(        0.8075499103u,        0.3333333333u)--
(        0.6666666667u,        0.1924500897u)--
(        0.8333333333u,        0.0962250449u)--
(        0.9037749551u,        0.1666666667u)--
cycle;
pols[10] = 
(        0.5000000000u,        0.2886751346u)--
(        0.5555555556u,        0.1924500897u)--
(        0.6666666667u,        0.1924500897u)--
(        0.8075499103u,        0.3333333333u)--
(        0.8075499103u,        0.4444444444u)--
(        0.7113248654u,        0.5000000000u)--
(        0.5000000000u,        0.2886751346u)--
cycle;
pols[11] = 
(        0.8075499103u,        0.6666666667u)--
(        0.9037749551u,        0.6111111111u)--
(        1.0000000000u,        0.6666666667u)--
(        1.0000000000u,        0.7777777778u)--
(        0.9037749551u,        0.8333333333u)--
(        0.8075499103u,        0.6666666667u)--
cycle;
pols[12] = 
(        0.9037749551u,        0.8333333333u)--
(        1.0000000000u,        0.8888888889u)--
(        1.0000000000u,        1.0000000000u)--
(        0.8888888889u,        1.0000000000u)--
(        0.8333333333u,        0.9037749551u)--
(        0.9037749551u,        0.8333333333u)--
cycle;
pols[13] = 
(        0.8333333333u,        0.9037749551u)--
(        0.7777777778u,        1.0000000000u)--
(        0.6666666667u,        1.0000000000u)--
(        0.6111111111u,        0.9037749551u)--
(        0.6666666667u,        0.8075499103u)--
(        0.8333333333u,        0.9037749551u)--
cycle;
pols[14] = 
(        0.8333333333u,        0.9037749551u)--
(        0.6666666667u,        0.8075499103u)--
(        0.8075499103u,        0.6666666667u)--
(        0.9037749551u,        0.8333333333u)--
(        0.8333333333u,        0.9037749551u)--
cycle;
pols[15] = 
(        0.8075499103u,        0.6666666667u)--
(        0.6666666667u,        0.8075499103u)--
(        0.5555555556u,        0.8075499103u)--
(        0.5000000000u,        0.7113248654u)--
(        0.7113248654u,        0.5000000000u)--
(        0.8075499103u,        0.5555555556u)--
(        0.8075499103u,        0.6666666667u)--
cycle;
pols[16] = 
(        0.3333333333u,        0.8075499103u)--
(        0.3888888889u,        0.9037749551u)--
(        0.3333333333u,        1.0000000000u)--
(        0.2222222222u,        1.0000000000u)--
(        0.1666666667u,        0.9037749551u)--
(        0.3333333333u,        0.8075499103u)--
cycle;
pols[17] = 
(        0.1666666667u,        0.9037749551u)--
(        0.1111111111u,        1.0000000000u)--
(        0.0000000000u,        1.0000000000u)--
(        0.0000000000u,        0.8888888889u)--
(        0.0962250449u,        0.8333333333u)--
(        0.1666666667u,        0.9037749551u)--
cycle;
pols[18] = 
(        0.0962250449u,        0.8333333333u)--
(        0.0000000000u,        0.7777777778u)--
(        0.0000000000u,        0.6666666667u)--
(        0.0962250449u,        0.6111111111u)--
(        0.1924500897u,        0.6666666667u)--
(        0.0962250449u,        0.8333333333u)--
cycle;
pols[19] = 
(        0.0962250449u,        0.8333333333u)--
(        0.1924500897u,        0.6666666667u)--
(        0.3333333333u,        0.8075499103u)--
(        0.1666666667u,        0.9037749551u)--
(        0.0962250449u,        0.8333333333u)--
cycle;
pols[20] = 
(        0.3333333333u,        0.8075499103u)--
(        0.1924500897u,        0.6666666667u)--
(        0.1924500897u,        0.5555555556u)--
(        0.2886751346u,        0.5000000000u)--
(        0.5000000000u,        0.7113248654u)--
(        0.4444444444u,        0.8075499103u)--
(        0.3333333333u,        0.8075499103u)--
cycle;
pols[21] = 
(        0.5000000000u,        0.7113248654u)--
(        0.2886751346u,        0.5000000000u)--
(        0.5000000000u,        0.2886751346u)--
(        0.7113248654u,        0.5000000000u)--
(        0.5000000000u,        0.7113248654u)--
cycle;
fill pols[1] withcolor 0.9Dandelion;
fill pols[2] withcolor 0.9Dandelion;
fill pols[3] withcolor 0.9Dandelion;
fill pols[4] withcolor 0.9Dandelion;
fill pols[5] withcolor 0.9Dandelion;
fill pols[6] withcolor 0.9Dandelion;
fill pols[7] withcolor 0.9Dandelion;
fill pols[8] withcolor 0.9Dandelion;
fill pols[9] withcolor 0.9Dandelion;
fill pols[10] withcolor 0.9Dandelion;
fill pols[11] withcolor 0.9Dandelion;
fill pols[12] withcolor 0.9Dandelion;
fill pols[13] withcolor 0.9Dandelion;
fill pols[14] withcolor 0.9Dandelion;
fill pols[15] withcolor 0.9Dandelion;
fill pols[16] withcolor 0.9Dandelion;
fill pols[17] withcolor 0.9Dandelion;
fill pols[18] withcolor 0.9Dandelion;
fill pols[19] withcolor 0.9Dandelion;
fill pols[20] withcolor 0.9Dandelion;
fill pols[21] withcolor 0.9Dandelion;
draw (        0.1666666667u,        0.0962250449u)--(        0.0962250449u,        0.1666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3333333333u,        0.1924500897u)--(        0.1666666667u,        0.0962250449u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.0962250449u,        0.1666666667u)--(        0.1924500897u,        0.3333333333u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.0962250449u,        0.1666666667u)--(        0.1666666667u,        0.0962250449u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1666666667u,        0.0962250449u)--(        0.3333333333u,        0.1924500897u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3333333333u,        0.1924500897u)--(        0.1924500897u,        0.3333333333u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1924500897u,        0.3333333333u)--(        0.0962250449u,        0.1666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1924500897u,        0.3333333333u)--(        0.3333333333u,        0.1924500897u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5000000000u,        0.2886751346u)--(        0.2886751346u,        0.5000000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.8333333333u,        0.0962250449u)--(        0.6666666667u,        0.1924500897u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.9037749551u,        0.1666666667u)--(        0.8333333333u,        0.0962250449u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.8075499103u,        0.3333333333u)--(        0.9037749551u,        0.1666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.9037749551u,        0.1666666667u)--(        0.8075499103u,        0.3333333333u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.8075499103u,        0.3333333333u)--(        0.6666666667u,        0.1924500897u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.6666666667u,        0.1924500897u)--(        0.8333333333u,        0.0962250449u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.8333333333u,        0.0962250449u)--(        0.9037749551u,        0.1666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.6666666667u,        0.1924500897u)--(        0.8075499103u,        0.3333333333u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.7113248654u,        0.5000000000u)--(        0.5000000000u,        0.2886751346u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.9037749551u,        0.8333333333u)--(        0.8075499103u,        0.6666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.8333333333u,        0.9037749551u)--(        0.9037749551u,        0.8333333333u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.6666666667u,        0.8075499103u)--(        0.8333333333u,        0.9037749551u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.8333333333u,        0.9037749551u)--(        0.6666666667u,        0.8075499103u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.6666666667u,        0.8075499103u)--(        0.8075499103u,        0.6666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.8075499103u,        0.6666666667u)--(        0.9037749551u,        0.8333333333u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.9037749551u,        0.8333333333u)--(        0.8333333333u,        0.9037749551u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.8075499103u,        0.6666666667u)--(        0.6666666667u,        0.8075499103u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5000000000u,        0.7113248654u)--(        0.7113248654u,        0.5000000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1666666667u,        0.9037749551u)--(        0.3333333333u,        0.8075499103u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.0962250449u,        0.8333333333u)--(        0.1666666667u,        0.9037749551u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1924500897u,        0.6666666667u)--(        0.0962250449u,        0.8333333333u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.0962250449u,        0.8333333333u)--(        0.1924500897u,        0.6666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1924500897u,        0.6666666667u)--(        0.3333333333u,        0.8075499103u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3333333333u,        0.8075499103u)--(        0.1666666667u,        0.9037749551u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1666666667u,        0.9037749551u)--(        0.0962250449u,        0.8333333333u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3333333333u,        0.8075499103u)--(        0.1924500897u,        0.6666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2886751346u,        0.5000000000u)--(        0.5000000000u,        0.7113248654u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5000000000u,        0.7113248654u)--(        0.2886751346u,        0.5000000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2886751346u,        0.5000000000u)--(        0.5000000000u,        0.2886751346u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5000000000u,        0.2886751346u)--(        0.7113248654u,        0.5000000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.7113248654u,        0.5000000000u)--(        0.5000000000u,        0.7113248654u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.0000000000u,        0.0000000000u)--(        0.1111111111u,        0.0000000000u) withpen pencircle scaled 1.0;
draw (        0.1111111111u,        0.0000000000u)--(        0.1666666667u,        0.0962250449u) withpen pencircle scaled 1.0;
draw (        0.0962250449u,        0.1666666667u)--(        0.0000000000u,        0.1111111111u) withpen pencircle scaled 1.0;
draw (        0.0000000000u,        0.1111111111u)--(        0.0000000000u,        0.0000000000u) withpen pencircle scaled 1.0;
draw (        0.1666666667u,        0.0962250449u)--(        0.2222222222u,        0.0000000000u) withpen pencircle scaled 1.0;
draw (        0.2222222222u,        0.0000000000u)--(        0.3333333333u,        0.0000000000u) withpen pencircle scaled 1.0;
draw (        0.3333333333u,        0.0000000000u)--(        0.3888888889u,        0.0962250449u) withpen pencircle scaled 1.0;
draw (        0.3888888889u,        0.0962250449u)--(        0.3333333333u,        0.1924500897u) withpen pencircle scaled 1.0;
draw (        0.1924500897u,        0.3333333333u)--(        0.0962250449u,        0.3888888889u) withpen pencircle scaled 1.0;
draw (        0.0962250449u,        0.3888888889u)--(        0.0000000000u,        0.3333333333u) withpen pencircle scaled 1.0;
draw (        0.0000000000u,        0.3333333333u)--(        0.0000000000u,        0.2222222222u) withpen pencircle scaled 1.0;
draw (        0.0000000000u,        0.2222222222u)--(        0.0962250449u,        0.1666666667u) withpen pencircle scaled 1.0;
draw (        0.3333333333u,        0.1924500897u)--(        0.4444444444u,        0.1924500897u) withpen pencircle scaled 1.0;
draw (        0.4444444444u,        0.1924500897u)--(        0.5000000000u,        0.2886751346u) withpen pencircle scaled 1.0;
draw (        0.2886751346u,        0.5000000000u)--(        0.1924500897u,        0.4444444444u) withpen pencircle scaled 1.0;
draw (        0.1924500897u,        0.4444444444u)--(        0.1924500897u,        0.3333333333u) withpen pencircle scaled 1.0;
draw (        0.6666666667u,        0.1924500897u)--(        0.6111111111u,        0.0962250449u) withpen pencircle scaled 1.0;
draw (        0.6111111111u,        0.0962250449u)--(        0.6666666667u,        0.0000000000u) withpen pencircle scaled 1.0;
draw (        0.6666666667u,        0.0000000000u)--(        0.7777777778u,        0.0000000000u) withpen pencircle scaled 1.0;
draw (        0.7777777778u,        0.0000000000u)--(        0.8333333333u,        0.0962250449u) withpen pencircle scaled 1.0;
draw (        0.8333333333u,        0.0962250449u)--(        0.8888888889u,        0.0000000000u) withpen pencircle scaled 1.0;
draw (        0.8888888889u,        0.0000000000u)--(        1.0000000000u,        0.0000000000u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,        0.0000000000u)--(        1.0000000000u,        0.1111111111u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,        0.1111111111u)--(        0.9037749551u,        0.1666666667u) withpen pencircle scaled 1.0;
draw (        0.9037749551u,        0.1666666667u)--(        1.0000000000u,        0.2222222222u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,        0.2222222222u)--(        1.0000000000u,        0.3333333333u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,        0.3333333333u)--(        0.9037749551u,        0.3888888889u) withpen pencircle scaled 1.0;
draw (        0.9037749551u,        0.3888888889u)--(        0.8075499103u,        0.3333333333u) withpen pencircle scaled 1.0;
draw (        0.5000000000u,        0.2886751346u)--(        0.5555555556u,        0.1924500897u) withpen pencircle scaled 1.0;
draw (        0.5555555556u,        0.1924500897u)--(        0.6666666667u,        0.1924500897u) withpen pencircle scaled 1.0;
draw (        0.8075499103u,        0.3333333333u)--(        0.8075499103u,        0.4444444444u) withpen pencircle scaled 1.0;
draw (        0.8075499103u,        0.4444444444u)--(        0.7113248654u,        0.5000000000u) withpen pencircle scaled 1.0;
draw (        0.8075499103u,        0.6666666667u)--(        0.9037749551u,        0.6111111111u) withpen pencircle scaled 1.0;
draw (        0.9037749551u,        0.6111111111u)--(        1.0000000000u,        0.6666666667u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,        0.6666666667u)--(        1.0000000000u,        0.7777777778u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,        0.7777777778u)--(        0.9037749551u,        0.8333333333u) withpen pencircle scaled 1.0;
draw (        0.9037749551u,        0.8333333333u)--(        1.0000000000u,        0.8888888889u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,        0.8888888889u)--(        1.0000000000u,        1.0000000000u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,        1.0000000000u)--(        0.8888888889u,        1.0000000000u) withpen pencircle scaled 1.0;
draw (        0.8888888889u,        1.0000000000u)--(        0.8333333333u,        0.9037749551u) withpen pencircle scaled 1.0;
draw (        0.8333333333u,        0.9037749551u)--(        0.7777777778u,        1.0000000000u) withpen pencircle scaled 1.0;
draw (        0.7777777778u,        1.0000000000u)--(        0.6666666667u,        1.0000000000u) withpen pencircle scaled 1.0;
draw (        0.6666666667u,        1.0000000000u)--(        0.6111111111u,        0.9037749551u) withpen pencircle scaled 1.0;
draw (        0.6111111111u,        0.9037749551u)--(        0.6666666667u,        0.8075499103u) withpen pencircle scaled 1.0;
draw (        0.6666666667u,        0.8075499103u)--(        0.5555555556u,        0.8075499103u) withpen pencircle scaled 1.0;
draw (        0.5555555556u,        0.8075499103u)--(        0.5000000000u,        0.7113248654u) withpen pencircle scaled 1.0;
draw (        0.7113248654u,        0.5000000000u)--(        0.8075499103u,        0.5555555556u) withpen pencircle scaled 1.0;
draw (        0.8075499103u,        0.5555555556u)--(        0.8075499103u,        0.6666666667u) withpen pencircle scaled 1.0;
draw (        0.3333333333u,        0.8075499103u)--(        0.3888888889u,        0.9037749551u) withpen pencircle scaled 1.0;
draw (        0.3888888889u,        0.9037749551u)--(        0.3333333333u,        1.0000000000u) withpen pencircle scaled 1.0;
draw (        0.3333333333u,        1.0000000000u)--(        0.2222222222u,        1.0000000000u) withpen pencircle scaled 1.0;
draw (        0.2222222222u,        1.0000000000u)--(        0.1666666667u,        0.9037749551u) withpen pencircle scaled 1.0;
draw (        0.1666666667u,        0.9037749551u)--(        0.1111111111u,        1.0000000000u) withpen pencircle scaled 1.0;
draw (        0.1111111111u,        1.0000000000u)--(        0.0000000000u,        1.0000000000u) withpen pencircle scaled 1.0;
draw (        0.0000000000u,        1.0000000000u)--(        0.0000000000u,        0.8888888889u) withpen pencircle scaled 1.0;
draw (        0.0000000000u,        0.8888888889u)--(        0.0962250449u,        0.8333333333u) withpen pencircle scaled 1.0;
draw (        0.0962250449u,        0.8333333333u)--(        0.0000000000u,        0.7777777778u) withpen pencircle scaled 1.0;
draw (        0.0000000000u,        0.7777777778u)--(        0.0000000000u,        0.6666666667u) withpen pencircle scaled 1.0;
draw (        0.0000000000u,        0.6666666667u)--(        0.0962250449u,        0.6111111111u) withpen pencircle scaled 1.0;
draw (        0.0962250449u,        0.6111111111u)--(        0.1924500897u,        0.6666666667u) withpen pencircle scaled 1.0;
draw (        0.1924500897u,        0.6666666667u)--(        0.1924500897u,        0.5555555556u) withpen pencircle scaled 1.0;
draw (        0.1924500897u,        0.5555555556u)--(        0.2886751346u,        0.5000000000u) withpen pencircle scaled 1.0;
draw (        0.5000000000u,        0.7113248654u)--(        0.4444444444u,        0.8075499103u) withpen pencircle scaled 1.0;
draw (        0.4444444444u,        0.8075499103u)--(        0.3333333333u,        0.8075499103u) withpen pencircle scaled 1.0;

x1 = .9353869397u;
y1 = .0613534634u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .7260182385u;
y2 = .0648210783u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .8499362567u;
y3 = .1593266819u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .9263386178u;
y4 = .2793382552u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .6784814142u;
y5 = .1396351547u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .7475827363u;
y6 = .2768450100u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .8426437848u;
y7 = .3975257297u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .5862608249u;
y8 = .2804590300u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .6823575569u;
y9 = .4174147200u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = .7577147553u;
y10 = .5462819576u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;

x11 = .8660592980u;
y11 = .6706959407u;
draw fullcircle scaled D shifted z11 withcolor RoyalBlue withpen pencircle scaled 1.0;

x12 = .9373971091u;
y12 = .7318952594u;
draw fullcircle scaled D shifted z12 withcolor RoyalBlue withpen pencircle scaled 1.0;

x13 = .2724947234u;
y13 = .0644554250u;
draw fullcircle scaled D shifted z13 withcolor RoyalBlue withpen pencircle scaled 1.0;

x14 = .3321166523u;
y14 = .1424405810u;
draw fullcircle scaled D shifted z14 withcolor RoyalBlue withpen pencircle scaled 1.0;

x15 = .4318511865u;
y15 = .2832467915u;
draw fullcircle scaled D shifted z15 withcolor RoyalBlue withpen pencircle scaled 1.0;

x16 = .5314075775u;
y16 = .4247824444u;
draw fullcircle scaled D shifted z16 withcolor RoyalBlue withpen pencircle scaled 1.0;

x17 = .6198527251u;
y17 = .5682496133u;
draw fullcircle scaled D shifted z17 withcolor RoyalBlue withpen pencircle scaled 1.0;

x18 = .7199429222u;
y18 = .7089804764u;
draw fullcircle scaled D shifted z18 withcolor RoyalBlue withpen pencircle scaled 1.0;

x19 = .8279366923u;
y19 = .8334785168u;
draw fullcircle scaled D shifted z19 withcolor RoyalBlue withpen pencircle scaled 1.0;

x20 = .9362979430u;
y20 = .9328542237u;
draw fullcircle scaled D shifted z20 withcolor RoyalBlue withpen pencircle scaled 1.0;

x21 = .0637020570u;
y21 = .0671457762u;
draw fullcircle scaled D shifted z21 withcolor RoyalBlue withpen pencircle scaled 1.0;

x22 = .1720633076u;
y22 = .1665214833u;
draw fullcircle scaled D shifted z22 withcolor RoyalBlue withpen pencircle scaled 1.0;

x23 = .2800570778u;
y23 = .2910195235u;
draw fullcircle scaled D shifted z23 withcolor RoyalBlue withpen pencircle scaled 1.0;

x24 = .3801472748u;
y24 = .4317503868u;
draw fullcircle scaled D shifted z24 withcolor RoyalBlue withpen pencircle scaled 1.0;

x25 = .4685924225u;
y25 = .5752175557u;
draw fullcircle scaled D shifted z25 withcolor RoyalBlue withpen pencircle scaled 1.0;

x26 = .5681488135u;
y26 = .7167532083u;
draw fullcircle scaled D shifted z26 withcolor RoyalBlue withpen pencircle scaled 1.0;

x27 = .6678833476u;
y27 = .8575594190u;
draw fullcircle scaled D shifted z27 withcolor RoyalBlue withpen pencircle scaled 1.0;

x28 = .7275052765u;
y28 = .9355445750u;
draw fullcircle scaled D shifted z28 withcolor RoyalBlue withpen pencircle scaled 1.0;

x29 = .0626028909u;
y29 = .2681047407u;
draw fullcircle scaled D shifted z29 withcolor RoyalBlue withpen pencircle scaled 1.0;

x30 = .1339407018u;
y30 = .3293040593u;
draw fullcircle scaled D shifted z30 withcolor RoyalBlue withpen pencircle scaled 1.0;

x31 = .2422852446u;
y31 = .4537180423u;
draw fullcircle scaled D shifted z31 withcolor RoyalBlue withpen pencircle scaled 1.0;

x32 = .3176424431u;
y32 = .5825852800u;
draw fullcircle scaled D shifted z32 withcolor RoyalBlue withpen pencircle scaled 1.0;

x33 = .4137391752u;
y33 = .7195409700u;
draw fullcircle scaled D shifted z33 withcolor RoyalBlue withpen pencircle scaled 1.0;

x34 = .1573562152u;
y34 = .6024742704u;
draw fullcircle scaled D shifted z34 withcolor RoyalBlue withpen pencircle scaled 1.0;

x35 = .2524172638u;
y35 = .7231549901u;
draw fullcircle scaled D shifted z35 withcolor RoyalBlue withpen pencircle scaled 1.0;

x36 = .3215185858u;
y36 = .8603648454u;
draw fullcircle scaled D shifted z36 withcolor RoyalBlue withpen pencircle scaled 1.0;

x37 = .0736613823u;
y37 = .7206617447u;
draw fullcircle scaled D shifted z37 withcolor RoyalBlue withpen pencircle scaled 1.0;

x38 = .1500637432u;
y38 = .8406733181u;
draw fullcircle scaled D shifted z38 withcolor RoyalBlue withpen pencircle scaled 1.0;

x39 = .2739817615u;
y39 = .9351789218u;
draw fullcircle scaled D shifted z39 withcolor RoyalBlue withpen pencircle scaled 1.0;

x40 = .0646130603u;
y40 = .9386465366u;
draw fullcircle scaled D shifted z40 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.7260182385u,        0.0648210783u)--(        0.8189925327u,        0.0714058929u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8499362567u,        0.1593266819u)--(        0.8189925327u,        0.0714058929u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7260182385u,        0.0648210783u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7260182385u,        0.0648210783u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8499362567u,        0.1593266819u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8499362567u,        0.1593266819u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8189925327u,        0.0714058929u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8189925327u,        0.0714058929u) withcolor Green;
draw (        0.7260182385u,        0.0648210783u)--(        0.6330426939u,        0.0582539414u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6784814142u,        0.1396351547u)--(        0.6330426939u,        0.0582539414u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7260182385u,        0.0648210783u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7260182385u,        0.0648210783u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6784814142u,        0.1396351547u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6784814142u,        0.1396351547u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6330426939u,        0.0582539414u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6330426939u,        0.0582539414u) withcolor Green;
draw (        0.8499362567u,        0.1593266819u)--(        0.8373304089u,        0.2516774901u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9263386178u,        0.2793382552u)--(        0.8373304089u,        0.2516774901u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8499362567u,        0.1593266819u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8499362567u,        0.1593266819u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9263386178u,        0.2793382552u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9263386178u,        0.2793382552u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8373304089u,        0.2516774901u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8373304089u,        0.2516774901u) withcolor MidnightBlue;
draw (        0.8499362567u,        0.1593266819u)--(        0.7601886267u,        0.1844850221u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6784814142u,        0.1396351547u)--(        0.7601886267u,        0.1844850221u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8499362567u,        0.1593266819u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8499362567u,        0.1593266819u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6784814142u,        0.1396351547u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6784814142u,        0.1396351547u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7601886267u,        0.1844850221u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7601886267u,        0.1844850221u) withcolor MidnightBlue;
draw (        0.8499362567u,        0.1593266819u)--(        0.7601913046u,        0.1844945732u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7475827363u,        0.2768450100u)--(        0.7601913046u,        0.1844945732u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8499362567u,        0.1593266819u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8499362567u,        0.1593266819u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7475827363u,        0.2768450100u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7475827363u,        0.2768450100u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7601913046u,        0.1844945732u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7601913046u,        0.1844945732u) withcolor MidnightBlue;
draw (        0.8499362567u,        0.1593266819u)--(        0.8373276884u,        0.2516771187u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7475827363u,        0.2768450100u)--(        0.8373276884u,        0.2516771187u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8499362567u,        0.1593266819u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8499362567u,        0.1593266819u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7475827363u,        0.2768450100u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7475827363u,        0.2768450100u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8373276884u,        0.2516771187u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8373276884u,        0.2516771187u) withcolor MidnightBlue;
draw (        0.9263386178u,        0.2793382552u)--(        1.0000000000u,        0.2222222222u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9263386178u,        0.2793382552u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9263386178u,        0.2793382552u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        0.2222222222u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        0.2222222222u) withcolor Yellow;
draw (        0.9263386178u,        0.2793382552u)--(        0.8373290336u,        0.2516819160u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7475827363u,        0.2768450100u)--(        0.8373290336u,        0.2516819160u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9263386178u,        0.2793382552u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9263386178u,        0.2793382552u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7475827363u,        0.2768450100u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7475827363u,        0.2768450100u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8373290336u,        0.2516819160u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8373290336u,        0.2516819160u) withcolor MidnightBlue;
draw (        0.9263386178u,        0.2793382552u)--(        0.8365923205u,        0.3045013492u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7475827363u,        0.2768450100u)--(        0.8365923205u,        0.3045013492u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9263386178u,        0.2793382552u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9263386178u,        0.2793382552u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7475827363u,        0.2768450100u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7475827363u,        0.2768450100u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8365923205u,        0.3045013492u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8365923205u,        0.3045013492u) withcolor MidnightBlue;
draw (        0.9263386178u,        0.2793382552u)--(        0.8365961350u,        0.3045149501u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8426437848u,        0.3975257297u)--(        0.8365961350u,        0.3045149501u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9263386178u,        0.2793382552u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9263386178u,        0.2793382552u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8426437848u,        0.3975257297u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8426437848u,        0.3975257297u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8365961350u,        0.3045149501u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8365961350u,        0.3045149501u) withcolor MidnightBlue;
draw (        0.6784814142u,        0.1396351547u)--(        0.6658802149u,        0.2319865974u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7475827363u,        0.2768450100u)--(        0.6658802149u,        0.2319865974u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6784814142u,        0.1396351547u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6784814142u,        0.1396351547u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7475827363u,        0.2768450100u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7475827363u,        0.2768450100u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6658802149u,        0.2319865974u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6658802149u,        0.2319865974u) withcolor MidnightBlue;
draw (        0.6784814142u,        0.1396351547u)--(        0.7601839355u,        0.1844935674u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7475827363u,        0.2768450100u)--(        0.7601839355u,        0.1844935674u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6784814142u,        0.1396351547u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6784814142u,        0.1396351547u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7475827363u,        0.2768450100u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7475827363u,        0.2768450100u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7601839355u,        0.1844935674u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7601839355u,        0.1844935674u) withcolor MidnightBlue;
draw (        0.6784814142u,        0.1396351547u)--(        0.6658716887u,        0.2319854336u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5862608249u,        0.2804590300u)--(        0.6658716887u,        0.2319854336u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6784814142u,        0.1396351547u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6784814142u,        0.1396351547u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5862608249u,        0.2804590300u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5862608249u,        0.2804590300u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6658716887u,        0.2319854336u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6658716887u,        0.2319854336u) withcolor MidnightBlue;
draw (        0.7475827363u,        0.2768450100u)--(        0.8365879274u,        0.3045154842u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8426437848u,        0.3975257297u)--(        0.8365879274u,        0.3045154842u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7475827363u,        0.2768450100u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7475827363u,        0.2768450100u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8426437848u,        0.3975257297u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8426437848u,        0.3975257297u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8365879274u,        0.3045154842u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8365879274u,        0.3045154842u) withcolor MidnightBlue;
draw (        0.7475827363u,        0.2768450100u)--(        0.6658765080u,        0.2319933496u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5862608249u,        0.2804590300u)--(        0.6658765080u,        0.2319933496u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7475827363u,        0.2768450100u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7475827363u,        0.2768450100u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5862608249u,        0.2804590300u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5862608249u,        0.2804590300u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6658765080u,        0.2319933496u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6658765080u,        0.2319933496u) withcolor MidnightBlue;
draw (        0.7475827363u,        0.2768450100u)--(        0.6679670531u,        0.3253106905u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5862608249u,        0.2804590300u)--(        0.6679670531u,        0.3253106905u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7475827363u,        0.2768450100u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7475827363u,        0.2768450100u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5862608249u,        0.2804590300u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5862608249u,        0.2804590300u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6679670531u,        0.3253106905u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6679670531u,        0.3253106905u) withcolor MidnightBlue;
draw (        0.7475827363u,        0.2768450100u)--(        0.6679750812u,        0.3253238759u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6823575569u,        0.4174147200u)--(        0.6679750812u,        0.3253238759u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7475827363u,        0.2768450100u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7475827363u,        0.2768450100u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6823575569u,        0.4174147200u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6823575569u,        0.4174147200u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6679750812u,        0.3253238759u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6679750812u,        0.3253238759u) withcolor MidnightBlue;
draw (        0.8426437848u,        0.3975257297u)--(        0.7682313293u,        0.4536538479u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6823575569u,        0.4174147200u)--(        0.7682313293u,        0.4536538479u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8426437848u,        0.3975257297u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8426437848u,        0.3975257297u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6823575569u,        0.4174147200u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6823575569u,        0.4174147200u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7682313293u,        0.4536538479u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7682313293u,        0.4536538479u) withcolor MidnightBlue;
draw (        0.8426437848u,        0.3975257297u)--(        0.7682446309u,        0.4536714782u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7577147553u,        0.5462819576u)--(        0.7682446309u,        0.4536714782u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8426437848u,        0.3975257297u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8426437848u,        0.3975257297u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7577147553u,        0.5462819576u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7577147553u,        0.5462819576u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7682446309u,        0.4536714782u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7682446309u,        0.4536714782u) withcolor MidnightBlue;
draw (        0.5862608249u,        0.2804590300u)--(        0.5555555556u,        0.1924500897u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5862608249u,        0.2804590300u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5862608249u,        0.2804590300u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5555555556u,        0.1924500897u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5555555556u,        0.1924500897u) withcolor Yellow;
draw (        0.5862608249u,        0.2804590300u)--(        0.6679583805u,        0.3253264857u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6823575569u,        0.4174147200u)--(        0.6679583805u,        0.3253264857u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5862608249u,        0.2804590300u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5862608249u,        0.2804590300u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6823575569u,        0.4174147200u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6823575569u,        0.4174147200u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6679583805u,        0.3253264857u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6679583805u,        0.3253264857u) withcolor MidnightBlue;
draw (        0.5862608249u,        0.2804590300u)--(        0.5099983245u,        0.3340464372u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4318511865u,        0.2832467915u)--(        0.5099983245u,        0.3340464372u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5862608249u,        0.2804590300u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5862608249u,        0.2804590300u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4318511865u,        0.2832467915u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4318511865u,        0.2832467915u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5099983245u,        0.3340464372u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5099983245u,        0.3340464372u) withcolor MidnightBlue;
draw (        0.5862608249u,        0.2804590300u)--(        0.5100109498u,        0.3340644004u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5314075775u,        0.4247824444u)--(        0.5100109498u,        0.3340644004u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5862608249u,        0.2804590300u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5862608249u,        0.2804590300u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5314075775u,        0.4247824444u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5314075775u,        0.4247824444u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5100109498u,        0.3340644004u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5100109498u,        0.3340644004u) withcolor MidnightBlue;
draw (        0.6823575569u,        0.4174147200u)--(        0.7682248365u,        0.4536692297u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7577147553u,        0.5462819576u)--(        0.7682248365u,        0.4536692297u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6823575569u,        0.4174147200u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6823575569u,        0.4174147200u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7577147553u,        0.5462819576u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7577147553u,        0.5462819576u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7682248365u,        0.4536692297u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7682248365u,        0.4536692297u) withcolor MidnightBlue;
draw (        0.6823575569u,        0.4174147200u)--(        0.6095427531u,        0.4756004875u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5314075775u,        0.4247824444u)--(        0.6095427531u,        0.4756004875u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6823575569u,        0.4174147200u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6823575569u,        0.4174147200u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5314075775u,        0.4247824444u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5314075775u,        0.4247824444u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6095427531u,        0.4756004875u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6095427531u,        0.4756004875u) withcolor MidnightBlue;
draw (        0.6823575569u,        0.4174147200u)--(        0.6095529610u,        0.4756132595u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6198527251u,        0.5682496133u)--(        0.6095529610u,        0.4756132595u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6823575569u,        0.4174147200u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6823575569u,        0.4174147200u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6198527251u,        0.5682496133u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6198527251u,        0.5682496133u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6095529610u,        0.4756132595u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6095529610u,        0.4756132595u) withcolor MidnightBlue;
draw (        0.7577147553u,        0.5462819576u)--(        0.7791606268u,        0.6369883729u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8660592980u,        0.6706959407u)--(        0.7791606268u,        0.6369883729u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7577147553u,        0.5462819576u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7577147553u,        0.5462819576u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8660592980u,        0.6706959407u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8660592980u,        0.6706959407u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7791606268u,        0.6369883729u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7791606268u,        0.6369883729u) withcolor MidnightBlue;
draw (        0.7577147553u,        0.5462819576u)--(        0.6985037609u,        0.6182655483u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6198527251u,        0.5682496133u)--(        0.6985037609u,        0.6182655483u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7577147553u,        0.5462819576u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7577147553u,        0.5462819576u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6198527251u,        0.5682496133u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6198527251u,        0.5682496133u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6985037609u,        0.6182655483u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6985037609u,        0.6182655483u) withcolor MidnightBlue;
draw (        0.7577147553u,        0.5462819576u)--(        0.6985103028u,        0.6182709289u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7199429222u,        0.7089804764u)--(        0.6985103028u,        0.6182709289u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7577147553u,        0.5462819576u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7577147553u,        0.5462819576u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7199429222u,        0.7089804764u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7199429222u,        0.7089804764u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6985103028u,        0.6182709289u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6985103028u,        0.6182709289u) withcolor MidnightBlue;
draw (        0.7577147553u,        0.5462819576u)--(        0.7791473746u,        0.6369915051u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7199429222u,        0.7089804764u)--(        0.7791473746u,        0.6369915051u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7577147553u,        0.5462819576u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7577147553u,        0.5462819576u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7199429222u,        0.7089804764u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7199429222u,        0.7089804764u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7791473746u,        0.6369915051u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7791473746u,        0.6369915051u) withcolor MidnightBlue;
draw (        0.8660592980u,        0.6706959407u)--(        0.7791565009u,        0.6369990116u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7199429222u,        0.7089804764u)--(        0.7791565009u,        0.6369990116u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8660592980u,        0.6706959407u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8660592980u,        0.6706959407u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7199429222u,        0.7089804764u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7199429222u,        0.7089804764u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7791565009u,        0.6369990116u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7791565009u,        0.6369990116u) withcolor MidnightBlue;
draw (        0.8660592980u,        0.6706959407u)--(        0.8068457193u,        0.7426774055u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7199429222u,        0.7089804764u)--(        0.8068457193u,        0.7426774055u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8660592980u,        0.6706959407u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8660592980u,        0.6706959407u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7199429222u,        0.7089804764u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7199429222u,        0.7089804764u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8068457193u,        0.7426774055u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8068457193u,        0.7426774055u) withcolor MidnightBlue;
draw (        0.8660592980u,        0.6706959407u)--(        0.8068566407u,        0.7426863882u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8279366923u,        0.8334785168u)--(        0.8068566407u,        0.7426863882u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8660592980u,        0.6706959407u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8660592980u,        0.6706959407u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8279366923u,        0.8334785168u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8279366923u,        0.8334785168u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8068566407u,        0.7426863882u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8068566407u,        0.7426863882u) withcolor MidnightBlue;
draw (        0.3321166523u,        0.1424405810u)--(        0.3531993654u,        0.2332320916u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4318511865u,        0.2832467915u)--(        0.3531993654u,        0.2332320916u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3321166523u,        0.1424405810u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3321166523u,        0.1424405810u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4318511865u,        0.2832467915u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4318511865u,        0.2832467915u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3531993654u,        0.2332320916u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3531993654u,        0.2332320916u) withcolor MidnightBlue;
draw (        0.3321166523u,        0.1424405810u)--(        0.4107684735u,        0.1924552810u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4318511865u,        0.2832467915u)--(        0.4107684735u,        0.1924552810u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3321166523u,        0.1424405810u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3321166523u,        0.1424405810u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4318511865u,        0.2832467915u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4318511865u,        0.2832467915u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4107684735u,        0.1924552810u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4107684735u,        0.1924552810u) withcolor Green;
draw (        0.3321166523u,        0.1424405810u)--(        0.2589699334u,        0.2002085366u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1720633076u,        0.1665214833u)--(        0.2589699334u,        0.2002085366u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3321166523u,        0.1424405810u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3321166523u,        0.1424405810u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1720633076u,        0.1665214833u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1720633076u,        0.1665214833u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2589699334u,        0.2002085366u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2589699334u,        0.2002085366u) withcolor MidnightBlue;
draw (        0.3321166523u,        0.1424405810u)--(        0.2589836212u,        0.2002258640u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2800570778u,        0.2910195235u)--(        0.2589836212u,        0.2002258640u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3321166523u,        0.1424405810u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3321166523u,        0.1424405810u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2800570778u,        0.2910195235u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2800570778u,        0.2910195235u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2589836212u,        0.2002258640u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2589836212u,        0.2002258640u) withcolor MidnightBlue;
draw (        0.3321166523u,        0.1424405810u)--(        0.3531901090u,        0.2332342405u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2800570778u,        0.2910195235u)--(        0.3531901090u,        0.2332342405u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3321166523u,        0.1424405810u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3321166523u,        0.1424405810u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2800570778u,        0.2910195235u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2800570778u,        0.2910195235u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3531901090u,        0.2332342405u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3531901090u,        0.2332342405u) withcolor MidnightBlue;
draw (        0.4318511865u,        0.2832467915u)--(        0.4532765275u,        0.3739580584u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5314075775u,        0.4247824444u)--(        0.4532765275u,        0.3739580584u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4318511865u,        0.2832467915u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4318511865u,        0.2832467915u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5314075775u,        0.4247824444u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5314075775u,        0.4247824444u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4532765275u,        0.3739580584u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4532765275u,        0.3739580584u) withcolor MidnightBlue;
draw (        0.4318511865u,        0.2832467915u)--(        0.5099822365u,        0.3340711775u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5314075775u,        0.4247824444u)--(        0.5099822365u,        0.3340711775u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4318511865u,        0.2832467915u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4318511865u,        0.2832467915u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5314075775u,        0.4247824444u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5314075775u,        0.4247824444u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5099822365u,        0.3340711775u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5099822365u,        0.3340711775u) withcolor MidnightBlue;
draw (        0.4318511865u,        0.2832467915u)--(        0.3531944808u,        0.2332397738u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2800570778u,        0.2910195235u)--(        0.3531944808u,        0.2332397738u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4318511865u,        0.2832467915u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4318511865u,        0.2832467915u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2800570778u,        0.2910195235u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2800570778u,        0.2910195235u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3531944808u,        0.2332397738u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3531944808u,        0.2332397738u) withcolor MidnightBlue;
draw (        0.4318511865u,        0.2832467915u)--(        0.3587137836u,        0.3410265412u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2800570778u,        0.2910195235u)--(        0.3587137836u,        0.3410265412u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4318511865u,        0.2832467915u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4318511865u,        0.2832467915u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2800570778u,        0.2910195235u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2800570778u,        0.2910195235u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3587137836u,        0.3410265412u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3587137836u,        0.3410265412u) withcolor MidnightBlue;
draw (        0.4318511865u,        0.2832467915u)--(        0.3587234414u,        0.3410387639u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3801472748u,        0.4317503868u)--(        0.3587234414u,        0.3410387639u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4318511865u,        0.2832467915u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4318511865u,        0.2832467915u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3801472748u,        0.4317503868u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3801472748u,        0.4317503868u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3587234414u,        0.3410387639u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3587234414u,        0.3410387639u) withcolor MidnightBlue;
draw (        0.4318511865u,        0.2832467915u)--(        0.4532750200u,        0.3739584145u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3801472748u,        0.4317503868u)--(        0.4532750200u,        0.3739584145u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4318511865u,        0.2832467915u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4318511865u,        0.2832467915u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3801472748u,        0.4317503868u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3801472748u,        0.4317503868u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4532750200u,        0.3739584145u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4532750200u,        0.3739584145u) withcolor MidnightBlue;
draw (        0.5314075775u,        0.4247824444u)--(        0.5417273032u,        0.5174165766u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6198527251u,        0.5682496133u)--(        0.5417273032u,        0.5174165766u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5314075775u,        0.4247824444u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5314075775u,        0.4247824444u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6198527251u,        0.5682496133u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6198527251u,        0.5682496133u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5417273032u,        0.5174165766u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5417273032u,        0.5174165766u) withcolor MidnightBlue;
draw (        0.5314075775u,        0.4247824444u)--(        0.6095329995u,        0.4756154811u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6198527251u,        0.5682496133u)--(        0.6095329995u,        0.4756154811u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5314075775u,        0.4247824444u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5314075775u,        0.4247824444u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6198527251u,        0.5682496133u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6198527251u,        0.5682496133u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6095329995u,        0.4756154811u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6095329995u,        0.4756154811u) withcolor MidnightBlue;
draw (        0.5314075775u,        0.4247824444u)--(        0.4532757198u,        0.3739593001u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3801472748u,        0.4317503868u)--(        0.4532757198u,        0.3739593001u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5314075775u,        0.4247824444u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5314075775u,        0.4247824444u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3801472748u,        0.4317503868u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3801472748u,        0.4317503868u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4532757198u,        0.3739593001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4532757198u,        0.3739593001u) withcolor MidnightBlue;
draw (        0.5314075775u,        0.4247824444u)--(        0.4582791325u,        0.4825735311u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3801472748u,        0.4317503868u)--(        0.4582791325u,        0.4825735311u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5314075775u,        0.4247824444u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5314075775u,        0.4247824444u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3801472748u,        0.4317503868u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3801472748u,        0.4317503868u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4582791325u,        0.4825735311u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4582791325u,        0.4825735311u) withcolor MidnightBlue;
draw (        0.5314075775u,        0.4247824444u)--(        0.4582857983u,        0.4825819650u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4685924225u,        0.5752175557u)--(        0.4582857983u,        0.4825819650u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5314075775u,        0.4247824444u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5314075775u,        0.4247824444u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4685924225u,        0.5752175557u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4685924225u,        0.5752175557u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4582857983u,        0.4825819650u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4582857983u,        0.4825819650u) withcolor MidnightBlue;
draw (        0.5314075775u,        0.4247824444u)--(        0.5417142018u,        0.5174180352u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4685924225u,        0.5752175557u)--(        0.5417142018u,        0.5174180352u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5314075775u,        0.4247824444u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5314075775u,        0.4247824444u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4685924225u,        0.5752175557u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4685924225u,        0.5752175557u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5417142018u,        0.5174180352u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5417142018u,        0.5174180352u) withcolor MidnightBlue;
draw (        0.6198527251u,        0.5682496133u)--(        0.6412970711u,        0.6589563892u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7199429222u,        0.7089804764u)--(        0.6412970711u,        0.6589563892u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6198527251u,        0.5682496133u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6198527251u,        0.5682496133u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7199429222u,        0.7089804764u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7199429222u,        0.7089804764u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6412970711u,        0.6589563892u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6412970711u,        0.6589563892u) withcolor MidnightBlue;
draw (        0.6198527251u,        0.5682496133u)--(        0.6984985762u,        0.6182737004u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7199429222u,        0.7089804764u)--(        0.6984985762u,        0.6182737004u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6198527251u,        0.5682496133u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6198527251u,        0.5682496133u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7199429222u,        0.7089804764u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7199429222u,        0.7089804764u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6984985762u,        0.6182737004u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6984985762u,        0.6182737004u) withcolor MidnightBlue;
draw (        0.6198527251u,        0.5682496133u)--(        0.5417208675u,        0.5174264689u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4685924225u,        0.5752175557u)--(        0.5417208675u,        0.5174264689u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6198527251u,        0.5682496133u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6198527251u,        0.5682496133u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4685924225u,        0.5752175557u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4685924225u,        0.5752175557u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5417208675u,        0.5174264689u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5417208675u,        0.5174264689u) withcolor MidnightBlue;
draw (        0.6198527251u,        0.5682496133u)--(        0.5467242802u,        0.6260407001u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4685924225u,        0.5752175557u)--(        0.5467242802u,        0.6260407001u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6198527251u,        0.5682496133u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6198527251u,        0.5682496133u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4685924225u,        0.5752175557u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4685924225u,        0.5752175557u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5467242802u,        0.6260407001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5467242802u,        0.6260407001u) withcolor MidnightBlue;
draw (        0.6198527251u,        0.5682496133u)--(        0.5467249798u,        0.6260415854u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5681488135u,        0.7167532083u)--(        0.5467249798u,        0.6260415854u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6198527251u,        0.5682496133u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6198527251u,        0.5682496133u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5681488135u,        0.7167532083u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5681488135u,        0.7167532083u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5467249798u,        0.6260415854u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5467249798u,        0.6260415854u) withcolor MidnightBlue;
draw (        0.6198527251u,        0.5682496133u)--(        0.6412765588u,        0.6589612362u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5681488135u,        0.7167532083u)--(        0.6412765588u,        0.6589612362u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6198527251u,        0.5682496133u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6198527251u,        0.5682496133u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5681488135u,        0.7167532083u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5681488135u,        0.7167532083u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6412765588u,        0.6589612362u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6412765588u,        0.6589612362u) withcolor MidnightBlue;
draw (        0.7199429222u,        0.7089804764u)--(        0.7410388052u,        0.7997689278u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8279366923u,        0.8334785168u)--(        0.7410388052u,        0.7997689278u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7199429222u,        0.7089804764u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7199429222u,        0.7089804764u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8279366923u,        0.8334785168u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8279366923u,        0.8334785168u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7410388052u,        0.7997689278u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7410388052u,        0.7997689278u) withcolor MidnightBlue;
draw (        0.7199429222u,        0.7089804764u)--(        0.8068408093u,        0.7426900654u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8279366923u,        0.8334785168u)--(        0.8068408093u,        0.7426900654u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7199429222u,        0.7089804764u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7199429222u,        0.7089804764u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8279366923u,        0.8334785168u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8279366923u,        0.8334785168u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8068408093u,        0.7426900654u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8068408093u,        0.7426900654u) withcolor MidnightBlue;
draw (        0.7199429222u,        0.7089804764u)--(        0.6412862164u,        0.6589734586u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5681488135u,        0.7167532083u)--(        0.6412862164u,        0.6589734586u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7199429222u,        0.7089804764u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7199429222u,        0.7089804764u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5681488135u,        0.7167532083u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5681488135u,        0.7167532083u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6412862164u,        0.6589734586u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6412862164u,        0.6589734586u) withcolor MidnightBlue;
draw (        0.7199429222u,        0.7089804764u)--(        0.6468055192u,        0.7667602261u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5681488135u,        0.7167532083u)--(        0.6468055192u,        0.7667602261u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7199429222u,        0.7089804764u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7199429222u,        0.7089804764u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5681488135u,        0.7167532083u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5681488135u,        0.7167532083u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6468055192u,        0.7667602261u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6468055192u,        0.7667602261u) withcolor MidnightBlue;
draw (        0.7199429222u,        0.7089804764u)--(        0.6468098910u,        0.7667657594u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6678833476u,        0.8575594190u)--(        0.6468098910u,        0.7667657594u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7199429222u,        0.7089804764u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7199429222u,        0.7089804764u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6678833476u,        0.8575594190u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6678833476u,        0.8575594190u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6468098910u,        0.7667657594u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6468098910u,        0.7667657594u) withcolor MidnightBlue;
draw (        0.7199429222u,        0.7089804764u)--(        0.7410163787u,        0.7997741359u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6678833476u,        0.8575594190u)--(        0.7410163787u,        0.7997741359u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7199429222u,        0.7089804764u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7199429222u,        0.7089804764u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6678833476u,        0.8575594190u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6678833476u,        0.8575594190u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7410163787u,        0.7997741359u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7410163787u,        0.7997741359u) withcolor MidnightBlue;
draw (        0.8279366923u,        0.8334785168u)--(        0.7410300665u,        0.7997914633u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6678833476u,        0.8575594190u)--(        0.7410300665u,        0.7997914633u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8279366923u,        0.8334785168u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8279366923u,        0.8334785168u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6678833476u,        0.8575594190u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6678833476u,        0.8575594190u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7410300665u,        0.7997914633u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7410300665u,        0.7997914633u) withcolor MidnightBlue;
draw (        0.1720633076u,        0.1665214833u)--(        0.1931591906u,        0.2573099346u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2800570778u,        0.2910195235u)--(        0.1931591906u,        0.2573099346u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1720633076u,        0.1665214833u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1720633076u,        0.1665214833u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2800570778u,        0.2910195235u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2800570778u,        0.2910195235u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1931591906u,        0.2573099346u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1931591906u,        0.2573099346u) withcolor MidnightBlue;
draw (        0.1720633076u,        0.1665214833u)--(        0.2589611948u,        0.2002310721u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2800570778u,        0.2910195235u)--(        0.2589611948u,        0.2002310721u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1720633076u,        0.1665214833u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1720633076u,        0.1665214833u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2800570778u,        0.2910195235u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2800570778u,        0.2910195235u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2589611948u,        0.2002310721u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2589611948u,        0.2002310721u) withcolor MidnightBlue;
draw (        0.1720633076u,        0.1665214833u)--(        0.1931433593u,        0.2573136118u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1339407018u,        0.3293040593u)--(        0.1931433593u,        0.2573136118u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1720633076u,        0.1665214833u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1720633076u,        0.1665214833u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1339407018u,        0.3293040593u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1339407018u,        0.3293040593u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1931433593u,        0.2573136118u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1931433593u,        0.2573136118u) withcolor MidnightBlue;
draw (        0.2800570778u,        0.2910195235u)--(        0.3015014239u,        0.3817262994u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3801472748u,        0.4317503868u)--(        0.3015014239u,        0.3817262994u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2800570778u,        0.2910195235u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2800570778u,        0.2910195235u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3801472748u,        0.4317503868u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3801472748u,        0.4317503868u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3015014239u,        0.3817262994u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3015014239u,        0.3817262994u) withcolor MidnightBlue;
draw (        0.2800570778u,        0.2910195235u)--(        0.3587029287u,        0.3410436109u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3801472748u,        0.4317503868u)--(        0.3587029287u,        0.3410436109u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2800570778u,        0.2910195235u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2800570778u,        0.2910195235u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3801472748u,        0.4317503868u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3801472748u,        0.4317503868u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3587029287u,        0.3410436109u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3587029287u,        0.3410436109u) withcolor MidnightBlue;
draw (        0.2800570778u,        0.2910195235u)--(        0.1931542806u,        0.2573225945u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1339407018u,        0.3293040593u)--(        0.1931542806u,        0.2573225945u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2800570778u,        0.2910195235u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2800570778u,        0.2910195235u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1339407018u,        0.3293040593u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1339407018u,        0.3293040593u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1931542806u,        0.2573225945u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1931542806u,        0.2573225945u) withcolor MidnightBlue;
draw (        0.2800570778u,        0.2910195235u)--(        0.2208434990u,        0.3630009882u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1339407018u,        0.3293040593u)--(        0.2208434990u,        0.3630009882u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2800570778u,        0.2910195235u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2800570778u,        0.2910195235u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1339407018u,        0.3293040593u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1339407018u,        0.3293040593u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2208434990u,        0.3630009882u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2208434990u,        0.3630009882u) withcolor MidnightBlue;
draw (        0.2800570778u,        0.2910195235u)--(        0.2208526254u,        0.3630084948u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2422852446u,        0.4537180423u)--(        0.2208526254u,        0.3630084948u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2800570778u,        0.2910195235u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2800570778u,        0.2910195235u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2422852446u,        0.4537180423u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2422852446u,        0.4537180423u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2208526254u,        0.3630084948u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2208526254u,        0.3630084948u) withcolor MidnightBlue;
draw (        0.2800570778u,        0.2910195235u)--(        0.3014896971u,        0.3817290710u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2422852446u,        0.4537180423u)--(        0.3014896971u,        0.3817290710u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2800570778u,        0.2910195235u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2800570778u,        0.2910195235u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2422852446u,        0.4537180423u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2422852446u,        0.4537180423u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3014896971u,        0.3817290710u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3014896971u,        0.3817290710u) withcolor MidnightBlue;
draw (        0.3801472748u,        0.4317503868u)--(        0.3904670006u,        0.5243845190u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4685924225u,        0.5752175557u)--(        0.3904670006u,        0.5243845190u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3801472748u,        0.4317503868u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3801472748u,        0.4317503868u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4685924225u,        0.5752175557u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4685924225u,        0.5752175557u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3904670006u,        0.5243845190u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3904670006u,        0.5243845190u) withcolor MidnightBlue;
draw (        0.3801472748u,        0.4317503868u)--(        0.4582726967u,        0.4825834236u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4685924225u,        0.5752175557u)--(        0.4582726967u,        0.4825834236u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3801472748u,        0.4317503868u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3801472748u,        0.4317503868u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4685924225u,        0.5752175557u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4685924225u,        0.5752175557u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4582726967u,        0.4825834236u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4582726967u,        0.4825834236u) withcolor MidnightBlue;
draw (        0.3801472748u,        0.4317503868u)--(        0.3014962391u,        0.3817344517u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2422852446u,        0.4537180423u)--(        0.3014962391u,        0.3817344517u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3801472748u,        0.4317503868u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3801472748u,        0.4317503868u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2422852446u,        0.4537180423u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2422852446u,        0.4537180423u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3014962391u,        0.3817344517u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3014962391u,        0.3817344517u) withcolor MidnightBlue;
draw (        0.3801472748u,        0.4317503868u)--(        0.3904470391u,        0.5243867406u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3176424431u,        0.5825852800u)--(        0.3904470391u,        0.5243867406u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3801472748u,        0.4317503868u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3801472748u,        0.4317503868u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3176424431u,        0.5825852800u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3176424431u,        0.5825852800u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3904470391u,        0.5243867406u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3904470391u,        0.5243867406u) withcolor MidnightBlue;
draw (        0.4685924225u,        0.5752175557u)--(        0.4900177633u,        0.6659288227u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5681488135u,        0.7167532083u)--(        0.4900177633u,        0.6659288227u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4685924225u,        0.5752175557u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4685924225u,        0.5752175557u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5681488135u,        0.7167532083u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5681488135u,        0.7167532083u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4900177633u,        0.6659288227u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4900177633u,        0.6659288227u) withcolor MidnightBlue;
draw (        0.4685924225u,        0.5752175557u)--(        0.5467234727u,        0.6260419414u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5681488135u,        0.7167532083u)--(        0.5467234727u,        0.6260419414u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4685924225u,        0.5752175557u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4685924225u,        0.5752175557u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5681488135u,        0.7167532083u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5681488135u,        0.7167532083u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5467234727u,        0.6260419414u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5467234727u,        0.6260419414u) withcolor MidnightBlue;
draw (        0.4685924225u,        0.5752175557u)--(        0.3904572469u,        0.5243995125u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3176424431u,        0.5825852800u)--(        0.3904572469u,        0.5243995125u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4685924225u,        0.5752175557u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4685924225u,        0.5752175557u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3176424431u,        0.5825852800u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3176424431u,        0.5825852800u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3904572469u,        0.5243995125u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3904572469u,        0.5243995125u) withcolor MidnightBlue;
draw (        0.4685924225u,        0.5752175557u)--(        0.4899890503u,        0.6659355997u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4137391752u,        0.7195409700u)--(        0.4899890503u,        0.6659355997u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4685924225u,        0.5752175557u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4685924225u,        0.5752175557u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4137391752u,        0.7195409700u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4137391752u,        0.7195409700u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4899890503u,        0.6659355997u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4899890503u,        0.6659355997u) withcolor MidnightBlue;
draw (        0.5681488135u,        0.7167532083u)--(        0.5892315265u,        0.8075447189u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6678833476u,        0.8575594190u)--(        0.5892315265u,        0.8075447189u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5681488135u,        0.7167532083u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5681488135u,        0.7167532083u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6678833476u,        0.8575594190u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6678833476u,        0.8575594190u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5892315265u,        0.8075447189u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5892315265u,        0.8075447189u) withcolor Green;
draw (        0.5681488135u,        0.7167532083u)--(        0.6468006346u,        0.7667679084u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6678833476u,        0.8575594190u)--(        0.6468006346u,        0.7667679084u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5681488135u,        0.7167532083u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5681488135u,        0.7167532083u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6678833476u,        0.8575594190u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6678833476u,        0.8575594190u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6468006346u,        0.7667679084u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6468006346u,        0.7667679084u) withcolor MidnightBlue;
draw (        0.5681488135u,        0.7167532083u)--(        0.4900016755u,        0.6659535627u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4137391752u,        0.7195409700u)--(        0.4900016755u,        0.6659535627u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5681488135u,        0.7167532083u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5681488135u,        0.7167532083u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4137391752u,        0.7195409700u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4137391752u,        0.7195409700u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4900016755u,        0.6659535627u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4900016755u,        0.6659535627u) withcolor MidnightBlue;
draw (        0.1339407018u,        0.3293040593u)--(        0.2208393730u,        0.3630116271u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2422852446u,        0.4537180423u)--(        0.2208393730u,        0.3630116271u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1339407018u,        0.3293040593u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1339407018u,        0.3293040593u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2422852446u,        0.4537180423u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2422852446u,        0.4537180423u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2208393730u,        0.3630116271u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2208393730u,        0.3630116271u) withcolor MidnightBlue;
draw (        0.2422852446u,        0.4537180423u)--(        0.2317751635u,        0.5463307703u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3176424431u,        0.5825852800u)--(        0.2317751635u,        0.5463307703u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2422852446u,        0.4537180423u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2422852446u,        0.4537180423u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3176424431u,        0.5825852800u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3176424431u,        0.5825852800u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2317751635u,        0.5463307703u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2317751635u,        0.5463307703u) withcolor MidnightBlue;
draw (        0.2422852446u,        0.4537180423u)--(        0.2317553691u,        0.5463285218u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1573562152u,        0.6024742704u)--(        0.2317553691u,        0.5463285218u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2422852446u,        0.4537180423u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2422852446u,        0.4537180423u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1573562152u,        0.6024742704u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1573562152u,        0.6024742704u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2317553691u,        0.5463285218u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2317553691u,        0.5463285218u) withcolor MidnightBlue;
draw (        0.3176424431u,        0.5825852800u)--(        0.3320416195u,        0.6746735143u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4137391752u,        0.7195409700u)--(        0.3320416195u,        0.6746735143u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3176424431u,        0.5825852800u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3176424431u,        0.5825852800u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4137391752u,        0.7195409700u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4137391752u,        0.7195409700u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3320416195u,        0.6746735143u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3320416195u,        0.6746735143u) withcolor MidnightBlue;
draw (        0.3176424431u,        0.5825852800u)--(        0.2317686707u,        0.5463461522u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1573562152u,        0.6024742704u)--(        0.2317686707u,        0.5463461522u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3176424431u,        0.5825852800u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3176424431u,        0.5825852800u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1573562152u,        0.6024742704u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1573562152u,        0.6024742704u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2317686707u,        0.5463461522u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2317686707u,        0.5463461522u) withcolor MidnightBlue;
draw (        0.3176424431u,        0.5825852800u)--(        0.3320249188u,        0.6746761241u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2524172638u,        0.7231549901u)--(        0.3320249188u,        0.6746761241u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3176424431u,        0.5825852800u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3176424431u,        0.5825852800u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2524172638u,        0.7231549901u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2524172638u,        0.7231549901u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3320249188u,        0.6746761241u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3320249188u,        0.6746761241u) withcolor MidnightBlue;
draw (        0.4137391752u,        0.7195409700u)--(        0.4444444444u,        0.8075499103u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4137391752u,        0.7195409700u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4137391752u,        0.7195409700u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4444444444u,        0.8075499103u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4444444444u,        0.8075499103u) withcolor Yellow;
draw (        0.4137391752u,        0.7195409700u)--(        0.3320329469u,        0.6746893095u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2524172638u,        0.7231549901u)--(        0.3320329469u,        0.6746893095u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4137391752u,        0.7195409700u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4137391752u,        0.7195409700u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2524172638u,        0.7231549901u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2524172638u,        0.7231549901u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3320329469u,        0.6746893095u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3320329469u,        0.6746893095u) withcolor MidnightBlue;
draw (        0.4137391752u,        0.7195409700u)--(        0.3341234921u,        0.7680066505u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2524172638u,        0.7231549901u)--(        0.3341234921u,        0.7680066505u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4137391752u,        0.7195409700u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4137391752u,        0.7195409700u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2524172638u,        0.7231549901u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2524172638u,        0.7231549901u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3341234921u,        0.7680066505u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3341234921u,        0.7680066505u) withcolor MidnightBlue;
draw (        0.4137391752u,        0.7195409700u)--(        0.3341283115u,        0.7680145666u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3215185858u,        0.8603648454u)--(        0.3341283115u,        0.7680145666u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4137391752u,        0.7195409700u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4137391752u,        0.7195409700u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3215185858u,        0.8603648454u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3215185858u,        0.8603648454u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3341283115u,        0.7680145666u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3341283115u,        0.7680145666u) withcolor MidnightBlue;
draw (        0.1573562152u,        0.6024742704u)--(        0.1634120727u,        0.6954845159u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2524172638u,        0.7231549901u)--(        0.1634120727u,        0.6954845159u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1573562152u,        0.6024742704u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1573562152u,        0.6024742704u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2524172638u,        0.7231549901u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2524172638u,        0.7231549901u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1634120727u,        0.6954845159u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1634120727u,        0.6954845159u) withcolor MidnightBlue;
draw (        0.1573562152u,        0.6024742704u)--(        0.1634038651u,        0.6954850499u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0736613823u,        0.7206617447u)--(        0.1634038651u,        0.6954850499u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1573562152u,        0.6024742704u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1573562152u,        0.6024742704u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0736613823u,        0.7206617447u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0736613823u,        0.7206617447u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1634038651u,        0.6954850499u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1634038651u,        0.6954850499u) withcolor MidnightBlue;
draw (        0.2524172638u,        0.7231549901u)--(        0.2398160645u,        0.8155064327u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3215185858u,        0.8603648454u)--(        0.2398160645u,        0.8155064327u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2524172638u,        0.7231549901u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2524172638u,        0.7231549901u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3215185858u,        0.8603648454u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3215185858u,        0.8603648454u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2398160645u,        0.8155064327u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2398160645u,        0.8155064327u) withcolor MidnightBlue;
draw (        0.2524172638u,        0.7231549901u)--(        0.3341197851u,        0.7680134028u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3215185858u,        0.8603648454u)--(        0.3341197851u,        0.7680134028u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2524172638u,        0.7231549901u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2524172638u,        0.7231549901u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3215185858u,        0.8603648454u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3215185858u,        0.8603648454u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3341197851u,        0.7680134028u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3341197851u,        0.7680134028u) withcolor MidnightBlue;
draw (        0.2524172638u,        0.7231549901u)--(        0.1634076796u,        0.6954986508u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0736613823u,        0.7206617447u)--(        0.1634076796u,        0.6954986508u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2524172638u,        0.7231549901u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2524172638u,        0.7231549901u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0736613823u,        0.7206617447u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0736613823u,        0.7206617447u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1634076796u,        0.6954986508u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1634076796u,        0.6954986508u) withcolor MidnightBlue;
draw (        0.2524172638u,        0.7231549901u)--(        0.1626709665u,        0.7483180840u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0736613823u,        0.7206617447u)--(        0.1626709665u,        0.7483180840u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2524172638u,        0.7231549901u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2524172638u,        0.7231549901u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0736613823u,        0.7206617447u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0736613823u,        0.7206617447u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1626709665u,        0.7483180840u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1626709665u,        0.7483180840u) withcolor MidnightBlue;
draw (        0.2524172638u,        0.7231549901u)--(        0.1626723117u,        0.7483228813u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1500637432u,        0.8406733181u)--(        0.1626723117u,        0.7483228813u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2524172638u,        0.7231549901u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2524172638u,        0.7231549901u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1500637432u,        0.8406733181u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1500637432u,        0.8406733181u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1626723117u,        0.7483228813u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1626723117u,        0.7483228813u) withcolor MidnightBlue;
draw (        0.2524172638u,        0.7231549901u)--(        0.2398086953u,        0.8155054269u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1500637432u,        0.8406733181u)--(        0.2398086953u,        0.8155054269u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2524172638u,        0.7231549901u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2524172638u,        0.7231549901u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1500637432u,        0.8406733181u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1500637432u,        0.8406733181u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2398086953u,        0.8155054269u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2398086953u,        0.8155054269u) withcolor MidnightBlue;
draw (        0.3215185858u,        0.8603648454u)--(        0.2398113733u,        0.8155149781u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1500637432u,        0.8406733181u)--(        0.2398113733u,        0.8155149781u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3215185858u,        0.8603648454u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3215185858u,        0.8603648454u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1500637432u,        0.8406733181u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1500637432u,        0.8406733181u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2398113733u,        0.8155149781u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2398113733u,        0.8155149781u) withcolor MidnightBlue;
draw (        0.3215185858u,        0.8603648454u)--(        0.3669573061u,        0.9417460587u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2739817615u,        0.9351789218u)--(        0.3669573061u,        0.9417460587u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3215185858u,        0.8603648454u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3215185858u,        0.8603648454u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2739817615u,        0.9351789218u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2739817615u,        0.9351789218u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3669573061u,        0.9417460587u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3669573061u,        0.9417460587u) withcolor Green;
draw (        0.0736613823u,        0.7206617447u)--(        0.0000000000u,        0.7777777778u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0736613823u,        0.7206617447u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0736613823u,        0.7206617447u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.7777777778u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.7777777778u) withcolor Yellow;
draw (        0.0736613823u,        0.7206617447u)--(        0.1626695912u,        0.7483225099u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1500637432u,        0.8406733181u)--(        0.1626695912u,        0.7483225099u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0736613823u,        0.7206617447u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0736613823u,        0.7206617447u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1500637432u,        0.8406733181u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1500637432u,        0.8406733181u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1626695912u,        0.7483225099u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1626695912u,        0.7483225099u) withcolor MidnightBlue;
draw (        0.1500637432u,        0.8406733181u)--(        0.1810074672u,        0.9285941070u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2739817615u,        0.9351789218u)--(        0.1810074672u,        0.9285941070u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1500637432u,        0.8406733181u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1500637432u,        0.8406733181u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2739817615u,        0.9351789218u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2739817615u,        0.9351789218u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1810074672u,        0.9285941070u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1810074672u,        0.9285941070u) withcolor Green;
% r =         0.0932071842
endfig;
end;
