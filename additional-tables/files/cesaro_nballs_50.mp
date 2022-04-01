prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(50);
u := 20cm;
path pols[];
pair p[];
R := .0829721218u;
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

x1 = .9444423212u;
y1 = .0616302625u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .8392771146u;
y2 = .0988355720u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .9274937486u;
y3 = .2197257381u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .9404835244u;
y4 = .3179047032u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .7142287775u;
y5 = .0499777522u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .7741497143u;
y6 = .2014099061u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .8202322263u;
y7 = .3191860203u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .9289555557u;
y8 = .6804681843u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .9422783475u;
y9 = .7829731243u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = .9444333216u;
y10 = .9383624146u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;

x11 = .6470408194u;
y11 = .1418429022u;
draw fullcircle scaled D shifted z11 withcolor RoyalBlue withpen pencircle scaled 1.0;

x12 = .6727106157u;
y12 = .2924489315u;
draw fullcircle scaled D shifted z12 withcolor RoyalBlue withpen pencircle scaled 1.0;

x13 = .7274447211u;
y13 = .4228157968u;
draw fullcircle scaled D shifted z13 withcolor RoyalBlue withpen pencircle scaled 1.0;

x14 = .7433895439u;
y14 = .5678408587u;
draw fullcircle scaled D shifted z14 withcolor RoyalBlue withpen pencircle scaled 1.0;

x15 = .8073944162u;
y15 = .7034284938u;
draw fullcircle scaled D shifted z15 withcolor RoyalBlue withpen pencircle scaled 1.0;

x16 = .8328491722u;
y16 = .8449558938u;
draw fullcircle scaled D shifted z16 withcolor RoyalBlue withpen pencircle scaled 1.0;

x17 = .7897360084u;
y17 = .9611210712u;
draw fullcircle scaled D shifted z17 withcolor RoyalBlue withpen pencircle scaled 1.0;

x18 = .5621819437u;
y18 = .2731929515u;
draw fullcircle scaled D shifted z18 withcolor RoyalBlue withpen pencircle scaled 1.0;

x19 = .5887818314u;
y19 = .4169959813u;
draw fullcircle scaled D shifted z19 withcolor RoyalBlue withpen pencircle scaled 1.0;

x20 = .6236325111u;
y20 = .5400368830u;
draw fullcircle scaled D shifted z20 withcolor RoyalBlue withpen pencircle scaled 1.0;

x21 = .6626011231u;
y21 = .6703710872u;
draw fullcircle scaled D shifted z21 withcolor RoyalBlue withpen pencircle scaled 1.0;

x22 = .6987225336u;
y22 = .8010836705u;
draw fullcircle scaled D shifted z22 withcolor RoyalBlue withpen pencircle scaled 1.0;

x23 = .6913423248u;
y23 = .9209027120u;
draw fullcircle scaled D shifted z23 withcolor RoyalBlue withpen pencircle scaled 1.0;

x24 = .3206286984u;
y24 = .0729004634u;
draw fullcircle scaled D shifted z24 withcolor RoyalBlue withpen pencircle scaled 1.0;

x25 = .4160894773u;
y25 = .2621257505u;
draw fullcircle scaled D shifted z25 withcolor RoyalBlue withpen pencircle scaled 1.0;

x26 = .4609295037u;
y26 = .3855477990u;
draw fullcircle scaled D shifted z26 withcolor RoyalBlue withpen pencircle scaled 1.0;

x27 = .4793352746u;
y27 = .5262465187u;
draw fullcircle scaled D shifted z27 withcolor RoyalBlue withpen pencircle scaled 1.0;

x28 = .5192180895u;
y28 = .6460893848u;
draw fullcircle scaled D shifted z28 withcolor RoyalBlue withpen pencircle scaled 1.0;

x29 = .5518470317u;
y29 = .7549271244u;
draw fullcircle scaled D shifted z29 withcolor RoyalBlue withpen pencircle scaled 1.0;

x30 = .2188782361u;
y30 = .0550328119u;
draw fullcircle scaled D shifted z30 withcolor RoyalBlue withpen pencircle scaled 1.0;

x31 = .2882406904u;
y31 = .1870191904u;
draw fullcircle scaled D shifted z31 withcolor RoyalBlue withpen pencircle scaled 1.0;

x32 = .3403610580u;
y32 = .3393666653u;
draw fullcircle scaled D shifted z32 withcolor RoyalBlue withpen pencircle scaled 1.0;

x33 = .3521242518u;
y33 = .4885502458u;
draw fullcircle scaled D shifted z33 withcolor RoyalBlue withpen pencircle scaled 1.0;

x34 = .3704458338u;
y34 = .6321023209u;
draw fullcircle scaled D shifted z34 withcolor RoyalBlue withpen pencircle scaled 1.0;

x35 = .3974132415u;
y35 = .7391899244u;
draw fullcircle scaled D shifted z35 withcolor RoyalBlue withpen pencircle scaled 1.0;

x36 = .0555660359u;
y36 = .0616416844u;
draw fullcircle scaled D shifted z36 withcolor RoyalBlue withpen pencircle scaled 1.0;

x37 = .1402378357u;
y37 = .1556180595u;
draw fullcircle scaled D shifted z37 withcolor RoyalBlue withpen pencircle scaled 1.0;

x38 = .2142319893u;
y38 = .2873409642u;
draw fullcircle scaled D shifted z38 withcolor RoyalBlue withpen pencircle scaled 1.0;

x39 = .2375129981u;
y39 = .4370800858u;
draw fullcircle scaled D shifted z39 withcolor RoyalBlue withpen pencircle scaled 1.0;

x40 = .2414982237u;
y40 = .5891898766u;
draw fullcircle scaled D shifted z40 withcolor RoyalBlue withpen pencircle scaled 1.0;

x41 = .2443466175u;
y41 = .7244836743u;
draw fullcircle scaled D shifted z41 withcolor RoyalBlue withpen pencircle scaled 1.0;

x42 = .2884127802u;
y42 = .8404480221u;
draw fullcircle scaled D shifted z42 withcolor RoyalBlue withpen pencircle scaled 1.0;

x43 = .3180818603u;
y43 = .9331627017u;
draw fullcircle scaled D shifted z43 withcolor RoyalBlue withpen pencircle scaled 1.0;

x44 = .0696073572u;
y44 = .2673796612u;
draw fullcircle scaled D shifted z44 withcolor RoyalBlue withpen pencircle scaled 1.0;

x45 = .0806788227u;
y45 = .3319110498u;
draw fullcircle scaled D shifted z45 withcolor RoyalBlue withpen pencircle scaled 1.0;

x46 = .1149882993u;
y46 = .6801789990u;
draw fullcircle scaled D shifted z46 withcolor RoyalBlue withpen pencircle scaled 1.0;

x47 = .1349475170u;
y47 = .8243861806u;
draw fullcircle scaled D shifted z47 withcolor RoyalBlue withpen pencircle scaled 1.0;

x48 = .2060965105u;
y48 = .9457282471u;
draw fullcircle scaled D shifted z48 withcolor RoyalBlue withpen pencircle scaled 1.0;

x49 = .0470830083u;
y49 = .7259044026u;
draw fullcircle scaled D shifted z49 withcolor RoyalBlue withpen pencircle scaled 1.0;

x50 = .0555721452u;
y50 = .9383563519u;
draw fullcircle scaled D shifted z50 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.9444423212u,        0.0616302625u)--(        0.8888888889u,        0.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9444423212u,        0.0616302625u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9444423212u,        0.0616302625u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8888888889u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8888888889u,        0.0000000000u) withcolor Yellow;
draw (        0.9444423212u,        0.0616302625u)--(        1.0000000000u,        0.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9444423212u,        0.0616302625u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9444423212u,        0.0616302625u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.9444423212u,        0.0616302625u)--(        0.9123473232u,        0.1381435561u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8392771146u,        0.0988355720u)--(        0.9123473232u,        0.1381435561u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9444423212u,        0.0616302625u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9444423212u,        0.0616302625u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8392771146u,        0.0988355720u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8392771146u,        0.0988355720u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9123473232u,        0.1381435561u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9123473232u,        0.1381435561u) withcolor MidnightBlue;
draw (        0.9444423212u,        0.0616302625u)--(        0.9123543494u,        0.1381465030u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9274937486u,        0.2197257381u)--(        0.9123543494u,        0.1381465030u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9444423212u,        0.0616302625u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9444423212u,        0.0616302625u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9274937486u,        0.2197257381u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9274937486u,        0.2197257381u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9123543494u,        0.1381465030u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9123543494u,        0.1381465030u) withcolor MidnightBlue;
draw (        0.8392771146u,        0.0988355720u)--(        0.8544260809u,        0.1804130310u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9274937486u,        0.2197257381u)--(        0.8544260809u,        0.1804130310u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8392771146u,        0.0988355720u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8392771146u,        0.0988355720u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9274937486u,        0.2197257381u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9274937486u,        0.2197257381u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8544260809u,        0.1804130310u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8544260809u,        0.1804130310u) withcolor MidnightBlue;
draw (        0.8392771146u,        0.0988355720u)--(        0.9123447823u,        0.1381482791u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9274937486u,        0.2197257381u)--(        0.9123447823u,        0.1381482791u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8392771146u,        0.0988355720u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8392771146u,        0.0988355720u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9274937486u,        0.2197257381u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9274937486u,        0.2197257381u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9123447823u,        0.1381482791u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9123447823u,        0.1381482791u) withcolor MidnightBlue;
draw (        0.8392771146u,        0.0988355720u)--(        0.7945006607u,        0.0289825679u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7142287775u,        0.0499777522u)--(        0.7945006607u,        0.0289825679u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8392771146u,        0.0988355720u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8392771146u,        0.0988355720u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7142287775u,        0.0499777522u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7142287775u,        0.0499777522u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7945006607u,        0.0289825679u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7945006607u,        0.0289825679u) withcolor Green;
draw (        0.8392771146u,        0.0988355720u)--(        0.7590052313u,        0.1198307563u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7142287775u,        0.0499777522u)--(        0.7590052313u,        0.1198307563u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8392771146u,        0.0988355720u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8392771146u,        0.0988355720u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7142287775u,        0.0499777522u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7142287775u,        0.0499777522u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7590052313u,        0.1198307563u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7590052313u,        0.1198307563u) withcolor MidnightBlue;
draw (        0.8392771146u,        0.0988355720u)--(        0.7590054454u,        0.1198315749u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7741497143u,        0.2014099061u)--(        0.7590054454u,        0.1198315749u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8392771146u,        0.0988355720u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8392771146u,        0.0988355720u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7741497143u,        0.2014099061u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7741497143u,        0.2014099061u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7590054454u,        0.1198315749u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7590054454u,        0.1198315749u) withcolor MidnightBlue;
draw (        0.8392771146u,        0.0988355720u)--(        0.8544213834u,        0.1804139032u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7741497143u,        0.2014099061u)--(        0.8544213834u,        0.1804139032u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8392771146u,        0.0988355720u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8392771146u,        0.0988355720u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7741497143u,        0.2014099061u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7741497143u,        0.2014099061u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8544213834u,        0.1804139032u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8544213834u,        0.1804139032u) withcolor MidnightBlue;
draw (        0.9274937486u,        0.2197257381u)--(        0.9999898248u,        0.2600827938u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9404835244u,        0.3179047032u)--(        0.9999898248u,        0.2600827938u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9274937486u,        0.2197257381u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9274937486u,        0.2197257381u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9404835244u,        0.3179047032u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9404835244u,        0.3179047032u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9999898248u,        0.2600827938u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9999898248u,        0.2600827938u) withcolor Green;
draw (        0.9274937486u,        0.2197257381u)--(        0.8544227665u,        0.1804191919u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7741497143u,        0.2014099061u)--(        0.8544227665u,        0.1804191919u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9274937486u,        0.2197257381u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9274937486u,        0.2197257381u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7741497143u,        0.2014099061u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7741497143u,        0.2014099061u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8544227665u,        0.1804191919u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8544227665u,        0.1804191919u) withcolor MidnightBlue;
draw (        0.9274937486u,        0.2197257381u)--(        0.8472206964u,        0.2407164524u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7741497143u,        0.2014099061u)--(        0.8472206964u,        0.2407164524u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9274937486u,        0.2197257381u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9274937486u,        0.2197257381u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7741497143u,        0.2014099061u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7741497143u,        0.2014099061u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8472206964u,        0.2407164524u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8472206964u,        0.2407164524u) withcolor MidnightBlue;
draw (        0.9274937486u,        0.2197257381u)--(        0.8472234026u,        0.2407267990u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8202322263u,        0.3191860203u)--(        0.8472234026u,        0.2407267990u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9274937486u,        0.2197257381u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9274937486u,        0.2197257381u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8202322263u,        0.3191860203u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8202322263u,        0.3191860203u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8472234026u,        0.2407267990u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8472234026u,        0.2407267990u) withcolor MidnightBlue;
draw (        0.9404835244u,        0.3179047032u)--(        0.8809670487u,        0.3757161392u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8202322263u,        0.3191860203u)--(        0.8809670487u,        0.3757161392u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9404835244u,        0.3179047032u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9404835244u,        0.3179047032u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8202322263u,        0.3191860203u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8202322263u,        0.3191860203u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8809670487u,        0.3757161392u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8809670487u,        0.3757161392u) withcolor Green;
draw (        0.7142287775u,        0.0499777522u)--(        0.7293751455u,        0.1315556937u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7741497143u,        0.2014099061u)--(        0.7293751455u,        0.1315556937u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7142287775u,        0.0499777522u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7142287775u,        0.0499777522u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7741497143u,        0.2014099061u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7741497143u,        0.2014099061u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7293751455u,        0.1315556937u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7293751455u,        0.1315556937u) withcolor MidnightBlue;
draw (        0.7142287775u,        0.0499777522u)--(        0.7590033463u,        0.1198319646u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7741497143u,        0.2014099061u)--(        0.7590033463u,        0.1198319646u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7142287775u,        0.0499777522u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7142287775u,        0.0499777522u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7741497143u,        0.2014099061u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7741497143u,        0.2014099061u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7590033463u,        0.1198319646u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7590033463u,        0.1198319646u) withcolor MidnightBlue;
draw (        0.7142287775u,        0.0499777522u)--(        0.7293728016u,        0.1315561288u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6470408194u,        0.1418429022u)--(        0.7293728016u,        0.1315561288u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7142287775u,        0.0499777522u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7142287775u,        0.0499777522u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6470408194u,        0.1418429022u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6470408194u,        0.1418429022u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7293728016u,        0.1315561288u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7293728016u,        0.1315561288u) withcolor MidnightBlue;
draw (        0.7741497143u,        0.2014099061u)--(        0.8472164218u,        0.2407243978u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8202322263u,        0.3191860203u)--(        0.8472164218u,        0.2407243978u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7741497143u,        0.2014099061u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7741497143u,        0.2014099061u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8202322263u,        0.3191860203u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8202322263u,        0.3191860203u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8472164218u,        0.2407243978u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8472164218u,        0.2407243978u) withcolor MidnightBlue;
draw (        0.7741497143u,        0.2014099061u)--(        0.7293729251u,        0.1315571169u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6470408194u,        0.1418429022u)--(        0.7293729251u,        0.1315571169u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7741497143u,        0.2014099061u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7741497143u,        0.2014099061u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6470408194u,        0.1418429022u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6470408194u,        0.1418429022u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7293729251u,        0.1315571169u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7293729251u,        0.1315571169u) withcolor MidnightBlue;
draw (        0.7741497143u,        0.2014099061u)--(        0.6918176086u,        0.2116956913u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6470408194u,        0.1418429022u)--(        0.6918176086u,        0.2116956913u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7741497143u,        0.2014099061u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7741497143u,        0.2014099061u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6470408194u,        0.1418429022u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6470408194u,        0.1418429022u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6918176086u,        0.2116956913u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6918176086u,        0.2116956913u) withcolor MidnightBlue;
draw (        0.7741497143u,        0.2014099061u)--(        0.6918190363u,        0.2117071130u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6727106157u,        0.2924489315u)--(        0.6918190363u,        0.2117071130u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7741497143u,        0.2014099061u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7741497143u,        0.2014099061u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6727106157u,        0.2924489315u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6727106157u,        0.2924489315u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6918190363u,        0.2117071130u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6918190363u,        0.2117071130u) withcolor MidnightBlue;
draw (        0.8202322263u,        0.3191860203u)--(        0.7401285604u,        0.3408141436u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6727106157u,        0.2924489315u)--(        0.7401285604u,        0.3408141436u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8202322263u,        0.3191860203u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8202322263u,        0.3191860203u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6727106157u,        0.2924489315u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6727106157u,        0.2924489315u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7401285604u,        0.3408141436u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7401285604u,        0.3408141436u) withcolor MidnightBlue;
draw (        0.8202322263u,        0.3191860203u)--(        0.7401298971u,        0.3408190937u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7274447211u,        0.4228157968u)--(        0.7401298971u,        0.3408190937u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8202322263u,        0.3191860203u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8202322263u,        0.3191860203u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7274447211u,        0.4228157968u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7274447211u,        0.4228157968u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7401298971u,        0.3408190937u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7401298971u,        0.3408190937u) withcolor MidnightBlue;
draw (        0.8202322263u,        0.3191860203u)--(        0.8075470503u,        0.4011827233u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7274447211u,        0.4228157968u)--(        0.8075470503u,        0.4011827233u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8202322263u,        0.3191860203u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8202322263u,        0.3191860203u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7274447211u,        0.4228157968u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7274447211u,        0.4228157968u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8075470503u,        0.4011827233u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8075470503u,        0.4011827233u) withcolor Green;
draw (        0.9422783475u,        0.7829731243u)--(        0.9142595922u,        0.8610712828u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9444333216u,        0.9383624146u)--(        0.9142595922u,        0.8610712828u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9422783475u,        0.7829731243u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9422783475u,        0.7829731243u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9444333216u,        0.9383624146u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9444333216u,        0.9383624146u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9142595922u,        0.8610712828u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9142595922u,        0.8610712828u) withcolor MidnightBlue;
draw (        0.9422783475u,        0.7829731243u)--(        0.8608868060u,        0.7668551225u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8073944162u,        0.7034284938u)--(        0.8608868060u,        0.7668551225u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9422783475u,        0.7829731243u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9422783475u,        0.7829731243u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8073944162u,        0.7034284938u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8073944162u,        0.7034284938u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8608868060u,        0.7668551225u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8608868060u,        0.7668551225u) withcolor MidnightBlue;
draw (        0.9422783475u,        0.7829731243u)--(        0.8608850707u,        0.7668638879u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8328491722u,        0.8449558938u)--(        0.8608850707u,        0.7668638879u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9422783475u,        0.7829731243u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9422783475u,        0.7829731243u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8328491722u,        0.8449558938u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8328491722u,        0.8449558938u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8608850707u,        0.7668638879u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8608850707u,        0.7668638879u) withcolor MidnightBlue;
draw (        0.9422783475u,        0.7829731243u)--(        0.9142424490u,        0.8610651303u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8328491722u,        0.8449558938u)--(        0.9142424490u,        0.8610651303u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9422783475u,        0.7829731243u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9422783475u,        0.7829731243u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8328491722u,        0.8449558938u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8328491722u,        0.8449558938u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9142424490u,        0.8610651303u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9142424490u,        0.8610651303u) withcolor MidnightBlue;
draw (        0.9444333216u,        0.9383624146u)--(        1.0000000000u,        1.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9444333216u,        0.9383624146u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9444333216u,        0.9383624146u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        1.0000000000u) withcolor Yellow;
draw (        0.9444333216u,        0.9383624146u)--(        0.8888888889u,        1.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9444333216u,        0.9383624146u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9444333216u,        0.9383624146u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8888888889u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8888888889u,        1.0000000000u) withcolor Yellow;
draw (        0.9444333216u,        0.9383624146u)--(        0.9142396918u,        0.8610790546u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8328491722u,        0.8449558938u)--(        0.9142396918u,        0.8610790546u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9444333216u,        0.9383624146u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9444333216u,        0.9383624146u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8328491722u,        0.8449558938u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8328491722u,        0.8449558938u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9142396918u,        0.8610790546u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9142396918u,        0.8610790546u) withcolor MidnightBlue;
draw (        0.9444333216u,        0.9383624146u)--(        0.8630428019u,        0.9222392538u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8328491722u,        0.8449558938u)--(        0.8630428019u,        0.9222392538u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9444333216u,        0.9383624146u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9444333216u,        0.9383624146u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8328491722u,        0.8449558938u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8328491722u,        0.8449558938u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8630428019u,        0.9222392538u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8630428019u,        0.9222392538u) withcolor MidnightBlue;
draw (        0.9444333216u,        0.9383624146u)--(        0.8630403751u,        0.9222515095u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7897360084u,        0.9611210712u)--(        0.8630403751u,        0.9222515095u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9444333216u,        0.9383624146u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9444333216u,        0.9383624146u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7897360084u,        0.9611210712u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7897360084u,        0.9611210712u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8630403751u,        0.9222515095u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8630403751u,        0.9222515095u) withcolor MidnightBlue;
draw (        0.6470408194u,        0.1418429022u)--(        0.6279464583u,        0.2225880468u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6727106157u,        0.2924489315u)--(        0.6279464583u,        0.2225880468u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6470408194u,        0.1418429022u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6470408194u,        0.1418429022u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6727106157u,        0.2924489315u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6727106157u,        0.2924489315u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6279464583u,        0.2225880468u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6279464583u,        0.2225880468u) withcolor MidnightBlue;
draw (        0.6470408194u,        0.1418429022u)--(        0.6918049768u,        0.2117037869u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6727106157u,        0.2924489315u)--(        0.6918049768u,        0.2117037869u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6470408194u,        0.1418429022u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6470408194u,        0.1418429022u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6727106157u,        0.2924489315u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6727106157u,        0.2924489315u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6918049768u,        0.2117037869u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6918049768u,        0.2117037869u) withcolor MidnightBlue;
draw (        0.6470408194u,        0.1418429022u)--(        0.5812898284u,        0.1924510062u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5621819437u,        0.2731929515u)--(        0.5812898284u,        0.1924510062u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6470408194u,        0.1418429022u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6470408194u,        0.1418429022u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5621819437u,        0.2731929515u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5621819437u,        0.2731929515u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5812898284u,        0.1924510062u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5812898284u,        0.1924510062u) withcolor Green;
draw (        0.6470408194u,        0.1418429022u)--(        0.6279329347u,        0.2225848476u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5621819437u,        0.2731929515u)--(        0.6279329347u,        0.2225848476u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6470408194u,        0.1418429022u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6470408194u,        0.1418429022u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5621819437u,        0.2731929515u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5621819437u,        0.2731929515u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6279329347u,        0.2225848476u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6279329347u,        0.2225848476u) withcolor MidnightBlue;
draw (        0.6727106157u,        0.2924489315u)--(        0.6600298396u,        0.3744463152u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7274447211u,        0.4228157968u)--(        0.6600298396u,        0.3744463152u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6727106157u,        0.2924489315u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6727106157u,        0.2924489315u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7274447211u,        0.4228157968u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7274447211u,        0.4228157968u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6600298396u,        0.3744463152u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6600298396u,        0.3744463152u) withcolor MidnightBlue;
draw (        0.6727106157u,        0.2924489315u)--(        0.7401254972u,        0.3408184132u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7274447211u,        0.4228157968u)--(        0.7401254972u,        0.3408184132u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6727106157u,        0.2924489315u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6727106157u,        0.2924489315u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7274447211u,        0.4228157968u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7274447211u,        0.4228157968u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7401254972u,        0.3408184132u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7401254972u,        0.3408184132u) withcolor MidnightBlue;
draw (        0.6727106157u,        0.2924489315u)--(        0.6279390503u,        0.2225927941u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5621819437u,        0.2731929515u)--(        0.6279390503u,        0.2225927941u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6727106157u,        0.2924489315u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6727106157u,        0.2924489315u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5621819437u,        0.2731929515u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5621819437u,        0.2731929515u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6279390503u,        0.2225927941u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6279390503u,        0.2225927941u) withcolor MidnightBlue;
draw (        0.6727106157u,        0.2924489315u)--(        0.6600114356u,        0.3744434669u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5887818314u,        0.4169959813u)--(        0.6600114356u,        0.3744434669u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6727106157u,        0.2924489315u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6727106157u,        0.2924489315u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5887818314u,        0.4169959813u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5887818314u,        0.4169959813u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6600114356u,        0.3744434669u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6600114356u,        0.3744434669u) withcolor MidnightBlue;
draw (        0.7274447211u,        0.4228157968u)--(        0.8075499103u,        0.4444444444u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7274447211u,        0.4228157968u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7274447211u,        0.4228157968u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8075499103u,        0.4444444444u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8075499103u,        0.4444444444u) withcolor Yellow;
draw (        0.7274447211u,        0.4228157968u)--(        0.6961223275u,        0.4996486066u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7433895439u,        0.5678408587u)--(        0.6961223275u,        0.4996486066u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7274447211u,        0.4228157968u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7274447211u,        0.4228157968u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7433895439u,        0.5678408587u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7433895439u,        0.5678408587u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6961223275u,        0.4996486066u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6961223275u,        0.4996486066u) withcolor MidnightBlue;
draw (        0.7274447211u,        0.4228157968u)--(        0.6600207252u,        0.3744590208u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5887818314u,        0.4169959813u)--(        0.6600207252u,        0.3744590208u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7274447211u,        0.4228157968u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7274447211u,        0.4228157968u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5887818314u,        0.4169959813u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5887818314u,        0.4169959813u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6600207252u,        0.3744590208u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6600207252u,        0.3744590208u) withcolor MidnightBlue;
draw (        0.7274447211u,        0.4228157968u)--(        0.6961079256u,        0.4996427338u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6236325111u,        0.5400368830u)--(        0.6961079256u,        0.4996427338u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7274447211u,        0.4228157968u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7274447211u,        0.4228157968u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6236325111u,        0.5400368830u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6236325111u,        0.5400368830u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6961079256u,        0.4996427338u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6961079256u,        0.4996427338u) withcolor MidnightBlue;
draw (        0.7433895439u,        0.5678408587u)--(        0.7432386189u,        0.6508128433u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8073944162u,        0.7034284938u)--(        0.7432386189u,        0.6508128433u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7433895439u,        0.5678408587u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7433895439u,        0.5678408587u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8073944162u,        0.7034284938u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8073944162u,        0.7034284938u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7432386189u,        0.6508128433u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7432386189u,        0.6508128433u) withcolor MidnightBlue;
draw (        0.7433895439u,        0.5678408587u)--(        0.8075453412u,        0.6204565092u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8073944162u,        0.7034284938u)--(        0.8075453412u,        0.6204565092u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7433895439u,        0.5678408587u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7433895439u,        0.5678408587u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8073944162u,        0.7034284938u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8073944162u,        0.7034284938u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8075453412u,        0.6204565092u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8075453412u,        0.6204565092u) withcolor Green;
draw (        0.7433895439u,        0.5678408587u)--(        0.6961142994u,        0.4996541718u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6236325111u,        0.5400368830u)--(        0.6961142994u,        0.4996541718u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7433895439u,        0.5678408587u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7433895439u,        0.5678408587u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6236325111u,        0.5400368830u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6236325111u,        0.5400368830u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6961142994u,        0.4996541718u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6961142994u,        0.4996541718u) withcolor MidnightBlue;
draw (        0.7433895439u,        0.5678408587u)--(        0.7432351597u,        0.6508128369u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6626011231u,        0.6703710872u)--(        0.7432351597u,        0.6508128369u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7433895439u,        0.5678408587u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7433895439u,        0.5678408587u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6626011231u,        0.6703710872u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6626011231u,        0.6703710872u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7432351597u,        0.6508128369u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7432351597u,        0.6508128369u) withcolor MidnightBlue;
draw (        0.8073944162u,        0.7034284938u)--(        0.7793645913u,        0.7815226800u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8328491722u,        0.8449558938u)--(        0.7793645913u,        0.7815226800u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8073944162u,        0.7034284938u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8073944162u,        0.7034284938u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8328491722u,        0.8449558938u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8328491722u,        0.8449558938u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7793645913u,        0.7815226800u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7793645913u,        0.7815226800u) withcolor MidnightBlue;
draw (        0.8073944162u,        0.7034284938u)--(        0.8608789970u,        0.7668617076u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8328491722u,        0.8449558938u)--(        0.8608789970u,        0.7668617076u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8073944162u,        0.7034284938u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8073944162u,        0.7034284938u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8328491722u,        0.8449558938u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8328491722u,        0.8449558938u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8608789970u,        0.7668617076u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8608789970u,        0.7668617076u) withcolor MidnightBlue;
draw (        0.8073944162u,        0.7034284938u)--(        0.7432359504u,        0.6508160972u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6626011231u,        0.6703710872u)--(        0.7432359504u,        0.6508160972u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8073944162u,        0.7034284938u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8073944162u,        0.7034284938u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6626011231u,        0.6703710872u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6626011231u,        0.6703710872u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7432359504u,        0.6508160972u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7432359504u,        0.6508160972u) withcolor MidnightBlue;
draw (        0.8073944162u,        0.7034284938u)--(        0.7267595889u,        0.7229834838u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6626011231u,        0.6703710872u)--(        0.7267595889u,        0.7229834838u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8073944162u,        0.7034284938u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8073944162u,        0.7034284938u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6626011231u,        0.6703710872u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6626011231u,        0.6703710872u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7267595889u,        0.7229834838u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7267595889u,        0.7229834838u) withcolor MidnightBlue;
draw (        0.8073944162u,        0.7034284938u)--(        0.7267618730u,        0.7229928999u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6987225336u,        0.8010836705u)--(        0.7267618730u,        0.7229928999u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8073944162u,        0.7034284938u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8073944162u,        0.7034284938u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6987225336u,        0.8010836705u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6987225336u,        0.8010836705u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7267618730u,        0.7229928999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7267618730u,        0.7229928999u) withcolor MidnightBlue;
draw (        0.8073944162u,        0.7034284938u)--(        0.7793550767u,        0.7815192644u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6987225336u,        0.8010836705u)--(        0.7793550767u,        0.7815192644u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8073944162u,        0.7034284938u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8073944162u,        0.7034284938u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6987225336u,        0.8010836705u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6987225336u,        0.8010836705u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7793550767u,        0.7815192644u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7793550767u,        0.7815192644u) withcolor MidnightBlue;
draw (        0.8328491722u,        0.8449558938u)--(        0.8630354074u,        0.9222421423u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7897360084u,        0.9611210712u)--(        0.8630354074u,        0.9222421423u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8328491722u,        0.8449558938u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8328491722u,        0.8449558938u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7897360084u,        0.9611210712u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7897360084u,        0.9611210712u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8630354074u,        0.9222421423u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8630354074u,        0.9222421423u) withcolor MidnightBlue;
draw (        0.8328491722u,        0.8449558938u)--(        0.7793573803u,        0.7815287608u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6987225336u,        0.8010836705u)--(        0.7793573803u,        0.7815287608u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8328491722u,        0.8449558938u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8328491722u,        0.8449558938u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6987225336u,        0.8010836705u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6987225336u,        0.8010836705u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7793573803u,        0.7815287608u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7793573803u,        0.7815287608u) withcolor MidnightBlue;
draw (        0.8328491722u,        0.8449558938u)--(        0.7522143254u,        0.8645108035u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6987225336u,        0.8010836705u)--(        0.7522143254u,        0.8645108035u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8328491722u,        0.8449558938u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8328491722u,        0.8449558938u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6987225336u,        0.8010836705u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6987225336u,        0.8010836705u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7522143254u,        0.8645108035u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7522143254u,        0.8645108035u) withcolor MidnightBlue;
draw (        0.8328491722u,        0.8449558938u)--(        0.7522174580u,        0.8645237163u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6913423248u,        0.9209027120u)--(        0.7522174580u,        0.8645237163u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8328491722u,        0.8449558938u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8328491722u,        0.8449558938u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6913423248u,        0.9209027120u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6913423248u,        0.9209027120u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7522174580u,        0.8645237163u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7522174580u,        0.8645237163u) withcolor MidnightBlue;
draw (        0.7897360084u,        0.9611210712u)--(        0.7164316416u,        0.9999906327u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6913423248u,        0.9209027120u)--(        0.7164316416u,        0.9999906327u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7897360084u,        0.9611210712u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7897360084u,        0.9611210712u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6913423248u,        0.9209027120u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6913423248u,        0.9209027120u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7164316416u,        0.9999906327u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7164316416u,        0.9999906327u) withcolor Green;
draw (        0.5621819437u,        0.2731929515u)--(        0.5369230303u,        0.3522268702u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5887818314u,        0.4169959813u)--(        0.5369230303u,        0.3522268702u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5621819437u,        0.2731929515u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5621819437u,        0.2731929515u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5887818314u,        0.4169959813u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5887818314u,        0.4169959813u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5369230303u,        0.3522268702u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5369230303u,        0.3522268702u) withcolor MidnightBlue;
draw (        0.5621819437u,        0.2731929515u)--(        0.4861926223u,        0.3065095546u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4160894773u,        0.2621257505u)--(        0.4861926223u,        0.3065095546u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5621819437u,        0.2731929515u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5621819437u,        0.2731929515u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4160894773u,        0.2621257505u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4160894773u,        0.2621257505u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4861926223u,        0.3065095546u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4861926223u,        0.3065095546u) withcolor MidnightBlue;
draw (        0.5621819437u,        0.2731929515u)--(        0.4861955138u,        0.3065161488u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4609295037u,        0.3855477990u)--(        0.4861955138u,        0.3065161488u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5621819437u,        0.2731929515u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5621819437u,        0.2731929515u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4609295037u,        0.3855477990u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4609295037u,        0.3855477990u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4861955138u,        0.3065161488u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4861955138u,        0.3065161488u) withcolor MidnightBlue;
draw (        0.5621819437u,        0.2731929515u)--(        0.5369159336u,        0.3522246018u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4609295037u,        0.3855477990u)--(        0.5369159336u,        0.3522246018u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5621819437u,        0.2731929515u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5621819437u,        0.2731929515u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4609295037u,        0.3855477990u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4609295037u,        0.3855477990u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5369159336u,        0.3522246018u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5369159336u,        0.3522246018u) withcolor MidnightBlue;
draw (        0.5887818314u,        0.4169959813u)--(        0.5553317268u,        0.4929266307u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6236325111u,        0.5400368830u)--(        0.5553317268u,        0.4929266307u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5887818314u,        0.4169959813u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5887818314u,        0.4169959813u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6236325111u,        0.5400368830u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6236325111u,        0.5400368830u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5553317268u,        0.4929266307u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5553317268u,        0.4929266307u) withcolor MidnightBlue;
draw (        0.5887818314u,        0.4169959813u)--(        0.5369185140u,        0.3522304866u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4609295037u,        0.3855477990u)--(        0.5369185140u,        0.3522304866u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5887818314u,        0.4169959813u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5887818314u,        0.4169959813u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4609295037u,        0.3855477990u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4609295037u,        0.3855477990u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5369185140u,        0.3522304866u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5369185140u,        0.3522304866u) withcolor MidnightBlue;
draw (        0.5887818314u,        0.4169959813u)--(        0.5127928210u,        0.4503132938u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4609295037u,        0.3855477990u)--(        0.5127928210u,        0.4503132938u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5887818314u,        0.4169959813u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5887818314u,        0.4169959813u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4609295037u,        0.3855477990u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4609295037u,        0.3855477990u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5127928210u,        0.4503132938u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5127928210u,        0.4503132938u) withcolor MidnightBlue;
draw (        0.5887818314u,        0.4169959813u)--(        0.5127960031u,        0.4503205504u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4793352746u,        0.5262465187u)--(        0.5127960031u,        0.4503205504u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5887818314u,        0.4169959813u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5887818314u,        0.4169959813u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4793352746u,        0.5262465187u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4793352746u,        0.5262465187u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5127960031u,        0.4503205504u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5127960031u,        0.4503205504u) withcolor MidnightBlue;
draw (        0.5887818314u,        0.4169959813u)--(        0.5553211029u,        0.4929219496u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4793352746u,        0.5262465187u)--(        0.5553211029u,        0.4929219496u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5887818314u,        0.4169959813u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5887818314u,        0.4169959813u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4793352746u,        0.5262465187u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4793352746u,        0.5262465187u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5553211029u,        0.4929219496u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5553211029u,        0.4929219496u) withcolor MidnightBlue;
draw (        0.6236325111u,        0.5400368830u)--(        0.5975899556u,        0.6188160573u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6626011231u,        0.6703710872u)--(        0.5975899556u,        0.6188160573u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6236325111u,        0.5400368830u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6236325111u,        0.5400368830u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6626011231u,        0.6703710872u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6626011231u,        0.6703710872u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5975899556u,        0.6188160573u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5975899556u,        0.6188160573u) withcolor MidnightBlue;
draw (        0.6236325111u,        0.5400368830u)--(        0.5553264868u,        0.4929342286u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4793352746u,        0.5262465187u)--(        0.5553264868u,        0.4929342286u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6236325111u,        0.5400368830u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6236325111u,        0.5400368830u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4793352746u,        0.5262465187u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4793352746u,        0.5262465187u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5553264868u,        0.4929342286u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5553264868u,        0.4929342286u) withcolor MidnightBlue;
draw (        0.6236325111u,        0.5400368830u)--(        0.5975783417u,        0.6188122171u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5192180895u,        0.6460893848u)--(        0.5975783417u,        0.6188122171u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6236325111u,        0.5400368830u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6236325111u,        0.5400368830u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5192180895u,        0.6460893848u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5192180895u,        0.6460893848u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5975783417u,        0.6188122171u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5975783417u,        0.6188122171u) withcolor MidnightBlue;
draw (        0.6626011231u,        0.6703710872u)--(        0.6345695141u,        0.7484646330u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6987225336u,        0.8010836705u)--(        0.6345695141u,        0.7484646330u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6626011231u,        0.6703710872u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6626011231u,        0.6703710872u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6987225336u,        0.8010836705u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6987225336u,        0.8010836705u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6345695141u,        0.7484646330u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6345695141u,        0.7484646330u) withcolor MidnightBlue;
draw (        0.6626011231u,        0.6703710872u)--(        0.7267541426u,        0.7229901246u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6987225336u,        0.8010836705u)--(        0.7267541426u,        0.7229901246u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6626011231u,        0.6703710872u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6626011231u,        0.6703710872u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6987225336u,        0.8010836705u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6987225336u,        0.8010836705u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7267541426u,        0.7229901246u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7267541426u,        0.7229901246u) withcolor MidnightBlue;
draw (        0.6626011231u,        0.6703710872u)--(        0.5975828130u,        0.6188250655u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5192180895u,        0.6460893848u)--(        0.5975828130u,        0.6188250655u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6626011231u,        0.6703710872u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6626011231u,        0.6703710872u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5192180895u,        0.6460893848u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5192180895u,        0.6460893848u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5975828130u,        0.6188250655u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5975828130u,        0.6188250655u) withcolor MidnightBlue;
draw (        0.6626011231u,        0.6703710872u)--(        0.6345670272u,        0.7484637403u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5518470317u,        0.7549271244u)--(        0.6345670272u,        0.7484637403u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6626011231u,        0.6703710872u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6626011231u,        0.6703710872u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5518470317u,        0.7549271244u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5518470317u,        0.7549271244u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6345670272u,        0.7484637403u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6345670272u,        0.7484637403u) withcolor MidnightBlue;
draw (        0.6987225336u,        0.8010836705u)--(        0.6378554895u,        0.8574713992u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6913423248u,        0.9209027120u)--(        0.6378554895u,        0.8574713992u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6987225336u,        0.8010836705u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6987225336u,        0.8010836705u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6913423248u,        0.9209027120u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6913423248u,        0.9209027120u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6378554895u,        0.8574713992u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6378554895u,        0.8574713992u) withcolor Green;
draw (        0.6987225336u,        0.8010836705u)--(        0.7522093689u,        0.8645149833u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6913423248u,        0.9209027120u)--(        0.7522093689u,        0.8645149833u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6987225336u,        0.8010836705u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6987225336u,        0.8010836705u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6913423248u,        0.9209027120u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6913423248u,        0.9209027120u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7522093689u,        0.8645149833u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7522093689u,        0.8645149833u) withcolor MidnightBlue;
draw (        0.6987225336u,        0.8010836705u)--(        0.6345673072u,        0.7484673239u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5518470317u,        0.7549271244u)--(        0.6345673072u,        0.7484673239u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6987225336u,        0.8010836705u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6987225336u,        0.8010836705u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5518470317u,        0.7549271244u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5518470317u,        0.7549271244u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6345673072u,        0.7484673239u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6345673072u,        0.7484673239u) withcolor MidnightBlue;
draw (        0.6987225336u,        0.8010836705u)--(        0.6160022581u,        0.8075434710u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5518470317u,        0.7549271244u)--(        0.6160022581u,        0.8075434710u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6987225336u,        0.8010836705u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6987225336u,        0.8010836705u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5518470317u,        0.7549271244u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5518470317u,        0.7549271244u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6160022581u,        0.8075434710u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6160022581u,        0.8075434710u) withcolor Green;
draw (        0.4160894773u,        0.2621257505u)--(        0.4861896366u,        0.3065142701u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4609295037u,        0.3855477990u)--(        0.4861896366u,        0.3065142701u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4160894773u,        0.2621257505u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4160894773u,        0.2621257505u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4609295037u,        0.3855477990u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4609295037u,        0.3855477990u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4861896366u,        0.3065142701u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4861896366u,        0.3065142701u) withcolor MidnightBlue;
draw (        0.4160894773u,        0.2621257505u)--(        0.3710347644u,        0.1924518971u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2882406904u,        0.1870191904u)--(        0.3710347644u,        0.1924518971u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4160894773u,        0.2621257505u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4160894773u,        0.2621257505u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2882406904u,        0.1870191904u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2882406904u,        0.1870191904u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3710347644u,        0.1924518971u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3710347644u,        0.1924518971u) withcolor Green;
draw (        0.4160894773u,        0.2621257505u)--(        0.3332954033u,        0.2566930438u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2882406904u,        0.1870191904u)--(        0.3332954033u,        0.2566930438u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4160894773u,        0.2621257505u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4160894773u,        0.2621257505u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2882406904u,        0.1870191904u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2882406904u,        0.1870191904u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3332954033u,        0.2566930438u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3332954033u,        0.2566930438u) withcolor MidnightBlue;
draw (        0.4160894773u,        0.2621257505u)--(        0.3332952125u,        0.2566959521u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3403610580u,        0.3393666653u)--(        0.3332952125u,        0.2566959521u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4160894773u,        0.2621257505u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4160894773u,        0.2621257505u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3403610580u,        0.3393666653u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3403610580u,        0.3393666653u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3332952125u,        0.2566959521u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3332952125u,        0.2566959521u) withcolor MidnightBlue;
draw (        0.4609295037u,        0.3855477990u)--(        0.5127881049u,        0.4503170702u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4793352746u,        0.5262465187u)--(        0.5127881049u,        0.4503170702u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4609295037u,        0.3855477990u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4609295037u,        0.3855477990u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4793352746u,        0.5262465187u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4793352746u,        0.5262465187u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5127881049u,        0.4503170702u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5127881049u,        0.4503170702u) withcolor MidnightBlue;
draw (        0.4609295037u,        0.3855477990u)--(        0.3820008253u,        0.4111336687u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3403610580u,        0.3393666653u)--(        0.3820008253u,        0.4111336687u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4609295037u,        0.3855477990u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4609295037u,        0.3855477990u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3403610580u,        0.3393666653u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3403610580u,        0.3393666653u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3820008253u,        0.4111336687u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3820008253u,        0.4111336687u) withcolor MidnightBlue;
draw (        0.4609295037u,        0.3855477990u)--(        0.3820045671u,        0.4111452087u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3521242518u,        0.4885502458u)--(        0.3820045671u,        0.4111452087u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4609295037u,        0.3855477990u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4609295037u,        0.3855477990u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3521242518u,        0.4885502458u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3521242518u,        0.4885502458u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3820045671u,        0.4111452087u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3820045671u,        0.4111452087u) withcolor MidnightBlue;
draw (        0.4793352746u,        0.5262465187u)--(        0.4482142994u,        0.6031611330u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5192180895u,        0.6460893848u)--(        0.4482142994u,        0.6031611330u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4793352746u,        0.5262465187u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4793352746u,        0.5262465187u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5192180895u,        0.6460893848u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5192180895u,        0.6460893848u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4482142994u,        0.6031611330u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4482142994u,        0.6031611330u) withcolor MidnightBlue;
draw (        0.4793352746u,        0.5262465187u)--(        0.4015709920u,        0.5551790151u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3521242518u,        0.4885502458u)--(        0.4015709920u,        0.5551790151u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4793352746u,        0.5262465187u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4793352746u,        0.5262465187u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3521242518u,        0.4885502458u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3521242518u,        0.4885502458u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4015709920u,        0.5551790151u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4015709920u,        0.5551790151u) withcolor MidnightBlue;
draw (        0.4793352746u,        0.5262465187u)--(        0.4015755415u,        0.5551912404u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3704458338u,        0.6321023209u)--(        0.4015755415u,        0.5551912404u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4793352746u,        0.5262465187u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4793352746u,        0.5262465187u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3704458338u,        0.6321023209u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3704458338u,        0.6321023209u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4015755415u,        0.5551912404u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4015755415u,        0.5551912404u) withcolor MidnightBlue;
draw (        0.4793352746u,        0.5262465187u)--(        0.4482055668u,        0.6031575991u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3704458338u,        0.6321023209u)--(        0.4482055668u,        0.6031575991u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4793352746u,        0.5262465187u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4793352746u,        0.5262465187u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3704458338u,        0.6321023209u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3704458338u,        0.6321023209u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4482055668u,        0.6031575991u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4482055668u,        0.6031575991u) withcolor MidnightBlue;
draw (        0.5192180895u,        0.6460893848u)--(        0.4776081282u,        0.7178736736u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5518470317u,        0.7549271244u)--(        0.4776081282u,        0.7178736736u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5192180895u,        0.6460893848u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5192180895u,        0.6460893848u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5518470317u,        0.7549271244u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5518470317u,        0.7549271244u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4776081282u,        0.7178736736u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4776081282u,        0.7178736736u) withcolor MidnightBlue;
draw (        0.5192180895u,        0.6460893848u)--(        0.4482097083u,        0.6031687276u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3704458338u,        0.6321023209u)--(        0.4482097083u,        0.6031687276u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5192180895u,        0.6460893848u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5192180895u,        0.6460893848u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3704458338u,        0.6321023209u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3704458338u,        0.6321023209u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4482097083u,        0.6031687276u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4482097083u,        0.6031687276u) withcolor MidnightBlue;
draw (        0.5192180895u,        0.6460893848u)--(        0.4775990660u,        0.7178684199u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3974132415u,        0.7391899244u)--(        0.4775990660u,        0.7178684199u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5192180895u,        0.6460893848u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5192180895u,        0.6460893848u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3974132415u,        0.7391899244u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3974132415u,        0.7391899244u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4775990660u,        0.7178684199u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4775990660u,        0.7178684199u) withcolor MidnightBlue;
draw (        0.5518470317u,        0.7549271244u)--(        0.4776031265u,        0.7178836965u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3974132415u,        0.7391899244u)--(        0.4776031265u,        0.7178836965u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5518470317u,        0.7549271244u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5518470317u,        0.7549271244u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3974132415u,        0.7391899244u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3974132415u,        0.7391899244u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4776031265u,        0.7178836965u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4776031265u,        0.7178836965u) withcolor MidnightBlue;
draw (        0.2188782361u,        0.0550328119u)--(        0.2213194446u,        0.1379690133u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2882406904u,        0.1870191904u)--(        0.2213194446u,        0.1379690133u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2188782361u,        0.0550328119u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2188782361u,        0.0550328119u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2882406904u,        0.1870191904u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2882406904u,        0.1870191904u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2213194446u,        0.1379690133u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2213194446u,        0.1379690133u) withcolor MidnightBlue;
draw (        0.2188782361u,        0.0550328119u)--(        0.1378021087u,        0.0726689841u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0555660359u,        0.0616416844u)--(        0.1378021087u,        0.0726689841u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2188782361u,        0.0550328119u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2188782361u,        0.0550328119u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0555660359u,        0.0616416844u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0555660359u,        0.0616416844u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1378021087u,        0.0726689841u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1378021087u,        0.0726689841u) withcolor MidnightBlue;
draw (        0.2188782361u,        0.0550328119u)--(        0.1378048575u,        0.0726816162u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1402378357u,        0.1556180595u)--(        0.1378048575u,        0.0726816162u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2188782361u,        0.0550328119u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2188782361u,        0.0550328119u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1402378357u,        0.1556180595u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1402378357u,        0.1556180595u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1378048575u,        0.0726816162u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1378048575u,        0.0726816162u) withcolor MidnightBlue;
draw (        0.2188782361u,        0.0550328119u)--(        0.2213112143u,        0.1379692551u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1402378357u,        0.1556180595u)--(        0.2213112143u,        0.1379692551u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2188782361u,        0.0550328119u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2188782361u,        0.0550328119u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1402378357u,        0.1556180595u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1402378357u,        0.1556180595u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2213112143u,        0.1379692551u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2213112143u,        0.1379692551u) withcolor MidnightBlue;
draw (        0.2882406904u,        0.1870191904u)--(        0.2953114992u,        0.2696894792u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3403610580u,        0.3393666653u)--(        0.2953114992u,        0.2696894792u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2882406904u,        0.1870191904u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2882406904u,        0.1870191904u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3403610580u,        0.3393666653u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3403610580u,        0.3393666653u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2953114992u,        0.2696894792u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2953114992u,        0.2696894792u) withcolor MidnightBlue;
draw (        0.2882406904u,        0.1870191904u)--(        0.3332902491u,        0.2566963764u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3403610580u,        0.3393666653u)--(        0.3332902491u,        0.2566963764u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2882406904u,        0.1870191904u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2882406904u,        0.1870191904u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3403610580u,        0.3393666653u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3403610580u,        0.3393666653u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3332902491u,        0.2566963764u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3332902491u,        0.2566963764u) withcolor MidnightBlue;
draw (        0.2882406904u,        0.1870191904u)--(        0.2213130582u,        0.1379777277u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1402378357u,        0.1556180595u)--(        0.2213130582u,        0.1379777277u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2882406904u,        0.1870191904u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2882406904u,        0.1870191904u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1402378357u,        0.1556180595u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1402378357u,        0.1556180595u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2213130582u,        0.1379777277u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2213130582u,        0.1379777277u) withcolor MidnightBlue;
draw (        0.2882406904u,        0.1870191904u)--(        0.2071654679u,        0.2046595221u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1402378357u,        0.1556180595u)--(        0.2071654679u,        0.2046595221u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2882406904u,        0.1870191904u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2882406904u,        0.1870191904u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1402378357u,        0.1556180595u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1402378357u,        0.1556180595u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2071654679u,        0.2046595221u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2071654679u,        0.2046595221u) withcolor MidnightBlue;
draw (        0.2882406904u,        0.1870191904u)--(        0.2071677727u,        0.2046701118u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2142319893u,        0.2873409642u)--(        0.2071677727u,        0.2046701118u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2882406904u,        0.1870191904u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2882406904u,        0.1870191904u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2142319893u,        0.2873409642u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2142319893u,        0.2873409642u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2071677727u,        0.2046701118u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2071677727u,        0.2046701118u) withcolor MidnightBlue;
draw (        0.2882406904u,        0.1870191904u)--(        0.2953049070u,        0.2696900428u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2142319893u,        0.2873409642u)--(        0.2953049070u,        0.2696900428u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2882406904u,        0.1870191904u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2882406904u,        0.1870191904u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2142319893u,        0.2873409642u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2142319893u,        0.2873409642u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2953049070u,        0.2696900428u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2953049070u,        0.2696900428u) withcolor MidnightBlue;
draw (        0.3403610580u,        0.3393666653u)--(        0.3819903743u,        0.4111397315u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3521242518u,        0.4885502458u)--(        0.3819903743u,        0.4111397315u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3403610580u,        0.3393666653u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3403610580u,        0.3393666653u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3521242518u,        0.4885502458u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3521242518u,        0.4885502458u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3819903743u,        0.4111397315u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3819903743u,        0.4111397315u) withcolor MidnightBlue;
draw (        0.3403610580u,        0.3393666653u)--(        0.2953056129u,        0.2696932854u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2142319893u,        0.2873409642u)--(        0.2953056129u,        0.2696932854u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3403610580u,        0.3393666653u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3403610580u,        0.3393666653u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2142319893u,        0.2873409642u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2142319893u,        0.2873409642u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2953056129u,        0.2696932854u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2953056129u,        0.2696932854u) withcolor MidnightBlue;
draw (        0.3403610580u,        0.3393666653u)--(        0.2592874344u,        0.3570143442u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2142319893u,        0.2873409642u)--(        0.2592874344u,        0.3570143442u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3403610580u,        0.3393666653u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3403610580u,        0.3393666653u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2142319893u,        0.2873409642u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2142319893u,        0.2873409642u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2592874344u,        0.3570143442u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2592874344u,        0.3570143442u) withcolor MidnightBlue;
draw (        0.3403610580u,        0.3393666653u)--(        0.2592878333u,        0.3570161765u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2375129981u,        0.4370800858u)--(        0.2592878333u,        0.3570161765u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3403610580u,        0.3393666653u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3403610580u,        0.3393666653u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2375129981u,        0.4370800858u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2375129981u,        0.4370800858u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2592878333u,        0.3570161765u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2592878333u,        0.3570161765u) withcolor MidnightBlue;
draw (        0.3521242518u,        0.4885502458u)--(        0.3210081498u,        0.5654668318u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3704458338u,        0.6321023209u)--(        0.3210081498u,        0.5654668318u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3521242518u,        0.4885502458u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3521242518u,        0.4885502458u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3704458338u,        0.6321023209u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3704458338u,        0.6321023209u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3210081498u,        0.5654668318u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3210081498u,        0.5654668318u) withcolor MidnightBlue;
draw (        0.3521242518u,        0.4885502458u)--(        0.4015619358u,        0.5551857349u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3704458338u,        0.6321023209u)--(        0.4015619358u,        0.5551857349u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3521242518u,        0.4885502458u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3521242518u,        0.4885502458u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3704458338u,        0.6321023209u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3704458338u,        0.6321023209u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4015619358u,        0.5551857349u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4015619358u,        0.5551857349u) withcolor MidnightBlue;
draw (        0.3521242518u,        0.4885502458u)--(        0.2726125965u,        0.5122624720u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2375129981u,        0.4370800858u)--(        0.2726125965u,        0.5122624720u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3521242518u,        0.4885502458u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3521242518u,        0.4885502458u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2375129981u,        0.4370800858u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2375129981u,        0.4370800858u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2726125965u,        0.5122624720u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2726125965u,        0.5122624720u) withcolor MidnightBlue;
draw (        0.3521242518u,        0.4885502458u)--(        0.2726160293u,        0.5122739799u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2414982237u,        0.5891898766u)--(        0.2726160293u,        0.5122739799u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3521242518u,        0.4885502458u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3521242518u,        0.4885502458u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2414982237u,        0.5891898766u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2414982237u,        0.5891898766u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2726160293u,        0.5122739799u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2726160293u,        0.5122739799u) withcolor MidnightBlue;
draw (        0.3521242518u,        0.4885502458u)--(        0.3210064462u,        0.5654661426u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2414982237u,        0.5891898766u)--(        0.3210064462u,        0.5654661426u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3521242518u,        0.4885502458u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3521242518u,        0.4885502458u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2414982237u,        0.5891898766u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2414982237u,        0.5891898766u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3210064462u,        0.5654661426u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3210064462u,        0.5654661426u) withcolor MidnightBlue;
draw (        0.3704458338u,        0.6321023209u)--(        0.3238721572u,        0.7007701129u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3974132415u,        0.7391899244u)--(        0.3238721572u,        0.7007701129u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3704458338u,        0.6321023209u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3704458338u,        0.6321023209u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3974132415u,        0.7391899244u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3974132415u,        0.7391899244u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3238721572u,        0.7007701129u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3238721572u,        0.7007701129u) withcolor MidnightBlue;
draw (        0.3704458338u,        0.6321023209u)--(        0.3210069233u,        0.5654677417u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2414982237u,        0.5891898766u)--(        0.3210069233u,        0.5654677417u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3704458338u,        0.6321023209u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3704458338u,        0.6321023209u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2414982237u,        0.5891898766u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2414982237u,        0.5891898766u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3210069233u,        0.5654677417u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3210069233u,        0.5654677417u) withcolor MidnightBlue;
draw (        0.3704458338u,        0.6321023209u)--(        0.2909371341u,        0.6558244558u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2414982237u,        0.5891898766u)--(        0.2909371341u,        0.6558244558u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3704458338u,        0.6321023209u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3704458338u,        0.6321023209u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2414982237u,        0.5891898766u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2414982237u,        0.5891898766u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2909371341u,        0.6558244558u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2909371341u,        0.6558244558u) withcolor MidnightBlue;
draw (        0.3704458338u,        0.6321023209u)--(        0.2909382009u,        0.6558280310u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2443466175u,        0.7244836743u)--(        0.2909382009u,        0.6558280310u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3704458338u,        0.6321023209u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3704458338u,        0.6321023209u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2443466175u,        0.7244836743u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2443466175u,        0.7244836743u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2909382009u,        0.6558280310u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2909382009u,        0.6558280310u) withcolor MidnightBlue;
draw (        0.3704458338u,        0.6321023209u)--(        0.3238542504u,        0.7007579642u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2443466175u,        0.7244836743u)--(        0.3238542504u,        0.7007579642u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3704458338u,        0.6321023209u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3704458338u,        0.6321023209u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2443466175u,        0.7244836743u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2443466175u,        0.7244836743u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3238542504u,        0.7007579642u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3238542504u,        0.7007579642u) withcolor MidnightBlue;
draw (        0.3974132415u,        0.7391899244u)--(        0.4444444444u,        0.8075499103u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3974132415u,        0.7391899244u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3974132415u,        0.7391899244u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4444444444u,        0.8075499103u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4444444444u,        0.8075499103u) withcolor Yellow;
draw (        0.3974132415u,        0.7391899244u)--(        0.3238630647u,        0.7007875222u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2443466175u,        0.7244836743u)--(        0.3238630647u,        0.7007875222u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3974132415u,        0.7391899244u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3974132415u,        0.7391899244u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2443466175u,        0.7244836743u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2443466175u,        0.7244836743u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3238630647u,        0.7007875222u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3238630647u,        0.7007875222u) withcolor MidnightBlue;
draw (        0.3974132415u,        0.7391899244u)--(        0.3178967943u,        0.7628860765u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2443466175u,        0.7244836743u)--(        0.3178967943u,        0.7628860765u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3974132415u,        0.7391899244u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3974132415u,        0.7391899244u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2443466175u,        0.7244836743u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2443466175u,        0.7244836743u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3178967943u,        0.7628860765u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3178967943u,        0.7628860765u) withcolor MidnightBlue;
draw (        0.3974132415u,        0.7391899244u)--(        0.3178985051u,        0.7628918165u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2884127802u,        0.8404480221u)--(        0.3178985051u,        0.7628918165u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3974132415u,        0.7391899244u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3974132415u,        0.7391899244u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2884127802u,        0.8404480221u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2884127802u,        0.8404480221u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3178985051u,        0.7628918165u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3178985051u,        0.7628918165u) withcolor MidnightBlue;
draw (        0.0555660359u,        0.0616416844u)--(        0.0000000000u,        0.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0555660359u,        0.0616416844u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0555660359u,        0.0616416844u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.0555660359u,        0.0616416844u)--(        0.1111111111u,        0.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0555660359u,        0.0616416844u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0555660359u,        0.0616416844u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1111111111u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1111111111u,        0.0000000000u) withcolor Yellow;
draw (        0.0555660359u,        0.0616416844u)--(        0.1378003962u,        0.0726817472u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1402378357u,        0.1556180595u)--(        0.1378003962u,        0.0726817472u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0555660359u,        0.0616416844u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0555660359u,        0.0616416844u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1402378357u,        0.1556180595u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1402378357u,        0.1556180595u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1378003962u,        0.0726817472u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1378003962u,        0.0726817472u) withcolor MidnightBlue;
draw (        0.1402378357u,        0.1556180595u)--(        0.1473127912u,        0.2382879936u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2142319893u,        0.2873409642u)--(        0.1473127912u,        0.2382879936u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1402378357u,        0.1556180595u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1402378357u,        0.1556180595u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2142319893u,        0.2873409642u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2142319893u,        0.2873409642u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1473127912u,        0.2382879936u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1473127912u,        0.2382879936u) withcolor MidnightBlue;
draw (        0.1402378357u,        0.1556180595u)--(        0.2071570338u,        0.2046710301u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2142319893u,        0.2873409642u)--(        0.2071570338u,        0.2046710301u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1402378357u,        0.1556180595u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1402378357u,        0.1556180595u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2142319893u,        0.2873409642u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2142319893u,        0.2873409642u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2071570338u,        0.2046710301u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2071570338u,        0.2046710301u) withcolor MidnightBlue;
draw (        0.1402378357u,        0.1556180595u)--(        0.1473122513u,        0.2382880398u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0696073572u,        0.2673796612u)--(        0.1473122513u,        0.2382880398u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1402378357u,        0.1556180595u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1402378357u,        0.1556180595u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0696073572u,        0.2673796612u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0696073572u,        0.2673796612u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1473122513u,        0.2382880398u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1473122513u,        0.2382880398u) withcolor MidnightBlue;
draw (        0.2142319893u,        0.2873409642u)--(        0.1924594297u,        0.3674054923u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2375129981u,        0.4370800858u)--(        0.1924594297u,        0.3674054923u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2142319893u,        0.2873409642u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2142319893u,        0.2873409642u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2375129981u,        0.4370800858u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2375129981u,        0.4370800858u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1924594297u,        0.3674054923u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1924594297u,        0.3674054923u) withcolor Green;
draw (        0.2142319893u,        0.2873409642u)--(        0.2592855578u,        0.3570155577u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2375129981u,        0.4370800858u)--(        0.2592855578u,        0.3570155577u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2142319893u,        0.2873409642u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2142319893u,        0.2873409642u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2375129981u,        0.4370800858u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2375129981u,        0.4370800858u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2592855578u,        0.3570155577u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2592855578u,        0.3570155577u) withcolor MidnightBlue;
draw (        0.2142319893u,        0.2873409642u)--(        0.1473124224u,        0.2382884968u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0696073572u,        0.2673796612u)--(        0.1473124224u,        0.2382884968u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2142319893u,        0.2873409642u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2142319893u,        0.2873409642u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0696073572u,        0.2673796612u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0696073572u,        0.2673796612u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1473124224u,        0.2382884968u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1473124224u,        0.2382884968u) withcolor MidnightBlue;
draw (        0.2142319893u,        0.2873409642u)--(        0.1613576454u,        0.3512837258u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0806788227u,        0.3319110498u)--(        0.1613576454u,        0.3512837258u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2142319893u,        0.2873409642u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2142319893u,        0.2873409642u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0806788227u,        0.3319110498u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0806788227u,        0.3319110498u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1613576454u,        0.3512837258u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1613576454u,        0.3512837258u) withcolor Green;
draw (        0.2375129981u,        0.4370800858u)--(        0.2726009812u,        0.5122678936u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2414982237u,        0.5891898766u)--(        0.2726009812u,        0.5122678936u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2375129981u,        0.4370800858u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2375129981u,        0.4370800858u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2414982237u,        0.5891898766u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2414982237u,        0.5891898766u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2726009812u,        0.5122678936u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2726009812u,        0.5122678936u) withcolor MidnightBlue;
draw (        0.2414982237u,        0.5891898766u)--(        0.1949097143u,        0.6578476060u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2443466175u,        0.7244836743u)--(        0.1949097143u,        0.6578476060u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2414982237u,        0.5891898766u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2414982237u,        0.5891898766u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2443466175u,        0.7244836743u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2443466175u,        0.7244836743u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1949097143u,        0.6578476060u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1949097143u,        0.6578476060u) withcolor MidnightBlue;
draw (        0.2414982237u,        0.5891898766u)--(        0.2909351269u,        0.6558259450u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2443466175u,        0.7244836743u)--(        0.2909351269u,        0.6558259450u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2414982237u,        0.5891898766u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2414982237u,        0.5891898766u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2443466175u,        0.7244836743u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2443466175u,        0.7244836743u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2909351269u,        0.6558259450u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2909351269u,        0.6558259450u) withcolor MidnightBlue;
draw (        0.2414982237u,        0.5891898766u)--(        0.1948961966u,        0.6578384314u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1149882993u,        0.6801789990u)--(        0.1948961966u,        0.6578384314u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2414982237u,        0.5891898766u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2414982237u,        0.5891898766u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1149882993u,        0.6801789990u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1149882993u,        0.6801789990u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1948961966u,        0.6578384314u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1948961966u,        0.6578384314u) withcolor MidnightBlue;
draw (        0.2443466175u,        0.7244836743u)--(        0.2148648208u,        0.8020413732u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2884127802u,        0.8404480221u)--(        0.2148648208u,        0.8020413732u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2443466175u,        0.7244836743u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2443466175u,        0.7244836743u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2884127802u,        0.8404480221u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2884127802u,        0.8404480221u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2148648208u,        0.8020413732u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2148648208u,        0.8020413732u) withcolor MidnightBlue;
draw (        0.2443466175u,        0.7244836743u)--(        0.3178945769u,        0.7628903232u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2884127802u,        0.8404480221u)--(        0.3178945769u,        0.7628903232u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2443466175u,        0.7244836743u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2443466175u,        0.7244836743u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2884127802u,        0.8404480221u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2884127802u,        0.8404480221u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3178945769u,        0.7628903232u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3178945769u,        0.7628903232u) withcolor MidnightBlue;
draw (        0.2443466175u,        0.7244836743u)--(        0.1949006419u,        0.6578543377u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1149882993u,        0.6801789990u)--(        0.1949006419u,        0.6578543377u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2443466175u,        0.7244836743u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2443466175u,        0.7244836743u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1149882993u,        0.6801789990u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1149882993u,        0.6801789990u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1949006419u,        0.6578543377u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1949006419u,        0.6578543377u) withcolor MidnightBlue;
draw (        0.2443466175u,        0.7244836743u)--(        0.1644342749u,        0.7468083356u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1149882993u,        0.6801789990u)--(        0.1644342749u,        0.7468083356u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2443466175u,        0.7244836743u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2443466175u,        0.7244836743u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1149882993u,        0.6801789990u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1149882993u,        0.6801789990u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1644342749u,        0.7468083356u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1644342749u,        0.7468083356u) withcolor MidnightBlue;
draw (        0.2443466175u,        0.7244836743u)--(        0.1644411660u,        0.7468329881u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1349475170u,        0.8243861806u)--(        0.1644411660u,        0.7468329881u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2443466175u,        0.7244836743u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2443466175u,        0.7244836743u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1349475170u,        0.8243861806u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1349475170u,        0.8243861806u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1644411660u,        0.7468329881u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1644411660u,        0.7468329881u) withcolor MidnightBlue;
draw (        0.2443466175u,        0.7244836743u)--(        0.2148529685u,        0.8020368668u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1349475170u,        0.8243861806u)--(        0.2148529685u,        0.8020368668u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2443466175u,        0.7244836743u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2443466175u,        0.7244836743u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1349475170u,        0.8243861806u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1349475170u,        0.8243861806u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2148529685u,        0.8020368668u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2148529685u,        0.8020368668u) withcolor MidnightBlue;
draw (        0.2884127802u,        0.8404480221u)--(        0.2148579224u,        0.8020545864u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1349475170u,        0.8243861806u)--(        0.2148579224u,        0.8020545864u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2884127802u,        0.8404480221u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2884127802u,        0.8404480221u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1349475170u,        0.8243861806u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1349475170u,        0.8243861806u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2148579224u,        0.8020545864u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2148579224u,        0.8020545864u) withcolor MidnightBlue;
draw (        0.2884127802u,        0.8404480221u)--(        0.2085023747u,        0.8627796164u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1349475170u,        0.8243861806u)--(        0.2085023747u,        0.8627796164u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2884127802u,        0.8404480221u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2884127802u,        0.8404480221u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1349475170u,        0.8243861806u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1349475170u,        0.8243861806u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2085023747u,        0.8627796164u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2085023747u,        0.8627796164u) withcolor MidnightBlue;
draw (        0.2884127802u,        0.8404480221u)--(        0.2085055865u,        0.8627911061u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2060965105u,        0.9457282471u)--(        0.2085055865u,        0.8627911061u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2884127802u,        0.8404480221u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2884127802u,        0.8404480221u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2060965105u,        0.9457282471u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2060965105u,        0.9457282471u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2085055865u,        0.8627911061u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2085055865u,        0.8627911061u) withcolor MidnightBlue;
draw (        0.0696073572u,        0.2673796612u)--(        0.0000000000u,        0.2222222222u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0696073572u,        0.2673796612u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0696073572u,        0.2673796612u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.2222222222u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.2222222222u) withcolor Yellow;
draw (        0.0696073572u,        0.2673796612u)--(        0.0000002339u,        0.3125374002u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0806788227u,        0.3319110498u)--(        0.0000002339u,        0.3125374002u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0696073572u,        0.2673796612u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0696073572u,        0.2673796612u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0806788227u,        0.3319110498u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0806788227u,        0.3319110498u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000002339u,        0.3125374002u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000002339u,        0.3125374002u) withcolor Green;
draw (        0.1149882993u,        0.6801789990u)--(        0.1644146372u,        0.7468229043u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1349475170u,        0.8243861806u)--(        0.1644146372u,        0.7468229043u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1149882993u,        0.6801789990u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1149882993u,        0.6801789990u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1349475170u,        0.8243861806u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1349475170u,        0.8243861806u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1644146372u,        0.7468229043u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1644146372u,        0.7468229043u) withcolor MidnightBlue;
draw (        0.1349475170u,        0.8243861806u)--(        0.1325475138u,        0.9073235847u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2060965105u,        0.9457282471u)--(        0.1325475138u,        0.9073235847u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1349475170u,        0.8243861806u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1349475170u,        0.8243861806u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2060965105u,        0.9457282471u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2060965105u,        0.9457282471u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1325475138u,        0.9073235847u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1325475138u,        0.9073235847u) withcolor MidnightBlue;
draw (        0.1349475170u,        0.8243861806u)--(        0.2084965136u,        0.8627908430u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2060965105u,        0.9457282471u)--(        0.2084965136u,        0.8627908430u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1349475170u,        0.8243861806u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1349475170u,        0.8243861806u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2060965105u,        0.9457282471u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2060965105u,        0.9457282471u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2084965136u,        0.8627908430u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2084965136u,        0.8627908430u) withcolor MidnightBlue;
draw (        0.1349475170u,        0.8243861806u)--(        0.1325220930u,        0.9073228452u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0555721452u,        0.9383563519u)--(        0.1325220930u,        0.9073228452u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1349475170u,        0.8243861806u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1349475170u,        0.8243861806u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0555721452u,        0.9383563519u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0555721452u,        0.9383563519u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1325220930u,        0.9073228452u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1325220930u,        0.9073228452u) withcolor MidnightBlue;
draw (        0.2060965105u,        0.9457282471u)--(        0.1325333386u,        0.9073507440u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0555721452u,        0.9383563519u)--(        0.1325333386u,        0.9073507440u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2060965105u,        0.9457282471u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2060965105u,        0.9457282471u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0555721452u,        0.9383563519u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0555721452u,        0.9383563519u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1325333386u,        0.9073507440u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1325333386u,        0.9073507440u) withcolor MidnightBlue;
draw (        0.0555721452u,        0.9383563519u)--(        0.1111111111u,        1.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0555721452u,        0.9383563519u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0555721452u,        0.9383563519u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1111111111u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1111111111u,        1.0000000000u) withcolor Yellow;
draw (        0.0555721452u,        0.9383563519u)--(        0.0000000000u,        1.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0555721452u,        0.9383563519u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0555721452u,        0.9383563519u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        1.0000000000u) withcolor Yellow;
% r =         0.0829721218
endfig;
end;
