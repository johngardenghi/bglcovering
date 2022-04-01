prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(70);
u := 20cm;
path pols[];
pair p[];
R := .0701584274u;
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

x1 = .9624239322u;
y1 = .0592474057u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .8949337732u;
y2 = .0634615700u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .9159231418u;
y3 = .1752043439u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .9608255455u;
y4 = .2773237940u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .7516790074u;
y5 = .0528604410u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .7957262199u;
y6 = .1379054632u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .8292647776u;
y7 = .2424842946u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .9121278713u;
y8 = .3279087609u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .6625664864u;
y9 = .0554509374u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = .6706025817u;
y10 = .1590622606u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;

x11 = .7126255460u;
y11 = .2425771671u;
draw fullcircle scaled D shifted z11 withcolor RoyalBlue withpen pencircle scaled 1.0;

x12 = .7857304858u;
y12 = .3501206032u;
draw fullcircle scaled D shifted z12 withcolor RoyalBlue withpen pencircle scaled 1.0;

x13 = .9616537462u;
y13 = .7214743040u;
draw fullcircle scaled D shifted z13 withcolor RoyalBlue withpen pencircle scaled 1.0;

x14 = .5811608303u;
y14 = .2536968300u;
draw fullcircle scaled D shifted z14 withcolor RoyalBlue withpen pencircle scaled 1.0;

x15 = .6636604619u;
y15 = .3405281145u;
draw fullcircle scaled D shifted z15 withcolor RoyalBlue withpen pencircle scaled 1.0;

x16 = .7399317444u;
y16 = .4355099169u;
draw fullcircle scaled D shifted z16 withcolor RoyalBlue withpen pencircle scaled 1.0;

x17 = .9073736550u;
y17 = .6675540618u;
draw fullcircle scaled D shifted z17 withcolor RoyalBlue withpen pencircle scaled 1.0;

x18 = .9657264479u;
y18 = .8576690571u;
draw fullcircle scaled D shifted z18 withcolor RoyalBlue withpen pencircle scaled 1.0;

x19 = .9427460627u;
y19 = .9594468408u;
draw fullcircle scaled D shifted z19 withcolor RoyalBlue withpen pencircle scaled 1.0;

x20 = .5395420993u;
y20 = .3319789412u;
draw fullcircle scaled D shifted z20 withcolor RoyalBlue withpen pencircle scaled 1.0;

x21 = .6103567289u;
y21 = .4379139369u;
draw fullcircle scaled D shifted z21 withcolor RoyalBlue withpen pencircle scaled 1.0;

x22 = .6857371258u;
y22 = .5330324251u;
draw fullcircle scaled D shifted z22 withcolor RoyalBlue withpen pencircle scaled 1.0;

x23 = .7505718912u;
y23 = .5954886892u;
draw fullcircle scaled D shifted z23 withcolor RoyalBlue withpen pencircle scaled 1.0;

x24 = .8010874372u;
y24 = .7062894847u;
draw fullcircle scaled D shifted z24 withcolor RoyalBlue withpen pencircle scaled 1.0;

x25 = .8769902523u;
y25 = .7955573154u;
draw fullcircle scaled D shifted z25 withcolor RoyalBlue withpen pencircle scaled 1.0;

x26 = .8445471863u;
y26 = .9110821019u;
draw fullcircle scaled D shifted z26 withcolor RoyalBlue withpen pencircle scaled 1.0;

x27 = .4382938293u;
y27 = .3070342841u;
draw fullcircle scaled D shifted z27 withcolor RoyalBlue withpen pencircle scaled 1.0;

x28 = .4854976856u;
y28 = .4276065471u;
draw fullcircle scaled D shifted z28 withcolor RoyalBlue withpen pencircle scaled 1.0;

x29 = .5581398423u;
y29 = .5332272012u;
draw fullcircle scaled D shifted z29 withcolor RoyalBlue withpen pencircle scaled 1.0;

x30 = .6254751065u;
y30 = .6323931453u;
draw fullcircle scaled D shifted z30 withcolor RoyalBlue withpen pencircle scaled 1.0;

x31 = .6962804319u;
y31 = .7089143801u;
draw fullcircle scaled D shifted z31 withcolor RoyalBlue withpen pencircle scaled 1.0;

x32 = .7583985616u;
y32 = .8238558723u;
draw fullcircle scaled D shifted z32 withcolor RoyalBlue withpen pencircle scaled 1.0;

x33 = .7330268173u;
y33 = .9459670691u;
draw fullcircle scaled D shifted z33 withcolor RoyalBlue withpen pencircle scaled 1.0;

x34 = .2800211366u;
y34 = .0290135471u;
draw fullcircle scaled D shifted z34 withcolor RoyalBlue withpen pencircle scaled 1.0;

x35 = .3344852680u;
y35 = .0976237023u;
draw fullcircle scaled D shifted z35 withcolor RoyalBlue withpen pencircle scaled 1.0;

x36 = .4074091382u;
y36 = .2373075126u;
draw fullcircle scaled D shifted z36 withcolor RoyalBlue withpen pencircle scaled 1.0;

x37 = .3726159544u;
y37 = .4056856317u;
draw fullcircle scaled D shifted z37 withcolor RoyalBlue withpen pencircle scaled 1.0;

x38 = .4334123363u;
y38 = .5232449880u;
draw fullcircle scaled D shifted z38 withcolor RoyalBlue withpen pencircle scaled 1.0;

x39 = .5083645275u;
y39 = .6283960881u;
draw fullcircle scaled D shifted z39 withcolor RoyalBlue withpen pencircle scaled 1.0;

x40 = .5703551299u;
y40 = .7389688965u;
draw fullcircle scaled D shifted z40 withcolor RoyalBlue withpen pencircle scaled 1.0;

x41 = .6543932654u;
y41 = .8188439505u;
draw fullcircle scaled D shifted z41 withcolor RoyalBlue withpen pencircle scaled 1.0;

x42 = .6734560538u;
y42 = .9314223362u;
draw fullcircle scaled D shifted z42 withcolor RoyalBlue withpen pencircle scaled 1.0;

x43 = .1428250833u;
y43 = .0395134359u;
draw fullcircle scaled D shifted z43 withcolor RoyalBlue withpen pencircle scaled 1.0;

x44 = .1996312069u;
y44 = .1169609445u;
draw fullcircle scaled D shifted z44 withcolor RoyalBlue withpen pencircle scaled 1.0;

x45 = .2833059342u;
y45 = .1927184758u;
draw fullcircle scaled D shifted z45 withcolor RoyalBlue withpen pencircle scaled 1.0;

x46 = .3120121119u;
y46 = .3030399763u;
draw fullcircle scaled D shifted z46 withcolor RoyalBlue withpen pencircle scaled 1.0;

x47 = .3171107118u;
y47 = .5010859598u;
draw fullcircle scaled D shifted z47 withcolor RoyalBlue withpen pencircle scaled 1.0;

x48 = .3793531751u;
y48 = .6034449238u;
draw fullcircle scaled D shifted z48 withcolor RoyalBlue withpen pencircle scaled 1.0;

x49 = .4431540309u;
y49 = .6931898212u;
draw fullcircle scaled D shifted z49 withcolor RoyalBlue withpen pencircle scaled 1.0;

x50 = .0423981266u;
y50 = .0559113870u;
draw fullcircle scaled D shifted z50 withcolor RoyalBlue withpen pencircle scaled 1.0;

x51 = .0795828341u;
y51 = .1573146636u;
draw fullcircle scaled D shifted z51 withcolor RoyalBlue withpen pencircle scaled 1.0;

x52 = .1553269137u;
y52 = .2249745269u;
draw fullcircle scaled D shifted z52 withcolor RoyalBlue withpen pencircle scaled 1.0;

x53 = .1999059948u;
y53 = .2974231983u;
draw fullcircle scaled D shifted z53 withcolor RoyalBlue withpen pencircle scaled 1.0;

x54 = .2495238225u;
y54 = .4079901991u;
draw fullcircle scaled D shifted z54 withcolor RoyalBlue withpen pencircle scaled 1.0;

x55 = .2501164158u;
y55 = .5923692347u;
draw fullcircle scaled D shifted z55 withcolor RoyalBlue withpen pencircle scaled 1.0;

x56 = .3089074226u;
y56 = .6946386864u;
draw fullcircle scaled D shifted z56 withcolor RoyalBlue withpen pencircle scaled 1.0;

x57 = .4040594714u;
y57 = .7597689986u;
draw fullcircle scaled D shifted z57 withcolor RoyalBlue withpen pencircle scaled 1.0;

x58 = .0343451404u;
y58 = .2755920259u;
draw fullcircle scaled D shifted z58 withcolor RoyalBlue withpen pencircle scaled 1.0;

x59 = .0884715503u;
y59 = .3416154546u;
draw fullcircle scaled D shifted z59 withcolor RoyalBlue withpen pencircle scaled 1.0;

x60 = .2001368279u;
y60 = .7020700066u;
draw fullcircle scaled D shifted z60 withcolor RoyalBlue withpen pencircle scaled 1.0;

x61 = .2825292671u;
y61 = .8079750725u;
draw fullcircle scaled D shifted z61 withcolor RoyalBlue withpen pencircle scaled 1.0;

x62 = .3348659826u;
y62 = .9019795847u;
draw fullcircle scaled D shifted z62 withcolor RoyalBlue withpen pencircle scaled 1.0;

x63 = .0881652236u;
y63 = .6578469932u;
draw fullcircle scaled D shifted z63 withcolor RoyalBlue withpen pencircle scaled 1.0;

x64 = .1552978314u;
y64 = .7744234362u;
draw fullcircle scaled D shifted z64 withcolor RoyalBlue withpen pencircle scaled 1.0;

x65 = .1995151346u;
y65 = .8836770498u;
draw fullcircle scaled D shifted z65 withcolor RoyalBlue withpen pencircle scaled 1.0;

x66 = .2802505534u;
y66 = .9715157229u;
draw fullcircle scaled D shifted z66 withcolor RoyalBlue withpen pencircle scaled 1.0;

x67 = .0342998121u;
y67 = .7245202987u;
draw fullcircle scaled D shifted z67 withcolor RoyalBlue withpen pencircle scaled 1.0;

x68 = .0797963749u;
y68 = .8428152589u;
draw fullcircle scaled D shifted z68 withcolor RoyalBlue withpen pencircle scaled 1.0;

x69 = .1430256671u;
y69 = .9609014283u;
draw fullcircle scaled D shifted z69 withcolor RoyalBlue withpen pencircle scaled 1.0;

x70 = .0423686090u;
y70 = .9440711029u;
draw fullcircle scaled D shifted z70 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.9624239322u,        0.0592474057u)--(        1.0000000000u,        0.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9624239322u,        0.0592474057u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9624239322u,        0.0592474057u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.9624239322u,        0.0592474057u)--(        0.9248478146u,        0.0000001245u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8949337732u,        0.0634615700u)--(        0.9248478146u,        0.0000001245u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9624239322u,        0.0592474057u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9624239322u,        0.0592474057u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8949337732u,        0.0634615700u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8949337732u,        0.0634615700u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9248478146u,        0.0000001245u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9248478146u,        0.0000001245u) withcolor Green;
draw (        0.9624239322u,        0.0592474057u)--(        0.9688180508u,        0.1291138510u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9159231418u,        0.1752043439u)--(        0.9688180508u,        0.1291138510u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9624239322u,        0.0592474057u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9624239322u,        0.0592474057u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9159231418u,        0.1752043439u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9159231418u,        0.1752043439u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9688180508u,        0.1291138510u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9688180508u,        0.1291138510u) withcolor Green;
draw (        0.8949337732u,        0.0634615700u)--(        0.8650199958u,        0.1269231399u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9159231418u,        0.1752043439u)--(        0.8650199958u,        0.1269231399u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8949337732u,        0.0634615700u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8949337732u,        0.0634615700u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9159231418u,        0.1752043439u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9159231418u,        0.1752043439u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8650199958u,        0.1269231399u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8650199958u,        0.1269231399u) withcolor MidnightBlue;
draw (        0.8949337732u,        0.0634615700u)--(        0.8650197316u,        0.1269230154u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7957262199u,        0.1379054632u)--(        0.8650197316u,        0.1269230154u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8949337732u,        0.0634615700u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8949337732u,        0.0634615700u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7957262199u,        0.1379054632u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7957262199u,        0.1379054632u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8650197316u,        0.1269230154u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8650197316u,        0.1269230154u) withcolor MidnightBlue;
draw (        0.9159231418u,        0.1752043439u)--(        0.8994175297u,        0.2433935638u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9608255455u,        0.2773237940u)--(        0.8994175297u,        0.2433935638u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9159231418u,        0.1752043439u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9159231418u,        0.1752043439u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9608255455u,        0.2773237940u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9608255455u,        0.2773237940u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8994175297u,        0.2433935638u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8994175297u,        0.2433935638u) withcolor MidnightBlue;
draw (        0.9159231418u,        0.1752043439u)--(        0.9773311577u,        0.2091345741u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9608255455u,        0.2773237940u)--(        0.9773311577u,        0.2091345741u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9159231418u,        0.1752043439u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9159231418u,        0.1752043439u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9608255455u,        0.2773237940u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9608255455u,        0.2773237940u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9773311577u,        0.2091345741u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9773311577u,        0.2091345741u) withcolor Green;
draw (        0.9159231418u,        0.1752043439u)--(        0.8650197864u,        0.1269233607u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7957262199u,        0.1379054632u)--(        0.8650197864u,        0.1269233607u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9159231418u,        0.1752043439u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9159231418u,        0.1752043439u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7957262199u,        0.1379054632u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7957262199u,        0.1379054632u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8650197864u,        0.1269233607u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8650197864u,        0.1269233607u) withcolor MidnightBlue;
draw (        0.9159231418u,        0.1752043439u)--(        0.8994173132u,        0.2433935114u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8292647776u,        0.2424842946u)--(        0.8994173132u,        0.2433935114u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9159231418u,        0.1752043439u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9159231418u,        0.1752043439u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8292647776u,        0.2424842946u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8292647776u,        0.2424842946u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8994173132u,        0.2433935114u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8994173132u,        0.2433935114u) withcolor MidnightBlue;
draw (        0.9608255455u,        0.2773237940u)--(        0.8994173073u,        0.2433939662u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8292647776u,        0.2424842946u)--(        0.8994173073u,        0.2433939662u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9608255455u,        0.2773237940u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9608255455u,        0.2773237940u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8292647776u,        0.2424842946u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8292647776u,        0.2424842946u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8994173073u,        0.2433939662u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8994173073u,        0.2433939662u) withcolor MidnightBlue;
draw (        0.9608255455u,        0.2773237940u)--(        0.9802365824u,        0.3447434987u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9121278713u,        0.3279087609u)--(        0.9802365824u,        0.3447434987u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9608255455u,        0.2773237940u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9608255455u,        0.2773237940u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9121278713u,        0.3279087609u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9121278713u,        0.3279087609u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9802365824u,        0.3447434987u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9802365824u,        0.3447434987u) withcolor Green;
draw (        0.7516790074u,        0.0528604410u)--(        0.7281730935u,        0.1189639730u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7957262199u,        0.1379054632u)--(        0.7281730935u,        0.1189639730u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7516790074u,        0.0528604410u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7516790074u,        0.0528604410u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7957262199u,        0.1379054632u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7957262199u,        0.1379054632u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7281730935u,        0.1189639730u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7281730935u,        0.1189639730u) withcolor MidnightBlue;
draw (        0.7516790074u,        0.0528604410u)--(        0.8192321338u,        0.0718019312u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7957262199u,        0.1379054632u)--(        0.8192321338u,        0.0718019312u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7516790074u,        0.0528604410u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7516790074u,        0.0528604410u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7957262199u,        0.1379054632u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7957262199u,        0.1379054632u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8192321338u,        0.0718019312u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8192321338u,        0.0718019312u) withcolor Green;
draw (        0.7516790074u,        0.0528604410u)--(        0.7055484573u,        0.0000004629u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6625664864u,        0.0554509374u)--(        0.7055484573u,        0.0000004629u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7516790074u,        0.0528604410u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7516790074u,        0.0528604410u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6625664864u,        0.0554509374u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6625664864u,        0.0554509374u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7055484573u,        0.0000004629u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7055484573u,        0.0000004629u) withcolor Green;
draw (        0.7516790074u,        0.0528604410u)--(        0.7281727291u,        0.1189638434u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6706025817u,        0.1590622606u)--(        0.7281727291u,        0.1189638434u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7516790074u,        0.0528604410u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7516790074u,        0.0528604410u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6706025817u,        0.1590622606u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6706025817u,        0.1590622606u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7281727291u,        0.1189638434u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7281727291u,        0.1189638434u) withcolor MidnightBlue;
draw (        0.7957262199u,        0.1379054632u)--(        0.7709143964u,        0.2035299887u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8292647776u,        0.2424842946u)--(        0.7709143964u,        0.2035299887u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7957262199u,        0.1379054632u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7957262199u,        0.1379054632u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8292647776u,        0.2424842946u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8292647776u,        0.2424842946u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7709143964u,        0.2035299887u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7709143964u,        0.2035299887u) withcolor MidnightBlue;
draw (        0.7957262199u,        0.1379054632u)--(        0.7281730148u,        0.1189642536u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6706025817u,        0.1590622606u)--(        0.7281730148u,        0.1189642536u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7957262199u,        0.1379054632u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7957262199u,        0.1379054632u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6706025817u,        0.1590622606u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6706025817u,        0.1590622606u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7281730148u,        0.1189642536u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7281730148u,        0.1189642536u) withcolor MidnightBlue;
draw (        0.7957262199u,        0.1379054632u)--(        0.7709136235u,        0.2035296965u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7126255460u,        0.2425771671u)--(        0.7709136235u,        0.2035296965u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7957262199u,        0.1379054632u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7957262199u,        0.1379054632u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7126255460u,        0.2425771671u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7126255460u,        0.2425771671u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7709136235u,        0.2035296965u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7709136235u,        0.2035296965u) withcolor MidnightBlue;
draw (        0.8292647776u,        0.2424842946u)--(        0.8440192319u,        0.3110737326u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9121278713u,        0.3279087609u)--(        0.8440192319u,        0.3110737326u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8292647776u,        0.2424842946u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8292647776u,        0.2424842946u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9121278713u,        0.3279087609u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9121278713u,        0.3279087609u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8440192319u,        0.3110737326u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8440192319u,        0.3110737326u) withcolor MidnightBlue;
draw (        0.8292647776u,        0.2424842946u)--(        0.7709141083u,        0.2035304202u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7126255460u,        0.2425771671u)--(        0.7709141083u,        0.2035304202u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8292647776u,        0.2424842946u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8292647776u,        0.2424842946u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7126255460u,        0.2425771671u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7126255460u,        0.2425771671u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7709141083u,        0.2035304202u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7709141083u,        0.2035304202u) withcolor MidnightBlue;
draw (        0.8292647776u,        0.2424842946u)--(        0.7709762153u,        0.2815310414u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7126255460u,        0.2425771671u)--(        0.7709762153u,        0.2815310414u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8292647776u,        0.2424842946u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8292647776u,        0.2424842946u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7126255460u,        0.2425771671u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7126255460u,        0.2425771671u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7709762153u,        0.2815310414u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7709762153u,        0.2815310414u) withcolor MidnightBlue;
draw (        0.8292647776u,        0.2424842946u)--(        0.7709762646u,        0.2815311150u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7857304858u,        0.3501206032u)--(        0.7709762646u,        0.2815311150u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8292647776u,        0.2424842946u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8292647776u,        0.2424842946u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7857304858u,        0.3501206032u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7857304858u,        0.3501206032u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7709762646u,        0.2815311150u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7709762646u,        0.2815311150u) withcolor MidnightBlue;
draw (        0.8292647776u,        0.2424842946u)--(        0.8440189988u,        0.3110737828u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7857304858u,        0.3501206032u)--(        0.8440189988u,        0.3110737828u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8292647776u,        0.2424842946u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8292647776u,        0.2424842946u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7857304858u,        0.3501206032u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7857304858u,        0.3501206032u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8440189988u,        0.3110737828u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8440189988u,        0.3110737828u) withcolor MidnightBlue;
draw (        0.9121278713u,        0.3279087609u)--(        0.8440191600u,        0.3110740234u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7857304858u,        0.3501206032u)--(        0.8440191600u,        0.3110740234u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9121278713u,        0.3279087609u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9121278713u,        0.3279087609u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7857304858u,        0.3501206032u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7857304858u,        0.3501206032u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8440191600u,        0.3110740234u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8440191600u,        0.3110740234u) withcolor MidnightBlue;
draw (        0.6625664864u,        0.0554509374u)--(        0.6195851126u,        0.1109018746u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6706025817u,        0.1590622606u)--(        0.6195851126u,        0.1109018746u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6625664864u,        0.0554509374u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6625664864u,        0.0554509374u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6706025817u,        0.1590622606u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6706025817u,        0.1590622606u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6195851126u,        0.1109018746u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6195851126u,        0.1109018746u) withcolor Green;
draw (        0.6706025817u,        0.1590622606u)--(        0.6448801617u,        0.2243352416u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7126255460u,        0.2425771671u)--(        0.6448801617u,        0.2243352416u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6706025817u,        0.1590622606u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6706025817u,        0.1590622606u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7126255460u,        0.2425771671u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7126255460u,        0.2425771671u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6448801617u,        0.2243352416u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6448801617u,        0.2243352416u) withcolor MidnightBlue;
draw (        0.6706025817u,        0.1590622606u)--(        0.6448796715u,        0.2243350484u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5811608303u,        0.2536968300u)--(        0.6448796715u,        0.2243350484u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6706025817u,        0.1590622606u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6706025817u,        0.1590622606u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5811608303u,        0.2536968300u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5811608303u,        0.2536968300u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6448796715u,        0.2243350484u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6448796715u,        0.2243350484u) withcolor MidnightBlue;
draw (        0.7126255460u,        0.2425771671u)--(        0.7273798821u,        0.3111666305u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7857304858u,        0.3501206032u)--(        0.7273798821u,        0.3111666305u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7126255460u,        0.2425771671u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7126255460u,        0.2425771671u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7857304858u,        0.3501206032u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7857304858u,        0.3501206032u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7273798821u,        0.3111666305u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7273798821u,        0.3111666305u) withcolor MidnightBlue;
draw (        0.7126255460u,        0.2425771671u)--(        0.7709761496u,        0.2815311397u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7857304858u,        0.3501206032u)--(        0.7709761496u,        0.2815311397u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7126255460u,        0.2425771671u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7126255460u,        0.2425771671u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7857304858u,        0.3501206032u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7857304858u,        0.3501206032u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7709761496u,        0.2815311397u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7709761496u,        0.2815311397u) withcolor MidnightBlue;
draw (        0.7126255460u,        0.2425771671u)--(        0.6448800137u,        0.2243357911u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5811608303u,        0.2536968300u)--(        0.6448800137u,        0.2243357911u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7126255460u,        0.2425771671u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7126255460u,        0.2425771671u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5811608303u,        0.2536968300u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5811608303u,        0.2536968300u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6448800137u,        0.2243357911u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6448800137u,        0.2243357911u) withcolor MidnightBlue;
draw (        0.7126255460u,        0.2425771671u)--(        0.6489063626u,        0.2719382059u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5811608303u,        0.2536968300u)--(        0.6489063626u,        0.2719382059u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7126255460u,        0.2425771671u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7126255460u,        0.2425771671u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5811608303u,        0.2536968300u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5811608303u,        0.2536968300u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6489063626u,        0.2719382059u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6489063626u,        0.2719382059u) withcolor MidnightBlue;
draw (        0.7126255460u,        0.2425771671u)--(        0.6489065278u,        0.2719385645u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6636604619u,        0.3405281145u)--(        0.6489065278u,        0.2719385645u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7126255460u,        0.2425771671u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7126255460u,        0.2425771671u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6636604619u,        0.3405281145u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6636604619u,        0.3405281145u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6489065278u,        0.2719385645u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6489065278u,        0.2719385645u) withcolor MidnightBlue;
draw (        0.7126255460u,        0.2425771671u)--(        0.7273794800u,        0.3111667170u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6636604619u,        0.3405281145u)--(        0.7273794800u,        0.3111667170u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7126255460u,        0.2425771671u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7126255460u,        0.2425771671u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6636604619u,        0.3405281145u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6636604619u,        0.3405281145u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7273794800u,        0.3111667170u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7273794800u,        0.3111667170u) withcolor MidnightBlue;
draw (        0.7857304858u,        0.3501206032u)--(        0.7273796207u,        0.3111670222u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6636604619u,        0.3405281145u)--(        0.7273796207u,        0.3111670222u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7857304858u,        0.3501206032u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7857304858u,        0.3501206032u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6636604619u,        0.3405281145u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6636604619u,        0.3405281145u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7273796207u,        0.3111670222u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7273796207u,        0.3111670222u) withcolor MidnightBlue;
draw (        0.7857304858u,        0.3501206032u)--(        0.8075493730u,        0.4167999935u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7399317444u,        0.4355099169u)--(        0.8075493730u,        0.4167999935u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7857304858u,        0.3501206032u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7857304858u,        0.3501206032u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7399317444u,        0.4355099169u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7399317444u,        0.4355099169u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8075493730u,        0.4167999935u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8075493730u,        0.4167999935u) withcolor Green;
draw (        0.9616537462u,        0.7214743040u)--(        0.9759610698u,        0.6527902053u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9073736550u,        0.6675540618u)--(        0.9759610698u,        0.6527902053u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9616537462u,        0.7214743040u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9616537462u,        0.7214743040u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9073736550u,        0.6675540618u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9073736550u,        0.6675540618u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9759610698u,        0.6527902053u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9759610698u,        0.6527902053u) withcolor Green;
draw (        0.9616537462u,        0.7214743040u)--(        0.9469403316u,        0.7900725574u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9657264479u,        0.8576690571u)--(        0.9469403316u,        0.7900725574u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9616537462u,        0.7214743040u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9616537462u,        0.7214743040u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9657264479u,        0.8576690571u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9657264479u,        0.8576690571u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9469403316u,        0.7900725574u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9469403316u,        0.7900725574u) withcolor MidnightBlue;
draw (        0.9616537462u,        0.7214743040u)--(        0.9804398625u,        0.7890708037u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9657264479u,        0.8576690571u)--(        0.9804398625u,        0.7890708037u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9616537462u,        0.7214743040u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9616537462u,        0.7214743040u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9657264479u,        0.8576690571u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9657264479u,        0.8576690571u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9804398625u,        0.7890708037u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9804398625u,        0.7890708037u) withcolor Green;
draw (        0.9616537462u,        0.7214743040u)--(        0.9469338518u,        0.7900711673u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8769902523u,        0.7955573154u)--(        0.9469338518u,        0.7900711673u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9616537462u,        0.7214743040u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9616537462u,        0.7214743040u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8769902523u,        0.7955573154u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8769902523u,        0.7955573154u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9469338518u,        0.7900711673u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9469338518u,        0.7900711673u) withcolor MidnightBlue;
draw (        0.5811608303u,        0.2536968300u)--(        0.6489062499u,        0.2719386243u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6636604619u,        0.3405281145u)--(        0.6489062499u,        0.2719386243u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5811608303u,        0.2536968300u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5811608303u,        0.2536968300u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6636604619u,        0.3405281145u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6636604619u,        0.3405281145u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6489062499u,        0.2719386243u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6489062499u,        0.2719386243u) withcolor MidnightBlue;
draw (        0.5811608303u,        0.2536968300u)--(        0.5123358586u,        0.2673103599u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5395420993u,        0.3319789412u)--(        0.5123358586u,        0.2673103599u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5811608303u,        0.2536968300u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5811608303u,        0.2536968300u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5395420993u,        0.3319789412u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5395420993u,        0.3319789412u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5123358586u,        0.2673103599u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5123358586u,        0.2673103599u) withcolor Green;
draw (        0.6636604619u,        0.3405281145u)--(        0.6746454848u,        0.4098212180u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7399317444u,        0.4355099169u)--(        0.6746454848u,        0.4098212180u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6636604619u,        0.3405281145u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6636604619u,        0.3405281145u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7399317444u,        0.4355099169u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7399317444u,        0.4355099169u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6746454848u,        0.4098212180u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6746454848u,        0.4098212180u) withcolor MidnightBlue;
draw (        0.6636604619u,        0.3405281145u)--(        0.5993718698u,        0.3686204954u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5395420993u,        0.3319789412u)--(        0.5993718698u,        0.3686204954u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6636604619u,        0.3405281145u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6636604619u,        0.3405281145u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5395420993u,        0.3319789412u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5395420993u,        0.3319789412u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5993718698u,        0.3686204954u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5993718698u,        0.3686204954u) withcolor MidnightBlue;
draw (        0.6636604619u,        0.3405281145u)--(        0.5993719973u,        0.3686207871u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6103567289u,        0.4379139369u)--(        0.5993719973u,        0.3686207871u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6636604619u,        0.3405281145u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6636604619u,        0.3405281145u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6103567289u,        0.4379139369u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6103567289u,        0.4379139369u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5993719973u,        0.3686207871u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5993719973u,        0.3686207871u) withcolor MidnightBlue;
draw (        0.6636604619u,        0.3405281145u)--(        0.6746451935u,        0.4098212642u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6103567289u,        0.4379139369u)--(        0.6746451935u,        0.4098212642u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6636604619u,        0.3405281145u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6636604619u,        0.3405281145u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6103567289u,        0.4379139369u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6103567289u,        0.4379139369u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6746451935u,        0.4098212642u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6746451935u,        0.4098212642u) withcolor MidnightBlue;
draw (        0.7399317444u,        0.4355099169u)--(        0.6746453372u,        0.4098215932u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6103567289u,        0.4379139369u)--(        0.6746453372u,        0.4098215932u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7399317444u,        0.4355099169u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7399317444u,        0.4355099169u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6103567289u,        0.4379139369u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6103567289u,        0.4379139369u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6746453372u,        0.4098215932u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6746453372u,        0.4098215932u) withcolor MidnightBlue;
draw (        0.7399317444u,        0.4355099169u)--(        0.6756431361u,        0.4636022606u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6103567289u,        0.4379139369u)--(        0.6756431361u,        0.4636022606u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7399317444u,        0.4355099169u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7399317444u,        0.4355099169u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6103567289u,        0.4379139369u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6103567289u,        0.4379139369u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6756431361u,        0.4636022606u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6756431361u,        0.4636022606u) withcolor MidnightBlue;
draw (        0.7399317444u,        0.4355099169u)--(        0.6756438209u,        0.4636038277u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6857371258u,        0.5330324251u)--(        0.6756438209u,        0.4636038277u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7399317444u,        0.4355099169u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7399317444u,        0.4355099169u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6857371258u,        0.5330324251u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6857371258u,        0.5330324251u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6756438209u,        0.4636038277u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6756438209u,        0.4636038277u) withcolor MidnightBlue;
draw (        0.9073736550u,        0.6675540618u)--(        0.8400176750u,        0.6479230430u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8010874372u,        0.7062894847u)--(        0.8400176750u,        0.6479230430u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9073736550u,        0.6675540618u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9073736550u,        0.6675540618u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8010874372u,        0.7062894847u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8010874372u,        0.7062894847u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8400176750u,        0.6479230430u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8400176750u,        0.6479230430u) withcolor Green;
draw (        0.9073736550u,        0.6675540618u)--(        0.8684434172u,        0.7259205035u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8010874372u,        0.7062894847u)--(        0.8684434172u,        0.7259205035u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9073736550u,        0.6675540618u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9073736550u,        0.6675540618u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8010874372u,        0.7062894847u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8010874372u,        0.7062894847u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8684434172u,        0.7259205035u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8684434172u,        0.7259205035u) withcolor MidnightBlue;
draw (        0.9073736550u,        0.6675540618u)--(        0.8684445903u,        0.7259212859u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8769902523u,        0.7955573154u)--(        0.8684445903u,        0.7259212859u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9073736550u,        0.6675540618u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9073736550u,        0.6675540618u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8769902523u,        0.7955573154u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8769902523u,        0.7955573154u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8684445903u,        0.7259212859u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8684445903u,        0.7259212859u) withcolor MidnightBlue;
draw (        0.9657264479u,        0.8576690571u)--(        0.9999938677u,        0.9188895510u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9427460627u,        0.9594468408u)--(        0.9999938677u,        0.9188895510u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9657264479u,        0.8576690571u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9657264479u,        0.8576690571u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9427460627u,        0.9594468408u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9427460627u,        0.9594468408u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9999938677u,        0.9188895510u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9999938677u,        0.9188895510u) withcolor Green;
draw (        0.9657264479u,        0.8576690571u)--(        0.9469340967u,        0.7900742905u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8769902523u,        0.7955573154u)--(        0.9469340967u,        0.7900742905u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9657264479u,        0.8576690571u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9657264479u,        0.8576690571u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8769902523u,        0.7955573154u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8769902523u,        0.7955573154u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9469340967u,        0.7900742905u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9469340967u,        0.7900742905u) withcolor MidnightBlue;
draw (        0.9657264479u,        0.8576690571u)--(        0.8957826035u,        0.8631520820u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8769902523u,        0.7955573154u)--(        0.8957826035u,        0.8631520820u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9657264479u,        0.8576690571u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9657264479u,        0.8576690571u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8769902523u,        0.7955573154u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8769902523u,        0.7955573154u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8957826035u,        0.8631520820u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8957826035u,        0.8631520820u) withcolor MidnightBlue;
draw (        0.9657264479u,        0.8576690571u)--(        0.8957827373u,        0.8631537882u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8445471863u,        0.9110821019u)--(        0.8957827373u,        0.8631537882u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9657264479u,        0.8576690571u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9657264479u,        0.8576690571u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8445471863u,        0.9110821019u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8445471863u,        0.9110821019u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8957827373u,        0.8631537882u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8957827373u,        0.8631537882u) withcolor MidnightBlue;
draw (        0.9427460627u,        0.9594468408u)--(        1.0000000000u,        1.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9427460627u,        0.9594468408u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9427460627u,        0.9594468408u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        1.0000000000u) withcolor Yellow;
draw (        0.9427460627u,        0.9594468408u)--(        0.8742527738u,        0.9746413892u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8445471863u,        0.9110821019u)--(        0.8742527738u,        0.9746413892u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9427460627u,        0.9594468408u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9427460627u,        0.9594468408u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8445471863u,        0.9110821019u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8445471863u,        0.9110821019u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8742527738u,        0.9746413892u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8742527738u,        0.9746413892u) withcolor Green;
draw (        0.5395420993u,        0.3319789412u)--(        0.5505271700u,        0.4012720372u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6103567289u,        0.4379139369u)--(        0.5505271700u,        0.4012720372u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5395420993u,        0.3319789412u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5395420993u,        0.3319789412u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6103567289u,        0.4379139369u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6103567289u,        0.4379139369u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5505271700u,        0.4012720372u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5505271700u,        0.4012720372u) withcolor MidnightBlue;
draw (        0.5395420993u,        0.3319789412u)--(        0.5993716583u,        0.3686208408u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6103567289u,        0.4379139369u)--(        0.5993716583u,        0.3686208408u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5395420993u,        0.3319789412u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5395420993u,        0.3319789412u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6103567289u,        0.4379139369u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6103567289u,        0.4379139369u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5993716583u,        0.3686208408u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5993716583u,        0.3686208408u) withcolor MidnightBlue;
draw (        0.5395420993u,        0.3319789412u)--(        0.5505261923u,        0.4012721922u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4854976856u,        0.4276065471u)--(        0.5505261923u,        0.4012721922u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5395420993u,        0.3319789412u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5395420993u,        0.3319789412u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4854976856u,        0.4276065471u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4854976856u,        0.4276065471u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5505261923u,        0.4012721922u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5505261923u,        0.4012721922u) withcolor MidnightBlue;
draw (        0.6103567289u,        0.4379139369u)--(        0.6756424404u,        0.4636040284u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6857371258u,        0.5330324251u)--(        0.6756424404u,        0.4636040284u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6103567289u,        0.4379139369u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6103567289u,        0.4379139369u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6857371258u,        0.5330324251u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6857371258u,        0.5330324251u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6756424404u,        0.4636040284u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6756424404u,        0.4636040284u) withcolor MidnightBlue;
draw (        0.6103567289u,        0.4379139369u)--(        0.5505265437u,        0.4012730599u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4854976856u,        0.4276065471u)--(        0.5505265437u,        0.4012730599u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6103567289u,        0.4379139369u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6103567289u,        0.4379139369u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4854976856u,        0.4276065471u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4854976856u,        0.4276065471u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5505265437u,        0.4012730599u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5505265437u,        0.4012730599u) withcolor MidnightBlue;
draw (        0.6103567289u,        0.4379139369u)--(        0.5453278708u,        0.4642474241u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4854976856u,        0.4276065471u)--(        0.5453278708u,        0.4642474241u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6103567289u,        0.4379139369u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6103567289u,        0.4379139369u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4854976856u,        0.4276065471u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4854976856u,        0.4276065471u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5453278708u,        0.4642474241u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5453278708u,        0.4642474241u) withcolor MidnightBlue;
draw (        0.6103567289u,        0.4379139369u)--(        0.5453282840u,        0.4642484443u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5581398423u,        0.5332272012u)--(        0.5453282840u,        0.4642484443u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6103567289u,        0.4379139369u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6103567289u,        0.4379139369u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5581398423u,        0.5332272012u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5581398423u,        0.5332272012u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5453282840u,        0.4642484443u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5453282840u,        0.4642484443u) withcolor MidnightBlue;
draw (        0.6857371258u,        0.5330324251u)--(        0.7554903397u,        0.5255028778u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7505718912u,        0.5954886892u)--(        0.7554903397u,        0.5255028778u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6857371258u,        0.5330324251u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6857371258u,        0.5330324251u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7505718912u,        0.5954886892u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7505718912u,        0.5954886892u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7554903397u,        0.5255028778u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7554903397u,        0.5255028778u) withcolor Green;
draw (        0.6857371258u,        0.5330324251u)--(        0.6219830390u,        0.5623175979u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5581398423u,        0.5332272012u)--(        0.6219830390u,        0.5623175979u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6857371258u,        0.5330324251u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6857371258u,        0.5330324251u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5581398423u,        0.5332272012u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5581398423u,        0.5332272012u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6219830390u,        0.5623175979u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6219830390u,        0.5623175979u) withcolor MidnightBlue;
draw (        0.6857371258u,        0.5330324251u)--(        0.6219848718u,        0.5623215876u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6254751065u,        0.6323931453u)--(        0.6219848718u,        0.5623215876u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6857371258u,        0.5330324251u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6857371258u,        0.5330324251u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6254751065u,        0.6323931453u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6254751065u,        0.6323931453u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6219848718u,        0.5623215876u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6219848718u,        0.5623215876u) withcolor MidnightBlue;
draw (        0.7505718912u,        0.5954886892u)--(        0.8075464073u,        0.6364290049u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8010874372u,        0.7062894847u)--(        0.8075464073u,        0.6364290049u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7505718912u,        0.5954886892u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7505718912u,        0.5954886892u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8010874372u,        0.7062894847u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8010874372u,        0.7062894847u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8075464073u,        0.6364290049u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8075464073u,        0.6364290049u) withcolor Green;
draw (        0.7505718912u,        0.5954886892u)--(        0.6953443404u,        0.6387567209u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6254751065u,        0.6323931453u)--(        0.6953443404u,        0.6387567209u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7505718912u,        0.5954886892u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7505718912u,        0.5954886892u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6254751065u,        0.6323931453u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6254751065u,        0.6323931453u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6953443404u,        0.6387567209u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6953443404u,        0.6387567209u) withcolor MidnightBlue;
draw (        0.7505718912u,        0.5954886892u)--(        0.6953485875u,        0.6387621413u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6962804319u,        0.7089143801u)--(        0.6953485875u,        0.6387621413u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7505718912u,        0.5954886892u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7505718912u,        0.5954886892u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6962804319u,        0.7089143801u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6962804319u,        0.7089143801u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6953485875u,        0.6387621413u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6953485875u,        0.6387621413u) withcolor MidnightBlue;
draw (        0.8010874372u,        0.7062894847u)--(        0.8096345524u,        0.7759253358u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8769902523u,        0.7955573154u)--(        0.8096345524u,        0.7759253358u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8010874372u,        0.7062894847u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8010874372u,        0.7062894847u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8769902523u,        0.7955573154u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8769902523u,        0.7955573154u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8096345524u,        0.7759253358u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8096345524u,        0.7759253358u) withcolor MidnightBlue;
draw (        0.8010874372u,        0.7062894847u)--(        0.8684431371u,        0.7259214643u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8769902523u,        0.7955573154u)--(        0.8684431371u,        0.7259214643u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8010874372u,        0.7062894847u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8010874372u,        0.7062894847u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8769902523u,        0.7955573154u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8769902523u,        0.7955573154u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8684431371u,        0.7259214643u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8684431371u,        0.7259214643u) withcolor MidnightBlue;
draw (        0.8010874372u,        0.7062894847u)--(        0.7498514182u,        0.7542172980u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6962804319u,        0.7089143801u)--(        0.7498514182u,        0.7542172980u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8010874372u,        0.7062894847u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8010874372u,        0.7062894847u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6962804319u,        0.7089143801u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6962804319u,        0.7089143801u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7498514182u,        0.7542172980u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7498514182u,        0.7542172980u) withcolor MidnightBlue;
draw (        0.8010874372u,        0.7062894847u)--(        0.7498537062u,        0.7542197438u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7583985616u,        0.8238558723u)--(        0.7498537062u,        0.7542197438u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8010874372u,        0.7062894847u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8010874372u,        0.7062894847u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7583985616u,        0.8238558723u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7583985616u,        0.8238558723u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7498537062u,        0.7542197438u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7498537062u,        0.7542197438u) withcolor MidnightBlue;
draw (        0.8010874372u,        0.7062894847u)--(        0.8096322926u,        0.7759256131u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7583985616u,        0.8238558723u)--(        0.8096322926u,        0.7759256131u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8010874372u,        0.7062894847u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8010874372u,        0.7062894847u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7583985616u,        0.8238558723u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7583985616u,        0.8238558723u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8096322926u,        0.7759256131u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8096322926u,        0.7759256131u) withcolor MidnightBlue;
draw (        0.8769902523u,        0.7955573154u)--(        0.8257559957u,        0.8434870126u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8445471863u,        0.9110821019u)--(        0.8257559957u,        0.8434870126u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8769902523u,        0.7955573154u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8769902523u,        0.7955573154u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8445471863u,        0.9110821019u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8445471863u,        0.9110821019u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8257559957u,        0.8434870126u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8257559957u,        0.8434870126u) withcolor MidnightBlue;
draw (        0.8769902523u,        0.7955573154u)--(        0.8957814430u,        0.8631524046u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8445471863u,        0.9110821019u)--(        0.8957814430u,        0.8631524046u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8769902523u,        0.7955573154u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8769902523u,        0.7955573154u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8445471863u,        0.9110821019u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8445471863u,        0.9110821019u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8957814430u,        0.8631524046u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8957814430u,        0.8631524046u) withcolor MidnightBlue;
draw (        0.8769902523u,        0.7955573154u)--(        0.8096339540u,        0.7759273890u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7583985616u,        0.8238558723u)--(        0.8096339540u,        0.7759273890u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8769902523u,        0.7955573154u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8769902523u,        0.7955573154u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7583985616u,        0.8238558723u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7583985616u,        0.8238558723u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8096339540u,        0.7759273890u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8096339540u,        0.7759273890u) withcolor MidnightBlue;
draw (        0.8769902523u,        0.7955573154u)--(        0.8257548599u,        0.8434857986u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7583985616u,        0.8238558723u)--(        0.8257548599u,        0.8434857986u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8769902523u,        0.7955573154u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8769902523u,        0.7955573154u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7583985616u,        0.8238558723u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7583985616u,        0.8238558723u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8257548599u,        0.8434857986u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8257548599u,        0.8434857986u) withcolor MidnightBlue;
draw (        0.8445471863u,        0.9110821019u)--(        0.8257543748u,        0.8434874633u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7583985616u,        0.8238558723u)--(        0.8257543748u,        0.8434874633u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8445471863u,        0.9110821019u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8445471863u,        0.9110821019u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7583985616u,        0.8238558723u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7583985616u,        0.8238558723u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8257543748u,        0.8434874633u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8257543748u,        0.8434874633u) withcolor MidnightBlue;
draw (        0.8445471863u,        0.9110821019u)--(        0.7771913731u,        0.8914505109u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7583985616u,        0.8238558723u)--(        0.7771913731u,        0.8914505109u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8445471863u,        0.9110821019u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8445471863u,        0.9110821019u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7583985616u,        0.8238558723u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7583985616u,        0.8238558723u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7771913731u,        0.8914505109u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7771913731u,        0.8914505109u) withcolor MidnightBlue;
draw (        0.8445471863u,        0.9110821019u)--(        0.7771905986u,        0.8914531685u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7330268173u,        0.9459670691u)--(        0.7771905986u,        0.8914531685u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8445471863u,        0.9110821019u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8445471863u,        0.9110821019u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7330268173u,        0.9459670691u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7330268173u,        0.9459670691u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7771905986u,        0.8914531685u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7771905986u,        0.8914531685u) withcolor MidnightBlue;
draw (        0.4382938293u,        0.3070342841u)--(        0.4367239363u,        0.3771751450u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4854976856u,        0.4276065471u)--(        0.4367239363u,        0.3771751450u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4382938293u,        0.3070342841u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4382938293u,        0.3070342841u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4854976856u,        0.4276065471u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4854976856u,        0.4276065471u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4367239363u,        0.3771751450u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4367239363u,        0.3771751450u) withcolor MidnightBlue;
draw (        0.4382938293u,        0.3070342841u)--(        0.4766978412u,        0.2483202585u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4074091382u,        0.2373075126u)--(        0.4766978412u,        0.2483202585u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4382938293u,        0.3070342841u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4382938293u,        0.3070342841u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4074091382u,        0.2373075126u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4074091382u,        0.2373075126u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4766978412u,        0.2483202585u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4766978412u,        0.2483202585u) withcolor Green;
draw (        0.4382938293u,        0.3070342841u)--(        0.3741895994u,        0.3355448549u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3726159544u,        0.4056856317u)--(        0.3741895994u,        0.3355448549u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4382938293u,        0.3070342841u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4382938293u,        0.3070342841u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3726159544u,        0.4056856317u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3726159544u,        0.4056856317u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3741895994u,        0.3355448549u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3741895994u,        0.3355448549u) withcolor MidnightBlue;
draw (        0.4382938293u,        0.3070342841u)--(        0.4367201843u,        0.3771750609u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3726159544u,        0.4056856317u)--(        0.4367201843u,        0.3771750609u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4382938293u,        0.3070342841u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4382938293u,        0.3070342841u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3726159544u,        0.4056856317u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3726159544u,        0.4056856317u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4367201843u,        0.3771750609u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4367201843u,        0.3771750609u) withcolor MidnightBlue;
draw (        0.4382938293u,        0.3070342841u)--(        0.3741881129u,        0.3355415122u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3120121119u,        0.3030399763u)--(        0.3741881129u,        0.3355415122u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4382938293u,        0.3070342841u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4382938293u,        0.3070342841u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3120121119u,        0.3030399763u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3120121119u,        0.3030399763u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3741881129u,        0.3355415122u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3741881129u,        0.3355415122u) withcolor MidnightBlue;
draw (        0.4854976856u,        0.4276065471u)--(        0.4983104151u,        0.4965850865u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5581398423u,        0.5332272012u)--(        0.4983104151u,        0.4965850865u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4854976856u,        0.4276065471u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4854976856u,        0.4276065471u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5581398423u,        0.5332272012u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5581398423u,        0.5332272012u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4983104151u,        0.4965850865u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4983104151u,        0.4965850865u) withcolor MidnightBlue;
draw (        0.4854976856u,        0.4276065471u)--(        0.5453271129u,        0.4642486618u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5581398423u,        0.5332272012u)--(        0.5453271129u,        0.4642486618u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4854976856u,        0.4276065471u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4854976856u,        0.4276065471u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5581398423u,        0.5332272012u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5581398423u,        0.5332272012u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5453271129u,        0.4642486618u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5453271129u,        0.4642486618u) withcolor MidnightBlue;
draw (        0.4854976856u,        0.4276065471u)--(        0.4367213389u,        0.3771776572u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3726159544u,        0.4056856317u)--(        0.4367213389u,        0.3771776572u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4854976856u,        0.4276065471u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4854976856u,        0.4276065471u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3726159544u,        0.4056856317u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3726159544u,        0.4056856317u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4367213389u,        0.3771776572u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4367213389u,        0.3771776572u) withcolor MidnightBlue;
draw (        0.4854976856u,        0.4276065471u)--(        0.4983082243u,        0.4965854934u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4334123363u,        0.5232449880u)--(        0.4983082243u,        0.4965854934u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4854976856u,        0.4276065471u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4854976856u,        0.4276065471u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4334123363u,        0.5232449880u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4334123363u,        0.5232449880u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4983082243u,        0.4965854934u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4983082243u,        0.4965854934u) withcolor MidnightBlue;
draw (        0.5581398423u,        0.5332272012u)--(        0.6219811360u,        0.5623217738u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6254751065u,        0.6323931453u)--(        0.6219811360u,        0.5623217738u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5581398423u,        0.5332272012u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5581398423u,        0.5332272012u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6254751065u,        0.6323931453u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6254751065u,        0.6323931453u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6219811360u,        0.5623217738u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6219811360u,        0.5623217738u) withcolor MidnightBlue;
draw (        0.5581398423u,        0.5332272012u)--(        0.4983090038u,        0.4965873910u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4334123363u,        0.5232449880u)--(        0.4983090038u,        0.4965873910u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5581398423u,        0.5332272012u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5581398423u,        0.5332272012u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4334123363u,        0.5232449880u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4334123363u,        0.5232449880u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4983090038u,        0.4965873910u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4983090038u,        0.4965873910u) withcolor MidnightBlue;
draw (        0.5581398423u,        0.5332272012u)--(        0.4932431748u,        0.5598847982u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4334123363u,        0.5232449880u)--(        0.4932431748u,        0.5598847982u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5581398423u,        0.5332272012u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5581398423u,        0.5332272012u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4334123363u,        0.5232449880u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4334123363u,        0.5232449880u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4932431748u,        0.5598847982u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4932431748u,        0.5598847982u) withcolor MidnightBlue;
draw (        0.5581398423u,        0.5332272012u)--(        0.4932438544u,        0.5598864524u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5083645275u,        0.6283960881u)--(        0.4932438544u,        0.5598864524u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5581398423u,        0.5332272012u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5581398423u,        0.5332272012u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5083645275u,        0.6283960881u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5083645275u,        0.6283960881u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4932438544u,        0.5598864524u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4932438544u,        0.5598864524u) withcolor MidnightBlue;
draw (        0.6254751065u,        0.6323931453u)--(        0.6953438407u,        0.6387622045u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6962804319u,        0.7089143801u)--(        0.6953438407u,        0.6387622045u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6254751065u,        0.6323931453u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6254751065u,        0.6323931453u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6962804319u,        0.7089143801u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6962804319u,        0.7089143801u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6953438407u,        0.6387622045u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6953438407u,        0.6387622045u) withcolor MidnightBlue;
draw (        0.6254751065u,        0.6323931453u)--(        0.5656033498u,        0.6689660544u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5083645275u,        0.6283960881u)--(        0.5656033498u,        0.6689660544u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6254751065u,        0.6323931453u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6254751065u,        0.6323931453u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5083645275u,        0.6283960881u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5083645275u,        0.6283960881u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5656033498u,        0.6689660544u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5656033498u,        0.6689660544u) withcolor MidnightBlue;
draw (        0.6254751065u,        0.6323931453u)--(        0.5656065860u,        0.6689713518u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5703551299u,        0.7389688965u)--(        0.5656065860u,        0.6689713518u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6254751065u,        0.6323931453u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6254751065u,        0.6323931453u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5703551299u,        0.7389688965u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5703551299u,        0.7389688965u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5656065860u,        0.6689713518u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5656065860u,        0.6689713518u) withcolor MidnightBlue;
draw (        0.6962804319u,        0.7089143801u)--(        0.7498488450u,        0.7542203405u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7583985616u,        0.8238558723u)--(        0.7498488450u,        0.7542203405u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6962804319u,        0.7089143801u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6962804319u,        0.7089143801u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7583985616u,        0.8238558723u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7583985616u,        0.8238558723u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7498488450u,        0.7542203405u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7498488450u,        0.7542203405u) withcolor MidnightBlue;
draw (        0.6962804319u,        0.7089143801u)--(        0.6395990786u,        0.7502596230u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5703551299u,        0.7389688965u)--(        0.6395990786u,        0.7502596230u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6962804319u,        0.7089143801u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6962804319u,        0.7089143801u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5703551299u,        0.7389688965u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5703551299u,        0.7389688965u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6395990786u,        0.7502596230u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6395990786u,        0.7502596230u) withcolor MidnightBlue;
draw (        0.6962804319u,        0.7089143801u)--(        0.6396012521u,        0.7502626026u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6543932654u,        0.8188439505u)--(        0.6396012521u,        0.7502626026u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6962804319u,        0.7089143801u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6962804319u,        0.7089143801u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6543932654u,        0.8188439505u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6543932654u,        0.8188439505u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6396012521u,        0.7502626026u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6396012521u,        0.7502626026u) withcolor MidnightBlue;
draw (        0.7583985616u,        0.8238558723u)--(        0.7771883540u,        0.8914513502u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7330268173u,        0.9459670691u)--(        0.7771883540u,        0.8914513502u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7583985616u,        0.8238558723u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7583985616u,        0.8238558723u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7330268173u,        0.9459670691u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7330268173u,        0.9459670691u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7771883540u,        0.8914513502u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7771883540u,        0.8914513502u) withcolor MidnightBlue;
draw (        0.7583985616u,        0.8238558723u)--(        0.7041322945u,        0.8683235818u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6543932654u,        0.8188439505u)--(        0.7041322945u,        0.8683235818u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7583985616u,        0.8238558723u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7583985616u,        0.8238558723u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6543932654u,        0.8188439505u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6543932654u,        0.8188439505u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7041322945u,        0.8683235818u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7041322945u,        0.8683235818u) withcolor MidnightBlue;
draw (        0.7583985616u,        0.8238558723u)--(        0.7041353130u,        0.8683272652u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6734560538u,        0.9314223362u)--(        0.7041353130u,        0.8683272652u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7583985616u,        0.8238558723u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7583985616u,        0.8238558723u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6734560538u,        0.9314223362u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6734560538u,        0.9314223362u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7041353130u,        0.8683272652u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7041353130u,        0.8683272652u) withcolor MidnightBlue;
draw (        0.7330268173u,        0.9459670691u)--(        0.7777777778u,        1.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7330268173u,        0.9459670691u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7330268173u,        0.9459670691u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7777777778u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7777777778u,        1.0000000000u) withcolor Yellow;
draw (        0.7330268173u,        0.9459670691u)--(        0.6882736497u,        0.9999981613u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6734560538u,        0.9314223362u)--(        0.6882736497u,        0.9999981613u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7330268173u,        0.9459670691u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7330268173u,        0.9459670691u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6734560538u,        0.9314223362u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6734560538u,        0.9314223362u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6882736497u,        0.9999981613u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6882736497u,        0.9999981613u) withcolor Green;
draw (        0.2800211366u,        0.0290135471u)--(        0.2124720549u,        0.0479694565u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1428250833u,        0.0395134359u)--(        0.2124720549u,        0.0479694565u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2800211366u,        0.0290135471u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2800211366u,        0.0290135471u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1428250833u,        0.0395134359u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1428250833u,        0.0395134359u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2124720549u,        0.0479694565u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2124720549u,        0.0479694565u) withcolor MidnightBlue;
draw (        0.2800211366u,        0.0290135471u)--(        0.2124774403u,        0.0479886368u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1996312069u,        0.1169609445u)--(        0.2124774403u,        0.0479886368u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2800211366u,        0.0290135471u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2800211366u,        0.0290135471u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1996312069u,        0.1169609445u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1996312069u,        0.1169609445u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2124774403u,        0.0479886368u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2124774403u,        0.0479886368u) withcolor MidnightBlue;
draw (        0.3344852680u,        0.0976237023u)--(        0.2694431932u,        0.1239245281u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1996312069u,        0.1169609445u)--(        0.2694431932u,        0.1239245281u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3344852680u,        0.0976237023u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3344852680u,        0.0976237023u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1996312069u,        0.1169609445u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1996312069u,        0.1169609445u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2694431932u,        0.1239245281u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2694431932u,        0.1239245281u) withcolor MidnightBlue;
draw (        0.3344852680u,        0.0976237023u)--(        0.2694502021u,        0.1239418546u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2833059342u,        0.1927184758u)--(        0.2694502021u,        0.1239418546u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3344852680u,        0.0976237023u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3344852680u,        0.0976237023u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2833059342u,        0.1927184758u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2833059342u,        0.1927184758u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2694502021u,        0.1239418546u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2694502021u,        0.1239418546u) withcolor MidnightBlue;
draw (        0.4074091382u,        0.2373075126u)--(        0.3534638516u,        0.1924509525u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2833059342u,        0.1927184758u)--(        0.3534638516u,        0.1924509525u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4074091382u,        0.2373075126u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4074091382u,        0.2373075126u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2833059342u,        0.1927184758u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2833059342u,        0.1927184758u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3534638516u,        0.1924509525u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3534638516u,        0.1924509525u) withcolor Green;
draw (        0.4074091382u,        0.2373075126u)--(        0.3372512208u,        0.2375750359u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2833059342u,        0.1927184758u)--(        0.3372512208u,        0.2375750359u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4074091382u,        0.2373075126u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4074091382u,        0.2373075126u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2833059342u,        0.1927184758u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2833059342u,        0.1927184758u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3372512208u,        0.2375750359u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3372512208u,        0.2375750359u) withcolor MidnightBlue;
draw (        0.4074091382u,        0.2373075126u)--(        0.3372512345u,        0.2375786013u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3120121119u,        0.3030399763u)--(        0.3372512345u,        0.2375786013u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4074091382u,        0.2373075126u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4074091382u,        0.2373075126u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3120121119u,        0.3030399763u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3120121119u,        0.3030399763u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3372512345u,        0.2375786013u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3372512345u,        0.2375786013u) withcolor MidnightBlue;
draw (        0.3726159544u,        0.4056856317u)--(        0.3823138573u,        0.4751705626u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4334123363u,        0.5232449880u)--(        0.3823138573u,        0.4751705626u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3726159544u,        0.4056856317u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3726159544u,        0.4056856317u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4334123363u,        0.5232449880u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4334123363u,        0.5232449880u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3823138573u,        0.4751705626u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3823138573u,        0.4751705626u) withcolor MidnightBlue;
draw (        0.3726159544u,        0.4056856317u)--(        0.3741864029u,        0.3355447833u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3120121119u,        0.3030399763u)--(        0.3741864029u,        0.3355447833u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3726159544u,        0.4056856317u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3726159544u,        0.4056856317u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3120121119u,        0.3030399763u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3120121119u,        0.3030399763u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3741864029u,        0.3355447833u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3741864029u,        0.3355447833u) withcolor MidnightBlue;
draw (        0.3726159544u,        0.4056856317u)--(        0.3104416634u,        0.3731808248u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3120121119u,        0.3030399763u)--(        0.3104416634u,        0.3731808248u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3726159544u,        0.4056856317u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3726159544u,        0.4056856317u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3120121119u,        0.3030399763u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3120121119u,        0.3030399763u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3104416634u,        0.3731808248u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3104416634u,        0.3731808248u) withcolor MidnightBlue;
draw (        0.3726159544u,        0.4056856317u)--(        0.3823076591u,        0.4751714274u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3171107118u,        0.5010859598u)--(        0.3823076591u,        0.4751714274u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3726159544u,        0.4056856317u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3726159544u,        0.4056856317u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3171107118u,        0.5010859598u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3171107118u,        0.5010859598u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3823076591u,        0.4751714274u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3823076591u,        0.4751714274u) withcolor MidnightBlue;
draw (        0.3726159544u,        0.4056856317u)--(        0.3104398166u,        0.3731843576u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2495238225u,        0.4079901991u)--(        0.3104398166u,        0.3731843576u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3726159544u,        0.4056856317u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3726159544u,        0.4056856317u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2495238225u,        0.4079901991u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2495238225u,        0.4079901991u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3104398166u,        0.3731843576u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3104398166u,        0.3731843576u) withcolor MidnightBlue;
draw (        0.4334123363u,        0.5232449880u)--(        0.4485349666u,        0.5917541917u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5083645275u,        0.6283960881u)--(        0.4485349666u,        0.5917541917u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4334123363u,        0.5232449880u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4334123363u,        0.5232449880u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5083645275u,        0.6283960881u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5083645275u,        0.6283960881u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4485349666u,        0.5917541917u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4485349666u,        0.5917541917u) withcolor MidnightBlue;
draw (        0.4334123363u,        0.5232449880u)--(        0.4932418972u,        0.5598868844u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5083645275u,        0.6283960881u)--(        0.4932418972u,        0.5598868844u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4334123363u,        0.5232449880u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4334123363u,        0.5232449880u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5083645275u,        0.6283960881u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5083645275u,        0.6283960881u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4932418972u,        0.5598868844u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4932418972u,        0.5598868844u) withcolor MidnightBlue;
draw (        0.4334123363u,        0.5232449880u)--(        0.3823092584u,        0.4751754513u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3171107118u,        0.5010859598u)--(        0.3823092584u,        0.4751754513u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4334123363u,        0.5232449880u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4334123363u,        0.5232449880u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3171107118u,        0.5010859598u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3171107118u,        0.5010859598u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3823092584u,        0.4751754513u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3823092584u,        0.4751754513u) withcolor MidnightBlue;
draw (        0.4334123363u,        0.5232449880u)--(        0.4485308658u,        0.5917550968u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3793531751u,        0.6034449238u)--(        0.4485308658u,        0.5917550968u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4334123363u,        0.5232449880u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4334123363u,        0.5232449880u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3793531751u,        0.6034449238u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3793531751u,        0.6034449238u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4485308658u,        0.5917550968u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4485308658u,        0.5917550968u) withcolor MidnightBlue;
draw (        0.5083645275u,        0.6283960881u)--(        0.5131204162u,        0.6983931341u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5703551299u,        0.7389688965u)--(        0.5131204162u,        0.6983931341u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5083645275u,        0.6283960881u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5083645275u,        0.6283960881u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5703551299u,        0.7389688965u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5703551299u,        0.7389688965u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5131204162u,        0.6983931341u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5131204162u,        0.6983931341u) withcolor MidnightBlue;
draw (        0.5083645275u,        0.6283960881u)--(        0.5655992412u,        0.6689718504u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5703551299u,        0.7389688965u)--(        0.5655992412u,        0.6689718504u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5083645275u,        0.6283960881u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5083645275u,        0.6283960881u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5703551299u,        0.7389688965u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5703551299u,        0.7389688965u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5655992412u,        0.6689718504u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5655992412u,        0.6689718504u) withcolor MidnightBlue;
draw (        0.5083645275u,        0.6283960881u)--(        0.4485316327u,        0.5917596360u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3793531751u,        0.6034449238u)--(        0.4485316327u,        0.5917596360u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5083645275u,        0.6283960881u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5083645275u,        0.6283960881u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3793531751u,        0.6034449238u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3793531751u,        0.6034449238u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4485316327u,        0.5917596360u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4485316327u,        0.5917596360u) withcolor MidnightBlue;
draw (        0.5083645275u,        0.6283960881u)--(        0.5131192340u,        0.6983932145u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4431540309u,        0.6931898212u)--(        0.5131192340u,        0.6983932145u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5083645275u,        0.6283960881u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5083645275u,        0.6283960881u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4431540309u,        0.6931898212u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4431540309u,        0.6931898212u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5131192340u,        0.6983932145u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5131192340u,        0.6983932145u) withcolor MidnightBlue;
draw (        0.5703551299u,        0.7389688965u)--(        0.5555555556u,        0.8075499103u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5703551299u,        0.7389688965u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5703551299u,        0.7389688965u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5555555556u,        0.8075499103u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5555555556u,        0.8075499103u) withcolor Yellow;
draw (        0.5703551299u,        0.7389688965u)--(        0.5851498996u,        0.8075496498u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6543932654u,        0.8188439505u)--(        0.5851498996u,        0.8075496498u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5703551299u,        0.7389688965u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5703551299u,        0.7389688965u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6543932654u,        0.8188439505u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6543932654u,        0.8188439505u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5851498996u,        0.8075496498u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5851498996u,        0.8075496498u) withcolor Green;
draw (        0.5703551299u,        0.7389688965u)--(        0.6395984957u,        0.7502631972u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6543932654u,        0.8188439505u)--(        0.6395984957u,        0.7502631972u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5703551299u,        0.7389688965u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5703551299u,        0.7389688965u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6543932654u,        0.8188439505u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6543932654u,        0.8188439505u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6395984957u,        0.7502631972u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6395984957u,        0.7502631972u) withcolor MidnightBlue;
draw (        0.5703551299u,        0.7389688965u)--(        0.5131191020u,        0.6983949880u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4431540309u,        0.6931898212u)--(        0.5131191020u,        0.6983949880u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5703551299u,        0.7389688965u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5703551299u,        0.7389688965u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4431540309u,        0.6931898212u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4431540309u,        0.6931898212u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5131191020u,        0.6983949880u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5131191020u,        0.6983949880u) withcolor MidnightBlue;
draw (        0.6543932654u,        0.8188439505u)--(        0.6237185104u,        0.8819412114u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6734560538u,        0.9314223362u)--(        0.6237185104u,        0.8819412114u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6543932654u,        0.8188439505u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6543932654u,        0.8188439505u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6734560538u,        0.9314223362u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6734560538u,        0.9314223362u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6237185104u,        0.8819412114u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6237185104u,        0.8819412114u) withcolor Green;
draw (        0.6543932654u,        0.8188439505u)--(        0.7041308088u,        0.8683250753u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6734560538u,        0.9314223362u)--(        0.7041308088u,        0.8683250753u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6543932654u,        0.8188439505u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6543932654u,        0.8188439505u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6734560538u,        0.9314223362u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6734560538u,        0.9314223362u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7041308088u,        0.8683250753u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7041308088u,        0.8683250753u) withcolor MidnightBlue;
draw (        0.1428250833u,        0.0395134359u)--(        0.1299863958u,        0.1084871485u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1996312069u,        0.1169609445u)--(        0.1299863958u,        0.1084871485u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1428250833u,        0.0395134359u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1428250833u,        0.0395134359u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1996312069u,        0.1169609445u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1996312069u,        0.1169609445u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1299863958u,        0.1084871485u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1299863958u,        0.1084871485u) withcolor MidnightBlue;
draw (        0.1428250833u,        0.0395134359u)--(        0.2124698944u,        0.0479872318u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1996312069u,        0.1169609445u)--(        0.2124698944u,        0.0479872318u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1428250833u,        0.0395134359u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1428250833u,        0.0395134359u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1996312069u,        0.1169609445u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1996312069u,        0.1169609445u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2124698944u,        0.0479872318u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2124698944u,        0.0479872318u) withcolor MidnightBlue;
draw (        0.1428250833u,        0.0395134359u)--(        0.1299568439u,        0.1084816412u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0795828341u,        0.1573146636u)--(        0.1299568439u,        0.1084816412u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1428250833u,        0.0395134359u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1428250833u,        0.0395134359u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0795828341u,        0.1573146636u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0795828341u,        0.1573146636u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1299568439u,        0.1084816412u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1299568439u,        0.1084816412u) withcolor MidnightBlue;
draw (        0.1996312069u,        0.1169609445u)--(        0.2134958581u,        0.1857357682u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2833059342u,        0.1927184758u)--(        0.2134958581u,        0.1857357682u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1996312069u,        0.1169609445u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1996312069u,        0.1169609445u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2833059342u,        0.1927184758u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2833059342u,        0.1927184758u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2134958581u,        0.1857357682u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2134958581u,        0.1857357682u) withcolor MidnightBlue;
draw (        0.1996312069u,        0.1169609445u)--(        0.2694412830u,        0.1239436520u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2833059342u,        0.1927184758u)--(        0.2694412830u,        0.1239436520u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1996312069u,        0.1169609445u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1996312069u,        0.1169609445u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2833059342u,        0.1927184758u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2833059342u,        0.1927184758u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2694412830u,        0.1239436520u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2694412830u,        0.1239436520u) withcolor MidnightBlue;
draw (        0.1996312069u,        0.1169609445u)--(        0.1299836972u,        0.1085093571u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0795828341u,        0.1573146636u)--(        0.1299836972u,        0.1085093571u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1996312069u,        0.1169609445u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1996312069u,        0.1169609445u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0795828341u,        0.1573146636u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0795828341u,        0.1573146636u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1299836972u,        0.1085093571u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1299836972u,        0.1085093571u) withcolor MidnightBlue;
draw (        0.1996312069u,        0.1169609445u)--(        0.2134875238u,        0.1857374478u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1553269137u,        0.2249745269u)--(        0.2134875238u,        0.1857374478u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1996312069u,        0.1169609445u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1996312069u,        0.1169609445u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1553269137u,        0.2249745269u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1553269137u,        0.2249745269u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2134875238u,        0.1857374478u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2134875238u,        0.1857374478u) withcolor MidnightBlue;
draw (        0.2833059342u,        0.1927184758u)--(        0.2580690670u,        0.2581807204u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3120121119u,        0.3030399763u)--(        0.2580690670u,        0.2581807204u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2833059342u,        0.1927184758u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2833059342u,        0.1927184758u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3120121119u,        0.3030399763u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3120121119u,        0.3030399763u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2580690670u,        0.2581807204u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2580690670u,        0.2581807204u) withcolor MidnightBlue;
draw (        0.2833059342u,        0.1927184758u)--(        0.3372489791u,        0.2375777317u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3120121119u,        0.3030399763u)--(        0.3372489791u,        0.2375777317u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2833059342u,        0.1927184758u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2833059342u,        0.1927184758u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3120121119u,        0.3030399763u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3120121119u,        0.3030399763u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3372489791u,        0.2375777317u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3372489791u,        0.2375777317u) withcolor MidnightBlue;
draw (        0.2833059342u,        0.1927184758u)--(        0.2134946364u,        0.1857479928u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1553269137u,        0.2249745269u)--(        0.2134946364u,        0.1857479928u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2833059342u,        0.1927184758u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2833059342u,        0.1927184758u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1553269137u,        0.2249745269u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1553269137u,        0.2249745269u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2134946364u,        0.1857479928u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2134946364u,        0.1857479928u) withcolor MidnightBlue;
draw (        0.2833059342u,        0.1927184758u)--(        0.2251382116u,        0.2319450098u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1553269137u,        0.2249745269u)--(        0.2251382116u,        0.2319450098u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2833059342u,        0.1927184758u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2833059342u,        0.1927184758u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1553269137u,        0.2249745269u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1553269137u,        0.2249745269u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2251382116u,        0.2319450098u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2251382116u,        0.2319450098u) withcolor MidnightBlue;
draw (        0.2833059342u,        0.1927184758u)--(        0.2251511139u,        0.2319641356u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1999059948u,        0.2974231983u)--(        0.2251511139u,        0.2319641356u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2833059342u,        0.1927184758u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2833059342u,        0.1927184758u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1999059948u,        0.2974231983u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1999059948u,        0.2974231983u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2251511139u,        0.2319641356u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2251511139u,        0.2319641356u) withcolor MidnightBlue;
draw (        0.2833059342u,        0.1927184758u)--(        0.2580608151u,        0.2581775385u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1999059948u,        0.2974231983u)--(        0.2580608151u,        0.2581775385u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2833059342u,        0.1927184758u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2833059342u,        0.1927184758u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1999059948u,        0.2974231983u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1999059948u,        0.2974231983u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2580608151u,        0.2581775385u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2580608151u,        0.2581775385u) withcolor MidnightBlue;
draw (        0.3120121119u,        0.3030399763u)--(        0.2580656968u,        0.2581847733u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1999059948u,        0.2974231983u)--(        0.2580656968u,        0.2581847733u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3120121119u,        0.3030399763u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3120121119u,        0.3030399763u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1999059948u,        0.2974231983u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1999059948u,        0.2974231983u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2580656968u,        0.2581847733u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2580656968u,        0.2581847733u) withcolor MidnightBlue;
draw (        0.3120121119u,        0.3030399763u)--(        0.3104377477u,        0.3731807370u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2495238225u,        0.4079901991u)--(        0.3104377477u,        0.3731807370u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3120121119u,        0.3030399763u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3120121119u,        0.3030399763u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2495238225u,        0.4079901991u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2495238225u,        0.4079901991u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3104377477u,        0.3731807370u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3104377477u,        0.3731807370u) withcolor MidnightBlue;
draw (        0.3171107118u,        0.5010859598u)--(        0.3170207866u,        0.5712443295u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3793531751u,        0.6034449238u)--(        0.3170207866u,        0.5712443295u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3171107118u,        0.5010859598u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3171107118u,        0.5010859598u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3793531751u,        0.6034449238u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3793531751u,        0.6034449238u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3170207866u,        0.5712443295u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3170207866u,        0.5712443295u) withcolor MidnightBlue;
draw (        0.3171107118u,        0.5010859598u)--(        0.2508118346u,        0.4781368024u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2495238225u,        0.4079901991u)--(        0.2508118346u,        0.4781368024u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3171107118u,        0.5010859598u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3171107118u,        0.5010859598u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2495238225u,        0.4079901991u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2495238225u,        0.4079901991u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2508118346u,        0.4781368024u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2508118346u,        0.4781368024u) withcolor Green;
draw (        0.3171107118u,        0.5010859598u)--(        0.2502082214u,        0.5222108673u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2501164158u,        0.5923692347u)--(        0.2502082214u,        0.5222108673u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3171107118u,        0.5010859598u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3171107118u,        0.5010859598u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2501164158u,        0.5923692347u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2501164158u,        0.5923692347u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2502082214u,        0.5222108673u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2502082214u,        0.5222108673u) withcolor Green;
draw (        0.3171107118u,        0.5010859598u)--(        0.3170189062u,        0.5712443271u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2501164158u,        0.5923692347u)--(        0.3170189062u,        0.5712443271u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3171107118u,        0.5010859598u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3171107118u,        0.5010859598u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2501164158u,        0.5923692347u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2501164158u,        0.5923692347u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3170189062u,        0.5712443271u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3170189062u,        0.5712443271u) withcolor MidnightBlue;
draw (        0.3793531751u,        0.6034449238u)--(        0.3758111591u,        0.6735138831u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4431540309u,        0.6931898212u)--(        0.3758111591u,        0.6735138831u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3793531751u,        0.6034449238u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3793531751u,        0.6034449238u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4431540309u,        0.6931898212u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4431540309u,        0.6931898212u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3758111591u,        0.6735138831u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3758111591u,        0.6735138831u) withcolor MidnightBlue;
draw (        0.3793531751u,        0.6034449238u)--(        0.3170196200u,        0.5712465879u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2501164158u,        0.5923692347u)--(        0.3170196200u,        0.5712465879u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3793531751u,        0.6034449238u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3793531751u,        0.6034449238u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2501164158u,        0.5923692347u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2501164158u,        0.5923692347u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3170196200u,        0.5712465879u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3170196200u,        0.5712465879u) withcolor MidnightBlue;
draw (        0.3793531751u,        0.6034449238u)--(        0.3124499710u,        0.6245675705u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2501164158u,        0.5923692347u)--(        0.3124499710u,        0.6245675705u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3793531751u,        0.6034449238u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3793531751u,        0.6034449238u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2501164158u,        0.5923692347u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2501164158u,        0.5923692347u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3124499710u,        0.6245675705u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3124499710u,        0.6245675705u) withcolor MidnightBlue;
draw (        0.3793531751u,        0.6034449238u)--(        0.3124506716u,        0.6245697895u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3089074226u,        0.6946386864u)--(        0.3124506716u,        0.6245697895u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3793531751u,        0.6034449238u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3793531751u,        0.6034449238u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3089074226u,        0.6946386864u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3089074226u,        0.6946386864u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3124506716u,        0.6245697895u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3124506716u,        0.6245697895u) withcolor MidnightBlue;
draw (        0.3793531751u,        0.6034449238u)--(        0.3758099261u,        0.6735138207u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3089074226u,        0.6946386864u)--(        0.3758099261u,        0.6735138207u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3793531751u,        0.6034449238u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3793531751u,        0.6034449238u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3089074226u,        0.6946386864u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3089074226u,        0.6946386864u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3758099261u,        0.6735138207u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3758099261u,        0.6735138207u) withcolor MidnightBlue;
draw (        0.4431540309u,        0.6931898212u)--(        0.3758105760u,        0.6735158789u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3089074226u,        0.6946386864u)--(        0.3758105760u,        0.6735158789u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4431540309u,        0.6931898212u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4431540309u,        0.6931898212u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3089074226u,        0.6946386864u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3089074226u,        0.6946386864u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3758105760u,        0.6735158789u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3758105760u,        0.6735158789u) withcolor MidnightBlue;
draw (        0.4431540309u,        0.6931898212u)--(        0.4741241167u,        0.7561426505u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4040594714u,        0.7597689986u)--(        0.4741241167u,        0.7561426505u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4431540309u,        0.6931898212u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4431540309u,        0.6931898212u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4040594714u,        0.7597689986u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4040594714u,        0.7597689986u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4741241167u,        0.7561426505u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4741241167u,        0.7561426505u) withcolor Green;
draw (        0.0423981266u,        0.0559113870u)--(        0.0000000000u,        0.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0423981266u,        0.0559113870u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0423981266u,        0.0559113870u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.0795828341u,        0.1573146636u)--(        0.0852053135u,        0.2272474360u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1553269137u,        0.2249745269u)--(        0.0852053135u,        0.2272474360u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0795828341u,        0.1573146636u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0795828341u,        0.1573146636u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1553269137u,        0.2249745269u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1553269137u,        0.2249745269u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0852053135u,        0.2272474360u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0852053135u,        0.2272474360u) withcolor MidnightBlue;
draw (        0.0795828341u,        0.1573146636u)--(        0.0851893119u,        0.2272487207u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0343451404u,        0.2755920259u)--(        0.0851893119u,        0.2272487207u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0795828341u,        0.1573146636u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0795828341u,        0.1573146636u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0343451404u,        0.2755920259u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0343451404u,        0.2755920259u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0851893119u,        0.2272487207u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0851893119u,        0.2272487207u) withcolor MidnightBlue;
draw (        0.1553269137u,        0.2249745269u)--(        0.2251368497u,        0.2319586361u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1999059948u,        0.2974231983u)--(        0.2251368497u,        0.2319586361u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1553269137u,        0.2249745269u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1553269137u,        0.2249745269u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1999059948u,        0.2974231983u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1999059948u,        0.2974231983u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2251368497u,        0.2319586361u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2251368497u,        0.2319586361u) withcolor MidnightBlue;
draw (        0.1553269137u,        0.2249745269u)--(        0.0852059242u,        0.2272661983u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0343451404u,        0.2755920259u)--(        0.0852059242u,        0.2272661983u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1553269137u,        0.2249745269u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1553269137u,        0.2249745269u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0343451404u,        0.2755920259u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0343451404u,        0.2755920259u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0852059242u,        0.2272661983u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0852059242u,        0.2272661983u) withcolor MidnightBlue;
draw (        0.1553269137u,        0.2249745269u)--(        0.1044661300u,        0.2733003544u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0343451404u,        0.2755920259u)--(        0.1044661300u,        0.2733003544u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1553269137u,        0.2249745269u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1553269137u,        0.2249745269u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0343451404u,        0.2755920259u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0343451404u,        0.2755920259u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1044661300u,        0.2733003544u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1044661300u,        0.2733003544u) withcolor MidnightBlue;
draw (        0.1553269137u,        0.2249745269u)--(        0.1044712700u,        0.2733057635u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0884715503u,        0.3416154546u)--(        0.1044712700u,        0.2733057635u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1553269137u,        0.2249745269u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1553269137u,        0.2249745269u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0884715503u,        0.3416154546u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0884715503u,        0.3416154546u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1044712700u,        0.2733057635u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1044712700u,        0.2733057635u) withcolor MidnightBlue;
draw (        0.1999059948u,        0.2974231983u)--(        0.1924527400u,        0.3671846054u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2495238225u,        0.4079901991u)--(        0.1924527400u,        0.3671846054u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1999059948u,        0.2974231983u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1999059948u,        0.2974231983u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2495238225u,        0.4079901991u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2495238225u,        0.4079901991u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1924527400u,        0.3671846054u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1924527400u,        0.3671846054u) withcolor Green;
draw (        0.2501164158u,        0.5923692347u)--(        0.3124488717u,        0.6245696985u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3089074226u,        0.6946386864u)--(        0.3124488717u,        0.6245696985u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2501164158u,        0.5923692347u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2501164158u,        0.5923692347u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3089074226u,        0.6946386864u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3089074226u,        0.6946386864u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3124488717u,        0.6245696985u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3124488717u,        0.6245696985u) withcolor MidnightBlue;
draw (        0.2501164158u,        0.5923692347u)--(        0.1924531339u,        0.6323336044u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2001368279u,        0.7020700066u)--(        0.1924531339u,        0.6323336044u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2501164158u,        0.5923692347u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2501164158u,        0.5923692347u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2001368279u,        0.7020700066u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2001368279u,        0.7020700066u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1924531339u,        0.6323336044u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1924531339u,        0.6323336044u) withcolor Green;
draw (        0.3089074226u,        0.6946386864u)--(        0.3339023294u,        0.7601936982u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4040594714u,        0.7597689986u)--(        0.3339023294u,        0.7601936982u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3089074226u,        0.6946386864u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3089074226u,        0.6946386864u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4040594714u,        0.7597689986u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4040594714u,        0.7597689986u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3339023294u,        0.7601936982u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3339023294u,        0.7601936982u) withcolor MidnightBlue;
draw (        0.3089074226u,        0.6946386864u)--(        0.2575325765u,        0.7424176580u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2001368279u,        0.7020700066u)--(        0.2575325765u,        0.7424176580u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3089074226u,        0.6946386864u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3089074226u,        0.6946386864u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2001368279u,        0.7020700066u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2001368279u,        0.7020700066u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2575325765u,        0.7424176580u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2575325765u,        0.7424176580u) withcolor MidnightBlue;
draw (        0.3089074226u,        0.6946386864u)--(        0.2575346931u,        0.7424199338u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2825292671u,        0.8079750725u)--(        0.2575346931u,        0.7424199338u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3089074226u,        0.6946386864u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3089074226u,        0.6946386864u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2825292671u,        0.8079750725u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2825292671u,        0.8079750725u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2575346931u,        0.7424199338u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2575346931u,        0.7424199338u) withcolor MidnightBlue;
draw (        0.3089074226u,        0.6946386864u)--(        0.3339019966u,        0.7601938251u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2825292671u,        0.8079750725u)--(        0.3339019966u,        0.7601938251u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3089074226u,        0.6946386864u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3089074226u,        0.6946386864u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2825292671u,        0.8079750725u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2825292671u,        0.8079750725u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3339019966u,        0.7601938251u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3339019966u,        0.7601938251u) withcolor MidnightBlue;
draw (        0.4040594714u,        0.7597689986u)--(        0.3339023324u,        0.7601941861u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2825292671u,        0.8079750725u)--(        0.3339023324u,        0.7601941861u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4040594714u,        0.7597689986u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4040594714u,        0.7597689986u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2825292671u,        0.8079750725u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2825292671u,        0.8079750725u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3339023324u,        0.7601941861u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3339023324u,        0.7601941861u) withcolor MidnightBlue;
draw (        0.4040594714u,        0.7597689986u)--(        0.3526864061u,        0.8075498849u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2825292671u,        0.8079750725u)--(        0.3526864061u,        0.8075498849u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4040594714u,        0.7597689986u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4040594714u,        0.7597689986u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2825292671u,        0.8079750725u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2825292671u,        0.8079750725u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3526864061u,        0.8075498849u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3526864061u,        0.8075498849u) withcolor Green;
draw (        0.0343451404u,        0.2755920259u)--(        0.1044662683u,        0.2733045922u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0884715503u,        0.3416154546u)--(        0.1044662683u,        0.2733045922u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0343451404u,        0.2755920259u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0343451404u,        0.2755920259u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0884715503u,        0.3416154546u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0884715503u,        0.3416154546u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1044662683u,        0.2733045922u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1044662683u,        0.2733045922u) withcolor MidnightBlue;
draw (        0.2001368279u,        0.7020700066u)--(        0.2251364797u,        0.7676232089u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2825292671u,        0.8079750725u)--(        0.2251364797u,        0.7676232089u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2001368279u,        0.7020700066u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2001368279u,        0.7020700066u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2825292671u,        0.8079750725u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2825292671u,        0.8079750725u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2251364797u,        0.7676232089u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2251364797u,        0.7676232089u) withcolor MidnightBlue;
draw (        0.2001368279u,        0.7020700066u)--(        0.2575296153u,        0.7424218701u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2825292671u,        0.8079750725u)--(        0.2575296153u,        0.7424218701u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2001368279u,        0.7020700066u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2001368279u,        0.7020700066u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2825292671u,        0.8079750725u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2825292671u,        0.8079750725u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2575296153u,        0.7424218701u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2575296153u,        0.7424218701u) withcolor MidnightBlue;
draw (        0.2001368279u,        0.7020700066u)--(        0.2251262971u,        0.7676270913u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1552978314u,        0.7744234362u)--(        0.2251262971u,        0.7676270913u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2001368279u,        0.7020700066u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2001368279u,        0.7020700066u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1552978314u,        0.7744234362u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1552978314u,        0.7744234362u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2251262971u,        0.7676270913u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2251262971u,        0.7676270913u) withcolor MidnightBlue;
draw (        0.2825292671u,        0.8079750725u)--(        0.2693492986u,        0.8768843851u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3348659826u,        0.9019795847u)--(        0.2693492986u,        0.8768843851u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2825292671u,        0.8079750725u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2825292671u,        0.8079750725u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3348659826u,        0.9019795847u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3348659826u,        0.9019795847u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2693492986u,        0.8768843851u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2693492986u,        0.8768843851u) withcolor MidnightBlue;
draw (        0.2825292671u,        0.8079750725u)--(        0.2251272031u,        0.7676364064u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1552978314u,        0.7744234362u)--(        0.2251272031u,        0.7676364064u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2825292671u,        0.8079750725u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2825292671u,        0.8079750725u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1552978314u,        0.7744234362u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1552978314u,        0.7744234362u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2251272031u,        0.7676364064u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2251272031u,        0.7676364064u) withcolor MidnightBlue;
draw (        0.2825292671u,        0.8079750725u)--(        0.2126998954u,        0.8147621024u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1552978314u,        0.7744234362u)--(        0.2126998954u,        0.8147621024u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2825292671u,        0.8079750725u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2825292671u,        0.8079750725u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1552978314u,        0.7744234362u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1552978314u,        0.7744234362u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2126998954u,        0.8147621024u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2126998954u,        0.8147621024u) withcolor MidnightBlue;
draw (        0.2825292671u,        0.8079750725u)--(        0.2127005445u,        0.8147687782u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1995151346u,        0.8836770498u)--(        0.2127005445u,        0.8147687782u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2825292671u,        0.8079750725u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2825292671u,        0.8079750725u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1995151346u,        0.8836770498u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1995151346u,        0.8836770498u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2127005445u,        0.8147687782u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2127005445u,        0.8147687782u) withcolor MidnightBlue;
draw (        0.2825292671u,        0.8079750725u)--(        0.2693438572u,        0.8768833441u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1995151346u,        0.8836770498u)--(        0.2693438572u,        0.8768833441u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2825292671u,        0.8079750725u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2825292671u,        0.8079750725u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1995151346u,        0.8836770498u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1995151346u,        0.8836770498u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2693438572u,        0.8768833441u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2693438572u,        0.8768833441u) withcolor MidnightBlue;
draw (        0.3348659826u,        0.9019795847u)--(        0.2693450395u,        0.8768955071u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1995151346u,        0.8836770498u)--(        0.2693450395u,        0.8768955071u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3348659826u,        0.9019795847u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3348659826u,        0.9019795847u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1995151346u,        0.8836770498u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1995151346u,        0.8836770498u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2693450395u,        0.8768955071u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2693450395u,        0.8768955071u) withcolor MidnightBlue;
draw (        0.0881652236u,        0.6578469932u)--(        0.1044434736u,        0.7260908465u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1552978314u,        0.7744234362u)--(        0.1044434736u,        0.7260908465u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0881652236u,        0.6578469932u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0881652236u,        0.6578469932u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1552978314u,        0.7744234362u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1552978314u,        0.7744234362u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1044434736u,        0.7260908465u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1044434736u,        0.7260908465u) withcolor MidnightBlue;
draw (        0.0881652236u,        0.6578469932u)--(        0.1044406432u,        0.7260915216u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0342998121u,        0.7245202987u)--(        0.1044406432u,        0.7260915216u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0881652236u,        0.6578469932u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0881652236u,        0.6578469932u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0342998121u,        0.7245202987u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0342998121u,        0.7245202987u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1044406432u,        0.7260915216u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1044406432u,        0.7260915216u) withcolor MidnightBlue;
draw (        0.1552978314u,        0.7744234362u)--(        0.2126958365u,        0.8147678775u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1995151346u,        0.8836770498u)--(        0.2126958365u,        0.8147678775u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1552978314u,        0.7744234362u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1552978314u,        0.7744234362u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1995151346u,        0.8836770498u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1995151346u,        0.8836770498u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2126958365u,        0.8147678775u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2126958365u,        0.8147678775u) withcolor MidnightBlue;
draw (        0.1552978314u,        0.7744234362u)--(        0.1044405903u,        0.7260938804u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0342998121u,        0.7245202987u)--(        0.1044405903u,        0.7260938804u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1552978314u,        0.7744234362u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1552978314u,        0.7744234362u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0342998121u,        0.7245202987u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0342998121u,        0.7245202987u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1044405903u,        0.7260938804u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1044405903u,        0.7260938804u) withcolor MidnightBlue;
draw (        0.1552978314u,        0.7744234362u)--(        0.0851570531u,        0.7728498545u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0342998121u,        0.7245202987u)--(        0.0851570531u,        0.7728498545u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1552978314u,        0.7744234362u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1552978314u,        0.7744234362u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0342998121u,        0.7245202987u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0342998121u,        0.7245202987u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0851570531u,        0.7728498545u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0851570531u,        0.7728498545u) withcolor MidnightBlue;
draw (        0.1552978314u,        0.7744234362u)--(        0.0851567837u,        0.7728619105u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0797963749u,        0.8428152589u)--(        0.0851567837u,        0.7728619105u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1552978314u,        0.7744234362u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1552978314u,        0.7744234362u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0797963749u,        0.8428152589u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0797963749u,        0.8428152589u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0851567837u,        0.7728619105u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0851567837u,        0.7728619105u) withcolor MidnightBlue;
draw (        0.1995151346u,        0.8836770498u)--(        0.2126939661u,        0.9525865799u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2802505534u,        0.9715157229u)--(        0.2126939661u,        0.9525865799u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1995151346u,        0.8836770498u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1995151346u,        0.8836770498u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2802505534u,        0.9715157229u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2802505534u,        0.9715157229u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2126939661u,        0.9525865799u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2126939661u,        0.9525865799u) withcolor MidnightBlue;
draw (        0.1995151346u,        0.8836770498u)--(        0.1298494171u,        0.8919772148u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0797963749u,        0.8428152589u)--(        0.1298494171u,        0.8919772148u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1995151346u,        0.8836770498u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1995151346u,        0.8836770498u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0797963749u,        0.8428152589u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0797963749u,        0.8428152589u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1298494171u,        0.8919772148u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1298494171u,        0.8919772148u) withcolor MidnightBlue;
draw (        0.1995151346u,        0.8836770498u)--(        0.1298510715u,        0.8919910883u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1430256671u,        0.9609014283u)--(        0.1298510715u,        0.8919910883u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1995151346u,        0.8836770498u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1995151346u,        0.8836770498u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1430256671u,        0.9609014283u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1430256671u,        0.9609014283u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1298510715u,        0.8919910883u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1298510715u,        0.8919910883u) withcolor MidnightBlue;
draw (        0.1995151346u,        0.8836770498u)--(        0.2126897303u,        0.9525873898u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1430256671u,        0.9609014283u)--(        0.2126897303u,        0.9525873898u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1995151346u,        0.8836770498u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1995151346u,        0.8836770498u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1430256671u,        0.9609014283u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1430256671u,        0.9609014283u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2126897303u,        0.9525873898u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2126897303u,        0.9525873898u) withcolor MidnightBlue;
draw (        0.2802505534u,        0.9715157229u)--(        0.2126909275u,        0.9525974278u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1430256671u,        0.9609014283u)--(        0.2126909275u,        0.9525974278u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2802505534u,        0.9715157229u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2802505534u,        0.9715157229u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1430256671u,        0.9609014283u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1430256671u,        0.9609014283u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2126909275u,        0.9525974278u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2126909275u,        0.9525974278u) withcolor MidnightBlue;
draw (        0.0342998121u,        0.7245202987u)--(        0.0851463528u,        0.7728611120u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0797963749u,        0.8428152589u)--(        0.0851463528u,        0.7728611120u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0342998121u,        0.7245202987u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0342998121u,        0.7245202987u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0797963749u,        0.8428152589u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0797963749u,        0.8428152589u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0851463528u,        0.7728611120u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0851463528u,        0.7728611120u) withcolor MidnightBlue;
draw (        0.0797963749u,        0.8428152589u)--(        0.1298322471u,        0.8919946899u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1430256671u,        0.9609014283u)--(        0.1298322471u,        0.8919946899u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0797963749u,        0.8428152589u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0797963749u,        0.8428152589u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1430256671u,        0.9609014283u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1430256671u,        0.9609014283u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1298322471u,        0.8919946899u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1298322471u,        0.8919946899u) withcolor MidnightBlue;
draw (        0.0797963749u,        0.8428152589u)--(        0.0190427854u,        0.8779038082u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0423686090u,        0.9440711029u)--(        0.0190427854u,        0.8779038082u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0797963749u,        0.8428152589u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0797963749u,        0.8428152589u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0423686090u,        0.9440711029u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0423686090u,        0.9440711029u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0190427854u,        0.8779038082u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0190427854u,        0.8779038082u) withcolor Green;
draw (        0.0423686090u,        0.9440711029u)--(        0.0000000000u,        1.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0423686090u,        0.9440711029u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0423686090u,        0.9440711029u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        1.0000000000u) withcolor Yellow;
% r =         0.0701584274
endfig;
end;
