prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(3);
u := 20cm;
path pols[];
pair p[];
R := .2886440620u;
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

x1 = .2500215274u;
y1 = .1443499960u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .7499784727u;
y2 = .1443499961u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .4999999999u;
y3 = .5773254115u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.2500215274u,        0.1443499960u)--(        0.0000000000u,        0.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2500215274u,        0.1443499960u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500215274u,        0.1443499960u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.2500215274u,        0.1443499960u)--(        0.5000000000u,        0.2886627042u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.7499784727u,        0.1443499961u)--(        0.5000000000u,        0.2886627042u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2500215274u,        0.1443499960u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500215274u,        0.1443499960u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7499784727u,        0.1443499961u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7499784727u,        0.1443499961u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.2886627042u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.2886627042u) withcolor MidnightBlue;
draw (        0.2500215274u,        0.1443499960u)--(        0.5000000001u,        0.0000372880u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.7499784727u,        0.1443499961u)--(        0.5000000001u,        0.0000372880u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2500215274u,        0.1443499960u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500215274u,        0.1443499960u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7499784727u,        0.1443499961u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7499784727u,        0.1443499961u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.0000372880u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.0000372880u) withcolor Green;
draw (        0.2500215274u,        0.1443499960u)--(        0.2500322924u,        0.4329940578u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.4999999999u,        0.5773254115u)--(        0.2500322924u,        0.4329940578u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2500215274u,        0.1443499960u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500215274u,        0.1443499960u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4999999999u,        0.5773254115u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999999u,        0.5773254115u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2500322924u,        0.4329940578u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500322924u,        0.4329940578u) withcolor Green;
draw (        0.2500215274u,        0.1443499960u)--(        0.4999892350u,        0.2886813498u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.4999999999u,        0.5773254115u)--(        0.4999892350u,        0.2886813498u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2500215274u,        0.1443499960u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500215274u,        0.1443499960u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4999999999u,        0.5773254115u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999999u,        0.5773254115u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4999892350u,        0.2886813498u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999892350u,        0.2886813498u) withcolor MidnightBlue;
draw (        0.7499784727u,        0.1443499961u)--(        1.0000000000u,        0.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7499784727u,        0.1443499961u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7499784727u,        0.1443499961u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.7499784727u,        0.1443499961u)--(        0.5000107651u,        0.2886813498u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.4999999999u,        0.5773254115u)--(        0.5000107651u,        0.2886813498u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7499784727u,        0.1443499961u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7499784727u,        0.1443499961u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4999999999u,        0.5773254115u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999999u,        0.5773254115u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000107651u,        0.2886813498u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000107651u,        0.2886813498u) withcolor MidnightBlue;
draw (        0.7499784727u,        0.1443499961u)--(        0.7499677076u,        0.4329940579u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.4999999999u,        0.5773254115u)--(        0.7499677076u,        0.4329940579u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7499784727u,        0.1443499961u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7499784727u,        0.1443499961u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4999999999u,        0.5773254115u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999999u,        0.5773254115u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7499677076u,        0.4329940579u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7499677076u,        0.4329940579u) withcolor Green;
draw (        0.4999999999u,        0.5773254115u)--(        0.5000000000u,        0.8660254038u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4999999999u,        0.5773254115u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999999u,        0.5773254115u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.8660254038u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.8660254038u) withcolor Yellow;
% r =         0.2886440620
endfig;
end;
