% function to plot geometry 
function plotgeom(geom, fignum, colorofpoint) 
figure(fignum) 
hold on 
for j=1:3:18 
plot3(geom(j), geom(j+1), geom(j+2), ... 
'o','MarkerSize',10,'MarkerEdgeColor','k', ... 
'MarkerFaceColor',colorofpoint) 
end