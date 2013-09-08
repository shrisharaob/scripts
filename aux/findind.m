
function ind = findind(f,per)
ind =[];f=f(:);
for ii=1:size(per,1)
   ind = [ind; find(f>per(ii,1) & f<per(ii,2))];
end
end
