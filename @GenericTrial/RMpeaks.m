function rmp = 
roi 
for i = 1 : 45
rm1 = gt.pfObject.rateMap{cluid(i)}; 
if ~isempty(rm1)
rmp1(i) = nanmax(rm1(:));            
end

end