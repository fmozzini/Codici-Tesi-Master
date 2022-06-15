function [medianvalue,maxvalue,M_tag] = maxmedian(tag,R)
 M = zeros(length(R), 2);
 c = 0;
 ln = 1;
 val = R(1,3);
 for i = 2:length(R)
     c = c + 1;
     cur_val = R(i,3);
     if cur_val ~= val
         M(ln, 1) = val;
         M(ln, 2) = c;
         c = 0;
         val = cur_val;
         ln = ln+1;
     end
 end
 M(ln, 1) = val;
 M(ln, 2) = c;
 M(ln:end,:) = [];

 posM_tag = find(M(:,1) == tag); 
 M_tag = M(posM_tag,:);
 medianvalue = median(M_tag(:,2)); 
 maxvalue = max(M_tag(:,2));
end