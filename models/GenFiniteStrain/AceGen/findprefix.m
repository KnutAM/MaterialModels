function found = findprefix(str,pref)
%Look for the pref in str
tmp = regexp(str,pref); 
%Check if it was found and that it was first
found = (~isempty(tmp))&&(tmp==1);
