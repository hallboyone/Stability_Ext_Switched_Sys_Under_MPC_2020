function terminalProgressBar(curIdx, endIdx)
if curIdx == 1
    fprintf('|+>%54c|', ' ');
elseif curIdx == endIdx
    fprintf([repmat('\b', 1, 3), '++|\n']); 
else
    digit = round(curIdx/endIdx*55);
    fprintf([repmat('\b', 1, 58-digit), '+>', repmat(' ', 1, 55-digit), '|']);
end
end