% Encompass rectangleA inside rectangleB.
% The rectangle format is [left top width height].
function [left,top,width,height] = encompass(rectangleA, rectangleB)
    left = rectangleA(1); top = rectangleA(2);
    width = rectangleA(3); height = rectangleA(4);
    
    if (left < rectangleB(1))
        left = rectangleB(1);
    end
    
    if (top < rectangleB(2))
        top = rectangleB(2);
    end
    
    if (left + width > rectangleB(1) + rectangleB(3))
        width = rectangleB(1) + rectangleB(3) - left - 1;
    end
    
    if (top + height > rectangleB(2) + rectangleB(4))
        height = rectangleB(2) + rectangleB(4) - top - 1;
    end   

end
