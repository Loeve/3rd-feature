function displayIter(i, b, weights, f, g)
    
if 1 % mod(i,10)==0
  fprintf('iter %d, batch %d, f=%5.3f\n',...
	  i-1,b,f);
end
  
