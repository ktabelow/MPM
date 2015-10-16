function D = sdiag(v)
D = spdiags(v(:),0,numel(v),numel(v));
end