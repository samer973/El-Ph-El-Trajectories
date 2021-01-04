function cutoff=sntrajcut(t,q)
if(q(2)+t>q(3)-t)
    cutoff=0;
else
    cutoff=1;
end