# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# Try cluster-based permutation test
function welch_statistic(x,y) # degrees of freedom depend on the variances...
    nx=length(x); ny=length(y);
    μx=mean(x); μy=mean(y);
    σx=std(x); σy=std(y);
    d = σx^2/nx + σy^2/ny
    tw = (μx-μy)/sqrt(d)

    sx = (σx^2/nx)^2/(nx-1) ; sy = (σy^2/ny)^2/(ny-1)
    df = d^2/(sx+sy)
    return tw, df
end


id=100; data="cortex_";
t, f, tfhm = timefreq(id, data*"pre");
idx, f_pre = movfilter(t, tfhm, "pre");

t, f, tfhm = timefreq(id, data*"post");
idx, f_post = movfilter(t, tfhm, "post");

plot([f_pre[100,:],f_post[100,:]])

l=2000
ps=[ pvalue(UnequalVarianceTTest(f_pre[i,:],f_post[i,:])) for i in 1:l]

pss = @. (ps<0.05) 

ncluster = 0
cluster  = 0
for i in 1:l
    if ps[i]<0.05 && cluster==0
        ncluster+=1
    end
    cluster = (ps[i]<0.05 ? ncluster : 0) 
    cluster_list[i] = cluster
end
# sum ts of each cluster
# randomize and compute p-value (use chuncks of hypothesis test packg)


# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
