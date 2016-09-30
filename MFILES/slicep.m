function slicedp=slicep(p,sub)

% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Sep 30 09:07:23 EDT 2016

slicedp=p;
N=length(p.targregions);
fn=fieldnames(p);
for www=1:length(fn)
    if size(p.(fn{www}),1)==N
        slicedp.(fn{www})=slicedp.(fn{www})(sub,:,:);
    end
    if size(p.(fn{www}),2)==N
        slicedp.(fn{www})=slicedp.(fn{www})(:,sub,:);
    end
end
