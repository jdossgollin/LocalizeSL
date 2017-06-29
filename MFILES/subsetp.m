function p=subsetp(p,targregions)

% Last updated by Robert Kopp, robert-dot-kopp-at-rutgers-dot-edu, Fri Sep 30 09:07:23 EDT 2016

sub=find(ismember(p.targregions,targregions));

fn=fieldnames(p);
oldn=length(p.targregions);
for sss=1:length(fn)
    if size(p.(fn{sss}),1)==oldn
        p.(fn{sss})=p.(fn{sss})(sub,:);
    elseif size(p.(fn{sss}),2)==oldn
        if ndims(p.(fn{sss}))==2
            p.(fn{sss})=p.(fn{sss})(:,sub);
        elseif ndims(p.(fn{sss}))==3
            p.(fn{sss})=p.(fn{sss})(:,sub,:);
        end
    end
end
