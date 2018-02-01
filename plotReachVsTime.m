% Author:   Nicole Chan
% Created:  1/28/18
% Description: Plots 1-dimensional projections of reachsets 
% (reach.projection(dim)) with respect to time (tVec) and optional label
% (labelString) for the y-axis
%
function figHand = plotReachVsTime(reach,dim,tVec,labelString)
    figHand = figure();
    hold on;
    grid on; grid minor;
    if length(dim)==1
        outp = zeros(length(tVec),2);
        proj = reach.projection(dim);
        % % WEIRD MPT PROJECTION METHOD BEHAVIOR: some Polyhedron in the projection will have duplicate vertices % %
        for i=1:length(outp)
            vert = proj(i).V;
            outp(i,:) = [max(vert),min(vert)];
        end
        fill([tVec;flipud(tVec)],[outp(:,1);flipud(outp(:,2))],'b','EdgeColor','none');
    elseif length(dim)==2
        plot(reach.projection(dim),'color','b','edgealpha',0);
    else
        error('Please specify up to 2 dimensions to plot.');
    end
    
    if nargin==4 && isa(labelString,'char')
        xlabel('Time (s)');
        ylabel(labelString);
    end
end