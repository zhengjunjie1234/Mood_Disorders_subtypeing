function [regional_nd_p,regional_nd_stats] = normative_deviation_hc(age_sex_edu_hc,hc_data,mean_age)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明

%%% normative model

b0 = ones(size(hc_data,1),1);
b1 = age_sex_edu_hc(:,2);
b1= b1-1;
b2 = age_sex_edu_hc(:,1)-mean_age;
b3 = b1.*b2;
b4 = b2.*b2;
x_hc = [b0,b1,b2,b3,b4];

regional_nd_p ={};
regional_nd_stats = {};

for i=1:68
    i
    y_hc = hc_data(:,i);
    
    [p1,stats1]=quantreg(x_hc,y_hc,0.95,[],1000);
    [p2,stats2]=quantreg(x_hc,y_hc,0.5,[],1000);
    [p3,stats3]=quantreg(x_hc,y_hc,0.05,[],1000);
    
    regional_nd_p{1,i} = p1;
    regional_nd_p{2,i} = p2;
    regional_nd_p{3,i} = p3;
    
    regional_nd_stats{1,i} = stats1;
    regional_nd_stats{2,i} = stats2;
    regional_nd_stats{3,i} = stats3;
    
end




