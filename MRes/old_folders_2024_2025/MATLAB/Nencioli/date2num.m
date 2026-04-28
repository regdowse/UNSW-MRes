function  day_num=date2num(day)  
% day_num .m (call function ) 
% date conversion  
 
% Convert date and time to string forma  
daytime=datestr(day,30);  
% Read the year month day in the date  
yy=str2num(daytime(1:4));  
mm=str2num(daytime(5:6));  
dd=str2num(daytime(7:8));  
% Convert date and time to serial date number  
day_num=datenum(yy,mm,dd);  
 
