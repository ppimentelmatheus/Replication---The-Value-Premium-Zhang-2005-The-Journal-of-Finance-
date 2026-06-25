
function lfv = TimeAggr(hfv, basef, targetf, sorf)

%-------------------------------------------------------------------------------
% TIMEAGGR: Perform time aggregation from high frequency to low frequency 
% 
% USAGE:
%   hfv     - matrix of variable with high frequency to be transformed
%   basef   - original high frequency:
%      'd' - daily  
%      'w' - weekly   
%      'm' - monthly
%      'q' - quarterly
%   targetf - target low frequency of interest: 
%      'm'  - monthly   
%      'q'  - quarterly   
%      'y'  - yearly    
%      '2y' - bi-annual
%      '4y' - 4-yearly
%   sorf    - whether hfv is of stock or flow variable
%      's'  - stock     
%      'f'  - flow
%
% RESULT:
%   lfv     - matrix of low frequency variable
%
% Example: Vf = TimeAggr(V_f, 'w', 'y', 's')
%          Rf = TimeAggr(R_f, 'w', 'm', 'f')
% 
% NOTE: Always use continuous-compound return!
%-----------------------------------------------------------------------------

[N, T] = size(hfv);
if basef == 'd'
   bfnum = 360; 
   if targetf == 'm'
      tfnum = 30;
   elseif targetf == 'q'
      tfnum = 90;
   elseif targetf == 'y'
      tfnum = 360;
   end
elseif basef == 'w'
   bfnum = 48;
   if targetf == 'm'
      tfnum = 4;
   elseif targetf == 'q'
      tfnum = 12;
   elseif targetf == 'y'
      tfnum = 48;
   elseif targetf == '2y'
      tfnum = 96;
   elseif targetf == '4y'
      tfnum = 96 * 2;
   end   
elseif basef == 'm'
   bfnum = 12;
   if targetf == 'm'
      tfnum = 1;
   elseif targetf == 'q'
      tfnum = 4;
   elseif targetf == 'y'
      tfnum = 12;
   elseif targetf == '2y'
      tfnum = 24;
   elseif targetf == '4y'
      tfnum = 48;
   end   
elseif basef == 'q'
   bfnum = 4;
   if targetf == 'm'
       error('Attempt to got from low frequency to high frequency data.')
   elseif targetf == 'q'
       lfv = hfv;
       return
   elseif targetf == 'y'
       tfnum = 4;
   elseif targetf == '2y'
       tfnum = 8;
   elseif targetf == '4y'
       tfnum = 16;
   end
else  
   error('Your input of high frequency is not supported by this routine!')
end

% time aggregation of a stock variable - take the beginning-of-period observation as that for the low-freq period observ
lfv = zeros(N, T/tfnum);
if sorf == 's'   
   col = zeros(1, T/tfnum + 1); 
   col(1) = 1;
   for j = 1 : T/tfnum
      col( j + 1)  = tfnum * j ;
   end
   lfv = hfv( :, col(2 : T/tfnum + 1) );
% time aggregation of a flow variable such as R_f 
elseif sorf == 'f'  
   mtmp = reshape (hfv, N, tfnum, T/tfnum);
   mtmp = sum(mtmp, 2);  
   lfv  = squeeze(mtmp);   
   % a small trick to ensure the column-wise dimension is time since when N = 1 squeeze command operates on N dimension as well.
   if N == 1, lfv = lfv'; end
end
