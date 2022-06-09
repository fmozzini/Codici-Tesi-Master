function [data] = openLargeDatafiles(filename)

blockSize = 500e6 ; % Choose large enough so not too many blocks.
 tailSize  = 100 ;   % Choose large enough so larger than one value representation.
 % - Open file.
 fId  = fopen(filename) ;
 % - Read first line, convert to double, determine #columns.
 line    = fgetl( fId ) ;
 data    = sscanf( line, '%f;' ) ;
 nCols   = numel( data ) ;
 lastBit = '' ;
 while ~feof( fId )
     % - Read and pre-process block.
     buffer = fread( fId, blockSize, '*char' ) ;
     isLast = length( buffer ) < blockSize ;
     buffer(buffer==10) = ';' ;
     buffer(buffer==13) = '' ;
     % - Pre-pend last bit of last block.  
     if ~isempty( lastBit )
         buffer = [lastBit; buffer] ; %#ok<AGROW>
     end
     % - Truncate to last ',' and keep last bit for next iteration.
     if ~isLast
         n = find( buffer(end-tailSize:end)==';', 1, 'last' ) ;
         cutAt   = length(buffer) - tailSize + n ;
         lastBit = buffer(cutAt:end) ;
         buffer(cutAt:end) = [] ;
     end
     % - Parse.
     data = [data; sscanf( buffer, '%f;' )] ; %#ok<AGROW>
 end
 % - Reshape data vector -> array.
 data = reshape( data, nCols, [] )' ;
 % - Close file.
 fclose( fId ) ;

 


end
