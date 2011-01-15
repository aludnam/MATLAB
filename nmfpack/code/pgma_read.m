function gray = pgma_read ( file_name )

% function gray = pgma_read ( file_name )
%
%  PGMA_READ opens an ASCII PGM file and reads the data.
%  Can also read binary data (replace P2 by P5) - which 
%  must come directly after the max value of 255, in the next line.
%
%  The PGMA file is assumed to have the format:
%
%    P2
%    # feep.pgma
%    ncol nrow
%    maxgray, the largest legal value for the data.
%    row 1, a list of ncol numbers between 0 and maxgray.
%    row 2, a list of ncol numbers between 0 and maxgray.
%    ...
%    row nrow, a list of ncol numbers between 0 and maxgray.
%
%    Lines beginning with '#' are comments.
%
%  Example:
%
%    P2
%    # feep.pgm
%    24 7
%    15
%    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
%    0  3  3  3  3  0  0  7  7  7  7  0  0 11 11 11 11  0  0 15 15 15 15  0
%    0  3  0  0  0  0  0  7  0  0  0  0  0 11  0  0  0  0  0 15  0  0 15  0
%    0  3  3  3  0  0  0  7  7  7  0  0  0 11 11 11  0  0  0 15 15 15 15  0
%    0  3  0  0  0  0  0  7  0  0  0  0  0 11  0  0  0  0  0 15  0  0  0  0
%    0  3  0  0  0  0  0  7  7  7  7  0  0 11 11 11 11  0  0 15  0  0  0  0
%    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
%
%  Modified:
%
%    21 February 2004
%
%    1. June 2005 by Fabian Theis (fabian@theis.name) to include binary pgms as well
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, character *FILE_NAME is the name of the PGMA file to read.
%
%    Output, integer GRAY[NROW,NCOL], is the gray scale data read from the file.
%
  FALSE = 0;
  TRUE = 1;

  gray = [];

  fid = fopen ( file_name, 'r' );

  if ( fid < 0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'PGMA_READ - Fatal error!\n' );
    fprintf ( 1, '  Could not open the input file.\n' );
    error ( 'PGMA_READ - Fatal error!' );
    return;
  end
%
%  Read the first line.
%
  line = fgets ( fid );
%
%  Verify that the first two characters are the "magic number".
%  Matlab strncmp returns 1 for equality, and 0 for inequality.
%
binary = 0;
  if ( strncmp ( line, 'P2', 2 ) == 0 )
      if (strncmp ( line, 'P5', 2 ) == 0 )
        error('P2 or P5 missing in header');
      else
          binary = 1;
      end
  end
%
%  Move to the next noncomment line.
%
  while ( 1 )
    line = fgets ( fid );
    if ( line(1) ~= '#' ) 
      break;
    end
  end
%
%  Extract NCOL and NROW.
%
  [ array, count ] = sscanf ( line, '%d' );
  ncol = array(1);
  nrow = array(2);
%
%  Move to the next noncomment line.
%
  while ( 1 )
    line = fgets ( fid );
    if ( line(1) ~= '#' ) 
      break;
    end
  end
%
%  Extract MAXGRAY, and ignore it.
%
  [ array, count ] = sscanf ( line, '%d' );
  maxgray = array(1);
  if (binary & maxgray>255)
      error('binary image - only byte values are accepted');      
  end
  
  % here binary part:
  if (binary)
%
%  Move to the next noncomment line.
%
%    while ( 1 )
%      line = fgets ( fid );
%      if ( line(1) ~= '#' ) 
%        break;
%      end
%    end
%%    fid = fopen ( file_name, 'r' );
%%    gray = fread(fid,[nrow,ncol],'uchar');
%    if (size(line,2)~=ncol*nrow)
%        error(['wrong binary line size! (image: ' file_name ')']);
%    end
%    gray=reshape(double(line),ncol,nrow)';    
    gray = fread(fid,[ncol,nrow],'uchar')';
    fclose(fid);
    return
  end
  
%
%  Set aside space for GRAY.
%
  gray = zeros ( nrow, ncol );

  i = 1;
  j = 0;
  done = FALSE;

  while ( done == FALSE )
%
%  Move to the next noncomment line.
%
    while ( 1 )
      line = fgets ( fid );
      if ( line(1) ~= '#' ) 
        break;
      end
    end
    [ array, count ] = sscanf ( line, '%d' );
%
%  Each value that you read goes into the "next" open entry in GRAY.
%
    for k = 1 : count

      j = j + 1;
      if ( ncol < j )
        j = 1;
        i = i + 1;
      end

      if ( i <= nrow )
        gray(i,j) = array(k);
      end

      if ( i == nrow & j == ncol )
        done = TRUE;
      end

    end
  
  end

  fclose ( fid );
