function s=dec2bin_array(d,numBits)
%DEC2BIN Convert decimal integer to its binary representation
%   DEC2BIN(D) returns the binary representation of D as a character
%   vector. D must be an integer. If D is greater than flintmax, DEC2BIN
%   might not return an exact representation of D.
%
%   DEC2BIN(D,numBits) produces a binary representation with at least
%   numBits bits.
%
%   Example
%      dec2bin(23) returns '10111'
%
%   See also BIN2DEC, DEC2HEX, DEC2BASE, FLINTMAX.

%   Copyright 1984-2019 The MathWorks, Inc.
% Modified @ 2022-09-25

narginchk(1,2);

if isempty(d)
    s = [];
    return;
end

if ~(isnumeric(d) || islogical(d) || ischar(d))
    error(message('MATLAB:dec2bin:InvalidDecimalArg'));
elseif ~isreal(d)
    error(message('MATLAB:dec2bin:MustBeReal'));
end

d = d(:); % Make sure d is a column vector.

if ~all(isfinite(d))
    error(message('MATLAB:dec2bin:MustBeFinite'));
end

if nargin<2
    numBits=1; % Need at least one digit even for 0.
else
    if ~(isnumeric(numBits) || ischar(numBits)) || ~isscalar(numBits) ...
            || numBits<0 || ~isfinite(numBits) || ~isreal(numBits)
        error(message('MATLAB:dec2bin:InvalidBitArg'));
    end
    numBits = double(numBits);
    numBits = round(numBits); % Make sure n is an integer.
end

if any(d<0)
    if any(d<intmin('int64'))
        error(message('MATLAB:dec2bin:NegativeValueOutOfRange'));
    end
    
    d = upcastNegatives(d);
end
    
if isfloat(d) || all(d < flintmax)
    d = double(d);

    [~, e] = log2(max(d)); % How many digits do we need to represent the numbers?
    numBits = max(e, numBits);
    powersOf2 = pow2(1-numBits:0);
    bits = rem(floor(d*powersOf2),2);
    s = bits;
    % s = char(bits + 48); % char(48) is '0'
else
    numOfHexDigitsForValue = getNumberOfBinaryDigits(max(d), numBits);
    numBits = max(numBits, numOfHexDigitsForValue);
    s = char(ones(numel(d),numBits)*48);

    for i=1:numel(d)
        s(i,:) = bigIntDec2bin(d(i),numBits);
    end
end
end

function minDigitsForGivenValue = ...
    getNumberOfBinaryDigits(biggestDec, inputWidth)

    maxInt = intmax('uint64');

    if biggestDec > maxInt
        [~, powersOfTwo] = log2(double(biggestDec));
        minDigitsForGivenValue = max(inputWidth, powersOfTwo);
    elseif double(biggestDec) + eps(double(biggestDec)) > maxInt
        minDigitsForGivenValue = 64;
    else
        biggestDec = double(biggestDec);
        minDigitsForGivenValue = 0;
        while (biggestDec >= 1)
            biggestDec = biggestDec / 2; 
            minDigitsForGivenValue = minDigitsForGivenValue+1; 
        end 
    end
    
end

function b = bigIntDec2bin(n, columnWidth)

    switch class(n)
    case 'uint64'
        two = uint64(2);
    case 'char'
        n = uint64(n);
        two = uint64(2);
    case 'int64' 
        two = int64(2);
    case 'logical'
        n = uint8(n);
        two = 2;
    otherwise
        two = 2;
    end
    
    binaryNum = false(1,columnWidth);
    bit = 1;
    
    while (n > 0)
        binaryNum(bit) = mod(n, 2);
        n = (n - rem(n,two)) ./ two;
        bit = bit+1; 
    end
    
    b = char(ones(1,columnWidth)*'0');
    b(binaryNum) = '1';
    b = flip(b,2);
    
end
