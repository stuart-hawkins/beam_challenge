% Incident beam object.
%
%  p = beam(theta,w,k) returns a beam object p with wavenumber k, beam
%  width w and direction exp(1i*theta).
%
% Also:
%
%   f = p.evaluate(z) returns the values f of the beam at points z.
%
%   f = u.evaluate(z,mask) returns the values f of the beam at
%   points z for which mask==1 and NaN elsewhere.
%
%   [dx,dy] = u.evaluateGradient(z) returns dx and dy the partial 
%   derivatives of the beam in the x and y directions respectively
%   at the points z.
%
%   [dx,dy] = u.evaluateGradient(z,mask) returns dx and dy the partial 
%   derivatives of the beam in the x and y directions respectively
%   at the points z for which mask==1 and NaN elsewhere.
%
%   cof = p.get_coefficients(x,n) returns the vector cof of regular
%   wavefunction expansion coefficients of the beam field with 
%   wavefunction origin x and order n.
%
% Note: in the functions above vectors in the plane are represented by
% complex numbers.
%
% See also: point_source, incident, plane_wave.
%
% Stuart C. Hawkins - 11 December 2022

classdef beam < incident
   
    properties
        direction
        kwave
        width
    end
    
    methods

        %-----------------------------------------------
        % constructor
        %-----------------------------------------------

        function self = beam(varargin)
            
            if nargin==1
                
                if ~isa(varargin{1},'beam')
                    
                    error('Single argument must be a beam')
                    
                end
                
                % set wavenumber
                self.kwave = varargin{1}.kwave;
                
                % set incident direction
                self.direction = varargin{1}.direction;

                % set beam width
                self.width = varargin{1}.width;

            else
            
                % set wavenumber
                self.kwave = varargin{3};
                
                % set incident direction
                self.direction = exp(1i*varargin{1});

                % set beam width
                self.width = varargin{2};
                
            end
            
        end
        
        %-----------------------------------------------
        % return vector of coefficients for scatterer at
        % centre
        %-----------------------------------------------

        function cof = get_coefficients(self,centre,nmax)
                     
            % get the angle of the incident direction
            incang = angle(self.direction);
            
            n=-nmax:nmax;
            n=n(:);
            
            % work out distance of centre from the midline of the beam by
            % projecting centre onto the line z = mu * self.direction
            d = abs( centre - self.direction * real(self.direction * conj(centre))  );

            if d > self.width/2

                % outside the beam... coefficients are zero
                cof = 0 * n;

            else

                % inside the beam... use plane wave coefficients

                % set coefficients cf (3.13) in Ganesh, Hawkins, Hiptmair
                % IMA Journal of Numerical Analysis doi:10.1093/imanum/drr041
                cof = 1i.^abs(n).*exp(-1i*n*incang);

                % adjust the coefficients if the centre of the scatterer is not
                % the origin
                if centre~=0
                    dp = real(conj(self.direction)*centre);
                    phase = exp(1i*self.kwave*dp);
                    cof = phase * cof;
                end
            
            end

        end

        %-----------------------------------------------
        % evaluate method
        %-----------------------------------------------

        function val = evaluate(self,points,mask)
            
            % intialize return array
            val=zeros(size(points));

            % apply mask if necessary
            if nargin>2
                points=points(mask);
            end

            % work out distance of points from the midline of the beam by
            % projecting points onto the line z = mu * self.direction
            d = abs( points - self.direction * real(self.direction * conj(points)) );

            % compute incident field
            v = (d<=self.width/2) .* exp(1i*self.kwave*real(conj(self.direction)*points));
            
            % insert values into the return array
            if nargin>2
                val(mask)=v;
            else
                val=v;
            end
            
        end
        
        %-----------------------------------------------
        % evaluate gradient
        %
        % This is useful for eg Neumann BCs. Note that we
        % cannot use complex numbers to represent vectors
        % in this case because the components of the vector
        % are, in general, complex.
        %-----------------------------------------------

        function [dx,dy] = evaluateGradient(self,points,mask)
            
            % initialize return values
            dx = zeros(size(points));
            dy = zeros(size(points));

            % apply mask if necessary
            if nargin>2
                points=points(mask);
            end

            % work out distance of points from the midline of the beam by
            % projecting points onto the line z = mu * self.direction
            d = abs( points - self.direction * real(self.direction * conj(points)) );

            % compute incident field
            scalarPart = (d<=self.width/2) .* 1i*self.kwave*exp(1i*self.kwave*real(conj(self.direction)*points));
            vx = real(self.direction) * scalarPart;
            vy = imag(self.direction) * scalarPart;
            
            % insert values into the return array
            if nargin>2
                dx(mask) = vx;
                dy(mask) = vy;
            else
                dx = vx;
                dy = vy;
            end
            
        end        
        
    end % end methods
        
end