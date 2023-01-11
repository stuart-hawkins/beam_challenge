% Sensor class for the INI MSW beam challenge
%
%   s = sensor(z,r) creates a circular sensor with radius r at z. Here 
%   z = x + 1i * y  is a complex number representing the point (x,y).
%
%   s.plot() plots the sensor.
%
%   s.plot(linetype) plots the sensor with specified linetype eg 'k-'
%
%   v = s.measure(u) measures the field u over the sensor. Here u should be
%   of class wavefunctionexpansion or incident.
%
%  See [1] for description of the sensor.
%
% References:
%
%   [1] https://gist.github.com/arturgower/0251e3d5696bf2f49dc8e311a5492bab
%
% Stuart C. Hawkins - 11 January 2023

classdef sensor < handle

    properties
        position
        radius
    end

    methods
    
        %----------------------------------
        % constructor
        %----------------------------------

        function self = sensor(x,r)

            % store position and radius
            self.position = x;
            self.radius = r;

        end

        %----------------------------------
        % plot
        %----------------------------------

        function varargout = plot(self,opts)

            % set default linetype
            if nargin<2
                opts = 'k-';
            end

            % create array of angles for plotting the circle
            tp = linspace(0,2*pi,100);

            % plot the circle
            h = plot(real(self.position)+self.radius*cos(tp),...
                imag(self.position)+self.radius*sin(tp),opts);

            % return figure handle if required
            if nargout>0
                varargout{1} = h;
            end

        end

        %----------------------------------
        % measure the field
        %----------------------------------

        % Note: we assume varargin (a cell array) holds
        % the components of the field

        function val = measure(self,varargin)

            % set resolution for grid inside the circle
            res = 20;

            % setup square grid of points with sides of lengt 2*self.radius
            [x,y] = meshgrid(linspace(-self.radius,self.radius,res));
            
            % get the points in complex format and centre it on
            % self.position
            z = self.position + x(:) + 1i*y(:);

            % throw away points outside the circular sensor
            z = z(find(abs(z-self.position)<self.radius));

            % evaluate the field from the first component
            u = varargin{1}.evaluate(z);

            % loop through the other components and add on
            % their contribution to the field
            for j=2:length(varargin)   
                u = u + varargin{j}.evaluate(z);
            end

            % compute the mean of the modulus of the field
            % over the points
            val = mean(abs(u));

        end

    end

end