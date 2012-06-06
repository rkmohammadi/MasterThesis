% a formula for Gabor filter is implemented according to: 
%J. Sung, S.Y. Bang, S. Choi,  
%A Bayesian network classifier and hierarchical Gabor features for handwritten numeral recognition, 
%Pattern recognition letters, 2006., 
%http://linkinghub.elsevier.com/retrieve/pii/S0167865505001935

classdef GaborKernel
   
  properties(GetAccess = 'public', SetAccess = 'private')
      KernelValues     
  end
    
  methods      

        function this = GaborKernel(matrixSize, scale, orientation, frequency, centerPoint)
        % matrixSize  - Size of matrix (odd values are desirable -> center point)
        % scale       - size of Gaussian envelope (bigger value -> bigger filter)
        % orientation - orientation angle (in radians) 
        % frequency   - circular frequency of  sin part (in spatial domain)
        % centerPoint - center point of filter [x y]
            
          pointX = centerPoint(1); 
          pointY = centerPoint(2);
          
          maxX = matrixSize - pointX;
          maxY = matrixSize - pointY;
          
          if(mod(matrixSize, 2) ==0) %if matrix size is even
              maxX = maxX - 1; %kernel size should be odd
              maxY = maxY - 1;
          end

          [x, y] = meshgrid(-pointX :1: maxX, -pointY :1: maxY);
          x = double(x); y = double(y); 
          
          kernel = GaborKernel.Create(x,y,scale, orientation, frequency); 
		  	
          this.KernelValues = GaborKernel.RemoveDC(kernel); %only real part has DC  
          %this.KernelValues = GaborKernel.ScaleValues(this.KernelValues, -1, 1);
        end
  end
  
  methods(Static)
       function amp = GetAmplitudes(kernelValues)
            realPart = real(kernelValues);
            imagPart = imag(kernelValues);
            
            amp = sqrt( realPart.^2 + imagPart.^2 );
       end

       function phase = GetPhases(kernelValues)
            realPart = real(kernelValues);
            imagPart = imag(kernelValues);
            
            phase = atan(imagPart ./ realPart);
       end 
       
       function realPart = GetRealParts(kernelValues)
           realPart = real(kernelValues);          
       end
        
       function imagPart = GetImagParts(kernelValues)
           imagPart = imag(kernelValues);          
       end
        
  end
  
  methods(Static, Access = 'private') %main method
      
        function kernel = Create(x, y, sigma, theta, omega)
		%sigma = scale
		%theta = orientation
		%omega = frequency
		
           e = exp(1); 
            
           r = 1; %constant (values ~=1 -> eliptical filter)

           R1 = x.*cos(theta) + y.*sin(theta);
           R2 =-x.*sin(theta) + y.*cos(theta);

           expFactor = -1/2 * ( (R1/sigma).^2 + (R2/(r*sigma)).^2  );

           gauss = 1 / ( sqrt(r*pi)*sigma) ;
           gauss =  gauss .* e.^expFactor;

           gaborReal = gauss .* cos(omega*R1);
           gaborImag = gauss .* sin(omega*R1);

           kernel = gaborReal + gaborImag*1i;
           kernel = kernel'; %we want kernel(y,x); y=row; x=column;
        end
  end
  
  methods(Static, Access = 'private') %helper methods
      
        function scaledVals = ScaleValues(matrix, newMin, newMax)
            %matrix = double(matrix); %unmark if matrix type is ~=double

            oldMin = min(matrix(:));
            oldMax = max(matrix(:));

            if oldMax-oldMin == 0
                k=0;
            else
                k = (newMax - newMin) / (oldMax - oldMin);
            end

            scaledVals = k.*(matrix - oldMin) + newMin;
        end
       
        function rez = RemoveDC(matrix)
            
            valueDC = mean( matrix(:) );
            rez = matrix - valueDC;
        end  
  end
  
end