function ImgF = threeD_to_twoD(ImgF, varargin);
%THREED_TO_TWOD convert a 3D datacube to a 2D one
%if validPixels are provided then subselect on those pixels

i_p = inputParser;
i_p.addRequired('ImgF')
i_p.addOptional('validPixels',nan);

i_p.parse(ImgF,  varargin{:});

ImgF = reshape(ImgF, prod([size(ImgF,1),size(ImgF,2)]), size(ImgF,3));

if ~isnan(i_p.Results.validPixels)
    ImgF = ImgF(i_p.Results.validPixels, :);
end
end

