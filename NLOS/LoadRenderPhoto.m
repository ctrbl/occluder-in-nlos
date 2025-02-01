function [test_im_crop_ds] = LoadRenderPhoto(fname,downsamp_factor,numPixels)

% test_im = imread([fname]);

test_im = load([fname]);
test_im = test_im.(fname);
test_im_crop = double(test_im(1:numPixels(2),1:numPixels(1),:));


test_im_crop_ds = cat(3, downsample_2(test_im_crop(:,:,1),downsamp_factor),...
                    downsample_2(test_im_crop(:,:,2),downsamp_factor), ...
                    downsample_2(test_im_crop(:,:,3),downsamp_factor));

end


function [image] = downsample_2(im,iterations)
    %Downsampling 2^iterations times

    for i=1:iterations
        im1 = im(1:2:end,1:2:end);
        im2 = im(1:2:end,2:2:end);
        im3 = im(2:2:end,1:2:end);
        im4 = im(2:2:end,2:2:end);

        image = (im1+im2+im3+im4)/4;
        im = image;
    end
    
    if iterations == 0 
        image = im;
    end
end