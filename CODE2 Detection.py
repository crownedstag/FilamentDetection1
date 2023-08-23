###Load the packages that may be used.
from scipy import ndimage
import numpy as np
import pandas as pd
from astropy.io import fits
from skimage import measure
from skimage import exposure
import math
from scipy.ndimage import rotate
from scipy.ndimage import binary_closing
import os



###
#####y-coordinate-latitude conversion function for later adding the latitude column
def y_coord_to_solar_latitude(y0,B0, image_center_y=1024, solar_radius=900):

    dy = y0 - image_center_y
    
    # Compute the angle between the y-coordinate and the center of the sun
    angle_from_center = np.arcsin(dy / solar_radius)

    # consideration of the B0 angle
    latitude = np.degrees(angle_from_center) - B0

    return latitude


######Angle correction function
def adjust_orientation(angle):
    if angle < 0:
        return 90 + angle
    else:
        return angle - 90


######Filter area function as commonly used to functionalize it.
####Find regions larger than the threshold area
def image_area_filter(threshold_image, threshold):
    regions_image1 = np.zeros(threshold_image.shape, dtype='uint8')
    label1 = measure.label(threshold_image)
    
    for region in measure.regionprops(label1):
        if region.area >= threshold:
            for coordinates in region.coords: 
                regions_image1[coordinates[0], coordinates[1]] = 255  # 或者使用任何你喜欢的值，如1或255
                
    return regions_image1




##########Main code for filament detection
def process_sun_image(filepath):
    with fits.open(filepath, mode='readonly', ignore_missing_end=True) as hdulist:
        ## Unzip the fz format file to get its header and image
        
        image_data = hdulist[1].data
        header = hdulist[1].header
        #Flip the image for subsequent image processing
        flipped_image_data = np.flipud(image_data)
        #Header information that will be used after collection, e.g. date, B0, etc.
        obs_date = header['DATE-OBS']
        date_only = obs_date.split('T')[0]
        B0 = header['SOLAR-B0']

        # Plot the image

        ######far filed correction

        image = flipped_image_data

       # Smoothing images using high intensity Gaussian filtering
        smoothed_image = ndimage.gaussian_filter(image, sigma=40)

        # Normalize the smoothed image 
        flat_field_image = smoothed_image / np.nanmean(smoothed_image)

         # Flat field corrected with 'flat_field_image'.
        corrected_image = image / flat_field_image



 ###############Check if the image is corrupted or incomplete
        ###############
        ###Image normalization by setting the average brightness of
        ### the center of the solar image to 200

        center_y, center_x = np.array(image.shape) // 2

        y, x = np.ogrid[:image.shape[0], :image.shape[1]]
        dist_from_center = np.sqrt((x - center_x)**2 + (y - center_y)**2)

        mask = dist_from_center <= 200

        average_brightness = image[mask].mean()

        scaling_factor = 200 / average_brightness

        normalized_sunimage = image * scaling_factor
        
        normalized_sunimage_ref = normalized_sunimage.copy()




        # Discard data with brightness below 50 
        ##(usually considered to be outside the solar disk)
        
        normalized_sunimage[normalized_sunimage < 50] = np.nan
        
        ##Create empty dataset df
        
        df = pd.DataFrame(columns=[    "length", 
                                            "area",
                                            "orientation", 
                                            "centroid_x", 
                                            "centroid_y", 
                                            "Latitude", 
                                            "exist",
                                            "date"])
                

        # Create a circular mask with a radius of 900, 
        ##centered on the center of the image
        
        radius = 900
        normalized_sunimage_ref[dist_from_center >= radius] = np.nan
        # Create binary versions of the images
        image1_binary = np.isfinite(normalized_sunimage)
        image2_binary = np.isfinite(normalized_sunimage_ref)

        # Calculate the intersection 
        intersection = np.sum(image1_binary & image2_binary)

        # Calculate the union 
        union = np.sum(image1_binary | image2_binary)

        # Calculate the IoU (Jaccard index)
        intersection_rate = intersection / union
        
        #If the overlap is less than 98.5%,
        #the code ends and is labeled as an incomplete solar image.
        if intersection_rate < 0.985 or math.isnan(intersection_rate):
                   print("Incomplete sun image")
                   df = df.append({"length": 0,
                                    "area":0,
                                    "orientation":0, 
                                   "centroid_x": 0, 
                                   "centroid_y": 0,
                                   "Latitude": 0,
                                   "exist": 2},ignore_index=True)
                   

        else:
            #######If greater, then continue
            print("complete sun image")

            #Remove regions in the corrected image that lie outside the sun's disk
            corrected_image1 = corrected_image.copy()
            corrected_image1[dist_from_center >= radius] = np.nan


           ###############
           ####Ensure that the median image is 200 to
           ###facilitate the processing of large numbers of images

            median = np.nanmedian(corrected_image1)

            scale_factor = 200 / median

            new_image = corrected_image1 * scale_factor

            new_median = np.nanmedian(new_image)


            image_ref = new_image.copy()



            # # Sift out regions suspected to be sunspots (brightness below 140)
            image_ref[image_ref < 140] = np.nan



            ######Perform histogram equalization

            mask = np.isnan(image_ref) 

            image_no_nan_no_zero = image_ref[~mask]

            equalized_img_no_nan_no_zero = exposure.equalize_hist(image_no_nan_no_zero)

            equalized_img = np.copy(image_ref)
            equalized_img[~mask] = equalized_img_no_nan_no_zero
            
            
            
            
            
            #######Apply a global threshold of 0.03 of the median of the equalized image
            new_median = np.nanmedian(equalized_img)

            threshold_image = equalized_img < new_median*0.03
            
            
            
            ####First area threshold to screen out small area noise
            regions_image1 = image_area_filter(threshold_image, 50)
            
            binary_image = regions_image1
            
            
            
            
            ####Create linear structure elements 
            ####and use them to process images with closed operations

            selem_orig = np.zeros((41, 41), dtype=np.uint8)
            selem_orig[20, :] = 1

            #  Define degrees
            angles = [0, 45, 90, 135, 67.5, 112.5, 22.5, 157.5]

            result_image = np.zeros_like(binary_image)

            # for each degree
            for i, angle in enumerate(angles):
                selem_rotated_large = rotate(selem_orig, angle, reshape=False)
                selem_rotated = selem_rotated_large[9:31, 9:31]
                
            
                ##Performing closing
                closed_image = binary_closing(binary_image, selem_rotated)
                

                ###Combining Images
                result_image = np.logical_or(result_image, closed_image)
                
                
                
            ###Second area thresholding process, retaining only large areas   
            filterd_image = image_area_filter(result_image, 150)


              #Perform the closing operation again
            closed_image1 = binary_closing(filterd_image, structure=np.ones((6, 6)))
            closed_image1[dist_from_center >= 895] = 0

            
                
                      
               #### # Processing of each region
            labels2 = measure.label(closed_image1)
            filament_image = np.zeros_like(closed_image1)
            props = measure.regionprops(labels2)
            
            
            #Use the flag variable to keep track of whether 
            #a region exists that meets the conditions.
            has_valid_region = False




           # Selection of regions based on conditions, illustrated in Section 2.2.3 in the disseration 
            for region in measure.regionprops(labels2):
               if region.eccentricity > 0.920 and region.minor_axis_length <= 200 and region.major_axis_length > 50 and region.area >= 500:
                  has_valid_region = True  # Region that meets the conditions
                  
                 
                  orientation_deg = np.rad2deg(region.orientation)
                  ellipse_major_axis = region.major_axis_length
                  ellipse_minor_axis = region.minor_axis_length

                  ellipse_orientation = np.rad2deg(region.orientation) 
                  y0, x0 = region.centroid
                  df = df.append({
                 "length": ellipse_major_axis,
                 "area": region.area,
                 "orientation": adjust_orientation(ellipse_orientation),
                 "centroid_x": x0,
                 "centroid_y": 2048-y0,
                 "Latitude": y_coord_to_solar_latitude(2048-y0, B0), 
                 "exist": 1}, ignore_index=True)

                    # Images that do not satisfy the condition, i.e. images without filaments
            if not has_valid_region:
                  print("No regions detected that satisfy the conditions.")
                  ##If no filament is detected, mark index as 0 
                  df = df.append({
                  "length": 0,
                  "area": 0,
                  "orientation": 0,
                  "centroid_x": 0,
                  "centroid_y": 0,
                  "Latitude": 0,
                  "exist": 0}, ignore_index=True)
            

            

         ###Assign a value to the date column based on the current date.
        df['date'] = date_only
            
        df = df.round(1)
        cols = df.columns.tolist()
        cols = [cols[-1]] + cols[:-1]
        df = df[cols]
            
            
    return df
                              




###Functions for processing sun image files and generating csv files
def process_directory(directory, output_directory, year_to_process):
      # Create output_directory if it doesn't exist
      os.makedirs(output_directory, exist_ok=True)

      for root, dirs, files in os.walk(directory):
          # Only proceed if root has a format like {directory}/{year}/{month}
          path_parts = os.path.normpath(root).split(os.sep)
          if len(path_parts) >= 3 and path_parts[-2] == str(year_to_process):
              year = path_parts[-2]
              month = path_parts[-1]

              dfs = []  # list to store dataframes
              for file in files:
                  filepath = os.path.join(root, file)
                  df = process_sun_image(filepath)
                  dfs.append(df)

              # If there were any .fts files in this directory, save them to a CSV
              if dfs:
                  result = pd.concat(dfs)
                  output_filename = os.path.join(output_directory, "{}_{}.csv".format(year, month))
                  result.to_csv(output_filename)
#####This function will generate a csv file with basic characterization information for 
#####each detected filament for that time period in year_month format   
##Next, by changing the value of year_to_process, 
#####can download a csv file of all filaments for each month from 2011-2023     
