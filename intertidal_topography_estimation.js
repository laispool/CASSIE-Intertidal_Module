//////////////////////////////////////////////////
//  Extract the geometry of the intertidal zone //
//  Find waterline and assign water level       //
//////////////////////////////////////////////////

//////////////////////////////////////////////////
///// Laís Pool  (lais.pool@gmail.com)       /////                        
///// Florianopolis, 24/12/2024              /////              
///// Code Editor - Earth Engine             /////                             
//////////////////////////////////////////////////

///// SOME DEFINITIONS
var img0 = ee.Image('COPERNICUS/S2_HARMONIZED/20230718T132239_20230718T132238_T23MNT');
var img1 = ee.Image('COPERNICUS/S2_HARMONIZED/20230723T132241_20230723T132237_T23MNT');
var img2 = ee.Image('COPERNICUS/S2_HARMONIZED/20230812T132241_20230812T132237_T23MNT');
var img3 = ee.Image('COPERNICUS/S2_HARMONIZED/20230827T132239_20230827T132236_T23MNT');
var img4 = ee.Image('COPERNICUS/S2_HARMONIZED/20230916T132239_20230916T132235_T23MNT');
var img5 = ee.Image('COPERNICUS/S2_HARMONIZED/20231031T132231_20231031T132235_T23MNT');
var img6 = ee.Image('COPERNICUS/S2_HARMONIZED/20231120T132231_20231120T132229_T23MNT');
var img7 = ee.Image('COPERNICUS/S2_HARMONIZED/20231125T132229_20231125T132228_T23MNT');

var imageCollection = ee.ImageCollection.fromImages([img0,img1, img2, img3, img4,img5,img6,img7]);  

var geometry = ee.Geometry.Polygon(
    [[[-44.35507322619092, -2.5298727414928885],
      [-44.35657526323926, -2.5323165296313523],
      [-44.3555023796333, -2.5330882512438886],
      [-44.351639998651855, -2.5324451499320286],
      [-44.34906507799756, -2.531630554478721],
      [-44.347391379572265, -2.530044235547647],
      [-44.348678839899414, -2.5282006732707076],
      [-44.35018087694775, -2.5283292939796906]]]);
Map.centerObject(geometry, 15);

var scale = 10;
var index = 'awei'; // ndwi, mndwi_sharp, awei
var thLand = -0.1;
var stdDev = 'awei_stdDev'

var vis_rgb = {max: 0.25,bands: ['B4', 'B3', 'B2']};
var singleBandVis = {'min': -0.5,'max': 1};
var aweiVis = {'max': 0.2,'min': -0.2};
var vis_intertidal_area = {'min': -0.5,'max': 1,'dimensions': 500,'palette': ['white','blue']};
              
///// SCALE FACTOR
// Function to apply a scale factor to optical bands in the image.
// Sentinel-2 bands are scaled by a factor of 0.0001 to convert digital numbers to reflectance values.
// The adjusted bands are added back to the image, replacing the original bands.
function reflec_corr (image){
  var opticalBands = image.select('B.*').multiply(0.0001); // Apply scale factor to optical bands.
    return image.addBands(opticalBands, null, true); // Add scaled bands back to the image, replacing originals.
}
  

///// CALCULATE INDICES
// Function to calculate various water indices and add them as new bands to the image.
// These indices are used for water body detection and analysis.
var indexFunction = function(image) { 
  // NDWI (Normalized Difference Water Index) using Green (B3) and NIR (B8) bands. (McFeeters, 1996)
 var ndwi = image.normalizedDifference(['B3', 'B8']).rename ('ndwi'); //Mc Feeters, 1996
 // MNDWI (Modified Normalized Difference Water Index) using Green (B3) and SWIR1 (B11) bands. (Xu, 2006)
 var mndwi = image.normalizedDifference(['B3', 'B11']).rename ('mndwi'); //Xu, 2006
 // AWEI (Automated Water Extraction Index) calculation using a custom expression. (Feyisa et al., 2014)
 var awei = image.expression('(B+ (2.5*G) -1.5*(N+S1) -(0.25*S2)) ',{ //Feyisa etal, 2014 
   B: image.select('B2'),
   G: image.select('B3'), 
   S1: image.select('B11'), 
   S2: image.select('B12'),
   N: image.select('B8'),
  }).rename('awei');
   
  return image.addBands([ndwi,mndwi,awei]); // Add the calculated indices as new bands to the image.
};

///// CLIP BY GEOMETRY
// Function to clip an image to the specified geometry.
// This restricts the analysis to the area of interest.
var clip_image = function (image){
  return image.clip(geometry); // Clip the image to the 'geometry' region.
};

///// MASK LAND WITH MNDWI
// Function to mask out land areas using the MNDWI_sharp index.
// Pixels with MNDWI_sharp values greater than or equal to the threshold (thLand) are retained.

var mask_land = function (image){
  var ndwi = image.select('mndwi'); // Select the MNDWI_sharp band.
  return image.updateMask(ndwi.gte(thLand)); // Apply the mask, keeping only water pixels.
};

///// OTSU THRESHOLDING (Otsu, 1979)
// Function to calculate the optimal threshold value for separating two classes (e.g., water and land) 
// based on the Otsu method. This method maximizes the between-class variance of pixel values in a histogram.

var otsu = function(histogram) {
  // Extract the histogram counts and bucket means (bin centers) from the input.
  var counts = ee.Array(ee.Dictionary(index_stdDev).get('histogram'));
  var means = ee.Array(ee.Dictionary(index_stdDev).get('bucketMeans'));
  
  // Compute the total number of bins in the histogram.
  var size = means.length().get([0]);
  
  // Calculate the total sum of counts and the weighted sum of pixel values.
  var total = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
  var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
  
  // Compute the global mean (mean pixel value across all bins).
  var mean = sum.divide(total);
  
  // Create a list of indices from 1 to the size of the histogram.
  var indices = ee.List.sequence(1, size);
  
  // Calculate between-class variance (BSS) for each possible threshold.
  var bss = indices.map(function(i) {
    // Slice the histogram to obtain the counts and means for the first (a) and second (b) classes.
    var aCounts = counts.slice(0, 0, i); // Class a: all bins up to the current threshold.
    var aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]); // Total count in class a.
    var aMeans = means.slice(0, 0, i);// Means in class a.
    var aMean = aMeans.multiply(aCounts)
        .reduce(ee.Reducer.sum(), [0]).get([0]).divide(aCount); // Mean pixel value in class a.
        
    var bCount = total.subtract(aCount); // Total count in class b.
    var bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount);// Mean pixel value in class b.
    
    // Compute between-class variance (BSS) for the current threshold.
    return aCount.multiply(aMean.subtract(mean).pow(2)).add(
          bCount.multiply(bMean.subtract(mean).pow(2)));
  });
  
  // Visualize the between-class variance as a chart for analysis.
  //print(ui.Chart.array.values(ee.Array(bss), 0, means));
  
  // Return the threshold corresponding to the maximum between-class variance. 
  return means.sort(bss).get([-1]);
};

///// APPLY PRE PROCESSING FUNCTIONS
imageCollection = imageCollection.map(reflec_corr);

var NWI = imageCollection.map(indexFunction);
var NWINoMask=NWI.map(clip_image); 
NWI=NWINoMask.map(mask_land);  

///// BEGIN OF CODE ///// 

// Create the Standard Deviation (STD) index image using all the images on the collection
// choosing the index band
var NWI_STD = NWI.select(index).reduce(ee.Reducer.stdDev());

// Create the histogram through the standard deviation index image
var histogram = NWI_STD.reduceRegion({
  reducer: ee.Reducer.histogram(),
  geometry: geometry, 
  scale: 10
});

var index_stdDev = histogram.get(stdDev); // Select the wanted histogram values

var threshold = otsu(histogram.get('histogram')); // Apply Otsu's thresholding

// Mask pixels in the NWI_STD (Normalized Water Index Standard Deviation) image that are below the threshold.
// Only pixels with values greater than or equal to the threshold are retained in the resulting masked image.
var stdMasked = NWI_STD.updateMask(NWI_STD.gte(threshold));

// Create a binary image (zones) where pixels greater than or equal to the threshold are assigned a value of 1 (true),
// and all other pixels are 0 (false).
var zones = NWI_STD.gte(threshold);

// Update the mask for the binary zones image, ensuring that only non-zero (true) pixels are displayed.
// This removes any remaining zero-value pixels from visualization or further processing.
var zones = zones.updateMask(zones.neq(0));

// Convert the raster `zones` image into vector polygons with additional properties.
// This transformation enables further spatial analysis and visualization of the masked regions.
var vectors = zones.addBands(NWI_STD).reduceToVectors({ 
  crs: 'EPSG:4326',
  scale: 10,
  geometryType:'polygon',
  labelProperty: 'stdMasked',
  eightConnected: false,
  geometry: geometry,
  maxPixels: 100e9,
  geometryInNativeProjection: true,
  reducer: ee.Reducer.mean()
});

var intertidal_zones = NWI.map(function(image) {return image.clip(vectors)}); // Clip all images on the collection by the intertidal zone vector

Map.addLayer(intertidal_zones.select(index), vis_intertidal_area, 'intertidal zone '+index); // Add the intertidal zone on map

///////////////// END FIRST CODE

///// COMBINE DATE AND TIME FROM TIDE FILE
// Function to combine the 'Date' (date) and 'Time' (time) properties from a feature into a single property.
// This creates a unified 'dateTime' property as a string for easier handling of temporal data.
var timeParser = function(feature) {
  // Retrieve 'Date' and 'Time' properties
  var dateStr = ee.String(feature.get('Date')); // e.g., "01/01/2023"
  var timeStr = ee.String(feature.get('Time')); // e.g., "00:00:00"

  // Combine 'Date' and 'Time' into a single string
  var dateTimeStr = dateStr.cat(' ').cat(timeStr); // e.g., "01/01/2023 00:00:00"

  // Parse the combined string into an Earth Engine Date object
  var dateTime = ee.Date.parse('dd/MM/yyyy HH:mm:ss', dateTimeStr);

  // Return the feature with the new 'dateTime' property
  return feature.set('dateTime', dateTime);
};

///// OTSU THRESHOLDING (Otsu, 1979)
// Function to calculate the optimal threshold value for separating two classes (e.g., water and land) 
// based on the Otsu method. This method maximizes the between-class variance of pixel values in a histogram.
var otsu = function(histogram) {
  // Extract the histogram counts and bucket means (bin centers) from the input.
  var counts = ee.Array(ee.Dictionary(index).get('histogram')); //RE DO // Frequency of pixel values in each bin.
  var means = ee.Array(ee.Dictionary(index).get('bucketMeans')); //RE DO // The center values of each bin.
  
  var size = means.length().get([0]);// Compute the total number of bins in the histogram.
  
  // Calculate the total sum of counts and the weighted sum of pixel values.
  var total = counts.reduce(ee.Reducer.sum(), [0]).get([0]); // Total number of pixels.
  var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]); // Weighted sum of pixel values.
  var mean = sum.divide(total); // Compute the global mean (mean pixel value across all bins).

  var indices = ee.List.sequence(1, size); // Create a list of indices from 1 to the size of the histogram.
  
  // Calculate between-class variance (BSS) for each possible threshold.
  var bss = indices.map(function(i) {
    // Slice the histogram to obtain the counts and means for the first (a) and second (b) classes.
    var aCounts = counts.slice(0, 0, i); // Class a: all bins up to the current threshold.
    var aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]); // Total count in class a.
    var aMeans = means.slice(0, 0, i);// Means in class a.
    var aMean = aMeans.multiply(aCounts)
        .reduce(ee.Reducer.sum(), [0]).get([0]).divide(aCount); // Mean pixel value in class a.
        
    var bCount = total.subtract(aCount); // Total count in class b.
    var bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount); // Mean pixel value in class b.
    
    // Compute between-class variance (BSS) for the current threshold.
    return aCount.multiply(aMean.subtract(mean).pow(2)).add(
          bCount.multiply(bMean.subtract(mean).pow(2)));
  });
  
  // Visualize the between-class variance as a chart for analysis.
  //print(ui.Chart.array.values(ee.Array(bss), 0, means));
  
  // Return the threshold corresponding to the maximum between-class variance.
  return means.sort(bss).get([-1]);
};

///// ADD DATETIME PROPERTY and CORRECT TIMEZONE and CACULATE TIME DIFFERENCE
// Function to add datetime property to each image and correct the timezone
// Function to add datetime property to each image and correct the timezone
var addDatetimeProperty = function(image) {
  var datetime = ee.Date(image.get('system:time_start')); // Retrieve the system start time as a date object.
  var updatedDatetime = datetime.advance(-3, 'hour'); // Adjust the datetime by subtracting 3 hours (correct for timezone).

  // Return the image with the added datetime property.
  return image.set('datetime', updatedDatetime.format('YYYY-MM-dd\'T\'HH:mm:ss')); // ISO 8601 format
};

// Combined function to calculate time difference and set properties
var calculateTimeDifferenceAndSetProperties = function(image, feature) {
  var imageDate = ee.Date(image.get('datetime'));
  var featureDate = ee.Date(feature.get('dateTime')); // Retrieve the 'dateTime' property from the feature.
  var timeDifference = (featureDate.difference(imageDate, 'second')).abs(); // Compute the absolute time difference in seconds.

  return feature.set({
    'timeDifference': timeDifference,
    'dateImage': imageDate,
    'dateGauge': featureDate
  });
};

///// ADD WATER LEVEL AND DATE TIME
// This function adds the water level ('Tide') and formatted date properties 
// ('dateImage' and 'dateGauge') from the closest image feature to the given feature.
var addTideValue = function(feature) {
  var tideValue = ee.Number(closestFeatureImage.get('Value')).format('%.2f'); // Extract the tide value property.
  var dateImage = closestFeatureImage.get('dateImage'); // Get the image date property.
  dateImage = ee.Date(dateImage).format('dd/MM/yyyy HH:mm:ss'); // Format the date to 'yyyy-MM-dd'.
  var dateGauge = closestFeatureImage.get('dateGauge'); // Get the tide gauge date property.
  dateGauge = ee.Date(dateGauge).format('dd/MM/yyyy HH:mm:ss'); // Format the date to 'yyyy-MM-dd'.
  
  // Add the 'Tide' value and formatted dates as new properties to the feature.
  return feature.set({
    'Tide': tideValue,
    'dateImage': dateImage,
    'dateGauge': dateGauge
  });
};

///// IDENTIFY WATERLINE POSITION
// This function identifies the position of the waterline based on a specified image band 
// and thresholding using the Otsu method. It outputs the geometry of the intertidal feature.
var identifyWaterFeature = function (image,geometry,scale,band){

  var internalBandName = 'ndwi'; // Placeholder for the band to be processed.
  // Compute the histogram of the selected band over the given geometry.
  var histogram = image.reduceRegion({
    reducer: ee.Reducer.histogram(),
    geometry: geometry, 
    scale: 10
  });
  var histNDWI = histogram.get(internalBandName); // Retrieve the histogram for the NDWI band.

  var threshold = otsu(histogram.get('histogram'));// Compute the threshold using the histogram.
  
  // Identify intertidal features by thresholding and vectorizing.
  var intertidalFeature = ee.FeatureCollection(image.clip(geometry) // Clip the NDWI image to the area of interest.
        .lte(threshold) // Apply the threshold to create a binary mask.
        .reduceToVectors({ scale: scale, maxPixels: 1e12, eightConnected: false}) // Convert the binary mask to vector features.
        .filter(ee.Filter.eq("label", 1))); // Filter to retain only the waterline features.

  var intertidalFeatureList = intertidalFeature.toList(intertidalFeature.size()); // Convert the feature collection to a list for further processing.
  
  // Handle cases where no intertidal feature is detected.
  if(intertidalFeature.getInfo() === null ){
    return null; // Return null if no feature is detected.
  }else{
    return intertidalFeature.geometry(); // Return the geometry of the detected feature.
  }
};

///// GAUSSIAN KERNEL CURVE
// This function generates a Gaussian kernel curve of a specified size, mean, and standard deviation (sigma).
// The resulting kernel is normalized so that the sum of all values equals 1.
var gaussianKernel = function (size, mean, sigma) {
  // Sub-function to calculate a single value of the Gaussian curve for a given x, mean, and sigma.
  var gaussianCurve = function (x, mean, sigma)  {
    var divider = ee.Number(sigma) // Divider term (normalizing factor for the Gaussian function)
                    .multiply(ee.Number(2).multiply(Math.PI).sqrt());
      
    var exponent = ee.Number(-1).multiply( // Exponent term (Gaussian function formula)
      ee.Number(x)
        .subtract(mean)
        .pow(2)
        .divide(ee.Number(2).multiply(ee.Number(sigma).pow(2)))
    );

    return ee.Number(1).divide(divider).multiply(exponent.exp()); // Gaussian value for the given x
  };
  
  // Calculate the range of x values for the kernel (centered around 0).
  var half = ee.Number(size).divide(2).floor();
  var begin = ee.Number(0).subtract(half),
    end = ee.Number(0).add(half);

  // Get the normal distribution Y value for each X in the interval
  var kernel = ee.List.sequence(begin, end).map(function (i) { 
    return gaussianCurve(i, mean, sigma); // Apply the Gaussian function to each x
  });
  
  var sum = kernel.reduce(ee.Reducer.sum()); // Calculate the sum of all Gaussian values to normalize the kernel.
  
  // Normalize each value, so that the sum of the list will be equal to one
  var normalizedKernel = kernel.map(function (val) { return ee.Number(val).divide(sum)});

  return normalizedKernel; // Return the normalized Gaussian kernel.
};

///// LINEAR GAUSSIAN FILTER
// This function applies a linear Gaussian filter to a list of 2D coordinates, smoothing the path 
// using a Gaussian kernel. The output is a smoothed version of the input coordinates.
var linearGaussianFilter = function (coordinates) {
  // Parameters for the Gaussian kernel
  var samples = 3; // Number of samples in the kernel
  var mean = 0; // Mean for the Gaussian distribution
  var sd = 1; // Standard deviation for the Gaussian distribution
  
  var coordinateList = ee.List(coordinates); // Convert the input coordinates to an Earth Engine list
  // Create the Gaussian kernel
  var kernelSize = ee.Number(samples);
  var kernelMean = ee.Number(mean);
  var kernelSd = ee.Number(sd);
  var kernel = gaussianKernel(kernelSize, kernelMean, kernelSd);
  // Extract the first and last points of the coordinate list
  var first = coordinateList.reduce(ee.Reducer.first()),
    last = coordinateList.reduce(ee.Reducer.last());
  // Generate a sequence to iterate through the coordinate list for kernel application
  var sequence = ee.List.sequence(
    ee.Number(0),
    coordinateList.length().subtract(kernelSize)
  );
  
  // Apply Gaussian smoothing to each segment of the coordinate list
  var path = sequence.map(function (index) {
    // Take interval of the kernel size to apply the smoothing
    // and zip it to the kernel, so each element in the new list
    // will be a pair of a 2d point and its weight
    var interval = coordinateList
      .slice(ee.Number(index), ee.Number(index).add(kernelSize))
      .zip(kernel);

    // Map the elements, multiplying their axis values by their weight
    var gaussian = interval.map(function (element) {
      // Each element contains a 2d point (0) and a kernel weight (1)
      var asList = ee.List(element); // Element contains [point, weight]
      var point = ee.List(asList.get(0)); // Extract the 2D point
      var weight = ee.Number(asList.get(1)); // Extract the kernel weight

       // Multiply each coordinate value (longitude and latitude) by the weight
      return point.map( function (value) { return ee.Number(value).multiply(weight)});
    });

    // Sum the weighted longitude and latitude separately
    var smoothenLong = gaussian
      .map(function (point){ return  ee.List(point).get(0)})
      .reduce(ee.Reducer.sum());
    var smoothenLat = gaussian
      .map(function (point) { return ee.List(point).get(1)})
      .reduce(ee.Reducer.sum());

    // Return the smoothed 2D point
    return ee.List([smoothenLong, smoothenLat]);
  });
  
  // Combine the smoothed path with the first and last points to retain the original endpoints
  var smoothen = ee.List([]).add(first).cat(path).add(last);

  // return original coordinates if the kernelSize is less than or equal to the length
  // of the given coordinates, otherwise return smoothen coordinates.
  return ee.Algorithms.If(
    coordinateList.size().lte(kernelSize),
    coordinateList, // Return original coordinates if kernel size is too large
    smoothen // Return smoothed coordinates otherwise
  );
};

///// APPLY TIME DIFFERENCE FUNCTION
// Map over the tide gauge collection and calculate the time difference
// between each tide gauge feature and the waterline image.
function getClosestFeatureImageAndValue(waterlineImage, tideGauge) {
  var closestCollection = tideGauge.map(function(feature) {
    return calculateTimeDifferenceAndSetProperties(waterlineImage, feature);
  });
  // Sort the tide gauge collection by the calculated time difference in ascending order,
  // so the closest feature appears first.
  var sortedCollection = closestCollection.sort('timeDifference');
  var closestFeatureImage = sortedCollection.first(); // Extract the feature with the smallest time difference (the closest feature).
  // Retrieve specific properties from the closest feature:
  var tideValue = closestFeatureImage.get('Value'); 
  var dateImage = closestFeatureImage.get('dateImage');
  var dateGauge = closestFeatureImage.get('dateGauge');
  
  // Return an object containing the closest feature and its relevant properties.
  return {
    closestFeatureImage: closestFeatureImage,
    tideValue: tideValue,
    dateImage: dateImage,
    dateGauge: dateGauge
  };
}

///// MASK LAND
// This function masks land areas in the image by filtering out pixels with values above a specified threshold.
var waterLandFilter = function (image){
  return image.updateMask(image.lt(0.2))}; // Retain only pixels with values less than 0.2, masking all others (land areas).
  
///// FILTER NOISY PIXELS
// This function filters out small, noisy polygons based on the number of vertices.
// It ensures that only polygons with more than 20 vertices are retained.
// Needed to overcome the limitations of GEE
var filterPolygons = function(polygon) {
  var numVertices = ee.List(ee.List(polygon).get(0)).length(); // Calculate the number of vertices in the polygon.
  return ee.Algorithms.If(numVertices.gt(20), polygon, null); // Return the polygon if it has more than 20 vertices; otherwise, return null.
};


///// SELECT MAIN WATERLINE FEATURE
// This function identifies the main waterline feature by selecting the polygon with the highest number of vertices.
// Needed to overcome the noise produced by waterline partitions
var getPolygonWithMostVertices = function(feature) {
  var polygons = feature.geometry().geometries(); // Extract all polygons (geometries) from the input feature's geometry.
  
  var polygonsWithVertexCount = polygons.map(function(polygon) { // Map through each polygon to calculate the number of vertices and associate it as a property.
    var numVertices = ee.Geometry(polygon).coordinates().flatten().length(); // Calculate the number of vertices in the current polygon.
    return ee.Feature(ee.Geometry(polygon)).set('numVertices', numVertices); // Create a feature from the polygon and set the vertex count as a property.
  });
  var polygonsCollection = ee.FeatureCollection(polygonsWithVertexCount); // Convert the list of polygons with vertex counts into a FeatureCollection.
  var sortedPolygons = polygonsCollection.sort('numVertices', false); // Sort the polygons by the 'numVertices' property in descending order.
  var largestPolygon = sortedPolygons.first().geometry(); // Select the polygon with the highest number of vertices (the first in the sorted collection).
  
  return ee.Feature(largestPolygon); // Return the largest polygon as a Feature.
};  

///// EXTRACT VALUES FROM LIST OF POLYGONS
// This function processes a polygon by checking its size and returning either the first item in the list or the polygon itself.
// Needed to handle nested lists correctly on GEE
var processPolygon = function(polygon) {
      var size = ee.List(polygon).size(); // Get the size of the polygon list.
      // If the size of the list is greater than 1, return the first item; otherwise, return the polygon as-is.
      // Needed to overcome the limitations of GEE
      return ee.Algorithms.If(size.gt(1), [ee.List(polygon).get(0)], polygon);
};

///// CREATE FINAL DICTIONARY
// This function creates a dictionary-like result by adding custom properties to the feature.
// Needed for exporting SHP with properties
var createDictionaryResult = function(feature) {
        return ee.Feature(feature.geometry(), feature.toDictionary())
          .set('id', 'WL_' + datetime + '_' + ii)
          .set('dateTime', datetime)
          .set('tideValue', tideValue)
          .set('dateGauge', ee.Date(dateGauge).format('dd/MM/yyyy HH:mm:ss'))
          .set('dateImage', ee.Date(dateImage).format('dd/MM/yyyy HH:mm:ss')); 
};

//////////////////// FOR LEGEND
var colors = ['#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#FFA500', '#800080', '#008080', '#000000'];

// Create a legend panel
var legend = ui.Panel({
style: {
  position: 'bottom-right',
  padding: '8px 15px'
}
});

// Create a title for the legend
var legendTitle = ui.Label({
  value: 'Waterline Legend',
  style: {fontWeight: 'bold', fontSize: '18px', margin: '0 0 10px 0', padding: '0'}
});

legend.add(legendTitle);

// Define a function to add a color and label to the legend
var addColorAndLabel = function(color, label) {

  // Create a label
  var colorBox = ui.Label({
    style: {
      backgroundColor: color,
      padding: '8px',
      margin: '0 0 4px 0'
    }
  });

  // Create a panel with the color box and the label
  var legendEntry = ui.Panel({
    widgets: [colorBox, ui.Label(label)],
    layout: ui.Panel.Layout.Flow('horizontal'),
    style: {margin: '0 0 4px 0'}
  });

  legend.add(legendEntry);
};
//////////////////// END LEGEND


///// SOME IMPORTS
var pathTideg = 'projects/ee-cassie-intetidal/assets/'; // Replace by your asset path
var nameTideg = 'pamor_tide_sig_jun23_dez23_dhn'; // Replace by the name of tide file
var tideGauge = ee.FeatureCollection(pathTideg + nameTideg);
tideGauge = tideGauge.map(timeParser);

///// APPLY PRE PROCESSING FUNCTIONS
var imageCollectionWithDatetime = intertidal_zones.map(addDatetimeProperty); // Apply function to add datetime property with specified name and format
var numImages = imageCollectionWithDatetime.size().getInfo(); //get the size of image collection


///// BEGIN OF MAIN CODE ///// 
for (var i = 0; i < numImages; i++) { 
// Run over all images on the image collection
  var ndwi = ee.Image(imageCollectionWithDatetime.toList(numImages) // Convert imageCollection into list to interact over each
    .get(i)) // Get image according to counter "i"
    .select('ndwi'); // Retrieve the NDWI band for the current image
    
  // Compute histogram for the NDWI values within the specified geometry
  var histogram = ndwi.reduceRegion({
    reducer: ee.Reducer.histogram(),
    geometry: geometry, 
    scale: 10
  });
  var index = histogram.get('ndwi'); // Extract histogram data
  
  var threshold = otsu(histogram.get('histogram')); // Determine the threshold using the Otsu method
  
  // Calculate statistical measures of NDWI values
  var minPixelValue = ndwi.reduceRegion({reducer: ee.Reducer.min(),geometry: geometry,scale: 10, maxPixels: 1e9});
  var maxPixelValue = ndwi.reduceRegion({reducer: ee.Reducer.max(),geometry: geometry,scale: 10, maxPixels: 1e9});
  
  // Calculate standard deviation of NDWI values
  var stdDevPixelValue = ndwi.reduceRegion({reducer: ee.Reducer.stdDev(),geometry: geometry,scale: 10, maxPixels: 1e9});
  stdDevPixelValue = stdDevPixelValue.get('ndwi').getInfo();

  // Adjust the threshold  by the standard deviation value.
  // Needed to adjust the value according to each image to improve waterline definition
  var THajustment = ee.Number(stdDevPixelValue).add(1);
  var thresholdValue = ee.Number(threshold).subtract(ee.Number(THajustment).multiply(ee.Number(threshold))); 

  // Generate water-land image features by thresholding
  var waterLandImageFeature = ee.FeatureCollection(ndwi.clip(geometry).lt(thresholdValue)
      .reduceToVectors({ scale: scale, maxPixels: 1e12, eightConnected: false})
      .filter(ee.Filter.eq("label", 1))
  );
  
  // Filter the water-land image to retain the polygon with the most vertices
  var filteredWaterLandImageFeature = ee.Geometry.Polygon(ee.List(getPolygonWithMostVertices(waterLandImageFeature)
    .geometry().coordinates()).get(0));
  
  // Get the closest tide gauge data and corresponding properties, by comparing the datetime of them
  var result = getClosestFeatureImageAndValue(ndwi, tideGauge);
  var datetime = ee.Date(ndwi.get("datetime"));

  var closestFeatureImage = result.closestFeatureImage; // Defining which tide gauge value is the closest to the date of image acquisition
  var tideValue = result.tideValue; // Get the tide value of the closest tide gauge time series
  var dateImage = result.dateImage; // Get the datetime of the image
  var dateGauge = result.dateGauge; // Get the datetime of the tide gauge data

   // Format the datetime and tide data for display
  var formattedDatetime = datetime.format('YYYY-MM-dd HH:mm:ss').getInfo();
  var tideString = formattedDatetime + ' - ' + ee.String(ee.Number(tideValue).format('%.2f')).getInfo() + 'm';
  
  // Identify the waterline feature by applying function "identifyWaterFeature"
  var band = "ndwi";
  var waterSegment = identifyWaterFeature(ndwi, geometry, scale, band); 

  if(waterSegment !== null){
    // Process and filter polygons to remove noise                 
    var polygons = waterSegment.coordinates(); // Get only the coordinates of the water segment identified above
    var filteredPolygons = polygons.map(filterPolygons).removeAll([null]); // Remove polygons with less than 20 vertices
    
    waterSegment = ee.Geometry.MultiPolygon(filteredPolygons); // Transform into multpolygon again
    polygons = waterSegment.coordinates(); // Get the coordinate of the polygons
    var processedPolygons = polygons.map(processPolygon); // Needed to overcome GEE limitations to handle nested lists correctly

    // Create a new MultiPolygon with the processed polygons
    var waterSegment = ee.Geometry.MultiPolygon(processedPolygons); 

    var waterline = ee.FeatureCollection([]); 
    // Iterate over waterline segments to handle nested lists correctly
    for (var ii=0; ii < waterSegment.coordinates().size().getInfo(); ii++) {
      // Iterate through each set of coordinates in the water segment 
      var coordinates = waterSegment.coordinates().get(ii);
      
      // Create a MultiLineString from the coordinates and intersect it with the geometry
      // This ensures the waterline is limited to the specified region of interest
      var wl = ee.FeatureCollection(ee.Geometry.MultiLineString(coordinates)
        .intersection(geometry)); 
        
      // Verify that the waterline has valid coordinates
      var checkCoordinates = wl.first().geometry().coordinates(); 
      if (ee.List(checkCoordinates).size().gt(0).getInfo()) {
        
        // Check the type of the first coordinate to handle nested lists correctly
        var type2test = ee.Algorithms.ObjectType(ee.List(ee.List(checkCoordinates).get(0)).get(0));
          if (type2test.equals('List').getInfo()) {
            // Extract the last set of coordinates if the structure is nested
            var coordsList = ee.List(ee.List(checkCoordinates)).get(-1);
            wl = ee.FeatureCollection(ee.Geometry.MultiLineString(coordsList));
          } else {
        } 
      
      // Apply a Gaussian smoothing filter to refine the waterline coordinates
      wl = ee.FeatureCollection(ee.Geometry.MultiLineString(
        linearGaussianFilter((wl).geometry().coordinates()))); 
      
      // Add metadata to each waterline feature, such as tide and date information,
      // and merge it into the final waterline collection
      var updatedWl = wl.map(createDictionaryResult);  //createDictionaryResult
      waterline = waterline.merge(updatedWl);
      
      // Assign a unique color to this waterline based on the index of the image
      // Add a label with tide and datetime information for visualization purposes
      var color = colors[i % colors.length];
      addColorAndLabel(color, tideString); // Add color and datetime to the legend
    
      var geometry_wl = wl.geometry();
      Map.addLayer(geometry_wl,{'color': color },'waterline: '+ datetime.format('YYYY-MM-dd').getInfo());
      
      // Extract only essential properties (date and tide information) from the waterline
      // This reduces the size of the exported feature collection
      var smallFeatures = waterline.map(function(feature) {
        return feature.select(['dateGauge', 'dateImage', 'tideValue']); 
      });
      
      ///// EXPORT RESULTED WATERLINES
      Export.table.toDrive({
            collection: ee.FeatureCollection(smallFeatures),
            description: 'waterlines' + datetime.format('YYYY-MM-dd').getInfo(),
            fileFormat: 'SHP',
          });
          
      } else {
      } 
    }

  }else{
    print("Error image: "+ i);
  }
}
Map.add(legend);