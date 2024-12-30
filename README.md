# Waterline Detection and Exporting GEE Code 

This repository contains a script to detect waterlines using Water Indices, based on [Mason et al (1995)](https://www.researchgate.net/publication/253023601_Construction_of_an_inter-tidal_digital_elevation_model_by_the_'Water-Line'_Method) Waterline Method,from satellite imagery and export the results to Google Drive. The code leverages the TPXO tide signal to match waterline images with tide data for accurate topographic mapping. 
This README provides instructions on how to run the code for the example area of interest and modify it for other areas of interest.

---

## Table of Contents
- [Prerequisites](#prerequisites)
- [Running the Code](#running-the-code)
- [Customizing the Code](#customizing-the-code)
- [Troubleshooting](#troubleshooting)
- [License](#license)

---

## Prerequisites

Before running the code, ensure the following:

1. **Google Earth Engine (GEE) Account**: You must have a Google Earth Engine account. [Sign up here](https://signup.earthengine.google.com/).
2. **TPXO Asset**: The script uses a tide signal dataset from TPXO. The asset is available under the "input_asset" folder. Ensure the correct asset is imported from Google Earth Engine.
3. **Code Editor**: Use the Google Earth Engine code editor to run the script.
4. **Tide Signal Files**: Make sure you have the necessary assets in your Google Earth Engine account. You may need to adjust the paths to these files.

---

## Running the Code

### 1. **Import Assets**

The code requires an asset to run correctly. An example is located in the `input_asset` folder. Specifically, the script uses:

- **TPXO Tide Signal**: The TPXO data is required for tide information. It’s already included in the `input_asset` folder, but you must import it into your Earth Engine account. You can read about TPXO here: [TPXO Reference](https://www.tpxo.net/home).

- **Area Geometry**: This variable defines the region of interest where the waterline detection will occur. You don't need to change this to run the example; only modify if you want to test another area.

- **Satellite Imagery**: The script processes satellite imagery from an image collection. If you would like to test with another set of images, follow the instructions below to modify the image ID.

---

### 2. **Running the Exporting Task**

Once the code is loaded in the Google Earth Engine by pasting the provided "intertidal_topography_estimation" script from this repository as a repository on Code Editor, follow these steps:

1. **Import the necessary assets**:
   Ensure you’ve imported the correct tide signal file. The asset for the TPXO tide signal is available in the `input_asset` folder.

2. **Press Run**:
   Press "Run" in the Earth Engine Code Editor. This will execute the script and initiate the exporting task. The script will process the waterline detection and return the results. This may take a while. 

---

### 3. **View the Result on the Map**

After running the code, the resulting waterlines will be displayed on the map with the appropriate color-coded layers for each time step.

[Insert image of the Exporting Task]


### 4. **Export the Result to Google Drive**

Also, the results will be ready to be exported to your Google Drive as a shapefile (.SHP). Below is an example of the exporting task displayed in the Google Earth Engine console.
Once the task is complete (which may take a while), the waterline shapefiles can be downloaded to your local computer. They will include properties for the tide and image datetime, allowing you to check the reliability of the values, as well as the tide elevation value.

[Insert image of the Result on the Map]

---

## Customizing the Code

### 1. **Change the Area of Interest**

If you want to run the script for a different area, replace the geometry coordinates. To do this, you will need to:

1. Replace the `geometry` variable in the script with a new geometry coodrinates or import a new geometry asset.
    If you’re using an external geometry file, upload it to Google Earth Engine and import it using its asset ID, or define the geometry by its coordinates as follows:

    ```javascript
    var geometry = ee.Geometry.Polygon([[[...], [...], [...], [...]]]);  // Replace with your coordinates
    ```

2. Change Image Collection (Images ID)
    To use a different satellite image collection, create a new file on GEE Code Editor and run following lines:

    ```javascript
    var imageCollectionWithDatetime = ee.ImageCollection('COPERNICUS/S2_HARMONIZED')
        .filterDate('2022-01-01', '2022-12-31') // Choose your period of interest
        .filterBounds(geometry); // Address your geometry file (you must define the geometry first)
        .filter('CLOUDY_PIXEL_PERCENTAGE < 5') // Choose the acceptable cloud coverage percentage
    ```
    If you wish to test with another image collection, refer to the corresponding code in the folder for retrieving the image IDs. Ensure the images are processed and filtered with respect to your region of interest.

3. Change Tide Signal Asset
    To use a different tide signal file, modify the line that imports the TPXO tide signal by the asset path. For instance:

    ```javascript
    var tideGauge = ee.FeatureCollection('users/yourusername/tpxo_tide_signal'); // Replace 'users/yourusername/tpxo_tide_signal' with the correct asset ID for your tide signal data.
    ```

### Handling the File Headers:

You need to modify the file headers before importing the tide signal file so the code runs correctly. Ensure the following headers are used:
    - **"Date"** (format: dd/MM/yyyy)
    - **"Time"** (format: HH:mm:ss)  
    - **"Value"** (decimal separator as a dot)
    - Comma-Separated File (.CSV)

    If your elevation value is defined in meters, the resulting shapefile data will use the same unit. The reference level of the resulting data will be the same as the input tide signal data. Pay attention to this. If you need a different reference level, you must preprocess the tide signal file before importing it into the GEE asset or handle it in a GIS software afterward.

---

## Troubleshooting

- **Error with asset path**: Double-check the paths to the assets (geometry, tide signal and image IDs). You can import assets to Google Earth Engine via the Earth Engine Asset Manager.
- **No results on the map**: Ensure the region you are analyzing has adequate satellite image coverage and the geometry file correctly bounds the area of interest.

---

## License

The code is public, open, and free, provided under the Creative Commons CC-BY-SA license. The source must be referenced in the following format:

    CASSIE - accessed on [date] via the link: https://cassiengine.org/

---