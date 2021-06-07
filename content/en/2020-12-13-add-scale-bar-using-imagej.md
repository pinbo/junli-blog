---
title: Add scale bar using ImageJ
author: Junli Zhang
date: '2020-12-13'
slug: add-scale-bar-using-imagej
categories:
  - Method
tags:
  - Tips
  - tutorial
---

If you take a picture comparing the differences between a wild type plant and its mutant, you need to add a scale bar for publication. [ImageJ](https://imagej.nih.gov/ij/) is commonly used software to annotate images for scientific publication. Here is just a simple tutorial on how to add a scale bar to your image.

1. Open ImageJ and open your image.
2. Set the scale: you should include an object with known distance in your images, such as a ruler or just a tape. 
  1. Select the "Line" tool on the toolbar and draw a straight line along with the known distance.
  ![select length](/images/20201213-imageJ-select-line-tool.png)
  ![select length](/images/20201213-imageJ-select-distance.png)
  2. Go to Analyze -> Set Scale : change the "Known distance" (2cm here) and the "Unit" (cm here).
  ![set scale](/images/20201213-imageJ-set-scale.png)
3. Now the scale is set, you can measure any distance: just draw a line or select an area.
4. Add scale bar: go to Analyze -> Tools -> Scale bar: you can set the scale length (1cm here) and the location of the scale bar. If you want to save it as a TIFF image, you should also uncheck the "Overlay" option, otherwise, it will possibly not shown in other image viewers.
![add scale bar](/images/20201213-imageJ-add-scale-bar.png)

For more information, you can check ImageJ's documents: https://imagej.nih.gov/ij/docs/index.html