---
title: 'Plasmid: ApE to edit and SnapGene Viewer to view'
author: Junli Zhang
date: '2024-03-13'
slug: plasmid-ape-to-edit-and-snapgene-viewer-to-view
categories:
  - Method
tags:
  - tutorial
---

Working on molecular cloning, I like to use [SnapGene Viewer](https://www.snapgene.com/snapgene-viewer) to view the plasmid map and features. It creates beautiful maps, and the layout of features and the enzyme is elegant. It always gives me a pleasant working experience. However, the free version has limited functions: view the sequences and add features. You cannot save features as a template and you cannot edit the sequences to make a new construct. So originally, I made the new construct in a text editor, then made a new file in SnapGene Viewer and re-annotated all the features. It is kind of a pain, especially when you find some mistake after all the annotations. Later I found the free and powerful software [ApE](https://jorgensen.biology.utah.edu/wayned/ape/). It takes me a little time to learn how to use it. The author of ApE made good tutorial videos on [YouTube](https://www.youtube.com/channel/UC_-pObWrnUZRhsO8YblX6gQ). Soon I know how to edit sequences, add new features, and annotate sequences. Because I am not satisfied with its map drawing, I still use SnapGene Viewer to view the final construct if you edited the .gb files with ApE. But there are some tricks before SnapGene Viewer can correctly display all the features. Here are my procedures to use SnapGene Viewer and ApE together.

1. I annotate new features in SnapGene Viewer, then export it to a .gbk file.
1. Open the .gbk file with ApE, select a feature or multiple features in the feature panel, then go to `Features -> Edit Feature Library -> New -> plasmid_name -> New Feature From Selected`
1. Then you can edit to change color etc.
1. You can also select a sequence and right-click to add a new feature.
1. For SnapGene Viewer to correctly show the color of the features, you also need to edit the feature file. To see where the library file is located, you can check `Features -> Open Feature Library`. From there, you can find the feature location. Go to the saving folder, open the file `Default_Features.txt`, and add a "`note {color: #ff9966}`" for example.
1. To annotate sequences, select part or the whole plasmid sequence, then go to `Features -> Annotate Features using Library` or `Ctrl + K` shortcut.
1. After you edit and save the sequence file, open it with a text editor, and replace all the `locus_tag` with `label`.
1. If the .gb or .gbk file is created from blank by ApE, please also add the reference keywords (author, title, and journal should have two space indent at the beginning).  


```perl
REFERENCE   1  (bases 1 to 100)
      AUTHORS   .
      TITLE     .
      JOURNAL   .
```

