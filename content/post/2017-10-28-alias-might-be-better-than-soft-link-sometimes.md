---
title: Alias might be better than soft link sometimes
author: Junli Zhang
date: '2017-10-28'
slug: alias-might-be-better-than-soft-link-sometimes
categories:
  - Linux
tags:
  - setting
---

I have been using Tassel and Haploview for a long time. Usually, I make a new folder under the folders of the two software for each project. However, this means that files in each project are distributed everywhere. Although I have a NOTE file in my project folder to record all the data and analysis locations, I still feel it is better to keep all files and analysis in the same project folder. To do this, I need let the command at least user-wide accessible.

To do this, I firstly make a soft link to the executable file of each software, but some software needs to load other files in the software folder. With a soft link, they cannot find other necessary files. Therefore, I made aliases for these files in the `.bashrc` file. This way I call the software with the full path and I can also give any options to the software, and I can give any names I like for them. For example, my alias for tassel is `alias tassel='~/Software/tassel-5-standalone/run_pipeline.pl'`, which is very easy to use and does not need to set PATH parameters.